#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <mpi.h>


using namespace std;


vector<vector<double> >
create_local_matrix(
	int n, 
	int m,
	int n_proc,
	int rank)
{
	int loc_n = rank < n % n_proc ? n / n_proc + 1 : n / n_proc;

	vector<vector<double> > F(loc_n, vector<double>(m));
	for (auto &i: F) {
		for (auto &j: i) {
			j = rand() % 1000;
		}
	}	

	return F;
}

vector<double> 
calc_max_in_col(const vector<vector<double> > &F)
{
	vector<double> F_max(F.size());
	for (int i = 0; i < F.size(); ++i) {
		F_max[i] = F[i][0];
		for (int j = 1; j < F[i].size(); ++j) {
			F_max[i] = max(F_max[i], F[i][j]);
		}
	}

	return F_max;
}


int
main(int argc, char *argv[])
{
	int n_proc;
	int rank;
	int n, m;



	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &n_proc);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);



//reading and sending n and m;
	if (rank == 0) {
		// this is main process
		cin >> n >> m;
		if (n < m) {
			swap(n, m);
		}
	}

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);

//Create and fill its part of matrices F and G
	srand(time(NULL));
	vector<vector<double> > F = create_local_matrix(n, m, n_proc, rank);
	vector<vector<double> > G = create_local_matrix(n, m, n_proc, rank);

//every process calculate max in every its rows of matrix G and columns of matrix F

	int loc_n = rank < n % n_proc ? n / n_proc + 1 : n / n_proc;

	//columns of F
	vector<double> loc_F_col_max(m);
	for (int i = 0; i < m; ++i) {
		loc_F_col_max[i] = F[0][i];
		for (int j = 1; j < loc_n; ++j) {
			loc_F_col_max[i] = max(loc_F_col_max[i], F[j][i]);
		}
	}

	//and rows of G
	vector<double> loc_G_row_max = calc_max_in_col(G);


//now we gather this local maximums in on common vector

	vector<double> F_max(m);
	MPI_Allreduce(loc_F_col_max.data(), F_max.data(), m, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

//find nash equilibriums

	vector<pair<int, int> > loc_ans;

	for (int i = 0; i < loc_n; ++i) {
		for (int j = 0; j < m; ++j) {
			if (F[i][j] == F_max[j] && G[i][j] == loc_G_row_max[i]) {
				loc_ans.push_back(make_pair(i, j));
			}
		}
	}


	MPI_Finalize();

}