#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <mpi.h>


using namespace std;


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

vector<double>
my_gather(
	const std::vector<double> &loc_max,
	int n,
	int n_proc)
{
	int s_loc_n = n / n_proc;
	int cnt_big = n % n_proc;

	vector<int> pos(n_proc);
	for (int i = 0; i < n_proc; ++i) {
		pos[i] = i <= cnt_big ? i * s_loc_n : i * s_loc_n + (i - cnt_big);
	}

	vector<int> rcount(n_proc);
	for (int i = 0; i < n_proc; ++i) {
		rcount[i] = i < cnt_big ? s_loc_n + 1 : s_loc_n;
	}

	vector<double> comm_max(n);	
	MPI_Allgatherv(loc_max.data(), loc_max.size(), MPI_DOUBLE,
				comm_max.data(), rcount.data(), pos.data(), MPI_DOUBLE, MPI_COMM_WORLD);

	return comm_max;
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
	}

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);


//getting local n and m

	int n_b = n % n_proc;
	int m_b = m % n_proc;

	int loc_n = n / n_proc;
	if (rank < n % n_proc) {
		++loc_n;
	}

	int loc_m = m / n_proc;
	if (rank < m * n_proc) {
		++loc_m;
	}	

//Create and fill its part of matrices F and G
	vector<vector<double> > F(loc_m, vector<double>(n));
	for (auto &i: F) {
		for (auto &j: i) {
			j = 10.0;
		}
	}


	vector<vector<double> > G(loc_n, vector<double>(m));
	for (auto &i: G) {
		for (auto &j: i) {
			j = 15.0;
		}
	}

//every process calculate max in every its rows of matrix G and columns of matrix F

	//columns of F
	vector<double> loc_F_col_max = calc_max_in_col(F);
	//and rows of G
	vector<double> loc_G_row_max = calc_max_in_col(G);


//now we gather this local maximums in on common vector
	vector<double> F_max = my_gather(loc_F_col_max, n, n_proc);

	vector<double> G_max = my_gather(loc_G_row_max, m, n_proc);

//



	MPI_Finalize();

}