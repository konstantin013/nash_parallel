#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <mpi.h>

void
get_local_dim()
{

}






void 
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
	}

	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);


//getting local n and m

	int n_b = n % n_proc;
	int m_b = m % n_proc;

	loc_n = n / n_proc;
	if (rank < n % n_proc) {
		++loc_n;
	}

	loc_m = m / n_proc;
	if (rank < m * n_proc) {
		++loc_m;
	}	

//Create and fill its part of matrices F and G
	vector<vector<double> > F(loc_m, vector<double>(n));
	for (auto &i: F) {
		for (auto &j: i) {
			j = 10.0
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
	vector<double> F_max(m);
	vector<int> pos(n_proc);
	for (int i = 0; i < n_proc; ++i) {
		pos[i] = i <= n_b ? (n / n_proc + 1) * i : ;
	} 

	vector<int> rcounts(n_proc);
	for (int i = 0; i < n_proc; ++i) {
		rcounts[i] = i < n_b ? n / n_proc + 1 : n / n_proc * n_b + (i - n_b) * (n / n_proc + 1);
	}
	MPI_Gatherv(loc_F_col_max.data(), loc_m, MPI_DOUBLE, F_max.data(), rcounts, pos, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(F_max.data(), m, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	vector<double> G_max(loc_n);

//








	if (rank == 0) {

	cout << "master hello" << fflush;


		int tag = 0;
		for (int i = 0; i < n; ++i) {
			int dest = i % (n_proc - 1);
			MPI_Send(&i, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
		}

		cout << "master send all assessments" << fflush;
/*
		tag = 1;
		for (int j = 0; j < m; ++j) {
			int dest = j % (n_proc - 1);
			MPI_Send(&j, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
		}
	*/

		for (int i = 1; i < n_proc; ++i) {
			int mes = -1;
			MPI_Send(&mes, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		}

	} else {
		vector<int> nums_rows;
		int mes;
		MPI_Status status;
		while (true) {
			MPI_Recv(&mes, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			if (mes == -1) {
				break;
			} else {
				nums_rows.push_back(mes);
			}
		}

		cout << "i am process " << rank << " i have got ";
		for (int i: nums_rows) {
			cout << i << ", ";
		}
		cout << endl;
	}



	MPI_Finalize();

}