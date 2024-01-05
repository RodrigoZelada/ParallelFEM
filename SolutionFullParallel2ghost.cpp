// mpic++ SolutionFullParallel2.cpp -o SolutionFullParallel2 -larmadillo
// mpirun -np 8 ./SolutionFullParallel2

#include <armadillo>
#include <vector>
#include <mpi.h>
#include <cmath>

using namespace std;
using namespace arma;

double Frhs(double x, double y) { return 1.;}

#include "mesh.h"
#include "linalg.h"
#include "fespace.h"
#include "laplace2DFullParallel.h"
#include "reduceAndReconstruct.h"
#include "medit.h"

int main(int argc, char * argv[]) {
	mesh ThGlobal("meshes/Th.mesh");
	int totalnodes, mynode;
	const double tol=1e-10;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	mesh Th("meshes/Th"+to_string(mynode)+".mesh");
	cout << "meshes/Th"+to_string(mynode)+".mesh" << endl;

	int N=ThGlobal.nv;
	int N_local=Th.nv;
	int D = Th.d;
	int s=0;
	for (int n = 0; n < N_local; n++){
		//Iterate each vertex on triangle T_n: v_ii
		for (int i = 0; i <= D; i++){
			//Iterate a pair of vertex: edges
			for (int j = 0; j <= D; j++){
				s++;
			}
		}
	}

	int* Is = new int[s];
	int* Js = new int[s];
	double* As = new double[s];

	for (int k = 0; k<s; k++){
		Is[k]=Js[k]=As[k]=0;
	}
	//We define b of dimension N instead of N_local because in that way we just have to reduce
	//We cannot do the same for A because we have the sparse structure
	double* b = new double[N_local];
	double* F = new double[N_local];
	for (int i = 0; i<N_local; i++){
		b[i] = 0.;
	}
	for (int i = 0; i<N_local; i++){
		F[i] = 1.;
	}

	laplace2DFullParallel(As,Is,Js, b, Th, F, totalnodes, mynode, ThGlobal, s);
	delete[] F;
	//Eliminate the extra memory used in the dynamic allocation, creating a new array
	int* I = new int[s];
	int* J = new int[s];
	double* A = new double[s];
	for (int k = 0; k<s; k++){
		I[k] = Is[k];
		J[k] = Js[k];
		A[k] = As[k];
	}
	delete[] Is;
	delete[] Js;
	delete[] As;

	//Reduce matrix
	int Nghost=0;
	int* listEdgesMem = new int[Th.d*Th.ne];
	for (int i=0; i<Th.ne; i++){
		for (int j=0; j<Th.d; j++){
			listEdgesMem[i + Th.ne*j] = Th.edges[i][j];
			if (Th.edges[i][Th.d] >= 0){Nghost++;}
			// if (Th.edges[i][Th.d] >= 0){
			// 	listEdgesMem[i + Th.ne*j] = Th.edges[i][j];
			// }
			// else{
			// 	listEdgesMem[i + Th.ne*j] = Th.edges[0][0];
			// }
		}
	}
	sort(listEdgesMem,listEdgesMem+Th.d*Th.ne);
	int le = 0;
	for (int i = 0; i < Th.d*Th.ne - 1; i++){
			if (listEdgesMem[i] != listEdgesMem[i + 1]){
					listEdgesMem[le++] = listEdgesMem[i];
			}
	}
	listEdgesMem[le++] = listEdgesMem[Th.d*Th.ne - 1];

	int* listEdges = new int[le];
	for (int i=0; i<le; i++){
		listEdges[i] = listEdgesMem[i];
	}

	delete[] listEdgesMem;


	//bijection between local to global mesh
	int* localToGlobalVertices = new int[N_local];
	for (int i=0;i<N_local;i++){
		int i_global = 0;
		while ( (Th.vertices[i][0] != ThGlobal.vertices[i_global][0]) || (Th.vertices[i][1] != ThGlobal.vertices[i_global][1])){
			i_global++;
		}
		localToGlobalVertices[i] = i_global;
	}

	int Nred=N_local-le;
	double* bR = new double[Nred];
	int* localToGlobalVerticesR = new int[Nred];
	int removed_rows=0;

	for (int i=0; i<N_local; i++){
		int k=0;
		while ((k<le) && (listEdges[k] != i)){
			k++;
		}
		if (k < le){
			removed_rows++;
		}
		else{//if i doesnt belong to listEdges, i add it to bR
	 		bR[i-removed_rows] = b[i];
			localToGlobalVerticesR[i-removed_rows] = localToGlobalVertices[i];
		}
	}
	delete[] b;
	delete[] localToGlobalVertices;

	int* index = new int[s];
	int* Isorted = new int[s];
	int* Jsorted = new int[s];
	double* Asorted = new double[s];

	for (int i=0; i<s; i++){index[i]=i;}
	sort(index, index+s, [&](const int& i, const int& j) {
												if (I[i] == I[j]){return (J[i]<J[j]);}
												else {return (I[i]<I[j]); } });
	for (int i=0; i<s; i++){
		Isorted[i] = I[index[i]];
		Jsorted[i] = J[index[i]];
		Asorted[i] = A[index[i]];
	}
	delete[] index;
	delete[] I;
	delete[] J;
	delete[] A;

	int nred=s-le;
	double* AR = new double[nred];
	int* IR = new int[nred];
	int* JR = new int[nred];
	int removed_cols, removed_k;
	removed_k = removed_rows = removed_cols = 0;

	for (int k = 0; k<s; k++){
		int i = Isorted[k];
		int j = Jsorted[k];

		int kk=0;
		while ((kk<le) && (listEdges[kk] != i) ){
			kk++;
		}
		int kkk=0;
		while ((kkk<le) && (listEdges[kkk] != j)){
			kkk++;
		}
		if ( (kk < le ) || (kkk < le) ){
			removed_k++;
		}
		else {
			kk=kkk=0;
			while ((kk<le) && (i > listEdges[kk] )){
				kk++;
			}
			while ((kkk<le) && (j > listEdges[kkk] )){
				kkk++;
			}
			IR[k-removed_k] = i-kk;
			JR[k-removed_k] = j-kkk;
			AR[k-removed_k] = Asorted[k];
			nred=k+1-removed_k;
		}
	}
	delete[] Asorted;
	delete[] Isorted;
	delete[] Jsorted;

	//Parallel conjugate gradient
	int* count = new int[totalnodes];
	int* displacements = new int[totalnodes];

	int total_s;
	MPI_Allreduce(&nred,&total_s,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	int N_global;
	MPI_Allreduce(&Nred,&N_global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

	//count and displacements
	MPI_Allgather(&Nred,1,MPI_INT,count,1,MPI_INT,MPI_COMM_WORLD);
	displacements[0] = 0;
	for (int i=1; i<totalnodes; i++){
		displacements[i] = 0;
		for (int j = 0; j < i; j++){
			displacements[i] += count[j];
		}
	}

	double *p = new double[N_global]; //direction
  double *z = new double[Nred]; //z = A*p
  double *r = new double[Nred]; //residual r -= alpha*z;
  double *mr = new double[Nred];

	double *x = new double[Nred];

	double* M = new double[Nred]; //Preconditioning: diag(ARed)

	for(int i=0;i<N_global;i++){
		p[i] = 0.;
	}

	for (int k=0; k<nred; k++){
		if (IR[k] == JR[k]){
			M[IR[k]] = 1./AR[k];
		}

	}
	for(int i=0;i<Nred;i++){
		x[i] = 0.; //1.;
    r[i] = bR[i];  //calculation of residual
  }

	//2 options: 1. Use local index in I, J, A
	//           2. Use size N in every vector
	for (int k=0; k<nred; k++){
		//j global is the error, i never reach N_Global
		int j_global = JR[k] + displacements[mynode];
		r[IR[k]] -= AR[k]*p[j_global]; //dot(A[i],p)
	}

	for(int i=0;i<Nred;i++){
		mr[i] = M[i]*r[i];   //calculation of modified residual: M the diagonal of A
	}

	double sum,local_sum,c,d,alpha,beta;
	local_sum = dotLA(Nred,mr,r);
  MPI_Allreduce(&local_sum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  c = sum;

	//I gather each mr in the different processors and i send it to every processor
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allgatherv(mr,Nred,MPI_DOUBLE,p,count,displacements,MPI_DOUBLE,MPI_COMM_WORLD); //p0 = mr0

	int nbitermax=1000;
	for (int iter = 0; iter < nbitermax; iter++){

		//z = A*p
		for (int i=0;i<Nred;i++){
			z[i] = 0.;
		}
		for (int k=0; k<nred; k++){
			int j_global = JR[k] + displacements[mynode];
			z[IR[k]] -= AR[k]*p[j_global]; //dot(A[i],p)
		}

		//parallel sum = dot(z,p)
		local_sum = 0;
		for (int i=0; i<Nred; i++){
			int i_global =  i + displacements[mynode];
			local_sum += z[i]*p[i_global];
		}
		MPI_Allreduce(&local_sum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

		//alpha = c/dot(p,z)
    alpha = c/sum;

		//parallel x = x + alpha*p and r = r - alpha*z
    for(int i=0;i<Nred;i++){
			int i_global =  i + displacements[mynode];
      x[i] += alpha*p[i_global];
      r[i] -= alpha*z[i];
    }

		//Preconditioning Stage
		//supposing M = Diag(A), M mr = r -> mr[i] = r[i] / diag(A)[i] = r[i]*M[i]
		for(int i=0;i<Nred;i++){
      mr[i] = M[i]*r[i];
		}

		// parallel sum = dot(mr,r)
		local_sum = dotLA(Nred,mr,r);
		MPI_Allreduce(&local_sum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

		//d = dot(mr,r)
		d = sum; //contains inner product of residual and modified residual

		//parallel sum = dot(r,r)
		local_sum = dotLA(Nred,r,r); //parallel sum = ||r||^2
    MPI_Allreduce(&local_sum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

		if(fabs(d) < tol) {
			if (mynode==0) cout << "mynode = " << mynode << ", iter = " << iter << "\t" << "dot(mr,r) = " << d << endl;
			break;
		}

		//beta = dot(mr,mr)/dot(old mr, old mr) //old mr: mr of previous iteration
    beta = d/c;

		//parallel p = mr + beta*p,
		for(int i=0;i<Nred;i++){
			int i_global =  i + displacements[mynode];
      z[i] = mr[i] + beta*p[i_global];
		}
    //i use z because has the dimensions of local rows, then i gathered it in p
		MPI_Allgatherv(z,Nred,MPI_DOUBLE,p,count,displacements,MPI_DOUBLE,MPI_COMM_WORLD); //p0 = mr0

		//save dot(old mr, old mr)
    c = d;
	}

	removed_rows = 0;
	double* uh = new double[N_local];
	for (int i = 0; i<N_local; i++){
		uh[i] = 0.;
		int k=0;
		while ((k<le) && (listEdges[k] != i) ){
			k++;
		}
		if (k < le){
			removed_rows++;
		}
		else{
			uh[i] = -x[i-removed_rows];
		}
	}
	cout << "output/T"+to_string(mynode)+".sol" << endl;
	saveSolution("output/T"+to_string(mynode)+".sol", N_local, uh);

	// int* indexGlobal = new int[N_global];
	// MPI_Gatherv(x,Nred,MPI_DOUBLE,p,count,displacements,MPI_DOUBLE,0,MPI_COMM_WORLD); //p0 = mr0
	// MPI_Gatherv(localToGlobalVerticesR,Nred,MPI_INT,indexGlobal,count,displacements,MPI_INT,0,MPI_COMM_WORLD); //p0 = mr0

	// if (mynode=0){
	// 	double* uhred = new double[N_global];
	// 	sort(indexGlobal, indexGlobal+N_global, [&](const int& i, const int& j) {return (indexGlobal[i]<indexGlobal[j]);});
	// 	for (int i=0; i<N_global; i++){
	// 		uhred[i] = p[indexGlobal[i]];
	// 	}
	// }

	MPI_Finalize();

	return 0;
}
