// mpic++ Solution.cpp -o Solution -larmadillo
// mpirun -np 8 ./Solution

#include <armadillo>
#include <vector>
#include <mpi.h>

using namespace std;
using namespace arma;

double Frhs(double x, double y) { return 1.;}

#include "mesh.h"
#include "linalg.h"
#include "fespace.h"
#include "laplace2D.h"
#include "reduceAndReconstruct.h"
#include "medit.h"

int main(int argc, char * argv[]) {
	mesh ThGlobal("meshes/Th.mesh");
	int totalnodes, mynode;

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
	double* b = new double[N];
	double* F = new double[N_local];
	for (int i = 0; i<N; i++){
		b[i] = 0.;
	}
	for (int i = 0; i<N_local; i++){
		F[i] = 1.;
	}
	int* localToGlobalVertices = new int[N_local];
	for (int i=0;i<N_local;i++){
		int i_global = 0;
		while ( (Th.vertices[i][0] != ThGlobal.vertices[i_global][0]) || (Th.vertices[i][1] != ThGlobal.vertices[i_global][1])){
			i_global++;
		}
		localToGlobalVertices[i] = i_global;
	}
	laplace2D(As,Is,Js, b, Th, F, totalnodes, mynode, ThGlobal, s, localToGlobalVertices);
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

	int* count = new int[totalnodes];
	int* displacements = new int[totalnodes];
	MPI_Allgather(&s,1,MPI_INT,count,1,MPI_INT,MPI_COMM_WORLD);

	int total_s;
	MPI_Allreduce(&s,&total_s,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

	displacements[0] = 0;
	for (int i=1; i<totalnodes; i++){
		displacements[i] = 0;
		for (int j = 0; j < i; j++){
			displacements[i] += count[j];
		}
	}

	int* IGlobal = new int[total_s];
	int* JGlobal = new int[total_s];
	double* AGlobal = new double[total_s];
	double* bGlobal = new double[N];

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gatherv(A, s, MPI_DOUBLE, AGlobal, count, displacements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Gatherv(I, s, MPI_INT, IGlobal, count, displacements, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gatherv(J, s, MPI_INT, JGlobal, count, displacements, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Reduce(b,bGlobal,N,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

	delete[] count;
	delete[] displacements;
	delete[] I;
	delete[] J;
	delete[] A;
	delete[] b;

	if (mynode == 0){
		//sort
		int* index = new int[total_s];
		int* Isorted = new int[total_s];
		int* Jsorted = new int[total_s];
		double* Asorted = new double[total_s];

		for (int i=0; i<total_s; i++){index[i]=i;}
		sort(index, index+total_s, [&](const int& i, const int& j) {
													if (IGlobal[i] == IGlobal[j]){return (JGlobal[i]<JGlobal[j]);}
													else {return (IGlobal[i]<IGlobal[j]); } });
    for (int i=0; i<total_s; i++){
			Isorted[i] = IGlobal[index[i]];
			Jsorted[i] = JGlobal[index[i]];
			Asorted[i] = AGlobal[index[i]];
		}
		delete[] index;
		delete[] IGlobal;
		delete[] JGlobal;
		delete[] AGlobal;

		//sum duplicated index
		int* IRed = new int[total_s];
		int* JRed = new int[total_s];
		double* ARed = new double[total_s];

		//Sum duplicated index
		for (int k=0; k<total_s; k++){
			ARed[k] = 0.;
		}

		int n=0;
		IRed[0] = Isorted[0];
		JRed[0] = Jsorted[0];
		ARed[0] = Asorted[0];
		n++;

		for (int k=1; k<total_s; k++){
			if ((Isorted[k] == IRed[n-1]) && (Jsorted[k] == JRed[n-1]) ){ //si el elemento AGlobal[k] est+a en ARed
				IRed[n-1] = Isorted[k];
				JRed[n-1] = Jsorted[k];
				ARed[n-1]+= Asorted[k];
			}
			else {
				IRed[n] = Isorted[k];
				JRed[n] = Jsorted[k];
				ARed[n] = Asorted[k];
				n++;
			}
		}

		delete[] Isorted;
		delete[] Jsorted;
		delete[] Asorted;

		//Reduce matrix->Dirichlet boundary condition
		int* listEdgesMem = new int[ThGlobal.d*ThGlobal.ne];
		for (int i=0; i<ThGlobal.ne; i++){
			for (int j=0; j<ThGlobal.d; j++){
				listEdgesMem[i + ThGlobal.ne*j] = ThGlobal.edges[i][j];
			}
		}
		sort(listEdgesMem,listEdgesMem+ThGlobal.d*ThGlobal.ne);
		int le = 0;
    for (int i = 0; i < ThGlobal.d*ThGlobal.ne - 1; i++)
        if (listEdgesMem[i] != listEdgesMem[i + 1]){
            listEdgesMem[le++] = listEdgesMem[i];
				}
    listEdgesMem[le++] = listEdgesMem[ThGlobal.d*ThGlobal.ne - 1];

		int* listEdges = new int[le];
		for (int i=0; i<le; i++){
			listEdges[i] = listEdgesMem[i];
		}
		delete[] listEdgesMem;
		//le is the size of listEdges

		int Nred=N-le;
		double* bR = new double[Nred];
		int removed_rows=0;
		for (int i=0; i<N; i++){
			int k=0;
			while ((k<le) && (listEdges[k] != i)){
				k++;
			}
			if (k < le){
				removed_rows++;
			}
			else{//if i doesnt belong to listEdges, i add it to bR
				bR[i-removed_rows] = bGlobal[i];
			}
		}
		delete[] bGlobal;

		int nred=n-le;
		cout << "n = " << n << ", le = " << le << ", nred = " << nred << ", Nred = " << Nred << endl;
 		double* AR = new double[nred];
		int* IR = new int[nred];
		int* JR = new int[nred];
		int removed_cols, removed_k;
		removed_k = removed_rows = removed_cols = 0;
		for (int k = 0; k<n; k++){
			int i = IRed[k];
			int j = JRed[k];

			int kk=0;
			while ((kk<le) && (listEdges[kk] != i)){
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
				AR[k-removed_k] = ARed[k];
				nred=k+1-removed_k;//k+1-removed_k;
			}
		}
		delete[] ARed;
		delete[] IRed;
		delete[] JRed;

		//Conjugate gradient
		double* C = new double[nred]; //Preconditioning: diag(ARed)
		double* x = new double[Nred];
		double* grad = new double[Nred];
		double* Cg = new double[Nred];
		double* h = new double[Nred];
		double* Ah = new double[Nred];

		for (int i=0; i<Nred; i++){
			x[i] = bR[i];
		}
		grad = MultMatSparseVec(Nred, nred, AR, IR, JR, x);

		for (int i=0; i<Nred; i++){
			grad[i] -= bR[i];
			Cg[i] = h[i] = 0.;
		}

		for (int k=0; k<nred; k++){
	    int i = IR[k];
	    int j = JR[k];
			if (i==j){
		  	C[k] = 1./AR[k];
			}
			else{
				C[k] = 0.;
			}
	  }
		Cg = MultMatSparseVec(Nred, nred, C, IR, JR, grad);

		for (int i = 0; i<Nred; i++){
			h[i] = -Cg[i];
		}

		double g2 = dotLA(Nred, Cg, grad);

		double eps=1e-10;
		int nbitermax=1000;
		for (int iter = 0; iter < nbitermax; iter++){
			Ah = MultMatSparseVec(Nred, nred, AR, IR, JR, h);

			double ro = -dotLA(Nred,grad,h)/dotLA(Nred,h,Ah);

			for (int i=0; i<Nred; i++){
				x[i] += ro*h[i];
				grad[i] += ro*Ah[i];
			}
			Cg = MultMatSparseVec(Nred, nred, C, IR, JR, grad);
			double g2p = g2;
			g2 = dotLA(Nred,Cg,grad);
			if (g2 < eps){
	      cout << "iter = " << iter << ", ||g|| = " << sqrt(g2)  << endl;
				break;
	    }
			double gamma = g2/g2p;
			for (int i=0; i<Nred; i++){
				h[i] = -Cg[i] + h[i]*gamma;
			}
		}
		delete[] bR;
		delete[] IR;
		delete[] JR;
		delete[] AR;
		delete[] C;
		delete[] grad;
		delete[] Cg;
		delete[] h;
		delete[] Ah;

		removed_rows = 0;
		double* uh = new double[N];
		for (int i = 0; i<N; i++){
			uh[i] = 0.;
			int k=0;
			while ((k<le) && (listEdges[k] != i)){
				k++;
			}
			if (k < le){
				removed_rows++;
			}
			else{
				uh[i] = x[i-removed_rows];
			}
		}
		saveSolution("output/T.sol", N, uh);

		delete[] listEdges;
		delete[] x;
		delete[] uh;
	}
	else{
		delete[] IGlobal;
		delete[] JGlobal;
		delete[] AGlobal;
		delete[] bGlobal;
	}
	MPI_Finalize();

	return 0;
}
