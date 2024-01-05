// mpic++ SolutionFullParallel.cpp -o SolutionFullParallel -larmadillo
// mpirun -np 8 ./SolutionFullParallel

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
	delete[] localToGlobalVertices;
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

	double *bR = new double[N];
	double *AR = new double[total_s];
	int *IR = new int[total_s];
	int *JR = new int[total_s];
	int* sendcounts=new int[totalnodes];
	int* disp=new int[totalnodes];
	int* sendcountsSparse=new int[totalnodes];
	int* dispSparse=new int[totalnodes];
	int Nred, N_locallast, nred, n_local, N_localdisp, n_localdisp;

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

		Nred=N-le;
		bR = new double[Nred];
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

		nred=n-le;

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
				nred=k+1-removed_k;//k-removed_k;
			}
		}
		delete[] ARed;
		delete[] IRed;
		delete[] JRed;
		cout << "nred = " << nred << ", Nred = " << Nred << endl;

		//send b;
		N_local = (int) floor((double)Nred/totalnodes);
		N_locallast = Nred - N_local*(totalnodes-1); //41 - 5*7 = 6

		for (int i=0; i<totalnodes; i++){
			sendcounts[i] = N_local;
			disp[i] = i*N_local;
			if (i==totalnodes-1) {
				sendcounts[i]=N_locallast;
			}
		}

		int c=0;
		for (int i=0;i<totalnodes;i++){
			dispSparse[i] = sendcountsSparse[i] = 0;
		}
		for (int k=0;k<nred;k++){
			if ( (IR[k] >= disp[c]) && (IR[k] < disp[c]+sendcounts[c] ) ){
				sendcountsSparse[c]+=1;
			}
			else{
				sendcountsSparse[c+1]+=1;
				for (int i=c+1;i<totalnodes;i++){
					dispSparse[i] += sendcountsSparse[c];
				}
				c++;
			}
		}

		n_local = sendcountsSparse[0];
		N_localdisp = disp[0];
		n_localdisp = dispSparse[0];

		for (int i=1; i<totalnodes; i++){
			// MPI_Send(sendcounts,totalnodes,MPI_INT,i,0,MPI_COMM_WORLD);
			// MPI_Send(disp,totalnodes,MPI_INT,i,2,MPI_COMM_WORLD);
			// MPI_Send(sendcountsSparse+i,1,MPI_INT,i,1,MPI_COMM_WORLD);
			// MPI_Send(dispSparse+i,1,MPI_INT,i,3,MPI_COMM_WORLD);
			MPI_Send(&Nred,1,MPI_INT,i,4,MPI_COMM_WORLD);
			MPI_Send(&nred,1,MPI_INT,i,5,MPI_COMM_WORLD);
		}
	}
	else{
		delete[] IGlobal;
		delete[] JGlobal;
		delete[] AGlobal;
		delete[] bGlobal;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if (mynode!=0){
		// MPI_Recv(&N_local,1,MPI_INT,0,0,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// MPI_Recv(&n_local,1,MPI_INT,0,1,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// MPI_Recv(&N_localdisp,1,MPI_INT,0,2,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		// MPI_Recv(&n_localdisp,1,MPI_INT,0,3,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&Nred,1,MPI_INT,0,4,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		MPI_Recv(&nred,1,MPI_INT,0,5,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	MPI_Bcast(sendcounts,totalnodes,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(sendcountsSparse,totalnodes,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(disp,totalnodes,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(dispSparse,totalnodes,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(bR,Nred,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(IR,nred,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(JR,nred,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(AR,nred,MPI_DOUBLE,0,MPI_COMM_WORLD);

	N_local = sendcounts[mynode];
	n_local = sendcountsSparse[mynode];
	N_localdisp = disp[mynode];
	n_localdisp = dispSparse[mynode];

	delete[] sendcountsSparse;
	delete[] dispSparse;

	double* bRLocal = new double[N_local];
	int* IRLocal = new int[n_local];
	int* JRLocal = new int[n_local];
	double* ARLocal = new double[n_local];

	for (int i=0; i<N_local; i++){
		bRLocal[i] = bR[i+N_localdisp];
	}

	for (int i=0; i<n_local; i++){
		IRLocal[i] = IR[i+n_localdisp]-N_localdisp;
		JRLocal[i] = JR[i+n_localdisp]-N_localdisp;
		ARLocal[i] = AR[i+n_localdisp];
	}

	delete[] IR;
	delete[] JR;
	delete[] AR;

	//Parallel conjugate gradient
	double *p = new double[Nred]; //direction
  double *z = new double[N_local]; //z = A*p
  double *r = new double[N_local]; //residual r -= alpha*z;
  double *mr = new double[N_local];
	double *x = new double[N_local];
	double* M = new double[N_local]; //Preconditioning: diag(ARed)

	for(int i=0;i<Nred;i++){
		p[i] = 0.; //bR[i];
	}
	delete[] bR;

	for (int k=0; k<n_local; k++){
		if (IRLocal[k] == JRLocal[k]){
			M[IRLocal[k]] = 1./ARLocal[k];
		}

	}
	for(int i=0;i<N_local;i++){
		x[i] = 0.; //bRLocal[i]; //1.;
    r[i] = bRLocal[i];  //calculation of residual
  }
	for (int k=0; k<n_local; k++){
		int j_global = JRLocal[k]+N_localdisp;
		r[IRLocal[k]] -= ARLocal[k]*p[j_global]; //dot(A[i],p)
	}
	for(int i=0;i<N_local;i++){
    mr[i] = M[i]*r[i];   //calculation of modified residual: M the diagonal of A
  }
	double sum,local_sum,c,d,alpha,beta;
	local_sum = dotLA(N_local,mr,r);
  MPI_Allreduce(&local_sum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  c = sum;

	cout << "sum = " << sum << endl;

	//I gather each mr in the different processors and i send it to every processor
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allgatherv(mr,N_local,MPI_DOUBLE,p,sendcounts,disp,MPI_DOUBLE,MPI_COMM_WORLD); //p0 = mr0

	double tol=1e-10;
	int nbitermax=1000;
	for (int iter = 0; iter < nbitermax; iter++){
		//z = A*p
		for (int i=0;i<N_local;i++){
			z[i] = 0.;
		}
		for (int k=0; k<n_local; k++){
			int j_global = JRLocal[k]+N_localdisp;
			z[IRLocal[k]] -= ARLocal[k]*p[j_global]; //dot(A[i],p)
		}

		//parallel sum = dot(z,p)
		local_sum = 0;
		for (int i=0; i<N_local; i++){
			local_sum += z[i]*p[i+N_localdisp];
		}
		MPI_Allreduce(&local_sum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

		//alpha = c/dot(p,z)
    alpha = c/sum;

		//parallel x = x + alpha*p and r = r - alpha*z
    for(int i=0;i<N_local;i++){
      x[i] += alpha*p[i+N_localdisp];
      r[i] -= alpha*z[i];
    }

		//Preconditioning Stage
		//supposing M = Diag(A), M mr = r -> mr[i] = r[i] / diag(A)[i] = r[i]*M[i]
		for(int i=0;i<N_local;i++){
      mr[i] = M[i]*r[i];
		}

		// parallel sum = dot(mr,r)
		local_sum = dotLA(N_local,mr,r);
		MPI_Allreduce(&local_sum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

		//d = dot(mr,r)
		d = sum; //contains inner product of residual and modified residual

    if(fabs(d) < tol) {
			if (mynode==0) cout << "mynode = " << mynode << ", iter = " << iter << "\t" << "dot(mr,r) = " << d << endl;
			break;
		}
		//beta = dot(mr,mr)/dot(old mr, old mr) //old mr: mr of previous iteration
    beta = d/c;

		//parallel p = mr + beta*p,
		for(int i=0;i<N_local;i++){
      z[i] = mr[i] + beta*p[i+N_localdisp];
		}
    //i use z because has the dimensions of local rows, then i gathered it in p
		MPI_Allgatherv(z,N_local,MPI_DOUBLE,p,sendcounts,disp,MPI_DOUBLE,MPI_COMM_WORLD); //p0 = mr0

		//save dot(old mr, old mr)
    c = d;

	}

	delete[] bRLocal;
	delete[] IRLocal;
	delete[] JRLocal;
	delete[] ARLocal;
	delete[] z;
	delete[] r;
	delete[] mr;
	delete[] M;
	MPI_Gatherv(x,N_local,MPI_DOUBLE,p,sendcounts,disp,MPI_DOUBLE,0,MPI_COMM_WORLD); //p0 = mr0
	delete[] sendcounts;
	delete[] disp;
	delete[] x;

	if(mynode==0){
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
		int removed_rows = 0;
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
				uh[i] = -p[i-removed_rows];
			}
		}
		saveSolution("output/TFullParallel.sol", N, uh);

		delete[] listEdges;
		delete[] uh;
	}
	delete[] p;

	MPI_Finalize();

	return 0;
}
