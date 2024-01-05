// g++ SolutionSequential.cpp -o SolutionSequential -larmadillo
// time ./SolutionSequential
//LA ESTRUCTURA DE MATRICES SPARSES QUE UTILICE ES MUY LENTA!
//LAS MATRICES SPARSES NO SE DEFINEN CON DYNAMIC ARRAYS, PUES HAY QUE RESERVAR MEMORIA FIJA
//SE HACE CON "UNORDER_MAP" QUE SON COMO DICCIONARIOS, PARA INSERTAR CADA VEZ UN NUEVO INDICE.

#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

double Frhs(double x, double y) { return 1.;}

#include "mesh.h"
#include "fespace.h"
#include "laplace2DSequantial.h"
#include "reduceAndReconstruct.h"
#include "medit.h"

int main(int argc, char * argv[]) {
	mesh Th("meshes/Th.mesh");

	int N=Th.nv;
	int D = Th.d;
	int s=0;
	for (int n = 0; n < N; n++){
		//Iterate each vertex on triangle T_n: v_ii
		for (int i = 0; i <= D; i++){
			//Iterate a pair of vertex: edges
			for (int j = 0; j <= D; j++){
				s++;
			}
		}
	}

	double* A = new double[s];
	int* I = new int[s];
	int* J = new int[s];

	for (int i = 0; i<s; i++){
		I[i]=J[i]=-1;
		A[i]=0;
	}
	double* b = new double[N];
	for (int i = 0; i<N; i++){
		b[i] = 0.;
	}
	vec F(N, fill::ones);

	laplace2D(A,I,J, b, Th, F,s);

	cout << "matrix assembled " << endl;

	// sp_mat AG(N,N);
	// vec bG(N, fill::zeros);
	// for (int k=0; k<s; k++){
	// 	AG(I[k],J[k]) = A[k];
	// }
	// for (int i=0; i<N; i++){
	// 	bG(i) = b[i];
	// }
	// vec uh = resolution(Th, AG, bG);
	// saveSolutionArma("output/Tseq.sol", N, uh);
	delete[] A;
	delete[] I;
	delete[] J;

	return 0;
}
