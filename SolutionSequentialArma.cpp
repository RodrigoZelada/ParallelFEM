// g++ SolutionSequential.cpp -o SolutionSequential -larmadillo
// time ./SolutionSequential

#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

double Frhs(double x, double y) { return 1.;}

#include "mesh.h"
#include "fespace.h"
#include "laplace2DSequantialArma.h"
#include "reduceAndReconstruct.h"
#include "medit.h"

int main(int argc, char * argv[]) {
	mesh Th("meshes/Th.mesh");

	int N=Th.nv;
	sp_mat AG(N,N);
	// cout << "AG = " << AG << endl;
	vec bG(N, fill::zeros);
	cout << "bG +1 = " << bG+1 << endl;
	vec F(N, fill::ones);
	laplace2D(AG, bG, Th, F);
	cout << "matrix assembled " << endl;
	vec uh = resolution(Th, AG, bG);
	saveSolutionArma("output/Tseq.sol", N, uh);

	return 0;
}
