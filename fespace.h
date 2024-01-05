
//interpolation coefficients.
template<typename M>
mat baseSpace(M &Th, int degree){

	uvec index = {0,1,2};
	uvec l(2);
	mat Aii(3,3);
	vec oneAtVertex = {1,0,0};
	vec coefii(3);

	int ii, kk1, kk2;
	double xii, yii, xkk1, ykk1, xkk2, ykk2;

	mat coef(Th.nt, (Th.d+1)*3, fill::zeros);

	for (int n = 0; n < Th.nt; n++){
		for (int i = 0; i <= Th.d; i++){
			l = find(index != i); //the other 2 vertex of triangle T_n

			/*Aii = { {1, Th.verticesInTriangle[i][0][n],    Th.verticesInTriangle[i][1][n]},
							{1, Th.verticesInTriangle[l(0)][0][n], Th.verticesInTriangle[l(0)][1][n]},
							{1, Th.verticesInTriangle[l(1)][0][n], Th.verticesInTriangle[l(1)][1][n]} };*/
			Aii = {
						{1,Th.vertices[Th.triangles[n][i]][0],Th.vertices[Th.triangles[n][i]][1]},
						{1,Th.vertices[Th.triangles[n][l(0)]][0],Th.vertices[Th.triangles[n][l(0)]][1]},
						{1,Th.vertices[Th.triangles[n][l(1)]][0],Th.vertices[Th.triangles[n][l(1)]][1]}
			};

			coefii = solve(Aii, oneAtVertex);
			coef(n, span(3*i,3*i+2)) = coefii.t();
		}
	}
	return coef;
}
