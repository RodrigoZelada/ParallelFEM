void laplace2D(double* A, int* I, int* J, double* b, mesh &Th, vec &F, int &s)
{
	int Nv = Th.nv;
	int Nt = Th.nt;

	int D = Th.d; //dimension 2

	//Integers
	int n, k, d;
	int i, j, ii, jj, ii_global, jj_global, kk;

	mat Pe(3, 3);
	mat v(3,2);
	mat vones(3, 1, fill::ones);
	double area;

	vec coefii(3), coefjj(3);

	mat coef = baseSpace<mesh>(Th, 1);
	int boundary = 0;
	double large;

	s=0;

	//Iterate each triangle: T_n
	for (n = 0; n < Nt; n++){
		for (k = 0; k < 3; k++) {
			for (d = 0; d < 2; d++) {
				v(k,d) = Th.vertices[Th.triangles[n][k]][d];
			}
		}

		Pe = join_rows(vones, v);
		area = 0.5*det(Pe);

		//Iterate each vertex on triangle T_n: v_ii
		for (i = 0; i <= D; i++){
			ii = Th.triangles[n][i]; //index ii of vertex v_ii belonging to triangle T_n
			coefii = coef(n, span(3*i,3*i+2)).t();

			//b(ii_global) += (area/3)*dot(coefii,Pe.row(i))*F(ii); //(aii + bii*xii + cii*yii);
			b[ii] += (area/3)*dot(coefii,Pe.row(i))*F(ii);

			//Iterate a pair of vertex: edges
			for (j = 0; j <= D; j++){
				jj = Th.triangles[n][j];
				coefjj = coef(n, span(3*j,3*j+2)).t();

				//a(u,v) = grad(u)*grad(v), u = aii + bii*x + cii*y
				//A(ii,jj) += area*dot(coefii(span(1,2)),coefjj(span(1,2)));
				/*k=0;
				if (s>0){
					while ((k<s) && ((I[k] != ii) || (J[k] != jj))){
						k++;
					}
					if (k==s){
						I[s] = ii;
						J[s] = jj;
						A[s] += area*dot(coefii(span(1,2)),coefjj(span(1,2)));
						s++;
					}
					else{
						I[k] = ii;
						J[k] = jj;
						A[k] += area*dot(coefii(span(1,2)),coefjj(span(1,2)));
					}

				}
				else {
					I[s] = ii;
					J[s] = jj;
					A[s] += area*dot(coefii(span(1,2)),coefjj(span(1,2)));
					s++;
				}*/

			}
		}
	}
}
