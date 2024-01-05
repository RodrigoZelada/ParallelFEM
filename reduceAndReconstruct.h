template<typename M>
vec resolution(M &Th, sp_mat &A, vec &b){

  umat matEdges(Th.ne,2) ;
  for (int i = 0; i < Th.ne; i++){
    for (int j = 0; j < 2; j++){
      matEdges(i,j) = Th.edges[i][j];
    }
  }

  uvec edgesU = unique(matEdges);
  b.shed_rows(edgesU);
  int acr = 0;
  int acc = 0;
  for (int i = 0; i < edgesU.n_rows; i++){
    A.shed_row(edgesU(i) - acr);
    acr ++;
    A.shed_col(edgesU(i) - acc);
    acc ++;
  }

  vec uhred = spsolve(A,b);
  uvec not_boundary = regspace<uvec>(0,Th.nv-1);
  not_boundary.shed_rows(edgesU);

  vec uh(Th.nv, fill::zeros);
  uh.elem(not_boundary) = uhred;

  return uh;
}
