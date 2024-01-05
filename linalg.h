double dotLA(int N, double *a, double *b){
  double sum = 0.0;

  for(int i=0;i<N;i++)
    sum += a[i]*b[i];

  return sum;
}

double* MultMatSparseVec(int N, int n, double *A, int* I, int *J, double *x){
  double* v = new double[N];
  for (int k=0; k<N; k++){
    v[k] = 0.;
  }

  for (int k=0; k<n; k++){
    int i = I[k];
    int j = J[k];
    v[i] += A[k]*x[j];
  }

  return v;
}
