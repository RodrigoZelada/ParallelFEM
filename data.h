void data(meshGamma &ThGamma, vec &TExact, vec &F, vec &f, vec &g, int hotregion, double (*FHot)(double, double), double (*FCold)(double, double), double (*THot)(double, double), double (*TCold)(double, double), double (*fgamma)(double, double), double (*ggamma)(double, double)){
  int ii, region;
  for (int n=0; n < ThGamma.nt; n++){
    for (int i=0; i < ThGamma.d+1; i++){
      ii = ThGamma.triangles[n][i];
      region = ThGamma.triangles[n][ThGamma.d+1];

      TExact(ii) = TCold(ThGamma.vertices[ii][0],ThGamma.vertices[ii][1]);
      F(ii) = FCold(ThGamma.vertices[ii][0],ThGamma.vertices[ii][1]);
      f(ii) = fgamma(ThGamma.vertices[ii][0],ThGamma.vertices[ii][1]);
      g(ii) = ggamma(ThGamma.vertices[ii][0],ThGamma.vertices[ii][1]);

      if (region==hotregion){
        TExact(ii) = THot(ThGamma.vertices[ii][0],ThGamma.vertices[ii][1]);
        F(ii) = FHot(ThGamma.vertices[ii][0],ThGamma.vertices[ii][1]);
      }


    }
  }
}
