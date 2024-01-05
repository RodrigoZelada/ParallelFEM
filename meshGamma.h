#include <bits/stdc++.h>

class meshGamma:public mesh {
public:
  int neGamma;
  int neGamma2;
  int nvGamma;
  int** Gamma;
  int** newGamma;
  int* VerticesGamma;
  int* VerticesGammaDuplicated;
  double* area;
  double*** verticesInTriangle;

  meshGamma(mesh &Th, int region_2, int label_interface, int dirichlet2[2]){
    duplicateNodesInterface(Th, region_2, label_interface, dirichlet2);
    computeArea();
  }

  ~meshGamma() {
    for (int i = 0; i < neGamma; i++) {
      delete Gamma[i];
    }
    for (int i = 0; i < neGamma2; i++) {
      delete newGamma[i];
    }
    delete[] VerticesGamma;
    delete[] VerticesGammaDuplicated;
    delete[] area;

    for (int k = 0; k<d+1; k++){
      for (int j = 0; j < d; j++){
        delete[] verticesInTriangle[k][j];
      }
      delete[] verticesInTriangle[k];
    }
	}

  void duplicateNodesInterface(mesh &Th, int region_2, int label_interface, int dirichlet2[2]){
    d = Th.d;
    nt = Th.nt;

    neGamma=0;
    for (int i = 0; i < Th.ne; i++){
      if (Th.edges[i][Th.d] == label_interface){
        neGamma++;
      }
    }

    nvGamma = neGamma + 1;

    Gamma = new int* [neGamma];
    for (int i = 0; i < neGamma; i++) {
      Gamma[i] = new int[Th.d-1];
    }

    neGamma = 0;
    for (int i = 0; i < Th.ne; i++){
      if (Th.edges[i][Th.d] == label_interface){
        Gamma[neGamma][0] = Th.edges[i][0];
        Gamma[neGamma][1] = Th.edges[i][1];
        neGamma++;
      }
    }

    VerticesGamma = new int[nvGamma];
    VerticesGammaDuplicated = new int[nvGamma];

    VerticesGamma[0] = Gamma[0][0];
    VerticesGammaDuplicated[0] = Th.nv;
    int k = 0;
    int nVG = 1;

    for (int i=0; i < neGamma; i++){
      for (int j = 0; j < Th.d; j++){
        k=0;
        while ((Gamma[i][j] != VerticesGamma[k]) && (k<nVG)){
          k++;
        }
        if (k==nVG){
          VerticesGamma[nVG] = Gamma[i][j];
          nVG++;
        }
      }
    }

    sort(VerticesGamma, VerticesGamma + nVG);
    for (int i=0; i<nVG; i++){
      VerticesGammaDuplicated[i] = Th.nv + i;
    }

    nv = Th.nv + nvGamma;
    vertices = new double* [nv];
    for (int i = 0; i < nv; i++) {
      vertices[i] = new double[Th.d+1];

      if (i < Th.nv){
        vertices[i][0] = Th.vertices[i][0];
        vertices[i][1] = Th.vertices[i][1];
        vertices[i][2] = Th.vertices[i][2];
      }
      else {
        vertices[i][0] = Th.vertices[VerticesGamma[i-Th.nv]][0];
        vertices[i][1] = Th.vertices[VerticesGamma[i-Th.nv]][1];
        vertices[i][2] = Th.vertices[VerticesGamma[i-Th.nv]][2];
      }
    }

    neGamma2=2*neGamma;
    newGamma = new int* [neGamma2];
    for (int i = 0; i < neGamma2; i++) {
      newGamma[i] = new int[Th.d-1];

      if (i < neGamma){
        for (int j = 0; j < Th.d; j++){
          newGamma[i][j] = Gamma[i][j];
        }
      }
    }

    ne = Th.ne + neGamma;
    edges = new int* [ne];
    for (int i = 0; i < ne; i++) {
      edges[i] = new int[Th.d + 1];

      if (i < Th.ne){

        for (int j = 0; j < Th.d+1; j++){
          edges[i][j] = Th.edges[i][j];
        }
      }
    }

    int ii, ii_tilde, jj, jj_tilde;
    int mi, mj;
    for (int l = 0; l<neGamma; l++){
      ii = Gamma[l][0];
      jj = Gamma[l][1];

      ii_tilde = -1;
      jj_tilde = -1;
      for (k=0; k<nvGamma; k++){
        if (VerticesGamma[k] == ii){ii_tilde = k; }
        if (VerticesGamma[k] == jj){jj_tilde = k; }
        //if ((ii_tilde >= 0) && (jj_tilde >= 0)){break;}
      }
      newGamma[neGamma+l][0] = Th.nv + ii_tilde;
      newGamma[neGamma+l][1] = Th.nv + jj_tilde;

      edges[Th.ne+l][0] = Th.nv + ii_tilde;
      edges[Th.ne+l][1] = Th.nv + jj_tilde;
      edges[Th.ne+l][2] = label_interface;

      //nodos de la interface conectados a un nodo con condición de borde dirichlet
      for (int m=0; m<ne; m++){
        if ((edges[m][0] == ii) || (edges[m][1] == ii)){
          mi = 0;
          if (edges[m][mi] == ii){
            mi=1;
          }
          if ((edges[m][Th.d] == dirichlet2[0]) || (edges[m][Th.d] == dirichlet2[1])){
              edges[m][mi] = Th.nv + ii_tilde;
          }
        }

        if ((edges[m][0] == jj) || (edges[m][1] == jj)){
          mj = 0;
          if (edges[m][mj] == jj){
            mj=1;
          }
          if ((edges[m][Th.d] == dirichlet2[0]) || (edges[m][Th.d] == dirichlet2[1])){
              edges[m][mj] = Th.nv + jj_tilde;
          }
        }
      }
    }

    triangles = new int* [Th.nt];
    for (int n = 0; n < Th.nt; n++) {
      triangles[n] = new int[Th.d + 2];

      for (int i=0; i<Th.d+2; i++){
        triangles[n][i] = Th.triangles[n][i];
      }
    }

    for (int n=0; n<Th.nt; n++){
      if (Th.triangles[n][Th.d+1] == region_2){
        for (int i=0; i<Th.d+1; i++){

          ii = Th.triangles[n][i];
          for (k=0; k<nvGamma; k++){
            if (VerticesGamma[k] == ii){
              ii_tilde = k;
              triangles[n][i] = Th.nv + ii_tilde;
            }
          }
        }
      }
    }

  };

  void computeArea(){
    mat Pe(3, 3);
    mat v(3,2);
    mat vones(3, 1, fill::ones);

    area = new double [nt];
    verticesInTriangle = new double** [d+1];

    for (int k = 0; k<d+1; k++){
      verticesInTriangle[k] = new double* [d];
      for (int j = 0; j < d; j++){
        verticesInTriangle[k][j] = new double [nt];
      }
    }

    for (int n = 0; n < nt; n++){
      for (int k = 0; k <= d; k++) {
        for (int j = 0; j < d; j++) {
          verticesInTriangle[k][j][n] = vertices[triangles[n][k]][j];
          v(k,j) = vertices[triangles[n][k]][j];
        }
      }
      Pe = join_rows(vones, v);
      area[n] = 0.5*abs(det(Pe));
    }
  };

};
