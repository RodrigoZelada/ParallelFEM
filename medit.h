//void saveSolution(string name, int Nv, vec &uh){
template <typename T>
void saveSolution(string name, int Nv, T* uh){
  ofstream solution(name);
  solution << "MeshVersionFormatted 1 \n ";
  solution << "\n ";
  solution << "Dimension 2 \n ";
  solution << "\n ";
  solution << "SolAtVertices \n ";
  solution << Nv << " \n ";
  solution << "1 1 \n" ;
  for (int i=0; i<Nv; i++){
    solution << uh[i] << " \n "; //uh(i)
  }
  solution << " \n " ;
  solution << "End ";
}

void saveSolutionArma(string name, int Nv, vec &uh){
  ofstream solution(name);
  solution << "MeshVersionFormatted 1 \n ";
  solution << "\n ";
  solution << "Dimension 2 \n ";
  solution << "\n ";
  solution << "SolAtVertices \n ";
  solution << Nv << " \n ";
  solution << "1 1 \n" ;
  for (int i=0; i<Nv; i++){
    solution << uh(i) << " \n "; //uh(i)
  }
  solution << " \n " ;
  solution << "End ";
}

vec readSolution(string name, int Nv){
  ifstream file(name, ios::in);
  string line;
  vec uh(Nv, fill::zeros);
  int N;

  while (getline(file, line)) {
    if (line == "SolAtVertices") {
      getline(file, line);
      N = stoi(line);
      /*if (Nv != N){
        cout << "error " << endl;
        cout << "Nv = " << Nv << ", stoi(line) = " << N << endl;
      }*/
    }

    if (line == "1 1") {
      getline(file, line);

      for (int i = 0; i < N; i++) {
        file >> uh(i) ;
      }
    }
  }
  return uh;
}
