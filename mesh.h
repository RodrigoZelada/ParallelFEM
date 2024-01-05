#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

class mesh {
public:
	int d; //dimension
	int nv; //number of vertices
	int nt; //number of triangles
	int ne; //number of edges

	double** vertices;
	int** edges;
	int** triangles;

	mesh(){};

	mesh(string meshFile) {
		loadMesh(meshFile);
	}

	~mesh() {
		for (int i = 0; i < nv; i++) {
			delete[] vertices[i];
		}
		for (int i = 0; i < ne; i++) {
			delete[] edges[i];
		}
		for (int i = 0; i < nt; i++) {
			delete[] triangles[i];
		}
	}

	template <typename T>
	//T **Array
	void readField(ifstream &file, T **Array, int n, int cols, T dif) {
		T tmp;
		string line;
		for (int i = 0; i < n; i++) {
			getline(file, line);
			istringstream iss(line);

			int j = 0;

			while (iss >> tmp) {
				if (j<cols-1){
					Array[i][j] = tmp-dif;
				}
				else {
					Array[i][j] = tmp;
				}

				j++;
			}
		}
	}

	void loadMesh(string meshFile) {

		ifstream file(meshFile, ios::in);
		string line;

		while (getline(file, line)) {
			if (line == "Dimension") {
				getline(file, line);
				d = stoi(line);
			}
			if (line == "Dimension 2") {
				d = 2;
			}

			if (line == "Vertices") {
				getline(file, line);
				nv = stoi(line);

				vertices = new double* [nv];
				for (int i = 0; i < nv; i++) {
					vertices[i] = new double[d+1];
				}

				readField<double>(file, vertices, nv, d+1, 0.);
			}
			if (line == "Edges") {
				getline(file, line);
				ne = stoi(line);

				edges = new int* [ne];
				for (int i = 0; i < ne; i++) {
					edges[i] = new int[d + 1];
				}

				readField<int>(file, edges, ne, d+1, 1);
			}
			if (line == "Triangles") {
				getline(file, line);
				nt = stoi(line);

				triangles = new int* [nt];
				for (int i = 0; i < nt; i++) {
					triangles[i] = new int[d + 2];
				}

				readField<int>(file, triangles, nt, d+2, 1);
			}
		}
	}
};
