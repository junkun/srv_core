
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class kNN{
	private:
		double **matriz;
		int instancias;
		int atributos;
		vector <double> promedios;
		vector <double> desviaciones;
	public:
		kNN(){
		}
		
		kNN(int i, int a){
			instancias = i;
			atributos = a;
			matriz = new double*[instancias];
			for(int i = 0; i < instancias; i++){
				matriz[i] = new double[atributos + 1];
			}
		}

		kNN(string nomArchivo){
			int i, j;
			ifstream archivo(nomArchivo.c_str());

			archivo >> instancias;
			archivo >> atributos;
			matriz = new double*[instancias];
			for(i = 0; i < instancias; i++){
				matriz[i] = new double[atributos + 1];
			}

			for(i = 0; i < instancias; i++){
				for(j = 0; j < atributos; j++){
					archivo >> matriz[i][j];
				}
			}

			for(i = 0; i < instancias; i++){
				for(j = i; j < i + 10 && j < instancias; j++){
					matriz[j][atributos] = i / 10 + 1;
				}
			}
		}

		~kNN(){
			for(int i = 0; i < atributos; i++){
				delete[] matriz[i];
			}
			delete[] matriz;
		}

		void normalizarInstancias(){
			int i, j;
			double desv_est, media;

			for(i = 0 ; i < atributos; i++){
				media = 0;
				desv_est = 0;
				for(j = 0; j < instancias; j++){
					media += matriz[j][i];
				}
				media /= instancias;
				for(j = 0; j < instancias; j++){
					desv_est += (matriz[j][i] - media) * (matriz[j][i] - media);
				}
				desv_est = sqrt(desv_est / instancias);
				for(j = 0; j < instancias; j++){
					matriz[j][i] = (matriz[j][i] - media) / desv_est;
				}
				promedios.push_back(media);
				desviaciones.push_back(desv_est);
			}
		}

		void localizarInstancia(int k){
			int i, tam;
			double *instancia;
			ifstream archivo("Instancia.txt");

			archivo >> tam;
			instancia = new double[tam];
			for(i = 0; i < tam; i++){
				archivo >> instancia[i];
			}
			localizarInstancia(instancia, k);
			delete instancia;
		}

		int localizarInstancia(double *instancia, int k){
			int i, j;
			double **resultados;

			resultados = new double*[instancias];
			for(i = 0; i < instancias; i++){
				resultados[i] = new double[2];
			}

			for(i = 0; i < promedios.size(); i++){
				instancia[i] = (instancia[i] - promedios[i]) / desviaciones[i];
			}

			for(i = 0; i < instancias; i++){
				double dist = 0; //distancia euclidea
				for(j = 0; j < atributos; j++){
					dist += (matriz[i][j] - instancia[j]) * (matriz[i][j] - instancia[j]);
				}
				dist = sqrt(dist);
				resultados[i][0] = dist;
				resultados[i][1] = matriz[i][atributos];
			}

			ordenarQuicksort(resultados, 0, instancias - 1);

			for(i = 0; i < k && i < instancias; i++){
				cout << resultados[i][1] << " ";
			}
			cout << endl;

			for(i = 0; i < instancias; i++){
				delete[] resultados[i];
			}
			delete[] resultados;
		}

		void realizarPruebas(){
			int f, c, nMat, dist;
			double *instancia;
			ifstream archivo("InstanciasPrueba.txt");

			archivo >> dist;
			archivo >> nMat;

			for(int k = 0; k < nMat; k++){
				archivo >> f;
				archivo >> c;
				instancia = new double[c];
				for(int i = 0; i < f; i++){
					for(int j = 0; j < c; j++){
						archivo >> instancia[j];
					}
					localizarInstancia(instancia, dist);
				}
				cout << endl;
				delete[] instancia;
			}
			archivo.close();
		}
	private:
		void ordenarQuicksort(double **vector, int primero, int ultimo){
			int i = primero, j = ultimo;
			double pivote = vector[(primero + ultimo) / 2][0];
		 
			do{
				while(vector[i][0] < pivote){
					i++;
				}
				while(vector[j][0] > pivote){
					j--;
				}
				if (i <= j){
					double *auxiliar = vector[j];
				        vector[j] = vector[i];
				        vector[i] = auxiliar;
				        i++;
				        j--;
				}
		 
			}while(i <= j);
		 
			if(primero < j){
				ordenarQuicksort(vector,primero, j);
			}
			if(ultimo > i){
				ordenarQuicksort(vector,i, ultimo);
			}
		}
};

int main(int argc, char* argv[]){
	string archivo = argv[1];
	kNN knn(archivo);

	//knn.localizarInstancia(3);
	knn.realizarPruebas();

	return 0;
}
