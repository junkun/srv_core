#include <stdlib.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <vector>

#include "eigen/Eigen/Core"
#include "eigen/Eigen/Eigenvalues"

using namespace std;
using namespace Eigen;

class PCA{
	private:
		MatrixXd matriz;
		MatrixXd correlaciones;
		MatrixXd eigenVectores;
		vector <double> promedios;
		vector <double> desviaciones;
	public:
		PCA(){
		}

		PCA(int f, int c){
			matriz = MatrixXd::Zero(f, c);
		}

		PCA(string nomArchivo, int modo){
			int i, j, f, c;
			ifstream archivo(nomArchivo.c_str());

			archivo >> f;
			archivo >> c;
			
			if(modo == 0){ //cargar matriz a reducir
				matriz = MatrixXd::Zero(f, c);
				for(i = 0; i < f; i++){
					for(j = 0; j < c; j++){
						archivo >> matriz(i, j);
					}
				}
			}else if(modo == 1){ //cargar matriz eigenVectores y vectores de promedios y desviaciones
				eigenVectores = MatrixXd::Zero(f, c);
				for(i = 0; i < f; i++){
					for(j = 0; j < c; j++){
						archivo >> eigenVectores(i, j);
					}
				}
				
				double v;
				for(i = 0; i < f; i++){	
					archivo >> v;
					promedios.push_back(v);
				}
				for(i = 0; i < f; i++){	
					archivo >> v;
					desviaciones.push_back(v);
				}
			}
			archivo.close();
		}

		~PCA(){
		}

		MatrixXd reducir(double pctObjetivo){
			normalizar(); //normaliza matriz de entrada
			correlacion();
			eigenvalores(pctObjetivo);

			MatrixXd reducida, aux = eigenVectores, aux2 = matriz; // Obtener la matriz de eigenvectores ordenados y reducidos
			aux.transposeInPlace(); // La operación para reducir las dimensiones: CP = EigenvectoresReducidos' x MatrizNormalizada'
			aux2.transposeInPlace();
			reducida = aux * aux2;
			reducida.transposeInPlace();

			return reducida;
		}

		MatrixXd transformarVector(){
			int c;
			ifstream archivo("vEntrada.txt");

			archivo >> c;
			MatrixXd v = MatrixXd::Zero(1, c);
			for(int i = 0; i < c; i++){
				archivo >> v(0, i);
			}

			return transformarVector(v);
		}

		MatrixXd transformarVector(MatrixXd v){
			for(int i = 0; i < promedios.size(); i++){
				v(0, i) = (v(0, i) - promedios[i]) / desviaciones[i];
			}

			MatrixXd res, aux = eigenVectores;
			aux.transposeInPlace();
			v.transposeInPlace();
			res = aux * v;
			res.transposeInPlace();
			return res;
		}

		void imprimirArchivo(string nomArchivo){
			ofstream archivo(nomArchivo.c_str());
			archivo << eigenVectores.rows() << " " << eigenVectores.cols() << endl;
			archivo << eigenVectores << endl << endl;
			
			int i = 0;
			for(i = 0; i < promedios.size(); i++){
				archivo << promedios[i] << " ";
			}
			archivo << endl;
			for(i = 0; i < desviaciones.size(); i++){
				archivo << desviaciones[i] << " ";
			}
			archivo.close();
		}
	private:
		void normalizar(){ //normalizar vectores para hacerlos comparables es decir z = (X(i) - X(Media)) / Desviacion_Estandar
			int i, j, filas = matriz.rows(), columnas = matriz.cols();
			double desv_est, media;

			for(i = 0 ; i < columnas; i++){
				media = 0;
				desv_est = 0;

				for(j = 0; j < filas; j++){
					media += matriz(j, i);
				}
				media /= filas;

				for(j = 0; j < filas; j++){
					desv_est += (matriz(j, i) - media) * (matriz(j, i) - media);
				}
				desv_est = sqrt(desv_est / filas);

				for(j = 0; j < filas; j++){
					matriz(j, i) = (matriz(j, i) - media) / desv_est;
				}
				promedios.push_back(media);
				desviaciones.push_back(desv_est);
			}
		}

		void correlacion(){
			MatrixXd aux = matriz;
			aux.transposeInPlace();
			correlaciones = aux * matriz;

			double coef = correlaciones(0, 0);
			for(int i = 0 ; i < correlaciones.rows(); i++){
				for(int j = 0 ; j < correlaciones.cols(); j++){
					correlaciones(i, j) = correlaciones(i, j) / coef;
				}
			}
        	}

		MatrixXd eigenvalores(double porcentaje){
			EigenSolver<MatrixXd> sol(correlaciones);
			VectorXcd e, auxV;
			MatrixXcd v;
			MatrixXd res;
			complex<double> auxC;
			double *modulo, menor, auxD;
			int i, j, *posOrig, filas = correlaciones.rows();

			e = sol.eigenvalues();
			v = sol.eigenvectors();

			/***** Obtención del módulo de los eigenvalores para hacerlos comparables *****/
			modulo = new double[filas];
			for(i = 0; i < filas; i++){
				modulo[i] = sqrt((e(i).imag() * e(i).imag()) + (e(i).real() * e(i).real()));
			}

			/**************** Ordenar datos de mayor a menor ****************/
			for(i = 0; i < filas; i++){
				for(j = i + 1; j < filas; j++){
					if(modulo[i] < modulo[j]){
						auxD = modulo[j];
						auxC = e(j);
						auxV = v.col(j);
						modulo[j] = modulo[i];
						modulo[i] = auxD;
						e(j) = e(i);						
						e(i)= auxC;
						v.col(j) = v.col(i);
						v.col(i) = auxV;
					}
				}
			}

			/**************** Obtencion del Porcentaje de Variabilidad ****************/
			auxD = 0;
			int importantes;
			for(i = 0; i < filas; i++){ // eigenvalor / # de eigenvalores * 100
				auxD += (modulo[i] / filas * 100);
				if(auxD >= porcentaje){
					importantes = i + 1;
					break;
				}
			}

			eigenVectores = MatrixXd::Zero(filas, importantes);
			for(i = 0; i < filas; i++){
				for(j = 0; j < importantes; j++){
					eigenVectores(i, j) = sqrt((v(i, j).imag() * v(i, j).imag()) + (v(i, j).real() * v(i, j).real()));
				}
			}
			delete[] modulo;

			return eigenVectores;
		}
};

int main(int argc, char* argv[]){

	if(argc < 3){
		cout << "argumentos insuficientes !" << endl;
	}

	string nomArchivo = argv[1];
	int modo = atoi(argv[2]);
	PCA pca(nomArchivo, modo);

	if(modo == 0){ //reducir
		cout << "archivo cargado. reduciendo matriz.." << endl;
		MatrixXd reducida = pca.reducir(95.0); //95% de informacion deseado
		pca.imprimirArchivo(nomArchivo);

		ofstream archivo("Reducida.txt");
		archivo << reducida.rows() << " " << reducida.cols() << endl;
		archivo << reducida << endl << endl;	
		archivo.close();
	}
	else if(modo == 1){ //solo cargar archivo
		cout << "archivo cargado" << endl;
	}

	
	string continuar;
	cout << "reducir vector? ";
	cin >> continuar;

	while(continuar == "si" || continuar == "s"){
		ofstream archivo("Patron.txt");
		MatrixXd vector = pca.transformarVector();

		archivo << vector.cols() << endl;
		archivo << vector << endl << endl;

		cout << "continuar reduciendo? ";
		cin >> continuar;
		archivo.close();
	}
	cout << endl;
	
	return 0;

}
