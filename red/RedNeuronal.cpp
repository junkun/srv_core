
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class RedNeuronal{
	private:
		int *ranpat; //arreglo para movimientos aleatorios de los patrones
		int numPatrones, numEntradas, numOcultas, numSalidas; //numeros de patrones, neuronas de entrada, neuronas capa oculta y neuronas de salida
		double **Objetivo, **Entrada, **Oculta, **Salida; //matrices de objetivos, entrada, capa oculta y salida
		double **PesosEO, **PesosOS; //matrices de pesos de entradaAoculta y ocultaAsalida
		double *SumRet, **SumO, **SumS; //sumatorias de retropropagacion, capa oculta y salida
		double **DeltaPesosEO, *DeltaO, **DeltaPesosOS, *DeltaS; //deltas de entradaAoculta, oculta, ocultaAsalida y salida
		double Error, velAp, alpha, valMax, errormin; //error global, razon de aprendizaje, alfa, valor maximo de pesos iniciales, error minimo de parada
		long maxEpoch; //maximo de epocas
		bool Entrenada; //indica si se ha entrenado satisfactoriamente la red
		vector <double> promedios, desviaciones;
	public:
		RedNeuronal(){
		}
		
		RedNeuronal(int nP, int nE, int nO, int nS){ //numero de patrones, numero neuronas de entrada, numero neuronas capa oculta, numero neuronas salida
			numPatrones = nP; //asignamos numeros de neuronas
			numEntradas = nE;
			numOcultas = nO;
			numSalidas = nS;

			ranpat = new int[numPatrones + 1]; //se asigna memoria para arreglos y matrices
			Entrada = crearMat(numPatrones + 1, numEntradas + 1);
			Objetivo = crearMat(numPatrones + 1, numSalidas + 1);
			SumO = crearMat(numPatrones + 1, numOcultas + 1);
			PesosEO = crearMat(numEntradas + 1, numOcultas + 1);
			Oculta = crearMat(numPatrones + 1, numOcultas + 1);
			SumS = crearMat(numPatrones + 1, numSalidas + 1);
			PesosOS = crearMat(numOcultas + 1, numSalidas + 1);
			Salida = crearMat(numPatrones + 1, numSalidas + 1);
			DeltaS = new double[numSalidas + 1];
			SumRet = new double[numOcultas + 1];
			DeltaO = new double[numOcultas + 1];
			DeltaPesosEO = crearMat(numEntradas + 1, numOcultas + 1);
			DeltaPesosOS = crearMat(numOcultas + 1, numSalidas + 1);

			velAp = 0.5; //valores por defecto
			alpha = 0.9;
			valMax = 0.5;
			errormin = 0.0004;
			maxEpoch = 10000;
			Entrenada = false;
		}

		RedNeuronal(string nomArchivo, int modo){
			int i, j;
			ifstream archivo(nomArchivo.c_str());

			archivo >> numPatrones;
			archivo >> numEntradas;
			archivo >> numOcultas;
			archivo >> numSalidas;
			archivo >> velAp;
			archivo >> errormin;
			archivo >> maxEpoch;

			ranpat = new int[numPatrones + 1]; //se asigna memoria para arreglos y matrices
			Entrada = crearMat(numPatrones + 1, numEntradas + 1);
			Objetivo = crearMat(numPatrones + 1, numSalidas + 1);
			SumO = crearMat(numPatrones + 1, numOcultas + 1);
			PesosEO = crearMat(numEntradas + 1, numOcultas + 1);
			Oculta = crearMat(numPatrones + 1, numOcultas + 1);
			SumS = crearMat(numPatrones + 1, numSalidas + 1);
			PesosOS = crearMat(numOcultas + 1, numSalidas + 1);
			Salida = crearMat(numPatrones + 1, numSalidas + 1);
			DeltaS = new double[numSalidas + 1];
			SumRet = new double[numOcultas + 1];
			DeltaO = new double[numOcultas + 1];
			DeltaPesosEO = crearMat(numEntradas + 1, numOcultas + 1);
			DeltaPesosOS = crearMat(numOcultas + 1, numSalidas + 1);

			for(i = 1; i <= numPatrones; i++){
				for(j = 1; j <= numEntradas; j++){
					archivo >> Entrada[i][j];
				}
			}

			for(i = 1; i <= numPatrones; i++){
				for(j = 1; j <= numSalidas; j++){
					archivo >> Objetivo[i][j];
				}
			}

			Entrenada = false;
			if(modo == 1 || modo == 2){ //cargar ademas matrices de pesos y umbrales
				for(i = 0; i <= numEntradas; i++){
					for(j = 1; j <= numOcultas; j++){
						archivo >> PesosEO[i][j];
					}
				}


				for(i = 0; i <= numOcultas; i++){
					for(j = 1; j <= numSalidas; j++){
						archivo >> PesosOS[i][j];
					}
				}

				for(i = 0; i < numEntradas; i++){
					double aux;					
					archivo >> aux;
					promedios.push_back(aux);
				}

				for(i = 0; i < numEntradas; i++){
					double aux;					
					archivo >> aux;
					desviaciones.push_back(aux);
				}				

				Entrenada = true;
			}
			
			alpha = 0.9;
			valMax = 1.5; //0.5

			archivo.close();
		}

		~RedNeuronal(){
			delete ranpat;
			destruirMat(Entrada, numPatrones + 1);
			destruirMat(Objetivo, numPatrones + 1);
			destruirMat(SumO, numPatrones + 1);
			destruirMat(PesosEO, numEntradas + 1);
			destruirMat(Oculta, numPatrones + 1);
			destruirMat(SumS, numPatrones + 1);
			destruirMat(PesosOS, numOcultas + 1);
			destruirMat(Salida, numPatrones + 1);
			delete DeltaS;
			delete SumRet;
			delete DeltaO;
			destruirMat(DeltaPesosEO, numEntradas + 1);
			destruirMat(DeltaPesosOS, numOcultas + 1);
		}

		void setEntradas(double *e){
			int i, j, k = 0;

			for(i = 0; i <= numPatrones; i++){
				Entrada[i][0] = 0;
			}

			for(i = 0; i <= numEntradas; i++){
				Entrada[0][i] = 0;
			}

			for(i = 1; i <= numPatrones; i++){
				for(j = 1; j <= numEntradas; j++){
					Entrada[i][j] = e[k++];
				}
			}
		}

		void setObjetivos(double *o){
			int i, j, k = 0;

			for(i = 0; i <= numPatrones; i++){
				Entrada[i][0] = 0;
			}

			for(i = 0; i <= numSalidas; i++){
				Entrada[0][i] = 0;
			}

			for(i = 1; i <= numPatrones; i++){
				for(j = 1; j <= numSalidas; j++){
					Objetivo[i][j] = o[k++];
				}
			}
		}

		void setRazonAprendizaje(double ra){
			velAp = ra;
		}

		void setEpocasMaximas(long em){
			maxEpoch = em;
		}

		void setErrorMinimo(double em){
			errormin = em;
		}

		bool entrenada(){
			return Entrenada;
		}

		void inicializarPesos(){
			int i, j;

			for(i = 0; i <= numEntradas; i++){ //inicializar pesos y deltas entradaAoculta
				for(j = 1; j <= numOcultas; j++){
					DeltaPesosEO[i][j] = 0.0;
					PesosEO[i][j] = 2.0 * (randomize() - 0.5) * valMax;
				}
			}

			for(i = 0; i <= numOcultas; i++){ //inicializar pesos y deltas ocultaAsalida
				for(j = 1; j <= numSalidas; j++){
					DeltaPesosOS[i][j] = 0.0 ;
					PesosOS[i][j] = 2.0 * (randomize() - 0.5) * valMax;
				}
			}
		}

		void normalizarPatrones(){
			int i, j;
			double desv_est, media;

			for(i = 1 ; i <= numEntradas; i++){
				media = 0;
				desv_est = 0;
				for(j = 1; j <= numPatrones; j++){
					media += Entrada[j][i];
				}
				media /= numPatrones;
				for(j = 1; j <= numPatrones; j++){
					desv_est += (Entrada[j][i] - media) * (Entrada[j][i] - media);
				}
				desv_est = sqrt(desv_est / numPatrones);
				for(j = 1; j <= numPatrones; j++){
					Entrada[j][i] = (Entrada[j][i] - media) / desv_est;
				}
				promedios.push_back(media);
				desviaciones.push_back(desv_est);
			}
		}

		void entrenar(){
			int i, j, k, p, np, epoch;

			for(epoch = 0; epoch < maxEpoch ; epoch++){
				for(p = 1; p <= numPatrones; p++){ //randomizacion de orden de entrada
					ranpat[p] = p;
				}
				for(p = 1; p <= numPatrones; p++){
					int aux;

					np = p + randomize() * (numPatrones + 1 - p);
					aux = ranpat[p];
					ranpat[p] = ranpat[np];
					ranpat[np] = aux;
				}

				Error = 0.0;
				for(np = 1; np <= numPatrones; np++){
					p = ranpat[np];
					for(j = 1; j <= numOcultas; j++){ //propagacion en capa oculta
						SumO[p][j] = PesosEO[0][j];
						for(i = 1; i <= numEntradas; i++){
							SumO[p][j] += Entrada[p][i] * PesosEO[i][j];
						}
						Oculta[p][j] = 1.0 / (1.0 + exp(-SumO[p][j]));
					}
					
					for(k = 1; k <= numSalidas; k++){ //propagacion en capa salida y calculo de errores
						SumS[p][k] = PesosOS[0][k];
						for(j = 1; j <= numOcultas; j++){
							SumS[p][k] += Oculta[p][j] * PesosOS[j][k];
						}

						//Salida[p][k] = SumS[p][k]; // salida lineal
						Salida[p][k] = 1.0 / (1.0 + exp(-SumS[p][k])); // salida sigmoidal
						
						//Error += 0.5 * (Objetivo[p][k] - Salida[p][k]) * (Objetivo[p][k] - Salida[p][k]); //error medio cuadratico
						Error -= (Objetivo[p][k] * log(Salida[p][k]) + (1.0 - Objetivo[p][k]) * log(1.0 - Salida[p][k])); //error de entropia cruzada

						//DeltaS[k] = (Objetivo[p][k] - Salida[p][k]) * Salida[p][k] * (1.0 - Salida[p][k]); //salidas sigmoidales, requiere ems
						//DeltaS[k] = Objetivo[p][k] - Salida[p][k]; //salidas lineales, requiere ems
						DeltaS[k] = Objetivo[p][k] - Salida[p][k]; //salidas sigmoidales, requiere eec
						
					}
					
					for(j = 1; j <= numOcultas; j++){ //retropropagacion de errores a capa oculta
						SumRet[j] = 0.0;
						for(k = 1; k <= numSalidas; k++){
							SumRet[j] += PesosOS[j][k] * DeltaS[k];
						}
						DeltaO[j] = SumRet[j] * Oculta[p][j] * (1.0 - Oculta[p][j]);
					}
					
					for(j = 1; j <= numOcultas; j++){ //actualizacion de pesos entradaAoculta
						DeltaPesosEO[0][j] = velAp * DeltaO[j] + alpha * DeltaPesosEO[0][j];
						PesosEO[0][j] += DeltaPesosEO[0][j];
						for(i = 1; i <= numEntradas; i++){
							DeltaPesosEO[i][j] = velAp * Entrada[p][i] * DeltaO[j] + alpha * DeltaPesosEO[i][j];
							PesosEO[i][j] += DeltaPesosEO[i][j];
						}
					}

					for(k = 1 ; k <= numSalidas; k++){ //actualizacion de pesos ocultaAsalida
						DeltaPesosOS[0][k] = velAp * DeltaS[k] + alpha * DeltaPesosOS[0][k];
						PesosOS[0][k] += DeltaPesosOS[0][k];
						for(j = 1; j <= numOcultas; j++){
							DeltaPesosOS[j][k] = velAp * Oculta[p][j] * DeltaS[k] + alpha * DeltaPesosOS[j][k];
							PesosOS[j][k] += DeltaPesosOS[j][k];
						}
					}
				}

				if(epoch % 2000 == 0){
					cout << "Epoca no." << epoch << " Error = " << Error << endl;
				}
				
				if(Error < errormin){ //detenerse cuando el error sea muy pequeño
					Entrenada = true;
					break;
				}
			}
		}

		void imprimirDatos(string nomArchivo){
			int i, j;
			ofstream archivo(nomArchivo.c_str());

			archivo << numPatrones << " " << numEntradas << " " << numOcultas << " " << numSalidas << " " << velAp << " " << errormin << " " << maxEpoch << endl << endl;
			for(i = 1; i <= numPatrones; i++){
				for(j = 1; j <= numEntradas; j++){
					archivo << Entrada[i][j] << " ";
				}
				archivo << endl;
			}
			archivo << endl;

			for(i = 1; i <= numPatrones; i++){
				for(j = 1; j <= numSalidas; j++){
					archivo << Objetivo[i][j] << " ";
				}
				archivo << endl;
			}
			archivo << endl;

			for(i = 0; i <= numEntradas; i++){
				for(j = 1; j <= numOcultas; j++){
					archivo << PesosEO[i][j] << " ";
				}
				archivo << endl;
			}
			archivo << endl;

			for(i = 0; i <= numOcultas; i++){
				for(j = 1; j <= numSalidas; j++){
					archivo << PesosOS[i][j] << " ";
				}
				archivo << endl;
			}
			archivo << endl;

			for(i = 0; i < promedios.size(); i++){
				archivo << promedios[i] << " ";
			}
			archivo << endl;
			
			for(i = 0; i < desviaciones.size(); i++){
				archivo << desviaciones[i] << " ";
			}
			archivo << endl;
			archivo.close();
		}

		void clasificar(){
			int i, tam;
			double *patron;
			ifstream archivo("Patron.txt");

			archivo >> tam;
			patron = new double[tam];
			for(i = 0; i < tam; i++){
				archivo >> patron[i];
			}
			clasificar(patron);
			delete patron;
		}

		void clasificar(double *patron){
			int i, j, k;
			double sum, *oculta, *salida;

			for(i = 0; i < promedios.size(); i++){
				patron[i] = (patron[i] - promedios[i]) / desviaciones[i];
			}

			oculta = new double[numOcultas + 1];
			salida = new double[numSalidas + 1];

			for(j = 1; j <= numOcultas; j++){ //propagacion en capa oculta
				sum = PesosEO[0][j]; //se comienza con el umbral
				for(i = 1; i <= numEntradas; i++){
					sum +=  patron[i - 1] * PesosEO[i][j];
				}
				oculta[j] = 1.0 / (1.0 + exp(-sum));
			}
					
			for(k = 1; k <= numSalidas; k++){ //propagacion en capa salida
				sum = PesosOS[0][k]; //se comienza con el umbral
				for(j = 1; j <= numOcultas; j++){
					sum += oculta[j] * PesosOS[j][k];
				}

				salida[k] = 1.0 / (1.0 + exp(-sum)); // salida sigmoidal
				//(salida[k] > .5) ? (salida[k] = 1) : (salida[k] = 0);
			}

			for(i = 1; i <= numSalidas; i++){
				//cout << salida[i] << " "; 
				printf("%.2lf%% ", salida[i] * 100);
			}
			cout << endl;

			delete[] oculta;
			delete[] salida;
		}

		void realizarPruebas(){
			int f, c, nMat;
			double *patron;
			ifstream archivo("PatronesPrueba.txt");

			archivo >> nMat;
			for(int k = 0; k < nMat; k++){
				archivo >> f;
				archivo >> c;
				patron = new double[c];
				for(int i = 0; i < f; i++){
					for(int j = 0; j < c; j++){
						archivo >> patron[j];
					}
					clasificar(patron);
				}
				cout << endl;
				delete[] patron;
			}
		}

	private:
		double** crearMat(int filas, int columnas){
			int i, j;
			double **matriz;

			matriz = new double*[filas];
			for(int i = 0; i < filas; i++){
				matriz[i] = new double[columnas];
			}
			
			for(i = 0; i < filas; i++){
				for(j = 0; j < columnas; j++){
					matriz[i][j] = 0;
				}
			}

			return matriz;
 		}
 		
 		void destruirMat(double **matriz, int filas){
			if(matriz != NULL){
				for(int i = 0; i < filas; i++){
					delete matriz[i];
				}
				delete[] matriz;
				matriz = NULL;
			}
		}

		double randomize(){ //regresa un numero aleatorio de [0, 1]
			srand(time(NULL) * rand());
			return (double)rand() / ((double)RAND_MAX + 1);
		}
};

int main(int argc, char* argv[]){
	
	if(argc < 3){
		cout << "argumentos insuficientes !" << endl;
	}

	string nomArchivo = argv[1];
	int modo = atoi(argv[2]);
	RedNeuronal Red(nomArchivo, modo);

	if(modo == 0){ //entrenar
		cout << "archivo cargado. entrenando red..." << endl;

		Red.inicializarPesos();
		Red.normalizarPatrones();
		Red.entrenar();

		if(Red.entrenada()){
			cout << "red entrenada con exito" << endl;
			Red.imprimirDatos(nomArchivo);
		}else{
			cout << "no se pudo entrenar la red" << endl;
			return -1;
		}
	}
	else if(modo == 1){ //solo clasificar
		cout << "archivo cargado" << endl;
	}
	else if(modo == 2){ //realizar pruebas
		Red.realizarPruebas();
		return 0;
	}

	string continuar;
	cout << "clasificar patron? ";
	cin >> continuar;

	while(continuar == "si" || continuar == "s"){
		Red.clasificar();
		cout << "continuar clasificando? ";
		cin >> continuar;
	}

	cout << endl;
	return 0;
}
