
#include <stdio.h>
#include <math.h>

#include <iostream>
#include <vector>
#include <fstream>

using namespace std;

class Complejo{
	private:
		double real;
		double imaginaria;
	public:
		Complejo(){
			real = 0;
			imaginaria = 0;
		}

		Complejo(double real_, double imaginaria_){
			real = real_;
			imaginaria = imaginaria_;
		}
		
		Complejo(double thetaRadianes){ //forma polar
			real = cos(thetaRadianes);
			imaginaria = sin(thetaRadianes);
		}

		void setReal(double real_){
			real = real_;
		}

		void setImaginaria(double imaginaria_){
			imaginaria = imaginaria_;
		}

		double getReal(){
			return real;
		}

		double getImaginaria(){
			return imaginaria;
		}

		double Magnitud(){
			return sqrt(real * real + imaginaria * imaginaria);
		}

		Complejo operator+(Complejo otro){ //suma de complejos
			Complejo resultado(real + otro.getReal(), imaginaria + otro.getImaginaria());
			return resultado;
		}
		
		Complejo operator-(Complejo otro){ //resta de complejos
			Complejo resultado(real - otro.getReal(), imaginaria - otro.getImaginaria());
			return resultado;
		}

		Complejo operator*(Complejo otro){ //multiplicacion de complejos
			double a = real, b = imaginaria, c = otro.getReal(), d = otro.getImaginaria();
			Complejo resultado(a * c - b * d, a * d + c * b);
			
			return resultado;
		}

		Complejo operator/(Complejo otro){ //division de complejos
			double a = real, b = imaginaria, c = otro.getReal(), d = otro.getImaginaria();
			Complejo resultado((a * c + b * d) / (c * c + d * d), (b * c - a * d) / (c * c + d * d));
			
			return resultado;
		}

		bool operator==(Complejo otro){ //comparacion de complejos
			if(real == otro.getReal() && imaginaria == otro.getImaginaria()){
				return true;
			}
			return false;
		}

		void imprimir(){ // imprime el numero en la forma n + mi
			cout << " " << real << " + " << imaginaria << "i " << endl;
		}
};

class WAV{
	private:
		char chunkID[5];
		int chunkSize;
		char format[5];
		char subChunk1ID[5];
		int subChunk1Size;
		int audioFormat;
		int numChannels;
		int sampleRate;
		int byteRate;
		int blockAlign;
		int bitsPerSample;
		char subChunk2ID[5];
		int subChunk2Size;
		vector <Complejo> data;
	public:
		WAV(){
		}
		
		WAV(string archivo){
			FILE *f = fopen(archivo.c_str(), "rb");
			
			chunkSize = 0; //inicializamos enteros
			subChunk1Size = 0;
			audioFormat = 0;
			numChannels = 0;
			sampleRate = 0;
			byteRate = 0;
			blockAlign = 0;
			bitsPerSample = 0;
			subChunk2Size = 0;
			
			fread(chunkID, 1, 4, f); //se leen y almacenan cada uno de los parametros del archivo wav
			fread(&chunkSize, 1, 4, f);
			fread(format, 1, 4, f);
			fread(subChunk1ID, 1, 4, f);
			fread(&subChunk1Size, 1, 4, f);
			fread(&audioFormat, 1, 2, f);
			fread(&numChannels, 1, 2, f);
			fread(&sampleRate, 1, 4, f);
			fread(&byteRate, 1, 4, f);
			fread(&blockAlign, 1, 2, f);
			fread(&bitsPerSample, 1, 2, f);
			fread(subChunk2ID, 1, 4, f);
			fread(&subChunk2Size, 1, 4, f);

			chunkID[4] = format[4] = subChunk1ID[4] = subChunk2ID[4] = '\0'; //por si acaso =P
			string cID = chunkID;

			double nMuestras = subChunk2Size / (numChannels * bitsPerSample / 8);
			for(int i = 0; i < nMuestras && !feof(f); i++){
				int muestra = 0;
				if(fread(&muestra, blockAlign, 1, f) != 1){ //si no se lee un bloque entero de informacion romper el ciclo
					break;
				}

				Complejo dato(muestra, 0);
				data.push_back(dato);
			}
			if(cID == "RIFX"){ //usa ordenamiento de bytes big-endian
				for(int i = 0; i < data.size(); i++){
					int muestra = data[i].getReal();

					voltearBytes((char*)&muestra, sizeof(int));
					data[i].setReal(muestra);
				}
			}
			fclose(f);
		}

		~WAV(){
		}
		
		int getSampleRate(){
			return sampleRate;
		}

		int getBitsPerSample(){
			return bitsPerSample;
		}

		double getFS(){
			double fs = 1;
			int i;
			for(i = 0; i < bitsPerSample; i++){
				fs *= 2;
			}
			if(bitsPerSample > 16){
				return fs / 2;
			}
			return fs;
		}

		vector <Complejo> getData(){ //regresar vector de muestras
			return data;
		}

		void imprimirDatos(){
			for(int i = 0; i < data.size(); i++){
				data[i].imprimir();
			}
		}
	private:
		void voltearBytes(char *array, int tam){
			int i;
			char aux;
			
			for(i = 0; i < (tam / 2); i++){
				aux = array[i];
				array[i] = array[tam - i - 1];
				array[tam - i - 1] = aux;
			}
		}    
};

class Segmentador{
	private:
		WAV wav;
	public:
		Segmentador(){
		}
		
		Segmentador(WAV archivo){
			wav = archivo;
		}

		~Segmentador(){
		}

		vector < vector <Complejo> > getVentanas(){ //regresa los segmentos de voz ventaneados con hamming
			vector <Complejo> datos = wav.getData();
			int i, j, sampleRate = wav.getSampleRate(), tiempoSegmento = (int)(sampleRate * 0.01), tam = datos.size() / tiempoSegmento;			
			int** segmentos = new int*[tam];
			
			double M = 0;
			for(i = 0; i < datos.size(); i++){ //obtenemos el valor de amplitud mayor
				double dato = datos[i].getReal();
				if(dato < 0){ //valor absoluto
					dato *= -1;
				}
				if(dato > M){
					M = dato;
				}
			}

			for(i = 0; i < datos.size(); i ++){ //se acoplan amplitudes para que solo caigan de rangos de |0-1|
				datos[i].setReal(datos[i].getReal() / (M + 1));
			}

			j = 0;
			for(i = 0; i < datos.size(); i += tiempoSegmento){ //se hacen divisiones de 10ms
				segmentos[j] = new int[2];
				
				segmentos[j][0] = i;
				if((i + tiempoSegmento) <= datos.size()){
					segmentos[j][1] = i + tiempoSegmento;
				}else{
					segmentos[j][1] = datos.size();
				}
				j++;
			}

			double *EnEs = new double[tam]; //arreglo para energias escaladas de cada segmento
			for(i = 0; i < tam; i++){
				double max = 0, prom = 0;
				for(j = segmentos[i][0]; j < segmentos[i][1]; j++){
					double dato = datos[j].getReal();

					if(dato < 0){ //valor absoluto
						dato *= -1;
					}
					if(dato > max){
						max = dato;
					}
					prom += dato;
				}
				prom /= segmentos[i][1] - segmentos[i][0];

				double b = 0.01, c = 1, d = 2, R;
				
				R = (b * max) / prom;
				EnEs[i] = (c * R) + (d * prom);
			}

			for(i = 0; i < tam; i++){
				double silencio = 0.039, sinVoz = 0.007843; //limites para clasificar segmentos

				if(EnEs[i] < silencio){
					EnEs[i] = 0;
				}else if(EnEs[i] < sinVoz){
					EnEs[i] = 0.5;
				}else{
					EnEs[i] = 1;
				}
			}

			/*for(i = 0; i < tam; i++){
				for(j = segmentos[i][0]; j < segmentos[i][1]; j++){
					cout << EnEs[i] << endl;
				}
			}*/

			vector < vector <Complejo> > ventanas;
			datos = wav.getData(); //regresamos a la informacion original
			for(i = 0; i < tam; i++){

				if(EnEs[i] == 1){ //segmento de voz
					int limInf, limSup;

					limInf = segmentos[i][0];
					while(EnEs[i] == 1 && i < tam){ //mientras siga siendo voz
						limSup = segmentos[i][1];
						i++;
					}

					vector <Complejo> ventana;
					for(j = limInf; j < limSup; j++){
						ventana.push_back(datos[j]);
					}
						ventanas.push_back(ventana);
				}
			}

			/*for(i = 0; i < tam; i++){
				delete[] segmentos[i];
			}
			delete[] segmentos;
			delete[] EnEs;*/

			return ventanas;
		}
};

class Parametros{
	private:
		double F0; //Frecuencia Fundamental
		double T0;
		double Fmax;
		double Fmin;
		double DEF0;
		double CVF0;
		double Jitter;
		double Jittfactor;
		double Jittradio;
		double PPQ;
		double DDP;
		double RAP;
		double Amp; //Amplitud/Intensidad
		double DEAmp;
		double CVAmp;
		double Shimmer;
		double Shimmfactor;
		double APQ;
	public:
		Parametros(){
		}
		
		~Parametros(){
		}

		void setF0(double n){
			F0 = n;
		}

		void setT0(double n){
			T0 = n;
		}

		void setFmax(double n){
			Fmax = n;
		}

		void setFmin(double n){
			Fmin = n;
		}

		void setDEF0(double n){
			DEF0 = n;
		}

		void setCVF0(double n){
			CVF0 = n;
		}

		void setJitter(double n){
			Jitter = n;
		}

		void setJittfactor(double n){
			Jittfactor = n;
		}

		void setJittradio(double n){
			Jittradio = n;
		}
		
		void setPPQ(double n){
			PPQ = n;
		}

		void setDDP(double n){
			DDP = n;
		}

		void setRAP(double n){
			RAP = n;
		}

		void setAmp(double n){
			Amp = n;
		}

		void setDEAmp(double n){
			DEAmp = n;
		}

		void setCVAmp(double n){
			CVAmp = n;
		}

		void setShimmer(double n){
			Shimmer = n;
		}

		void setShimmfactor(double n){
			Shimmfactor = n;
		}

		void setAPQ(double n){
			APQ = n;
		}

		void imprimirDatos(){
			cout << "-- Frecuencia Fundamental --" << endl;
			cout << "f0 promedio = " << F0 << endl;
			cout << "t0 promedio = " << T0 << endl;
			cout << "fmax = " << Fmax << endl;
			cout << "fmin = " << Fmin << endl;
			cout << "desviacion estandar = " << DEF0 << endl;
			cout << "coeficiente de variacion = " << CVF0 << endl;
			cout << "perturbacion ciclo a ciclo (jitter prom.) = " << Jitter << endl;
			cout << "jitter factor = " << Jittfactor << endl;
			cout << "jitter radio = " << Jittradio << endl;
			cout << "PPQ = " << PPQ << endl;
			cout << "DDP = " << DDP << endl;
			cout << "RAP = " << RAP << endl;
			cout << "-- Amplitud/Intensidad --" << endl;
			cout << "amplitud promedio = " << Amp << endl;
			cout << "desviacion estandar = " << DEAmp << endl;
			cout << "coeficiente de variacion = " << CVAmp << endl;
			cout << "perturbacion ciclo a ciclo (shimmer prom.) = " << Shimmer << endl;
			cout << "shimmer factor = " << Shimmfactor << endl;
			cout << "APQ = " << APQ << endl;
		}

		void imprimirVector(){
			cout << F0 << " ";
			cout << T0 << " ";
			cout << Fmax << " ";
			cout << Fmin << " ";
			cout << DEF0 << " ";
			cout << CVF0 << " ";
			cout << Jitter << " ";
			cout << Jittfactor << " ";
			cout << Jittradio << " ";
			cout << PPQ << " ";
			cout << DDP << " ";
			cout << RAP << " ";
			cout << Amp << " ";
			cout << DEAmp << " ";
			cout << CVAmp << " ";
			cout << Shimmer << " ";
			cout << Shimmfactor << " ";
			cout << APQ << " ";
			cout << endl;
		}
};

class Analizador{
	private:
		int sampleRate;
		double FS;
		Parametros parametros;
		vector <Complejo> muestras;
	public:
		Analizador(){
		}
	
		void setSampleRate(int sampleRate_){
			sampleRate = sampleRate_;
		}

		void setFS(double fs){
			FS = fs;
		}

		void setMuestras(vector <Complejo> muestras_){
			muestras = muestras_;
		}

		Parametros getParametros(){
			return parametros;
		}

		///////////////////////
		void analizar(){
			int i, j, tiempoSegmento = (int)(sampleRate * 0.0312), tam;
			int** segmentos;
			int lim = 1;

			while((lim * 2) <= tiempoSegmento){ //para computar un numero de muestras potencia de 2 para la fft
				lim *= 2;
			}
			tam = muestras.size() / lim;

			if(tam < 11){
				return;
			}

			double f0 = 0, t0 = 0, fmin = 1500, fmax = 0;
			vector <double> frecuencias, periodos, amplitudes; //Hertz, segundos, dbFS (decibelFullScale)

			for(i = 0; i < tam; i++){ //obtener frecuencias fundamentales (f0's)
				if((i * (lim + 1)) <= muestras.size()){
					vector <Complejo> segmento;
					for(j = 0; j < lim; j++){
						segmento.push_back(muestras[(i * lim) + j]);
					}
					hamming(segmento); //se aplica ventaneo hamming

					double f = cepstrum(segmento);
					if(f < fmin){
		            			fmin = f;
					}
					if(f > fmax){
						fmax = f;
					}

					f0 += f;
					frecuencias.push_back(f);
					periodos.push_back(1 / f);
					
					double rms = 0;
					int datosXperiodo = (int)(sampleRate * f);
					for(j = 0; j < datosXperiodo && j < lim; j++){
						double muestra = (muestras[(i * lim) + j]).getReal() / FS;
						rms += muestra * muestra;
					}
					rms /= j;
					rms = sqrt(rms);
					amplitudes.push_back(-20 * log10(rms)); //para que se guarden en unidades positivas
				}
			}
			f0 /= frecuencias.size();

			t0 = 0;
			for(i = 0; i < periodos.size(); i++){
				t0 += periodos[i];
			}
			t0 /= periodos.size();

			parametros.setF0(f0);
			parametros.setT0(t0);
			parametros.setFmax(fmax);
			parametros.setFmin(fmin);
			
			////desviacion estandar
			double desEst = 0;
			for(i = 0; i < frecuencias.size(); i++){
				double aux = frecuencias[i] - f0;
				desEst += aux * aux;
			}
			desEst = sqrt(desEst / frecuencias.size());

			parametros.setDEF0(desEst);
			parametros.setCVF0(100 * desEst / f0);
			//

			////jitter
			double jitt = 0;
			for(i = 0; i < frecuencias.size() - 1; i++){
				double aux = frecuencias[i] - frecuencias[i + 1];
				
				if(aux < 0){
					aux *= -1;
				}
				jitt += aux;
			}
			jitt /= frecuencias.size() - 1;

			parametros.setJitter(jitt);
			parametros.setJittfactor(jitt / f0);
			parametros.setJittradio((1 / jitt) / (t0));
			//

			////period perturbation quotient
			double ppq = 0;
			for(i = 2; i < periodos.size() - 2; i++){
				double aux = ((periodos[i - 2] + periodos[i - 1] + periodos[i] + periodos[i + 1] + periodos[i + 2]) / 5) - periodos[i];
				
				if(aux < 0){
					aux *= -1;
				}
				ppq += aux;
			}
			ppq /= periodos.size() - 4;
			ppq /= t0;

			parametros.setPPQ(ppq);
			//

			////second-order difference of the interval process
			double ddp = 0;
			for(i = 1; i < periodos.size() - 1; i++){
				double aux = (periodos[i + 1] - periodos[i]) - (periodos[i] - periodos[i - 1]);
				
				if(aux < 0){
					aux *= -1;
				}
				ddp += aux;
			}
			ddp /= periodos.size() - 2;
			ddp /= t0;

			parametros.setDDP(ddp);
			//

			////relative average perturbation
			double rap = 0;
			for(i = 1; i < periodos.size() - 1; i++){
				double aux = ((periodos[i - 1] + periodos[i] + periodos[i + 1]) / 3) - periodos[i];
				
				if(aux < 0){
					aux *= -1;
				}
				rap += aux;
			}
			rap /= periodos.size() - 2;
			rap /= t0;

			parametros.setRAP(rap);
			//

			// -- amplitudes 
			double amp = 0;
			for(i = 0; i < amplitudes.size(); i++){
				amp += amplitudes[i];
			}
			amp /= amplitudes.size();
			parametros.setAmp(amp);

			////desviacion estandar
			double desEst2 = 0;
			for(i = 0; i < amplitudes.size(); i++){
				double aux = amplitudes[i] - amp;
				desEst2 += aux * aux;
			}
			desEst2 = sqrt(desEst2 / amplitudes.size());
			
			parametros.setDEAmp(desEst2);
			parametros.setCVAmp(100 * desEst2 / amp);
			//
			
			////shimmer
			double shimm = 0;
			for(i = 0; i < amplitudes.size() - 1; i++){
				double aux = 20 * log10(amplitudes[i] / amplitudes[i + 1]);

				if(aux < 0){
					aux *= -1;
				}
				shimm += aux;
			}
			shimm /= amplitudes.size() - 1;

			parametros.setShimmer(shimm);
			parametros.setShimmfactor(100 * shimm / amp);
			//

			////amplitude perturbation quotient
			double apq = 0;
			for(i = 5; i < amplitudes.size() - 5; i++){
				double aux = 0;
				for(int j = -5; j <= 5; j++){
					aux += amplitudes[i + j] / 11;
				}
				aux -= amplitudes[i];

				if(aux < 0){
					aux *= -1;
				}
				apq += aux;
			}
			apq /= amplitudes.size() - 10;
			apq /= amp;

			parametros.setAPQ(apq);
			//

			parametros.imprimirVector();
		}
		//////////////
	private:
		void hamming(vector <Complejo> &entrada){
			for(int i = 0; i < entrada.size(); i++){
				Complejo hamm(0.54  - (0.46 * cos((2 * 3.14159265 * i) / (entrada.size() - 1))), 0);
				entrada[i] = entrada[i] * hamm;
			}
		}

		double cepstrum(vector <Complejo> muestras){
			int i, muestra, inicio, fin;
			double pico = 0;
			vector <Complejo> fft, ifft;

			fft = FFT(muestras);
			for(i = 0; i < fft.size(); i++){
				fft[i].setReal(log10(fft[i].Magnitud()));
				fft[i].setImaginaria(0);
			}
			ifft = IFFT(fft);
			inicio = (int)(0.0009091 * sampleRate); //rango de frecuencias en las que se encuentra la voz humana
			fin = (int)(0.0125 * sampleRate);

			for(i = inicio; i < fin; i++){
				if(ifft[i].getReal() > pico){
					pico = ifft[i].getReal();
					muestra = i;
				}
			}
			
			return 1 / ((double)muestra / sampleRate);
		}

		vector <Complejo> FFT(vector <Complejo> entrada){
			int i, lim = entrada.size();
			Complejo *X;
			vector <Complejo> fft;
			
			X = new Complejo[lim];
			for(int i = 0; i < lim; i++){
				X[i] = entrada[i];
			}

			X = FFT(X, lim, 1);
			for(int i = 0; i < lim; i++){
				fft.push_back(X[i]);
			}
			delete X;
			
			return fft;
		}

		vector <Complejo> IFFT(vector <Complejo> entrada){
			int i, lim = entrada.size();
			Complejo *X;
			vector <Complejo> ifft;

			X = new Complejo[lim];
			for(int i = 0; i < lim; i++){
				X[i] = entrada[i];
			}

			X = FFT(X, lim, -1);
			Complejo ajPer(1.0 / lim, 0); //para ajustar el periodo al regresar al dominio del tiempo
			for(int i = 0; i < lim; i++){
				ifft.push_back(X[i] * ajPer);
			}
			delete X;

			return ifft;
		}

		Complejo* FFT(Complejo *x, int N, int direccion){ //N debe ser potencia de 2
			Complejo *X = new Complejo[N];
			Complejo *d, *e, *D, *E;
			int k;

			if (N == 1){
				X[0] = x[0];
				return X;
			}

			e = new Complejo[N / 2];
			d = new Complejo[N / 2];
			for(k = 0; k < N / 2; k++){
				e[k] = x[2 * k];
				d[k] = x[2 * k + 1];
			}
			E = FFT(e, N / 2, direccion);
			D = FFT(d, N / 2, direccion);
			delete e;
			delete d;

			for(k = 0; k < N / 2; k++){ //Multipliclar D por los factores e^(-2 * pi * i / N * k)
				Complejo aux(direccion * -2.0 * 3.14159265 * k / N);
				D[k] = aux * D[k];
			}
			for(k = 0; k < N / 2; k++){
				X[k] = E[k] + D[k];
				X[k + N / 2] = E[k] - D[k];
			}

			delete D;
			delete E;

			return X;
		}
};

int main(int argc, char* argv[]){
	if(argc != 2){
		return -1;
	}

	WAV archivo(argv[1]);
	Segmentador segmentador(archivo);	
	Analizador analizador;
	
	vector < vector <Complejo> > ventanas = segmentador.getVentanas();
	analizador.setSampleRate(archivo.getSampleRate());
	analizador.setFS(archivo.getFS());
	for(int i = 0; i < ventanas.size(); i++){
		analizador.setMuestras(ventanas[i]);
		analizador.analizar();
		//analizador.getParametros().imprimirVector();
	}

	return 0;
}
