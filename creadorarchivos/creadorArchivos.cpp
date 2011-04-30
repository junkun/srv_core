
#include <stdio.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class creadorArchivos{
	private:
		vector <string> archivos;
	public:
		creadorArchivos(string nomArchivo){
			int i;
			ifstream archivo(nomArchivo.c_str());

			while(!archivo.eof()){
				string aux;
				archivo >> aux;
				archivos.push_back(aux);
				cout << aux;
			}
			archivo.close();

			vector < vector <string> > datos;

			for(i = 0; i < archivos.size(); i++){
				char c;
				vector <string> aux;
				FILE *f = fopen(archivos[i].c_str(), "r");

				while(!feof(f)){
					char buffer[1000];
					
					fgets(buffer, 500, f);
					string aux2;
					aux2 = buffer;
					aux.push_back(aux2);
					buffer[0] = 0;
				}
				datos.push_back(aux);
				//fclose(f);
			}

			for(i = 0; i < datos.size(); i++){
				vector <string> aux = datos[i];
				for(int j = 0; j < aux.size(); j++){
					cout << aux[i];
				}
				cout << endl;
			}

			for(i = 0; i < datos.size(); i++){
				vector <string> aux = datos[i];
				for(int j = 0; j < aux.size(); j++){
					for(int k = 0; k < datos.size(); k++){
						if(k == i){
							cout << "1 ";
						}else{
							cout << "0 ";
						}
					}
					cout << endl;
				}
				cout << endl;
			}
		}
};

int main(int argc, char* argv[]){
	string archivo = argv[1];
	creadorArchivos cA(archivo);

	//cA.hacerArchivo();

	return 0;
}
