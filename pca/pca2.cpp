#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <fstream>

#include "eigen/Eigen/Core"
#include "eigen/Eigen/Eigenvalues"

using namespace std;
using namespace Eigen;

class Matriz
{
    private:
        MatrixXd mat;
        int filas;
        int columnas;

    public:

        Matriz(int f, int c) // Constructor
        {
            mat = MatrixXd::Random(f, c);

            filas = f;
            columnas = c;
            
            for(int i = 0 ; i < f ; i++)
                for(int j = 0 ; j < c ; j++)
                    mat(i, j) = int(mat(i, j)*100) % 11;
        }
        
        void impMat()
        {
           /* ofstream salida;
            salida.open("salida.txt", fstream::app);
            salida << mat << endl << endl;
            salida.close(); */
            cout << mat << endl << endl;
        }
        
        MatrixXd regresaMatriz()
        {
            return mat;
        }
        
        void llenaMatriz()
        {
            double var;
            
            ifstream entrada;
            entrada.open("Valores3.txt");
            
            for(int i = 0 ; i < filas ; i++)
                for(int j = 0 ; j < columnas ; j++)
                {
                    entrada >> var;
                    mat(i, j) = var;
                }
            entrada.close();
        }
        
        float position(int f, int c)
        {
            return mat(f, c);
        }

        void setPosition(int f, int c, float val)
        {
            mat(f, c) = val;
        }

        void normalizar() // NORMALIZAR VECTORES PARA HACERLOS COMPARABLES
        {                 // es decir z = (X(i) - X(Media)) / Desviacion_Estandar
            int i = 0, j = 0;
            double desv_est, media;

            for(i = 0 ; i < columnas ; i++)
            {
                j = 0;
                media = 0;
                desv_est = 0;

                for(j = 0 ; j < filas ; j++)
                    media += mat(j, i);

                media /= filas;

                for(j = 0 ; j < filas ; j++)
                    desv_est += (mat(j, i) - media) * (mat(j, i) - media);


                desv_est = sqrt(desv_est/(filas-1));

                for(j = 0 ; j < filas ; j++)
                    mat(j, i) = (mat(j, i) - media) / desv_est;
            }
        }

        void transpuestaDe(Matriz *m, int f, int c)
        {
            for(int i = 0 ; i < f ; i++)
                for(int j = 0 ; j < c ; j++)
                    mat(j, i) = m->position(i, j);
        }

        void correlacion(Matriz *m1, Matriz *m2, int num_vars, int num_m)
        {
            double res, v1, v2, coef;
            int i, j, k;

            for(i = 0 ; i < num_vars ; i++) // Multiplicacion de la transpuesta por la normal
            {
                for(j = 0 ; j < num_vars ; j++)
                {
                    res = 0;
                    for(k = 0 ; k < num_m ; k++)
                    {
                        v1 = m1->position(i, k);
                        v2 = m2->position(k, j);
                        res += v1 * v2;
                    }
                    mat(i, j) = res;
                }
            }
            
            // Dividir entre un coeficiente para hacer unos la diagonal
            
            coef = mat(0, 0);

            for(i = 0 ; i < num_vars ; i++) 
                for(j = 0 ; j < num_vars ; j++)
                    mat(i, j) = mat(i, j) / coef;
        }
        
        MatrixXd eigenvalores(double porcentaje)
        {
            EigenSolver<MatrixXd> sol(mat);
            VectorXcd e, auxV;
            MatrixXcd v;
            MatrixXd res;
            complex<double> auxC;
            
            double *modulo, menor, auxD;
            int i, j, *posOrig, auxI;
            
            e = sol.eigenvalues();
            v = sol.eigenvectors();
            
            modulo = (double*)malloc(filas*sizeof(double));
            posOrig = (int*)malloc(filas*sizeof(int));
            
            for(i = 0 ; i < filas ; i++)
                posOrig[i] = i+1;
            
            /***** Obtención del módulo de los eigenvalores para hacerlos comparables *****/
            
            for(i = 0 ; i < filas ; i++)
                modulo[i] = sqrt((e(i).imag() * e(i).imag()) + (e(i).real() * e(i).real()));
            
            /**************** Ordenar datos de mayor a menor ****************/
            
            for(i = 0 ; i < filas ; i++)
                for(j = i + 1 ; j < filas ; j++)
            		if(modulo[i] < modulo[j])
            		{
            			auxD = modulo[j];
            			auxC = e(j);
            			auxV = v.col(j);
            			auxI = posOrig[j];
            			
            			modulo[j] = modulo[i];
            			e(j) = e(i);
            			v.col(j) = v.col(i);
            			posOrig[j] = posOrig[i];
            			
            			modulo[i] = auxD;
            			e(i)= auxC;
            			v.col(i) = auxV;
            			posOrig[i] = auxI;
            		}
                
            cout << "Eigenvalores Ordenados: " << endl << endl << e << endl << endl;
            cout << "Eigenvectores Ordenados: " << endl << endl << v << endl << endl;
            
            /**************** Obtencion del Porcentaje de Variabilidad Explicado ****************/
            // eigenvalor / # de eigenvalores * 100
            
            cout << "Porcentajes de Variabilidad Explicados" << endl;
            
            auxD = 0;
            int importantes;
            bool bandera = false;
            
            for(i = 0 ; i < filas ; i++)
            {
                auxD += (modulo[i] / filas * 100);
                cout << endl << e(i) << " = " << (modulo[i] / filas * 100) << "%";
                if(auxD >= porcentaje && !bandera)
                {
                    importantes = i+1;
                    bandera = true;
                }
            }
            
            cout << endl << endl << "Variables que representan el " << porcentaje << "%: " << importantes;
            cout << endl << endl << "Posicion Original de las Variables ya ordenadas: ";
            
            for(i = 0 ; i < filas ; i ++)
                cout << endl << i+1 << ": " << posOrig[i];
                
            res = MatrixXd::Zero(filas, importantes);
            
            for(i = 0 ; i < filas ; i++)
                  for(j = 0 ; j < importantes ; j++)
                        res(i, j) = sqrt((v(i, j).imag() * v(i, j).imag()) + (v(i, j).real() * v(i, j).real()));
            
            free(modulo);
            free(posOrig);
            
            return res;
        }
};

int main()
{
    /*****************************************************************************************************/
    int nMuestras = 5, nVariables = 5; // En orden: numero de muestras y numero de variables.
    double pctObjetivo = 95.00; // porcentaje de informacion que se desea conservar despues de la reducción
    /*****************************************************************************************************/
    
    MatrixXd m4, m5, matrizReducida;
    
    /*cin >> nMuestras;
    cin >> nVariables;*/
    
    Matriz *m1 = new Matriz(nMuestras, nVariables); // Matriz Original Normalizada
    Matriz *m2 = new Matriz(nVariables, nMuestras); // Matriz Transpuesta Normalizada
    Matriz *m3 = new Matriz(nVariables, nVariables); // Matriz de Correlaciones
    
    cout << "Matriz Original: " << endl << endl;
    m1->llenaMatriz();
    m1->impMat();
    
    cout << "Matriz Normalizada: " << endl << endl;
    m1->normalizar();
    m1->impMat();

    cout << "Matriz Transpuesta: " << endl << endl;
    m2->transpuestaDe(m1, nMuestras, nVariables);
    m2->impMat();
    
    cout << "Matriz de Correlaciones: " << endl << endl;
    m3->correlacion(m2, m1, nVariables, nMuestras);
    m3->impMat();
    
    m4 = m3->eigenvalores(pctObjetivo); // Obtener la matriz de eigenvectores ordenados y reducidos
    
    // La operación para reducir las dimensiones: CP = EigenvectoresReducidos' x MatrizNormalizada'
    
    m4.transposeInPlace();
    m5 = m2->regresaMatriz();
    
    cout << endl << endl << "Multiplicar (EigenvectoresReducidos')" << endl << endl;
    cout << m4 << endl << endl << "Por (MatrizNormalizada')" << endl << m5;
    
    cout << endl << endl << m4 * m5;
    matrizReducida = m4 * m5;
    matrizReducida.transposeInPlace();
    
    cout << endl << endl << "Resultado Final: " << endl << endl << matrizReducida;
    
    ofstream salida;
    salida.open("salida.txt", fstream::app);
    salida << matrizReducida << endl << endl;
    
    salida.close();
    
    free(m1);
    free(m2);
    free(m3);
    
    return 0;
}
