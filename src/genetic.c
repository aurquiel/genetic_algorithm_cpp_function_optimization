#include "genetic.h"
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <iostream>
#include <string>
#include "exprtk.hpp" //Header que maneja la entrada string de la forumla

genetic::genetic(std::string f_fitness_variable_string,bool f_maximizar_o_minimizar, double long f_limite_superior,
                 double long f_limite_inferior, long long int f_presicion, double long f_probabilidad_de_cruce,
                 double long f_probabilidad_de_mutacion,unsigned long long int f_tamanio_de_poblacion,
                 unsigned long long int f_numero_de_generaciones)
{
    //Pasando los valores externos a variables privadas de la clase
    fitness_variable_string=f_fitness_variable_string;
    maximizar_o_minimizar=f_maximizar_o_minimizar;
    limite_superior=f_limite_superior;
    limite_inferior=f_limite_inferior;
    presicion=f_presicion;
    probabilidad_de_cruce=f_probabilidad_de_cruce;
    probabilidad_de_mutacion=f_probabilidad_de_mutacion;
    tamanio_de_poblacion=f_tamanio_de_poblacion;
    numero_de_generaciones=f_numero_de_generaciones;
    //Final

    //Incialización de variable internas
    i=0; j=0;
    media=0;
    bits=1;
    while((pow(2,(bits-1)))<((limite_superior-limite_inferior)*pow(10,presicion))){
        bits++;
    }
    bits=bits-1;
    //Final

    generacion_de_poblacion_inicial(matriz); //se creaa toda la poblacioin incial, fitness, numero de copias, suma de numero de copias.
    ordenamiento_burbuja(matriz,elite);
    suma_copias(matriz);
    ruleta(matriz,matriz2);
    codificacion(matriz2);
    cruce(matriz2);
}

void genetic::generacion_de_poblacion_inicial(double long p[][4]) //Función que genera la poblacion incial
{
    srand48(time(NULL));  //Valores aleatorios con semilla de tiempo, rand48() solo funciona en linux

    for (i=0;i<tamanio_de_poblacion;i++) //desde i=0 hasta que i sea menor que el numero de inviduos, esto de debe a que el arreglo de los individuos de la poblacion comienza en cero
    {
        p[i][0]=limite_inferior+drand48()*(limite_superior-limite_inferior); //se le agina un numero al azar entre el numero infierior y el superior
        p[i][1]=fitness<double>(p[i][0]);  //Se genera el fitness de los individuos
        media=media+p[i][1];  //se suman todos los fitness
    }

    media=media/tamanio_de_poblacion; //La media de los fitness


    for (i=0;i<tamanio_de_poblacion;i++) //Ciclo para definir el numero de copias
    {
        p[i][2]=numero_copias(p[i][1]);
    }

}

template <typename T> double long genetic::fitness(double long &valor)
{
   typedef exprtk::symbol_table<T> symbol_table_t;
   typedef exprtk::expression<T>     expression_t;
   typedef exprtk::parser<T>             parser_t;

   T  x = T(valor);

   symbol_table_t symbol_table;
   symbol_table.add_variable("x",x);

   expression_t expression;
   expression.register_symbol_table(symbol_table);

   parser_t parser;
   parser.compile(fitness_variable_string,expression);

   if (maximizar_o_minimizar==true) //si es verdadero maximizo
   {
        return expression.value();
   }
   else //minimiza
   {
        return 1/(expression.value());
   }
}

double long genetic::numero_copias(double long &valor) //Funcion que devuelvel el numero de copias
{
    return valor/media; //Devuelve el numero de copias, fitness entre media de fitness
}


inline void genetic::ordenamiento_burbuja(double long p[][4], double long *ptr_elite)
{
    double long aux=0,aux1=0,aux2=0,aux3=0; //variables auxiliares

    for (j=0;j<=tamanio_de_poblacion;j++)
    {
        for (i=0;i<tamanio_de_poblacion-1;i++)
        {
            if(p[i][1]>p[i+1][1])
            {
                //Se guardan los que van a ser sobre escritos
                aux=p[i+1][0];
                aux1=p[i+1][1];
                aux2=p[i+1][2];
                aux3=p[i+1][3];
                //se intercambian
                p[i+1][0]=p[i][0];
                p[i+1][1]=p[i][1];
                p[i+1][2]=p[i][2];
                p[i+1][3]=p[i][3];
                //se recuperan los auxiliares
                p[i][0]=aux;
                p[i][1]=aux1;
                p[i][2]=aux2;
                p[i][3]=aux3;
            }
        }
    }

    if (ptr_elite[1]<p[tamanio_de_poblacion-1][1])
    {
        ptr_elite[0]=p[tamanio_de_poblacion-1][0];
        ptr_elite[1]=p[tamanio_de_poblacion-1][1];
        ptr_elite[2]=p[tamanio_de_poblacion-1][2];
    }
}

inline void genetic::suma_copias(double long p[][4]) //Funcion que realiza la suma de copias
{
    double suma=0;

    for (i=0;i<tamanio_de_poblacion;i++)
    {
        suma=p[i][2]+suma;
        p[i][3]=suma;
    }
}

void genetic::ruleta(double long p[][4],double long p2[][4])
{
    double r=0;

    for (j=0;j<tamanio_de_poblacion;j++)
    {
	   r=drand48()*(tamanio_de_poblacion);
        for (i=0;i<tamanio_de_poblacion;i++)
        {
            if(r>p[i][3])
            {
                //Se coloca la población y el fitness de la población deseada
                p2[j][0]=p[i][0];
                p2[j][1]=p[i][1];
            }
        }
    }
}

inline void genetic::codificacion(double long p2[][4])
{
    for (i=0;i<tamanio_de_poblacion;i++)
    {
        p2[i][2]=(p2[i][0]-limite_inferior)*((pow(2,bits)-1)/(limite_superior-limite_inferior));
    }
}

void genetic::cruce(double long p2[][4])
{
    unsigned long long int mascara_unos=pow(2,bits)-1;

    for (unsigned long int g=0;g<numero_de_generaciones;g++)
    {

        for (i=0;i<tamanio_de_poblacion/2;i++)
        {

            if(probabilidad_de_cruce>drand48())
            {

                unsigned long int punto_cruce=round((1) + ((bits-1)-(1))*drand48());
                unsigned long long int mascara_cola=pow(2,punto_cruce)-1;
                unsigned long long int mascara_cabeza=(mascara_unos^mascara_cola);

                unsigned long long int cabeza=((unsigned long long int)p2[i][2]&mascara_cabeza);
                unsigned long long int cola=((unsigned long long int)p2[i][2]&mascara_cola);
                unsigned long long int cabeza2=((unsigned long long int)p2[i+tamanio_de_poblacion/2][2]&mascara_cabeza);
                unsigned long long int cola2=((unsigned long long int)p2[i+tamanio_de_poblacion/2][2]&mascara_cola);

                p2[i][2]=(cabeza2|cola);
                p2[i+tamanio_de_poblacion/2][2]=(cabeza|cola2);

                mutacion(i,mascara_unos,matriz2);

                //Regreso el individuo nuevo decodificado
                p2[i][0]=limite_inferior+(p2[i][2]*((limite_superior-limite_inferior)) / (pow(2,bits)-1));
                p2[i+tamanio_de_poblacion/2][0]=limite_inferior+(p2[i+tamanio_de_poblacion/2][2]*((limite_superior-limite_inferior)) / (pow(2,bits)-1));

                //Regreso el nuevo fitness nuevo
                p2[i][1]=fitness<double>(p2[i][0]);
                p2[i+tamanio_de_poblacion/2][1]=fitness<double>(p2[i+tamanio_de_poblacion/2][0]);
            }
        }

        if (p2[tamanio_de_poblacion][1]<elite[1])
        {

            p2[tamanio_de_poblacion-1][0]=elite[0];
            p2[tamanio_de_poblacion-1][1]=elite[1];
            p2[tamanio_de_poblacion-1][2]=elite[2];

        }

        ordenamiento_burbuja(matriz2,elite);
    }
}


inline void genetic::mutacion(unsigned long long int i, unsigned long long int mascara_unos,double long p2[][4])
{

    unsigned long long int counter=1;

    while(counter<=bits)
    {
        if (probabilidad_de_mutacion>drand48())
        {
            if (0.5>drand48())
            {
                p2[i][2]=((unsigned long long int)p2[i][2]|((unsigned long long  int)pow(2,(counter-1))));
            }
            else
            {
                p2[i][2]=((unsigned long long int)p2[i][2]&((unsigned long long int)mascara_unos^((unsigned long long int)pow(2,(counter-1)))));
            }
        }

        if (probabilidad_de_mutacion>drand48())
        {
            if (0.5>drand48())
            {
                p2[i+tamanio_de_poblacion/2][2]=((unsigned long long int)p2[i+tamanio_de_poblacion/2][2]|((unsigned long long int)pow(2,(counter-1))));
            }
            else
            {
                p2[i+tamanio_de_poblacion/2][2]=((unsigned long long int)p2[i+tamanio_de_poblacion/2][2]&((unsigned long long int)mascara_unos^((unsigned long long int)pow(2,(counter-1)))));
            }
        }
        counter++;
    }
}

double long genetic::mejor_individuo()
{
    return matriz2[tamanio_de_poblacion-1][0];
}

genetic::~genetic()
{
    //dtor
}
