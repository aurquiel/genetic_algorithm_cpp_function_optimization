#ifndef GENETIC_H
#define GENETIC_H
#include <string>

class genetic
{
    public:
        genetic(std::string f_fitness_variable_string,bool f_maximizar_o_minimiza, double long f_limite_superior,
                 double long f_limite_inferior, long long int f_presicion, double long f_probabilidad_de_cruce,
                 double long f_probabilidad_de_mutacion,unsigned long long int f_tamanio_de_poblacion,
                 unsigned long long int f_numero_de_generaciones);

        double long mejor_individuo(); //Funcion que devuelve el mejor individuo

        virtual ~genetic();//Desctructor

    protected:

    private:
        //Variables externos del programa los toma en el construcctor
        bool maximizar_o_minimizar;
        unsigned long long int presicion;
        double long limite_superior;
        double long limite_inferior;
        double long probabilidad_de_cruce;
        double long probabilidad_de_mutacion;
        unsigned long long int tamanio_de_poblacion;
        unsigned long long int numero_de_generaciones;
        //Final

        //Variables internas del programa
        unsigned long long int i,j; //variable contadora
        unsigned long long int bits; //Varaible del numero de bits necesarios para el intervalo de numeros mas presicion
        std::string fitness_variable_string=" ";
        double long media;  //Variable para calcular la media suma de fitness entre fitness
        double long matriz[50000][4];  //Matriz primera poblacion inicial; individuos,fitness,copias, suma de copias
        double long matriz2[50000][4]; //Matriz de ruleta y operacion
        double long elite[3]={0}; //Arreglo para guardar al individuo elite
        //Final

        //Funciones privadas
        template <typename T> double long fitness(double long &valor);
        void generacion_de_poblacion_inicial(double long p[][4]); //Funcion que genera la probalción incial y se obtiene el fitness
        double long numero_copias(double long &valor);//Funcion para obtner el numero de copias
        void ordenamiento_burbuja(double long p[][4], double long *ptr_elite);
        void suma_copias(double long p[][4]);  //Función para obtner la suma de copias
        void ruleta(double long p[][4],double long p2[][4]); //Funcion que realiza la ruleta de seleccion
        void cruce(double long p2[][4]); //Funcion donde se relaiza el cruce de indivisuos
        void codificacion(double long p2[][4]); //Funcion que realiza la codificacion de los individuos
        void mutacion(unsigned long long int i, unsigned long long int mascara_unos,double long p2[][4]); //Funcion que realiza la mutacion
        //Final

};

#endif // GENETIC_
