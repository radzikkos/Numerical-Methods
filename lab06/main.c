#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fun.c"


#define delta 0.1


int main(){
    int nx = 4;
    int ny = 4;
    double epsilon_1 = 1;
    double epsilon_2 = 1;
    int V1 = 10;
    int V3 = 10;
    int V2 = -10;
    int V4 = -10;
    double xmax = 0;
    double ymax = 0;
    double ro = 0;
    int N = (nx+1)*(ny+1);

    FILE * vector = NULL;
    FILE * matrix = NULL;
    FILE * map = NULL;

    vector = fopen("vector_check.dat", "w");
    matrix = fopen("matrix_check.dat", "w");
    check_poisson(matrix,vector,delta,xmax,ymax,nx,ny,epsilon_1,epsilon_2,V1,V2,V3, V4,N);
    fclose(vector);
    fclose(matrix);

    /*Mapy potencjalow*/
    //1.a
    puts("podpunkt 1.a");
    nx = 50;
    ny = 50;

    map = fopen("map1a.dat", "w");
    algebraization(map, delta, xmax, ymax, nx, ny, epsilon_1, epsilon_2, V1, V2, V3, V4);
    fclose(map);

    //1.b
    puts("podpunkt 1.b");
    nx = 100;
    ny = 100;
   
    map = fopen("map1b.dat", "w");
    algebraization(map, delta, xmax, ymax, nx, ny, epsilon_1, epsilon_2, V1, V2, V3, V4);
    fclose(map);

    //1.c
    puts("podpunkt 1.c");
    nx = 200;
    ny = 200;
    
    map = fopen("map1c.dat", "w");
    algebraization(map, delta, xmax, ymax, nx, ny, epsilon_1, epsilon_2, V1, V2, V3, V4);
    fclose(map);


    /** Rozklady potencjalow **/
    
    nx = ny = 100;
	V1=V2=V3=V4=0;
	xmax = delta*nx;
	ymax = delta*ny;

    //2.a
    puts("podpunkt 2.a");

    epsilon_1 = 1;
    epsilon_2 = 1;

    map = fopen("map2a.dat", "w");
    algebraization(map, delta, xmax, ymax, nx, ny, epsilon_1, epsilon_2, V1, V2, V3, V4);
    fclose(map);

    //2.b
    puts("podpunkt 2.b");
    epsilon_1 = 1;
    epsilon_2 = 2;

    map = fopen("map2b.dat", "w");
    algebraization(map, delta, xmax, ymax, nx, ny, epsilon_1, epsilon_2, V1, V2, V3, V4);
    fclose(map);

    //2.c
    puts("podpunkt 2.c");
    epsilon_1 = 1;
    epsilon_2 = 10;

    map = fopen("map2c.dat", "w");
    algebraization(map, delta, xmax, ymax, nx, ny, epsilon_1, epsilon_2, V1, V2, V3, V4);
    fclose(map);

}