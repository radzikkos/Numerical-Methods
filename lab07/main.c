#include <stdio.h>
#include "fun.c"
#include <stdlib.h>

int main(){

    //Q = -1000
    FILE * fp = fopen("case_-1000.dat", "w");
    double Q = -1000;
    navier_stokes(Q, fp);
    fclose(fp);

    
    //Q = -4000
    fp = fopen("case_-4000.dat", "w");
    Q = -4000;
    navier_stokes(Q, fp);
    fclose(fp);

    //Q = 4000
    fp = fopen("case_4000.dat", "w");
    Q = 4000;
    navier_stokes(Q, fp);
    fclose(fp);
    
    return 0;
}
