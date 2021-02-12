#include "fun.h"

int main(){
    double omegaG;

    omegaG = 1.0;
    FILE * integral = fopen("integral_1.0_relglob.dat", "w");
    FILE * map = fopen("map_1.0_relglob.dat", "w");
    FILE * error = fopen("error_1.0_relglob.dat", "w");
    
    global_relaxation(omegaG, integral, map, error);

    fclose(integral);
    fclose(map);
    fclose(error);


    omegaG=0.6;
    integral = fopen("integral_0.6_relglob.dat", "w");
    map = fopen("map_0.6_relglob.dat", "w");
    error = fopen("error_0.6_relglob.dat", "w");

    global_relaxation(omegaG, integral, map, error);

    fclose(integral);
    fclose(map);
    fclose(error);


    double omegaL;

    omegaL = 1.0;
    integral = fopen("integral_1.0_relloc.dat", "w");
    local_relaxation(omegaL, integral);
    fclose(integral);

    omegaL = 1.4;
    integral = fopen("integral_1.4_relloc.dat", "w");
    local_relaxation(omegaL, integral);
    fclose(integral);

    omegaL = 1.8;
    integral = fopen("integral_1.8_relloc.dat", "w");
    local_relaxation(omegaL, integral);
    fclose(integral);

    omegaL = 1.9;
    integral = fopen("integral_1.9_relloc.dat", "w");
    local_relaxation(omegaL, integral);
    fclose(integral);

    return 0;
}