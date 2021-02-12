#include "fun.h"
#include <math.h>

/* Tablica gestosci */
double ro[nx + 1][ny + 1];


/*Zerowanie tablicy*/
void zeros(double V[0][ny + 1]){
    for(int i = 0 ; i <= nx; i++){
        for(int j = 0; j <= ny; j++){
            V[i][j] = 0.0;
        }
    }
}

/*Wypelnienie tablicy gestosci*/
void fill_ro(){
    for(int i = 0; i <= nx; i++){
        for(int j = 0; j <= ny; j++){
			ro[i][j] =  exp( -pow(delta*i - 0.35 * xMax, 2) / (pow(sigmaX,2)) - pow(delta*j - 0.5 * yMax, 2) / (pow(sigmaY,2)) );
            ro[i][j] += (-1.0)*exp( -pow(delta*i - 0.65 * xMax, 2) / (pow(sigmaX,2)) - pow(delta*j - 0.5 * yMax, 2) / (pow(sigmaY,2)) ) ;
            }
        }
}

// /*Warunek stopu*/
double stop_condition(double V[][ny + 1]){
    double S = 0;

    for(int i = 0; i < nx; i++) 
		for(int j = 0; j < ny; j++) 
				S += pow(delta,2) * ( 0.5*pow((V[i+1][j] - V[i][j])/delta, 2) +  0.5*pow((V[i][j+1] - V[i][j])/delta, 2) - ro[i][j]*V[i][j]);
                	
    return S;
}

/*Warunki brzegowe*/
void eadge_cases(double V[][ny + 1]){
    for(int i = 0; i <= nx; i++){
        V[i][0] = v1;
        V[i][ny] = v2;
    }
}

/*Kopiowanie tablicy V_n do V_s*/
void copy(double V_s[][ny + 1], double V_n[][ny + 1]){
    for(int i = 0 ; i <= nx; i++)
            for(int j = 0; j <= ny; j++)
                V_s[i][j] = V_n[i][j];
}

/*Kroki w relaksacji globalnej*/
void global_relaxation_core(double omegaG, double V_s[][ny + 1], double V_n[][ny + 1]){
    
    for(int i = 1; i < nx; i++)
        for(int j = 1; j < ny; j++)
            V_n[i][j] = 0.25 * ( V_s[i+1][j] + V_s[i-1][j] + V_s[i][j+1] + V_s[i][j-1] +  pow(delta,2)/epsilon * ro[i][j] );

    for(int j = 1; j < ny; j++){
        V_n[0][j] = V_n[1][j];
        V_n[nx][j] = V_n[nx - 1][j];
    }

    for(int i = 0 ; i <= nx ; i++)
        for(int j = 0; j <= ny; j++)
            V_s[i][j] = (1-omegaG)*V_s[i][j] + omegaG*V_n[i][j];  
}

/*Kroki w relaksacji lokalnej*/
void local_relaxation_core(double omegaL, double V[][ny + 1]){
    for(int i = 1; i < nx; i++)
            for(int j = 1; j < ny; j++)
                V[i][j] = (1 - omegaL)* V[i][j] +(omegaL/4.0) * ( V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1] +  pow(delta,2)/epsilon * ro[i][j] );

    for(int j = 1; j < ny; j++){
        V[0][j] = V[1][j];
        V[nx][j] = V[nx - 1][j];
    }
}

double err(double V[][ny + 1], int i, int j){
    return  (V[i+1][j] - 2*V[i][j] + V[i-1][j])/(pow(delta,2)) + (V[i][j+1] - 2*V[i][j] + V[i][j-1])/(pow(delta,2)) + ro[i][j]/epsilon ;
}

/*Relaksacja globalna dla podanej omegi, wypelnianie plikow integral, map, error*/
void global_relaxation(double omegaG, FILE * integral, FILE * map, FILE * error){

    double sIt[2] = {0.0, 0.0};

    double V_n[nx+1][ny+1];
    double V_s[nx+1][ny+1];

    fill_ro();

    zeros(V_n);
    
    eadge_cases(V_n);

    copy(V_s, V_n);

    sIt[1] = stop_condition(V_n);

    int it = 0;

    do{
    
        it+=1;

        global_relaxation_core(omegaG, V_s, V_n);

        sIt[1] = sIt[0];
        sIt[0] = stop_condition(V_s);

        
        fprintf(integral, "%d %f \n", it, sIt[0]);

    }while(fabs( (sIt[0] - sIt[1])/sIt[1] ) > TOL);


    for(int i = 0; i <= nx; i++){
        for(int j = 0; j <= ny; j++){
            
            fprintf(map, "%f %f %f \n", delta*i, delta*j, V_s[i][j]);

            
            if(i > 0 && j > 0 && i < nx && j < ny)
                fprintf(error, "%f %f %f \n", delta*i, delta*j, err(V_s, i, j) );
             else
                fprintf(error, "%f %f %f \n", delta*i, delta*j, 0.0);
        }
        fprintf(map, "\n");
        fprintf(error, "\n");

    }
}


/*Relaksacja lokalna dla podanej omegi, wypelnianie pliku integral*/
void local_relaxation(double omegaL, FILE * integral){

    double sIt[2] = {0.0, 0.0};

    double V[nx+1][ny+1];

    zeros(V);
    
    eadge_cases(V);

    sIt[1] = stop_condition(V);

    int it = 0;
    do{
        
        it+=1;

        local_relaxation_core(omegaL, V);        

        sIt[1] = sIt[0];
        sIt[0] = stop_condition(V);

        
        fprintf(integral, "%d %f \n", it, sIt[0]);
        
    }while(fabs( (sIt[0] - sIt[1])/sIt[1] ) > TOL);

}