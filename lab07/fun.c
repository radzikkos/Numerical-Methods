#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//Stałe
#define delta  0.01
#define ro  1.0
#define mi  1.0
#define n_x  200
#define n_y  90
#define i_1  50
#define j_1  55
#define it_max  20000


void navier_stokes(double, FILE * );
double q_out(double, double*);
double calc_gamma(int , double psi[][n_y+1], double zeta[][n_y+1]);
void eage_case_psi(double, double,double *, double psi[][n_y+1]);
void eadge_case_zeta(double, double,double *, double psi[][n_y+1], double zeta[][n_y+1]);
int check_if_border(int, int);
int check_if_body(int, int);
void relaxation_algorithm(double, double, double *, double psi[][n_y+1],double zeta[][n_y+1]);




void navier_stokes(double q_we, FILE * fp){

    double psi[n_x+1][n_y+1];
    double zeta[n_x+1][n_y+1];
    for( int i =0; i<=n_x; i++ ){
        for( int j=0; j<=n_y; j++ ){
            psi[i][j] = 0.0;
            zeta[i][j] = 0.0;
        }
    }
    double y[n_y+1];
    for( int i=0; i<=n_y; i++ )
        y[i] = delta*i;

    double q_wy = q_out(q_we, y);

    
    eage_case_psi(q_we,q_wy,y,psi);
    relaxation_algorithm( q_we, q_wy, y, psi, zeta);


    double u = 0;
    double v = 0;
    for(int i = 1; i <= n_x - 1; i++){

        for(int j = 1; j <= n_y - 1; j++){

            if(!check_if_border(i,j) && check_if_body(i,j)){

                u = (psi[i][j + 1] - psi[i][j - 1]) / (2 * delta);
                v = -(psi[i + 1][j] - psi[i - 1][j]) / (2 * delta);
            }
            else{

                u = 0;
                v = 0;
            }

            fprintf(fp,"%f %f %f %f %f %f\n", delta*i, delta*j, psi[i][j], zeta[i][j], u, v );
    
        }
        fprintf(fp, "\n");
    }

    
}


double q_out(double q_we, double *y){
    return  q_we * (pow(y[n_y], 3) - pow(y[j_1], 3) - 3 * pow(y[n_y], 2) * y[j_1] + 3 * pow(y[j_1], 2) * y[n_y]) / (pow(y[n_y], 3));
}


double calc_gamma(int j_2, double psi[][n_y+1], double zeta[][n_y+1]){
    double gamma = 0.0;

    for(int i = 1; i <= n_x - 1; i++)
        gamma += (psi[i + 1][j_2] + psi[i - 1][j_2] + psi[i][j_2 + 1] + psi[i][j_2 - 1] - 4 * psi[i][j_2] - pow(delta, 2) * zeta[i][j_2]);
        
    return gamma;
}


void eage_case_psi(double q_we, double q_wy, double *y, double psi[][n_y+1]){

    //A
    for(int j = j_1; j <= n_y; j++)    
        psi[0][j] = q_we / (2 * mi) * (y[j] * y[j] * y[j] / 3 - y[j] * y[j] / 2 * (y[j_1] + y[n_y]) + y[j] * y[j_1] * y[n_y]);
    
    //C
    for(int j = 0; j <= n_y; j++)   
        psi[n_x][j] = q_wy / (2 * mi) * (y[j] * y[j] * y[j] / 3 - y[j] * y[j] / 2 * y[n_y]) +  q_we * y[j_1] * y[j_1] * (-y[j_1] + 3 * y[n_y]) / (12 * mi);
    
    //B
    for(int i = 1; i <= n_x - 1; i++)  
        psi[i][n_y] = psi[0][n_y];
    
    //D
    for(int i = i_1; i <= n_x - 1; i++)
        psi[i][0] = psi[0][j_1];
    
    //E
    for(int j = 1; j <= j_1; j++)
        psi[i_1][j] = psi[0][j_1];
    
     //F
    for(int i = 1; i <= i_1; i++)
        psi[i][j_1] = psi[0][j_1];
    
}


void eadge_case_zeta(double q_we, double q_wy, double *y, double psi[][n_y+1], double zeta[][n_y+1]){

    //A
    for(int j = j_1; j <= n_y; j++)
        zeta[0][j] = q_we / (2 * mi) * (y[j] * 2 - y[j_1] - y[n_y]);
    
    //C
    for(int j = 0; j <= n_y; j++)  
        zeta[n_x][j] = q_wy / (2 * mi) * (2 * y[j] - y[n_y]);
    
    //B
    for(int i = 1; i <= n_x - 1; i++)
        zeta[i][n_y] = 2 / (pow(delta, 2)) * (psi[i][n_y - 1] - psi[i][n_y]);
    
    //D
    for(int i = i_1 + 1; i <= n_x - 1; i++)  
        zeta[i][0] = 2 / (pow(delta, 2)) * (psi[i][1] - psi[i][0]);
    
    //E
    for(int j = 1; j <= j_1 - 1; j++)
        zeta[i_1][j] = 2 / (pow(delta, 2)) * (psi[i_1 + 1][j] - psi[i_1][j]);
    
    //F
    for(int i = 1; i <= i_1; i++)
        zeta[i][j_1] = 2 / (delta * delta) * (psi[i][j_1 + 1] - psi[i][j_1]);
    
    //E/F
    zeta[i_1][j_1] = 0.5 * (zeta[i_1 - 1][j_1] + zeta[i_1][j_1 - 1]);
}

void relaxation_algorithm(double q_we, double q_wy, double *y, double psi[][n_y+1], double zeta[][n_y+1] ){
    int j_2 = j_1 + 2;
    double omega;

    for(int it = 1; it < it_max; it++)
    {
        
        if(it < 2000)
            omega = 0.0;
        
        else   
            omega = 1.0;
        
        for(int i = 1; i <= n_x - 1; i++)
        {
            for(int j = 1; j <= n_y - 1; j++)
            {
                if(!check_if_border(i,j) && check_if_body(i,j))
                {
                    psi[i][j] = 0.25 * (psi[i + 1][j] + psi[i - 1][j] + psi[i][j + 1] + psi[i][j - 1] - delta * delta * zeta[i][j]);

                    zeta[i][j] = 0.25 * (zeta[i + 1][j] + zeta[i - 1][j] + zeta[i][j + 1] + zeta[i][j - 1]) -
                        omega * ro / (16 * mi) * 
                        ((psi[i][j + 1] - psi[i][j - 1]) * (zeta[i + 1][j] - zeta[i - 1][j]) - (psi[i + 1][j] - psi[i - 1][j]) * (zeta[i][j + 1] - zeta[i][j - 1]));
                }
            }
        }

        eadge_case_zeta(q_we, q_wy, y, psi, zeta);
        double gamma = calc_gamma(j_2, psi, zeta);
        printf("Gamma: %f \n", gamma);
    }
}


int check_if_border(int i, int j){

    if(i==0 || i == n_x)
        return 1;
    if(j==0 || j == n_y)
        return 1;
    if(i == i_1 && j <=j_1)
        return 1;
    if(i <= i_1 && j == j_1)
        return 1;
    else
        return 0;
}


int check_if_body(int i, int j){

    if(i <= i_1 && j < j_1)
        return 0;
    else
        return 1;
    
}