#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define delta  0.2
#define nx  128
#define ny  128
const double xMax = nx * delta;
const double yMax = ny * delta;
const double TOL = std::pow(10,-8);

double stop_condition(double** V,int k){
    double S = 0.0;

    for(int i = 0; i <= nx-k; i += k){
		for(int j = 0; j <= ny-k; j += k){
            double part1 = 0.5*std::pow(k*delta, 2);
            double part2 = std::pow((V[i+k][j] - V[i][j])/(2*k*delta) + (V[i+k][j+k] - V[i][j+k])/(2*k*delta), 2);
            double part3 = std::pow((V[i][j+k] - V[i][j])/(2*k*delta) + (V[i+k][j+k] - V[i+k][j])/(2*k*delta), 2);
			S += part1 * (part2 + part3);
		}
	}
    return S;
}

void edge_cases(double** V){
    for(int y=0; y <= ny;y++)
        V[0][y] = sin(M_PI * delta*y / yMax);
    
    for(int x=0; x <= nx;x++)
        V[x][ny] = (-1.0) * sin(2*M_PI * delta*x / xMax);

    for(int y=0;y <= ny;y++)
        V[nx][y] = sin(M_PI * delta*y / yMax);

    for(int x=0;x <= nx;x++)
        V[x][0] =  sin(2*M_PI* delta*x/xMax);

}

void condensing(double** V,int k){

    for(int i = 0; i <= nx-k; i += k){
		for(int j = 0; j <= ny-k; j += k){

            V[i + k/2][j + k/2] = 0.25*(V[i][j] + V[i+k][j] + V[i][j+k] + V[i+k][j+k]);

            if(i!=nx-k)
                V[i + k][j + k/2] = 0.5*(V[i+k][j] + V[i+k][j+k]);
            if(j!=ny-k)
                V[i + k/2][j + k] = 0.5*(V[i][j+k] + V[i+k][j+k]);
            if(j!=0)
                V[i + k/2][j] = 0.5*(V[i][j] + V[i+k][j]);
            if(i!=0)
                V[i][j + k/2] = 0.5*(V[i][j] + V[i][j+k]);
        }
    }
}

double** zeros(){
    double** V = new double*[nx+1];
    for(int i = 0 ; i < nx+1; i++){
        V[i] = new double[ny+1];
        for(int j = 0; j < ny+1; j++){
            V[i][j] = 0.;
        }
    }
        
    return V;
}
void fr(double **V){
    for(int i = 0 ; i < nx+1; i++)
            delete[] V[i];
    delete[] V;
}
void discretization(double** V, int k){
    for(int i = k; i <= nx-k; i += k)
                for(int j = k; j <= ny-k; j += k)
                    V[i][j] = 0.25*(V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k]);

}
int main(){
    FILE * map = fopen("maps.dat" , "w");
    FILE * integral = fopen("integral.dat", "w");

    double s[2] = {0.0, 0.0};
    int it= 0;
    double** V = zeros();
    
    edge_cases(V);

    for(int k = 16; k>0; k=k/2){

        s[1] = stop_condition(V,k);

        do{
            discretization(V,k);
            s[0] = s[1];
            s[1] = stop_condition(V, k);
            fprintf(integral, "%d %d %f \n", k, it, s[1]);
            it++;

        }while(fabs((s[1] - s[0])/s[0]) > TOL);

        fprintf(integral, "\n \n");

        //Mapping
        for(int i=0;i <= nx;i+=k){
            for(int j=0;j <= ny;j+=k){
                fprintf(map, "%f %f %f \n", delta*i, delta*j, V[i][j]);
            }
            fprintf(map,"\n");
        }
        fprintf(map, "\n \n");
        condensing(V, k);

    }

    fclose(map);
    fclose(integral);
    fr(V);
    return 0;
}


