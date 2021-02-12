#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mgmres.c"

void fill_minus_one(int*,int);
void check_poisson(FILE *, FILE *, double, double, double, int, int, int, int, int, int, int, int, int);
int fun_j(int, int);
int fun_i(int, int);
double fun_ro1(double, double, double, double, double);
double fun_ro2(double, double, double, double, double);
double calc_epsl(int, int, int, int);
int check_Dirichlet(FILE * , FILE *, double, double, double, double *, int *, int *, double *,int, int, int, int, double, double, double, double, int);
int Eadge_cases(double, double, double, double *, int *, int *, double *,int, int, int, int, double, double, double, double,int);
void algebraization(FILE * , double, double, double, int, int, int, int, int, int, int, int);



void fill_minus_one(int* tab,int n){
    for(int i = 0; i < n ; i++){
        tab[i] = -1;
    }
}
void check_poisson(FILE * matrix, FILE * vector, double delta, double xmax, double ymax, int nx, int ny, int eps1, int eps2, int V1, int V2, int V3, int V4, int N){
    double a[5 * N];
    int ja[5 * N];
    int ia[5 * N];
    fill_minus_one(ia, 5*N);

    double b[N];
    double V[N];

    int nz_num = check_Dirichlet( matrix, vector, delta, xmax, ymax, a, ia, ja, b, nx, ny, eps1, eps2, V1, V2, V3, V4, N);
    int itr_max = 500;
    int mr = 500;
    double tol_abs = pow(10,-8);
    double tol_rel = pow(10,-8);

    pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);
}


/*Wzory 12 i 13 na indeksy*/
int fun_j(int nx, int l){
    return l / (nx + 1);
}
int fun_i(int nx, int l){
    return l - fun_j(nx,l) * (nx + 1);
}

//Funkcja obliczajca ro(1)  wzor (25)
double fun_ro1(double x, double y, double xmax, double ymax, double sigma){

    return exp(-1*pow(x - 0.25*xmax, 2)/pow(sigma, 2) - pow(y - 0.5*ymax, 2)/pow(sigma, 2));

}

//Funkcja obliczajaca ro(2) wzor (26)
double fun_ro2(double x, double y, double xmax, double ymax, double sigma){
    
    return -1*exp(-1*pow(x - 0.75*xmax, 2)/pow(sigma, 2) - pow(y - 0.5*ymax, 2)/pow(sigma, 2));

}
//Funkcja obliczajaca epsilon_l wzor (21)
double calc_epsl(int eps1, int eps2, int nx, int l){

    if(fun_i(nx,l) <= nx/2)
        return eps1;
    else
        return eps2;

}

int check_Dirichlet(FILE * matrix, FILE * vector, double delta, double xmax, double ymax, double *a, int *ia, int *ja, double *b,int nx, int ny, int eps1, int eps2, double V1, double V2, double V3, double V4, int N){

    int k = -1;

    for(int l = 0; l < N; l++){
        int brzeg =0; // wskaźnik położenia : 0 - środek obszaru ; 1 - brzeg
        double vb =0; // potencjal na brzegu

        if(fun_i(nx, l) == 0){
			brzeg = 1;
			vb = V1;
		}

        else if(fun_i(nx, l) == nx){
			brzeg = 1;
			vb = V3;
		}

        else if(fun_j(nx, l) == ny){
			brzeg = 1;
			vb = V2;
		}

        else if(fun_j(nx, l) == 0){
			brzeg = 1;
			vb = V4;
		}

        b[l] =  (-1)*(fun_ro1(delta*fun_i(nx,l), delta*fun_j(nx,l),xmax,ymax,xmax/10) + fun_ro2(delta*fun_i(nx,l), delta*fun_j(nx,l),xmax,ymax,xmax/10)); 
    
        if(brzeg == 1)
			b[l] = vb;

        ia[l] = -1;

        if(l - nx - 1 > 0 && brzeg == 0){
			k++;
			if(ia[l] < 0)
                ia[l] = k;

			a[k] = calc_epsl(nx,l,eps1,eps2)/(delta*delta);
			ja[k] = l - nx - 1;
		}

        if(l-1 > 0 && brzeg == 0) {
			k++;
			if(ia[l] < 0)
                ia[l] = k;
			a[k] = calc_epsl(nx,l,eps1,eps2)/(delta*delta);
			ja[k] = l-1;
		}

        k++;
		if(ia[l] < 0)
            ia[l] = k;

		if(brzeg == 0)
			a[k] = -(2*calc_epsl(nx,l,eps1,eps2)+calc_epsl(nx,l+1,eps1,eps2)+calc_epsl(nx,l+nx+1,eps1,eps2))/(delta*delta);
		else
			a[k] = 1;

		ja[k] = l;

		if(l < N && brzeg == 0){
			k++;
			a[k] = calc_epsl(nx,l+1,eps1,eps2)/(delta*delta);
			ja[k] = l + 1;
		}

		if(l < N-nx-1 && brzeg == 0){
			k++;
			a[k] = calc_epsl(nx,l+nx+1,eps1,eps2)/(delta*delta);
			ja[k] = l + nx + 1;
		}
        if(l%5 == 0 && l != 0)
            fprintf(vector, "\n");
        fprintf(vector, "%d %d %d %f \n", l, fun_i(nx,l), fun_j(nx,l), b[l]);

    }


	ia[N] = k+1;
    for(int i = 0; i < 5*N; i++)
        fprintf(matrix,"%d %0.f \n", i, a[i]);

    return k+1;
}

//Warunki Brzegowe Dirichleta 
int Eadge_cases(double delta, double xmax, double ymax, double *a, int *ia, int *ja, double *b,int nx, int ny, int eps1, int eps2, double V1, double V2, double V3, double V4, int N){
    int k = -1;

    for(int l = 0; l < N; l++){
        int brzeg =0; // wskaźnik położenia : 0 - środek obszaru ; 1 - brzeg
        double vb =0; // potencjal na brzegu

        if(fun_i(nx, l) == 0){
			brzeg = 1;
			vb = V1;
		}

        else if(fun_i(nx, l) == nx){
			brzeg = 1;
			vb = V3;
		}

        else if(fun_j(nx, l) == ny){
			brzeg = 1;
			vb = V2;
		}

        else if(fun_j(nx, l) == 0){
			brzeg = 1;
			vb = V4;
		}

        b[l] =  (-1)*(fun_ro1(delta*fun_i(nx,l), delta*fun_j(nx,l),xmax,ymax,xmax/10) + fun_ro2(delta*fun_i(nx,l), delta*fun_j(nx,l),xmax,ymax,xmax/10)); 
    
        if(brzeg == 1)
			b[l] = vb;

        ia[l] = -1;

        if(l - nx - 1 > 0 && brzeg == 0){
			k++;
			if(ia[l] < 0)
                ia[l] = k;

			a[k] = calc_epsl(nx,l,eps1,eps2)/(delta*delta);
			ja[k] = l - nx - 1;
		}

        if(l-1 > 0 && brzeg == 0) {
			k++;
			if(ia[l] < 0)
                ia[l] = k;
			a[k] = calc_epsl(nx,l,eps1,eps2)/(delta*delta);
			ja[k] = l-1;
		}

        k++;
		if(ia[l] < 0)
            ia[l] = k;

		if(brzeg == 0)
			a[k] = -(2*calc_epsl(nx,l,eps1,eps2)+calc_epsl(nx,l+1,eps1,eps2)+calc_epsl(nx,l+nx+1,eps1,eps2))/(delta*delta);
		else
			a[k] = 1;

		ja[k] = l;

		if(l < N && brzeg == 0){
			k++;
			a[k] = calc_epsl(nx,l+1,eps1,eps2)/(delta*delta);
			ja[k] = l + 1;
		}

		if(l < N-nx-1 && brzeg == 0){
			k++;
			a[k] = calc_epsl(nx,l+nx+1,eps1,eps2)/(delta*delta);
			ja[k] = l + nx + 1;
		}
        
    }


	ia[N] = k+1;
    return k+1;
}

//Algebrazacja rownania
void algebraization(FILE * map, double delta, double xmax, double ymax, int nx, int ny, int eps1, int eps2, int V1, int V2, int V3, int V4){

    
    int N = (nx+1)*(ny+1);   
    
    double *a = malloc((5*N) * sizeof(double));
    int *ja = malloc((5*N) * sizeof(int));
    int *ia = malloc((N + 1) * sizeof(int)); 
    double *b = malloc(N * sizeof(double));
    double *V = malloc(N * sizeof(double));
   

    int nz_num = Eadge_cases(delta, xmax, ymax, a, ia, ja, b, nx, ny, eps1, eps2, V1, V2, V3, V4,N);

    int itr_max = 500;
    int mr = 500;
    double tol_abs = pow(10,-8);
    double tol_rel = pow(10,-8);

    pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);


    double temp = 0.0;

	for(int l = 0; l < N; ++l){

		if( delta*fun_j(nx,l) > temp)
			fprintf(map,"\n");

		fprintf(map,"%f %f %f \n", delta*fun_j(nx,l), delta*fun_i(nx,l), V[l]);
		temp = delta*fun_j(nx,l);
	}

    free(ja);
    free(ia);
    free(a);
    free(b);
    free(V);
}