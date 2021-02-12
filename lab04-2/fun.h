#pragma once
#include <stdlib.h>
#include <stdio.h>
#define epsilon 1.
#define delta 0.1
#define nx 150
#define ny 100
#define v1 10
#define v2  0
#define xMax (delta * nx)
#define yMax (delta * ny)
#define TOL pow(10, -8)
#define sigmaX (0.1 * xMax)
#define sigmaY (0.1 * yMax)

void zeros(double V[0][ny + 1]);
void fill_ro();
double stop_condition(double V[][ny + 1]);
void eadge_cases(double V[][ny + 1]);
void copy(double V_s[][ny + 1], double V_n[][ny + 1]);
void global_relaxation_core(double omegaG, double V_s[][ny + 1], double V_n[][ny + 1]);
void local_relaxation_core(double omegaL, double V[][ny + 1]);
double err(double V[][ny + 1], int i, int j);
void global_relaxation(double omegaG, FILE * integral, FILE * map, FILE * error);
void local_relaxation(double omegaL, FILE * integral);
