/**
 * \file main.c
 * \brief The main function of oblique shock reflection simulation(serial version).
 * \author Lin li
 * \version 2.0
 * \date 2021-07-19
 */

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "solve.h"

double*** a_memory(int a, int b, int c);
void f_memory(double*** U, int k, int i);
int output(double ***U , double dx, double dy, double gam);

int main(void)
{
    //calculation area
    const double x = 4.0, y = 1.0;

    //physical constants
    const double sf = 0.8, gam = 1.4; 

    //time of simulation
    const double total_time = 3.0;

    clock_t start,finish;
    double dt, runtime;
    double current_time = 0;
    const double dx = x / gx;
    const double dy = y / gy;

    double ***U = a_memory(4, (gx + 2), (gy + 2));
    double ***U_half = a_memory(4, (gx + 2), (gy + 2));
    double ***F = a_memory(4, (gx + 2), (gy + 2));
    double ***G = a_memory(4, (gx + 2), (gy + 2));
    double ***temp = a_memory(4, (gx + 2), (gy + 2));

	start = clock();

    init(U, dx, dy, gam);  //initialization

    cfl(U, dx, dy, &dt, sf, gam);  //calculate dt using the cfl number

    while (current_time < total_time) 
	{
        solver_2d(U, U_half, F, G, temp, dx, dy, dt, gam);  //use a second order Lax-Wendroff two-step difference scheme

        current_time = current_time + dt; 
        //printf("time = %10g\n", current_time);  //the simulation process should be commented out when testing the runing time
    }

	finish = clock();

	runtime = (finish - start) / CLOCKS_PER_SEC;
	printf("Program done, the serial version use %.3lf seconds.\n", runtime);

    output(U, dx, dy, gam);

    f_memory(U, 4, (gx+2));
    f_memory(U_half, 4, (gx+2));
    f_memory(F, 4, (gx+2));
    f_memory(G, 4, (gx+2));
    f_memory(temp, 4, (gx+2));

    return 0;
}

//allocate memory
double*** a_memory(int a, int b, int c)
{
    double ***U = malloc(a * sizeof(double **));

    for(int k = 0; k < a; k++)
    {
        U[k] = malloc(b * sizeof(double *));
        for(int i = 0; i < b; i++)
        {
            U[k][i] = malloc(c * sizeof(double));
        }
    }

    return U;
}

//free memory
void f_memory(double*** U, int k, int i)
{
    for(int k = 0; k < 4; k++)
    {
        for(int i = 0; i< (gx + 2); i++)
        {
            free(U[k][i]);
        }
    }

    for(int k = 0; k < 4; k++)
    {
        free(U[k]);
    }

    free(U);
}

int output(double ***U , double dx, double dy, double gam)
{
    FILE* fp;
    double rou, u, v, p;
    
    //for density contour map
    fp = fopen("result.txt", "w");
    if(fp == NULL)
    {
        perror("Error opening result.txt for writing"); 
        return(-1);
    } 
    
    fprintf(fp, "variables = x, y, rou, u, v, p, E\n");
    fprintf(fp, "gx = %d, gy = %d\n", gx, gy);

    for (int i = 1; i <= gx; i++)
	{
        for (int j = 1; j <= gy; j++) 
		{
            rou = U[0][i][j];
            u = U[1][i][j] / rou;
            v = U[2][i][j] / rou;
            p = (gam - 1) * (U[3][i][j] - 0.5 * U[0][i][j] * (u * u + v * v));
            fprintf(fp, "%10lf%10lf%10lf%10lf%10lf%10lf%10lf\n", i * dx, j * dy, rou, u, v, p, U[3][i][j]);
        }
	}
    fclose(fp);

	//for pressure map at position y / 2
    fp = fopen("result_0.5y.txt", "w");
    if(fp == NULL)
    {
        perror("Error opening result_0.5y.txt for writing"); 
        return(-1);
    } 

    fprintf(fp, "variables = x, rou, u, v, p, E\n");
    fprintf(fp, "gx = %d, gy = %d, 0.5y position\n", gx, gy);

    int j = gy / 2;
    for (int i = 1; i <= gx; i++)
    {
        rou = U[0][i][j];
        u = U[1][i][j] / rou;
        v = U[2][i][j] / rou;
        p = (gam - 1) * (U[3][i][j] - 0.5 * U[0][i][j] * (u * u + v * v));
        fprintf(fp, "%10lf%10lf%10lf%10lf%10lf%10lf\n", i * dx, rou, u, v, p, U[3][i][j]);
    }
    fclose(fp);

    return 0;
}
