/**
 * \file main.c
 * \brief The main function of oblique shock reflection simulation(using rma).
 * \author Lin li
 * \version 3.0
 * \date 2021-07-23
 */

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <mpi.h>
#include "decomp1d.h"
#include "solve.h"

double*** a_memory(int k, int i, int j);
void f_memory(double*** U);
int output(double ***U, double dx, double dy, double gam);

int main(int argc, char *argv[])
{
    //calculation area
    const double x = 4.0, y = 1.0;

    //physical constants
    const double sf = 0.8, gam = 1.4; 

    //time of simulation
    const double total_time = 3.0;

    double dt, start, finish;
    double current_time = 0;
    const double dx = x / gx;
    const double dy = y / gy;
    int rank, size;
    int left, right, s, e;

    double ***U = a_memory(4, (gx + 2), (gy + 2));
    double ***U_half = a_memory(4, (gx + 2), (gy + 2));
    double ***F = a_memory(4, (gx + 2), (gy + 2));
    double ***G = a_memory(4, (gx + 2), (gy + 2));
    double ***temp = a_memory(4, (gx + 2), (gy + 2));

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

    left  = rank - 1;
    right = rank + 1;

    if( rank == 0 ){
        left = MPI_PROC_NULL;
    }

    if( rank == size -1 ){
        right  = MPI_PROC_NULL;
    } 

    decomp1d(rank, size, &s, &e); 
    printf("(rank: %d) gx: %d s: %3d; e: %3d; left: %2d; right: %2d\n", rank, gx, s, e, left, right);

    MPI_Win win;
    
    MPI_Win_create(&U[0][s - 1][0], ((e - s + 3) * (gy + 2) + 3 * (gx + 2) * (gy + 2))* sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);

	start = MPI_Wtime();

    init(U, dx, dy, gam, rank, s, e);  //initialization

    if(rank == 0)
    {
        cfl(U, dx, dy, &dt, sf, gam);  //calculate dt using the cfl number
    }
    MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //printf("(myid: %d) dt = %lf\n", rank, dt);
    
    while (current_time < total_time) 
	{
        solver_2d(U, U_half, F, G, temp, dx, dy, dt, gam, s, e, left, right, rank, size, win);  //use a second order Lax-Wendroff two-step difference scheme
        
        current_time = current_time + dt;
        //printf("time = %10g\n", current_time);  //the simulation process should be commented out when testing the runing time
    }    

	finish = MPI_Wtime();

    gathergrid(U, rank, size, s, e, left, right);

    if(rank == 0)
    {
        printf("Program done, the parallel version use %.3lf seconds.\n", finish -start);
        output(U, dx, dy, gam);
    }

    f_memory(U);
    f_memory(U_half);
    f_memory(F);
    f_memory(G);
    f_memory(temp);

    MPI_Finalize();
    return 0;
}

double*** a_memory(int k, int i, int j)
{	
	double*** ppp = (double***)malloc(k * sizeof(double**));
	double** pp = (double**)malloc(k * i * sizeof(double*));
	double* p = (double*)malloc(k * i * j * sizeof(double));
	
	for(int m = 0; m < k; m++)
	{
		ppp[m] = pp + m * i;
		for(int n = 0; n < i; n++)
		{
			ppp[m][n] = p + n * j;
		}
		p = p + i * j;
	}
	
    return ppp;
}

void f_memory(double*** U)
{
    if(*U != NULL)
    {
        if(**U != NULL)
        {
            free(**U);
            **U = NULL;
        }
        free(*U);
        *U = NULL;
    }
    
    free(U);
    U = NULL;
}

int output(double ***U, double dx, double dy, double gam)
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
