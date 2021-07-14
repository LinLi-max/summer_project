/**
 * \file main.c
 * \brief The main function of oblique shock reflection simulation.
 * \author Lin li
 * \version 1.0
 * \date 2021-07-12
 */

#include <stdio.h>
#include <time.h>
#include <mpi.h>
#include "decomp1d.h"
#include "solve.h"

int output(double U[gx + 2][gy + 2][4], double dx, double dy, double gam);

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
    double U[gx + 2][gy + 2][4], U_half[gx + 2][gy + 2][4], F_G[gx + 2][gy + 2][4];

	int rank, size;
    int left, right, s, e;

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

    decomp1d(rank, size, &s, &e); // x向，将gx平分给各核心
    printf("(rank: %d) gx: %d s: %3d; e: %3d; left: %2d; right: %2d\n", rank, gx, s, e, left, right);

	start = MPI_Wtime();

    init(U, dx, dy, gam);  //initialization
    if(rank == 0)
    {
        cfl(U, dx, dy, &dt, sf, gam);  //calculate dt using the cfl number
    }
    MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    printf("(myid: %d) dt = %lf\n", rank, dt);
    
    while (current_time < total_time) 
	{
        solver_2d(U, U_half, F_G, dx, dy, dt, gam, s, e, left, right, MPI_COMM_WORLD);  //use a second order Lax-Wendroff two-step difference scheme
        printf("time = %10g\n", current_time);  //the simulation process should be commented out when testing the runing time
        current_time = current_time + dt; 
    }    

    gathergrid(U, rank, size, s, e, left, right);

    MPI_Barrier(MPI_COMM_WORLD);

	finish = MPI_Wtime();

    if(rank == 0)
    {
        printf("Program done, the parallel version use %.3lf seconds.\n", finish -start);
        output(U, dx, dy, gam);
    }
	
    MPI_Finalize();
    return 0;
}

int output(double U[gx + 2][gy + 2][4], double dx, double dy, double gam)
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
            rou = U[i][j][0];
            u = U[i][j][1] / rou;
            v = U[i][j][2] / rou;
            p = (gam - 1) * (U[i][j][3] - 0.5 * U[i][j][0] * (u * u + v * v));
            fprintf(fp, "%10lf%10lf%10lf%10lf%10lf%10lf%10lf\n", i * dx, j * dy, rou, u, v, p, U[i][j][3]);
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
        rou = U[i][j][0];
        u = U[i][j][1] / rou;
        v = U[i][j][2] / rou;
        p = (gam - 1) * (U[i][j][3] - 0.5 * U[i][j][0] * (u * u + v * v));
        fprintf(fp, "%10lf%10lf%10lf%10lf%10lf%10lf\n", i * dx, rou, u, v, p, U[i][j][3]);
    }
    fclose(fp);

    return 0;
}
