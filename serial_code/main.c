/**
 * \file main.c
 * \brief The main function of oblique shock reflection simulation.
 * \author Lin li
 * \version 2.0
 * \date 2021-06-29
 */

#include <stdio.h>
#include <time.h>
#include "solve.h"

void output(double U[gx + 2][gy + 2][4], double dx, double dy, double gam);

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
    double U[gx + 2][gy + 2][4], U_half[gx + 2][gy + 2][4], F_G[gx + 2][gy + 2][4];

	start=clock();

    init(U, dx, dy, gam);  //initialization
	
    while (current_time < total_time) 
	{
        cfl(U, dx, dy, &dt, sf, gam);  //calculate dt using the cfl number
        solver_2d(U, U_half, F_G, dx, dy, dt, gam);  //use a second order Lax-Wendroff two-step difference scheme

        current_time = current_time + dt;  
        printf("time=%10g dt=%10g\n", current_time, dt);  //the simulation process should be commented out when testing the runing time
    }

	finish=clock();

	runtime=(finish-start)/CLOCKS_PER_SEC;
	printf("Program done, the serial version use %.3lf seconds.\n", runtime);

    output(U, dx, dy, gam);

    return 0;
}

void output(double U[gx + 2][gy + 2][4], double dx, double dy, double gam)
{
    FILE* fp;
    double rou, u, v, p;
    
    //for density contour map
    fp = fopen("result.txt", "w");
    fprintf(fp, "variables = x, y, rou, u, v, p, E\n");
    fprintf(fp, "gx=%d, gy=%d\n", gx, gy);

    for (int i = 0; i <= gx + 1; i++)
	{
        for (int j = 0; j <= gy + 1; j++) 
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
    int j = gy / 2;

    for (int i = 0; i <= gx + 1; i++)
    {
        rou = U[i][j][0];
        u = U[i][j][1] / rou;
        v = U[i][j][2] / rou;
        p = (gam - 1) * (U[i][j][3] - 0.5 * U[i][j][0] * (u * u + v * v));
        fprintf(fp, "%10lf%10lf%10lf%10lf%10lf%10lf\n", i * dx, rou, u, v, p, U[i][j][3]);
    }
    fclose(fp);
}
