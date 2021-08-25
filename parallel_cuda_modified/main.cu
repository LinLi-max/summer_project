/**
 * \file main.c
 * \brief The main function of oblique shock reflection simulation(cuda version).
 * \author Lin li
 * \version 2.0
 * \date 2021-07-29
 */

#include <stdio.h>
#include <time.h>
#include "kernel.h"

int output(double* d_U, double dx, double dy, double gam, int col, int row);

int main(int argc, char* argv[])
{
    //calculation area
    const double x = 4.0, y = 1.0;

    //physical constants
    const double sf = 0.8, gam = 1.4; 
    
    //time of simulation
    const double total_time = 3.0;

    clock_t start, finish;
    double runtime;
    double current_time = 0;
    const double dx = x / gx;
    const double dy = y / gy;
    int TPB = 128;
    int col = gx + 2;
    int row = gy + 2;
    double* dt, * d_temp;
    double* d_U, * d_U_half, * d_FG;
    double* d_U_cpu, *dt_cpu;

    if (argc != 2)
    {
        printf("Please input TPB\n");
        return -1;
    }

    TPB = atoi(argv[1]);
    printf("TPB = %d\n", TPB);
    printf("gx = %d, gy = %d\n", gx, gy);
    printf("Simulation time = %.3lf seconds.\n", total_time);

    //cudaMallocManaged((void**)&d_U, sizeof(double) * col * row * 4);
    //cudaMallocManaged((void**)&d_U_half, sizeof(double) * col * row * 4);
    //cudaMallocManaged((void**)&d_FG, sizeof(double) * col * row * 4);
    //cudaMallocManaged((void**)&d_temp, sizeof(double) * col * row * 4);
    //cudaMallocManaged((void**)&dt, sizeof(double));

    cudaMalloc((void**)&d_U, sizeof(double) * col * row * 4);
	d_U_cpu = (double*)malloc(sizeof(double) * col * row * 4);

    cudaMalloc((void**)&d_U_half, sizeof(double) * col * row * 4);
    cudaMalloc((void**)&d_FG, sizeof(double) * col * row * 4);
    cudaMalloc((void**)&d_temp, sizeof(double) * col * row * 4);

    cudaMalloc((void**)&dt, sizeof(double));
	dt_cpu = (double*)malloc(sizeof(double));

    start = clock();

    call_init(d_U, dx, dy, gam, col, row, TPB);  //initialization on the GPU

	// Need to bring d_U and dt back to the host before calling call_cfl
	cudaMemcpy(d_U_cpu, d_U, sizeof(double) * col * row * 4, cudaMemcpyDeviceToHost);
	cudaMemcpy(dt_cpu, dt, sizeof(double), cudaMemcpyDeviceToHost);

    call_cfl(d_U_cpu, dx, dy, dt_cpu, sf, gam, col, row);  //calculate dt using the cfl number - on the CPU

	// Need to bring dt to the GPU before calling the main loop
	cudaMemcpy(dt, dt_cpu, sizeof(double), cudaMemcpyHostToDevice);

    while (current_time < total_time)
    {
        call_solver_2d(d_U, d_U_half, d_FG, d_temp, dx, dy, dt_cpu[0], gam, col, row, TPB);  //use a second order Lax-Wendroff two-step difference scheme
        current_time = current_time + *dt_cpu;
        //printf("time = %10g\n", current_time);  //the simulation process should be commented out when testing the runing time  
    }

    finish = clock();

    runtime = (finish - start) / CLOCKS_PER_SEC;
    printf("Program done, the cuda version use %.3lf seconds.\n", runtime);

	// Need to bring d_U back to the host before calling output
	cudaMemcpy(d_U_cpu, d_U, sizeof(double) * col * row * 4, cudaMemcpyDeviceToHost);
    output(d_U_cpu, dx, dy, gam, col, row);

    cudaFree(d_U);
    cudaFree(d_U_half);
    cudaFree(d_FG);
    cudaFree(d_temp);
    cudaFree(dt);

	free(d_U_cpu);
	free(dt_cpu);

    return 0;
}

int output(double* d_U, double dx, double dy, double gam, int col, int row)
{
    FILE* fp;
    double rou, u, v, p;

    //for density contour map
    fp = fopen("result.txt", "w");
    if (fp == NULL)
    {
        perror("Error opening result.txt for writing");
        return(-1);
    }

    fprintf(fp, "variables = x, y, rou, u, v, p, E\n");
    fprintf(fp, "gx = %d, gy = %d\n", gx, gy);

    long int size = col * row;
    for (int i = 1; i <= gx; i++)
    {
        for (int j = 1; j <= gy; j++)
        {
            rou = d_U[i + col * j];
            u = d_U[i + col * j + size] / rou;
            v = d_U[i + col * j + size * 2] / rou;
            p = (gam - 1) * (d_U[i + col * j + size * 3] - 0.5 * d_U[i + col * j] * (u * u + v * v));
            fprintf(fp, "%10lf%10lf%10lf%10lf%10lf%10lf%10lf\n", i * dx, j * dy, rou, u, v, p, d_U[i + col * j + size * 3]);
        }
    }
    fclose(fp);

    //for pressure map at position y / 2
    fp = fopen("result_0.5y.txt", "w");
    if (fp == NULL)
    {
        perror("Error opening result_0.5y.txt for writing");
        return(-1);
    }

    fprintf(fp, "variables = x, rou, u, v, p, E\n");
    fprintf(fp, "gx = %d, gy = %d, 0.5y position\n", gx, gy);

    int j = gy / 2;
    for (int i = 1; i <= gx; i++)
    {
        rou = d_U[i + col * j];
        u = d_U[i + col * j + size] / rou;
        v = d_U[i + col * j + size * 2] / rou;
        p = (gam - 1) * (d_U[i + col * j + size * 3] - 0.5 * d_U[i + col * j] * (u * u + v * v));
        fprintf(fp, "%10lf%10lf%10lf%10lf%10lf%10lf\n", i * dx, rou, u, v, p, d_U[i + col * j + size * 3]);
    }
    fclose(fp);

    return 0;
}
