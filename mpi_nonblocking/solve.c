#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "solve.h"

//Lax-Wendroff 2d sovler
void solver_2d(double ***U, double ***U_half, double ***F, double ***G, double ***temp, double dx, double dy, double dt, double gam, int s, int e, int left, int right, int rank, int size, MPI_Comm comm)
{
    exchang1(U, temp, dx, dt / 2.0, s, e, left, right, comm);
    exchang2(U, U_half, F, dx, dt / 2.0, gam, s, e, left, right, comm);
    bound(U, dx, dy, gam, s, e, rank, size);
    
    solve_y(U, U_half, G, temp, dy, dt / 2.0, gam, s, e);
    bound(U, dx, dy, gam, s, e, rank, size);

    solve_y(U, U_half, G, temp, dy, dt / 2.0, gam, s, e);
    bound(U, dx, dy, gam, s, e, rank, size);

    exchang1(U, temp, dx, dt / 2.0, s, e, left, right, comm);
    exchang2(U, U_half, F, dx, dt / 2.0, gam, s, e, left, right, comm);
    bound(U, dx, dy, gam, s, e, rank, size);
}

//update U after artificial viscosity filtering by switching function
void exchang1(double ***U, double ***temp, double dx, double dt, int s, int e, int left, int right, MPI_Comm comm)
{
    //calculate theta with U[0]
    //because theta is equal in U[0] -U[3], it is only calculated once using U[0]
    const int a = 3.0;  //speed of sound not exceeding 3
    double eta = (a * dt / dx) * (1 - a * dt / dx); 

    //allocate memory for theta[gx+2][gy +2]
    double** theta = malloc((gx + 2) * sizeof(double *));
	theta[0] = malloc((gx + 2) * (gy + 2) * sizeof(double));
	for(int i = 0; i < (gx + 2); i++)
    {
		theta[i] = theta[0] + i * (gy+2);
	}

    //non-blocking data communication
    MPI_Request reqs[4];
    MPI_Irecv(&U[0][s-1][0], (gy + 2), MPI_DOUBLE, left, 0, comm, &reqs[0]);
    MPI_Irecv(&U[0][e+1][0], (gy + 2), MPI_DOUBLE, right, 1, comm, &reqs[1]);
    MPI_Isend(&U[0][e][0], (gy + 2), MPI_DOUBLE, right, 0, comm, &reqs[2]);
    MPI_Isend(&U[0][s][0], (gy + 2), MPI_DOUBLE, left, 1, comm, &reqs[3]);

    //calculate local theta
    for (int i = s + 1; i <= e - 1; i++)
	{
        for (int j = 0; j <= gy + 1; j++) 
		{	
            theta[i][j] = fabs(fabs(U[0][i + 1][j] - U[0][i][j]) - fabs(U[0][i][j] - U[0][i - 1][j]))
            / (fabs(U[0][i + 1][j] - U[0][i][j]) + fabs(U[0][i][j] - U[0][i - 1][j]) + 1e-100);
            temp[0][i][j] = U[0][i][j] + 0.5 * eta * theta[i][j] * (U[0][i + 1][j] - 2 * U[0][i][j] + U[0][i - 1][j]); 
        }
	}

    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

    //calculate the boundary of theta
    for (int j = 0; j <= gy + 1; j++) 
	{	
        int i = s; 
        theta[i][j] = fabs(fabs(U[0][i + 1][j] - U[0][i][j]) - fabs(U[0][i][j] - U[0][i - 1][j]))
            / (fabs(U[0][i + 1][j] - U[0][i][j]) + fabs(U[0][i][j] - U[0][i - 1][j]) + 1e-100);
        temp[0][i][j] = U[0][i][j] + 0.5 * eta * theta[i][j] * (U[0][i + 1][j] - 2 * U[0][i][j] + U[0][i - 1][j]);

        i = e; 
        theta[i][j] = fabs(fabs(U[0][i + 1][j] - U[0][i][j]) - fabs(U[0][i][j] - U[0][i - 1][j]))
            / (fabs(U[0][i + 1][j] - U[0][i][j]) + fabs(U[0][i][j] - U[0][i - 1][j]) + 1e-100);
        temp[0][i][j] = U[0][i][j] + 0.5 * eta * theta[i][j] * (U[0][i + 1][j] - 2 * U[0][i][j] + U[0][i - 1][j]);  
    }

    //update U[0]
    for (int i = s; i <= e; i++)
	{
        for (int j = 0; j <= gy + 1; j++)
		{
            U[0][i][j] = temp[0][i][j];
		}
	}

    //update U[1],U[2],U[3]
    for(int k = 1; k < 4; k++)
    {
        //non-blocking data communication
        MPI_Request reqs[4];
        MPI_Irecv(&U[k][s-1][0], (gy + 2), MPI_DOUBLE, left, 0, comm, &reqs[0]);
        MPI_Irecv(&U[k][e+1][0], (gy + 2), MPI_DOUBLE, right, 1, comm, &reqs[1]);
        MPI_Isend(&U[k][e][0], (gy + 2), MPI_DOUBLE, right, 0, comm, &reqs[2]);
        MPI_Isend(&U[k][s][0], (gy + 2), MPI_DOUBLE, left, 1, comm, &reqs[3]);

        //calculate local temp
        for (int i = s + 1; i <= e - 1; i++)
	    {
            for (int j = 0; j <= gy + 1; j++)
            {
                temp[k][i][j] = U[k][i][j] + 0.5 * eta * theta[i][j] * (U[k][i + 1][j] - 2 * U[k][i][j] + U[k][i - 1][j]);
            }
        }

        MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

        //calculate the boundary of temp
        for (int j = 0; j <= gy + 1; j++)
        {
            int i = s;
            {
                temp[k][i][j] = U[k][i][j] + 0.5 * eta * theta[i][j] * (U[k][i + 1][j] - 2 * U[k][i][j] + U[k][i - 1][j]);
            }

            i = e;
            {
                temp[k][i][j] = U[k][i][j] + 0.5 * eta * theta[i][j] * (U[k][i + 1][j] - 2 * U[k][i][j] + U[k][i - 1][j]);
            }
        }

        //update U
        for (int i = s; i <= e; i++)
		{
            for (int j = 0; j <= gy + 1; j++)
			{
                U[k][i][j] = temp[k][i][j];
			}
		}
    }

    free(theta[0]);
	free(theta);
}

//update U with Lax-Wendroff method
void exchang2(double ***U, double ***U_half, double ***F, double dx, double dt, double gam, int s, int e, int left, int right, MPI_Comm comm)
{
    //calculate the local F and put it outside the loop of k to save computation
    for (int i = s ; i <= e; i++)
	{
        for (int j = 0; j <= gy + 1; j++)
        {
            double u = U[1][i][j] / U[0][i][j];
            double v = U[2][i][j] / U[0][i][j];
            double p = (gam - 1) * (U[3][i][j] - 0.5 * U[0][i][j] *(u * u + v * v));

            F[0][i][j] = U[1][i][j];
            F[1][i][j] = U[0][i][j] * u * u + p;
            F[2][i][j] = U[0][i][j] * u * v;
            F[3][i][j] = (U[3][i][j] + p) * u;
        }
	}

    //non-blocking data communication
    for(int k = 0; k < 4; k++)
    {
        MPI_Request reqs[4];
        MPI_Irecv(&U[k][s-1][0], (gy + 2), MPI_DOUBLE, left, 0, comm, &reqs[0]);
        MPI_Irecv(&U[k][e+1][0], (gy + 2), MPI_DOUBLE, right, 1, comm, &reqs[1]);
        MPI_Isend(&U[k][e][0], (gy + 2), MPI_DOUBLE, right, 0, comm, &reqs[2]);
        MPI_Isend(&U[k][s][0], (gy + 2), MPI_DOUBLE, left, 1, comm, &reqs[3]);

        //update local U_half
        for (int i = s ; i <= e-1; i++)
		{
            for (int j = 0; j <= gy + 1; j++)
			{
                U_half[k][i][j] = 0.5 * (U[k][i + 1][j] + U[k][i][j]) - 0.5 * dt / dx * (F[k][i + 1][j] - F[k][i][j]);
			}
		}

        MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);
    }
    
    //calculate the boundry of F
    for (int j = 0; j <= gy + 1; j++)
    {
        int i = s - 1;

        double u = U[1][i][j] / U[0][i][j];
        double v = U[2][i][j] / U[0][i][j];
        double p = (gam - 1) * (U[3][i][j] - 0.5 * U[0][i][j] *(u * u + v * v));

        F[0][i][j] = U[1][i][j];
        F[1][i][j] = U[0][i][j] * u * u + p;
        F[2][i][j] = U[0][i][j] * u * v;
        F[3][i][j] = (U[3][i][j] + p) * u;

        i= e + 1;
                
        u = U[1][i][j] / U[0][i][j];
        v = U[2][i][j] / U[0][i][j];
        p = (gam - 1) * (U[3][i][j] - 0.5 * U[0][i][j] *(u * u + v * v));

        F[0][i][j] = U[1][i][j];
        F[1][i][j] = U[0][i][j] * u * u + p;
        F[2][i][j] = U[0][i][j] * u * v;
        F[3][i][j] = (U[3][i][j] + p) * u;
    }

    //update the boundry of U_half
    for (int k = 0; k < 4; k++)
	{
        for (int j = 0; j <= gy + 1; j++)
		{   
            int i = s - 1;
            U_half[k][i][j] = 0.5 * (U[k][i + 1][j] + U[k][i][j]) - 0.5 * dt / dx * (F[k][i + 1][j] - F[k][i][j]);

            i = e;
            U_half[k][i][j] = 0.5 * (U[k][i + 1][j] + U[k][i][j]) - 0.5 * dt / dx * (F[k][i + 1][j] - F[k][i][j]);
		}
	}

    //calculate the F in the second difference with U_half
    for (int i = s - 1; i <= e; i++)
	{
        for (int j = 0; j <= gy + 1; j++)
		{
            double u = U_half[1][i][j] / U_half[0][i][j];
            double v = U_half[2][i][j] / U_half[0][i][j];
            double p = (gam - 1) * (U_half[3][i][j] - 0.5 * U_half[0][i][j] *(u * u + v * v));

            F[0][i][j] = U_half[1][i][j];
            F[1][i][j] = U_half[0][i][j] * u * u + p;
            F[2][i][j] = U_half[0][i][j] * u * v;
            F[3][i][j] = (U_half[3][i][j] + p) * u;
		}
	}

    //update the final U
    for (int k = 0; k < 4; k++)
    {
        for (int i = s; i <= e; i++)
		{
            for (int j = 0; j <= gy + 1; j++)
			{
                U[k][i][j] = U[k][i][j] - dt / dx * (F[k][i][j] - F[k][i - 1][j]);
			}
		}
    }
} 

//differential in y-direction
void solve_y(double ***U, double ***U_half, double ***G, double ***temp, double dy, double dt, double gam, int s, int e)
{
    const int a = 3.0;  //speed of sound not exceeding 3
    double eta = (a * dt / dy) * (1 - a * dt / dy); 
 
    //find the U after artificial viscous filtering with a switching function prior
    for (int i = s - 1; i <= e + 1; i++)
	{
        for (int j = 1; j <= gy; j++) 
		{
			//switching function
            double theta = fabs(fabs(U[0][i][j + 1] - U[0][i][j]) - fabs(U[0][i][j] - U[0][i][j - 1]))
                / (fabs(U[0][i][j + 1] - U[0][i][j]) + fabs(U[0][i][j] - U[0][i][j - 1]) + 1e-100);
			
			//sticky items
            for (int k = 0; k < 4; k++)
			{
                temp[k][i][j] = U[k][i][j] + 0.5 * eta * theta * (U[k][i][j + 1] - 2 * U[k][i][j] + U[k][i][j - 1]); 
			}
        }
	}

    for (int k = 0; k < 4; k++)
	{
        for (int i = s - 1; i <= e + 1; i++)
		{
            for (int j = 1; j <= gy; j++)
			{
                U[k][i][j] = temp[k][i][j];
			}
		}
	}
    
    //find G for the first difference in terms of U
    for (int i = 0; i <= gx + 1; i++)
	{
        for (int j = 0; j <= gy + 1; j++)
		{
            double u = U[1][i][j] / U[0][i][j];
            double v = U[2][i][j] / U[0][i][j];
            double p = (gam - 1) * (U[3][i][j] - 0.5 * U[0][i][j] *(u * u + v * v));

            G[0][i][j] = U[2][i][j];
            G[1][i][j] = U[0][i][j] * u * v;
            G[2][i][j] = U[0][i][j] * v * v + p;
            G[3][i][j] = (U[3][i][j] + p) * v;
		}
	}

    //find the U of the first differential step according to G, U_half
    for (int k = 0; k < 4; k++)
	{
        for (int i = s; i <= e + 1; i++)
		{
            for (int j = 0; j <= gy; j++)
			{
                U_half[k][i][j] = 0.5 * (U[k][i][j + 1] + U[k][i][j]) - 0.5 * dt / dy * (G[k][i][j + 1] - G[k][i][j]); 
			}
		}
	}

    //according to U_half, find the G required by the second difference
    for (int i = s; i <= e + 1; i++)
	{
        for (int j = 0; j <= gy; j++)
		{
            double u = U_half[1][i][j] / U_half[0][i][j];
            double v = U_half[2][i][j] / U_half[0][i][j];
            double p = (gam - 1) * (U_half[3][i][j] - 0.5 * U_half[0][i][j] *(u * u + v * v));

            G[0][i][j] = U_half[2][i][j];
            G[1][i][j] = U_half[0][i][j] * u * v;
            G[2][i][j] = U_half[0][i][j] * v * v + p;
            G[3][i][j] = (U_half[3][i][j] + p) * v;
		}
	}

    //a second step of differencing with G gives the final result U
    for (int k = 0; k < 4; k++)
	{
        for (int i = s; i <= e + 1; i++)
		{
            for (int j = 1; j <= gy; j++)
			{
                U[k][i][j] = U[k][i][j] - dt / dy * (G[k][i][j] - G[k][i][j - 1]); 
			}
		}
	}
}

//gather data to rank 0
void gathergrid(double ***U, int rank, int size, int s, int e, int left, int right)
{
    for(int k = 0; k < 4; k++)
    {
        for (int r = size - 1; r > 0; r--) 
        {
            if (rank == r - 1) 
            {
                MPI_Recv(&U[k][e + 1][0], (gx + 1 - e) * (gy + 2), MPI_DOUBLE, right, right, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (rank == r) 
            {
                MPI_Send(&U[k][s][0], (gx + 2 - s) * (gy + 2), MPI_DOUBLE, left, rank, MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}
