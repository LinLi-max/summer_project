#include <mpi.h>
#include "solve.h"

//Lax-Wendroff 2d sovler
void solver_2d(double U[4][gx + 2][gy + 2], double U_half[4][gx + 2][gy + 2], double F_G[4][gx + 2][gy + 2], double dx, double dy, double dt, double gam, int s, int e, int left, int right, MPI_Comm comm)
{
    exchang(U, s, e, left, right, comm);
    solve_x(U, U_half, F_G, dx, dt / 2.0, gam, s, e);
    bound(U, dx, dy, gam);

    solve_y(U, U_half, F_G, dy, dt / 2.0, gam, s, e);
    bound(U, dx, dy, gam);

    solve_y(U, U_half, F_G, dy, dt / 2.0, gam, s, e);
    bound(U, dx, dy, gam);

    exchang(U, s, e, left, right, comm);
    solve_x(U, U_half, F_G, dx, dt / 2.0, gam, s, e);
    bound(U, dx, dy, gam);
}

//differential in x-direction
void solve_x(double U[4][gx + 2][gy + 2], double U_half[4][gx + 2][gy + 2], double F[4][gx + 2][gy + 2], double dx, double dt, double gam, int s, int e)
{
    const int a = 3.0;  //speed of sound not exceeding 3
    double eta = (a * dt / dx) * (1 - a * dt / dx); 
    double temp[4][gx + 2][gy + 2];

    //find the U after artificial viscous filtering with a switching function prior
    for (int i = s; i <= e; i++)
	{
        for (int j = 0; j <= gy + 1; j++) 
		{	
			//switching function
            double theta = fabs(fabs(U[0][i + 1][j] - U[0][i][j]) - fabs(U[0][i][j] - U[0][i - 1][j]))
                / (fabs(U[0][i + 1][j] - U[0][i][j]) + fabs(U[0][i][j] - U[0][i - 1][j]) + 1e-100);

			//sticky items
            for (int k = 0; k < 4; k++)
			{
                temp[k][i][j] = U[k][i][j] + 0.5 * eta * theta * (U[k][i + 1][j] - 2 * U[k][i][j] + U[k][i - 1][j]); 
			}
        }
	}

    for (int j = 0; j <= gy + 1; j++)
	{   
        for (int i = s; i <= e; i++)
		{
            for (int k = 0; k < 4; k++)
			{
                U[k][i][j] = temp[k][i][j];
			}
		}
	}

    //find F for the first difference in terms of U
    for (int i = s - 1 ; i <= e + 1; i++)
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
            //for (int k = 0; k <= 4; k++)
		    //{
                //U_F(U[][i][j], F[][i][j], gam);
		    //}
        }
        
	}

    //find the U of the first differential step according to F, U_half
    for (int i = s - 1; i <= e; i++)
	{
        for (int j = 0; j <= gy + 1; j++)
		{
            for (int k = 0; k < 4; k++)
			{
                U_half[k][i][j] = 0.5 * (U[k][i + 1][j] + U[k][i][j]) - 0.5 * dt / dx * (F[k][i + 1][j] - F[k][i][j]);
			}
		}
	}

    //according to U_half, find the F required by the second difference
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
            //U_F(U_half[i][j], F[i][j], gam); 
		}
	}

    //a second step of differencing with F gives the final result U
    for (int i = s; i <= e; i++)
	{
        for (int j = 0; j <= gy + 1; j++)
		{
            for (int k = 0; k < 4; k++)
			{
                U[k][i][j] = U[k][i][j] - dt / dx * (F[k][i][j] - F[k][i - 1][j]);
			}
		}
	}
}

//differential in y-direction
void solve_y(double U[4][gx + 2][gy + 2], double U_half[4][gx + 2][gy + 2], double G[4][gx + 2][gy + 2], double dy, double dt, double gam, int s, int e)
{
    const int a = 3.0;  //speed of sound not exceeding 3
    double eta = (a * dt / dy) * (1 - a * dt / dy); 
    double temp[4][gx + 2][gy + 2];
 
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

    for (int i = s - 1; i <= e + 1; i++)
	{
        for (int j = 1; j <= gy; j++)
		{
            for (int k = 0; k < 4; k++)
			{
                U[k][i][j] = temp[k][i][j];
			}
		}
	}
    
    //find G for the first difference in terms of U
    for (int i = s; i <= e + 1; i++)
	{
        for (int j = 0; j <= gy + 1; j++)
		{
            //U_G(U[i][j], G[i][j], gam);
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
    for (int i = s; i <= e + 1; i++)
	{
        for (int j = 0; j <= gy; j++)
		{
            for (int k = 0; k < 4; k++)
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
            //U_G(U_half[i][j], G[i][j], gam); 
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
    for (int i = s; i <= e + 1; i++)
	{
        for (int j = 1; j <= gy; j++)
		{
            for (int k = 0; k < 4; k++)
			{
                U[k][i][j] = U[k][i][j] - dt / dy * (G[k][i][j] - G[k][i][j - 1]); 
			}
		}
	}
}

void exchang(double U[4][gx + 2][gy + 2], int s, int e, int left, int right, MPI_Comm comm)
{   
    for(int i = 0; i < 4; i++)
    {
        MPI_Send(&U[i][e][0], (gy + 2), MPI_DOUBLE, right, 0, comm);
        MPI_Recv(&U[i][s-1][0], (gy + 2), MPI_DOUBLE, left, 0, comm, MPI_STATUS_IGNORE);

        MPI_Send(&U[i][s][0], (gy + 2), MPI_DOUBLE, left, 1, comm);
        MPI_Recv(&U[i][e+1][0], (gy + 2), MPI_DOUBLE, right, 1, comm, MPI_STATUS_IGNORE);
    }
}

void gathergrid(double U[4][gx + 2][gy + 2], int rank, int size, int s, int e, int left, int right)
{
    for(int m = 0; m < 4; m++)
    {
        for (int i = size - 1; i > 0; i--) 
        {
            if (rank == i - 1) 
            {
                //printf("%d recieve size = %d\n", rank, (gx + 1 - e));
                MPI_Recv(&U[m][e + 1][0], (gx + 1 - e) * (gy + 2), MPI_DOUBLE, right, right, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (rank == i) 
            {
                //printf("%d send size = %d\n", rank, (gx + 2 - s));
                MPI_Send(&U[m][s][0], (gx + 2 - s) * (gy + 2), MPI_DOUBLE, left, rank, MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}
