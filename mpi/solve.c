#include <mpi.h>
#include "solve.h"

//calculate the f term in the euler equation
void U_F(double U[4], double F[4], double gam)
{
    double u = U[1] / U[0];
    double v = U[2] / U[0];
    double p = (gam - 1) * (U[3] - 0.5 * U[0] *(u * u + v * v));

    F[0] = U[1];
    F[1] = U[0] * u * u + p;
    F[2] = U[0] * u * v;
    F[3] = (U[3] + p) * u;
}

//calculate the g term in the euler equation
void U_G(double U[4], double G[4], double gam) 
{
    double u = U[1] / U[0];
    double v = U[2] / U[0];
    double p = (gam - 1) * (U[3] - 0.5 * U[0] *(u * u + v * v));

    G[0] = U[2];
    G[1] = U[0] * u * v;
    G[2] = U[0] * v * v + p;
    G[3] = (U[3] + p) * v;
}

/* //differential in x-direction
void solve_x(double U[gx + 2][gy + 2][4], double U_half[gx + 2][gy + 2][4], double F[gx + 2][gy + 2][4], double dx, double dt, double gam)
{
    const int a = 3.0;  //speed of sound not exceeding 3
    double eta = (a * dt / dx) * (1 - a * dt / dx); 
    double temp[gx + 2][gy + 2][4];

    //find the U after artificial viscous filtering with a switching function prior
    for (int i = 1; i <= gx; i++)
	{
        for (int j = 0; j <= gy + 1; j++) 
		{	
			//switching function
            double theta = fabs(fabs(U[i + 1][j][0] - U[i][j][0]) - fabs(U[i][j][0] - U[i - 1][j][0]))
                / (fabs(U[i + 1][j][0] - U[i][j][0]) + fabs(U[i][j][0] - U[i - 1][j][0]) + 1e-100);

			//sticky items
            for (int k = 0; k < 4; k++)
			{
                temp[i][j][k] = U[i][j][k] + 0.5 * eta * theta * (U[i + 1][j][k] - 2 * U[i][j][k] + U[i - 1][j][k]); 
			}
        }
	}

    for (int j = 0; j <= gy + 1; j++)
	{   
        for (int i = 1; i <= gx; i++)
		{
            for (int k = 0; k < 4; k++)
			{
                U[i][j][k] = temp[i][j][k];
			}
		}
	}

    //find F for the first difference in terms of U
    for (int i = 0; i <= gx + 1; i++)
	{
        for (int j = 0; j <= gy + 1; j++)
		{
            U_F(U[i][j], F[i][j], gam);
		}
	}

    //find the U of the first differential step according to F, U_half
    for (int i = 0; i <= gx; i++)
	{
        for (int j = 0; j <= gy + 1; j++)
		{
            for (int k = 0; k < 4; k++)
			{
                U_half[i][j][k] = 0.5 * (U[i + 1][j][k] + U[i][j][k]) - 0.5 * dt / dx * (F[i + 1][j][k] - F[i][j][k]);
			}
		}
	}

    //according to U_half, find the F required by the second difference
    for (int i = 0; i <= gx; i++)
	{
        for (int j = 0; j <= gy + 1; j++)
		{
            U_F(U_half[i][j], F[i][j], gam); 
		}
	}

    //a second step of differencing with F gives the final result U
    for (int i = 1; i <= gx; i++)
	{
        for (int j = 0; j <= gy + 1; j++)
		{
            for (int k = 0; k < 4; k++)
			{
                U[i][j][k] = U[i][j][k] - dt / dx * (F[i][j][k] - F[i - 1][j][k]);
			}
		}
	}
} */

/* //differential in y-direction
void solve_y(double U[gx + 2][gy + 2][4], double U_half[gx + 2][gy + 2][4], double G[gx + 2][gy + 2][4], double dy, double dt, double gam)
{
    const int a = 3.0;  //speed of sound not exceeding 3
    double eta = (a * dt / dy) * (1 - a * dt / dy); 
    double temp[gx + 2][gy + 2][4];
 
    //find the U after artificial viscous filtering with a switching function prior
    for (int i = 0; i <= gx + 1; i++)
	{
        for (int j = 1; j <= gy; j++) 
		{
			//switching function
            double theta = fabs(fabs(U[i][j + 1][0] - U[i][j][0]) - fabs(U[i][j][0] - U[i][j - 1][0]))
                / (fabs(U[i][j + 1][0] - U[i][j][0]) + fabs(U[i][j][0] - U[i][j - 1][0]) + 1e-100);
			
			//sticky items
            for (int k = 0; k < 4; k++)
			{
                temp[i][j][k] = U[i][j][k] + 0.5 * eta * theta * (U[i][j + 1][k] - 2 * U[i][j][k] + U[i][j - 1][k]); 
			}
        }
	}

    for (int i = 0; i <= gx + 1; i++)
	{
        for (int j = 1; j <= gy; j++)
		{
            for (int k = 0; k < 4; k++)
			{
                U[i][j][k] = temp[i][j][k];
			}
		}
	}
    
    //find G for the first difference in terms of U
    for (int i = 0; i <= gx + 1; i++)
	{
        for (int j = 0; j <= gy + 1; j++)
		{
            U_G(U[i][j], G[i][j], gam);
		}
	}

    //find the U of the first differential step according to G, U_half
    for (int i = 0; i <= gx + 1; i++)
	{
        for (int j = 0; j <= gy; j++)
		{
            for (int k = 0; k < 4; k++)
			{
                U_half[i][j][k] = 0.5 * (U[i][j + 1][k] + U[i][j][k]) - 0.5 * dt / dy * (G[i][j + 1][k] - G[i][j][k]); 
			}
		}
	}

    //according to U_half, find the G required by the second difference
    for (int i = 0; i <= gx + 1; i++)
	{
        for (int j = 0; j <= gy; j++)
		{
            U_G(U_half[i][j], G[i][j], gam); 
		}
	}

    //a second step of differencing with G gives the final result U
    for (int i = 0; i <= gx + 1; i++)
	{
        for (int j = 1; j <= gy; j++)
		{
            for (int k = 0; k < 4; k++)
			{
                U[i][j][k] = U[i][j][k] - dt / dy * (G[i][j][k] - G[i][j - 1][k]); 
			}
		}
	}
} */

/* //Lax-Wendroff 2d sovler
void solver_2d(double U[gx + 2][gy + 2][4], double U_half[gx + 2][gy + 2][4], double F_G[gx + 2][gy + 2][4], double dx, double dy, double dt, double gam)
{
    solve_x(U, U_half, F_G, dx, dt / 2.0, gam);
    bound(U, dx, dy, gam);

    solve_y(U, U_half, F_G, dy, dt / 2.0, gam);
    bound(U, dx, dy, gam);

    solve_y(U, U_half, F_G, dy, dt / 2.0, gam);
    bound(U, dx, dy, gam);

    solve_x(U, U_half, F_G, dx, dt / 2.0, gam);
    bound(U, dx, dy, gam);
} */

//Lax-Wendroff 2d sovler
void solver_2d(double U[gx + 2][gy + 2][4], double U_half[gx + 2][gy + 2][4], double F_G[gx + 2][gy + 2][4], double dx, double dy, double dt, double gam, int s, int e, int left, int right, MPI_Comm comm)
{
    exchang(U, s, e, left, right, comm);
    solve_x(U, U_half, F_G, dx, dt / 2.0, gam, s, e);
    //solve_x(U, U_half, F_G, dx, dt / 2.0, gam);
    bound(U, dx, dy, gam);

    solve_y(U, U_half, F_G, dy, dt / 2.0, gam, s, e);
    bound(U, dx, dy, gam);

    solve_y(U, U_half, F_G, dy, dt / 2.0, gam, s, e);
    bound(U, dx, dy, gam);

    exchang(U, s, e, left, right, comm);
    solve_x(U, U_half, F_G, dx, dt / 2.0, gam, s, e);
    //solve_x(U, U_half, F_G, dx, dt / 2.0, gam);
    bound(U, dx, dy, gam);
}

//differential in x-direction
void solve_x(double U[gx + 2][gy + 2][4], double U_half[gx + 2][gy + 2][4], double F[gx + 2][gy + 2][4], double dx, double dt, double gam, int s, int e)
{
    const int a = 3.0;  //speed of sound not exceeding 3
    double eta = (a * dt / dx) * (1 - a * dt / dx); 
    double temp[gx + 2][gy + 2][4];

    //find the U after artificial viscous filtering with a switching function prior
    for (int i = s; i <= e; i++)
	{
        for (int j = 0; j <= gy + 1; j++) 
		{	
			//switching function
            double theta = fabs(fabs(U[i + 1][j][0] - U[i][j][0]) - fabs(U[i][j][0] - U[i - 1][j][0]))
                / (fabs(U[i + 1][j][0] - U[i][j][0]) + fabs(U[i][j][0] - U[i - 1][j][0]) + 1e-100);

			//sticky items
            for (int k = 0; k < 4; k++)
			{
                temp[i][j][k] = U[i][j][k] + 0.5 * eta * theta * (U[i + 1][j][k] - 2 * U[i][j][k] + U[i - 1][j][k]); 
			}
        }
	}

    for (int j = 0; j <= gy + 1; j++)
	{   
        for (int i = s; i <= e; i++)
		{
            for (int k = 0; k < 4; k++)
			{
                U[i][j][k] = temp[i][j][k];
			}
		}
	}

    //find F for the first difference in terms of U
    for (int i = s - 1 ; i <= e + 1; i++)
	{
        for (int j = 0; j <= gy + 1; j++)
		{
            U_F(U[i][j], F[i][j], gam);
		}
	}

    //find the U of the first differential step according to F, U_half
    for (int i = s - 1; i <= e; i++)
	{
        for (int j = 0; j <= gy + 1; j++)
		{
            for (int k = 0; k < 4; k++)
			{
                U_half[i][j][k] = 0.5 * (U[i + 1][j][k] + U[i][j][k]) - 0.5 * dt / dx * (F[i + 1][j][k] - F[i][j][k]);
			}
		}
	}

    //according to U_half, find the F required by the second difference
    for (int i = s - 1; i <= e; i++)
	{
        for (int j = 0; j <= gy + 1; j++)
		{
            U_F(U_half[i][j], F[i][j], gam); 
		}
	}

    //a second step of differencing with F gives the final result U
    for (int i = s; i <= e; i++)
	{
        for (int j = 0; j <= gy + 1; j++)
		{
            for (int k = 0; k < 4; k++)
			{
                U[i][j][k] = U[i][j][k] - dt / dx * (F[i][j][k] - F[i - 1][j][k]);
			}
		}
	}
}

//differential in y-direction
void solve_y(double U[gx + 2][gy + 2][4], double U_half[gx + 2][gy + 2][4], double G[gx + 2][gy + 2][4], double dy, double dt, double gam, int s, int e)
{
    const int a = 3.0;  //speed of sound not exceeding 3
    double eta = (a * dt / dy) * (1 - a * dt / dy); 
    double temp[gx + 2][gy + 2][4];
 
    //find the U after artificial viscous filtering with a switching function prior
    for (int i = 0; i <= gx + 1; i++)
	{
        for (int j = 1; j <= gy; j++) 
		{
			//switching function
            double theta = fabs(fabs(U[i][j + 1][0] - U[i][j][0]) - fabs(U[i][j][0] - U[i][j - 1][0]))
                / (fabs(U[i][j + 1][0] - U[i][j][0]) + fabs(U[i][j][0] - U[i][j - 1][0]) + 1e-100);
			
			//sticky items
            for (int k = 0; k < 4; k++)
			{
                temp[i][j][k] = U[i][j][k] + 0.5 * eta * theta * (U[i][j + 1][k] - 2 * U[i][j][k] + U[i][j - 1][k]); 
			}
        }
	}

    for (int i = 0; i <= gx + 1; i++)
	{
        for (int j = 1; j <= gy; j++)
		{
            for (int k = 0; k < 4; k++)
			{
                U[i][j][k] = temp[i][j][k];
			}
		}
	}
    
    //find G for the first difference in terms of U
    for (int i = 0; i <= gx + 1; i++)
	{
        for (int j = 0; j <= gy + 1; j++)
		{
            U_G(U[i][j], G[i][j], gam);
		}
	}

    //find the U of the first differential step according to G, U_half
    for (int i = 0; i <= gx + 1; i++)
	{
        for (int j = 0; j <= gy; j++)
		{
            for (int k = 0; k < 4; k++)
			{
                U_half[i][j][k] = 0.5 * (U[i][j + 1][k] + U[i][j][k]) - 0.5 * dt / dy * (G[i][j + 1][k] - G[i][j][k]); 
			}
		}
	}

    //according to U_half, find the G required by the second difference
    for (int i = 0; i <= gx + 1; i++)
	{
        for (int j = 0; j <= gy; j++)
		{
            U_G(U_half[i][j], G[i][j], gam); 
		}
	}

    //a second step of differencing with G gives the final result U
    for (int i = 0; i <= gx + 1; i++)
	{
        for (int j = 1; j <= gy; j++)
		{
            for (int k = 0; k < 4; k++)
			{
                U[i][j][k] = U[i][j][k] - dt / dy * (G[i][j][k] - G[i][j - 1][k]); 
			}
		}
	}
}

void exchang(double U[gx + 2][gy + 2][4], int s, int e, int left, int right, MPI_Comm comm)
{   
    MPI_Datatype ctype;
    MPI_Type_vector((gy+2), 1, (gx+2), MPI_DOUBLE, &ctype);
    MPI_Type_commit(&ctype);

    for(int i = 0; i < 4; i++)
    {
        MPI_Send(&U[e][0][i], 1, ctype, right, 0, comm);
        MPI_Recv(&U[s-1][0][i], 1, ctype, left, 0, comm, MPI_STATUS_IGNORE);

        MPI_Send(&U[s][0][i], 1, ctype, left, 1, comm);
        MPI_Recv(&U[e+1][0][i], 1, ctype, right, 1, comm, MPI_STATUS_IGNORE);
    }
    
    MPI_Type_free(&ctype);
}

void gathergrid(double U[gx + 2][gy + 2][4], int rank, int size, int s, int e, int left, int right)
{
    for(int m = 0; m < 4; m++)
    {
        for (int i = size - 1; i > 0; i--) 
        {
            if (rank == i - 1) 
            {
                //printf("%d recieve size = %d\n", rank, (gx + 1 - e));
                MPI_Recv(&U[e + 1][0][m], (gx + 1 - e) * (gy + 2), MPI_DOUBLE, right, right, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            if (rank == i) 
            {
                //printf("%d send size = %d\n", rank, (gx + 2 - s));
                MPI_Send(&U[s][0][m], (gx + 2 - s) * (gy + 2), MPI_DOUBLE, left, rank, MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}
