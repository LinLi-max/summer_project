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

//differential in x-direction
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
        for (int k = 0; k < 4; k++)
		{
            for (int i = 1; i <= gx; i++)
			{
                U[i][j][k] = temp[i][j][k];
			}
		}
	}
    //memcpy(U,F,(gy + 2)*gx*4*sizeof(double)); 

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
}

//differential in y-direction
void solve_y(double U[gx + 2][gy + 2][4], double U_half[gx + 2][gy + 2][4], double G[gx + 2][gy + 2][4], double dy, double dt, double gam)
{
    const int a = 3.0;  //speed of sound not exceeding 3
    double eta = (a * dt / dy) * (1 - a * dt / dy); //C.5b
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
        for (int k = 0; k < 4; k++)
		{
            for (int j = 1; j <= gy; j++)
			{
                U[i][j][k] = temp[i][j][k];
			}
		}
	}
    //memcpy(U,G,(gx + 2)*gy*4*sizeof(double)); 
    
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

//Lax-Wendroff 2d sovler
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
}
