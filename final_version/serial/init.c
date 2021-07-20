#include "init.h"

//initialization
void init(double ***U, double dx, double dy, double gam) 
{
    const double rou1 = 1.0, u1 = 2.9, v1 = 0.0, p1 = 0.71429;  //uniform inlet condition at left boundary
    const double rou2 = 1.69997, u2 = 2.61934, v2 = -0.50632, p2 = 1.52819;  //uniform inlet condition at up boundary
    const double pi = 3.141592654, alpha = 29 * pi / 180;  //incidence radian of oblique shock wave

    for (int i = 0; i <= gx + 1; i++)
	{
        for (int j = 0; j <= gy + 1; j++) 
		{
            double x = (1 - j * dy) / tan(alpha);
            if (i * dx <= x) 
			{
                U[0][i][j] = rou1; 
                U[1][i][j] = rou1 * u1; 
                U[2][i][j] = rou1 * v1;
                U[3][i][j] = p1 / (gam - 1) + rou1 * (u1 * u1 + v1 * v1) / 2; 
            } 
			else
			{
                U[0][i][j] = rou2;
                U[1][i][j] = rou2 * u2;
                U[2][i][j] = rou2 * v2;
                U[3][i][j] = p2 / (gam - 1) + rou2 * (u2 * u2 + v2 * v2) / 2;
            }
        }
	}
}

//deal with the boundary
void bound(double ***U, double dx, double dy, double gam)
{
    const double rou1 = 1.0, u1 = 2.9, v1 = 0.0, p1 = 0.71429;  //uniform inlet condition at left boundary
    const double rou2 = 1.69997, u2 = 2.61934, v2 = -0.50632, p2 = 1.52819;  //uniform inlet condition at up boundary

    //left
    for (int j = 0; j <= gy + 1; j++) 
	{
        U[0][0][j] = rou1;
        U[1][0][j] = rou1 * u1;
        U[2][0][j] = rou1 * v1;
        U[3][0][j] = p1 / (gam - 1) + rou1 * (u1 * u1 + v1 * v1) / 2;
    }

    //right, free boundary
    //the last column of parameters equal to the penultimate column
    for (int k = 0; k < 4; k++) 
	{
        for (int j = 0; j <= gy + 1; j++)
		{
            U[k][gx + 1][j] = U[k][gx][j]; 
        }
	}

    //up
    for (int i = 0; i <= gx + 1; i++) 
	{
        U[0][i][gy + 1] = rou2;
        U[1][i][gy + 1] = rou2 * u2;
        U[2][i][gy + 1] = rou2 * v2;
        U[3][i][gy + 1] = p2 / (gam - 1) + rou2 * (u2 * u2 + v2 * v2) / 2;
    }

    //down, horizontal rigid wall
    //the y-directional velocity perpendicular to the wall is 0
    //the other terms have values equal to the previous one
    for (int i = 0; i <= gx + 1; i++) 
	{
        U[0][i][0] = U[0][i][1];
        U[1][i][0] = U[1][i][1];
        U[2][i][0] = 0;
        U[3][i][0] = U[3][i][1];
    }
}

//calculation of dt based on the algorithm stability parameter cfl
void cfl(double ***U, double dx, double dy, double* dt, double sf, double gam)
{
	double maxvel = 0;
    
    for (int i = 1; i <= gx; i++)
	{
		for (int j = 1; j <= gy; j++) 
		{
            double u = U[1][i][j] / U[0][i][j];
            double v = U[2][i][j] / U[0][i][j];
            double p = (gam - 1) * (U[3][i][j] - 0.5 * U[0][i][j] * (u * u + v * v));
            double velocity = sqrt(gam * p / U[0][i][j]) + sqrt(u * u + v * v);  //feature speed
            if (velocity > maxvel)
			{
				maxvel = velocity;
			}    
        }
	}

    *dt = sf * min(dx, dy) / maxvel; 
}
