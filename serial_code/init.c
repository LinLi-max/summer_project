#include "init.h"

//calculation of dt based on the algorithm stability parameter cfl
void cfl(double U[gx + 2][gy + 2][4], double dx, double dy, double* dt, double sf, double gam)
{
	double maxvel = 0;
    
    for (int i = 1; i <= gx; i++)
	{
		for (int j = 1; j <= gy; j++) 
		{
            double u = U[i][j][1] / U[i][j][0];
            double v = U[i][j][2] / U[i][j][0];
            double p = (gam - 1) * (U[i][j][3] - 0.5 * U[i][j][0] * (u * u + v * v));
            double velocity = sqrt(gam * p / U[i][j][0]) + sqrt(u * u + v * v);  //feature speed
            if (velocity > maxvel)
			{
				maxvel = velocity;
			}    
        }
	}

    *dt = sf * min(dx, dy) / maxvel; 
}

//initialization
void init(double U[gx + 2][gy + 2][4], double dx, double dy, double gam) 
{
    const double rou1 = 1.0, u1 = 2.9, v1 = 0.0, p1 = 0.71429; 
    const double rou2 = 1.69997, u2 = 2.61934, v2 = -0.50632, p2 = 1.52819;
    const double pi = 3.141592654, alpha = 29 * pi / 180;  //incidence radian of oblique shock wave

    for (int j = 0; j <= gy + 1; j++)
	{
        for (int i = 0; i <= gx + 1; i++) 
		{
            double x = (1 - j * dy) / tan(alpha);
            if (i * dx <= x) 
			{
                U[i][j][0] = rou1; 
                U[i][j][1] = rou1 * u1; 
                U[i][j][2] = rou1 * v1;
                U[i][j][3] = p1 / (gam - 1) + rou1 * (u1 * u1 + v1 * v1) / 2; 
            } 
			else
			{
                U[i][j][0] = rou2;
                U[i][j][1] = rou2 * u2;
                U[i][j][2] = rou2 * v2;
                U[i][j][3] = p2 / (gam - 1) + rou2 * (u2 * u2 + v2 * v2) / 2;
            }
        }
	}
}

//deal with the boundary
void bound(double U[gx + 2][gy + 2][4], double dx, double dy, double gam)
{
    const double rou1 = 1.0, u1 = 2.9, v1 = 0.0, p1 = 0.71429; 
    const double rou2 = 1.69997, u2 = 2.61934, v2 = -0.50632, p2 = 1.52819;

    //left
    for (int j = 0; j <= gy + 1; j++) 
	{
        U[0][j][0] = rou1;
        U[0][j][1] = rou1 * u1;
        U[0][j][2] = rou1 * v1;
        U[0][j][3] = p1 / (gam - 1) + rou1 * (u1 * u1 + v1 * v1) / 2;
    }

    //right, free boundary
    //the last column of parameters equal to the penultimate column
    for (int j = 0; j <= gy + 1; j++)
	{
        for (int k = 0; k < 4; k++) 
		{
            U[gx + 1][j][k] = U[gx][j][k]; 
        }
	}

    //up
    for (int i = 0; i <= gx + 1; i++) 
	{
        U[i][gy + 1][0] = rou2;
        U[i][gy + 1][1] = rou2 * u2;
        U[i][gy + 1][2] = rou2 * v2;
        U[i][gy + 1][3] = p2 / (gam - 1) + rou2 * (u2 * u2 + v2 * v2) / 2;
    }

    //down, horizontal rigid wall
    //the y-directional velocity perpendicular to the wall is 0
    //the other terms have values equal to the previous one
    for (int i = 0; i <= gx + 1; i++) 
	{
        U[i][0][0] = U[i][1][0];
        U[i][0][1] = U[i][1][1];
        U[i][0][2] = 0;
        U[i][0][3] = U[i][1][3];
    }
}
