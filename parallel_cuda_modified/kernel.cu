#include "kernel.h"

//initialization
__global__ void init_kernel(double* d_U, double dx, double dy, double gam, int col, int row)
{
    unsigned long long tid = (blockIdx.x * blockDim.x + threadIdx.x);

    const double rou1 = 1.0, u1 = 2.9, v1 = 0.0, p1 = 0.71429;  //uniform inlet condition at left boundary
    const double rou2 = 1.69997, u2 = 2.61934, v2 = -0.50632, p2 = 1.52819;  //uniform inlet condition at up boundary
    const double pi = 3.141592654, alpha = 29 * pi / 180;  //incidence radian of oblique shock wave

    if (tid >= col * row)
    {
        return;
    }
        
    int size = col * row;
    int  xx = tid % col;
    int  yy = tid / col;

    double x = (1 - yy * dy) / tan(alpha);

    if (xx * dx <= x)
    {
        d_U[tid] = rou1;
        d_U[tid + size] = rou1 * u1;
        d_U[tid + size * 2] = rou1 * v1;
        d_U[tid + size * 3] = p1 / (gam - 1) + rou1 * (u1 * u1 + v1 * v1) / 2;
    }
    else
    {
        d_U[tid] = rou2;
        d_U[tid + size] = rou2 * u2;
        d_U[tid + size * 2] = rou2 * v2;
        d_U[tid + size * 3] = p2 / (gam - 1) + rou2 * (u2 * u2 + v2 * v2) / 2;
    }

}

void call_init(double* d_U, double dx, double dy, double gam, int col, int row, int TPB)
{
    int BSIZE = (col * row + (TPB - 1)) / TPB;
    init_kernel << <BSIZE, TPB >> > (d_U, dx, dy, gam, col, row);
    cudaDeviceSynchronize();
}

//calculation of dt based on the algorithm stability parameter cfl
__global__ void cfl_kernel(double* d_U, double dx, double dy, double* dt, double gam, double sf, int col, int row)
{
    unsigned long long tid = (blockIdx.x * blockDim.x + threadIdx.x);
    double maxvel = 0;

    int xx = tid % col;
    int yy = tid / col;
    int size = col * row;

    if (xx < 1 || xx >= col || yy < 1 || yy >= row)
    {
        return;
    }
        
    double u = d_U[tid + size] / d_U[tid];
    double v = d_U[tid + size * 2] / d_U[tid];
    double p = (gam - 1) * (d_U[tid + size * 3] - 0.5 * d_U[tid] * (u * u + v * v));
    double velocity = sqrt(gam * p / d_U[tid]) + sqrt(u * u + v * v);  //feature speed
    if (velocity > maxvel)
    {
        maxvel = velocity;
    }

    dt[0] = sf * min(dx, dy) / maxvel;
}

void call_cfl(double* d_U, double dx, double dy, double* dt, double sf, double gam, int col, int row)
{

    double maxvel = 0;
    long int size = col * row;

    for (int i = 1; i <= gx; i++)
    {
        for (int j = 1; j <= gy; j++)
        {
            double u0 = d_U[i + j * col];
            double u1 = d_U[i + j * col + size];
            double u2 = d_U[i + j * col + size * 2];
            double u3 = d_U[i + j * col + size * 3];

            double u = u1 / u0;
            double v = u2 / u0;
            double p = (gam - 1) * (u3 - 0.5 * u0 * (u * u + v * v));
            double velocity = sqrt(gam * p / u0) + sqrt(u * u + v * v);  //feature speed
            if (velocity > maxvel)
            {
                maxvel = velocity;
            }
        }
    }

    *dt = sf * min(dx, dy) / maxvel;

}

//deal with the boundary
__global__ void bound_kernel(double* d_U, double gam, int col, int row)
{
    const double rou1 = 1.0, u1 = 2.9, v1 = 0.0, p1 = 0.71429;  //uniform inlet condition at left boundary
    const double rou2 = 1.69997, u2 = 2.61934, v2 = -0.50632, p2 = 1.52819;  //uniform inlet condition at up boundary

    unsigned long long tid = (blockIdx.x * blockDim.x + threadIdx.x);
    int xx = tid % col;
    int yy = tid / col;
    int size = col * row;

    if (xx >= col || yy >= row)
    {
        return;
    }
        
    //left
    if (xx == 0 && yy <= gy + 1)
    {
        d_U[tid] = rou1;
        d_U[tid + size] = rou1 * u1;
        d_U[tid + size * 2] = rou1 * v1;
        d_U[tid + size * 3] = p1 / (gam - 1) + rou1 * (u1 * u1 + v1 * v1) / 2;
    }

    //right
    if (xx == gx + 1 && yy <= gy + 1)
    {
        for (int k = 0; k < 4; k++)
        {
            d_U[tid + k * size] = d_U[gx + yy * col + k * size];
        }
    }

    //up 
    if (yy == gy + 1 && xx <= gx + 1)
    {
        d_U[tid] = rou2;
        d_U[tid + size] = rou2 * u2;
        d_U[tid + size * 2] = rou2 * v2;
        d_U[tid + size * 3] = p2 / (gam - 1) + rou2 * (u2 * u2 + v2 * v2) / 2;
    }

    //down
    if (yy == 0 && xx <= gx + 1)
    {
        d_U[tid] = d_U[xx + 1 * col];
        d_U[tid + size] = d_U[xx + 1 * col + size];
        d_U[tid + size * 2] = 0;
        d_U[tid + size * 3] = d_U[xx + 1 * col + size * 3];
    }
}

void call_bound(double* d_U, double gam, int col, int row, int TPB)
{
    int BSIZE = (col * row + (TPB - 1)) / TPB;
    bound_kernel << <BSIZE, TPB >> > (d_U, gam, col, row);
    cudaDeviceSynchronize();
}

//differential in x-direction
__global__ void updataU_kernel(double* d_U, double* temp, double dx, double dy, double dt, double gam, int col, int row)
{
    const int a = 3.0;  //speed of sound not exceeding 3
    double eta = (a * dt / dx) * (1 - a * dt / dx);

    unsigned long long tid = (blockIdx.x * blockDim.x + threadIdx.x);

    int xx = tid % col;
    int yy = tid / col;
    int size = col * row;
 
    if (xx < 1 || xx > gx || yy < 0 || yy > gy + 1)
    {
        return;
    }
        
    //switching function
    double theta = fabs(fabs(d_U[xx + 1 + yy * col] - d_U[tid]) - fabs(d_U[tid] - d_U[xx - 1 + yy * col]))
        / (fabs(d_U[xx + 1 + yy * col] - d_U[tid]) + fabs(d_U[tid] - d_U[xx - 1 + yy * col]) + 1e-100);

    for (int k = 0; k < 4; k++)
    {
        temp[tid + k * size] = d_U[tid + k * size] + 0.5 * eta * theta * (d_U[xx + 1 + yy * col + k * size]
            - 2 * d_U[tid + k * size] + d_U[xx - 1 + yy * col + k * size]);
    }
}

__global__ void updataU2_kernel(double* d_U, double* d_F, double dx, double dy, double dt, double gam, int col, int row)
{
    unsigned long long tid = (blockIdx.x * blockDim.x + threadIdx.x);

    int xx = tid % col;
    int yy = tid / col;
    int size = col * row;

    if (xx < 1 || xx > gx || yy < 0 || yy > gy + 1)
    {
        return;
    }
        
    for (int i = 0; i < 4; i++)
    {
        d_U[tid + i * size] = d_U[tid + i * size] - (dt / dx) *
            (d_F[tid + i * size] - d_F[xx - 1 + yy * col + i * size]);

    }
}

__global__ void temp2U_kernel(double* d_U, double* temp, double dx, double dy, double dt, double gam, int col, int row)
{
    unsigned long long tid = (blockIdx.x * blockDim.x + threadIdx.x);

    int xx = tid % col;
    int yy = tid / col;
    int size = col * row;

    if (xx < 1 || xx > gx || yy < 0 || yy > gy + 1)
    {
        return;
    }
        
    for (int k = 0; k < 4; k++)
    {
        d_U[tid + k * size] = temp[tid + k * size];
    }
}

__global__ void updataF_kernel(double* d_U, double* d_F, double dx, double dy, double dt, double gam, int col, int row)
{
    unsigned long long tid = (blockIdx.x * blockDim.x + threadIdx.x);

    int xx = tid % col;
    int yy = tid / col;
    int size = col * row;

    if (xx < 0 || xx > gx + 1 || yy < 0 || yy > gy + 1)
    {
        return;
    }
        
    double* u0 = d_U;
    double* u1 = d_U + size;
    double* u2 = d_U + 2 * size;
    double* u3 = d_U + 3 * size;

    double u = u1[tid] / u0[tid];
    double v = u2[tid] / u0[tid];
    double p = (gam - 1) * (u3[tid] - 0.5 * u0[tid] * (u * u + v * v));

    d_F[tid] = u1[tid];
    d_F[tid + size] = u0[tid] * u * u + p;
    d_F[tid + size * 2] = u0[tid] * u * v;
    d_F[tid + size * 3] = (u3[tid] + p) * u;
}

__global__ void updataF_kernel2(double* d_U, double* d_F, double dx, double dy, double dt, double gam, int col, int row)
{
    unsigned long long tid = (blockIdx.x * blockDim.x + threadIdx.x);

    int xx = tid % col;
    int yy = tid / col;
    int size = col * row;

    if (xx < 0 || xx > gx || yy < 0 || yy > gy + 1)
    {
        return;
    }
        
    double* u0 = d_U;
    double* u1 = d_U + size;
    double* u2 = d_U + 2 * size;
    double* u3 = d_U + 3 * size;

    double u = u1[tid] / u0[tid];
    double v = u2[tid] / u0[tid];
    double p = (gam - 1) * (u3[tid] - 0.5 * u0[tid] * (u * u + v * v));

    d_F[tid] = u1[tid];
    d_F[tid + size] = u0[tid] * u * u + p;
    d_F[tid + size * 2] = u0[tid] * u * v;
    d_F[tid + size * 3] = (u3[tid] + p) * u;
}

__global__ void updataUhalf_kernel(double* d_U, double* d_U_half, double* d_F, double dx, double dy, double dt, double gam, int col, int row)
{
    unsigned long long tid = (blockIdx.x * blockDim.x + threadIdx.x);

    int xx = tid % col;
    int yy = tid / col;
    int size = col * row;

    if (xx < 0 || xx > gx || yy < 0 || yy > gy + 1)
    {
        return;
    }
       
    for (int i = 0; i < 4; i++)
    {
        d_U_half[tid + i * size] = 0.5 * (d_U[xx + 1 + yy * col + i * size]
            + d_U[tid + i * size]) - 0.5 * dt / dx * (d_F[xx + 1 + yy * col + i * size]
                - d_F[tid + i * size]);

    }
}

//differential in y-direction
__global__ void updataU_kernel_y(double* d_U, double* temp, double dx, double dy, double dt, double gam, int col, int row)
{
    const int a = 3.0;  //speed of sound not exceeding 3
    double eta = (a * dt / dx) * (1 - a * dt / dx);

    unsigned long long tid = (blockIdx.x * blockDim.x + threadIdx.x);

    int xx = tid % col;
    int yy = tid / col;
    int size = col * row;

    if (xx < 0 || xx > gx + 1 || yy < 1 || yy > gy)
    {
        return;
    }
       
    //switching function
    double theta = fabs(fabs(d_U[xx + (yy + 1) * col] - d_U[tid]) - fabs(d_U[tid] - d_U[xx + (yy - 1) * col]))
        / (fabs(d_U[xx + (yy + 1) * col] - d_U[tid]) + fabs(d_U[tid] - d_U[xx + (yy - 1) * col]) + 1e-100);

    for (int k = 0; k < 4; k++)
    {
        temp[tid + k * size] = d_U[tid + k * size] + 0.5 * eta * theta * (d_U[xx + (yy + 1) * col + k * size]
            - 2 * d_U[tid + k * size] + d_U[xx + (yy - 1) * col + k * size]);
    }
}

__global__ void updataU2_kernel_y(double* d_U, double* d_G, double dx, double dy, double dt, double gam, int col, int row)
{
    unsigned long long tid = (blockIdx.x * blockDim.x + threadIdx.x);

    int xx = tid % col;
    int yy = tid / col;
    int size = col * row;

    if (xx < 0 || xx > gx + 1 || yy < 1 || yy > gy)
    {
        return;
    }
        
    for (int i = 0; i < 4; i++)
    {
        d_U[tid + i * size] = d_U[tid + i * size] - (dt / dx) *
            (d_G[tid + i * size] - d_G[xx + (yy - 1) * col + i * size]);

    }
}

__global__ void temp2U_kernel_y(double* d_U, double* temp, double dx, double dy, double dt, double gam, int col, int row)
{
    unsigned long long tid = (blockIdx.x * blockDim.x + threadIdx.x);

    int xx = tid % col;
    int yy = tid / col;
    int size = col * row;

    if (xx < 0 || xx > gx + 1 || yy < 1 || yy > gy)
    {
        return;
    }
        
    for (int k = 0; k < 4; k++)
    {
        d_U[tid + k * size] = temp[tid + k * size];
    }
}

__global__ void updataG_kernel_y(double* d_U, double* d_G, double dx, double dy, double dt, double gam, int col, int row)
{
    unsigned long long tid = (blockIdx.x * blockDim.x + threadIdx.x);

    int xx = tid % col;
    int yy = tid / col;
    int size = col * row;

    if (xx < 0 || xx > gx + 1 || yy < 0 || yy > gy + 1)
    {
        return;
    }
        
    double* u0 = d_U;
    double* u1 = d_U + size;
    double* u2 = d_U + 2 * size;
    double* u3 = d_U + 3 * size;

    double u = u1[tid] / u0[tid];
    double v = u2[tid] / u0[tid];
    double p = (gam - 1) * (u3[tid] - 0.5 * u0[tid] * (u * u + v * v));

    d_G[tid] = u2[tid];
    d_G[tid + size] = u0[tid] * u * v;
    d_G[tid + size * 2] = u0[tid] * v * v + p;
    d_G[tid + size * 3] = (u3[tid] + p) * v;
}

__global__ void updataG_kernel2_y(double* d_U, double* d_G, double dx, double dy, double dt, double gam, int col, int row)
{
    unsigned long long tid = (blockIdx.x * blockDim.x + threadIdx.x);

    int xx = tid % col;
    int yy = tid / col;
    int size = col * row;

    if (xx < 0 || xx > gx + 1 || yy < 0 || yy > gy)
    {
        return;
    }
        
    double* u0 = d_U;
    double* u1 = d_U + size;
    double* u2 = d_U + 2 * size;
    double* u3 = d_U + 3 * size;

    double u = u1[tid] / u0[tid];
    double v = u2[tid] / u0[tid];
    double p = (gam - 1) * (u3[tid] - 0.5 * u0[tid] * (u * u + v * v));

    d_G[tid] = u2[tid];
    d_G[tid + size] = u0[tid] * u * v;
    d_G[tid + size * 2] = u0[tid] * v * v + p;
    d_G[tid + size * 3] = (u3[tid] + p) * v;
}

__global__ void updataUhalf_kernel_y(double* d_U, double* d_U_half, double* d_G, double dx, double dy, double dt, double gam, int col, int row)
{
    unsigned long long tid = (blockIdx.x * blockDim.x + threadIdx.x);

    int xx = tid % col;
    int yy = tid / col;
    int size = col * row;

    if (xx < 0 || xx > gx + 1 || yy < 0 || yy > gy)
    {
        return;
    }
        
    for (int i = 0; i < 4; i++)
    {
        d_U_half[tid + i * size] = 0.5 * (d_U[xx + (yy + 1) * col + i * size]
            + d_U[tid + i * size]) - 0.5 * dt / dx * (d_G[xx + (yy + 1) * col + i * size]
                - d_G[tid + i * size]);
    }
}

//Lax-Wendroff 2d sovler
void call_solve_x(double* d_U, double* d_U_half, double* d_FG, double* d_temp, double dx, double dy, double dt, double gam, int col, int row, int TPB)
{
    int BSIZE = (col * row + (TPB - 1)) / TPB;

    cudaMemset(d_temp, 0, sizeof(double) * col * row * 4);

    updataU_kernel << <BSIZE, TPB >> > (d_U, d_temp, dx, dy, dt, gam, col, row);

    temp2U_kernel << <BSIZE, TPB >> > (d_U, d_temp, dx, dy, dt, gam, col, row);

    updataF_kernel << <BSIZE, TPB >> > (d_U, d_FG, dx, dy, dt, gam, col, row);

    updataUhalf_kernel << <BSIZE, TPB >> > (d_U, d_U_half, d_FG, dx, dy, dt, gam, col, row);

    updataF_kernel2 << <BSIZE, TPB >> > (d_U_half, d_FG, dx, dy, dt, gam, col, row);

    updataU2_kernel << <BSIZE, TPB >> > (d_U, d_FG, dx, dy, dt, gam, col, row);
}

void call_solve_y(double* d_U, double* d_U_half, double* d_FG, double* d_temp, double dx, double dy, double dt, double gam, int col, int row, int TPB)
{
    int BSIZE = (col * row + (TPB - 1)) / TPB;

    cudaMemset(d_temp, 0, sizeof(double) * col * row * 4);

    updataU_kernel_y << <BSIZE, TPB >> > (d_U, d_temp, dx, dy, dt, gam, col, row);

    temp2U_kernel_y << <BSIZE, TPB >> > (d_U, d_temp, dx, dy, dt, gam, col, row);

    updataG_kernel_y << <BSIZE, TPB >> > (d_U, d_FG, dx, dy, dt, gam, col, row);

    updataUhalf_kernel_y << <BSIZE, TPB >> > (d_U, d_U_half, d_FG, dx, dy, dt, gam, col, row);

    updataG_kernel2_y << <BSIZE, TPB >> > (d_U_half, d_FG, dx, dy, dt, gam, col, row);

    updataU2_kernel_y << <BSIZE, TPB >> > (d_U, d_FG, dx, dy, dt, gam, col, row);
}

void call_solver_2d(double* d_U, double* d_U_half, double* d_FG, double* d_temp, double dx, double dy, double dt, double gam, int col, int row, int TPB)
{
    call_solve_x(d_U, d_U_half, d_FG, d_temp, dx, dy, dt / 2.0, gam, col, row, TPB);
    call_bound(d_U, gam, col, row, TPB);

    call_solve_y(d_U, d_U_half, d_FG, d_temp, dx, dy, dt / 2.0, gam, col, row, TPB);
    call_bound(d_U, gam, col, row, TPB);

    call_solve_y(d_U, d_U_half, d_FG, d_temp, dx, dy, dt / 2.0, gam, col, row, TPB);
    call_bound(d_U, gam, col, row, TPB);

    call_solve_x(d_U, d_U_half, d_FG, d_temp, dx, dy, dt / 2.0, gam, col, row, TPB);
    call_bound(d_U, gam, col, row, TPB);
}
