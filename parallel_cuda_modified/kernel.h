#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define gx 1600
#define gy 400

void call_solver_2d(double* d_U, double* d_U_half, double* d_FG, double* d_temp, double dx, double dy, double dt, double gam, int col, int row, int TPB);
void call_init(double* d_U, double dx, double dy, double gam, int col, int row, int TPB);
void call_cfl(double* d_U, double dx, double dy, double* dt, double sf, double gam, int col, int row);
