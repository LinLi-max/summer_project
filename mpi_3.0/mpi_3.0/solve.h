#include "init.h"

void solver_2d(double U[4][gx + 2][gy + 2], double U_half[4][gx + 2][gy + 2], double F_G[4][gx + 2][gy + 2], double dx, double dy, double dt, double gam, int s, int e, int left, int right, MPI_Comm comm);
void solve_x(double U[4][gx + 2][gy + 2], double U_half[4][gx + 2][gy + 2], double F[4][gx + 2][gy + 2], double dx, double dt, double gam, int s, int e);
void solve_y(double U[4][gx + 2][gy + 2], double U_half[4][gx + 2][gy + 2], double G[4][gx + 2][gy + 2], double dy, double dt, double gam, int s, int e);
void exchang(double U[4][gx + 2][gy + 2], int s, int e, int left, int right, MPI_Comm comm);
void gathergrid(double U[4][gx + 2][gy + 2], int rank, int size, int s, int e, int left, int right);