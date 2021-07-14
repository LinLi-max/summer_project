#include "init.h"

void U_F(double U[4], double F[4], double gam);
void U_G(double U[4], double G[4], double gam);
//void solve_x(double U[gx + 2][gy + 2][4], double Uf[gx + 2][gy + 2][4], double Ff[gx + 2][gy + 2][4], double dx, double dt, double gam);
//void solve_y(double U[gx + 2][gy + 2][4], double Uf[gx + 2][gy + 2][4], double Gf[gx + 2][gy + 2][4], double dy, double dt, double gam);
//void solver_2d(double U[gx + 2][gy + 2][4], double Uf[gx + 2][gy + 2][4], double FGf[gx + 2][gy + 2][4], double dx, double dy, double dt, double gam);
void solver_2d(double U[gx + 2][gy + 2][4], double U_half[gx + 2][gy + 2][4], double F_G[gx + 2][gy + 2][4], double dx, double dy, double dt, double gam, int s, int e, int left, int right, MPI_Comm comm);
void solve_x(double U[gx + 2][gy + 2][4], double U_half[gx + 2][gy + 2][4], double F[gx + 2][gy + 2][4], double dx, double dt, double gam, int s, int e);
void solve_y(double U[gx + 2][gy + 2][4], double U_half[gx + 2][gy + 2][4], double G[gx + 2][gy + 2][4], double dy, double dt, double gam, int s, int e);
void exchang(double U[gx + 2][gy + 2][4], int s, int e, int left, int right, MPI_Comm comm);
void gathergrid(double U[gx + 2][gy + 2][4], int rank, int size, int s, int e, int left, int right);