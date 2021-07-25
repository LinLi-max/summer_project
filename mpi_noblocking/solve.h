#include "init.h"

void solver_2d(double ***U, double ***U_half, double ***F, double ***G, double ***temp, double dx, double dy, double dt, double gam, int s, int e, int left, int right, int rank, int size, MPI_Comm comm);
void exchang1(double ***U, double ***temp, double dx, double dt, int s, int e, int left, int right, MPI_Comm comm);
void exchang2(double ***U, double ***U_half, double ***F, double dx, double dt, double gam, int s, int e, int left, int right, MPI_Comm comm);
void solve_y(double ***U, double ***U_half, double ***G, double ***temp, double dy, double dt, double gam, int s, int e);
void gathergrid(double ***U, int rank, int size, int s, int e, int left, int right);
