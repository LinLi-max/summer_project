#include "init.h"

void solver_2d(double ***U, double ***U_half, double ***F, double ***G, double ***temp, double dx, double dy, double dt, double gam, int s, int e, int left, int right, int rank, int size, MPI_Win win);
void solve_x(double ***U, double ***U_half, double ***F, double ***temp, double dx, double dt, double gam, int s, int e, int left, int right, MPI_Win win);
void solve_y(double ***U, double ***U_half, double ***G, double ***temp, double dy, double dt, double gam, int s, int e);
void exchange(double ***U, int s, int e, int left, int right, MPI_Win win);
void gathergrid(double ***U, int rank, int size, int s, int e, int left, int right);
