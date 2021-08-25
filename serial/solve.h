#include "init.h"

void solver_2d(double ***U, double ***U_half, double ***F, double ***G, double ***temp, double dx, double dy, double dt, double gam);
void solve_x(double ***U, double ***U_half, double ***F, double ***temp, double dx, double dt, double gam);
void solve_y(double ***U, double ***U_half, double ***G, double ***temp, double dy, double dt, double gam);
