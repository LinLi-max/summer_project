#include <math.h>
#include "grid.h"

#define min(x, y) (((x) < (y)) ? (x) : (y))

void cfl(double ***U, double dx, double dy, double* dt, double sf, double gam);
void init(double ***U, double dx, double dy, double gam);
void bound(double ***U, double dx, double dy, double gam);
