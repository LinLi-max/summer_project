#include "decomp1d.h"

int decomp1d(int rank, int size, int *s, int *e)
{
	int length = gx / size;
	int remainder = gx % size;
	/* when p divides into gx without remainder */
	if (remainder == 0) {
		*s = 1 + rank * length;
		*e = *s + length - 1;
	} else {
		/* when p divides into gx with remainder */
		if (rank < remainder) {
			*s = rank * (length + 1) + 1;
			*e = *s + length - 1 + 1;
		} else {
			*s = rank * length + remainder + 1;
			*e = *s + length - 1;
		}
	}
	return 0;
}
