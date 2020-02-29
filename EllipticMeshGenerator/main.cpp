#include <stdio.h>
#include <valarray>

typedef std::valarray<std::valarray<double> > GridType;

const double eps = 1e-5;

int nX = 0, nY = 0; //mesh size
int nXMin = 0, nYMin = 0; //cut
int nXMax = 0, nYMax = 0; //cut
GridType x, y;
double alpha; //attack angle
double xc, yc; //foil's center of rotation
int iRib, jRib; //sharp edge position
double alphaLast;


int main(int argc, char * argv[]) {

    return 0;
}
