#include <stdio.h>
#include "constants.h"
#include "Profile.h"
#include "crvLin.h"

int nX = 0, nY = 0; //mesh size
int nXMin = 0, nYMin = 0; //cut
int nXMax = 0, nYMax = 0; //cut
double ** x, ** y;
double alpha; //attack angle
double xc, yc; //foil's center of rotation
int iRib, jRib; //sharp edge position
double alphaLast;

void alloc(double *** a) {
    *a = new double*[nMax];
    for (int i = 0; i < nMax; i++) {
        (*a)[i] = new double[nMax];
    }
}

void free(double ** a) {
    for (int i = 0; i < nMax; i++) {
        delete [] a[i];
    }
    delete [] a;
}

int main(int argc, char * argv[]) {
    alloc(&x);
    alloc(&y);

    nX = nMax;
    nY = nMax;
    for (int i = 0; i < nX; i++) {
        for (int j = 0; j < nY; j++) {
            x[i][j] = 0.0;
            y[i][j] = 0.0;
        }
    }

    for (int i = 0; i < nX; i++) {
        x[i][0] = (double)i / (double)(nX - 1);
        x[i][nY - 1] = x[i][0];
        y[i][nY - 1] = 1.0;
    }

    for (int j = 0; j < nY; j++) {
        y[0][j] = (double)j / (double)(nY - 1);
        y[nX - 1][j] = y[0][j];
        x[nX - 1][j] = 1.0;
    }

    printf("start reading\n");
    setProfile(nX, nY, nXMin, nYMin, nXMax, nYMax, x, y, xc, yc, iRib, jRib, "profile.dat");
    printf("end reading\n");
    //turnProfile(nXMin, nYMin, nXMax, nYMa, alpha - alphaLast, xc, yc, x, y);
    computeGrid(nX, nY, nXMin, nYMin, nXMax, nYMax, x, y, eps);

    //FILE * f = fopen("mesh.csv", "w");
    /*for (int i = 0; i < nX; i++) {
        for (int j = 0; j < nY; j++) {
            if ((i >= nXMin) && (i <= nXMax) && (j >= nYMin) && (j <= nYMax)) continue;
            fprintf(f, "%.15lg\t%.15lg\n", x[i][j], y[i][j]);
        }
    }*/
    FILE * f = fopen("x.csv", "w");
    for (int j = 0; j < nY; j++) {
        for (int i = 0; i < nX; i++) {
            fprintf(f, "%.15lg%c", x[i][j], ((i < nX - 1) ? '\t' : '\n'));
        }
    }
    fclose(f);
    f = fopen("y.csv", "w");
    for (int j = 0; j < nY; j++) {
        for (int i = 0; i < nX; i++) {
            fprintf(f, "%.15lg%c", y[i][j], ((i < nX - 1) ? '\t' : '\n'));
        }
    }
    fclose(f);
    free(x);
    free(y);
    return 0;
}
