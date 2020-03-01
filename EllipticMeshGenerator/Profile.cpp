#include "Profile.h"
#include <stdio.h>
#include <stdlib.h>

void setProfile(int nX, int nY, int & nXMin, int & nYMin, int & nXMax, int & nYMax, double ** x, double ** y, 
                double & xc, double & yc, int & iRib, int & jRib, const char * fileName) {
    nXMin = 5;
    nYMin = 8;
    nXMax = 13;
    nYMax = 12;

    iRib = nXMax - 1;
    jRib = 10;

    FILE * f = fopen(fileName, "r");
    for (int j = nYMin; j < nYMax; j++) {
        if (fscanf(f, "%lf%lf", &x[nXMin][j], &y[nXMin][j]) != 2) exit(1);
    }
    for (int i = nXMin; i < nXMax; i++) {
        if (fscanf(f, "%lf%lf", &x[i][nYMax], &y[i][nYMax]) != 2) exit(1);
    }
    for (int j = nYMax; j > nYMin; j--) {
        if (fscanf(f, "%lf%lf", &x[nXMax][j], &y[nXMax][j]) != 2) exit(1);
    }
    for (int i = nXMax; i > nXMin; i--) {
        if (fscanf(f, "%lf%lf", &x[i][nYMin], &y[i][nYMin]) != 2) exit(1);
    }
    xc = 0.0;
    yc = 0.0;
    fclose(f);
}
