#include "crvLin.h"
#include <stdio.h>
#include <math.h>

struct MetricCoefficients {
    double G11, G12, G22, G;
};

void computeMetricCoefficients(int i, int j, int nX, int nY, double ** x, double ** y, MetricCoefficients & metrCof) {
    double hKsi, hEta;
    double dXdKsi, dXdEta, dYdKsi, dYdEta;

    hKsi = 1.0 / (nX - 1);
    hEta = 1.0 / (nY - 1);

    dXdKsi = (x[i + 1][j] - x[i - 1][j]) / (2.0 * hKsi);
    dYdKsi = (y[i + 1][j] - y[i - 1][j]) / (2.0 * hKsi);
    dXdEta = (x[i][j + 1] - x[i][j - 1]) / (2.0 * hEta);
    dYdEta = (y[i][j + 1] - y[i][j - 1]) / (2.0 * hEta);

    metrCof.G11 = dXdKsi * dXdKsi + dYdKsi * dYdKsi;
    metrCof.G12 = dXdKsi * dXdEta + dYdKsi * dYdEta;
    metrCof.G22 = dXdEta * dXdEta + dYdEta * dYdEta;
    metrCof.G   = dXdKsi * dYdEta - dXdEta * dYdKsi;
}

double fLaplace(int i, int j, int nX, int nY, MetricCoefficients & metrCof, double ** f) {
    double hKsi = 1.0 / (nX - 1);
    double hEta = 1.0 / (nY - 1);

    return (metrCof.G22 * (f[i + 1][j] + f[i - 1][j]) +
            metrCof.G11 * (f[i][j + 1] + f[i][j - 1]) * pow(hKsi / hEta, 2) -
            metrCof.G12 * (f[i + 1][j + 1] - f[i - 1][j + 1] - f[i + 1][j - 1] + f[i - 1][j - 1]) * 0.5 * hKsi / hEta) /
            (2.0 * (metrCof.G22 + metrCof.G11 * pow(hKsi / hEta, 2)));
}

double relaxation(double fOld, double fNew, double alpha) {
    return fOld * (1.0 - alpha) + fNew * alpha;
}

void computeGrid(int nX, int nY, int nXMin, int nYMin, int nXMax, int nYMax, double ** x, double ** y, double eps) {
    const double alpha = 1.8;
    double epsIter = 0.0;
    MetricCoefficients metrCofij;
    double xNew, yNew;
    double epsXij, epsYij;
    int kIter = 0;
    do {
        kIter++;
        epsIter = 0.0;
        for (int j = 1; j < nY - 1; j++) {
            for (int i = 1; i < nX - 1; i++) {
                if ((i >= nXMin) && (i <= nXMax) && (j >= nYMin) && (j <= nYMax)) continue;

                computeMetricCoefficients(i, j, nX, nY, x, y, metrCofij);

                xNew = fLaplace(i, j, nX, nY, metrCofij, x);
                yNew = fLaplace(i, j, nX, nY, metrCofij, y);

                epsXij = fabs(x[i][j] - xNew);
                epsYij = fabs(y[i][j] - yNew);

                if (epsIter < epsXij) epsIter = epsXij;
                if (epsIter < epsYij) epsIter = epsYij;

                x[i][j] = relaxation(x[i][j], xNew, alpha);
                y[i][j] = relaxation(y[i][j], yNew, alpha);
                
            }
        }
        printf("%d | %le\n", kIter, epsIter);
        if (kIter == 1000) break;
    } while (epsIter > eps);
}
