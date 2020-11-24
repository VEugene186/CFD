#include <stdio.h>
#include <math.h>
#include "boundaryConditions.h"
#include "source.h"
#include <omp.h>

const int N = 200;
const int M = 400;
const double tol = 1e-8;

double ap[N][M], aw[N][M], ae[N][M], as[N][M], an[N][M], b[N][M];
double xc[N], yc[M];
double dx, dy;

double T[N][M], R[N][M], p[N][M], AP[N][M];

void makeCoefficients();
double residual();
void init();
void overRelaxationStd();
void overRelaxationRedBlack();
void conjugateGradients();
void write();

int main(int argc, char * argv[]) {
    makeCoefficients();
    init();
    int num = -1;
    printf("1 - overRelaxationRedBlack\n");
    printf("2 - overRelaxationStd\n");
    printf("3 - conjugateGradient\n");
    printf("Select: ");
    scanf("%d", &num);
    switch (num) {
    case 1:
        overRelaxationRedBlack();
        break;
    case 2:
        overRelaxationStd();
        break;
    case 3:
        conjugateGradients();
        break;
    }
    write();
    return 0;
}

void makeCoefficients() {
    dx = 1.0 / N;
    dy = 1.0 / M;
    //Centers of cells
    xc[0] = 0.5 * dx;
    for (int i = 1; i < N; i++) {
        xc[i] = xc[i - 1] + dx;
    }
    yc[0] = 0.5 * dy;
    for (int j = 1; j < M; j++) {
        yc[j] = yc[j - 1] + dy;
    }
    //Coefficients
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            b[i][j] = heatSource(xc[i], yc[j]) * dx * dy;
            
            if (i > 0) {
                aw[i][j] = dy / dx;
            }
            else {
                aw[i][j] = 2.0 * dy / dx;
                b[i][j] += aw[i][j] * leftVal(yc[j]);
            }

            if (i < N - 1) {
                ae[i][j] = dy / dx;
            }
            else {
                ae[i][j] = 2.0 * dy / dx;
                b[i][j] += ae[i][j] * rightVal(yc[j]);
            }

            if (j > 0) {
                as[i][j] = dx / dy;
            }
            else {
                as[i][j] = 0.0;
                b[i][j] -= bottomDer(xc[i]) * dx;
            }

            if (j < M - 1) {
                an[i][j] = dx / dy;
            }
            else {
                an[i][j] = 2.0 * dx / dy;
                b[i][j] += an[i][j] * topVal(xc[i]);
            }

            ap[i][j] = aw[i][j] + ae[i][j] + as[i][j] + an[i][j];
        }
    }
}

double residual() {
    double maxR = 0.0, absR;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            R[i][j] = b[i][j] - ap[i][j] * T[i][j];
            if (i > 0) R[i][j] += aw[i][j] * T[i - 1][j];
            if (i < N - 1) R[i][j] += ae[i][j] * T[i + 1][j];
            if (j > 0) R[i][j] += as[i][j] * T[i][j - 1];
            if (j < M - 1) R[i][j] += an[i][j] * T[i][j + 1];
            absR = fabs(R[i][j]);
            if (maxR < absR) maxR = absR;
        }
    }
    return maxR / (dx * dy);
}

void init() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            T[i][j] = 0.0;
        }
    }
}

void overRelaxationStd() {
    double w = 1.0;
    double maxR = 0.0;
    double tmp;
    int nIter = 0;
    printf("w = ");
    scanf("%lf", &w);
    
    double st = omp_get_wtime();
    while (true) {
        maxR = residual();
        if (nIter % 10 == 0) printf("nIter = %6d | maxR = %.4le\n", nIter, maxR);
        if (maxR < tol) {
            break;
        }
        
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                tmp = b[i][j];
                if (i > 0) tmp += aw[i][j] * T[i - 1][j];
                if (i < N - 1) tmp += ae[i][j] * T[i + 1][j];
                if (j > 0) tmp += as[i][j] * T[i][j - 1];
                if (j < M - 1) tmp += an[i][j] * T[i][j + 1];
                T[i][j] = w * tmp / ap[i][j] + (1.0 - w) * T[i][j];
            }
        }
        nIter++; 
    }
    double fin = omp_get_wtime();
    printf("===Final===\n");
    printf("nIter = %6d | maxR = %.4le\n", nIter, maxR);
    printf("Time = %lg\n", fin - st);
}

void overRelaxationRedBlack() {
    double w = 1.0;
    double maxR = 0.0;
    double tmp;
    int nIter = 0;
    printf("w = ");
    scanf("%lf", &w);
    
    double st = omp_get_wtime();
    while (true) {
        maxR = residual();
        if (nIter % 10 == 0) printf("nIter = %6d | maxR = %.4le\n", nIter, maxR);
        if (maxR < tol) {
            break;
        }
        
        for (int i = 0; i < N; i++) {
            for (int j = i % 2; j < M; j += 2) {
                tmp = b[i][j];
                if (i > 0) tmp += aw[i][j] * T[i - 1][j];
                if (i < N - 1) tmp += ae[i][j] * T[i + 1][j];
                if (j > 0) tmp += as[i][j] * T[i][j - 1];
                if (j < M - 1) tmp += an[i][j] * T[i][j + 1];
                T[i][j] = w * tmp / ap[i][j] + (1.0 - w) * T[i][j];
            }
        }
 
        for (int i = 0; i < N; i++) {
            for (int j = (i + 1) % 2; j < M; j += 2) {
                tmp = b[i][j];
                if (i > 0) tmp += aw[i][j] * T[i - 1][j];
                if (i < N - 1) tmp += ae[i][j] * T[i + 1][j];
                if (j > 0) tmp += as[i][j] * T[i][j - 1];
                if (j < M - 1) tmp += an[i][j] * T[i][j + 1];
                T[i][j] = w * tmp / ap[i][j] + (1.0 - w) * T[i][j];
            }
        }
        nIter++; 
    }
    double fin = omp_get_wtime();
    printf("===Final===\n");
    printf("nIter = %6d | maxR = %.4le\n", nIter, maxR);
    printf("Time = %lg\n", fin - st);
}

void conjugateGradients() {
    double st = omp_get_wtime();
    int k = 0;
    double maxR = residual(), absR;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            p[i][j] = R[i][j];
        }
    }
    printf("nIter = %6d | maxR = %.4le\n", k, maxR);

    double App = 0.0, RR1 = 0.0, RR2 = 0.0, alpha = 0.0, beta = 0.0;
    while (true) {
        App = 0.0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                AP[i][j] = ap[i][j] * p[i][j];
                if (i > 0) AP[i][j] -= aw[i][j] * p[i - 1][j];
                if (i < N - 1) AP[i][j] -= ae[i][j] * p[i + 1][j];
                if (j > 0) AP[i][j] -= as[i][j] * p[i][j - 1];
                if (j < M - 1) AP[i][j] -= an[i][j] * p[i][j + 1];
                
                App += AP[i][j] * p[i][j];
                if (k == 0) RR1 += R[i][j] * R[i][j];
            }
        }
        alpha = RR1 / App;
        maxR = 0.0;
        RR2 = 0.0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                T[i][j] += alpha * p[i][j];
                R[i][j] -= alpha * AP[i][j];
                absR = fabs(R[i][j]);
                if (maxR < absR) maxR = absR;
                RR2 += R[i][j] * R[i][j];
            }
        }
        maxR /= (dx * dy);
        k++;
        if (k % 10 == 0) printf("nIter = %6d | maxR = %.4le\n", k, maxR);
        if (maxR < tol) break;
        beta = RR2 / RR1;
        RR1 = RR2;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                p[i][j] = p[i][j] * beta + R[i][j];
            }
        }
    }
    double fin = omp_get_wtime();
    printf("===Final===\n");
    printf("nIter = %6d | maxR = %.4le\n", k, maxR);
    printf("Time = %lg\n", fin - st);
}

void write() {
    FILE * fX = fopen("x.csv", "w");
    FILE * fY = fopen("y.csv", "w");
    FILE * fT = fopen("T.csv", "w");
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            char delim = (j < M - 1) ? '\t' : '\n';
            fprintf(fX, "%.15lg%c", xc[i], delim);
            fprintf(fY, "%.15lg%c", yc[j], delim);
            fprintf(fT, "%.15lg%c", T[i][j], delim);
        }
    }

    fclose(fX);
    fclose(fY);
    fclose(fT);
}
