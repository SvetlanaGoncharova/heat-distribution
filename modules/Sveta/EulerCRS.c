
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include "crs.h"

int main() {

    double xStart = 0.0, xEnd = 0.0;
    double sigma = 0.0;
    int nX = 1;
    double tStart = 0.0, tFinal = 0.0;
    double dt = 0.0;
    int check = 0;
    double **U;
    double step;
    int sizeTime = 0;
    int curTime, prevTime;


    FILE *fp;
    if ((fp = fopen("C:\\Users\\ролд\\Documents\\Visual Studio 2013\\Projects\\Eiler\\Eiler\\INPUT.txt", "r")) == NULL) {
        printf("File cant open!\n");
        exit(-1);
    }

    // Заполнение конфигурационных настроек
    fscanf(fp, "XSTART=%lf\n", &xStart);
    fscanf(fp, "XEND=%lf\n", &xEnd);
    fscanf(fp, "SIGMA=%lf\n", &sigma);
    fscanf(fp, "NX=%d\n", &nX);
    fscanf(fp, "TSTART=%lf\n", &tStart);
    fscanf(fp, "TFINISH=%lf\n", &tFinal);
    fscanf(fp, "dt=%lf\n", &dt);
    fscanf(fp, "BC=%d\n", &check);

    *U = (double*)malloc(sizeof(double) * (nX + 2));

    // Заполнение функции в нулевой момент времени
    for (int i = 1; i < nX + 1; i++) {
        fscanf(fp, "%lf", &(*U)[i]);
    }
    (*U)[0] = (*U)[nX + 1] = 0.0;

    fclose(fp);

    //  Вычисление шага по x
    step = fabs(xStart - xEnd) / nX;

    if (2 * dt*sigma > step*step) {
        printf("Not correct step!");
        exit(-1);
    }
    sizeTime = (int)((tFinal - tStart) / dt);
    int backstep = 1 / (step*step);

    //создание матрицы
    CRSMatrix coeff;
    double h = dt*backstep;
    initCRSMartix(nX + 2, 3 * nX + 2, &coeff);
    coeff.Value[0] = 1;
    for (int i = 1; i < nX + 1; i += 3){
        coeff.Value[i] = h;
        coeff.Value[i + 1] = 1 - 2 * h;
        coeff.Value[i + 2] = h;
    }
    coeff.Value[nX + 1] = 1;
    coeff.col[0] = 0;
    coeff.col[nX + 1] = nX + 1;
    int j = 0;
    for (int i = 1; i < nX + 1; i += 3){
        coeff.col[i] = j;
        coeff.col[i + 1] = j + 1;
        coeff.col[i + 2] = j + 2;
        j++;
    }
    coeff.rowindex[0] = 0;
    coeff.rowindex[1] = 1;
    for (int i = 2; i < coeff.nRows; i++)
        coeff.rowindex[i] = coeff.rowindex[i - 1] + 3;
    coeff.rowindex[coeff.nRows + 1] = coeff.rowindex[coeff.nRows] + 1;

    //вычисления
    double* Unext = (double*)malloc(sizeof(double) * (nX + 2));
    double *tmp;
    for (int i = 1; i <= sizeTime; i++) {
        multCRSMatrix(&Unext, coeff, U);

        tmp = U;
        U = Unext;
        Unext = tmp;
    }

    fp = 0;
    fp = fopen("C:\\Users\\ролд\\Desktop\\heat-distribution\\result\\svetaEulerCRS.txt", "w");
    if (!fp)
        return -1;
    for (int x = 1; x < nX + 1; x++)
        fprintf(fp, "%.15le\n", U[sizeTime % 2][x]);
    fclose(fp);
    free(Unext);
    return 0;
}