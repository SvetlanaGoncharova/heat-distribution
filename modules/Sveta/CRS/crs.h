#include <stdio.h>
#include <stdlib.h>

typedef struct {
    double* Value; //�������� �������
    int* col;//������ �������� ��� ��������� ���������
    int* rowindex;//������ ������
    int nz;//���������� ��������� ���������
    int nRows;
} CRSMatrix;

void initCRSMartix(int nRows, int nz, CRSMatrix *crsm);
void freeCRSMatrix(CRSMatrix *);
void multCRSMatrix(double** result, CRSMatrix crsm, double* vec);