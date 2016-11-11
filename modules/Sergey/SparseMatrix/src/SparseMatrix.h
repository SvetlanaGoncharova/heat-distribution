//
// Created by lenferd on 27.10.16.
//

#ifndef SPARSEMATRIX_SPARSEMATRIX_H
#define SPARSEMATRIX_SPARSEMATRIX_H
#include <vector>
#include <algorithm>
#include <omp.h>

const int ENABLE_PARALLEL = 0;
using std::vector;

class SparseMatrix {
private:
    vector<double> values;
    vector<int> columns;
    vector<int> pointerB;
    vector<int> pointerE;

public:
    SparseMatrix() {};
    void fillMatrix(double** &matrix, int widthSize, int heightSize);
    void printVectors();
    double* multiplicateVector(vector<double> vect);
    double* multiplicateVector(double* &vect);
    void multiplicateVector(double* &vect, double* &result, int size);
    void testEuler(int size, double expr);
    void fillMatrix2Expr(int size, double expr1, double expr2);
    void Rungek2(int size, double expr);
};


#endif //SPARSEMATRIX_SPARSEMATRIX_H
