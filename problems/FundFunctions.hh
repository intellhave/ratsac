#ifndef FTOOLS_H
#define FTOOLS_H

#include <vector>
#include "MathFunctions.hh"

namespace FTools
{
void normalizePoints(double* inputPointsUnNorm, double* inputPoints, int numPoints,
                     double* T1, double* T2);
void computeDataMatrix(double* data_matrix, int num_points, double* points);
int nullspaceQR7x9(const double* A, double* N);
int nullspace(double* matrix, double* nullspace, int n, int* buffer);
void makePolynomial(double* A, double* B, double* p);
int rroots3 (double* po, double* r);
void formCovMat(double* Cv, const double* A, int len, int siz);
void singulF(double* F);
void computeEpipole(double* e, const double* F);
double getOriSign(double* F, double* e, double* pt);
void computeHFromF(const std::vector<int>& sample, double* u, double* ep, double* F, double* H);
int getHError(std::vector<double>& errs, double* u, int nr, double* H, double threshold, const std::vector<int> *test=NULL);
int computeHFromCorrs(const std::vector<int>& sample, int numPoints,
                               int numDataPoints, double* u, double* H);
int computeHFromMinCorrs(const std::vector<int>& sample, int numPoints,
                                  int numDataPoints, double* u, double* H);
}
#endif
