#ifndef __MatBand_h__
#define __MatBand_h__

#include "Mat.h"

struct MatBand : Mat<double> {// banded (and square) matrix
    int u,l;// upper and lower width
    MatBand();
    MatBand(int m, int u, int l); // square matrix size is (m,m)
    MatBand(int m, int u, int l, double a);// assign a to every element
    MatBand(const MatBand& a);// copy constructor
    MatBand& operator=(const MatBand& a);
    void resize(int m, int u, int l);// data remain unchanged
    void resize(int m, int u, int l, double a);// assign a to every element
    virtual inline double& operator()(int i, int j) {
        return (*this)[j][i-j+u]; // i-th row, j-th column
    }
    virtual inline const double& operator()(int i, int j) const {
        return (*this)[j][i-j+u]; // i-th row, j-th column
    }
    inline int upper_width() const { return u; }
    inline int lower_width() const { return l; }
};

void MulAddTo(vec_double& y, const MatBand& a, const vec_double& x);
// y += ax; &y==&x is allowed
// assume y.size == x.size == a.nrows

#endif // __MatBand_h__
