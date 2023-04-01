// uses GNU Scientfic Library
//   http://www.gnu.org/software/gsl

#ifndef __MatBandSym_h__
#define __MatBandSym_h__

#include "MatBand.h"

struct MatBandSym : MatBand {// symmetric banded (and square) matrix
    MatBandSym();// square matrix size is (m,m)
    MatBandSym(int m, int w);// band width is 2w+1
    MatBandSym(int m, int w, double a);// assign a to every element
    void resize(int m, int w);// data remain unchanged
    void resize(int m, int w, double a);// assign a to every element
    inline double& operator()(int i, int j) {// i-th row, j-th column
        return i<j ? (*this)[i][j-i] : (*this)[j][i-j];
    }
    inline const double& operator()(int i, int j) const {
        return i<j ? (*this)[i][j-i] : (*this)[j][i-j];
    };
    inline int half_width() const { return u; } // u=l=w
};

void conv(MatBand& a, const MatBandSym& b); // data conversion

void LDLTdecomp(MatBandSym& x, const MatBandSym& a);
// x = LDLT decomposition of a; &x==&a is allowed
// see GSL reference section 14.22.5

void solve(vec_double& x, const MatBandSym& a, const vec_double& b,
           bool decomp=true);
// solve ax = b for x; &x==&b is allowed
// if decomp is true, a is matrix before LDLT decomposition
// else, a must be output of LDLTdecomp

#endif // __MatBand_h__
