#include "MatBand.h"
#include "util.h"

MatBand::MatBand() : Mat<double>(), u(0), l(0) {;}

MatBand::MatBand(int m, int u_, int l_) // square matrix size is (m,m)
: Mat<double>(m, u_+l_+1), u(u_), l(l_) {;}

MatBand::MatBand(int m, int u_, int l_, double a) // assign a to every element
: Mat<double>(m, u_+l_+1, a), u(u_), l(l_) {;}

MatBand::MatBand(const MatBand& a) // resized to the size of a
: Mat<double>(a), u(a.u), l(a.l) {;}

MatBand& MatBand::operator=(const MatBand& a)
{
    Mat<double>::operator=(a);
    u = a.u;
    l = a.l;
    return *this;
}

void MatBand::resize(int m, int u_, int l_)
// data remain unchanged
{
    u = u_;
    l = l_;
    Mat<double>::resize(m, u+l+1);
}

void MatBand::resize(int m, int u, int l, double a)
// assign a to every element
{
    resize(m, u, l);
    Mat<double>::operator=(a);
}

void MulAddTo(vec_double& y, const MatBand& a, const vec_double& x)
// y += Ax
{
    if(&y==&x) {
        vec_double z(x);
        MulAddTo(y,a,z);
        return;
    }
    int i,j,k,n,m(a.nrows());
    int u(a.upper_width()), l(a.lower_width());
    for(i=0; i<m; i++) {
        k = max(i-l, 0);
        n = min(i+u+1, m);
        for(j=k; j<n; j++) y[i] += a(i,j)*x[j];
    }
}
