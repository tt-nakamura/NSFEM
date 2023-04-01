// uses GNU Scientfic Library
//   http://www.gnu.org/software/gsl

#include "MatBandSym.h"
#include "util.h"
#include<gsl/gsl_linalg.h>

MatBandSym::MatBandSym() : MatBand() {;}

MatBandSym::MatBandSym(int m, int w) // square matrix size is (m,m)
: MatBand(m,w,0) { l=w;} // band width is 2w+1

MatBandSym::MatBandSym(int m, int w, double a)
: MatBand(m,w,0,a) { l=w;} // assign a to every element

void MatBandSym::resize(int m, int w)
// data remain unchanged
{
    u = w;
    l = w;
    Mat<double>::resize(m, w+1);
}

void MatBandSym::resize(int m, int w, double a)
// assign a to every element
{
    resize(m,w);
    Mat<double>::operator=(a);
}

void conv(MatBand& a, const MatBandSym& b)
// data conversion
{
    int i,j,k,m(b.nrows()),w(b.half_width());
    a.resize(m,w,w);
    for(i=0; i<m; i++)
        a(i,i) = b(i,i);
    for(k=1; k<=w; k++)
        for(i=k, j=0; i<m; i++, j++)
            a(i,j) = a(j,i) = b(i,j);
}

void LDLTdecomp(MatBandSym& x, const MatBandSym& a)
// x = LDLT decomposition of a
{
    int i,j,m(a.nrows()),n(a.ncols());
    gsl_matrix *p;
    p = gsl_matrix_alloc(m,n);
    for(i=0; i<m; i++)
        for(j=0; j<n; j++)
            gsl_matrix_set(p,i,j,a[i][j]);
    gsl_linalg_ldlt_band_decomp(p);
    x.resize(m, a.half_width());
    for(i=0; i<m; i++)
        for(j=0; j<n; j++)
            x[i][j] = gsl_matrix_get(p,i,j);
    gsl_matrix_free(p);
}

void solve(vec_double& x, const MatBandSym& a, const vec_double& b,
           bool decomp)
// solve ax = b for x
{
    int i,j,m(a.nrows()),n(a.ncols());
    gsl_matrix *p;
    gsl_vector *v;
    p = gsl_matrix_alloc(m,n);
    v = gsl_vector_alloc(m);
    for(i=0; i<m; i++)
        gsl_vector_set(v, i, b[i]);
    for(i=0; i<m; i++)
        for(j=0; j<n; j++)
            gsl_matrix_set(p, i, j, a[i][j]);
    if(decomp) gsl_linalg_ldlt_band_decomp(p);
    gsl_linalg_ldlt_band_svx(p,v);
    if(x.size() != m) x.resize(m);
    for(i=0; i<m; i++)
        x[i] = gsl_vector_get(v,i);
    gsl_matrix_free(p);
    gsl_vector_free(v);
}
