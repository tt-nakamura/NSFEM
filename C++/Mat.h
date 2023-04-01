#ifndef __Mat_h__
#define __Mat_h__

#include<vector>

template<class T>
struct Mat : std::vector<T*> {
    int n; // number of columns
    std::vector<T> v; // data buffer
    Mat() : n(0) {;}
    Mat(int m, int n);
    Mat(int m, int n, const T& a); //assign a to every element
    Mat(const Mat& a); // copy constructor
    Mat& operator=(const Mat& a); //assignment
    Mat& operator=(const T& a); //assign a to every element
    void resize(int m, int n); // data remain unchanged
    void resize(int m, int n, const T& a); //assign a to every element
    inline int nrows() const { return this->size(); }
    inline int ncols() const { return n; }
};

template<class T>
Mat<T>::Mat(int m, int n_)
: std::vector<T*>(m), n(n_), v(m*n)
{
    int i,j;
    for(i=j=0; i<m; i++, j+=n)
        (*this)[i] = &v[j];
}

template<class T>
Mat<T>::Mat(int m, int n_, const T& a)
: std::vector<T*>(m), n(n_), v(m*n, a)
{
    int i,j;
    for(i=j=0; i<m; i++, j+=n)
        (*this)[i] = &v[j];
}

template<class T>
Mat<T>::Mat(const Mat& a)
: std::vector<T*>(a.nrows()), n(a.n), v(a.v)
{
    int i,j;
    for(i=j=0; i<nrows(); i++, j+=n)
        (*this)[i] = &v[j];
}

template<class T>
Mat<T>& Mat<T>::operator=(const Mat<T>& a)
// if this matrix and a have different sizes,
// this matrix is resized to match the size of a
{
    if(this == &a) return *this;
    v = a.v;
    if(nrows() == a.nrows() &&
       ncols() == a.ncols()) return *this;
    int i,j;
    n = a.ncols();
    std::vector<T*>::resize(a.nrows());
    for(i=j=0; i<nrows(); i++, j+=n)
        (*this)[i] = &v[j];
    return *this;
}

template<class T>
Mat<T>& Mat<T>::operator=(const T& a)
//assign a to every element
{
    for(int i=0; i<v.size(); i++)
        v[i] = a;
    return *this;
}

template<class T>
void Mat<T>::resize(int m, int n_)
// data remain unchanged
{
    if(m==nrows() && n_==ncols()) return;
    int i,j;
    n = n_;
    std::vector<T*>::resize(m);
    v.resize(m*n);
    for(i=j=0; i<m; i++, j+=n)
        (*this)[i] = &v[j];
}

template<class T>
void Mat<T>::resize(int m, int n, const T& a)
//assign a to every element
{
    resize(m,n);
    for(int i=0; i<v.size(); i++)
        v[i] = a;
}

typedef Mat<double> mat_double;
typedef Mat<int> mat_int;
typedef std::vector<double> vec_double;
typedef std::vector<int> vec_int;

#endif // __Mat_h__
