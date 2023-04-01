#include<fstream>
#include<sstream>
#include<string>
#include "Mat.h"

void loadtxt(mat_double& A, std::ifstream& f)
// read data from text file f and write to matrix A
// data must be rectangular (m rows and n columns)
// A is automatically resized to shape (m,n)
// if f is not open, A is resized to (0,0)
{
    int m(0),n(0),i;
    char c;
    std::string s;
    std::istringstream ss;
    A.resize(0,0);
    if(!f.is_open()) return;
    while(!f.eof()) {
        getline(f,s);
        ss.str(s);
        if(!isdigit(c = ss.get()) && c!='-' && c!='.')
            continue;// skip line not beginning with number
        ss.unget();
        if(m==0) while(!ss.eof()) {
            if(!isdigit(c = ss.get()) && c!='-' && c!='.')
                continue;
            ss.unget();
            A.resize(1,n+1);
            ss >> A[0][n++];
        }
        else {
            A.resize(m+1,n);
            for(i=0; i<n; i++) ss >> A[m][i];
        }
        ss.clear();
        m++;
    }
}