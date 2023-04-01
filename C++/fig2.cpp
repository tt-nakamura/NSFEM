#include<iostream>
#include<fstream>
#include "NSFEM.h"
using namespace NSFEM;

main() {
    int i;
    double x,y;
    std::ofstream f;

    init("../obstacle/node.txt",
         "../obstacle/element.txt",
         "../obstacle/stream_fix.txt",
         "../obstacle/vortex_fix.txt",
         "../obstacle/vortex_wall.txt");
    i = run(0.2, 0.05, 0.01);
    std::cout << i << std::endl;

    f.open("fig2uv.txt");
    for(i=0; i<element.nrows(); i++) {
        elem_center(x, y, i);
        f << x << ' ';
        f << y << ' ';
        f << velocity[i][0] << ' ';
        f << velocity[i][1] << '\n';
    }
    f.close();

    f.open("fig2ps.txt");
    for(i=0; i<node.nrows(); i++) {
        f << node[i][0] << ' ';
        f << node[i][1] << ' ';
        f << stream[i] << '\n';
    }
    f.close();

    f.open("fig2om.txt");
    for(i=0; i<node.nrows(); i++) {
        f << node[i][0] << ' ';
        f << node[i][1] << ' ';
        f << vortex[i] << '\n';
    }
    f.close();
}
