#include<iostream>
#include<fstream>
#include "NSFEM.h"
using namespace NSFEM;

main() {
    int i;
    double x,y;
    std::ofstream f;

    init("../cavity/node.txt",
         "../cavity/element.txt",
         "../cavity/stream_fix.txt",
         "../cavity/vortex_fix.txt",
         "../cavity/vortex_wall.txt");
    i = run(0.1, 0.01, 0.01);
    std::cout << i << std::endl;

    f.open("fig6uv.txt");
    for(i=0; i<element.nrows(); i++) {
        elem_center(x, y, i);
        f << x << ' ';
        f << y << ' ';
        f << velocity[i][0] << ' ';
        f << velocity[i][1] << '\n';
    }
    f.close();

    f.open("fig6ps.txt");
    for(i=0; i<node.nrows(); i++) {
        f << node[i][0] << ' ';
        f << node[i][1] << ' ';
        f << stream[i] << '\n';
    }
    f.close();

    f.open("fig6om.txt");
    for(i=0; i<node.nrows(); i++) {
        f << node[i][0] << ' ';
        f << node[i][1] << ' ';
        f << vortex[i] << '\n';
    }
    f.close();
}
