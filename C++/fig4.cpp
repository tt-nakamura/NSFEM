#include<iostream>
#include<fstream>
#include "NSFEM.h"
using namespace NSFEM;

main() {
    int i;
    double x,y;
    std::ofstream f;

    init("../cylinder/node.txt",
         "../cylinder/element.txt",
         "../cylinder/stream_fix.txt",
         "../cylinder/vortex_fix.txt",
         "../cylinder/vortex_wall.txt");
    // Karman vortex appears behind the cylinder
    i = run(0.5, 0.001, 0.003, 440);
    std::cout << i << std::endl;

    f.open("fig4uv.txt");
    for(i=0; i<element.nrows(); i++) {
        elem_center(x, y, i);
        f << x << ' ';
        f << y << ' ';
        f << velocity[i][0] << ' ';
        f << velocity[i][1] << '\n';
    }
    f.close();

    f.open("fig4ps.txt");
    for(i=0; i<node.nrows(); i++) {
        f << node[i][0] << ' ';
        f << node[i][1] << ' ';
        f << stream[i] << '\n';
    }
    f.close();

    f.open("fig4om.txt");
    for(i=0; i<node.nrows(); i++) {
        f << node[i][0] << ' ';
        f << node[i][1] << ' ';
        f << vortex[i] << '\n';
    }
    f.close();
}
