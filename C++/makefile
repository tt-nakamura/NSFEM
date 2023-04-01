OBJ = NSFEM.o MatBand.o MatBandSym.o loadtxt.o
GSL = -lgsl -lcblas -L/usr/local/lib

fig2: fig2.o $(OBJ)
	g++ fig2.o $(OBJ) $(GSL)
fig4: fig4.o $(OBJ)
	g++ fig4.o $(OBJ) $(GSL)
fig6: fig6.o $(OBJ)
	g++ fig6.o $(OBJ) $(GSL)
