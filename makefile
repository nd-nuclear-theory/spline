CC = g++
CFLAGS = -std=c++11 -g -c -Wall
LIBS = -lgsl -lgslcblas -lm
LDIR = -L"C:\MingW\lib"
LIBH = -I"C:\MingW\include\gsl" 
Eigen = -I"C:\MingW\include\Eigen"
OBJ = h2utils_input.o spline.o wavefunction_basis.o wavefunction_class.o

h2_input: $(OBJ)
	$(CC) -std=c++11 -o h2_input $(OBJ) $(LDIR) $(LIBS)

h2utils_input.o: h2utils_input.cpp
	$(CC) $(CFLAGS) h2utils_input.cpp $(Eigen)

spline.o: spline.cpp spline.h
	$(CC) $(CFLAGS) spline.cpp $(LIBH)

wavefunction_basis.o: wavefunction_basis.cpp wavefunction_basis.h
	$(CC) $(CFLAGS) wavefunction_basis.cpp

wavefunction_class.o: wavefunction_class.cpp wavefunction_class.h
	$(CC) $(CFLAGS) wavefunction_class.cpp $(Eigen)

.PHONY : clean
clean :
	-rm edit $(OBJ)