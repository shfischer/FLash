HEADERS=-I../../../inst/include
CXXFLAGS = -std=c++11

DRVOBJS=drivers.o driversf.o odedrivers.o odedriversf.o taylor.o

%.o: %.c
	@echo CC  $<
	@${CC} -fPIC ${HEADERS} ${ADOLC_CFLAGS} -c $<

library: $(DRVOBJS)

clean:
	@-rm -f *.o
