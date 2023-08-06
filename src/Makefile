CPPFLAGS=-g
CC=g++


SIM_SRC=main.C FourierSim.C
SIM_OBJ=$(SIM_SRC:.C=.o)

AMJFOURIER_SRC=FourierSim.C
AMJFOURIER_OBJ=$(AMJFOURIER_SRC:.C=.o)

all: libamjFourier.so sim

sim: $(SIM_OBJ)
	$(CC) -o $@ $^

libamjFourier.so: $(AMJFOURIER_OBJ)
	$(CPP) -shared -Wl,-soname,$@ -o $@ $^

install: libamjFourier.so ../include/amjFourier.H
	mkdir -p /opt/amj/lib
	mkdir -p /opt/amj/include
	cp libamjFourier.so /opt/amj/lib
	cp ../include/amjFourier.H /opt/amj/include

uninstall:
	- rm /opt/amj/lib/libamjFourier.so
	- rm /opt/amj/include/amjFourier.H

clean:
	- rm $(SIM_OBJ) sim $(AMJFOURIER_OBJ) libamjFourier.so

