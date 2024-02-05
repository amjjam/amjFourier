FLAGS=-Wall -g -O3 -fPIC -I/opt/amjCom/include -L/opt/amjCom/lib -I/opt/amj/include -L/opt/amj/lib
CPPFLAGS=$(FLAGS)


SIM_SRC=main.C FourierSim.C
SIM_OBJ=$(SIM_SRC:.C=.o)

AMJFOURIER_SRC=FourierSim.C Frame.C FourierCompute.C
AMJFOURIER_OBJ=$(AMJFOURIER_SRC:.C=.o)

CPP=g++ $(CPPFLAGS)

prefix=/opt/amj

all: libamjFourier.so sim

sim: $(SIM_OBJ)
	$(CPP) -o $@ $^ -lamjCom

libamjFourier.so: $(AMJFOURIER_OBJ)
	$(CPP) -shared -Wl,-soname,$@ -o $@ $^ -lamjCom -lamjInterferometry

install: libamjFourier.so ../include/amjFourier.H
	mkdir -p $(prefix)/lib
	mkdir -p $(prefix)/include
	cp libamjFourier.so $(prefix)/lib
	cp ../include/amjFourier.H $(prefix)/include

uninstall:
	- rm $(prefix)/lib/libamjFourier.so
	- rm $(prefix)/include/amjFourier.H
	- rmdir $(prefix)/lib
	- rmdir $(prefix)/include

clean:
	- rm $(SIM_OBJ) sim $(AMJFOURIER_OBJ) libamjFourier.so

