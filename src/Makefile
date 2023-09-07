FLAGS=-Wall -g -fPIC -I/opt/amjCom/include -L/opt/amjCom/lib -I/opt/amj/include -L/opt/amj/lib
CPPFLAGS=$(FLAGS)


SIM_SRC=main.C FourierSim.C
SIM_OBJ=$(SIM_SRC:.C=.o)

AMJFOURIER_SRC=FourierSim.C Frame.C FourierCompute.C
AMJFOURIER_OBJ=$(AMJFOURIER_SRC:.C=.o)

CPP=g++ $(CPPFLAGS)

all: libamjFourier.so sim

sim: $(SIM_OBJ)
	$(CPP) -o $@ $^ -lamjCom

libamjFourier.so: $(AMJFOURIER_OBJ)
	$(CPP) -shared -Wl,-soname,$@ -o $@ $^ -lamjCom -lamjInterferometry

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

