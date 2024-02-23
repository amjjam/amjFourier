ifndef prefix
	prefix=/opt/amj
endif
ifndef prefix_amjCom
	prefix_amjCom=$(prefix)
endif

ifndef prefix_amjInterferometry
	prefix_amjInterferometry=$(prefix)
endif

LIBS=-lamjCom -lamjInterferometry
LDFLAGS=-L$(prefix_amjCom)/lib -L$(prefix_amjInterferometry)/lib

FLAGS=-Wall -g -O3 -fPIC -I$(prefix_amjCom)/include -I$(prefix_amjInterferometry)/include
CPPFLAGS=$(FLAGS)


SIM_SRC=main.C FourierSim.C
SIM_OBJ=$(SIM_SRC:.C=.o)

AMJFOURIER_SRC=FourierSim.C Frame.C FourierCompute.C
AMJFOURIER_OBJ=$(AMJFOURIER_SRC:.C=.o)

CPP=g++ $(CPPFLAGS)


all: libamjFourier.so sim

sim: $(SIM_OBJ)
	$(CPP) -o $@ $^ $(LDFLAGS) $(LIBS)

libamjFourier.so: $(AMJFOURIER_OBJ)
	$(CPP) -shared -Wl,-soname,$@ -o $@ $^ $(LDFLAGS) $(LIBS)

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

