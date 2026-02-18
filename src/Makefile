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

FLAGS=-Wall -g -fPIC -I$(prefix_amjCom)/include -I$(prefix_amjCom)/include/amjCom -I$(prefix_amjInterferometry)/include
#-g -fsanitize=address,undefined -fno-omit-frame-pointer
CPPFLAGS=$(FLAGS)


SIM_SRC=sim.C amjFourier.C
SIM_OBJ=$(SIM_SRC:.C=.o)

SIM2_SRC=sim2.C amjFourier.C
SIM2_OBJ=$(SIM2_SRC:.C=.o)

VALIDATION_SRC=validation.C amjFourier.C
VALIDATION_OBJ=$(VALIDATION_SRC:.C=.o)

AMJFOURIER_SRC=amjFourier.C
AMJFOURIER_OBJ=$(AMJFOURIER_SRC:.C=.o)

CPP=g++ $(CPPFLAGS)


all: libamjFourier.so sim sim2 validation

sim: $(SIM_OBJ)
	$(CPP) -o $@ $^ $(LDFLAGS) $(LIBS) -lamjArg

sim2: $(SIM2_OBJ)
	$(CPP) -o $@ $^ $(LDFLAGS) $(LIBS) -lamjArg

validation: $(VALIDATION_OBJ)
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
	- rm $(SIM_OBJ) sim $(SIM2_OBJ) sim2 $(VALIDATION_OBJ) validation $(AMJFOURIER_OBJ) \
	libamjFourier.so
