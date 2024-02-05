FLAGS=-I/opt/amjCom/include -L/opt/amjCom/lib -lamjCom -I/opt/amj/include \
	-L/opt/amj/lib -lamjInterferometry
CPPFLAGS=$(FLAGS)
CPP=g++ $(CPPFLAGS)

frametest: frametest.o ../../src/FourierSim.o
	$(CPP) -o frametest frametest.o ../../src/FourierSim.o -lamjCom -lamjInterferometry -lamjFourier

