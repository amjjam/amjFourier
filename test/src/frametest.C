
#include <time.h>
#include <iostream>

#include "../../include/amjFourier.H"

double tdiff(timespec t1, timespec t2);

//int nL=64,nF=64;
int nL=256,nF=320;

int main(int argc, char *argv[]){
  
  std::vector<double> positions(6,0);
  
  // Initialize the frame simulator
  std::vector<Beam> beams{
    Beam(0,13,[&positions](double t)->double{return positions[0];},
	 [](double L)->double{return 20;}),
    Beam(26,13,[&positions](double t)->double{return positions[1];},
	 [](double L)->double{return 20;}),
    Beam(78,13,[&positions](double t)->double{return positions[2];},
	 [](double L)->double{return 20;})};
  
  std::vector<Baseline> baselines{
    Baseline(beams[0],beams[1],
	     [](double L)->std::complex<double>{return 1;}),
    Baseline(beams[1],beams[2],
	     [](double L)->std::complex<double>{return 1;}),
    Baseline(beams[2],beams[0],
	     [](double L)->std::complex<double>{return 1;})};
  
  FourierSim f(beams,baselines,nL,nF);

  Frame<double> framed(nL,nF);
  Frame<uint16_t> frameu(nL,nF);

  long seed=1;

  timespec t0,t1,t2;
  
  clock_gettime(CLOCK_MONOTONIC,&t0);

  for(int i=0;i<1000;i++)
    f.frame(0,framed);

  clock_gettime(CLOCK_MONOTONIC,&t1);

  for(int i=0;i<1000;i++)
    poisson<double,uint16_t>(framed,frameu,seed);

  clock_gettime(CLOCK_MONOTONIC,&t2);

  std::cout << tdiff(t0,t1) << std::endl;
  std::cout << tdiff(t1,t2) << std::endl;
  
  return 0;
}

double tdiff(timespec t1, timespec t2){
  return ((t2.tv_sec-t1.tv_sec)*1e9+(t2.tv_nsec-t1.tv_nsec))/1e9;
}
