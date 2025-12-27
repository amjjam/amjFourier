#include <amjArg.H>

amjArg::Help

help({"sim",
    "--illumination A B C - Illumination for the three beams. Default 1",
    "--visibility A B C - visibilities for the three baselines. Default 1",
    "--window L0 F0 wL wF first corner, width, each in wavelength and fringe direction"});

#include "../include/amjFourier.H"

#include <stdio.h>
#include <math.h>
#include <vector>


std::vector<double> illumination({1,1,1});
std::vector<double> visibility({1,1,1});

void parse_args(int argc, char *argv[]);

int main(int argc, char *argv[]){
  parse_args(argc,argv);

  std::vector<amjFourier::Beam> beams{
    amjFourier::Beam(0,13,[](double t)->double{return 0;},
		     [](double L)->double{return illumination[0];}),
    amjFourier::Beam(26,13,[](double t)->double{return 2;},
		     [](double L)->double{return illumination[1];}),
    amjFourier::Beam(78,13,[](double t)->double{return -2;},
		     [](double L)->double{return illumination[2];})};
  
  std::vector<amjFourier::Baseline> baselines{
    amjFourier::Baseline("baseline1",beams[0],beams[1],
			 [](double L)->std::complex<double>{return visibility[0];}),
    amjFourier::Baseline("baseline2",beams[1],beams[2],
			 [](double L)->std::complex<double>{return visibility[1];}),
    amjFourier::Baseline("baseline3",beams[0],beams[2],
			 [](double L)->std::complex<double>{return visibility[2];})};

  amjFourier::Sim sim(beams,baselines);
  amjFourier::Frame<double> frame;
  sim(0,frame);
  
  //double nn;
  //std::vector<double> n;
  //std::vector<double> nv;
  //std::vector<double> nv2;
  
  //f.frame(0,frame,nn,n,nv,nv2);
  
  //double **frame=f.frame(0);
  
  int nL=frame.nL(),nF=frame.nF();
  //f.dim(nL,nF);
  
  FILE *fp=fopen("frame.dat","w");
  fwrite(&nL,sizeof(int),1,fp);
  fwrite(&nF,sizeof(int),1,fp);
  for(int iL=0;iL<nL;iL++)
    fwrite(frame[iL].data(),sizeof(double),nF,fp);
  fclose(fp);
  
  return 0;
}

void parse_args(int argc, char *argv[]){

}
