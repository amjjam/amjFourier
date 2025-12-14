#include "../include/amjFourier.H"
#include "../include/Frame.H"

#include <stdio.h>
#include <math.h>

int main(int argc, char *argv[]){
  amjFourier::Beam beam1(0,13,[](double t)->double{return 0;},[](double L)->double{return 1/*exp(-2*L+2)*/;});
  amjFourier::Beam beam2(26,13,[](double t)->double{return 2;},[](double L)->double{return 1/*exp(-2*L+2)*/;});
  amjFourier::Beam beam3(78,13,[](double t)->double{return -2;},[](double L)->double{return 1/*exp(-2*L+2)*/;});
  std::vector<amjFourier::Beam> beams={beam1,beam2,beam3};
  
  amjFourier::Baseline baseline1("baseline1",beam1,beam2,[](double L)->std::complex<double>{return 1/*exp(-2*L+2)*/;});
  amjFourier::Baseline baseline2("baseline2",beam2,beam3,[](double L)->std::complex<double>{return 1/*exp(-2*L+2)*/;});
  amjFourier::Baseline baseline3("baseline3",beam1,beam3,[](double L)->std::complex<double>{return 1/*exp(-2*L+2)*/;});
  std::vector<amjFourier::Baseline> baselines={baseline1,baseline2,baseline3};

  amjFourier::Sim f(beams,baselines);
  
  amjFourier::Frame<double> frame(256,320);
  double nn;
  std::vector<double> n;
  std::vector<double> nv;
  std::vector<double> nv2;
  
  f.frame(0,frame,nn,n,nv,nv2);
  
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

