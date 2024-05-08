#include "../include/amjFourier.H"
#include "../include/Frame.H"

#include <stdio.h>
#include <math.h>

int main(int argc, char *argv[]){
  Beam beam1(0,13,[](double t)->double{return 0;},[](double L)->double{return 1/*exp(-2*L+2)*/;});
  Beam beam2(26,13,[](double t)->double{return 2;},[](double L)->double{return 1/*exp(-2*L+2)*/;});
  Beam beam3(78,13,[](double t)->double{return -2;},[](double L)->double{return 1/*exp(-2*L+2)*/;});
  std::vector<Beam> beams={beam1,beam2,beam3};
  
  Baseline baseline1("baseline1",beam1,beam2,[](double L)->std::complex<double>{return 1/*exp(-2*L+2)*/;});
  Baseline baseline2("baseline2",beam2,beam3,[](double L)->std::complex<double>{return 1/*exp(-2*L+2)*/;});
  Baseline baseline3("baseline3",beam1,beam3,[](double L)->std::complex<double>{return 1/*exp(-2*L+2)*/;});
  std::vector<Baseline> baselines={baseline1,baseline2,baseline3};

  FourierSim f(beams,baselines);

  Frame<double> frame(256,320);
  std::vector<double> nv2;
  
  f.frame(0,frame,nv2);
  
  //double **frame=f.frame(0);
  
  int nL=frame.nL(),nF=frame.nF();
  //f.dim(nL,nF);
  
  FILE *fp=fopen("frame.dat","w");
  fwrite(&nL,sizeof(int),1,fp);
  fwrite(&nF,sizeof(int),1,fp);
  for(int iL=0;iL<nL;iL++)
    fwrite(frame[iL].data(),sizeof(double),nF,fp);
  fclose(fp);
  
  // std::complex<double> **vis=f.vis();
  
  // fp=fopen("vis.dat","w");
  // int nB=baselines.size();
  // fwrite(&nB,sizeof(int),1,fp);
  // fwrite(&nL,sizeof(int),1,fp);
  // double x,y;
  // for(int iB=0;iB<nB;iB++)
  //   for(int iL=0;iL<nL;iL++){
  //     x=vis[iB][iL].real();
  //     y=vis[iB][iL].imag();
  //     fwrite(&x,sizeof(double),1,fp);
  //     fwrite(&y,sizeof(double),1,fp);
  //   }
  // fclose(fp);
      
  return 0;
}

