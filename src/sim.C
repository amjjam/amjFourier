#include <amjArg.H>

amjArg::Help

help({"sim",
    "--illumination A B C - Illumination for the three beams. Default 1",
    "--visibility A B C - visibilities for the three baselines. Default 1",
    "--window L0 F0 wL wF first corner, width, each in wavelength and fringe direction",
    "--diameter A B C - beam diameters in mm (default 13 mm)"});

#include "../include/amjFourier.H"

#include <stdio.h>
#include <math.h>
#include <vector>
#include <cstdint>

std::vector<double> illumination({1,1,1});
std::vector<double> visibility({1,1,1});
std::vector<double> diameters({13,13,13});

void parse_args(int argc, char *argv[]);
unsigned int L0=0,F0=0,wL=256,wF=320;

int main(int argc, char *argv[]){
  parse_args(argc,argv);

  std::vector<amjFourier::Beam> beams{
    amjFourier::Beam(0,diameters[0],[](double t)->double{return 0;},
		     [](double L)->double{return illumination[0];}),
    amjFourier::Beam(26,diameters[1],[](double t)->double{return 2;},
		     [](double L)->double{return illumination[1];}),
    amjFourier::Beam(78,diameters[2],[](double t)->double{return -2;},
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
  sim.set_window(L0,F0,wL,wF);
  sim(frame);
  
  int nL=frame.nL(),nF=frame.nF();
  
  FILE *fp=fopen("frame.dat","w");
  std::uint16_t ui16=0;
  std::uint8_t ui8=0;
  std::uint32_t ui32=0;
  fwrite(&ui16,sizeof(std::uint16_t),1,fp);
  for(int i=0;i<5;i++)
    fwrite(&ui8,sizeof(std::uint8_t),1,fp);
  fwrite(&ui32,sizeof(std::uint32_t),1,fp);
  fwrite(&nL,sizeof(int),1,fp);
  fwrite(&nF,sizeof(int),1,fp);
  for(int iL=0;iL<nL;iL++)
    fwrite(frame[iL].data(),sizeof(double),nF,fp);
  fclose(fp);
  
  return 0;
}

void parse_args(int argc, char *argv[]){
  help(argc,argv);
  for(int i=1;i<argc;i++)
    if(strcmp(argv[i],"--illumination")==0){
      for(int j=0;j<3;j++){
	i++;
	illumination[j]=atof(argv[i]);
      }
    }
    else if(strcmp(argv[i],"--visibility")==0){
      for(int j=0;j<3;j++){
	i++;
	visibility[j]=atof(argv[i]);
      }
    }
    else if(strcmp(argv[i],"--window")==0){
      i++;
      L0=atoi(argv[i]);
      i++;
      F0=atoi(argv[i]);
      i++;
      wL=atoi(argv[i]);
      i++;
      wF=atoi(argv[i]);
    }
    else if(strcmp(argv[i],"--diameters")==0){
      for(int j=0;j<3;j++){
	i++;
	diameters[j]=atof(argv[i]);
      }
    }
    else{
      std::cout << "unknow parameter: " << argv[i] << std::endl;
      exit(1);
    }
}
