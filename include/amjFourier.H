#ifndef _FOURIERSIM_H_
#define _FOURIERSIM_H_

#include "Frame.H"

#include <vector>
#include <tuple>
#include <complex>

class Beam{
public:
  Beam(double x, double D, double (*delay)(double t)=[](double t)->double{return 0;},
       double (illumination)(double L)=[](double L)->double{return 1;}):
    _x(x),_D(D),_delay(delay),_illumination(illumination){};
  ~Beam(){};
  double x(){return _x;}
  double D(){return _D;}
  double delay(double t){return _delay(t);}
  double illumination(double L){return _illumination(L);}
private:
  double _x; // Beam position in mm
  double _D; // Beam diameter in mm
  double (*_delay)(double t); // Beam delay in microns
  double (*_illumination)(double L); // Beam brightness as a function of wavelength in microns
};

class Baseline{
public:
  Baseline(Beam &beam1, Beam &beam2,std::complex<double>(*visibility)(double)
	   =[](double L)->std::complex<double>{return 1;})
    :_beam1(beam1),_beam2(beam2),_visibility(visibility){};
  ~Baseline(){};
  Beam &beam1(){return _beam1;}
  Beam &beam2(){return _beam2;}
  std::complex<double> visibility(double L){return _visibility(L);}
private:
  Beam _beam1,_beam2;
  std::complex<double> (*_visibility)(double);
};

class FourierSim{
public:
  FourierSim(int nL=256, int nF=320, double d=24, double f1=1203, double d2=52.5, double f2=44,
	     double (*wavelength)(int)=[](int i)->double{return 1/(1+(1/2.5-1)/256*i);},
	     double (*bandpass)(int)=[](int i)->double{return (1/(1+(1/2.5-1)/256*(i+1))
							       -1/(1+(1/2.5-1)/256*(i-1)));});
  ~FourierSim();
  void beams(std::vector<Beam> beams){_beams=beams;};
  void baselines(std::vector<Baseline> baselines);
  void dim(int &_nL, int &_nF){_nL=nL;_nF=nF;};
  double **frame(double t);
  template<class T>
  int frame(Frame<T> &) const;
  std::complex<double> **vis();
  // Fringe direction position as a function of pixel number. Pixel 0 is center of left-most pixel,
  // pixel -0.5 is left side of left-most pixel, pixel NF-1 is center of right-most pixel, and
  // pixel NF-0.5 is right side of right-most pixel. etc... returned position is relative to center
  // of row which is pixel number (-0.5+NF-0.5)/2=(NF-1)/2
  double x(double iF) const {return ((double)iF-((double)nF-1)/2)*d;};
  double delay(double iF, double t, Baseline &b) const; // Total delay for a pixel for a baseline
private:
  int nL,nF; // Number of pixels in wavelength and fringe direction
  double d; // Pixel size in wavelength direction in microns
  double f1; // Focal length of M1 in mm
  double d2; // Distance from M1 focus to M2 in mm
  double f2; // Focal length of M2 in mm
  double m2; // Magnification in fringe direction by M2
  double (*wavelength)(int); // Wavelength in microns
  double (*bandpass)(int); // Bandpass in microns
  double sinc(double x) const {if(fabs(x)<1e-6) return 1-x*x/6; return sin(x)/x;};
  double bessj1(double x) const;
  double airy(double x, double L, Beam &b) const;
  double envelope(double delay, int iL) const;
  double ** _frame; // [nlambda][nfringe]
  std::complex<double> ** _vis; // [nBaselines][nLambda]
  std::vector<Beam> _beams;
  std::vector<Baseline> _baselines;
};

#endif /* _FOURIERSIM_H_ */