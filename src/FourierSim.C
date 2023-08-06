#include "../include/amjFourier.H"

FourierSim::FourierSim(int nL, int nF, double d, double f1, double d2, double f2,
			  double (*wavelength)(int),
			  double (*bandpass)(int)):
  nL(nL),nF(nF),d(d),f1(f1),d2(d2),f2(f2),wavelength(wavelength),bandpass(bandpass){
  _frame=(double**)malloc(nL*sizeof(double*));
  _frame[0]=(double*)malloc(nL*nF*sizeof(double));
  for(int i=1;i<nL;i++)
    _frame[i]=_frame[0]+i*nF;
  m2=f2/(f2-d2); // Magnification in fringe direction by M2 (5.18)
}

FourierSim::~FourierSim(){
  free(_frame[0]);
  free(_frame);
  /* Also free _vis, only if it was allocated */
}

void FourierSim::baselines(std::vector<Baseline> baselines){
  _baselines=baselines;
  _vis=(std::complex<double> **)malloc(_baselines.size()*sizeof(std::complex<double> *));
  _vis[0]=(std::complex<double> *)malloc(_baselines.size()*nL*sizeof(std::complex<double>));
  for(int i=1;i<nL;i++)
    _vis[i]=_vis[0]+i*nL;
}

#include<iostream>
double **FourierSim::frame(double t){
  // Clear the frame
  for(int iL=0;iL<nL;iL++)
    for(int iF=0;iF<nF;iF++)
      _frame[iL][iF]=0;
  // Frame illumination
  for(int iB=0;iB<_beams.size();iB++)
    for(int iL=0;iL<nL;iL++)
      for(int iF=0;iF<nF;iF++)
	_frame[iL][iF]+=airy(x(iF),wavelength(iL),_beams[iB]);
  
  // Fringes
  for(int iB=0;iB<_baselines.size();iB++)
    for(int iL=0;iL<nL;iL++){
      std::complex<double> V0=_baselines[iB].visibility(wavelength(iL));
      for(int iF=0;iF<nF;iF++){
	double I1=airy(x(iF),wavelength(iL),_baselines[iB].beam1());
	double I2=airy(x(iF),wavelength(iL),_baselines[iB].beam2());
	std::complex V=V0*envelope(delay(iF,t,_baselines[iB]),iL);
	//std::cout << iB << " " << iL << " " << iF << " " << delay(iF,t,_baselines[iB]) << " " << envelope(delay(iF,t,_baselines[iB]),iL) << " " << I1 << " " << I2 << std::endl;
	_frame[iL][iF]+=2*sqrt(I1*I2)*fabs(V)*cos(2*M_PI*delay(iF,t,_baselines[iB])/wavelength(iL)
						  +std::arg(V));
      }
    }
  return _frame;
}

template<class T>
int FourierSim::frame(Frame<T> &frame) const{
  // Clear the frame
  for(int iL=0;iL<nL;iL++)
    for(int iF=0;iF<nF;iF++)
      frame[iL][iF]=0;
  // Frame illumination
  for(int iB=0;iB<_beams.size();iB++)
    for(int iL=0;iL<nL;iL++)
      for(int iF=0;iF<nF;iF++)
	frame[iL][iF]+=airy(x(iF),wavelength(iL),_beams[iB]);
  
  // Fringes
  for(int iB=0;iB<_baselines.size();iB++)
    for(int iL=0;iL<nL;iL++){
      std::complex<double> V0=_baselines[iB].visibility(wavelength(iL));
      for(int iF=0;iF<nF;iF++){
	double I1=airy(x(iF),wavelength(iL),_baselines[iB].beam1());
	double I2=airy(x(iF),wavelength(iL),_baselines[iB].beam2());
	std::complex V=V0*envelope(delay(iF,t,_baselines[iB]),iL);
	frame[iL][iF]+=2*sqrt(I1*I2)*fabs(V)*cos(2*M_PI*delay(iF,t,_baselines[iB])/wavelength(iL)
						 +std::arg(V));
      }
    }
  return 0;      
}

// std::complex<double> **FourierSim::vis(){
//   for(int iB=0;iB<_baselines.size();iB++){
//     std::cout << iB << std::endl;
//     for(int iL=0;iL<nL;iL++){
//       double xx=0,yy=0,nn=0;
//       double xwavelength=wavelength(iL)*f1*m2/(_baselines[iB].beam1().x()-_baselines[iB].beam2().x());
//       double N=4*xwavelength/d;
//       int iF1=nF/2-N;
//       if(iF1<0) iF1=0;
//       int iF2=nF/2+N;
//       if(iF2>nF-1) iF2=nF-1;
//       for(int iF=iF1;iF<iF2;iF++){
// 	xx+=_frame[iL][iF]*cos(2*M_PI*x(iF)/xwavelength);
// 	yy+=_frame[iL][iF]*sin(2*M_PI*x(iF)/xwavelength);
// 	nn+=_frame[iL][iF];
//       }
//       _vis[iB][iL].real(xx);
//       _vis[iB][iL].imag(yy);
//     }
//   }
//   return _vis;
// }

double FourierSim::delay(double iF, double t, Baseline &b){
  return (b.beam1().x()-b.beam2().x())*x(iF)/f1/m2+b.beam2().delay(t)-b.beam1().delay(t);
}

/* Borrowed from Press et al. Numerical Recipes */
double FourierSim::bessj1(double x){
  double ax,z;
  double xx,y,ans,ans1,ans2;
  
  if ((ax=fabs(x)) < 8.0) {
    y=x*x;
    ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
					      +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
    ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
					   +y*(99447.43394+y*(376.9991397+y*1.0))));
    ans=ans1/ans2;
  } else {
    z=8.0/ax;
    y=z*z;
    xx=ax-2.356194491;
    ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
			       +y*(0.2457520174e-5+y*(-0.240337019e-6))));
    ans2=0.04687499995+y*(-0.2002690873e-3
			  +y*(0.8449199096e-5+y*(-0.88228987e-6
						 +y*0.105787412e-6)));
    ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
    if (x < 0.0) ans = -ans;
  }
  return ans;
}


double FourierSim::airy(double x, double L, Beam &b){
  double Ic=b.illumination(L);
  double y=M_PI*x*b.D()/f1/m2/L;
  if(fabs(y)<1e-6)
    return Ic;
  double tmp=2*bessj1(y)/y;
  return Ic*tmp*tmp;
}

double FourierSim::envelope(double delay, int iL){
  double y=M_PI*delay*bandpass(iL)/wavelength(iL)/wavelength(iL);
  return sinc(y);
}

