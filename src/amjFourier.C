#include "../include/amjFourier.H"
#include<iostream>

namespace amjFourier{
  Sim::Sim(const std::vector<Beam> &beams,
	   const std::vector<Baseline> &baselines):
    beams(beams),baselines(baselines),
    nL(NL),nF(NF),L0(0),F0(0),wL(NL),wF(NF),max_wavelength(2.2),
    wavelength(amjFourier::wavelength),bandpass(amjFourier::bandpass),
    center(amjFourier::center),d(24),f1(1203),d2(52.5),f2(44){
    m2=f2/(f2-d2); // Magnification in fringe direction by M2 (5.18)
    initialize();
  }

  Sim::~Sim(){
  }

  void Sim::set_wavelength(std::function<double(double)> w){
    wavelength=w;
    initialize();
  }

  void Sim::set_bandpass(std::function<double(double)> b){
    bandpass=b;
    initialize();
  }

  void Sim::set_center(std::function<double(double)> c){
    center=c;
    initialize();
  }

  void Sim::set_window(unsigned int L0_, unsigned int F0_, unsigned int wL_,
		       unsigned int wF_){
    L0=L0_; F0=F0_; wL=wL_; wF=wF_;
    initialize();
  }

  // void Sim::set_bias(std::function<double(int,int)> b){
  //   bias=b;
  //   initialize();
  // }
  
  void Sim::initialize(){
    //std::cout << "wL=" << wL << ", wF=" << wF << std::endl;
    // Compute airys for beams
    airys.resize(beams.size());
    unsigned int iB;
    unsigned int iL,iF; // pixel location in the window
    unsigned int jL,jF; // pixel location in the frame
    int i;
    double yyy,Ic,yy,y,L,tmp;
    for(iB=0;iB<beams.size();iB++){
      airys[iB].resize(wL*wF);
      yyy=M_PI*beams[iB].D()/f1/m2;
      //std::cout << "Beam " << iB << ":" << std::endl;
      for(iL=0,jL=L0;iL<wL;iL++,jL++){
	i=iL*wF;
	L=wavelength(jL);
	//std::cout << "L(" << jL << ")=" << L << std::endl;
	Ic=beams[iB].illumination(L);
	yy=yyy/L;
	for(iF=0,jF=F0;iF<wF;iF++,jF++){
	  if(L>max_wavelength)
	    airys[iB][i+iF]=0;
	  else{
	    y=yy*x(jL,jF);
	    if(fabs(y)<1e-6)
	      airys[iB][i+iF]=Ic;
	    else{
	      tmp=2*bessj1(y)/y;
	      airys[iB][i+iF]=Ic*tmp*tmp;
	    }
	  }
	}
      }
      //std::cout << std::endl << std::endl;
    }
    
    // Compute total illumination as sum of beam airys
    illumination.resize(wL*wF);
    for(iL=0;iL<wL;iL++){
      i=iL*wF;      
      for(iF=0;iF<wF;iF++)
	illumination[i+iF]=0;
    }
    for(iB=0;iB<beams.size();iB++)
      for(iL=0,jL=0;iL<wL;iL++,jL++){
	i=iL*wF;
	for(iF=0;iF<wF;iF++)
	  illumination[i+iF]+=airys[iB][i+iF];
      }
    
    // Determine beam indexes
    unsigned int iBeam;
    bi1.resize(baselines.size());
    bi2.resize(baselines.size());
    for(iB=0;iB<baselines.size();iB++){
      for(iBeam=0;iBeam<beams.size();iBeam++){
	if(&beams[iBeam]==&baselines[iB].beam1()){
	  bi1[iB]=iBeam;
	  break;
	}
      }
      if(iBeam==beams.size()){
	std::cerr << "Could not find beam1 of baseline " << iB << std::endl;
	abort();
      }
      for(iBeam=0;iBeam<beams.size();iBeam++)
	if(&beams[iBeam]==&baselines[iB].beam2()){
	  bi2[iB]=iBeam;
	  break;
	}
      if(iBeam==beams.size()){
	std::cerr << "Could not find beam2 of baseline " << iB << std::endl;
	abort();
      }
    }
  }
  
  /* Borrowed from Press et al. Numerical Recipes */
  double Sim::bessj1(double x) const{
    double ax,z;
    double xx,y,ans,ans1,ans2;
    
    if ((ax=fabs(x)) < 8.0) {
      y=x*x;
      ans1=x*(72362614232.0
	      +y*(-7895059235.0
		  +y*(242396853.1
		      +y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
      ans2=144725228442.0
	+y*(2300535178.0+y*(18583304.74+y*(99447.43394+y*(376.9991397+y*1.0))));
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

#include <math.h>
#define PI 3.141592654
  
  float poidev(float xm, long *idum){
    float gammln(float xx);
    float ran1(long *idum);
    static float sq,alxm,g,oldm=(-1.0);
    float em,t,y;
    
    if (xm < 12.0) {
      if (xm != oldm) {
	oldm=xm;
	g=exp(-xm);
      }
      em = -1;
      t=1.0;
      do {
	++em;
	t *= ran1(idum);
      } while (t > g);
    } else {
      if (xm != oldm) {
	oldm=xm;
	sq=sqrt(2.0*xm);
	alxm=log(xm);
	g=xm*alxm-gammln(xm+1.0);
      }
      do {
	do {
	  y=tan(PI*ran1(idum));
	  em=sq*y+xm;
	} while (em < 0.0);
	em=floor(em);
	t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
      } while (ran1(idum) > t);
    }
    return em;
  }
#undef PI
  /* (C) Copr. 1986-92 Numerical Recipes Software 1)0. */
  
#include <math.h>
  
  float gammln(float xx){
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
			  24.01409824083091,-1.231739572450155,
			  0.1208650973866179e-2,-0.5395239384953e-5};
    int j;
    
    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for (j=0;j<=5;j++) ser += cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
  }
  /* (C) Copr. 1986-92 Numerical Recipes Software 1)0. */
  
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
  
  float ran1(long *idum){
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;
    
    if (*idum <= 0 || !iy) {
      if (-(*idum) < 1) *idum=1;
      else *idum = -(*idum);
      for (j=NTAB+7;j>=0;j--) {
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	if (j < NTAB) iv[j] = *idum;
      }
      iy=iv[0];
    }
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0) *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
  }

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
  /* (C) Copr. 1986-92 Numerical Recipes Software 1)0. */

  float Sim::envelope(float delay, float L, float B) const{
    float y=M_PI*delay*B/L/L;
    return sinc(y);
  }

  CalcPhasors::CalcPhasors():nL(NL),nF(NF),L0(0),F0(0),wL(nL),wF(nF),
			     C0(0),nC(NF),cAvg(1),
			     wavelength(amjFourier::wavelength),
			     center(amjFourier::center),periods({3.34}){
    initialize();
  }
  
  void CalcPhasors::set_wavelength(std::function<double(double)> w){
    wavelength=w;
    initialize();
  }
  
  void CalcPhasors::set_center(std::function<double(double)> c){
    center=c;
    initialize();
  }
  
  void CalcPhasors::set_channels(unsigned int C0_, unsigned int nC_){
    C0=C0_; nC=nC_;
    initialize();
  }
  
  void CalcPhasors::set_periods(std::vector<double> p){
    periods=p;
    initialize();
  }

  void CalcPhasors::set_channel_sum(unsigned int cAvg_){
    cAvg=cAvg_;
    initialize();
  }  
  
  void CalcPhasors::initialize(){
    c.resize(periods.size()*wL*wF);
    s.resize(periods.size()*wL*wF);
    int i1,i2,i3;
    float period,f;
    for(unsigned int iB=0;iB<periods.size();iB++){
      i1=iB*wL*wF;
      for(unsigned int iL=0,jL=L0;iL<wL;iL++,jL++){
	i2=i1+iL*wF;
	period=periods[iB]/wavelength(NL-1)*wavelength(jL);
	for(unsigned int iF=0,jF=F0;iF<wF;iF++,jF++){
	  i3=i2+iF;
	  f=2*M_PI*(jF-center(jL))/period;
	  c[i3]=cos(f);
	  s[i3]=sin(f);
	}
      }
    }
    assert(c.size() == periods.size() * nL * nF);
assert(s.size() == periods.size() * nL * nF);
  }
  
  double wavelength(double i){ // returns wavelength in um for channel i
    // Coefficients from Paolo e-mail 9 October 2025
    //static const double c[5]={6.95694016e-06,-5.20602460e-03, 1.46461339e+00,
    //			      -1.96798243e+02,1.28465801e+04};
    // Coefficients from Paolo e-mail 19 December 2025
    static const double c[5]={-3.63299018e-06*1e-3,3.46705393e-03*1e-3,
			      -1.17335426e+00*1e-3,1.56363989e+02*1e-3,
			      -4.71516112e+03*1e-3};
    if(i<140)
      return 2.5;
    double v=c[4];
    double p=i;
    for(int j=3;j>=0;j--){
      v+=c[j]*p;
      p*=i;
    }
    return v;
  }

  double wavelength_lin(double i){
    // returns wavelength in channel i with a model which is linear in
    // wavenumber. The wavelength in channel 140 is 2.5 um, and the
    // wavelength in channel 255 is 1 um.
    double dk=(1.0-1.0/2.5)/(140-256);
    double k=1.0/2.5+i*dk;
    return 1/k;      
  }
  
  double bandpass(double i){ // Return bandpass in um for channel i.
    // Current model is that bandpass for channel i is difference in
    // wavelength between channels i-1 and i+1
    return fabs(wavelength(i+1)-wavelength(i-1));
  }

  double center(double i){
    return 125;
  }
}
