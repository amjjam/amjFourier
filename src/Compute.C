#include "../include/amjFourier.H" 

namespace amjFourier{
// nL and nF are frame size in wavleength and fringe direction
// periods are the fringe periods at wavelength channel 0.
// dimension of periods is number of phasor sets (baselines) to compute
  Phasors::Phasors(int nL, int nF,
		   const std::vector<float> &periods,
		   std::function<double(int)> wavelength)
    :nL(nL),nF(nF),wavelength(wavelength),periods(periods){
    c.resize(periods.size()*nL*nF);
    s.resize(periods.size()*nL*nF);
    int i1,i2,i3;
    float period,f;
    for(unsigned int iB=0;iB<periods.size();iB++){
      i1=iB*nL*nF;
      for(int iL=0;iL<nL;iL++){
	i2=i1+iL*nF;
	period=periods[iB]/wavelength(0)*wavelength(iL);
	for(int iF=0;iF<nF;iF++){
	  i3=i2+iF;
	  f=2*M_PI*(iF-(float)(nF-1)/2)/period;
	  c[i3]=cos(f);
	  s[i3]=sin(f);
	}
      }
    }
  }
  
  // 
  void Phasors::operator()(const Frame<double> &frame,PhasorSets &phasors){
    phasors.resize(periods.size());
    int i1,iL,i2,iF,i3;
    float real,imag;
    for(unsigned int iB=0;iB<periods.size();iB++){
      i1=iB*nL*nF;
      phasors[iB].resize(frame.nL());
      for(iL=0;iL<nL;iL++){
	i2=i1+iL*nF;
	real=0;
	imag=0;
	for(iF=0;iF<nF;iF++){
	  i3=i2+iF;
	  real+=frame[iL][iF]*c[i3];
	  imag+=frame[iL][iF]*s[i3];
	}
	phasors[iB][iL].real(real);
	phasors[iB][iL].imag(imag);
      }
    }
  }
  
  double wavelength(int i){ // returns wavelength in um for channel i
    static const double c[5]={6.95694016e-06,-5.20602460e-03, 1.46461339e+00,
			      -1.96798243e+02,1.28465801e+04};
    double v=c[4];
    double p=i;
    for(int j=3;j>=0;j--){
      v+=c[j]*p;
      p*=i;
    }
    return v;
  }
}
