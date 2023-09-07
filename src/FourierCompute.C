#include "../include/amjFourier.H" 

// nL and nF are frame size in wavleength and fringe direction
// periods are the fringe periods at wavelength channel 0.
// dimension of periods is number of phasor sets (baselines) to compute
FourierCompute::FourierCompute(int nL, int nF,
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

void FourierCompute::phasors(const Frame<uint16_t> &frame,PhasorSets &phasors){
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
