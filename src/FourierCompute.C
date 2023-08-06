#include "../include/FourierCompute.H" 

template<class T, class U>
void FourierCompute::FourierCompute(double d=24, double f1=1203, double d2=52.5, double f2=44,
				    double (*wavelength)):d(d),f1(f1),d2(d2),f2(f2),
							  wavelength(wavelength){}

template<class T, class U>
void FourierCompute::coherentFlux(Frame<T> &frame,
				  std::vector<Baseline> &baselines,
				  std::vector<std::vector<std::complex<U> > > &cf){
  for(int iB=0;iB<baselines.size();iB++)
    for(iL=0;iL<frame.nL();iL++){
      cf[iB][iL]=0;
      double xwavelength=wavelength(iL)*f1*m2/(baselines[iB].beam1.x()-baselines[iB].beam2.x());
      double N=4*xwavelength/d;
      int iF1=nF/2-N;
      if(iF1<0) iF1=0;
      int iF2=nF/2+N;
      if(iF2>NF-1) iF2=nF-1;
      for(int iF=iF1;iF<iF2;iF++)
	cf[iB][iL]+=frame[iL][iF]*std::polar(1,2*M_PI*x(iF)/xwavelength);
    }
}
