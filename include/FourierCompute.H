#ifndef _FOURIERCOMPUTE_H_
#define _FOURIERCOMPUTE_H_

template<class T, class U>
class FourierCompute{
public:
  FourierCompute(double d=24, double f1=1203, double d2=52.5, double f2=44,
		 double (*wavelength)(int)=
		 [](int i)->double{return 1/(1+(1/2.5-1)/256*i);});
  ~FourierCompute();
  void coherentflux(Frame<T> &frame, std::vector<Baseline> &baselines,
		    std::vector<std::vector<std::complex<U> > > &cf);
private:
  
};

#endif /* _FOURIERCOMPUTE_H_ */

