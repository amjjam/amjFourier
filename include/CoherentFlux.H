#ifndef _COHERENTFLUX_H_
#define _COHERENTFLUX_H_

#include <vector>
#include <complex>

template<class T>
class BaselineCoherentFlux{
public:
  BaselineCoherentFlux(int nL){_d.resize(nL);};
  ~BaselineCoherentFlux();
  void rotate(double a);
  int write(char *) const;
  int read(const char *);
  const std::vector<std::complex<T> >& data(){return _d;};
private:
  std::vector<std::complex<T> > _d;
};

template<class T>
class FrameCoherentFlux{
public:
  FrameCoherentFlux(int nB, int nL);
  ~FrameCoherentFlux();
  int write(char *) const;
  int read(const char *);
private:
  std::vector<BaselineCoherentFlux<T> > _d;
};



template<class T>
FrameCoherentFlux::FrameCoherentFlux(int nB, int nL){
  _d.resize(nB);
  for(int i=0;i<

#endif /* _COHERENTFLUX_H_ */
