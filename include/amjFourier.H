#ifndef _AMJFOURIER_H_
#define _AMJFOURIER_H_

#include <vector>
#include <cstdint>
#include <iostream>

#include <complex>
#include <tuple>
#include <functional>
  

#include <amjCom/amjPacketRW.H>
#include <amjCom/amjPacket.H>

#include <amjInterferometry.H>

namespace amjFourier{
  template<class T>
  class Frame:public amjPacketRW{
  public:
    Frame(int nL=1,int nF=1);
    template<class U>
    Frame(const Frame<U> &);
    ~Frame(){};
    void resize(int nL, int nF);
    unsigned int nL() const {return _nL;}
    unsigned int nF() const {return _nF;}
    void clear();
    template<class U>
    Frame<T> &operator=(const Frame<U> &);
    Frame<T> &operator*=(double);
    template<class U>
    Frame<T> &operator+=(const Frame<U> &);
    std::vector<T> &operator[](int iL){return _frame[iL];};
    const std::vector<T> &operator[](int iL) const {return _frame[iL];};
    virtual void write(amjPacket &p) const;
    int size(){return 2*sizeof(uint32_t)+_nL*_nF*sizeof(T);};
    virtual void write(uint8_t *) const;
    virtual void read(amjPacket &p);
    int size1(){return 2*sizeof(uint32_t);};
    void read1(uint8_t *);
    int size2(){return _nL*_nF*sizeof(T);};
    void read2(uint8_t *);
  private:
    uint32_t _nL,_nF;
    std::vector<std::vector<T> > _frame;
  };
  
  template<class T>
  Frame<T>::Frame(int nL, int nF):_nL(nL),_nF(nF){
    resize(nL,nF);
  }
  
  template<class T>
  template<class U>
  Frame<T>::Frame(const Frame<U> &f){
    operator=(f);
  }
  
  template<class T>
  void Frame<T>::resize(int nL, int nF){
    _nL=nL;
    _nF=nF;
    _frame.resize(nL);
    for(int iL=0;iL<nL;iL++)
      _frame[iL].resize(nF);
  }
  
  template<class T>
  void Frame<T>::clear(){
    unsigned int iF;
    for(unsigned int iL=0;iL<_nL;iL++)
      for(iF=0;iF<_nF;iF++)
	_frame[iL][iF]=0;
  }
  
  template<class T>
  template<class U>
  Frame<T> &Frame<T>::operator=(const Frame<U> &f){
    _nL=f.nL();
    _nF=f.nF();
    _frame.resize(_nL);
    for(unsigned int iL=0;iL<_nL;iL++){
      _frame[iL].resize(_nF);
      for(unsigned int iF=0;iF<_nF;iF++)
	_frame[iL][iF]=f[iL][iF];
    }
    return *this;
  }
  
  template<class T>
  Frame<T> &Frame<T>::operator*=(double f){
    for(int iL=0;iL<_nL;iL++)
      for(int iF=0;iF<_nF;iF++)
	_frame[iL][iF]*=f;
    return *this;
  }
  
  template<class T>
  template<class U>
  Frame<T> &Frame<T>::operator+=(const Frame<U> &f){
    if(_nL!=f.nL()||_nF!=f.nF()){
      std::cout << "Frame<T>: frame dimensions do not match: (_nL,_nF)=(" << _nL
		<< "," << _nF << ") and (f.nL,f.nF)=(" << f.nL() << "," << f.nF()
		<< ")" << std::endl;
      abort();
    }
    for(unsigned int iL=0;iL<_nL;iL++)
      for(unsigned int iF=0;iF<_nF;iF++)
	_frame[iL][iF]+=f[iL][iF];
    return *this;
  }
  
  template<class T>
  Frame<T> operator/(const Frame<T> &f, double d){
    Frame<T> r=f;
    for(unsigned int iL=0;iL<r.nL();iL++)
      for(unsigned int iF=0;iF<r.nF();iF++)
	r[iL][iF]=f[iL][iF]/d;
    return r;
  }
  
  template<class T, class U>
  Frame<T> operator-(const Frame<T> &p, const Frame<U> &m){
    Frame<T> r=p;
    for(unsigned int iL=0;iL<r.nL();iL++)
      for(unsigned int iF=0;iF<r.nF();iF++)
	r[iL][iF]=p[iL][iF]-m[iL][iF];
    return r;
  }
  
  template<class T>
  void Frame<T>::write(amjPacket &p) const{
    p.resize(2*sizeof(uint32_t)+_nL*_nF*sizeof(T));
    memcpy(p._data(),&_nL,sizeof(uint32_t));
    memcpy(p._data()+sizeof(uint32_t),&_nF,sizeof(uint32_t));
    for(uint32_t iL=0;iL<_nL;iL++)
      memcpy(p._data()+2*sizeof(uint32_t)+iL*_nF*sizeof(T),_frame[iL].data(),_nF*sizeof(T));
  }
  
  template<class T>
  void Frame<T>::write(uint8_t *p) const{
    memcpy(p,&_nL,sizeof(uint32_t));
    memcpy(p+sizeof(uint32_t),&_nF,sizeof(uint32_t));
    for(uint32_t iL=0;iL<_nL;iL++)
      memcpy(p+2*sizeof(uint32_t)+iL*_nF*sizeof(T),_frame[iL].data(),_nF*sizeof(T));
  }
  
  template<class T>
  void Frame<T>::read(amjPacket &p){
    memcpy(&_nL,p._data(),sizeof(uint32_t));
    memcpy(&_nF,p._data()+sizeof(uint32_t),sizeof(uint32_t));
    _frame.resize(_nL);
    for(uint32_t iL=0;iL<_nL;iL++){
      _frame[iL].resize(_nF);
      memcpy(_frame[iL].data(),p._data()+2*sizeof(uint32_t)+iL*_nF*sizeof(T),_nF*sizeof(T));
    }
  }
  
  template<class T>
  void Frame<T>::read1(uint8_t *p){
    memcpy(&_nL,p,sizeof(uint32_t));
    memcpy(&_nF,p+sizeof(uint32_t),sizeof(uint32_t));
  }
  
  template<class T>
  void Frame<T>::read2(uint8_t *p){
    _frame.resize(_nL);
    for(uint32_t iL=0;iL<_nL;iL++){
      _frame[iL].resize(_nF);
      memcpy(_frame[iL].data(),p+iL*_nF*sizeof(T),_nF*sizeof(T));
    }
  }
  
  float poidev(float xm, long *idum);
  
  template<class T, class U>
  void poisson(const Frame<T> &in, Frame<U> &out, long &seed){
    for(unsigned int iL=0;iL<in.nL();iL++)
      for(unsigned int iF=0;iF<in.nF();iF++)
	out[iL][iF]=poidev(in[iL][iF],&seed);
  }
    
  class Beam{
  public:
    Beam(double x=0, double D=0,
	 std::function<double(double)> delay=[](double t)->double{return 0;},
	 std::function<double(double)> illumination=[](double L)->double{return 1;}):
      _x(x),_D(D),_delay(delay),_illumination(illumination){};
    ~Beam(){};
    double x() const {return _x;}
    double D() const {return _D;}
    double delay(double t) const {return _delay(t);}
    double illumination(double L) const {return _illumination(L);}
  private:
    double _x; // Beam position in mm
    double _D; // Beam diameter in mm
    std::function<double(double)> _delay; // Beam delay in microns
    std::function<double(double)> _illumination; // Beam brightness as a function of wavelength in microns
  };
  
  class Baseline{
  public:
    Baseline(const std::string &name, Beam &beam1, Beam &beam2,std::complex<double>(*visibility)(double)
	     =[](double L)->std::complex<double>{return 1;})
      :_name(name),_beam1(beam1),_beam2(beam2),_visibility(visibility){};
    ~Baseline(){};
    std::string name(){return _name;};
    const Beam &beam1() const {return _beam1;}
    const Beam &beam2() const {return _beam2;}
    std::complex<double> visibility(double L) const {return _visibility(L);}
  private:
    std::string _name;
    Beam &_beam1,&_beam2;
    std::complex<double> (*_visibility)(double);
  };
  
  class Sim{
  public:
    Sim(const std::vector<Beam> &beams,const std::vector<Baseline> &baselines,
	       int nL=256, int nF=320, 
	       std::function<double(int)> wavelength=[](int i)->double{return 1/(1+(1/2.5-1)/256*i);},
	       std::function<double(int)> bandpass=[](int i)->double{return (1/(1+(1/2.5-1)/256*(i+1))-1/(1+(1/2.5-1)/256*(i-1)));},
	       double d=24, double f1=1203,double d2=52.5, double f2=44);
    ~Sim();
    void dim(int &_nL, int &_nF){_nL=nL;_nF=nF;};
    int frame(double t, Frame<double> &, double &, std::vector<double> &, std::vector<double> &, std::vector<double> &) const;
    // Fringe direction position as a function of pixel number. Pixel 0 is center of left-most pixel,
    // pixel -0.5 is left side of left-most pixel, pixel NF-1 is center of right-most pixel, and
    // pixel NF-0.5 is right side of right-most pixel. etc... returned position is relative to center
    // of row which is pixel number (-0.5+NF-0.5)/2=(NF-1)/2
    double x(double iF) const {return ((double)iF-((double)nF-1)/2)*d;};
  private:
    const std::vector<Beam> &beams; // reference to vector of beams
    const std::vector<Baseline> &baselines; // reference to vector of baselines
    std::vector<int> bi1; // index into beams of beam1 for each baseline
    std::vector<int> bi2; // index into beams of beam2 for each baseline
    int nL,nF; // Number of pixels in wavelength and fringe direction
    std::function<double(int)> wavelength; // Wavelength in microns
    std::function<double(int)> bandpass; // Bandpass in microns
    double d; // Pixel size in wavelength direction in microns
    double f1; // Focal length of M1 in mm
    double d2; // Distance from M1 focus to M2 in mm
    double f2; // Focal length of M2 in mm
    double m2; // Magnification in fringe direction by M2
    std::vector<std::vector<float> > airys; // Airy patterns from beams [nBeam][nL*nF]
    std::vector<float> illumination; // Sum of airy patterns [nL*nF];
    std::vector<std::vector<float> > airys2; // Airy patterns from baselines [nBaselines][2][nL*nF]
    void initialize();
    double bessj1(double x) const;  
    float envelope(float delay, float B, float L) const;
    double sinc(double x) const {if(fabs(x)<1e-6) return 1-x*x/6; return sin(x)/x;};
  };
  
  double wavelength(int i);

  class Phasors{
  public:
    Phasors(int nL, int nF, const std::vector<float> &periods={3.34},
	    std::function<double(int)> wavelength=amjFourier::wavelength);
    ~Phasors(){};
    //void resizephasors(PhasorSets &phasors);
    void operator()(const Frame<double> &frame, PhasorSets &);
  private:
    int nL; // number of wavelengths
    int nF; // number of pixels in fringe direction
    std::function<double(int)> wavelength;
    std::vector<float> periods; // [nB] - nB is number of baselines
    std::vector<float> c; // [iB][iL][iF]
    std::vector<float> s; // [iB][iL][iF]
  };
  
  
}

#endif /* _FOURIERSIM_H_ */
