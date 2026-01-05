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

// Default number of pixels in the camera
#define NL 256
#define NF 320

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
    Sim(const std::vector<Beam> &beams,const std::vector<Baseline> &baselines);
    ~Sim();
    void set_wavelength(std::function<double(double)>);
    void set_bandpass(std::function<double(double)>);
    void set_center(std::function<double(double)>);
    void set_window(unsigned int L0, unsigned int F0, unsigned int wL,
		    unsigned int wF);
    void set_max_wavelength(double);
    template<typename U>
    //void set_bias(std::function<double(int,int)> &);
    void dim(int &_nL, int &_nF){_nL=nL;_nF=nF;};
    template<typename T>
    void operator()(double t, Frame<T> &);
    template<typename T>
    int frame(double t, Frame<T> &, double &, std::vector<double> &, std::vector<double> &, std::vector<double> &) const;
    // Fringe direction position as a function of pixel number. Pixel
    // 0 is center of left-most pixel, pixel -0.5 is left side of
    // left-most pixel, pixel NF-1 is center of right-most pixel, and
    // pixel NF-0.5 is right side of right-most pixel. etc... returned
    // position is relative to center of row which is pixel number
    // (-0.5+NF-0.5)/2=(NF-1)/2
    double x(double iL, double iF) const {return ((double)iF-center(iL))*d;};
  private:
    const std::vector<Beam> &beams; // reference to vector of beams
    const std::vector<Baseline> &baselines; // reference to vector of baselines
    std::vector<int> bi1; // index into beams of beam1 for each baseline
    std::vector<int> bi2; // index into beams of beam2 for each baseline
    unsigned int nL,nF; // Detector pixels in wavelength and fringe direction
    unsigned int L0,F0; // Window start corner
    unsigned int wL,wF; // Window pixels in wavelength and fringe direction
    double max_wavelength; // Longest wavelength visible
    std::function<double(double)> wavelength; // Wavelength in microns for ch
    std::function<double(double)> bandpass; // Bandpass in microns for ch
    std::function<double(double)> center; // Center of fringes for ch
    //std::function<double(int,int)> &bias; // makes bias frame
    double d; // Pixel size in wavelength direction in microns
    double f1; // Focal length of M1 in mm
    double d2; // Distance from M1 focus to M2 in mm
    double f2; // Focal length of M2 in mm
    double m2; // Magnification in fringe direction by M2
    std::vector<std::vector<float> > airys; // Airy patterns from beams [nBeam][wL*wF]
    std::vector<float> illumination; // Sum of airy patterns [wL*wF];
    std::vector<std::vector<float> > airys2; // Airy patterns from baselines [nBaselines][2][wL*wF]
    void initialize();
    double bessj1(double x) const;  
    float envelope(float delay, float B, float L) const;
    double sinc(double x) const {if(fabs(x)<1e-6) return 1-x*x/6; return sin(x)/x;};
  };
  
  template<typename T>
  void Sim::operator()(double t, Frame<T> &frame){
    if(frame.nL()!=wL||frame.nF()!=wF)
      frame.resize(wL,wF);
    
    unsigned int iB1,iB2;
    unsigned int iL,iF; // pixel location in the window
    unsigned int jL,jF; // pixel location in the frame
    float L,B,fdelay1,delay2,delay;
    std::complex<float> V0,V;
    int i,j;
    
    // Clear the frame
    frame.clear();
    
    // Frame illumination
    for(iL=0;iL<wL;iL++)
      for(iF=0;iF<wF;iF++)
	frame[iL][iF]=illumination[iL*wF+iF];

    // Fringes
    for(unsigned int iB=0;iB<baselines.size();iB++){
      iB1=bi1[iB];
      iB2=bi2[iB];
      fdelay1=(beams[iB1].x()-beams[iB2].x())/f1/m2;
      delay2=beams[iB2].delay(0)-beams[iB1].delay(0);
      for(iL=0,jL=L0;iL<wL;iL++,jL++){
	L=wavelength(jL);
	B=bandpass(jL);
	V0=baselines[iB].visibility(L);
	i=iL*wF;
	for(iF=0,jF=F0;iF<wF;iF++,jF++){
	  j=i+iF;
	  delay=fdelay1*x(jL,jF)+delay2;
	  V=V0*envelope(delay,L,B);
	  frame[iL][iF]+=2*sqrt(airys[iB1][j]*airys[iB2][j])*fabs(V)
	    *cos(2*M_PI*delay/L+std::arg(V));
	}
      }
    }     
  }

  template<typename T>
  class Bias{
  public:
    Bias();
    ~Bias();
    void seed(long);
    T operator()(int,int);
  private:
    
  };
  
  class Phasors{
  public:
    Phasors();
    // Phasors(int nL, int nF, const std::vector<float> &periods={3.34},
    // 	    std::function<double(int)> wavelength=amjFourier::wavelength);
    ~Phasors(){};
    // Wavelength in um from channel in camera frame coordinates
    void set_wavelength(std::function<double(double)>);
    // Center of fringe packet in camera frame coordinates
    void set_center(std::function<double(double)>);
    // Channels to use first, last, in window coordinates
    void set_channels(unsigned int, unsigned int);
    // Fringe periods at shortest wavelength = highest numbered
    // channel, in camera frame
    void set_periods(std::vector<double>);
    // Number of consecutive channels to average
    void set_channel_averaging(unsigned int);
    template<typename T>
    void operator()(int L0, int F0, const Frame<T> &frame, PhasorSets &);
  private:
    unsigned int nL,nF; // number of pixels in camera in wavelength and fringe direction
    unsigned int L0,F0; // location of (0,0) pixel of window in camera frame
    unsigned int wL,wF; // number of pixels in window in wavelength and fringe direction
    unsigned int c0,c1; // range of channels in window to use from c0 to c1, inclusive
    unsigned int cAvg; // number of channels to average
    std::function<double(double)> wavelength; // wavelength of camera row
    std::function<double(double)> center; // center of fringe packet
					  // in camera coordinates
    std::vector<double> periods; // [nB] - nB is number of baselines
    std::vector<float> c; // [iB][iL][iF]
    std::vector<float> s; // [iB][iL][iF]
    void initialize();
  };

  template<typename T>
  void Phasors::operator()(int L0_, int F0_, const Frame<T> &frame,
			   PhasorSets &phasors){
    if(L0!=L0_||F0!=F0_||wL!=frame.nL()||wF!=frame.nF()){
      L0=L0_; F0=F0_; wL=frame.nL(); wF=frame.nF();
      initialize();
    }
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
 
  double wavelength(double i); // wavelength for channel i in um
  double bandpass(double i); // bandpass for channel i in um
  double center(double i); // center of fringes in channel i, in pixels
}

#endif /* _FOURIERSIM_H_ */
