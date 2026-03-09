#ifndef _AMJFOURIER_H_
#define _AMJFOURIER_H_

#include <vector>
#include <cstdint>
#include <iostream>

#include <complex>
#include <tuple>
#include <functional>
  
#include <amjTime.H>
#include <amjInterferometry.H>

// Default number of pixels in the camera
#define NL 256
#define NF 320

namespace amjFourier{
  template<typename T>
  class Frame{
  public:
    Frame(size_t wL=1,size_t wF=1);
    template<typename U>
    Frame(const Frame<U> &);
    ~Frame(){};
    void resize(size_t wL, size_t wF){_wL=wL; _wF=wF; resize();}
    size_t nL() const {return _nL;}
    void nL(size_t x) {_nL=x;}
    size_t nF() const {return _nF;}
    void nF(size_t x) {_nF=x;}
    size_t L0() const {return _L0;}
    void L0(size_t x) {_L0=x;}
    size_t F0() const {return _F0;}
    void F0(size_t x) {_F0=x;}
    size_t wL() const {return _wL;}
    void wL(size_t x){_wL=x; resize();}
    size_t wF() const {return _wF;}
    void wF(size_t x){_wF=x; resize();}
    amjTime &time(){return _t;}
    const amjTime &time() const {return _t;}
    void clear();
    template<typename U>
    Frame<T> &operator=(const Frame<U> &);
    Frame<T> &operator*=(double);
    Frame<T> &operator/=(double);
    template<typename U>
    Frame<T> &operator+=(const Frame<U> &);
    template<typename U>
    Frame<T> &operator-=(const Frame<U> &);
    std::vector<T> &operator[](size_t iL){return _frame[iL];}
    const std::vector<T> &operator[](size_t iL) const {return _frame[iL];}
    size_t size(){return size1()+size2();}
    size_t write(uint8_t *) const;
    size_t size1(){return 6*sizeof(uint32_t)+_t.size();}
    size_t read1(uint8_t *);
    size_t size2(){return _wL*_wF*sizeof(T);}
    size_t read2(uint8_t *);
  private:
    uint32_t _nL,_nF,_L0,_F0,_wL,_wF;
    amjTime _t;
    std::vector<std::vector<T> > _frame; /* [_wL][_wF] */
    void resize();
  };
  
  template<typename T>
  Frame<T>::Frame(size_t wL, size_t wF):_nL(NL),_nF(NF),_L0(0),_F0(0),
					_wL(wL),_wF(wF){
    resize();
  }
  
  template<typename T>
  template<typename U>
  Frame<T>::Frame(const Frame<U> &f){
    operator=(f);
  }
  
  template<typename T>
  void Frame<T>::clear(){
    for(uint32_t iL=0;iL<_nL;iL++)
      for(uint32_t iF=0;iF<_nF;iF++)
	_frame[iL][iF]=0;
  }
  
  template<typename T>
  template<typename U>
  Frame<T> &Frame<T>::operator=(const Frame<U> &f){
    _nL=f.nL();
    _nF=f.nF();
    _L0=f.L0();
    _F0=f.F0();
    _wL=f.wL();
    _wF=f.wF();
    _t=f.time();
    _frame.resize(_wL);
    for(uint32_t iL=0;iL<_wL;iL++){
      _frame[iL].resize(_wF);
      for(uint32_t iF=0;iF<_wF;iF++)
	_frame[iL][iF]=f[iL][iF];
    }
    return *this;
  }
  
  template<typename T>
  Frame<T> &Frame<T>::operator*=(double f){
    for(uint32_t iL=0;iL<_wL;iL++)
      for(uint32_t iF=0;iF<_wF;iF++)
	_frame[iL][iF]*=f;
    return *this;
  }
  
  template<typename T>
  Frame<T> &Frame<T>::operator/=(double f){
    for(uint32_t iL=0;iL<_wL;iL++)
      for(uint32_t iF=0;iF<_wF;iF++)
	_frame[iL][iF]/=f;
    return *this;
  }
  
  template<typename T>
  template<typename U>
  Frame<T> &Frame<T>::operator+=(const Frame<U> &f){
    if(_wL!=f.wL()||_wF!=f.wF()){
      std::cerr << "Frame<T>: frame dimensions do not match: (_wL,_wF)=("
		<< _wL << "," << _wF << ") and (f.wL,f.wF)=("
		<< f.wL() << "," << f.wF() << ")" << std::endl;
      abort();
    }
    for(uint32_t iL=0;iL<_wL;iL++)
      for(uint32_t iF=0;iF<_wF;iF++)
	_frame[iL][iF]+=f[iL][iF];
    return *this;
  }
  
  template<typename T>
  template<typename U>
  Frame<T> &Frame<T>::operator-=(const Frame<U> &f){
    if(_wL!=f._wL||_wF!=f._wF){
      std::cerr << "Frame<T>: frame dimensions do not match: (_wL,_wF)=("
		<< _wL << "," << _wF << ") and (f.wL,f.wF)=("
		<< f._wL << "," << f._wF << ")" << std::endl;
      abort();
    }
    for(uint32_t iL=0;iL<_wL;iL++)
      for(uint32_t iF=0;iF<_wF;iF++)
	_frame[iL][iF]-=f[iL][iF];
    return *this;
  }
  
  template<typename T>
  Frame<T> operator/(const Frame<T> &f, double d){
    Frame<T> r=f;
    for(uint32_t iL=0;iL<r.wL();iL++)
      for(uint32_t iF=0;iF<r.wF();iF++)
	r[iL][iF]=f[iL][iF]/d;
    return r;
  }
  
  template<typename T, typename U>
  Frame<T> operator-(const Frame<T> &p, const Frame<U> &m){
    Frame<T> r=p;
    for(uint32_t iL=0;iL<r.wL();iL++)
      for(uint32_t iF=0;iF<r.wF();iF++)
	r[iL][iF]=p[iL][iF]-m[iL][iF];
    return r;
  }

#include <vector>
#include <cstddef>
  
  template <typename T>
  void boxcar2D(Frame<T> &frame, size_t bL, size_t bF){
    const size_t wL = frame.wL();
    const size_t wF = frame.wF();
    
    // summed-area table
    Frame<T> sat(wL+1,wF+1);
    sat.clear();
    
    // build integral image
    for(size_t iL=1;iL<=wL;++iL){
      double row_sum=0.0;
      for(size_t iF=1;iF<=wF;++iF){
	row_sum+=frame[iL-1][iF-1];
	sat[iL][iF]=sat[iL-1][iF]+row_sum;
      }
    }
    
    const double norm=1.0/(bL*bF);
    double sum;
    
    // apply boxcar only where the window fully fits
    // edges are left untouched
    for(size_t iL=bL;iL<=wL;++iL){
      for(size_t iF=bF;iF<=wF;++iF){
	sum=sat[iL][iF]-sat[iL-bL][iF]-sat[iL][iF-bF]+sat[iL-bL][iF-bF];
	frame[iL-1][iF-1]=static_cast<T>(sum*norm);
      }
    }
  }
      
  template<typename T>
  size_t Frame<T>::write(uint8_t *pp) const{
    uint8_t *p=pp;
    memcpy(p,&_nL,sizeof(uint32_t)); p+=sizeof(uint32_t);
    memcpy(p,&_nF,sizeof(uint32_t)); p+=sizeof(uint32_t);
    memcpy(p,&_L0,sizeof(uint32_t)); p+=sizeof(uint32_t);
    memcpy(p,&_F0,sizeof(uint32_t)); p+=sizeof(uint32_t);
    memcpy(p,&_wL,sizeof(uint32_t)); p+=sizeof(uint32_t);
    memcpy(p,&_wF,sizeof(uint32_t)); p+=sizeof(uint32_t);
    p+=_t.write(p);
    for(uint32_t iL=0;iL<_wL;iL++){
      memcpy(p,_frame[iL].data(),_wF*sizeof(T));
      p+=_wF*sizeof(T);
    }
    return p-pp;
  }
  
  template<typename T>
  size_t Frame<T>::read1(uint8_t *pp){
    uint8_t *p=pp;
    memcpy(&_nL,p,sizeof(uint32_t)); p+=sizeof(uint32_t);
    memcpy(&_nF,p,sizeof(uint32_t)); p+=sizeof(uint32_t);
    memcpy(&_L0,p,sizeof(uint32_t)); p+=sizeof(uint32_t);
    memcpy(&_F0,p,sizeof(uint32_t)); p+=sizeof(uint32_t);
    memcpy(&_wL,p,sizeof(uint32_t)); p+=sizeof(uint32_t);
    memcpy(&_wF,p,sizeof(uint32_t)); p+=sizeof(uint32_t);
    p+=_t.read(p);
    resize();
    return p-pp;
  }
  
  template<typename T>
  size_t Frame<T>::read2(uint8_t *pp){
    uint8_t *p=pp;
    for(uint32_t iL=0;iL<_wL;iL++){
      memcpy(_frame[iL].data(),p,_wF*sizeof(T));
      p+=_wF*sizeof(T);
    }
    return p-pp;
  }
  
  template<typename T>
  void Frame<T>::resize(){
    _frame.resize(_wL);
    for(unsigned int iL=0;iL<_wL;iL++)
      _frame[iL].resize(_wF);
  }
  
  float poidev(float xm, long *idum);
  
  template<typename T, typename U>
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
    Baseline(const std::string &name, Beam &beam1, Beam &beam2,
	     std::function<std::complex<double>(double)> visibility
	     =[](double L)->std::complex<double>{return 1;})
      :_name(name),_beam1(beam1),_beam2(beam2),_visibility(std::move(visibility)){};
    ~Baseline(){};
    std::string name(){return _name;};
    const Beam &beam1() const {return _beam1;}
    const Beam &beam2() const {return _beam2;}
    std::complex<double> visibility(double L) const {return _visibility(L);}
  private:
    std::string _name;
    Beam &_beam1,&_beam2;
    std::function<std::complex<double>(double)> _visibility;
  };
  
  class Sim{
  public:
    Sim(const std::vector<Beam> &beams,const std::vector<Baseline> &baselines);
    ~Sim();
    void set_wavelength(std::function<double(double)>);
    void set_bandpass(std::function<double(double)>);
    void set_center(std::function<double(double)>);
    void set_window(size_t L0, size_t F0, size_t wL, size_t wF);
    void set_max_wavelength(double);
    template<typename U>
    //void set_bias(std::function<double(int,int)> &);
    void dim(int &_nL, int &_nF){_nL=nL;_nF=nF;};
    template<typename T>
    void operator()(Frame<T> &);
    //template<typename T>
    //int frame(double t, Frame<T> &, double &, std::vector<double> &, std::vector<double> &, std::vector<double> &) const;
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
    size_t nL,nF; // Detector pixels in wavelength and fringe direction
    size_t L0,F0; // Window start corner
    size_t wL,wF; // Window pixels in wavelength and fringe direction
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
  void Sim::operator()(Frame<T> &frame){
    frame.time().now();
    frame.nL(nL);
    frame.nF(nF);
    frame.L0(L0);
    frame.F0(F0);
    if(frame.wL()!=wL||frame.wF()!=wF)
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
  
  class CalcPhasors{
  public:
    CalcPhasors();
    ~CalcPhasors(){};
    // Wavelength in um from channel in camera frame coordinates
    void set_wavelength(std::function<double(double)>);
    // Center of fringe packet in camera frame coordinates
    void set_center(std::function<double(double)>);
    // Channels to use: first (lowest number), number of channels, in
    // camera frame coordinates
    void set_channels(unsigned int, unsigned int);
    // Fringe periods at shortest wavelength = highest numbered
    // channel, in camera frame
    void set_periods(std::vector<double>);
    // Names of baselines associated with fringe periods
    void set_names(std::vector<std::string>);
    // Number of consecutive channels to sum to produce fewer phasors
    void set_channel_sum(unsigned int);
    template<typename T, typename U>
    void operator()(const Frame<T> &frame,
		    std::vector<amjInterferometry::Phasors<U> > &);
  private:
    unsigned int nL,nF; // number of pixels in camera in wavelength and fringe direction
    unsigned int L0,F0; // location of (0,0) pixel of window in camera frame
    unsigned int wL,wF; // number of pixels in window in wavelength and fringe direction
    unsigned int C0,nC; // c0 is first channel to use for phasors, and
			// nC is number of chanels from there. c0 is
			// in camera frame
    //unsigned int c0,c1; // range of channels in window to use from c0 to c1, inclusive
    unsigned int cAvg; // number of channels to average
    std::function<double(double)> wavelength; // wavelength of camera row
    std::function<double(double)> center; // center of fringe packet
					  // in camera coordinates
    std::vector<double> periods; // [nB] - nB is number of baselines
    std::vector<std::string> names;
    std::vector<float> c; // [iB][iL][iF]
    std::vector<float> s; // [iB][iL][iF]
    void names_sizing();
    void initialize();
  };

  template<typename T, typename U>
  void CalcPhasors::operator()(const Frame<T> &frame,
			       std::vector<amjInterferometry::Phasors<U> > &phasors){
    if(L0!=frame.L0()||F0!=frame.F0()||wL!=frame.wL()||wF!=frame.wF()){
      L0=frame.L0(); F0=frame.F0(); wL=frame.wL(); wF=frame.wF();
      initialize();
    }
    phasors.resize(periods.size());
    unsigned int i1,iL,i2,iF,i3;
    float real,imag;
    unsigned int c0,c1,nP;
    c0=(L0>C0?L0:C0)-L0; // use the larger
    c1=(L0+wL<C0+nC?L0+wL:C0+nC)-L0; // use the smaller 
    nP=(c1-c0)/cAvg;
    for(unsigned int iB=0;iB<periods.size();iB++){
      phasors[iB].name(names[iB]);
      phasors[iB].time(frame.time());
      //i1=iB*nL*nF;
      i1=iB*wL*wF;
      phasors[iB].resize(nP);
      for(unsigned int iP=0;iP<nP;iP++){
	real=imag=0;
	float wsum=0;
	for(iL=c0+iP*cAvg;iL<c0+(iP+1)*cAvg;iL++){
	  //i2=i1+iL*nF;
	  i2=i1+iL*wF;
	  for(iF=0;iF<wF;iF++){
	    i3=i2+iF;
	    real+=frame[iL][iF]*c[i3];
	    imag+=frame[iL][iF]*s[i3];
	  }
	  wsum+=wavelength(iL);
	}
	phasors[iB][iP].real(real);
	phasors[iB][iP].imag(imag);
	phasors[iB](iP)=wsum/cAvg;
      }
    }
  }
  
  double wavelength(double i); // wavelength for channel i in um
  double wavelength_lin(double i); // wavelength for channel i, linear in wavenumber 
  double bandpass(double i); // bandpass for channel i in um
  double center(double i); // center of fringes in channel i, in pixels
  double window(double i); // window function (now Gaussian with sigma=20)
}

#endif /* _FOURIERSIM_H_ */
