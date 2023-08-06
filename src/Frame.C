#include "../include/Frame.H"

template<class T>
Frame::Frame(int nL, int nF):_nL(nL),_nF(nF){
  _frame.resize(nL);
  for(int iL=0;iL<nL;iL++)
    _frame[iL].resize(nF);
}

template<class T>
Frame<T> &Frame::operator*=(double f){
  for(int iL=0;iL<nL;iL++)
    for(int iF=0;iF<nF;iF++)
      _frame[iL][iF]*=f;
  return *this;
}

template<class T>
Frame<T> &Frame::operator+=(const Frame<T> &f){
  for(int iL=0;iL<nL;iL++)
    for(int iF=0;iF<nF;iF++)
      _frame[iL][iF]+=f[iL][iF];
  return *this;
}

template<class T>
Frame<T> &Frame::operator=(const Frame<T> &f){
  _nL=f._nL;
  _nF=f._nF;
  _frame.resize(_nL);
  for(int iL=0;i<nL;iL++)
    _frame[iL]=f._frame[iL];
}

template<class T, class U>
Frame<T> &Frame::poisson(const Frame
