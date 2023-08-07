#ifndef _FRAME_H_
#define _FRAME_H_

#include <vector>

template<class T>
class Frame{
public:
  Frame(int nL,int nF);
  ~Frame(){};
  int nL(){return _nL;}
  int nF(){return _nF;}
  template<class U>
  Frame<T> &operator=(const Frame<U> &);
  Frame<T> &operator*=(double);
  Frame<T> &operator+=(const Frame<T> &);
  std::vector<T> &operator[](int iL){return _frame[iL];};
private:
  int _nL,_nF;
  std::vector<std::vector<T> > _frame;
};

template<class T>
Frame<T>::Frame(int nL, int nF):_nL(nL),_nF(nF){
  _frame.resize(nL);
  for(int iL=0;iL<nL;iL++)
    _frame[iL].resize(nF);
}

template<class T>
template<class U>
Frame<T> &Frame<T>::operator=(const Frame<U> &f){
  _nL=f._nL;
  _nF=f._nF;
  _frame.resize(_nL);
  for(int iL=0;iL<nL;iL++)
    _frame[iL]=f._frame[iL];
}

template<class T>
Frame<T> &Frame<T>::operator*=(double f){
  for(int iL=0;iL<nL;iL++)
    for(int iF=0;iF<nF;iF++)
      _frame[iL][iF]*=f;
  return *this;
}

template<class T>
Frame<T> &Frame<T>::operator+=(const Frame<T> &f){
  for(int iL=0;iL<nL;iL++)
    for(int iF=0;iF<nF;iF++)
      _frame[iL][iF]+=f[iL][iF];
  return *this;
}

// template<class T, class U>
// Frame<T> &operator=(Frame<T> &f, Frame<U> &g){
//   f._nL=g._nL;
//   f._nF=g._nF;
//   f._frame.resize(_nL);
//   for(int iL=0;iL<f._nL;iL++){
//     _frame[iL].resize(f._nF);
//     for(int iF=0;iF<f._nF;iF++)
//       f._frame[iL][iF]=g._frame[iL][iF];
//   }
//   return *this;
// }

// template<class T, class U>
// Frame<T> &poisson(const Frame<U> &f){
//   f._nL=g._nL;
//   f._nF=g._nF;
//   f._frame.resize(nL);
//   for(int iL=0;iL,f._nL;iL++){
//     _frame[iL].resize(f._nF);
//     for(int iF=0;iF<f._nF;iF++)
//       f._frame[iL][iF]=poidev(g._frame[iL][iF]);
//   }
// }
  
#endif /* _FRAME_H_ */
