// This is to test whether std::vectro<std::vector<std::vector<float>
// > > or a single std::vector with index computation is faster to
// extra values from and put values into.

#include <vector>
#include <iostream>

#include <time.h>

double tdiff(timespec t1, timespec t2);

int n0=1000,n1=3,n2=256,n3=320;
int i,j,k,l,m,n,o;
std::vector<float> a(n1*n2*n3),aa(n1*n2*n3);
std::vector<std::vector<std::vector<float> > > b;
float c;
timespec t1,t2,t3,t4,t5;

int main(int argc, char *argv[]){
  b.resize(n1);
  for(i=0;i<n1;i++){
    b[i].resize(n2);
    for(j=0;j<n2;j++)
      b[i][j].resize(n3);
  }

  c=1;

  clock_gettime(CLOCK_MONOTONIC,&t1);
  // Write with index
  for(i=0;i<n0;i++)
    for(j=0;j<n1;j++){
      m=j*n2*n3;
      for(k=0;k<n2;k++){
	n=m+k*n3;
	for(l=0;l<n3;l++){
	  o=n+l;
	  a[o]=c;
	  aa[o]=c;
	}
      }
    }

  clock_gettime(CLOCK_MONOTONIC,&t2);
  // Read with index
  for(i=0;i<n0;i++)
    for(j=0;j<n1;j++){
      m=j*n2*n3;
      for(k=0;k<n2;k++){
	n=m+k*n3;
	for(l=0;l<n3;l++){
	  o=n+l;
	  c=a[o]*aa[o];
	}
      }
    }

  
  clock_gettime(CLOCK_MONOTONIC,&t3);
  // Write with pointers
  c=1;
  for(i=0;i<n0;i++)
    for(j=0;j<n1;j++)
      for(k=0;k<n2;k++)
	for(l=0;l<n3;l++)
	  b[j][k][l]=c;

  clock_gettime(CLOCK_MONOTONIC,&t4);
  // Read with index
  for(i=0;i<n0;i++)
    for(j=0;j<n1;j++)
      for(k=0;k<n2;k++)
	for(l=0;l<n3;l++)
	  c=b[j][k][l]*b[j][k][l];
  
  clock_gettime(CLOCK_MONOTONIC,&t5);

  std::cout << tdiff(t1,t2) << std::endl;

  std::cout << tdiff(t2,t3) << std::endl;

  std::cout << tdiff(t3,t4) << std::endl;

  std::cout << tdiff(t4,t5) << std::endl;

  
  return 0;
}

double tdiff(timespec t1, timespec t2){
  return ((t2.tv_sec-t1.tv_sec)*1e9+(t2.tv_nsec-t1.tv_nsec))/1e9;
}
