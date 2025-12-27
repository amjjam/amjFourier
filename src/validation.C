/* This is a test program which outputs values from several amjFourier
   functions for use in validation */

#include "../include/amjFourier.H"
#include <stdio.h>

int main(int argc, char *argv[]){
  // Wavelengths for channels 140 to 255
  FILE *fp=fopen("wavelength.dat","w");
  double w;
  for(int i=140;i<256;i++){
    fwrite(&i,sizeof(int),1,fp);
    w=amjFourier::wavelength(i);
    fwrite(&w,sizeof(double),1,fp);
  }
  
  fclose(fp);

  return 0;
}

