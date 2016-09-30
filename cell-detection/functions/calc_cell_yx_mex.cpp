#include"mex.h"
#include <stdlib.h>

void mexFunction( int Nreturned,  mxArray *returned[], int Noperand, const  mxArray *operand[] ){
 int  *celllabel = (int *)mxGetData(operand[0]);
 int  *totcell = (int *)mxGetData(operand[1]);
 double  *yx = (double *)mxGetData(operand[2]);
 int  *len=  (int *)mxGetData(operand[3]);
 double  *im = (double *)mxGetData(operand[4]);
  
  double *y,*x,*area;         
  returned[0] = mxCreateNumericMatrix(*totcell,1,mxDOUBLE_CLASS,mxREAL);
  returned[1] = mxCreateNumericMatrix(*totcell,1,mxDOUBLE_CLASS,mxREAL);    
  returned[2] = mxCreateNumericMatrix(*totcell,1,mxDOUBLE_CLASS,mxREAL);  
  returned[3] = mxCreateNumericMatrix(*totcell,1,mxDOUBLE_CLASS,mxREAL); 
  (double *)y =    (double *)mxGetData(returned[0]); 
  (double *)x =    (double *)mxGetData(returned[1]);
  (double *)area = (double *)mxGetData(returned[2]);
  (double *)mean = (double *)mxGetData(returned[3]);
 
  int t;
  
  int i=*len; 
  while (i--){
     t=celllabel[i];
     area[t-1] +=1;
     y[t-1] +=yx[i];
     x[t-1] +=yx[i+(*len)];
  }
  
  i=*totcell; 
  while (i--){
     y[i] /=area[i];
     x[i] /=area[i];
  }
}

