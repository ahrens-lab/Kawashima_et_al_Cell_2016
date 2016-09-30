#include"mex.h"


void localrank_assist_pc(int *Y1,const float *X1, const int *O, int *inds, int total, int *moveinds,int total_move){
    
  int i, j,index,mindex;  
  j=total;
  while (j--){
      index = *inds-1;
      i=total_move;      
      while (i--){
          mindex=index + moveinds[i];
          *Y1 +=(O[mindex])*(X1[index]>=X1[mindex]);
         
      }
      inds++;
      Y1++;
  }
  
}

void mexFunction( int Nreturned,  mxArray *returned[], int Noperand, const  mxArray *operand[] ){
    
  float *X1 = (float *)mxGetData(operand[0]);
  int *O = (int *)mxGetData(operand[1]);
  int *inds = (int *) mxGetData(operand[2]);
  int *moveinds = (int *) mxGetData(operand[3]);
  
  int total      = mxGetM(operand[2])*mxGetN(operand[2]);
  int total_move = mxGetM(operand[3])*mxGetN(operand[3]);
  
  int *Y1;
  returned[0] = mxCreateNumericMatrix(1,total, mxINT32_CLASS,0);  
  Y1 = (int *)mxGetData(returned[0]);
  
  localrank_assist_pc(Y1,X1,O,inds,total,moveinds,total_move);
}

