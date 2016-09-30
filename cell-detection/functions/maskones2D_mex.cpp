
#include "mex.h"
#include <stdlib.h>

void make_template(int *Y, int *X,  int *inds,  int total,  int *moveinds, int total_move ){
    
  int i, j,l;
  for (j=0;j<total;j++){
    Y[j]=0;
    l=inds[j];
    i=total_move;
    while (i--){
            Y[j] +=X[l+moveinds[i]];
    }
  }
  
}


void mexFunction(int nargout, mxArray *returned[], int nargin, const mxArray *operand[]){ 
   
  int *dim = (int *)mxGetData(operand[0]);
  int  *mask = (int *) mxGetData(operand[1]);
  int  *mask_dim = (int *) mxGetData(operand[2]);
 
  
  int maskp[2],dimp[2],margin[2];
  
  maskp[0]=mask_dim[0]*2-1;
  maskp[1]=mask_dim[1]*2-1;
  
  dimp[0]=dim[0]+mask_dim[1]-1;
  dimp[1]=dim[1]+mask_dim[1]-1;
  
  margin[0]=(mask_dim[0]-1)/2;
  margin[1]=(mask_dim[1]-1)/2;
  
 
  int total;
  total=dim[0]*dim[1];  
  
          
  int *Y,*T;
  returned[0] = mxCreateNumericMatrix(total,1,mxINT32_CLASS,mxREAL);
  returned[1] = mxCreateNumericMatrix(mask_dim[1]*mask_dim[0],1,mxINT32_CLASS,mxREAL);
  (int *)Y = (int *)mxGetData(returned[0]);
  (int *)T = (int *)mxGetData(returned[1]);
  
  // make move index ///////////////////////////////////////////////////////////////////////////

 int i0,i1, d0,d1, dp0,dp1, dt0,dt1;
 int *moveinds1;
 int total_move=0;
 
 moveinds1=(int *)calloc(mask_dim[1]*mask_dim[0],sizeof(int));
 
 for (i1=0;i1<mask_dim[1];i1++){
     dt1=maskp[0]*(i1-margin[1]);
     d1=mask_dim[0]*i1;

     for (i0=0;i0<mask_dim[0];i0++){
         dt0=i0-margin[0];             
         if (mask[d1+i0]==1){
                 moveinds1[total_move]=dt1+dt0;
                 total_move +=1;
         }
     }
 } 

  // make template  ///////////////////////////////////////////////////////////////////////////
  
  
  int total_template=mask_dim[0]*mask_dim[1];
  int total_template_padding=maskp[0]*maskp[1];
  
 int *templ;
 templ=(int *)calloc(total_template_padding,sizeof(int));
 int *temp_inds;
 temp_inds=(int *)calloc(total_template,sizeof(int));
 
 for (i1=0;i1<mask_dim[1];i1++){    
     dp1=maskp[0]*(margin[1]+i1);
     d1 =mask_dim[0]*i1;    

     for (i0=0;i0<mask_dim[0];i0++){
         d0=d1+i0;
         dp0=dp1+(margin[0]+i0);
         templ[dp0]=1;
         temp_inds[d0]=dp0;
     }  
 }

 make_template(T,templ, temp_inds,  total_template, moveinds1, total_move );
  
 // make move index ///////////////////////////////////////////////////////////////////////////
  
 int t1,t0; 
     
 for (i1=0;i1<dim[1];i1++){
     if (i1<margin[1]){t1=i1;}else{
     if (i1>=dim[1]-margin[1]){t1=mask_dim[1]-(dim[1]-i1);}else{
     t1=margin[1];}} 

     dp1=mask_dim[0]*t1;
     d1 =dim[0]*i1;

     for (i0=0;i0<dim[0];i0++){             
         if (i0<margin[0]){t0=i0;}else{
         if (i0>=dim[0]-margin[0]){t0=mask_dim[0]-(dim[0]-i0);}else{
         t0=margin[0];}}                      

         d0=i0;
         dp0=t0;
         Y[d1+d0]=T[dp1+dp0];             

     }  
 }
 
 
}
