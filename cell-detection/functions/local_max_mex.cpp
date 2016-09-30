#include"mex.h"   
#include <stdio.h>
#include <stdlib.h>


void mexFunction( int Nreturned, mxArray *returned[], int Noperand, const mxArray *operand[] ){
        
    float *img;       
    (float *)img =  (float *)mxGetData(operand[0]);
    
    int *moveY;       
    (int *)moveY =  (int *)mxGetData(operand[1]);  // [y x];
    
    int *moveX;       
    (int *)moveX =  (int *)mxGetData(operand[2]);  // [y x];
    
    int *mask;       
    (int *)mask =  (int *)mxGetData(operand[3]);  // [y x];
    
    int ylim=mxGetM(operand[0]);
    int xlim=mxGetN(operand[0]);
    int imlen=ylim*xlim;
    int movelen=mxGetM(operand[1])*mxGetN(operand[1]);
    int masklen=mxGetM(operand[3])*mxGetN(operand[3]);
    
    unsigned char *output;       
    returned[0] = mxCreateNumericMatrix(ylim, xlim, mxUINT8_CLASS, mxREAL);
    output =  (unsigned char *)mxGetData(returned[0]);
    
    int before_inds, move_inds, minds;
    unsigned char tmp;
    float cent;
    
    
    for (int i=0; i<masklen; i++){
        tmp=1;
        minds=mask[i]-1;
        cent=img[minds];
        for (int j=0; j<movelen & tmp==1 ; j++){
            move_inds   = moveX[j]*ylim+moveY[j];
            before_inds = ((minds-move_inds)+imlen) % imlen;
            tmp *= (cent>=img[before_inds]);
        }
        output[minds]=tmp;
    }
            
    mexUnlock();  
}
