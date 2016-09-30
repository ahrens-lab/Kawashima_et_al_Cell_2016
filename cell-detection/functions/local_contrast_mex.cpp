#include"mex.h"   
#include <stdio.h>
#include <stdlib.h>


void mexFunction( int Nreturned, mxArray *returned[], int Noperand, const mxArray *operand[] ){
        
    float *img;       
    (float *)img =  (float *)mxGetData(operand[0]);
    
    int *grid;       
    (int *)grid =  (int *)mxGetData(operand[1]);  
    
    float *thre;       
    (float *)thre =  (float *)mxGetData(operand[2]); 
    
    int ylim=mxGetM(operand[0]);
    int xlim=mxGetN(operand[0]);
    int imlen=ylim*xlim;
    int gridn=grid[0];
    int gridlen=grid[0]*grid[0];
    
    int ygrids=ylim/ gridn;
    int xgrids=xlim/ gridn;
    int xmod=xlim % gridn;
    int xgrids2=xgrids+(xmod>0);
    int ngrids=xgrids2*ygrids;
    int gridlen2=xmod*gridn;
    float vthre=*thre;
    
    
    
    unsigned char *output;       
    returned[0] = mxCreateNumericMatrix(ylim, xlim, mxUINT8_CLASS, mxREAL);
    output =  (unsigned char *)mxGetData(returned[0]);
    
    
    // calculating index within each grids
    
    int *gridplus;
    gridplus=(int *)mxCalloc(gridlen, sizeof (int));
    
    
    int ii=0;
    for (int xx=0;xx<gridn;xx++)            
    {
        for (int yy=0;yy<gridn;yy++)
        {
            gridplus[ii]=xx*ylim+yy;      
            ii++;
        }
    }
    
    // calculating minimum value for each grid
    float vmin,tempv;
    unsigned char v=0;
    int sourceinds,inds;   
    
    float *minarray;
    minarray=(float *)mxCalloc(ngrids,sizeof(float));
    
    for (int xx=0;xx<xgrids;xx++)            
    {
        for (int yy=0;yy<ygrids;yy++)
        {
            sourceinds=xx*gridn*ylim + yy*gridn;
            vmin=100000;
            for (int zz=0; zz<gridlen; zz++)
            {
                tempv=img[sourceinds+gridplus[zz]];
                if (vmin > tempv){vmin=tempv;};
            }
            minarray[ygrids*xx+yy]=vmin;
        }
    }
    
    if (xmod>0){
                
        for (int yy=0;yy<ygrids;yy++)
        {
            sourceinds=gridn*xgrids*ylim + yy*gridn;
            vmin=100000;
            for (int zz=0; zz<gridlen2; zz++)
            {
                tempv=img[sourceinds+gridplus[zz]];
                if (vmin > tempv){vmin=tempv;};
            }
            minarray[ygrids*xgrids+yy]=vmin;
        }
    }
    
    // local averaging minimum value for each grid
    
    int xave[3]={-1,0,1};
    int yave[3]={-1,0,1};
    int ave_inds[9];
    float *minarray2;
    float tmp,tt;
    
    minarray2=(float *)mxCalloc(ngrids,sizeof(float));
    
    for (int iii=0;iii<3;iii++){
        for (int jjj=0;jjj<3;jjj++){
            ave_inds[iii*3+jjj]=xave[iii]*ygrids+yave[jjj];
        }
    }
    
    
    for (int xx=0;xx<xgrids2;xx++)            
    {
        for (int yy=0;yy<ygrids;yy++)
        {
            sourceinds=xx*ygrids+ yy;
            
            tt=0; tmp=0;
            for (int zz=0; zz<9; zz++)
            {
                inds=sourceinds+ave_inds[zz];
                if (inds>=0 && inds <ngrids)
                {
                    tt += 1;
                    tmp +=minarray[inds];
                }                    
            }
            minarray2[sourceinds]=tmp/tt;
        }
    }
        
    
    
    
    
    // thresholding each pixel on entire image
    
    
    
    for (int xx=0;xx<xgrids;xx++)            
    {
        for (int yy=0;yy<ygrids;yy++)
        {
            sourceinds=xx*gridn*ylim + yy*gridn;
         
            vmin=minarray2[ygrids*xx+yy];
            
            for (int zz=0; zz<gridlen ; zz++){
                inds=sourceinds+gridplus[zz];
                output[inds]=((img[inds]-vmin)>vthre);
            }
        }
    }
    
    if (xmod>0){
        
        for (int yy=0;yy<ygrids;yy++)
        {
            sourceinds=gridn*xgrids*ylim + yy*gridn;
            
            vmin=minarray2[ygrids*xgrids+yy];
            
            for (int zz=0; zz<gridlen2 ; zz++){
                inds=sourceinds+gridplus[zz];
                output[inds]=((img[inds]-vmin)>vthre);
            }
        }
    }
            
    mexUnlock();  
}
