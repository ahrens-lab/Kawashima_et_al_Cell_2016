#include"mex.h"  
#include"matrix.h"   
#include <stdio.h>
#include <stdlib.h>
#include "tiffio.h"
#include "tiff.h"
#define UINT32  int
// This file needs libtiff library to build.
// http://www.libtiff.org/.


void mexFunction( int Nreturned, mxArray *returned[], int Noperand, const mxArray *operand[] ){

    // get the image information from the provided DIB
    
    
    char *input_buf;
    int buflen,status;
    unsigned char *imgdata;
    int* imgsize;
    
    buflen = (mxGetM(operand[0]) * mxGetN(operand[0])) + 1;
    input_buf = (char *)mxCalloc(buflen, sizeof(char));
    status = mxGetString(operand[0], input_buf, buflen);
        
    (unsigned char *)imgdata =  (unsigned char *)mxGetData(operand[1]);  // [y x];
    (int *)imgsize =  (int *)mxGetData(operand[2]);  // [y x];
    
    int sizelen = (mxGetM(operand[2]) * mxGetN(operand[2]));
    
    
    UINT32 w  = imgsize[1];
    UINT32 h  = imgsize[0];
    UINT32 ww = imgsize[1]*3;
    TIFF * tif;

    if ((w > 0) && (h > 0))
    {
       // open the output TIFF image for writing
        if ((tif = TIFFOpen(input_buf, "w")) == NULL){ return ;}
    }
    
        
    
    unsigned char* pdst;
    pdst=(unsigned char*) mxCalloc((size_t)ww, sizeof (unsigned char));
    unsigned char* pdst_ori=pdst;
    
    unsigned char *r_point,*g_point,*b_point;
    
    if (sizelen<3)
    {
        mexPrintf("Wrong format matrix: Y*X*3 or Y*X*Z*3 UINT8 matrix \n");
    }
    else if(sizelen==3)
    {
        if (imgsize[2]!=3)
        {
            mexPrintf("Wrong format matrix: Y*X*3 or Y*X*Z*3 UINT8 matrix \n");            
        }
        
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, h);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, 1);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 3);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE);
        TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
        TIFFSetField(tif, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);

        __int64 colorstep = (__int64) w*h;
        // now go line by line to write out the image data

        pdst=pdst_ori;
        for (int row = 0; row < h; row++ )
        {
            r_point = imgdata+(__int64)row;
            g_point = imgdata+(__int64)row+colorstep;
            b_point = imgdata+(__int64)row+colorstep*2;
            
            for (int col = 0; col < w; col++){
                *pdst++ = *r_point;r_point +=(__int64)h;
                *pdst++ = *g_point;g_point +=(__int64)h;
                *pdst++ = *b_point;b_point +=(__int64)h;
            }

            pdst=pdst_ori;
            TIFFWriteScanline(tif, pdst, row, 0);
        }
            
        
    }
    else if (sizelen==4)
    {        
        if (imgsize[3]!=3)
        {
            mexPrintf("Wrong format matrix: Y*X*3 or Y*X*Z*3 UINT8 matrix \n");            
        } 
        
        for (int pp=0;pp<imgsize[2];pp++)
        {
            TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w);
            TIFFSetField(tif, TIFFTAG_IMAGELENGTH, h);
            TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 8);
            TIFFSetField(tif, TIFFTAG_COMPRESSION, 1);
            TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
            TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 3);
            TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
            TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
            TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE);
            TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
            TIFFSetField(tif, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
            TIFFSetField(tif, TIFFTAG_PAGENUMBER, (unsigned short) pp, (unsigned short) imgsize[2]);

            __int64 zstart    = (__int64)pp*w*h;
            __int64 colorstep = (__int64)imgsize[2]*w*h;
            // now go line by line to write out the image data

            pdst=pdst_ori;
            for (int row = 0; row < h; row++ )
            {
                r_point = imgdata+zstart+(__int64)row;
                g_point = imgdata+zstart+(__int64)row+colorstep;
                b_point = imgdata+zstart+(__int64)row+colorstep*2;
                
                for (int col = 0; col < w; col++){
                    *pdst++ = *r_point;r_point +=(__int64)h;
                    *pdst++ = *g_point;g_point +=(__int64)h;
                    *pdst++ = *b_point;b_point +=(__int64)h;
                }
                pdst=pdst_ori;
                TIFFWriteScanline(tif, pdst, row, 0);
            }

            TIFFWriteDirectory(tif);
        }
    }

    TIFFClose(tif);
 
}