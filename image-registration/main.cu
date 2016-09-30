// A C/CUDA code for registering XY translations of timeseries.
// Registration will be separately performed for individual Z-planes.

// This program assumes microscope output where binary 3D stacks are written sequentially 
//(T00000.stack,T00001.stack,...).  This extracts time series of a single z-plane from these
// 3D stacks, register time series to their average, and write out registered 
// time series into a new binary file.  This process is repeated for all the z-planes, and 
// separate binary files for each z-plane time series are created (Plane00.stack, Plane01.stack,...).
// 2D registration of individual time series is based on phase correlation algorithm.
// Calculation is accelerated by GPU computation implemented by CUDA library.
//
// For general instruction of building .cu file please refer to online resources.
// Work with Windows Server 2012, Visual Studio 2010, Tesla K20 board and CUDA toolkit v5.5.
// Please compile as x64 software.  Libtiff library is required (http://www.libtiff.org/)
//
// required data format
// Image Files		-> Series of 3D stack binary files (UINT16). T00000.stack,T00001.stack,...
// Image Dimension File -> "StackDimensions.bin". 
//                         This uint32 binary contains 3 values (width, height, number of z-planes).
//
// input 
// "the source directory"		-> The data directory. Space is not allowed as part of the name.
// "the max timepoint"			-> Maximum time point. For the attached sample enter 59.
// "the number of digits in file name"	-> If the file name is like T00001.stack, enter 5.
//
// output ("registered" folder in the data directory)
// ave.tif		-> Averaged 3D stack after registration.
// PlaneXX.stack	-> Registered time series of individual Z-planes.
// stackdim.txt		-> Dimensions of the stack (width, height,number of z-planes, time length).

// Developed by Takashi Kawashima, HHMI Janelia Research Campus
// Sept 19, 2016

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cufft.h"
#include "cuda.h"
#include "cublas_v2.h"

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <windows.h>
#include <process.h>
#include <io.h>
#include <ctime>
#include <string>
#include <sstream>
#include <iomanip>
#include <locale>
#include <direct.h>
#include <omp.h>
#include <math.h>
#include "tiffio.h"
#include "tiff.h"


using namespace std;

typedef struct ST{ 
	
	unsigned short* outstack;
	int* zlist;
	int* status;
	int zrepeat;
	int imglen;
	unsigned short **outpoints;
	string fullpatho;

}thestruct; 

string convert_int_to_fname(int framenum, string fname_head, string fname_tail,int numdigit);
string convert_int_to_oname(int framenum, string fname_head, string fname_tail,int numdigit,int tt);
char*   str2char(string inputname);
void write_output(void* s);
void write_tiff_file(string path, float* stack, int w, int h, int zlen);
string ReplaceString(std::string subject, const std::string& search,const std::string& replace);




__global__ void MakeComp1( cufftDoubleComplex *a, cufftDoubleReal *b, int imglen)
{
	int id = (int) gridDim.x*blockDim.x*blockIdx.y+blockIdx.x *  blockDim.x +  threadIdx.x;
	if(id <imglen){
	
	    a[id].x  =  b[id];
	    a[id].y  =  0;  
	}

}

__global__ void MakeComp2( cufftDoubleComplex *a, unsigned short *b, unsigned short *c,int start,int imglen)
{
        int id = (int) gridDim.x*blockDim.x*blockIdx.y+blockIdx.x *  blockDim.x +  threadIdx.x;
        if(id <imglen){
            int id2= start+id;

            c[id]    =  b[id2];
            a[id].x  =  (double) b[id2];
            a[id].y  =  0;
        }
}


__global__ void ShiftPix( unsigned short *source, unsigned short *source2, int start, int shift,int imglen)
{
        int id = (int) gridDim.x*blockDim.x*blockIdx.y+blockIdx.x *  blockDim.x +  threadIdx.x ;
        if(id <imglen){    
            int id2 = (id+shift) % imglen;
            id2  += (id2 < 0 )*imglen;        
            source[start+id]  =  source2[id2];
        }
}
    

__global__ void KernelMult( cufftDoubleComplex *a, cufftDoubleComplex *b, cufftDoubleComplex *c, int imglen)
{
        int id = (int) gridDim.x*blockDim.x*blockIdx.y+blockIdx.x *  blockDim.x +  threadIdx.x;
        if(id <imglen){
            c[id].x  =  a[id].x * b[id].x+a[id].y * b[id].y;
            c[id].y  =  a[id].x * b[id].y-a[id].y * b[id].x;
        }
}

int main()
{
		

	
	string str1="\\";
	string str2="\\\\";
	
	string sdir, tdir;
	int maxtime, reftime,digits;

	cout << "Please enter the source directory \n>>";
	cin  >> sdir;

	tdir=sdir;

	cout << "Please enter the max timepoint\n>>";
	cin  >> maxtime;

	cout << "Please specify the number of digits in file name\n>>";
	cin  >> digits;
	
	reftime=maxtime/2;

	string indir = ReplaceString(sdir,str1, str2);
	indir.append("\\\\");

	string outdir= ReplaceString(tdir,str1, str2);
	outdir.append("\\\\\\Registered\\\\");

	string fname_head = "T";
	string fname_tail = ".stack";
	string oname_head = "Plane";
	string oname_tail = ".stack";
	string avetif    = "ave.tif";
	string stackdim    = "stackdim.txt";
	
	__int64 max_memory=3200000000; // this depends on the GPU board
	

	int numthreads=2;
	time_t tstart, tend; 
	tstart = time(0);
	int passedtime;
	
	cout<<"Starting...\n";
	string fullpathi;		
	string fullpatho;  


	

	string dimname = "StackDimensions.bin";

	string dimpath = indir;
	dimpath.append(dimname);
	unsigned int tmp[3];

	FILE *fdim = fopen(str2char(dimpath),"rb");
		fread(tmp, sizeof(unsigned int),3,fdim);
	fclose(fdim);
	
	


	// copy log files////////////////////////////////////////////////////

	
	_mkdir(str2char(outdir));




	
	printf("width=%d, height=%d, z=%d\n",tmp[0],tmp[1],tmp[2]);	

	int width  = (int)tmp[0];
	int height = (int)tmp[1];
	int zmin=1;
	int zmax   = (int)tmp[2];
	
			
	int imglen= width*height;
	int stacklen= width*height*(zmax-zmin+1);
	__int64 singlelen = imglen*sizeof(unsigned short);	
	


	//define input and output array/////////////////////////////////////////////////////////////////////////////////
	
				

	unsigned short *stack;
	stack =  (unsigned short *)calloc((__int64)imglen*(__int64)(maxtime+1),sizeof(unsigned short));
	unsigned short *stack_ori =  stack;	
	
	unsigned short **stackpoints;
    	stackpoints=(unsigned short**)malloc((maxtime+1)*sizeof(unsigned short*));

	
	for (int ii=0;ii<maxtime+1;ii++){
		stackpoints[ii]=stack+(__int64)imglen*(__int64)ii;
	}
				
	unsigned short *outstack;
	outstack=  (unsigned short *)calloc((__int64)imglen*(__int64)(maxtime+1),sizeof(unsigned short));
	unsigned short *outstack_ori =  outstack;	

	
	//define reference array/////////////////////////////////////////////////////////////////////////


	unsigned short *ref_i;
	(unsigned short *)ref_i=  (unsigned short *) calloc((__int64)stacklen ,sizeof(unsigned short));
	unsigned short *ref_i_ori=ref_i;
	double *ref_d;
	(double*)ref_d =  (double *) calloc((__int64)stacklen, sizeof(double));
	double* ref_d_ori=ref_d;


	
	//define calculation array/////////////////////////////////////////////////////////////////////////////////

	float *avestack;
	avestack=  (float *)calloc((__int64)stacklen, sizeof(float));
	float *avestack_ori=avestack;

	double *ave_d;
	ave_d=  (double *)calloc((__int64)imglen, sizeof(double));
	double *ave_d_ori=ave_d;
	


	/// Setup CUDA //////////////////////////////////////////////////////////////////////////////////////////////////////
	

	cudaSetDevice(1);
	cudaDeviceReset();
    
	int  zstep= (int) (max_memory/((__int64)(imglen*sizeof(unsigned short))));
	int  zrepeat= (maxtime+1)/zstep;

	int  amari=(maxtime+1)-zstep*zrepeat;
	if (amari!=0){
		zrepeat +=1;
	}
        
	int* zlist;
	zlist=(int*)calloc(zrepeat,sizeof(int));
    
	for (int i=0;i<zrepeat; i++){
		if(amari != 0 && i == zrepeat-1){
			zlist[i]=amari;
		}
		else{
			zlist[i]=zstep;
		}
	}

	
	unsigned short **outpoints;
	outpoints = (unsigned short**)malloc(zrepeat * sizeof (unsigned short*));

	outpoints[0]=outstack_ori;
	if(zrepeat>1)
	{
		for(int ii=1;ii<zrepeat;ii++){
			outpoints[ii]=outpoints[ii-1]+zlist[ii-1]*imglen;
		}
	}
	
    
	unsigned short  *source00,*source0;
	cufftDoubleReal *target0;
	cufftDoubleComplex *target1, *source1, *target2, *source2, *mult1, *mult2;
	cufftHandle fftPlan;
    
	cudaMalloc((void**)&target0,  sizeof(cufftDoubleReal)*imglen);
	cudaMalloc((void**)&source0,  sizeof(unsigned short)*((__int64)imglen*(__int64)zlist[0]));
	cudaMalloc((void**)&source00, sizeof(unsigned short)*imglen);

	cudaMalloc((void**)&target1, sizeof(cufftDoubleComplex)*imglen);
	cudaMalloc((void**)&source1, sizeof(cufftDoubleComplex)*imglen);       
	cudaMalloc((void**)&target2, sizeof(cufftDoubleComplex)*imglen);
	cudaMalloc((void**)&source2, sizeof(cufftDoubleComplex)*imglen);
	cudaMalloc((void**)&mult1,   sizeof(cufftDoubleComplex)*imglen);
	cudaMalloc((void**)&mult2,   sizeof(cufftDoubleComplex)*imglen);
    
    
    
	cufftPlan2d(&fftPlan, height, width, CUFFT_Z2Z);        
    
	cublasHandle_t handle;     
	cublasCreate(&handle);    
	cudaError_t ErrorHandle;
		
    
	
	int threads_num = 64;
    
	int g2= (imglen /(64*threads_num))+1;
	dim3 grids(64,g2,1);   
    
    
	int row_shift, col_shift;
	int shift,start;
	int peakind ;
    
	
	//create reference  images///////////////////////////////////////////////////////////////////////////////////////////
	
	printf("Creating reference image for alignment...\n") ;
	double refc=0;
	int rs=0;
	size_t readlen;

	for(int r=-15;r<15; r++)
	{
		ref_d=ref_d_ori;
		ref_i=ref_i_ori;
		int rtime=reftime+r;
		rs=0;

		if (rtime>0 && rtime<maxtime);
		{
			refc +=1;
			string fname= convert_int_to_fname(rtime, fname_head, fname_tail, digits);
			fullpathi=indir;
			fullpathi.append(fname);
		
	
			FILE *fi=fopen(str2char(fullpathi),"rb");   
			if(fi==NULL)
			{
				cout << "Invalid FilePath";
				return 0;
			}

			while(rs==0){
				_lseeki64(_fileno(fi), 0, SEEK_SET);  
				readlen=fread(ref_i,sizeof(unsigned short),(__int64)stacklen,fi);
				if (readlen==(__int64)stacklen){rs=1;}
			}

			fclose(fi);

			for (int i=0;i<stacklen;i++){
				*ref_d += (double)*ref_i;
				ref_d++;ref_i++;
			}	
		}
	}


	ref_d=ref_d_ori;

	for (int i=0;i<stacklen;i++){
		*ref_d /=refc;
		ref_d++;
	}	

	ref_d=ref_d_ori;
	
	tend = time(0); 
	passedtime=(int)difftime(tend,tstart);
	printf("%d sec elapsed\n",passedtime) ;

	//Process images///////////////////////////////////////////////////////////////////////////////////////////

	int zplane;
	for (zplane=zmin-1;zplane< zmax; zplane++){

		string oname= convert_int_to_fname(zplane+1, oname_head, oname_tail, 2);
		fullpatho=outdir;
		fullpatho.append(oname);
		FILE *fo=fopen(str2char(fullpatho),"rb");   
		if(fo!=NULL)
		{
			fclose(fo);
			remove(str2char(fullpatho));
		}
		tend = time(0); 		
		passedtime=(int)difftime(tend,tstart);

	

		//Acquire reference image//////////////////////////////////////////////////////////////////////////////


		ref_d=ref_d_ori+(__int64)imglen*(__int64)zplane;

	
		//read stack image///////////////////////////////////////////////////////////////////////////
	
		
	    
		omp_set_num_threads(numthreads);

		#pragma omp parallel default(none) shared() firstprivate(imglen,zplane,fname_head,fname_tail,indir,stackpoints)
		{
			string fullpathi2;	
			string fname2;
			FILE *fi2 = NULL;
			size_t result;
		    int i;

			#pragma omp for 
			for (int tt=0;tt<maxtime+1;tt++){
				fname2= convert_int_to_fname(tt, fname_head, fname_tail, digits);
				fullpathi2=indir;
				fullpathi2.append(fname2);
				i=0;
				fi2=fopen(str2char(fullpathi2),"rb");   

				while(i==0){
					_lseeki64(_fileno(fi2), singlelen*(__int64) (zplane), SEEK_SET);  
					result=fread(stackpoints[tt],sizeof(unsigned short),imglen,fi2);
					if (result==imglen){i=1;}
				}

				fclose(fi2);	
			
			}
		}
		
		omp_set_num_threads(1);

		
		
		
		tend = time(0); 
		passedtime=(int)difftime(tend,tstart);
		printf("%d sec elapsed\n",passedtime) ;
		
			
		//compute withGPU/////////////////////////////////////////////////////////////////////////////////

		stack    =stack_ori;   
		outstack =outstack_ori;   
    
    
		cudaMemcpy( target0,ref_d, sizeof(double)*imglen,cudaMemcpyHostToDevice);
		MakeComp1 <<< grids,  threads_num >>> (target1, target0, imglen);
		cufftExecZ2Z(fftPlan, target1, target2, CUFFT_FORWARD);

	
		HANDLE myhandleA;
		
		int* status;
		status=(int *)calloc(1,sizeof(int));
		status[0]=-1;
		
		thestruct st;

		st.outstack = outstack;
		st.zlist    = zlist;
		st.status   = status;
		st.zrepeat  = zrepeat;
		st.imglen    = imglen;
		st.fullpatho   = fullpatho;
		st.outpoints = outpoints;
		

		myhandleA = (HANDLE)_beginthread(write_output, 0, (void *)&st);

	


        
		for (int ii=0;ii<zrepeat;ii++){
        
			start  = 0;

			cudaMemcpyAsync(source0, stack, sizeof(unsigned short)*((__int64)imglen*(__int64)zlist[ii]),cudaMemcpyHostToDevice);
			status[0]=status[0]+1;
				
			
			for (int zz=0;zz<zlist[ii]; zz++){

				MakeComp2 <<< grids,  threads_num >>> (source1, source0,source00,start,imglen);
				cufftExecZ2Z(fftPlan, source1, source2, CUFFT_FORWARD);
				KernelMult <<< grids,  threads_num >>> (target2, source2, mult1,imglen);

				cufftExecZ2Z(fftPlan, mult1, mult2, CUFFT_INVERSE);
				cublasIzamax(handle,imglen,mult2,1,&peakind);

				row_shift = (peakind-1) % width;
				col_shift = (peakind-1) / width ;

				if (row_shift > (width/2)){
					row_shift -=  width;
				}

				if (col_shift > (height/2)){
					col_shift -=  height;
				}
				shift = width*col_shift+row_shift;
				ShiftPix <<< grids, threads_num >>> (source0, source00,start, shift, imglen);
				start += imglen;

			}
			

  			cudaMemcpyAsync(outpoints[ii], source0, sizeof(unsigned short)*imglen*zlist[ii],cudaMemcpyDeviceToHost);
			stack          += imglen*zlist[ii];

		}

		ErrorHandle = cudaGetLastError();
		status[0]=status[0]+1;

		
		
		
		tend = time(0); 
		passedtime=(int)difftime(tend,tstart);
		printf("%d sec elapsed\n",passedtime) ;
		
		WaitForSingleObject(myhandleA, INFINITE);


		// creating average image;///////////////////////////////////////////////////////

		outstack=outstack_ori;//outstack_ori;
		ave_d=ave_d_ori;

		for (int tt=0;tt<maxtime;tt++){

			for(int jj=0;jj<imglen;jj++){
				*ave_d += (double)*outstack;
				ave_d++;outstack++;
			}
			ave_d=ave_d_ori;
		}	
		outstack=outstack_ori;//outstack_o
		
		for(int jj=0;jj<imglen;jj++){
			*ave_d /= (double)maxtime;
			ave_d++;
		}
		ave_d=ave_d_ori;

		
		avestack=avestack_ori+imglen*zplane;

	
		for (int i=0;i<imglen;i++){
			*avestack = (float) *ave_d;
			avestack++; ave_d++;
		}	

			

		///disp time////////////////////////////////////////////////////////////////////
		

		tend = time(0); 
		passedtime=(int)difftime(tend,tstart);
		printf("Plane %d: Total %d sec elapsed\n",zplane+1,passedtime) ;
	}
		


	string avepath=outdir;
	avepath.append(avetif);

	avestack=avestack_ori;
	write_tiff_file(avepath, avestack, height, width, zmax);



	tend = time(0); 
	passedtime=(int)difftime(tend,tstart);
	printf("Total %d sec elapsed\n",passedtime) ;
	
	string dimmpath=outdir;
	dimmpath.append(stackdim);

	FILE *_fd=fopen(str2char(dimmpath),"w");	
	fprintf(_fd,"Y=%d \n",width);
	fprintf(_fd,"X=%d \n",height);
	fprintf(_fd,"Z=%d \n",maxtime);
	fclose(_fd);

	Sleep(5000);


	stack=stack_ori;
	outstack=outstack_ori;
	ref_i=ref_i_ori;
	ref_d=ref_d_ori;
	avestack=avestack_ori;
	ave_d=ave_d_ori;

	free(ref_i);
	free(ref_d);
	free(stack);
	free(outstack);
	free(avestack);
	free(ave_d);

    	cufftDestroy(fftPlan);
    	cudaFree(target0);
    	cudaFree(target1);
    	cudaFree(target2);
    	cudaFree(source0);
	cudaFree(source00);
    	cudaFree(source1);
    	cudaFree(source2);
    	cudaFree(mult1);
    	cudaFree(mult2);
    	cublasDestroy(handle);  

    	return 0;
}

void write_output(void* s)
{	
	thestruct* t = (thestruct *)s;
	FILE *_fo=fopen(str2char(t -> fullpatho),"ab+");

	size_t result;
	size_t position=0;
	__int64 slen;  
	int ok=0;

	for (int z=0;z<( t -> zrepeat); z++)
	{
		while( t -> status[0] <= z){
			Sleep(10);
		}
		slen=(__int64)(t -> imglen)*(__int64)(t -> zlist[z]);
		result=fwrite(t -> outpoints[z], sizeof(unsigned short),slen,_fo);

		if(result != slen)
		{

			while(ok==0)
			{
				_lseeki64(_fileno(_fo), position, SEEK_SET); 
				result=fwrite(t -> outpoints[z], sizeof(unsigned short),slen,_fo); 
				if(result ==slen){ok=1;};
			}
			
		}

		position=position+slen;
		ok=0;
	}

	fclose(_fo);
	_endthread();
}

string  convert_int_to_fname(int framenum, string fname_head, string fname_tail,int numdigit){


	
	ostringstream Convert;

	string Result, out;
	Convert << setw(numdigit) << setfill('0') << framenum;
	Result=Convert.str();

	out=fname_head;
	out.append(Result);
	out.append(fname_tail);

	return out;
}

string  convert_int_to_oname(int framenum, string fname_head, string fname_tail,int numdigit, int num){
	
	ostringstream Convert;

	string Result,out;
	Convert << setw(numdigit) << setfill('0') << framenum;
	Convert << setw(2) << setfill('0') << num;

	Result=Convert.str();

	out=fname_head;
	out.append(Result);
	out.append(fname_tail);

	return out;
}


char* str2char(string input){

	string search="\\";
	string replace="\\\\";
	size_t pos = 0;
	while ((pos = input.find(search, pos)) != std::string::npos) {
			input.replace(pos, search.length(), replace);
			pos += replace.length();
	}

	char *out=new char[input.length() +1];
	strcpy(out,input.c_str());

	return out;
}



string ReplaceString(std::string subject, const std::string& search,const std::string& replace) {
	size_t pos = 0;
	while ((pos = subject.find(search, pos)) != std::string::npos) {
			subject.replace(pos, search.length(), replace);
			pos += replace.length();
	}
	return subject;


}

void write_tiff_file(string path, float* stack, int w, int h, int zlen){
	
	TIFF *tif;

    unsigned short *pdst, *pdst_ori;
    pdst=(unsigned short*) calloc((size_t)w, sizeof (unsigned short));
	pdst_ori=pdst;

	
	tif=TIFFOpen(str2char(path),"wb");
	for (int pp=0;pp<zlen;pp++)
    {
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, h);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, 1);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
        TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, RESUNIT_NONE);
        TIFFSetField(tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
        TIFFSetField(tif, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
        TIFFSetField(tif, TIFFTAG_PAGENUMBER, (unsigned short) pp, (unsigned short) zlen);

        __int64 zstart    = (__int64)pp*w*h;
		int cstart;
        // now go line by line to write out the image data

        pdst=pdst_ori;
        for (int row = 0; row < h; row++ )
        {
			cstart=row+zstart;

            for (int col = 0; col < w; col++){
                *pdst++ = (unsigned short) stack[cstart+h*col];
            }
            pdst=pdst_ori;
            TIFFWriteScanline(tif, pdst, row, 0);
        }

        TIFFWriteDirectory(tif);
    }

    TIFFClose(tif);
}













