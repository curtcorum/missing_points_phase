#include "mex.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h> 

#define IDX2D(i1,i2,n1,n2) i1 + i2*n1
#define IDX3D(i1,i2,i3,n1,n2,n3) i1 + i2*n1 + i3*n1*n2

/* MATLAB syntax:  *****************************/
/* [outdata dcf] = cGrid(indata, kmap, kernel, opts) */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int ip, ix1, ix2, iy1, iy2, iz1, iz2, ix, iy, iz, x0, y0, z0;
	float cx, cy, cz, dz2, dy2, gx, lookupscale;
	mwSize nD;
	
	char errmsg[16384];
	strcpy(errmsg,"Invalid parameter configuration.\noutdata = cGrid(indata, kmap, kernel, opts);");
	strcat(errmsg,"\n  indata: 1D array, complex single precision, np");
	strcat(errmsg,"\n  kmap: 2D matrix, single precision, np x 3");
	strcat(errmsg,"\n  kernel: 1D array, single precision");
	strcat(errmsg,"\n  opts: 4 element struct");
	strcat(errmsg,"\n    opts.kernelradius");
	strcat(errmsg,"\n    opts.nx");
	strcat(errmsg,"\n    opts.ny");
	strcat(errmsg,"\n    opts.nz");

	/* Check for proper number of arguments */
	if (nrhs != 4 || nlhs > 2) mexErrMsgTxt(errmsg);

	/* Parse indata */
	nD = mxGetNumberOfDimensions(prhs[0]);
	const mwSize *indims1 = mxGetDimensions(prhs[0]);
	if(nD > 2)
	{
		mexErrMsgTxt(errmsg);
	}
	if(!mxIsSingle(prhs[0]))
	{
		mexErrMsgTxt(errmsg);
	}
	int np = indims1[0];
	float *rindata = (float*)mxGetData(prhs[0]);
	float *iindata = (float*)mxGetImagData(prhs[0]);
	/* mexPrintf("indata length = %d\n",np); */
	
	/* Parse kmap */
	nD = mxGetNumberOfDimensions(prhs[1]);
	const mwSize *indims2 = mxGetDimensions(prhs[1]);
	if(nD > 2)
	{
		mexErrMsgTxt(errmsg);
	}
	if(!mxIsSingle(prhs[1]))
	{
		mexErrMsgTxt(errmsg);
	}
	if(np != indims2[0])
	{
		mexErrMsgTxt(errmsg);
	}
	float *kmap = (float*)mxGetData(prhs[1]);

	/* Parse kernel */
	nD = mxGetNumberOfDimensions(prhs[2]);
	const mwSize *indims3 = mxGetDimensions(prhs[2]);
	if(nD > 2)
	{
		mexErrMsgTxt(errmsg);
	}
	if(!mxIsSingle(prhs[2]))
	{
		mexErrMsgTxt(errmsg);
	}
	int nk = indims3[1];
	float *kernel = (float*)mxGetData(prhs[2]);
	/* mexPrintf("kernel length = %d\n",nk); */

	/* Parse opts */
	if(!mxIsStruct(prhs[3]))
	{
		mexErrMsgTxt(errmsg);
	}
	if(mxGetNumberOfFields(prhs[3]) != 4)
	{
		mexErrMsgTxt(errmsg);
	}
	mxArray *field;
	field = mxGetField(prhs[3], 0, "kernelradius");
	int kr = (int)mxGetScalar(field);
	/* mexPrintf("kernel radius = %d\n",kr); */
	field = mxGetField(prhs[3], 0, "nx");
	int nx = (int)mxGetScalar(field);
	/* mexPrintf("nx = %d\n",nx); */
	field = mxGetField(prhs[3], 0, "ny");
	int ny = (int)mxGetScalar(field);
	/* mexPrintf("ny = %d\n",ny); */
	field = mxGetField(prhs[3], 0, "nz");
	int nz = (int)mxGetScalar(field);
	/* mexPrintf("nz = %d\n",nz); */
	x0 = floor(nx/2.0);
	y0 = floor(ny/2.0);
	z0 = floor(nz/2.0);
	
	/* Construct outdata */
	mwSize outdims[3];
	outdims[0] = nx;
	outdims[1] = ny;
	outdims[2] = nz;
	float *outdatar = (float*)mxCalloc((size_t)nx*ny*nz, sizeof(float));
	float *outdatai = (float*)mxCalloc((size_t)nx*ny*nz, sizeof(float));
	float *dcf = (float*)mxCalloc((size_t)nx*ny*nz, sizeof(float));
	plhs[0] = mxCreateNumericArray(3, outdims, mxSINGLE_CLASS, mxCOMPLEX);
	if(nlhs == 2)
		plhs[1] = mxCreateNumericArray(3, outdims, mxSINGLE_CLASS, mxREAL);

	/******* GRIDDING **************************/
	lookupscale = (nk-1)/kr;
	/* comment out pragmas for omp disabled debugging */
	#pragma omp parallel default(shared) private(ip, cx, cy, cz, ix1, ix2, iy1, iy2, iz1, iz2, iz, iy, ix, dz2, dy2, gx) num_threads(omp_get_num_procs())
	{
		#pragma omp for
		for(ip = 0; ip < np; ip++)
		{
			cx = kmap[IDX2D(ip,0,np,3)] + x0;
			/*if(cx > nx + kr) continue;*/
            if(cx > nx-1) continue;
			if(cx < 0) continue;
			cy = kmap[IDX2D(ip,1,np,3)] + y0;
			/*if(cy > ny + kr) continue;*/
            if(cy > ny-1) continue;
			if(cy < 0) continue;
			cz = kmap[IDX2D(ip,2,np,3)] + z0;
			/*if(cz > nz + kr) continue;*/
            if(cz > nz-1) continue;
			if(cz < 0) continue;
			ix1 = cx > kr ? ceil(cx - kr) : 0;
			ix2 = cx < nx - kr - 1? floor(cx + kr) : nx - 1;
			iy1 = cy > kr ? ceil(cy - kr) : 0;
			iy2 = cy < ny - kr - 1? floor(cy + kr) : ny - 1;
			iz1 = cz > kr ? ceil(cz - kr) : 0;
			iz2 = cz < nz - kr - 1? floor(cz + kr) : nz - 1;

			for(iz = iz1; iz <= iz2; iz++)
			{
				dz2 = (iz - cz)*(iz - cz);
				for(iy = iy1; iy <= iy2; iy++)
				{
					dy2 = (iy - cy)*(iy - cy);
					for(ix = ix1; ix <= ix2; ix++)
					{
						gx = sqrt(dz2 + dy2 + (ix - cx)*(ix - cx));
						if(gx < kr)
						{
							/* lookupscale = (nk-1)/kr; */
							outdatar[IDX3D(ix,iy,iz,nx,ny,nz)] += kernel[(int)(gx*lookupscale)]*rindata[ip]; 
							outdatai[IDX3D(ix,iy,iz,nx,ny,nz)] += kernel[(int)(gx*lookupscale)]*iindata[ip];
							dcf[IDX3D(ix,iy,iz,nx,ny,nz)] += kernel[(int)(gx*lookupscale)];
						}
					}
				}
			} 
		} 
	}
	mxSetData(plhs[0], (void *)outdatar);
	mxSetImagData(plhs[0], (void *)outdatai);
	if(nlhs == 2)
		mxSetData(plhs[1], (void *)dcf);
	else
		mxFree(dcf);
}
