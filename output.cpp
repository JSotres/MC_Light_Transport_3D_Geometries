/******************************************************************
*	This file contains all the functions needed for scoring 
*	and scaling physical quantities and for transfer them to an 
*	output file.
*******************************************************************/
#include "transport.h"



		 	


/***********************************************************
 *      Allocate an array with index from nl to nh inclusive.
 *
 *      Original matrix and vector from Numerical Recipes in C
 *      don't initialize the elements to zero. This will
 *      be accomplished by the following functions. 
 ****/
double* output::AllocVector(int nl, int nh)
{
	double *v;
  	int i;
  
  	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
  	if (!v) nrerror("allocation failure in vector()");
  
  	v -= nl;
  	for(i=nl;i<=nh;i++) 
		v[i] = 0.0;       /* init. */
	return v;
}





/***********************************************************
 *      Allocate a matrix with row index from nrl to nrh 
 *      inclusive, and column index from ncl to nch
 *      inclusive.
 ****/
double** output::AllocMatrix(int nrl,int nrh,int ncl,int nch)
{
	int i,j;
  	double **m;
  
  	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
  	if (!m) nrerror("allocation failure 1 in matrix()");
  		m -= nrl;
  
  	for(i=nrl;i<=nrh;i++) 
	{
    		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
    		if (!m[i]) nrerror("allocation failure 2 in matrix()");
    		m[i] -= ncl;
  	}
  
  	for(i=nrl;i<=nrh;i++)
    		for(j=ncl;j<=nch;j++) 
			m[i][j] = 0.0;
  	return m;
}





/***************************************************************************
*    Allocate a 3D matrix with row index from nrl to nrh
*    inclusive, column index from ncl to nch inclusive,
*    and depth index from ndl to ndh inclusive
**************************************************************************/

double*** output::Alloc3DMatrix( unsigned int nrl, unsigned int nrh, unsigned int ncl, 
				unsigned int nch, unsigned int ndl, unsigned int ndh)
{
	int i, j, k, nrow=nrh-nrl+1, ncol=nch-ncl+1, ndep=ndh-ndl+1;
	double*** t;
	
	/*Allocate pointers to rows. */
	
	t=(double ***) malloc((size_t)((nrow)*sizeof(double**)));
	if (!t)
		nrerror("allocation failure 1 in Alloc3DMatrix().\n");
	t -=nrl;
	
	/*Allocate pointers to rows and set pointers to them*/
	
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+1)*sizeof(double*)));
	if (!t[nrl])
	{
		cout<<"nrow: "<<nrow<<".\tncol: "<<ncol<<endl;
		nrerror("allocation failure in Alloc3DMatrix().\n");
	}
	t[nrl]-=ncl;
	
	/*Allocate rows and set pointers to them*/
	t[nrl][ncl]=(double*)malloc((size_t)((nrow*ncol*ndep+1)*sizeof(double)));
	if ( !t[nrl][ncl])
		nrerror("allocation failure 3 in Alloc3DMatrix().\n");
	t[nrl][ncl]-=ndl;
	
	for ( j=ncl+1; j<=nch; j++)
		t[nrl][j]=t[nrl][j-1]+ndep;
	for ( i=nrl+1; i<=nrh; i++)
	{
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch; j++)
			t[i][j]=t[i][j-1]+ndep;
	}
	
	
	/*Ensure that all elements are initialized to zero*/
	
	for (i=nrl; i<=nrh; i++)
	{
		for ( j=ncl; j<=nch; j++)
			for ( k=ndl; k<=ndh; k++)
				t[i][j][k]=0.0;
	}
	
	/*Return pointer to array of pointers to rows. */
	
	return t;
}



/***********************************************************
 *      Release the memory.
 ****/
void output::FreeVector(double *v,int nl,int nh)
{
	free((char*) (v+nl));
}



/***********************************************************
 *      Release the memory.
 ****/
void output::FreeMatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
	int i;
  
  	for(i=nrh;i>=nrl;i--) 
  		free((char*) (m[i]+ncl));
  	free((char*) (m+nrl));
}



/***********************************************************
 *      Release the memory.
 ****/

void output::Free3DMatrix(double***t,unsigned int nrl, unsigned int nrh, unsigned int ncl,
 			unsigned int nch, unsigned int ndl, unsigned int ndh)
{
	free((char*)(t[nrl][ncl]+ndl));
	free((char*)(t[nrl]+ncl));
	free((char*)(t+nrl));
} 





/***********************************************************
 *      Report error message to stderr, then exit the program
 *      with signal 1.
 ****/
void output::nrerror(char error_text[])
     
{
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}



void output::Set_number_grid_elements(input &rin)
{
	its_nx=rin.Get_nx();
	its_ny=rin.Get_ny();
	its_nz=rin.Get_nz();
	
	
	
	
}
	



void output::InitOutputData( input &rin )
{
	int nlayers;
	nlayers=rin.Get_nlayers()-2;
	
	Set_number_grid_elements(rin);
	
	int nx=its_nx;
	int ny=its_ny;
	int nz=its_nz;
	
	
	
	its_Rsra=0.0;
	its_Rdxy=AllocMatrix(0,nx-1,0,ny-1);
	its_R=0;
	its_Tdxy=AllocMatrix(0,nx-1,0,ny-1);
	its_T=0;
	its_Axyz=Alloc3DMatrix(0, nx-1, 0, ny-1, 0, nz-1 );
	its_Axy=AllocMatrix(0,nx-1,0,ny-1);
	its_Az=AllocVector(0,nz-1);
	its_Al=AllocVector(0,nlayers);
	its_A=0;
	its_Out=0.0;
	its_AOut=0.0;
	its_ROut=0.0;
	its_TOut=0.0;
	
	
}


void output::Add_to_Rdra_element(double extraR,int ir, int ia)
{
	its_Rdra[ir][ia] += extraR;
}

void output::Add_to_Rdxy_element(double extraR,int ix, int iy)
{
	its_Rdxy[ix][iy] += extraR;
}
void output::Add_to_Tdxy_element(double extraT,int ix, int iy)
{
	its_Tdxy[ix][iy] += extraT;
}

void output::Add_to_Tdra_element(double extraT,int ir, int ia)
{
	its_Tdra[ir][ia] += extraT;
}

void output::Add_to_Axyz_element(double extraA,unsigned int ix,unsigned int iy,unsigned int iz)
{
	its_Axyz[ix][iy][iz] += extraA;
}



void output::Add_to_Out( double i)
{
	if ( i<0.0 )
		nrerror("\nweight is negative!!!");
	its_Out+=i;
}

void output::Add_to_ROut( double i)
{
	if ( i<0.0 )
		nrerror("\nweight is negative!!!");
	its_ROut+=i;
}

void output::Add_to_TOut( double i)
{
	if ( i<0.0 )
		nrerror("\nweight is negative!!!");
	its_TOut+=i;
}

void output::Add_to_AOut( double i)
{
	if ( i<0.0 )
		nrerror("\nweight is negative!!!");
	its_AOut+=i;
}






void output::FreeMemory(input &rin)
{
	FreeMatrix(its_Rdxy,0,Get_nx()-1,0,Get_ny()-1);
	FreeMatrix(its_Tdxy,0,Get_nx()-1,0,Get_ny()-1);
	Free3DMatrix(its_Axyz,0, Get_nx(), 0,Get_ny(),0,Get_nz());
	FreeMatrix(its_Axy, 0, Get_nx(), 0,Get_ny() );
	FreeVector(its_Az, 0, Get_nz());
	FreeVector(its_Al, 0, rin.Get_nlayers()-2);	// is 2 because in the start you add two layers
							// the first one and the last one.
				
	
	
	
	
}







/************************************************************
*      Get 1D array elements by summing the 2D array elements.
*************************************************************/
void output::Sum2DRd()
{
	int nx = its_nx;
  	int ny = its_ny;
  	int ix,iy;
  	double sum;
  
  	
	sum = 0.0;
  	for ( ix=0; ix<nx; ix++)
		for ( iy=0; iy<ny; iy++)
			sum+=its_Rdxy[ix][iy];
	
  
  		
  	its_R = sum+its_ROut;
}



/***********************************************************
*      Return the index to the layer according to the index
*      to the grid line system in z direction (Iz).
*
*      Use the center of box.
************************************************************/
short output::IzToLayer(int Iz, input &rin)
{
	int i=0;    /* index to layer. */
  	int num_layers = rin.Get_nlayers();
  	double dz = rin.Get_gridz();
  
  	while( (Iz+0.5)*dz > rin.Get_zlpos(i+1) && i<(num_layers-1)) 
		i++;
  
  	return(i-1);
}






/***************************************************************
*      Get 1D array elements by summing the 2D array elements.
****************************************************************/
void output::Sum2DA(input &rin)
{
	int nz=Get_nz();
	int nx=Get_nx();
	int ny=Get_ny();
  	
	
  	int iz,ix, iy;
  	double sum;
  
  	for(ix=0; ix<nx; ix++)  
	{
		for(iy=0; iy<ny; iy++)
		{
    			sum = 0.0;
    			for(iz=0; iz<nz+1; iz++) 
				sum += its_Axyz[ix][iy][iz];
    			its_Axy[ix][iy] = sum;
		}
  	}
	
	for(iz=0; iz<nz; iz++)
	{
  		sum = 0.0;
  		for(ix=0; ix<nx; ix++) 
		{
			for(iy=0; iy<ny; iy++)
				sum += its_Axyz[ix][iy][iz];
		}	
		its_Az[iz] = sum;
	}
    		
		
	sum = 0.0;	
	for(iz=0; iz<nz; iz++)
	{
  		sum += its_Az[iz];
  		its_Al[IzToLayer(iz, rin)] += its_Az[iz];
  	}
  	its_A = sum+its_AOut;
}






/***********************************************************
*      Get 1D array elements by summing the 2D array elements.
**************************************************************/
void output::Sum2DTd()
{
	int nx = Get_nx();
  	int ny = Get_ny();
  	int ix,iy;
  	double sum;
  
  	
  	sum = 0.0;
	for ( ix=0; ix<nx; ix++)
		for ( iy=0; iy<ny; iy++)
			sum+=its_Tdxy[ix][iy];
	
	
  		
  	its_T = sum+its_TOut;
}







/***********************************************************
 *      Scale Rd and Td properly.
 *
 *      "a" stands for angle alpha.
 ****
 *      Scale Rd(r,a) and Tt(r,a) by
 *      (area perpendicular to photon direction)
 *              x(solid angle)x(No. of photons).
 *      or
 *              [2*PI*r*dr*cos(a)]x[2*PI*sin(a)*da]x[No. of photons]
 *      or
 *              [2*PI*PI*dr*da*r*sin(2a)]x[No. of photons]
 ****
 *      Scale Rd(r) and Tt(r) by
 *              (area on the surface)x(No. of photons).
 ****
 *      Scale Rd(a) and Tt(a) by
 *              (solid angle)x(No. of photons).
 ****/
void output::ScaleRdTd(input &rin)
{
	double scale1, scale2;
  	int nx = its_nx;
	int ny = its_ny;
  	double dx = rin.Get_gridx();
	double dy = rin.Get_gridy();
  	int ix, iy;
  	
	its_Out/=rin.Get_nphotons();
	
	scale1 = 1.0/(dx*dy*rin.Get_nphotons());
	
	for(ix=0; ix<nx; ix++)  
    		for(iy=0; iy<ny; iy++)
		{ 
			its_Rdxy[ix][iy] *= scale1;
			its_Tdxy[ix][iy] *= scale1;
		}
			
    		
  
  	
  	scale2 = 1.0/(double)rin.Get_nphotons();
  	its_R *= scale2;
  	its_T *= scale2;
	cout<<"its_Rsra:\t"<<its_Rsra<<endl;
	cout<<"its_R:\t\t"<<its_R<<endl;
	cout<<"its_T:\t\t"<<its_T<<endl;
	cout<<"its_Out:\t\t"<<its_Out<<endl;
	
}





/***********************************************************
 *      Scale absorption arrays properly.
 ****/
void output::ScaleA(input &rin)
{
	int nz = Get_nz();
  	int nx = Get_nx();
	int ny = Get_ny();
  	double dz = rin.Get_gridz();
	double dx = rin.Get_gridx();
	double dy = rin.Get_gridy();
  	int nl = rin.Get_nlayers()-2;
  	int iz,ix, iy;
  	int il;
  	double scale1;
  
  	/* Scale Axyz. */
  	scale1 = 1.0/dx*dy*dz*rin.Get_nphotons(); 
        	/* volume is dx*dy*dz.*/ 
        	
  	for(iz=0; iz<nz; iz++) 
    		for(ix=0; ix<nx; ix++)
			for(iy=0; iy<ny; iy++) 
      				its_Axyz[ix][iy][iz] *= scale1;
  
  	
	// Scale Axy
	scale1 = 1.0/(dx*dy*rin.Get_nphotons());
	for(ix=0; ix<nx; ix++)
			for(iy=0; iy<ny; iy++) 
      				its_Axy[ix][iy] *= scale1;
	
	
	/* Scale Az. */
  	scale1 = 1.0/(dz*rin.Get_nphotons());
  	for(iz=0; iz<nz; iz++) 
    		its_Az[iz] *= scale1;
  
  	
	
	
	/* Scale Al. Avoid int/int. */
  	scale1 = 1.0/(double)rin.Get_nphotons();     
  	for(il=0; il<nl; il++)
    		its_Al[il] *= scale1;
  
  	its_A *=scale1;
	
	cout<<"its_A:\t\t"<<its_A<<endl;
}

ofstream fout;

/***********************************************************
 *      Sum and scale results of current run.
 ****/
void output::SumScaleResult(input &rin)
{
	/* Get 1D & 0D results. */
  	Sum2DRd();
  	Sum2DA(rin);
  	Sum2DTd();
  
  	ScaleRdTd(rin);
  	ScaleA(rin);
	Total(rin);
}










	





/**********************************************
*	Sets the output file.
***********************************************/
void output::SetFile()
{
	cout<<"\nEnter name of output file( or . to exit ):\t";
	cin>>file;
	if(strlen(file) == 1 && file[0] == '.') 
      		exit(1);                  /* exit if no filename entered. */
 
}


/**********************************************
*	Gets the output file.
***********************************************/
char* output::GetFile()
{
	return (file);
}





/**********************************************
*	Opens the output file.
***********************************************/
void output::OpenFile(ofstream &rfout)
{
	rfout.open(file);
	if (! rfout.is_open())
		nrerror(" Could not open defined output file.\n");
}




/**********************************************
*	Closes the output file.
***********************************************/
void output::CloseFile(ofstream &rfout)
{
	rfout.close();
}


/**********************************************
*	Writes a title in the output file.
***********************************************/
void output::WriteHeadFile(ofstream &rfout)
{
	rfout<<"#This is the output file for the program.\n"<<endl;
}


/**********************************************
*	Writes the input data in the output 
*	file.
***********************************************/
void output::InputData(ofstream &rfout, input &rin)
{
	if ( rin.Get_specbeam() )
	{
		rfout<<"#x and y positions of the gaussian beam (cm)."<<endl;
		rfout<<rin.Get_xo()<<"\t"<<rin.Get_yo()<<endl;
		rfout<<"#standard desviation of the gaussian beam (cm)."<<endl;
		rfout<<rin.Get_stand()<<endl;
	}
	rfout<<"#Number of photon packets."<<endl;
	rfout<<rin.Get_nphotons()<<endl;
	rfout<<"#Maximum and minimum values of x and y where to represent the data"<<endl;
	rfout<<"#xmax	xmin	ymax	ymin"<<endl;
	rfout<<rin.Get_xmax()<<"\t"<<rin.Get_xmin()<<"\t"<<rin.Get_ymax()<<"\t"<<rin.Get_ymin()<<endl;
	rfout<<"# Maximun and minimun values of  x and y voxel index for each layer"<<endl;
	rfout<<"# ixmax     ixmin   iymax  iymin"<<endl;
	for ( int i=1; i<rin.Get_nlayers()-1; i++)
	{
		rfout<<rin.Get_bixmax(i)<<"\t"<<rin.Get_bixmin(i)<<"\t"<<rin.Get_biymax(i)<<"\t"<<rin.Get_biymin(i)<<endl;
	}
	rfout<<"# Number of grid elements for the output structure"<<endl;
	rfout<<"# ( for x, y, z )"<<endl;
	rfout<<rin.Get_nx()<<"\t"<<rin.Get_ny()<<"\t"<<rin.Get_nz()<<"\t"<<endl;
	rfout<<"# refraction index of the first layer"<<endl;
	rfout<<rin.Get_n(0)<<endl;
	rfout<<"# properties ( zlpos, dx, dy, dz, n, g, mua, mus) of every layer of the probe"<<endl;
	for ( int i=1; i<rin.Get_nlayers()-1; i++)
	{
		rfout<<rin.Get_zlpos(i)<<"\t"<<rin.Get_dx(i)<<"\t"<<rin.Get_dy(i)<<"\t"<<rin.Get_dz(i)<<"\t";
		rfout<<rin.Get_n(i)<<"\t"<<rin.Get_g(i)<<"\t"<<rin.Get_mua(i)<<"\t"<<rin.Get_mus(i)<<endl;
	}
	if ( rin.Get_specdata() )
	{
	rfout<<"#Cubes with different optical properties than those of the rest of the layer"<<endl;
	rfout<<"#layer, ix, iy, iz, n, g, mua, mus"<<endl;
	for ( int j=0; j<rin.Get_nspeclayers(); j++)
	{
	for ( int k=0; k<rin.Get_nspeccubes(j); k++)
	{
	rfout<<rin.Get_speclayer(j)<<"\t"<<rin.Get_specix(j,k)<<"\t"<<rin.Get_speciy(j,k)<<"\t"<<rin.Get_speciz(j,k);
	rfout<<"\t"<<rin.Get_specn(j,k)<<"\t"<<rin.Get_specg(j,k)<<"\t"<<rin.Get_specmua(j,k)<<"\t"<<rin.Get_specmus(j,k)<<endl;
	}
	}
	}
}
	







/**********************************************
*	Writes the total quantities scored
*	in the output 
*	file.
***********************************************/
void output::WriteGlobalResults(ofstream &rfout)
{
	rfout<<"\n# Reflectance, absorption, transmission."<<endl;
	rfout<<its_Rsra<<"\t\t#Specular reflectance [-]"<<endl;
	rfout<<its_R<<"\t\t#Diffuse reflectance [-]"<<endl;
	rfout<<its_A<<"\t\t#Absorbed fraction [-]"<<endl;
	rfout<<its_T<<"\t\t#Transmittance [-]\n"<<endl;
	rfout<<its_Out<<"\t\t#Lost in Boundaries [-]\n"<<endl;
	
}

/**********************************************
*	Writes the the absorption of every layer
*	in the output 
*	file.
***********************************************/
void output::WriteLayerAbsorption(ofstream &rfout, input &rin)
{
	rfout<<"#Absorption as a function of layer."<<endl;
	rfout<<"#layer\tA_l"<<endl;
	for ( int i=0; i<rin.Get_nlayers()-2; i++)
		rfout<<i+1<<"\t"<<its_Al[i]<<endl;
}

/**********************************************
*	Writes the the absorption of all the 
*	output grid elements 
*	in the output 
*	file.
***********************************************/
void output::WriteXYZAbsorption(ofstream &rfout, input &rin)
{
	rfout<<"Absorption in every grid element of the probe."<<endl;
	rfout<<"z1[cm]"<<endl;
	rfout<<"x1[cm]\t\ty1[cm]\t\tA(x1, y1, z1)"<<endl;
	rfout<<"x1[cm]\t\ty2[cm]\t\tA(x1, y2, z1)"<<endl;
	rfout<<"........."<<endl;
	rfout<<"z2[cm]"<<endl;
	rfout<<"x1[cm]\t\ty1[cm]\t\tA(x1, y1, z2)"<<endl;
	rfout<<"x1[cm]\t\ty2[cm]\t\tA(x1, y2, z2)"<<endl;
	rfout<<"........."<<endl;
	
	for ( int iz=0; iz<=Get_nz()-1; iz++)
	{
	rfout<<(iz+0.5)*rin.Get_gridz()<<endl;
	for ( int ix=0; ix<=Get_nx()-1; ix++)
	for ( int iy=0; iy<=Get_ny()-1; iy++)

	rfout<<rin.Get_xmin()+((ix+0.5)*rin.Get_gridx())<<"\t\t"<<rin.Get_ymin()+((iy+0.5)*rin.Get_gridy())
	<<"\t\t"<<its_Axyz[ix][iy][iz]<<endl;
	}

}	
	
/**********************************************
*	Writes the the absorption vs depth
*	in the output 
*	file.
***********************************************/
void output::WriteZetaAbsorption(ofstream &rfout, input &rin)
{
	rfout<<"#Absorption as a function of depth."<<endl;
	rfout<<"#z[cm]\tA_z[1/cm]"<<endl;
	for ( int iz=0; iz<=Get_nz()-1; iz++)
		rfout<<(iz+0.5)*rin.Get_gridz()<<"\t"<<its_Az[iz]<<endl;
}

/**********************************************
*	Writes the the diffuse reflectance
*	for all the grid elements
*	in the output file.
***********************************************/
void output::WriteXYReflectance(ofstream &rfout, input &rin)
{
	rfout<<"#Diffuse reflectance as a function of x and y coordinates."<<endl;
	rfout<<"#x[cm]\t\ty[cm]\t\tR(1/cm2)"<<endl;
	for ( int ix=0; ix<=Get_nx()-1; ix++)
		for ( int iy=0; iy<=Get_ny()-1; iy++)
			rfout<<rin.Get_xmin()+((ix+0.5)*rin.Get_gridx())<<"\t\t"<<rin.Get_ymin()+((iy+0.5)*rin.Get_gridy())
			<<"\t\t"<<its_Rdxy[ix][iy]<<endl;
}

/**********************************************
*	Writes the the diffuse reflectance
*	for all the grid elements
*	in the output file.
***********************************************/
void output::WriteXYTransmittance(ofstream &rfout, input &rin)
{
	rfout<<"#Diffuse transmittance as a function of x and y coordinates."<<endl;
	rfout<<"#x[cm]\t\ty[cm]\t\tT(1/cm2)"<<endl;
	for ( int ix=0; ix<=Get_nx()-1; ix++)
		for ( int iy=0; iy<=Get_ny()-1; iy++)
			rfout<<rin.Get_xmin()+((ix+0.5)*rin.Get_gridx())<<"\t\t"<<rin.Get_ymin()+((iy+0.5)*rin.Get_gridy())
			<<"\t\t"<<its_Tdxy[ix][iy]<<endl;
}


	

/**********************************************
*	Coordinates all the data transfer 
*	into the output file.
***********************************************/
void output::Total(input &rin)
{
	OpenFile(fout);
	WriteHeadFile(fout);
	InputData(fout, rin);
	WriteGlobalResults(fout);
	WriteLayerAbsorption(fout, rin);
	WriteZetaAbsorption(fout, rin);
	WriteXYZAbsorption(fout, rin);
	WriteXYReflectance(fout, rin);
	WriteXYTransmittance(fout, rin);
	CloseFile(fout);
}









