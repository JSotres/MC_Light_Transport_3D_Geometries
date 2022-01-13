/*****************************************************************
*	This file include all the function needed to read
*	data from the input file and to compute all the 
*	general properties of the probe.
******************************************************************/

#include "transport.h"

ifstream fin;



/***********************************************************
 *      Report error message to stderr, then exit the program
 *      with signal 1.
 ****/
void input::nrerror(char error_text[])
     
{
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

/**********************************************************************
*	Asks for the name of the input file and stores it  in the
*	input variable file.
***********************************************************************/
void input::SetFile()
{
	cout<<"\nEnter name of input file( or . to exit ):\t";
	cin>>file;
	if(strlen(file) == 1 && file[0] == '.') 
      		exit(1);                  /* exit if no filename entered. */
 
}


/***********************************************************************
*	Gets the name of the input file.
************************************************************************/
char* input::GetFile()
{
	return (file);
}




/**************************************************************************
*	Opens the input file.
***************************************************************************/
void input::OpenFile(ifstream &rfin)
{
	rfin.open(file);
	if (! rfin.is_open())
		nrerror(" Could not open defined file.\n");
}


/***************************************************************************
*	Closes the input file.
****************************************************************************/
void input::CloseFile(ifstream &rfin)
{
	rfin.close();
}


/************************************************************************
*	Finds the next data line.
*************************************************************************/
void input::FindData(ifstream &rfin)
{
	do
	{
		rfin.getline(buf,255);
	}while(CommentLine());
	
	
	
}

/************************************************************************
*	Gets the value of the character array buf.
*************************************************************************/
char* input::Getbuf()
{
	return(buf);
}






/***********************************************************
*      Return 1 if this line is a comment line in which the 
*      first non-space character is "#".
*
*      Also return 1 if this line is space line.
***********************************************************/
Boolean input::CommentLine()
{
  size_t spn, cspn;
  
  spn = strspn(buf, " \t");     
  /* length spanned by space or tab chars. */

  cspn = strcspn(buf, "#\n");
  /* length before the 1st # or return. */

  if(spn == cspn)       /* comment line or space line. */
        return(1);
  else                          /* the line has data. */         
        return(0);              
}






/****************************************************************
*	Reads the initial average positions of the photons.
******************************************************************/
void input::Read_initpos()
{
	double xo, yo;
	FindData(fin);
	
  	sscanf(buf, "%lf%lf",&xo,&yo);
	
	its_xo=xo;
	its_yo=yo;
}
	 

/****************************************************************
*	Reads the standard desviation of the gaussian beam.
******************************************************************/
void input::Read_stand()
{
	double stand;
	FindData(fin);
	
  	sscanf(buf, "%lf",&stand);
	its_stand=stand;
}
	



/*****************************************************************
*	Read the number of initial photons.
*****************************************************************/
void input::Read_nphotons()
{
	long l=0;
  
  	
	FindData(fin);
	
  	sscanf(buf, "%ld",&l);
	
	its_nphotons=l;
}

/******************************************************************
*	Gets the number of initial photons.
*******************************************************************/
long input::Get_nphotons()
{
	return its_nphotons;
}





/***************************************************************
*	Reads the number of overall layers.
****************************************************************/
void input::Read_nlayers()
{
	int l=0;
  
  	
	FindData(fin);
	
  	sscanf(buf, "%ld",&l);
	
	its_nlayers=l;
}

/******************************************************************
*	Gets the number of overall layers.
*******************************************************************/
int input::Get_nlayers()
{
	return (its_nlayers+2);
}




/******************************************************************
*	Reads the refractive index of the initial layer.
******************************************************************/
void input::Read_firstlayer()
{
	double b;
	FindData(fin);
	sscanf(buf, "%lf",&b);
	its_n[0]=b;
}


/******************************************************************
*	Reads the different properties of the layers of the probe.
*******************************************************************/
void input::Read_layers()
{
	double zlpos;
	double dx;
	double dy;
	double dz;
	double n;
	double g;
	double mua;
	double mus;
	
	
	for ( int i=1; i<Get_nlayers()-1; i++)
	{
		
		
		FindData(fin);
	
		sscanf(buf, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",&zlpos,&dx,&dy,&dz,&n,&g,&mua,&mus);
		
		its_zlpos[i]=zlpos;
		its_dx[i]=dx;
		its_dy[i]=dy;
		its_dz[i]=dz;
		its_n[i]=n;
		its_g[i]=g;
		its_mua[i]=mua;
		its_mus[i]=mus;
		
		
	}
	
	FindData(fin);
	
	sscanf(buf, "%lf%lf",&zlpos,&n);
	
	its_zlpos[Get_nlayers()-1]=zlpos;
	
	its_n[Get_nlayers()-1]=n;	
		
	
}



/***************************************************
*	Set the minimum values that the x and y 
*	dimensions can take.
*
*	For the moment they are set to 0, but these 
*	can be changed in future versions:
*****************************************************/
void input::Set_min()
{
	its_xmin=0.0;
	its_ymin=0.0;
	its_zmin=0.0;
	its_rmin=0.0;
	its_amin=0.0;
}
	
/***************************************************
*	Reads the dimensions in the x and y
*	directions.
****************************************************/
void input::Read_inc()
{
	double xmax, xmin, ymax, ymin;
	
	FindData(fin);
	
	sscanf(buf,"%lf%lf%lf%lf",&xmax, &xmin, &ymax,&ymin);
	
	its_xmax=xmax;
	its_xmin=xmin;
	its_ymax=ymax;
	its_ymin=ymin;
	
	
	its_incx=xmax-xmin;
	
	its_incy=ymax-ymin;
	its_zmin=0.0;
	
	
}



/***************************************************************
*	Reads from the input file the labels of the maximum and 
*	minimum voxels that will have a pointer with the 
*	optical properties assigned to them.
****************************************************************/
void input::Read_maxminind()
{
	int ixmax, ixmin, iymax, iymin;
	
	for ( int i=1; i<=its_nlayers; i++)
	{
		FindData(fin);
	
		sscanf(buf,"%ld%ld%ld%ld",&ixmax,&ixmin, &iymax, &iymin);
		
		its_bixmax[i]=ixmax;
		its_bixmin[i]=ixmin;
		its_biymax[i]=iymax;
		its_biymin[i]=iymin;
	}
}
		
		


/***************************************************
*	Read the number of grid elements.
****************************************************/
void input::Read_outgrid()
{
	int x, y, z;
	
	FindData(fin);
	
	sscanf(buf,"%ld%ld%ld",&x,&y, &z);
	
	its_nx=x;
	
	its_ny=y;
	
	its_nz=z;
}


/*****************************************************************
*	Sets the number of cubes in the x direction
*	for every layer.
******************************************************************/
void input::Set_ncubesx()
{
	double a;
	for ( int i=1; i<(Get_nlayers()-1); i++ )
	{
		a=(its_incx/Get_dx(i));
		its_ncubesx[i]=int(a+0.1);
	}	
		
}



/*****************************************************************
*	Sets the number of cubes in the y direction
*	for every layer.
******************************************************************/
void input::Set_ncubesy()
{
	double a;
	for ( int i=1; i<(Get_nlayers()-1); i++ )
	{
		a=(its_incy/Get_dy(i));
		its_ncubesy[i]=int(a+0.1);
	}
}


/*****************************************************************
*	Sets the number of cubes in the z direction
*	for every layer.
******************************************************************/
void input::Set_ncubesz()
{
	double a;
	double b;
	double c;
	double d;
	for ( int i=1; i<(Get_nlayers()-1); i++ )
	{
		a=((Get_zlpos(i+1)-Get_zlpos(i))/Get_dz(i));
		its_ncubesz[i]=int(a+0.1);
		
	}
}


/*******************************************************************
*	This three functions get the number of cubes in the three
*	cartesian coordinates for every layer.
********************************************************************/
long input::Get_ncubesx(long i){ return its_ncubesx[i];}


long input::Get_ncubesy(long i){ return its_ncubesy[i];}


long input::Get_ncubesz(long i){ return its_ncubesz[i];}




/***********************************************************
*      Compute the critical angles for total internal
*      reflection according to the relative refractive index
*      of the layer.
*      All layers are processed.
*************************************************************/
void input::CriticalAngle()
{
	short i=0;
  	double n1, n2;
  
  	for(i=1; i<=Get_nlayers()-1; i++)  
	{
    		n1 = Get_n(i);
    		n2 =Get_n(i-1);
    		its_coscrit0[i] = n1>n2 ?sqrt(1.0 - n2*n2/(n1*n1)) : 0.0;
    
    		n2 =Get_n(i+1);
    		its_coscrit1[i] = n1>n2 ?sqrt(1.0 - n2*n2/(n1*n1)) : 0.0;
  	}
}



/******************************************************************
*	Sets the number of grid elements.
*******************************************************************/
void input::Set_grid()
{
	its_gridx=its_incx/its_nx;
	its_gridz=its_zmax/its_nz;
	its_gridy=its_incy/its_ny;
	
}
	
/*****************************************************************
*	Assign dimensions to the general input variables.
*****************************************************************/
void input::AssignDimensions()
{
	its_zlpos=new double[Get_nlayers()];
	its_dx=new double[Get_nlayers()];
	its_dy=new double[Get_nlayers()];
	its_dz=new double[Get_nlayers()];
	its_n=new double[Get_nlayers()];
	its_g=new double[Get_nlayers()];
	its_mua=new double[Get_nlayers()];
	its_mus=new double[Get_nlayers()];
	its_ncubesx=new long[Get_nlayers()];
	its_ncubesy=new long[Get_nlayers()];
	its_ncubesz=new long[Get_nlayers()];
	its_coscrit0=new double[Get_nlayers()];
	its_coscrit1=new double[Get_nlayers()];
	its_bixmax=new int[Get_nlayers()];
	its_bixmin=new int[Get_nlayers()];
	its_biymax=new int[Get_nlayers()];
	its_biymin=new int[Get_nlayers()];
	
}

/*******************************************************************
*	Delete the input pointers. Used in the end of the program.
********************************************************************/
void input::DeleteInputPointers()
{
	delete[]its_zlpos;
	delete[]its_dx;
	delete[]its_dy;
	delete[]its_dz;
	delete[]its_n;
	delete[]its_g;
	delete[]its_mua;
	delete[]its_mus;
	delete[]its_ncubesx;
	delete[]its_ncubesy;
	delete[]its_ncubesz;
	delete[]its_coscrit1;
	delete[]its_coscrit0;
	delete[]its_bixmax;
	delete[]its_bixmin;
	delete[]its_biymax;
	delete[]its_biymin;
	
}
	


/********************************************************
*	Asks if there are any special input data that
*	the program should read from the input file.
*********************************************************/
void input::AskSpecInput()
{
	cout<<"\nAre there any special input data( y or n ):\t";
	cin>>special;
	if(strlen(special) == 1 && special[0] == 'y') 
      		specdata= 1;
	else
		if(strlen(special) == 1 && special[0] == 'n')
		{
			//cout<<"\nhere"<<endl;
			specdata= 0;
		}
		else
			nrerror("\nNo valid input. Remember to enter lowercase y or n.");
	
 
}
	


/********************************************************
*	Asks if there the user wants a infinite narrow 
*	beam or a "gaussian beam".
*********************************************************/
void input::AskGaussianBeam()
{
	cout<<"\nAre there any specifications in order to use a gaussian beam?( y or n ):\t";
	cin>>gaussian;
	if(strlen(gaussian) == 1 && gaussian[0] == 'y') 
      		specbeam= 1;
	else
		if(strlen(gaussian) == 1 && gaussian[0] == 'n')
			specbeam=0;
		else
			nrerror("\nNo valid input. Remember to enter lowercase y or n.");
	
 
}




/****************************************************
*	Reads the number of layers with cubes with 
*	different optical properties than the
*	surrounding medium.
*****************************************************/
void input::ReadNSpecLayers()
{
	int l=0;
	
	
	FindData(fin);
	
  	sscanf(buf,"%ld",&l);
	
	
	its_nspeclayers=l;
	cout<<"\nits_nspeclayers:\t"<<its_nspeclayers<<endl;
	
	its_nspeccubes=new int[its_nspeclayers];
	its_speclayer=new int[its_nspeclayers];
	//cout<<"\naqui"<<endl;
	
	
}

/********************************************************
*	Reads the number of cubes with different optical
*	properties in each of the layers. The layers are 
*	not yet specified specified.
**********************************************************/
void input::ReadNSpecCubes()
{
	int nspeccubes=0;
	
	for ( int i=0; i<its_nspeclayers; i++)
	{
		FindData(fin);
		sscanf(buf,"%ld",&nspeccubes);
		its_nspeccubes[i]=nspeccubes;
		
	}
	
	its_maxnspeccubes=0;
	for ( int i=0; i<its_nspeclayers; i++)
		its_maxnspeccubes=max(its_maxnspeccubes,its_nspeccubes[i]);
	
}



/*********************************************************************
*	Get the maximum number of cubes.
**********************************************************************/
void input::SetMaxNumCubes()
{
	for ( int i=1; i<its_nlayers+1; i++ )
	{
		its_maxncubes=max(its_maxncubes,its_ncubesz[i]);
		its_maxncubes=max(its_ncubesx[i],its_maxncubes);
		its_maxncubes=max(its_maxncubes,its_ncubesy[i]);
		
	}
	
}

/*********************************************************************
*	Gives dimensions to all the arrays needed for giving special 
*	optical properties to certain cubes.
***********************************************************************/
void input::GiveDimensions()	
{
	int nlayers=its_nlayers;
	
	its_specix=new int*[its_nspeclayers];
	for ( int b=0; b<its_nspeclayers; b++)
		its_specix[b]=new int[its_maxnspeccubes];
	
	its_speciy=new int*[its_nspeclayers];
	for ( int b=0; b<its_nspeclayers; b++)
		its_speciy[b]=new int[its_maxnspeccubes];
	
	its_speciz=new int*[its_nspeclayers];
	for ( int b=0; b<its_nspeclayers; b++)
		its_speciz[b]=new int[its_maxnspeccubes];
	
	
	
	its_specn= new double*[its_nspeclayers];
	for ( int b=0; b<its_nspeclayers; b++)
		its_specn[b]=new double[its_maxnspeccubes];
	
	
	its_specg= new double*[its_nspeclayers];
	for ( int b=0; b<its_nspeclayers; b++)
		its_specg[b]=new double[its_maxnspeccubes];
	
	
	its_specmua= new double*[its_nspeclayers];
	for ( int b=0; b<its_nspeclayers; b++)
		its_specmua[b]=new double[its_maxnspeccubes];
	
	
	its_specmus= new double*[its_nspeclayers];
	for ( int b=0; b<its_nspeclayers; b++)
		its_specmus[b]=new double[its_maxnspeccubes];
	
	
	
}
	
	
	
/*********************************************************************
*	Deletes all the arrays needed for giving special 
*	optical properties to certain cubes.
***********************************************************************/
void input::DeleteSpecVarDimensions()
{
	int nlayers=Get_nlayers();
	delete[]its_nspeccubes;
	delete[]its_speclayer;	
	
	for ( long b=0; b<its_maxnspeccubes; b++)
		delete its_specix[b];
	delete[]its_specix;
	
	for ( long b=0; b<its_maxnspeccubes; b++)
		delete its_speciy[b];
	delete[]its_speciy;
	
	for ( long b=0; b<its_maxnspeccubes; b++)
		delete its_speciz[b];
	delete[]its_speciz;
	
	
	for ( long b=0; b<its_maxnspeccubes; b++)
		delete its_specn[b];
	delete[]its_specn;
	
	
	for ( long b=0; b<its_maxnspeccubes; b++)
		delete its_specg[b];
	delete[]its_specg;
	
	
	for ( long b=0; b<its_maxnspeccubes; b++)
		delete its_specmua[b];
	delete[]its_specmua;
	
	
	for ( long b=0; b<its_maxnspeccubes; b++)
		delete its_specmus[b];
	delete[]its_specmus;
	
	
	

}
	

/***********************************************************************
*	Reads the labels and properties of the cubes with different
*	optical properties than those that surround them.
************************************************************************/
void input::ReadSpecialCubes()
{
	for ( int i=0; i<its_nspeclayers; i++)
	{
		FindData(fin);
		int speclayer;
		sscanf(buf,"%ld",&speclayer);
		its_speclayer[i]=speclayer;
		
		for ( int j=0; j<its_nspeccubes[i]; j++ )
		{
			
			int ix, iy, iz;
			double n, mua, mus, g;
			FindData(fin);
			sscanf(buf,"%ld%ld%ld",&ix,&iy,&iz);
			its_specix[i][j]=ix;
			its_speciy[i][j]=iy;
			its_speciz[i][j]=iz;
			
			FindData(fin);
			sscanf(buf,"%lf%lf%lf%lf",&n,&g,&mua,&mus);
			its_specn[i][j]=n;
			its_specmua[i][j]=mua;
			its_specmus[i][j]=mus;
			its_specg[i][j]=g;
		
			
		}
	}
	
}
			

/*****************************************************************
*	Reads all the input parameters.
*****************************************************************/
void input::ReadParam()
{
	SetFile();
	OpenFile(fin);
	AskGaussianBeam();
	if ( specbeam )
	{
		Read_initpos();
		Read_stand();
		cout<<"\nxo:\t"<<its_xo<<endl;
		cout<<"\nyo:\t"<<its_yo<<endl;
		cout<<"\nstandard desviation:\t"<<its_stand<<endl;
	}
	else
	{
		Read_initpos();
		cout<<"\nxo:\t"<<its_xo<<endl;
		cout<<"\nyo:\t"<<its_yo<<endl;
	}
	Read_nphotons();
	Read_nlayers();
	AssignDimensions();
	Read_inc();
	Read_maxminind();
	Read_outgrid();
	Read_firstlayer();
	Read_layers();
	Set_ncubesx();
	Set_ncubesy();
	Set_ncubesz();
	Set_zmax();
	Set_grid();
	CriticalAngle();
	AskSpecInput(); 
	if ( specdata )
	{
		ReadNSpecLayers();
		ReadNSpecCubes();
		GiveDimensions();
		ReadSpecialCubes();
		CloseFile(fin);
		CheckParam();
	}
	else
	{
		its_nspeclayers=0;
		CloseFile(fin);
		CheckParam();
	}
	
	
	
	
}
/*****************************************************************
*	Check input parameters.
******************************************************************/
void input::CheckParam()
{
	if ( Get_nlayers()<=0 )
		nrerror( "No positive number of layers.\n");
	if ( Get_nphotons()<=0)
		nrerror(" No positive number of photons packets.\n");
	for ( int i=1; i<Get_nlayers()-1; i++ )
		if ( Get_n(i)<0 )
			nrerror(" index of refraction can't be negative.\n");
	for ( int i=1; i<=its_nlayers; i++ )
		if ( Get_g(i)<0 || Get_g(i)>1 )
		{
			cout<<"\ng("<<i<<"):\t"<<Get_g(i)<<endl;
			cout<<"\nnumber of layers:\t"<<its_nlayers<<endl;
			nrerror(" the anisotropy can't be smaller than zero or bigger than one.\n");
		}
	for ( int i=1; i<Get_nlayers(); i++ )
		if ( Get_mua(i)<0 || Get_mus(i)<0 )
			nrerror("Non valid scattering or absorption coefficients.\n");
			
	
	
	
	cout<<"\nits_nphotons:\t"<<its_nphotons<<endl;
	cout<<"\nits_nlayers:\t"<<its_nlayers<<endl;
	cout<<"\nits_incx:\t"<<its_incx<<endl;
	cout<<"\nits_incy:\t"<<its_incy<<endl;
	for ( int i=1; i<Get_nlayers()-1; i++ )
	{
	cout<<"\nits_zlpos["<<i<<"]:\t"<<its_zlpos[i]<<endl;
	cout<<"\nits_dx["<<i<<"]:\t"<<its_dx[i]<<endl;
	cout<<"\nits_dy["<<i<<"]:\t"<<its_dy[i]<<endl;
	cout<<"\nits_dz["<<i<<"]:\t"<<its_dz[i]<<endl;
	cout<<"\nits_n["<<i<<"]:\t"<<its_n[i]<<endl;
	cout<<"\nits_g["<<i<<"]:\t"<<its_g[i]<<endl;
	cout<<"\nits_mua["<<i<<"]:\t"<<its_mua[i]<<endl;
	cout<<"\nits_mus["<<i<<"]:\t"<<its_mus[i]<<endl;
	cout<<"\nits_ncubesx["<<i<<"]:\t"<<its_ncubesx[i]<<endl;
	cout<<"\nits_ncubesy["<<i<<"]:\t"<<its_ncubesy[i]<<endl;
	cout<<"\nits_ncubesz["<<i<<"]:\t"<<its_ncubesz[i]<<endl;
	cout<<"\nits_bixmax["<<i<<"]:\t"<<its_bixmax[i]<<endl;
	cout<<"\nits_bixmin["<<i<<"]:\t"<<its_bixmin[i]<<endl;
	cout<<"\nits_biymax["<<i<<"]:\t"<<its_biymax[i]<<endl;
	cout<<"\nits_biymin["<<i<<"]:\t"<<its_biymin[i]<<endl;
	
	}
	cout<<"\nits_xmax:\t"<<its_xmax<<endl;
	cout<<"\nits_ymax:\t"<<its_ymax<<endl;
	cout<<"\nits_zmax:\t"<<its_zmax<<endl;
	cout<<"\nits_xmin:\t"<<its_xmin<<endl;
	cout<<"\nits_ymin:\t"<<its_ymin<<endl;
	cout<<"\nits_zmin:\t"<<its_zmin<<endl;
	cout<<"\nits_gridx:\t"<<its_gridx<<endl;
	cout<<"\nits_gridy:\t"<<its_gridy<<endl;
	cout<<"\nits_gridz:\t"<<its_gridz<<endl;
	cout<<"\nits_nx:\t"<<its_nx<<endl;
	cout<<"\nits_ny:\t"<<its_ny<<endl;
	cout<<"\nits_nz:\t"<<its_nz<<endl;
	cout<<"\n\n***Please, check the input data!******\n\n"<<endl;
	
	
}
	










