/*******************************************************************
* This is the head file of the program, it defines all the classes, 
* with its methods, that will be used in the program.
********************************************************************/
#include <fstream.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include <limits.h>

#define PI 3.1415926
#define WEIGHT 1E-4             /* Critical weight for roulette. */
#define CHANCE 0.1              /* Chance of roulette survival. */
#define STRLEN 256              /* String length. */

#define Boolean char

#define SIGN(x) ((x)>=0 ? 1:-1)

# define max(A,B) ((A)>(B)?(A):(B))     /* Definición macro de Max */
# define min(A,B) ((A)>(B)?(B):(A))     /* Definición macro de Min */
#define round(x) ((x)>=0?(int)((x)+0.5):(int)((x)-0.5))

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

#define COSZERO (1.0-1.0E-12)   /* cosine of about 1e-6 rad. */

#define COS90D  1.0E-6          /* cosine of about 1.57 - 1e-6 rad. */

#define PARTIALREFLECTION 0     
  /* 1=split photon, 0=statistical reflection. */






/******************************************************************
* This class has all the methods for reading the input information
* from the input file and for returning it to the rest of the 
* classes in the program.
*******************************************************************/
class input
{
	public:
		void nrerror(char error_text[]);
		void SetFile();
		char* GetFile();
		char* Getbuf();
		void OpenFile(ifstream &);
		void CloseFile(ifstream &);
		void FindData(ifstream &);
		Boolean CommentLine();
		void ReadParam();
		void Read_initpos();		// reads the mean values of the initial positions
		void Read_stand();		// reads the standard desviation of the gaussian beam.
		double Get_stand(){ return its_stand;}
		void Read_nphotons();
		long Get_nphotons();
		void Read_nlayers();
		int Get_nlayers();
		void Read_inc();
		void Read_maxminind();
		void Read_layers();
		double Get_xo(){ return its_xo;}
		double Get_yo(){ return its_yo;}
		void Set_ncubesx();
		void Set_ncubesy();
		void Set_ncubesz();
		long Get_ncubesx(long i);
		long Get_ncubesy(long i);
		long Get_ncubesz(long i);
		double Get_incx(){ return its_incx;}	// don't know if i will continue using these ones
		double Get_incy(){ return its_incy;}
		int Get_bixmax(int i){ return its_bixmax[i];}	// these functions return the boundary limits
		int Get_bixmin(int i){ return its_bixmin[i];}	// in the x and y directions for every
		int Get_biymax(int i){ return its_biymax[i];}	// layer.
		int Get_biymin(int i){ return its_biymin[i];}
		double Get_zlpos( int i){ return ( its_zlpos[i]);}
		double Get_dx( int i ){ return its_dx[i]; }
		double Get_dy( int i ){ return its_dy[i]; }
		double Get_dz( int i ){ return its_dz[i]; }
		double Get_n(int i){ return its_n[i];}
		double Get_g(int i){ return its_g[i];}
		double Get_mua(int i){ return its_mua[i];}
		double Get_mus(int i){ return its_mus[i];}
		void Set_layer(int i){ its_layer=i;}
		int Get_layer(){ return its_layer;}
		void Read_firstlayer();
		void CheckParam();
		void Set_xmax(){ its_xmax=its_incx;}
		void Set_ymax(){ its_ymax=its_incy;}
		void Set_zmax(){ its_zmax=Get_zlpos(Get_nlayers()-1);}
		double Get_xmax(){ return its_xmax;}
		double Get_ymax(){ return its_ymax;}
		double Get_zmax(){ return its_zmax;}
		double Get_xmin(){ return its_xmin;}
		double Get_ymin(){ return its_ymin;}
		double Get_zmin(){ return its_zmin;}
		void Set_min();	
		void Read_outgrid();
		double Get_gridx(){ return its_gridx;}
		double Get_gridy(){ return its_gridy;}
		double Get_gridz(){ return its_gridz;}
		double Get_gridr(){ return its_gridr;}
		double Get_grida(){ return its_grida;}
		long Get_nx(){return its_nx;}	
		long Get_ny(){return its_ny;}	
		long Get_nz(){return its_nz;}
		void Set_grid();			// sets the dimensions of the grid elements.
		void AssignDimensions();
		void DeleteInputPointers();
		Boolean Get_specdata(){ return specdata;}
			
		
		
		double Get_coscrit0(int i){return its_coscrit0[i];}
		double Get_coscrit1(int i){return its_coscrit1[i];}
		void CriticalAngle();
		
		// Functions for getting special cubes properties
		int Get_maxnspeccubes(){ return its_maxnspeccubes;}
		int Get_nspeclayers(){ return its_nspeclayers;}
		int Get_nspeccubes( int i ){ return its_nspeccubes[i];}
		int Get_speclayer( int i ){ return its_speclayer[i];}
		int Get_specix( int i, int j){ return its_specix[i][j];}
		int Get_speciy( int i, int j){ return its_speciy[i][j];}
		int Get_speciz( int i, int j){ return its_speciz[i][j];}
		double Get_specn(int i, int j){ return its_specn[i][j];}
		double Get_specg(int i, int j){ return its_specg[i][j];}
		double Get_specmua(int i, int j){ return its_specmua[i][j];}
		double Get_specmus(int i, int j){ return its_specmus[i][j];}
		void AskGaussianBeam();
		Boolean Get_specbeam(){ return specbeam;}
		
		//Functions for reading special cubes properties
		
		void GiveDimensions();
		void DeleteSpecVarDimensions();
		void ReadNSpecLayers();
		void ReadNSpecCubes();
		void ReadSpecialCubes();
		void AskSpecInput();
		void Set_maxncubes(int i){ its_maxncubes=i;}
		long Get_maxncubes(){ return its_maxncubes;}
		void SetMaxNumCubes();
		
			
		
		
	private:
		char file[STRLEN];
		char special[STRLEN];
		char gaussian[STRLEN];
		Boolean specdata;
		Boolean specbeam;
		char buf[STRLEN];
		char buffer[STRLEN];
		long its_nphotons;
		int its_nlayers;
		double its_incx;	
		double its_incy;
		int *its_bixmax;	// voxel index for the boundary positions of eacg layer
		int *its_bixmin;
		int *its_biymax;
		int *its_biymin;
		double *its_zlpos;
		double *its_dx;
		double *its_dy;
		double *its_dz;
		double *its_n;
		double *its_g;
		double *its_mua;
		double *its_mus;
		long *its_ncubesx;
		long *its_ncubesy;
		long *its_ncubesz;
		int its_layer;
		double its_xmax;		// maximum value of x
		double its_ymax;		// maximum value of y
		double its_zmax;		// maximum value of z
		double its_rmax;		// maximum value of the radius ( change the origin of x and y to the center)
		double its_amax;		// maximum angular value
		double its_xmin;		// minimum value of x
		double its_ymin;		// minimum value of y
		double its_zmin;		// minimum value of z
		double its_rmin;		// minimum value of the radius ( change the origin of x and y to the center)
		double its_amin;		// minimum angular value
		double its_gridx;			// dimension of the grid separations in x 
		double its_gridy;			// dimension of the grid separations in y 
		double its_gridz;			// dimension of the grid separations in z 
		double its_gridr;			// dimension of the grid separations in r 
		double its_grida;			// dimension of the grid separations in a
		int its_nx;			// number of x grid elements for the output.
		int its_ny;			// number of y grid elements for the output.
		int its_nz;			// number of z grid elements for the output.
		int its_nr;			// number of radial grid elements for the output.
		int its_na;			// number of angular grid elements for the output.
		
		/*critical angles for reflection*/
		
		double *its_coscrit0;
		double *its_coscrit1;
		
		// variables used in order to get specific optical
		// properties of some cubes.
		
		int its_maxnspeccubes;
		int its_maxncubes;
		int its_nspeclayers;
		int *its_nspeccubes;
		int *its_speclayer;
		int **its_specix;
		int **its_speciy;
		int **its_speciz;
		double **its_specn;
		double **its_specg;
		double **its_specmua;
		double **its_specmus;
		
		
		double its_xo;		// these two variables define where the gaussian beam is centerd.
		double its_yo;
		double its_stand;	// standard desviation of the gaussian beam.
		
		
		
		
};








/*************************************************************
* This class has all the methods for setting the properties of
* the voxels of the layer the photon is in and the methods for
* returning this properties to the rest of the classes.
**************************************************************/
class cube
{
	public:
		void nrerror(char error_text[]);
		void Set_properties( input &);
		double Get_dx(){ return its_dx;}
		double Get_dy(){ return its_dy;}
		double Get_dz(){ return its_dz;}
		double Get_n(){return its_n;}
		double Get_g(){ return its_g;}
		double Get_mua(){ return its_mua;}
		double Get_mus(){ return its_mus;}
		void Set_n(double i){ its_n=i;}
		void Set_g(double i){ its_g=i;}
		void Set_mua(double i){ its_mua=i;}
		void Set_mus(double i){ its_mus=i;}
		
		
		 
	private:
		double its_dx;
		double its_dy;
		double its_dz;
		double its_n;
		double its_g;
		double its_mua;
		double its_mus;
		
		 
};





/*************************************************************
* This class has all the methods for scoring and scaling the 
* physical quantities obtained and for writing them to the
* output file.
**************************************************************/ 
class output
{
	public:
		void Inc_R( ); 
		void Inc_T( ); 
		double** Get_Rdra(){ return its_Rdra;}
		double** Get_Tdra(){ return its_Tdra;}
		double*** Get_Axyz(){ return its_Axyz;}
		double Get_Rsra(){ return its_Rsra;}
		void Set_Rsra(double Rsra){its_Rsra=Rsra;}
		void Set_Rdra( double **Rdra){its_Rdra=Rdra;}  
		void Set_Tdra( double **Tdra){its_Tdra=Tdra;}  
		void Set_Axyz( double ***Axyz){its_Axyz=Axyz;}  
		double* AllocVector(int , int );
		double** AllocMatrix(int ,int ,int ,int );
 		void FreeVector(double *,int ,int );
		void FreeMatrix(double **,int ,int ,int ,int );
		double*** Alloc3DMatrix( unsigned int , unsigned int , unsigned int ,
 				unsigned int , unsigned int , unsigned int );
		void nrerror(char error_text[]);
		void Free3DMatrix(double***,unsigned int , unsigned int , unsigned int ,
 			unsigned int , unsigned int , unsigned int );

		void InitOutputData( input &);
		void Add_to_Rdra_element(double ,int, int);
		void Add_to_Tdra_element(double ,int, int);
		void Add_to_Tdxy_element(double ,int , int );
		void Add_to_Axyz_element(double ,unsigned int ,unsigned int ,unsigned int);
		void Add_to_Rdxy_element(double ,int , int );
		void Add_to_Out(double );
		void Add_to_AOut(double );
		void Add_to_ROut(double );
		void Add_to_TOut(double );
		void Set_number_grid_elements(input &);
		short Get_nr(){ return its_nr;}
		short Get_na(){ return its_na;}
		long Get_nx(){ return its_nx;}
		long Get_ny(){ return its_ny;}
		long Get_nz(){ return its_nz;}
		void FreeMemory(input &);
		void Sum2DRd();
		short IzToLayer(int , input &);
		void Sum2DA(input &);
		void Sum2DTd();
		void ScaleRdTd(input &);
		void ScaleA(input &);
		void SumScaleResult(input &);
		void SetFile();
		char* GetFile();
		void CloseFile(ofstream &);
		void OpenFile(ofstream &);
		void WriteHeadFile(ofstream &);
		void WriteGlobalResults(ofstream &);
		void WriteLayerAbsorption(ofstream &, input &);
		void WriteZetaAbsorption(ofstream &, input &);
		void WriteRadialReflectance(ofstream &, input &);
		void WriteAngularReflectance(ofstream &, input &);
		void WriteRadialTransmittance(ofstream &, input &);
		void WriteAngularTransmittance(ofstream &, input &);
		void WriteXYReflectance(ofstream &, input &);
		void WriteXYTransmittance(ofstream &, input &);
		void WriteXYZAbsorption(ofstream &, input &);
		void InputData(ofstream &, input &);
		void Total(input &);
		
		
	private:
		double its_Rsra;		// Specular Reflectance [-]	
		double** its_Rdra;
		double **its_Rdxy;		// 2D distribution of diffuse reflectance
		double* its_Rda;
		double* its_Rdr;
		double its_R;
		double** its_Tdra;		// 2D distribution of total transmittance.
		double* its_Tdr;
		double* its_Tda;
		double **its_Tdxy;		// 2D distribution of diffuse reflectance
		double its_T;
		double*** its_Axyz;		// 3D probability density in turbid media over x,y,z.
		double** its_Axy;		// 2D probability density in turbid media over x,y.
		double* its_Az;
		double* its_Al;
		double its_A;
		double its_Out;			// portion of weight lost out of the region of interest. 
		double its_AOut;		// portion of absorbed weight lost out of the region of interest.
		double its_ROut;		// portion of reflected weight lost out of the region of interest.
		double its_TOut;		// portion of transmitted weight lost out of the region of interest.
		double its_resx;		// Resolution in the x direction
		double its_resy;		// Resolution in the y direction
		double its_resz;		// Resolution in the z direction
		double its_resr;		// Radial resolution
		double its_resa;		// Angular resolution
		int its_nr;			// number of grids radius elements
		int its_na;			// number of grids angular elements
		int its_nx;			// number of grids elements in the x direction
		int its_ny;			// number of grids elements in the y direction
		int its_nz;			// number of grids elements in the z direction
		char file[STRLEN];
		char buf[STRLEN];
		
		
		
		 
		 
};



/********************************************************************
* This class has all the methods for propagation of the photon
* in the scattering field.
*********************************************************************/
class photon
{
	public:
		void nrerror(char error_text[]);
		void SetCubeArray(cube***&, input &);
		void SetCubeProperties(cube***&, input &);
		void DeleteCubeArray(cube***&, input &);
		void Set_weight( double i){ its_weight= i;}
		void Mod_weight();
		double Get_weight(){ return its_weight;}
		void Set_layer( int layer){ its_layer=layer;}
		int Get_layer(){ return its_layer;}
		void Set_max(input &);
		void Set_min(input &);
		void Set_ux( double i){ its_ux=i;}
		void Set_uy( double i){ its_uy=i;}
		void Set_uz( double i){ its_uz=i;}
		double Get_ux(){ return its_ux;}
		double Get_uy(){ return its_uy;}
		double Get_uz(){ return its_uz;}
		double Get_xnew(){return its_xnew;}
		double Get_ynew(){return its_ynew;}
		double Get_znew(){return its_znew;}
		void Move(input &);
		Boolean CrossOrNot(cube ***&, input &);
		void MoveBack();
		Boolean HitBoundary(input &);
		double RandomNum(void);
		float ran3(int*);
		double Rspecular(input &);
		void Set_dead(int i){ its_dead=i;}
		Boolean Get_dead(){ return its_dead;}
		void Set_s(double s){ its_s=s;}
		double Get_s(){ return its_s;}
		void LaunchPhoton(input &, cube ***&, output &);
		double SpinTheta(double );
		void Spin(input &, cube ***&);
		double Get_r(input &);
		long Get_ixmax(){ return its_ixmax;}
		long Get_iymax(){ return its_iymax;}
		long Get_izmax(){ return its_izmax;}
		int Get_ixmin(){ return its_ixmin;}
		int Get_iymin(){ return its_iymin;}
		int Get_izmin(){ return its_izmin;}
		void Set_ix();
		void Set_iy();
		void Set_iz();
		int Get_ix(){return its_ix;}
		int Get_iy(){return its_iy;}
		int Get_iz(){return its_iz;}
		double Get_x(){ return its_x;}
		double Get_y(){ return its_y;}
		double Get_z(){ return its_z;}
		void CrossLayerOrNot(input &, output &, cube ***&);
		void CrossUpOrNot(input &, output &, cube ***&);
		void CrossDnOrNot(input &, output &, cube ***&);
		double RFresnel(double, double, double, double* );
		void RecordR (double, input &, output &);
		void RecordT (double, input &, output &);
		void Set_ixnew(input &);
		void Set_iynew(input &);
		void Set_iznew(input &);
		int Get_ixnew(){ return its_ixnew;}
		int Get_iynew(){ return its_iynew;}
		int Get_iznew(){ return its_iznew;}
		void Set_CubeLabels();
		void Set_LateralLabels(input &);
		void Set_NewCubeLabels(input &);
	        Boolean HitDifferentCube(input &, cube ***&);
		void CrossDfCbOrNot(input &, output &, cube ***&);
		void Hop(input &);
		void CrossDfCbDnOrNot(input &, output &,cube ***&, int );
		void CrossDfCbUpOrNot(input &, output &, cube ***&, int );
		void HopMoveCheck(input &, output &, cube ***&);
		void Roulette(cube ***&, input &);
		void Drop(input &, output &, cube ***&);
		void HopMoveBackDropSpin(input &, output &, cube ***&);
		void DoOneRun(input &, output &, cube ***&);
		void Set_LaunchLabels(input &);
		void ListPhValues();
		double Get_sx(){ return its_sx;}
		double Get_sy(){ return its_sy;}
		double Get_sz(){ return its_sz;}
		Boolean OutOfBoundaries(input &, cube ***&, output &);
		double Get_stotal(){ return its_stotal;}
		void Set_stotal(double i){ its_stotal=i;}
		void StepSize(cube ***&, input &);
		void StepMoveHopeCheck(input &, output &, cube ***&);
		void StepMoveHopeCheckDropSpin(input &, output &, cube ***&);
		double normal(double, input &  );
		double uniform ( double , double );
	private:
		double its_x;
		double its_y;
		double its_z;
		double its_xnew;
		double its_ynew;
		double its_znew;
		double its_xold;
		double its_yold;
		double its_zold;
		double its_r;
		int its_ix;
		int its_iy;
		int its_iz;
		int its_ixnew;
		int its_iynew;
		int its_iznew;
		int its_ixold;
		int its_iyold;
		int its_izold;
		long its_ixmax;
		long its_iymax;
		long its_izmax;
		int its_ixmin;
		int its_iymin;
		int its_izmin;
		double its_weight;
		int its_nlayers;
		int its_layer;
		double its_ux;
		double its_uy;
		double its_uz;
		double its_s;
		double its_sx;
		double its_sy;
		double its_sz;
		double its_Rspecular;
		int its_dead;
		double its_stotal;
		
		
	
	
};






