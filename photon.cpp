/*********************************************************
*	This file includes the functions of the
*	photon class.
**********************************************************/
#include "transport.h"




/***********************************************************
*      Report error message to stderr, then exit the program
*      with signal 1.
*************************************************************/
void photon::nrerror(char error_text[])
     
{
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}





/******************************************************
*	set the array of cubes every time the photon 
*	enters a different layer
*******************************************************/

void photon::SetCubeArray (cube ***&pca, input &rin)

{
	
	pca=new cube**[rin.Get_bixmax(its_layer)-rin.Get_bixmin(its_layer)+1];
	for ( int a=0; a<=rin.Get_bixmax(its_layer)-rin.Get_bixmin(its_layer); a++)
	{
		pca[a]=new cube*[rin.Get_biymax(its_layer)-rin.Get_biymin(its_layer)+1];
	}
		
	for ( int a=0; a<=rin.Get_bixmax(its_layer)-rin.Get_bixmin(its_layer); a++)
		for ( int b=0; b<=rin.Get_biymax(its_layer)-rin.Get_biymin(its_layer); b++)
			{
			  	pca[a][b]=new cube[rin.Get_ncubesz(its_layer)+1];
			 
			}
	
	
	
			
}



/******************************************************
*	delete the array of cubes every time the photon 
*	enters a different layer
*******************************************************/


void photon::DeleteCubeArray( cube***&pca, input &rin)
{
	for ( int a=0; a<=rin.Get_bixmax(its_layer)-rin.Get_bixmin(its_layer); a++)
		for ( int b=0; b<=rin.Get_biymax(its_layer)-rin.Get_biymin(its_layer); b++)
			delete pca[a][b];
	for ( int a=0; a<=rin.Get_bixmax(its_layer)-rin.Get_bixmin(its_layer); a++)
		delete pca[a];
	delete[]pca;
	
}




/***************************************************************
*	set the values of the optical properties of the
*	voxels that conform a layer.
*
*	subject of further modifications.
****************************************************************/

void photon::SetCubeProperties( cube***&pca, input &rin)
{
	for ( int ix=rin.Get_bixmin(its_layer); ix<=rin.Get_bixmax(its_layer); ix++)
	{
	for ( int iy=rin.Get_biymin(its_layer); iy<=rin.Get_biymax(its_layer); iy++)
	{
	for ( int iz=0; iz<rin.Get_ncubesz(its_layer); iz++)
	{
	pca[ix-rin.Get_bixmin(its_layer)][iy-rin.Get_biymin(its_layer)][iz].Set_properties(rin);
	
	
	if (rin.Get_nspeclayers()!=0)
	{
	for ( int i=0; i<rin.Get_nspeclayers(); i++)
	{
	if ( its_layer==rin.Get_speclayer(i) )
	{
	for ( int j=0; j<rin.Get_nspeccubes(i); j++)
	{
	if ( ix==rin.Get_specix(i,j) )
	{
	if ( iy==rin.Get_speciy(i,j) )
	{
	if ( iz==rin.Get_speciz(i,j) )
	{
	
	
	pca[ix-rin.Get_bixmin(its_layer)][iy-rin.Get_biymin(its_layer)][iz].Set_n(rin.Get_specn(i,j));
	pca[ix-rin.Get_bixmin(its_layer)][iy-rin.Get_biymin(its_layer)][iz].Set_g(rin.Get_specg(i,j));
	pca[ix-rin.Get_bixmin(its_layer)][iy-rin.Get_biymin(its_layer)][iz].Set_mua(rin.Get_specmua(i, j));
	pca[ix-rin.Get_bixmin(its_layer)][iy-rin.Get_biymin(its_layer)][iz].Set_mus(rin.Get_specmus(i, j));
	}
	}
	}
	}
	}
	}
	}
	
	
	
	
	}				
	}
	}	
	
	
	

}




/*************************************************************************
* Set the labels of the boxes where the photon starts.
*************************************************************************/
void photon::Set_LaunchLabels(input &rin)
{
	double ix=(its_x/rin.Get_dx(its_layer));
	double iy=(its_y/rin.Get_dy(its_layer));
	double iz=((its_z-rin.Get_zlpos(its_layer))/rin.Get_dz(its_layer));
	its_ix=static_cast<int>(ix); 
  	its_iy=static_cast<int>(iy);
	its_iz=static_cast<int>(iz);
	if ( its_ix>its_ixmax )
		its_ix=its_ixmax;
	if ( its_ix<its_ixmin )
		its_ix=its_ixmin;
	if ( its_iy>its_iymax )
		its_iy=its_iymax;
	if ( its_iy<its_iymin )
		its_iy=its_iymin;
	

} 


/*************************************************************************
* Set the labels of the boxes where the photon cross a boundary.
*************************************************************************/
void photon::Set_LateralLabels(input &rin)
{
	double ix=(its_x/rin.Get_dx(its_layer));
	double iy=(its_y/rin.Get_dy(its_layer));
	
	its_ix=static_cast<int>(ix); 
  	its_iy=static_cast<int>(iy);
	
	if ( its_ix>its_ixmax )
		its_ix=its_ixmax;
	if ( its_ix<its_ixmin )
		its_ix=its_ixmin;
	if ( its_iy>its_iymax )
		its_iy=its_iymax;
	if ( its_iy<its_iymin )
		its_iy=its_iymin;
}
	
	

/*******************************************************************
*	These function generates a normal(gaussian) distribution of
*	random numbers.
********************************************************************/
double photon::normal(double mu, input &rin )
{
	double sigma=rin.Get_stand();
	assert ( sigma>0.);
	double p, p1, p2;
	
	do
	{
		p1=uniform(-1., 1. );
		p2=uniform(-1., 1. );
		p=p1*p1+p2*p2;
	}while(p>=1.);
	
	
	return mu+sigma*p1*sqrt(-2 * log(p)/p);
}

  
/*******************************************************************
*	These function generates a uniform distribution of
*	random numbers.
********************************************************************/

double photon::uniform ( double xMin=0., double xMax=1.)
{
	assert (xMin<xMax);
	return xMin+(xMax-xMin)*RandomNum();
}



/*********************************************************************
*	the following functions set the identification labels for 
*	the voxel the photon intends to go in, so they can be used
*	to decide if the photon enters a voxel with different optical 
*	properties
**********************************************************************/

void photon::Set_ixnew(input &rin)
{
	
	if(  its_ux<0.0 && its_ix!=its_ixmin)
	
		its_ix -= 1;
		
	else
	
		if (  its_ux>0.0 && its_ix!=its_ixmax)
		
			its_ix += 1;
			
		
	
	
}
void photon::Set_iynew(input &rin)
{
	
	if( its_uy<0.0 && its_iy!=its_iymin)
	
		its_iy-= 1;
		
	else
	
		if ( its_uy>0.0 && its_iy!=its_iymax)
		
			its_iy+=1;
		
}

void photon::Set_iznew(input &rin)
{
	
	if(  its_uz<0.0 )
			
		its_iz-=1;
		
	else
	{
		if (  its_uz>0.0 )
			its_iz+=1;
		else
			nrerror("\n The photon was supposed to change its iz label and it hasn't.");
		
	}
	

			
}
void photon::Set_NewCubeLabels( input &rin)
{
	if ( its_s==its_sx )
		Set_ixnew(rin);
	if ( its_s==its_sy )
		Set_iynew(rin);
	if ( its_s==its_sz )
		Set_iznew(rin);
	

}



/*****************************************************************
*	the next function returns the value of the radial 
*	position of the photon
*****************************************************************/
double photon::Get_r(input &rin)
{
	its_r=sqrt (  pow(its_x-(rin.Get_xmax()/2.0),2.0)+pow(its_y-(rin.Get_ymax()/2.0),2.0));
	
	return its_r;
}



/***********************************************************
 *      Compute the specular reflection. 
 *
 *      If the first layer is a turbid medium, use the Fresnel
 *      reflection from the boundary of the first layer as the 
 *      specular reflectance.
 *
 *      If the first layer is glass, multiple reflections in
 *      the first layer is considered to get the specular
 *      reflectance.
 *	
 *	It doesn't allow any special geometries on the cubes
 *	delimiting the layer.
 *
 *      
 ****/
double photon::Rspecular(input &rin)
{
	double r1, r2;	/* direct reflections from the 1st and 2nd layers. */
  	double temp;
  
  	temp =(rin.Get_n(0) - rin.Get_n(1))/(rin.Get_n(0) + rin.Get_n(1));
  	r1 = temp*temp;
  
  	if((rin.Get_mua(1) == 0.0) && (rin.Get_mus(1) == 0.0))   /* glass layer. */  
	{ 
    		temp = (rin.Get_n(1) - rin.Get_n(2))/(rin.Get_n(1) + rin.Get_n(2));
    		r2 = temp*temp;
    		r1 = r1 + (1-r1)*(1-r1)*r2/(1-r1*r2);
  	}
  
  	its_Rspecular=r1;
	return its_Rspecular; 
}





/***********************************************************
*	Initialize a photon packet.
***********************************************************/
void photon::LaunchPhoton(input &rin, cube ***&pca, output &rout)
{
	its_dead=0;
  	double a;
	a=Rspecular(rin);
	rout.Set_Rsra(a);
	Set_weight(1.0 - a);      
  	Set_layer(1);
	rin.Set_layer(1);
	Set_max(rin);
	Set_min(rin);
  	SetCubeArray( pca, rin);
	SetCubeProperties(pca, rin);
	its_s=0.0;
	its_sx=0.0;
	its_sy=0.0;
	its_sz=0.0;
	its_stotal=0.0;
	  
  	if (!rin.Get_specbeam())
	{
  		its_x= rin.Get_xo(); 
  		its_y= rin.Get_yo();
		
	}
	else
	{
		its_x=normal(rin.Get_xo(), rin);
		its_y=normal(rin.Get_yo(), rin);
		if ( its_x<rin.Get_xo()-(2*rin.Get_stand()) || its_x>rin.Get_xo()+(2*rin.Get_stand()) ||
			its_y<rin.Get_yo()-(2*rin.Get_stand()) || its_y>rin.Get_yo()+(2*rin.Get_stand()))
		{
			its_dead=1;
			DeleteCubeArray(pca, rin);
		}
	}
	
	
		
	
		 
  	its_z= rin.Get_zlpos(its_layer);
	its_ux= 0.0;  
  	its_uy= 0;  
  	its_uz= 1.0;
	Set_LaunchLabels(rin);
	
}


/***********************************************************
*      Pick a step size for a photon packet when it is in 
*      tissue.
*      If the member sleft is zero, make a new step size 
*      with: -log(rnd)/(mua+mus).
*      Otherwise, pick up the leftover in sleft.
*
*      Layer is the index to layer.
*      In_Ptr is the input parameters.
****/
void photon::StepSize(cube ***&pca, input &rin)
{
	double mua = pca[its_ix-its_ixmin][its_iy-its_iymin][its_iz].Get_mua();
  	double mus = pca[its_ix-its_ixmin][its_iy-its_iymin][its_iz].Get_mus();
  
  	if(its_stotal == 0.0) 
	{  /* make a new step. */
    		double rnd;

    		do
		{
			rnd = RandomNum();
		} 
      		while( rnd <= 0.0 );    /* avoid zero. */
        		its_stotal = -log(rnd)/(mua+mus);
  	}
  	else 
		nrerror("\nAsk to take a new step size when the anterior is not still zero.");
}



/***********************************************************
 *      Choose (sample) a new theta angle for photon propagation
 *      according to the anisotropy.
 *
 *      If anisotropy g is 0, then
 *              cos(theta) = 2*rand-1.
 *      otherwise
 *              sample according to the Henyey-Greenstein function.
 *
 *      Returns the cosine of the polar deflection angle theta.
 ****/
double photon::SpinTheta(double g)
{
	double cost;
  
  	if(g == 0.0) 
    		cost = 2*RandomNum() -1;
  	else 
	{
    		double temp = (1-g*g)/(1-g+2*g*RandomNum());
    		cost = (1+g*g - temp*temp)/(2*g);
        	if(cost < -1) 
			cost = -1;
        	else if(cost > 1) cost = 1;
  	}
	
  	return(cost);
}






/***********************************************************
 *      Choose a new direction for photon propagation by 
 *      sampling the polar deflection angle theta and the 
 *      azimuthal angle psi.
 *
 *      Note:
 *      theta: 0 - pi so sin(theta) is always positive 
 *      feel free to use sqrt() for cos(theta).
 * 
 *      psi:   0 - 2pi 
 *      for 0-pi  sin(psi) is + 
 *      for pi-2pi sin(psi) is - 
 ****/
void photon::Spin(input &rin, cube ***&pca)
{
	double cost, sint;    /* cosine and sine of the polar deflection angle theta. */
  	double cosp, sinp;    /* cosine and sine of the azimuthal angle psi. */
  	double ux = Get_ux();
  	double uy = Get_uy();
  	double uz = Get_uz();
  	double psi;
	
  	cost = SpinTheta(pca[its_ix-its_ixmin][its_iy-its_iymin][its_iz].Get_g());
  	sint = sqrt(1.0 - cost*cost); /* sqrt() is faster than sin(). */

  	psi = 2.0*PI*RandomNum(); /* spin psi 0-2pi. */
  	cosp = cos(psi);
  	if(psi<PI)
    		sinp = sqrt(1.0 - cosp*cosp);  /* sqrt() is faster than sin(). */
  	else
    		sinp = - sqrt(1.0 - cosp*cosp);     
  
  	if(fabs(uz) > COSZERO)           /* normal incident. */ 
	{     
    		its_ux = sint*cosp;
    		its_uy = sint*sinp;
    		its_uz = cost*SIGN(uz);     /* SIGN() is faster than division. */
  	}
  	else  		/* regular incident. */
	{               
    		double temp = sqrt(1.0 - uz*uz);
    		its_ux = sint*(ux*uz*cosp - uy*sinp)/temp + ux*cost;
    		its_uy = sint*(uy*uz*cosp + ux*sinp)/temp + uy*cost;
    		its_uz = -sint*cosp*temp + uz*cost;
  	}
	
	
}




/*********************************************************
*	Set the maximum values that the index to voxels
*	in the actual layer can take.
**********************************************************/
void photon::Set_max(input &rin)
{
	its_ixmax= rin.Get_bixmax(its_layer);
	
	its_iymax= rin.Get_biymax(its_layer);
	
	its_izmax= rin.Get_ncubesz(Get_layer())-1;
}




/*********************************************************
*	Set the minimum values that the index to voxels
*	in the actual layer can take.
**********************************************************/
void photon::Set_min(input &rin)
{
	its_ixmin=rin.Get_bixmin(its_layer);
	its_iymin=rin.Get_biymin(its_layer);
	its_izmin=0;
}


	

/***********************************************************
*	Move the photon to the end of the actual voxel
************************************************************/
void photon::Move( input &rin)
{
	double sx, sy, sz, s;
	
	if (its_ux!=0.0 )
	{
		if ( its_ix!=its_ixmax && its_ix!=its_ixmin )
		{
			if ( its_ux>0 )
				sx=(((its_ix+1)*rin.Get_dx(its_layer))-its_x)/its_ux;
			else
				sx=((its_ix*rin.Get_dx(its_layer))-its_x)/its_ux;
		}
		if ( its_ix==its_ixmax  &&  its_ux>0 )
			sx=its_stotal;
		if ( its_ix==its_ixmax  &&  its_ux<0 )
			sx=((its_ix*rin.Get_dx(its_layer))-its_x)/its_ux;
		if ( its_ix==its_ixmin  &&  its_ux<0 )
			sx=its_stotal;
		if ( its_ix==its_ixmin  &&  its_ux>0 )
			sx=(((its_ix+1)*rin.Get_dx(its_layer))-its_x)/its_ux;
	}
		
	if(its_ux==0) 
		sx=100000;
	
	sx=sx*SIGN(sx);
	
	its_sx=sx;
	
	if ( its_uy!=0.0)
	{
		if ( its_iy!=its_iymax && its_iy!=its_iymin )
		{
			if ( its_uy>0 )
				sy=(((its_iy+1)*rin.Get_dy(its_layer))-its_y)/its_uy;
			else
				sy=((its_iy*rin.Get_dy(its_layer))-its_y)/its_uy;
		}
		if ( its_iy==its_iymax  &&  its_uy>0 )
			sy=its_stotal;
		if ( its_iy==its_iymax  &&  its_uy<0 )
			sy=((its_iy*rin.Get_dy(its_layer))-its_y)/its_uy;
		if ( its_iy==its_iymin  &&  its_uy<0 )
			sy=its_stotal;
		if ( its_iy==its_iymin  &&  its_uy>0 )
			sy=(((its_iy+1)*rin.Get_dy(its_layer))-its_y)/its_uy;
	}
		
	if ( its_uy==0)
		sy=100000;
	
	sy=sy*SIGN(sy);
	
	its_sy=sy;
	
	if ( its_uz!=0.0)
		if ( its_uz>0 )
			sz=(((its_iz+1)*rin.Get_dz(its_layer))-its_z+rin.Get_zlpos(its_layer))/its_uz;
		else
			sz=((its_iz*rin.Get_dz(its_layer))+rin.Get_zlpos(its_layer)-its_z)/its_uz;
	else
		sz=100000;
	
	sz=sz*SIGN(sz);
	
	its_sz=sz;
	
	s=min(its_sx, its_sy);
	
	s=min(its_sz, s);
	
	
	its_s=s;
	
	if ( its_s>=its_stotal )
	{
		its_s=its_stotal;
		its_stotal=0.0;
	}
	else
		its_stotal-=its_s;
	
		
}









/*********************************************************************
*	Boolean HitBoundary: returns 1 if the photon hits boundary and 
*	0 otherwise
*********************************************************************/

Boolean photon::HitBoundary( input &rin)
{
	double uz=its_uz;
	double iz=its_iz;
	
	if ( uz>0.0 )
		if ( iz>its_izmax )
		{
			return 1;
		}
		else
			return 0;
	else
		if ( uz<0.0 )
			if ( iz<its_izmin )
				return 1;
			else
			{
				return 0;
			}
		else
		{
			return 0;	// in this case uz=0 and it will
					// never cross the surface
		}
	

}






/***********************************************************
 *      Decide whether the photon will be transmitted or 
 *      reflected on the upper boundary (uz<0) of the current 
 *      layer.
 *
 *      If "layer" is the first layer, the photon packet will 
 *      be partially transmitted and partially reflected if 
 *      PARTIALREFLECTION is set to 1,
 *      or the photon packet will be either transmitted or 
 *      reflected determined statistically if PARTIALREFLECTION 
 *      is set to 0.
 *
 *      Record the transmitted photon weight as reflection.  
 *
 *      If the "layer" is not the first layer and the photon 
 *      packet is transmitted, move the photon to "layer-1".
 *
 *      Update the photon parmameters.
 ****/
void photon::CrossUpOrNot(input &rin, output &rout, cube ***&pca)
{
	double uz = its_uz;	/* z directional cosine. */
	double ux=its_ux;	/* x directional cosine. */
	double uy=its_uy;	/* y directional cosine. */
  	double uz1;   /* cosines of transmission alpha. always */
                                /* positive. */
  	double r=0.0; /* reflectance */
  	short  layer = its_layer;
  	double ni = rin.Get_n(layer);
  	double nt = rin.Get_n(layer-1);
  	double mua1=rin.Get_mua(layer);
	double mus1=rin.Get_mus(layer);
	double mua2=rin.Get_mua(layer-1);
	double mus2=rin.Get_mus(layer-1);
	double mut1=mua1+mus1;		// interaction coefficient of the old layer.
	double mut2=mua2+mus2;		// interaction coefficient of the new layer.
  
  	
  	/* Get r. */
  	if( - uz <= rin.Get_coscrit0(layer))
	{
		r=1.0;
		
    		
	}                    /* total internal reflection. */
  	else
	{ 
		r = RFresnel(ni, nt, -uz, &uz1);
		
	}
  
	#if PARTIALREFLECTION
  	if(layer == 1 && r<1.0) 
	{     /* partially transmitted. */
    		Set_uz(-uz1);      /* transmitted photon. */
    		RecordR(r, rin, rout);
    		Set_uz(-uz);       /* reflected photon. */
		its_iznew=its_iz;
  	}             
  	else 
		if(RandomNum() > r) 
		{/* transmitted to layer-1. */
			DeleteCubeArray(pca, rin);
			layer--;
    			Set_layer(layer);
			rin.Set_layer(layer);
			Set_max(rin);
			Set_min(rin);
			SetCubeArray( pca, rin);
			SetCubeProperties(pca, rin);
			its_iz=Get_izmax();		// set new cube labels.
			Set_LateralLabels(rin);
			ux *= ni/nt;
			uy *= ni/nt;
    			Set_ux(ux);
    			Set_uy(uy);
    			Set_uz(-uz1);
			its_stotal*=(mut1/mut2);
			
  		}
  		else                                  /* reflected. */
		{
    			Set_uz (-uz);
			its_iz=its_izold;
			
		}
	#else
  	if(RandomNum() > r) 
	{         /* transmitted to layer-1. */
    		if(layer==1)  
		{
      			Set_uz( -uz1);
      			RecordR(0.0, rin, rout);	// add 0 because is totally transmitted.
      			Set_dead(1);
			DeleteCubeArray(pca, rin);
			
    		}
    		else 
		{
      			
			DeleteCubeArray(pca, rin);
			layer--;
			Set_layer(layer);
			rin.Set_layer(layer);
			Set_max(rin);
			Set_min(rin);
			SetCubeArray( pca, rin);
			SetCubeProperties(pca, rin);
			its_iz=Get_izmax();		// set new cube labels.
			Set_LateralLabels(rin);
			ux *= ni/nt;
			uy *= ni/nt;
      			Set_ux(ux );
      			Set_uy(uy );
      			Set_uz(-uz1);
			its_stotal*=(mut1/mut2);
    		}
  	}
  	else
	{                                          /* reflected. */
    		
		Set_uz(-uz);
		its_iz=its_izold;
		
	}
	#endif
	
}




/***********************************************************
 *      Decide whether the photon will be transmitted  or be 
 *      reflected on the bottom boundary (uz>0) of the current 
 *      layer.
 *
 *      If the photon is transmitted, move the photon to 
 *      "layer+1". If "layer" is the last layer, record the 
 *      transmitted weight as transmittance. See comments for 
 *      CrossUpOrNot.
 *
 *      Update the photon parmameters.
 ****/
void photon::CrossDnOrNot(input &rin, output &rout, cube ***&pca)
{
	double uz = its_uz; 	/* z directional cosine. */
	double ux=its_ux;	/* x directional cosine. */
	double uy=its_uy;	/* y directional cosine. */
  	
  	double uz1;   /* cosines of transmission alpha. */
  	double r=0.0; /* reflectance */
  	short  layer = Get_layer();
  	double ni = rin.Get_n(layer);
  	double nt = rin.Get_n(layer+1);
	double mua1=rin.Get_mua(layer);
	double mus1=rin.Get_mus(layer);
	double mua2=rin.Get_mua(layer+1);
	double mus2=rin.Get_mus(layer+1);
	double mut1=mua1+mus1;		// interaction coefficient of the old layer.
	double mut2=mua2+mus2;		// interaction coefficient of the new layer.
  	/* Get r. */
  	if( uz <= rin.Get_coscrit1(layer)) 
    		r=1.0;              /* total internal reflection. */
  	else 
		r = RFresnel(ni, nt, uz, &uz1);
		
	#if PARTIALREFLECTION   
  	if(layer == rin.Get_nlayers()-2 && r<1.0) 
	{
    		Set_uz(uz1);
    		RecordT(r, rin, rout);
    		Set_uz(-uz);
		
  	}
  	else 
		if(RandomNum() > r) 
		{	/* transmitted to layer+1. */
			DeleteCubeArray(pca, rin);
			layer++;
    			Set_layer(layer);
    			rin.Set_layer(layer);
			Set_max(rin);
			Set_min(rin);
			SetCubeArray( pca, rin);
			SetCubeProperties(pca, rin);
			its_iz=Get_izmin();		// set new cube labels.
			Set_LateralLabels(rin);
			ux *= ni/nt;
			uy *= ni/nt;
    			Set_ux(ux);
    			Set_uy(uy );
    			Set_uz(uz1);
			its_stotal*=(mut1/mut2);
  		}
  		else                                          /* reflected. */
		{
    			Set_uz(-uz);
			its_iz=its_izold;
			
		}
	#else
	if(RandomNum() > r) 
	{         /* transmitted to layer+1. */
		if(layer == rin.Get_nlayers()-2) 
		{
      			Set_uz(uz1);
      			RecordT(0.0,rin, rout);		// add 0 because is totally transmitted.
			Set_dead(1);
			DeleteCubeArray(pca, rin);
			
    		}
    		else 
		{
      			DeleteCubeArray(pca, rin);
			layer++;
			Set_layer(layer);
      			rin.Set_layer(layer);
			Set_max(rin);
			Set_min(rin);
			SetCubeArray( pca, rin);
			SetCubeProperties(pca, rin);
			its_iz=Get_izmin();		// set new cube labels.
			Set_LateralLabels(rin);
			ux *= ni/nt;
			uy *= ni/nt;
      			Set_ux(ux );
			Set_uy(uy );
      			Set_uz(uz1);
			its_stotal*=(mut1/mut2);
			
    		}
  	}
  	else
	{                                          /* reflected. */
    		Set_uz(-uz);
		its_iz=its_izold;
		
		
	}
	#endif
	
}




/***********************************************************
*	Once the photon has hit a boundary interface
*	between layers, it differenciates if the photon
*	is moving forward or backwards and call the
*	appropiate function.
*************************************************************/
void photon::CrossLayerOrNot(input &rin, output &rout, cube ***&pca)
{
	
	double uz=its_uz;
	
	if(uz < 0.0)
	{
		CrossUpOrNot(rin, rout, pca);
	}
  	else
		if ( uz>0.0 )
    			CrossDnOrNot(rin, rout, pca);
		else
			nrerror("Error, trying to cross layer with zero zeta directional cosine.\n");
	
}










/***********************************************************
 *      Compute the Fresnel reflectance.
 *
 *      Make sure that the cosine of the incident angle a1
 *      is positive, and the case when the angle is greater 
 *      than the critical angle is ruled out.
 *
 *      Avoid trigonometric function operations as much as
 *      possible, because they are computation-intensive.
 ****/
double photon::RFresnel(double n1,/* incident refractive index.*/double n2,/* transmit refractive index.*/
                                double ca1,/* cosine of the incident *//* angle. 0<a1<90 degrees. */
                                double * ca2_Ptr)/* pointer to the */
                                                        /* cosine of the transmission */
                                                        /* angle. a2>0. */
{
	double r;
  
  	if(n1==n2) 
	{                          /** matched boundary. **/
    		*ca2_Ptr = ca1;
    		r = 0.0;
  	}
  	else 
		if(ca1>COSZERO) 
		{        /** normal incident. **/
    			*ca2_Ptr = ca1;
    			r = (n2-n1)/(n2+n1);
    			r *= r;
  		}
  		else 
			if(ca1<COS90D)  
			{        /** very slant. **/
    				*ca2_Ptr = 0.0;
    				r = 1.0;
  			}
  			else  
			{                                       /** general. **/
    				double sa1, sa2;    
          			/* sine of the incident and transmission angles. */
    				double ca2;
    
    				sa1 = sqrt(1-ca1*ca1);
    				sa2 = n1*sa1/n2;
    				if(sa2>=1.0) 
				{      
          				/* double check for total internal reflection. */
      					*ca2_Ptr = 0.0;
      					r = 1.0;
    				}
    				else  
				{
      					double cap, cam;  /* cosines of the sum ap or */
                                                /* difference am of the two */
                                                /* angles. ap = a1+a2 */
                                                /* am = a1 - a2. */
      					double sap, sam;  /* sines. */
      
      					*ca2_Ptr = ca2 = sqrt(1-sa2*sa2);
      
      					cap = ca1*ca2 - sa1*sa2; /* c+ = cc - ss. */
      					cam = ca1*ca2 + sa1*sa2; /* c- = cc + ss. */
      					sap = sa1*ca2 + ca1*sa2; /* s+ = sc + cs. */
      					sam = sa1*ca2 - ca1*sa2; /* s- = sc - cs. */
      					r = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam); 
                			/* rearranged for speed. */
    				}
  			}
  	return(r);
}












/***********************************************************
 *      Record the photon weight exiting the first layer(uz<0), 
 *      no matter whether the layer is glass or not, to the 
 *      reflection array.
 *
 *      Update the photon weight as well.
 ****/
void photon::RecordR(double Refl,input &rin, output &rout)
{
	double x = its_x;
  	double y = its_y;
  	int  ir, ia, ix, iy;        /* index to r & angle. */
	
	double w=Get_weight();	//weight of the photon packet
  
  	
	ix=static_cast<int>(x/rin.Get_gridx());
	if ( ix>rout.Get_nx()-1)
		ix=rout.Get_nx()-1;
	
	iy=static_cast<int>(y/rin.Get_gridy());
	if ( iy>rout.Get_ny()-1)
		iy=rout.Get_ny()-1;
	
	
	
	
  	/* assign photon to the reflection array element. */
	
	if ( its_x>=rin.Get_xmin()  &&  its_x<=rin.Get_xmax() && its_y>=rin.Get_ymin()  &&  its_y<=rin.Get_ymax() )
	{
		ix=static_cast<int>((x-rin.Get_xmin())/rin.Get_gridx());
		if ( ix>rout.Get_nx()-1)
			ix=rout.Get_nx()-1;
	
		iy=static_cast<int>((y-rin.Get_ymin())/rin.Get_gridy());
		if ( iy>rout.Get_ny()-1)
			iy=rout.Get_ny()-1;
	
	
  		rout.Add_to_Rdxy_element(w*(1.0-Refl),ix,iy);
	}
	else
	{
		rout.Add_to_Out(w*(1.0-Refl));
		rout.Add_to_ROut(w*(1.0-Refl));
	}
  	
	
	
  	w *= Refl;			// actualize weight
	Set_weight(w);
  	
	
	
}





/***********************************************************
 *      Record the photon weight exiting the last layer(uz>0), 
 *      no matter whether the layer is glass or not, to the 
 *      transmittance array.
 *
 *      Update the photon weight as well.
 ****/
void photon::RecordT(double Refl,input &rin, output &rout)
{
	int  ir, ia, ix, iy;        /* index to r & angle. */
	
	double w=its_weight;	// weight of the photon packet
  
  	ix=static_cast<int>(its_x/rin.Get_gridx());
	if ( ix>rout.Get_nx()-1)
		ix=rout.Get_nx()-1;
	
	iy=static_cast<int>(its_y/rin.Get_gridy());
	if ( iy>rout.Get_ny()-1)
		iy=rout.Get_ny()-1;
	
	
  	/* assign photon to the transmittance array element. */
	
	if ( its_x>=rin.Get_xmin()  &&  its_x<=rin.Get_xmax() && its_y>=rin.Get_ymin()  &&  its_y<=rin.Get_ymax() )
		rout.Add_to_Tdxy_element(w*(1.0-Refl),ix,iy);
	else
	{
		rout.Add_to_Out(w*(1.0-Refl));
  		rout.Add_to_TOut(w*(1.0-Refl));
	}
	
	
  	w *= Refl;		// actualize weight.
	Set_weight(w);
  	
	
	
	
}





/********************************************************************
*	Determines if the photon finds in its way a voxel
*	with different refraction index.
*	Return 1 if is different and 0 otherwise.
*********************************************************************/
Boolean photon::HitDifferentCube(input &rin, cube ***&pca)
{
	
	int ix=its_ix;
	int iy=its_iy;
	int iz=its_iz;
	int ixold=its_ixold;
	int iyold=its_iyold;
	int izold=its_izold;
	if ( (pca[ix-its_ixmin][iy-its_iymin][iz].Get_n() !=
	 	pca[ixold-its_ixmin][iyold-its_iymin][izold].Get_n()) ||
		 (pca[ix-its_ixmin][iy-its_iymin][iz].Get_g() != 
		 pca[ixold-its_ixmin][iyold-its_iymin][izold].Get_g()) ||
		 (pca[ix-its_ixmin][iy-its_iymin][iz].Get_mua() != 
		 pca[ixold-its_ixmin][iyold-its_iymin][izold].Get_mua()) ||
		 (pca[ix-its_ixmin][iy-its_iymin][iz].Get_mus() != 
		 pca[ixold-its_ixmin][iyold-its_iymin][izold].Get_mus()))
		return 1;
	else
		return 0;
}



	


/***********************************************************
*      Decide whether the photon will be transmitted back or 
*      reflected forward on the boundary with a voxel with 
*      different refraction index.
*
*      Update the photon parmameters.
****/
void photon::CrossDfCbUpOrNot(input &rin, output &rout, cube ***&pca, int a)
{
	double u;	/* directional cosine to choose between ux, uy, uz. */
	double u1, u2;	// the other two directional cosines.
	if ( a==1 )	// the x direction is the crossing one.
		u=its_ux;
		u1=its_uy;	
		u2=its_uz;	
  	if ( a==2 )	// the y direction is the crossing one.
		u=its_uy;
		u1=its_uz;	
		u2=its_ux;	
  	if ( a==3 )	// the z direction is the crossing one.
		u=its_uz;
		u1=its_ux;	
		u2=its_uy;	
	
	double uu;   /* cosines of transmission alpha. always */
                                /* positive. */
  	double r=0.0; /* reflectance */
  	short  layer = Get_layer();
  	double ni = pca[its_ixold-its_ixmin][its_iyold-its_iymin][its_izold].Get_n();
  	double nt = pca[its_ix-its_ixmin][its_iy-its_iymin][its_iz].Get_n();
  
  	/* Get r. */
  	/*if( - u <= rin.Get_coscrit0(Get_layer()))	// have to modificate this 
    		r=1.0; */                   /* total internal reflection. */
  	//else 
		r = RFresnel(ni, nt, -u, &uu);
  
	if(RandomNum() > r) 
	{         
    		if ( a==3 )
		{
			u1 *= ni/nt;
			u2 *= ni/nt;
      			Set_ux(u1 );
      			Set_uy(u2 );
      			Set_uz(-uu);
		}
		if ( a== 1)
		{
			u1 *= ni/nt;
			u2 *= ni/nt;
			Set_uy(u1);
			Set_uz(u2);
			Set_ux(-uu);
		}
		if ( a==2 )
		{
			u1 *= ni/nt;
			u2 *= ni/nt;
			Set_uz(u1);
			Set_ux(u2);
			Set_uy(-uu);
		}
			
    		
  	}
  	else 
	{				/* reflected. */
		if ( a==3 )
		{
			its_iz=its_izold;                                         
    			Set_uz(-u);
		}
		if ( a==1 )
		{
			its_ix=its_ixold;                                         
    			Set_ux(-u);
		}
		if ( a==2 )
		{
			its_iy=its_iyold;  
			Set_uy(-u);
		}
	}			
	
}




/***********************************************************
*      Decide whether the photon will be transmitted forward or 
*      reflected back on the boundary with a voxel with 
*      different refraction index.
*
*      Update the photon parmameters.
***********************************************************/
void photon::CrossDfCbDnOrNot(input &rin, output &rout,cube ***&pca, int a)
{
	double u;	/* directional cosine to choose between ux, uy, uz. */
	double u1, u2;	// the other two directional cosines.
	if ( a==1 )	// the x direction is the crossing one.
		u=its_ux;
		u1=its_uy;	
		u2=its_uz;	
  	if ( a==2 )	// the y direction is the crossing one.
		u=its_uy;
		u1=its_uz;	
		u2=its_ux;	
  	if ( a==3 )	// the z direction is the crossing one.
		u=its_uz;
		u1=its_ux;	
		u2=its_uy;	
	
	
	
  	double uu;   /* cosines of transmission alpha. */
  	double r=0.0; /* reflectance */
  	
	double ni = pca[its_ixold-its_ixmin][its_iyold-its_iymin][its_izold].Get_n();
  	double nt = pca[its_ix-its_ixmin][its_iy-its_iymin][its_iz].Get_n();
  
	
  	/* Get r. */
  	/*if( u <= rin.Get_coscrit1(its_layer))	// have to modificate this 
    		r=1.0;*/              /* total internal reflection. */
  	//else 
	
	
		r = RFresnel(ni, nt, u, &uu);
  
 	if(RandomNum() > r) 
	{ 
		if ( a==3 )
		{
      			u1 *= ni/nt;
			u2 *= ni/nt;
      			Set_ux(u1);
			Set_uy(u2);
      			Set_uz(uu);
    		}
		if ( a==1 )
		{
			u1 *= ni/nt;
			u2 *= ni/nt;
			Set_uy(u1);
			Set_uz(u2);
			Set_ux(uu);
		}
		if ( a==2 )
		{
			u1 *= ni/nt;
			u2 *= ni/nt;
			Set_uz(u1);
			Set_ux(u2);
			Set_uy(uu);
		}
		
  	}
  	else
	{                                          /* reflected. */
    		if ( a==3 )
		{
			Set_uz(-u);
			its_iz=its_izold;
		}
		if ( a==1 )
		{
			Set_ux(-u);
			its_ix=its_ixold;
		}
		if ( a==2 )
		{
			Set_uy(-u);
			its_iy=its_iyold;
		}
		
	}	


}









/***********************************************************
*	Selects in which direction the photon is crossing
*	and calls the appropiate function.
***********************************************************/
void photon::CrossDfCbOrNot(input &rin, output &rout, cube ***&pca)
{
	
	if ( its_ix!=its_ixold )
	{
		
		if(its_ux < 0.0)
    			CrossDfCbUpOrNot(rin, rout, pca, 1);
  		else
    			CrossDfCbDnOrNot(rin, rout, pca, 1);
	}
	else
	{
	
		if ( its_iy!=its_iyold )
		{
		
			if(its_uy < 0.0)
    				CrossDfCbUpOrNot(rin, rout, pca, 2);
  			else
    				CrossDfCbDnOrNot(rin, rout, pca, 2);
		}
		else
		{
	
			if ( its_iz!=its_izold )
			{
		
				if(its_uz < 0.0)
    					CrossDfCbUpOrNot(rin, rout, pca, 3);
  				else
    					CrossDfCbDnOrNot(rin, rout, pca, 3);
			}
			else
				nrerror("\nIt has crossed to a cube with different opt prop but something is wrong.");
		}
	}
}



/***************************************************************
*	Actualize the photon positions.
****************************************************************/
void photon::Hop(input &rin)
{
	its_xold=its_x;
	its_yold=its_y;
	its_zold=its_z;
	its_ixold=its_ix;
	its_iyold=its_iy;
	its_izold=its_iz;
	
	its_x += (its_s*its_ux);
	its_y += (its_s*its_uy);
	its_z += (its_s*its_uz);
	
	if ( its_stotal>=0.0 )
		Set_NewCubeLabels(rin);
		
	its_s=0.0;
	its_sx=0.0;
	its_sy=0.0;
	its_sz=0.0;
	
}
	




/********************************************************************
*	Selects a step size for the photon and keeps moving
*	it until it is dead or it covers the assigned
*	path length.
*********************************************************************/

void photon::StepMoveHopeCheck(input &rin, output &rout, cube ***&pca)
{
	StepSize(pca, rin);
	do
	{
		Move(rin);
		Hop(rin);
		if ( !OutOfBoundaries(rin, pca, rout)  && its_dead==0 )
			if ( HitBoundary(rin) )
				CrossLayerOrNot(rin, rout, pca );
			else
				if (HitDifferentCube(rin, pca) )
					CrossDfCbOrNot(rin, rout, pca);
	}while( its_dead==0  &&  its_stotal>0 );	
	
}



/******************************************************************************
*	Combines the movement of the photon with the possibility of
*	scattering events.
*
*	If the photon weight is too small it has to survive a 
*	roulette to be still alive.
*******************************************************************************/
void photon::StepMoveHopeCheckDropSpin(input &rin, output &rout, cube ***&pca)
{
	do
	{
		StepMoveHopeCheck(rin, rout, pca);
		if ( its_dead==0 )
		{
			Drop(rin, rout, pca);
			Spin(rin, pca);
	
			if ( its_weight<WEIGHT && its_dead==0 )
				Roulette(pca, rin);
		}
	}while(its_dead==0);
}








/***********************************************************
 *      Drop photon weight inside the tissue.
 *
 *  	The photon is assumed not dead. 
 *
 *      The weight drop is dw = w*mua/(mua+mus).
 *
 *      The dropped weight is assigned to the absorption array 
 *      elements.
 ****/
void photon::Drop(input &rin, output &rout, cube ***&pca)
{
	double dwa;           /* absorbed weight.*/
  	double x=its_x;
  	double y=its_y;
	double z=its_z;
	int ix=its_ix;
	int iy=its_iy;
	int iz=its_iz;
	int  layer ;
  	double mua, mus;              
  

  	/* update photon weight. */
  	mua = pca[ix-its_ixmin][iy-its_iymin][iz].Get_mua();
	
  	mus = pca[ix-its_ixmin][iy-its_iymin][iz].Get_mus();
	
  	dwa = its_weight * mua/(mua+mus);
	
  	its_weight -= dwa;
	
	unsigned int ixx; 
	unsigned int iyy; 
	unsigned int izz; 
	
		
  	/* assign dwa to the absorption array element. */
	
	if ( its_x>=rin.Get_xmin()  &&  its_x<=rin.Get_xmax() && its_y>=rin.Get_ymin()  &&  its_y<=rin.Get_ymax() )
	{
		ixx=static_cast<int>( (x-rin.Get_xmin())/rin.Get_gridx() );
  		iyy=static_cast<int>( (y-rin.Get_ymin())/rin.Get_gridy() );
   		izz=static_cast<int>( z/rin.Get_gridz() );
	
  		if ( ixx>rout.Get_nx()-1)
			ixx=rout.Get_nx()-1;
		
		if ( iyy>rout.Get_ny()-1)
			iyy=rout.Get_ny()-1;
		
		if ( izz>rout.Get_nz()-1)
			izz=rout.Get_nz()-1;
	
		rout.Add_to_Axyz_element(dwa ,ixx ,iyy ,izz);
	}
	else
	{
		rout.Add_to_AOut(dwa);
		rout.Add_to_Out(dwa);
	} 
}







/*****************************************************************
 *      The photon weight is small, and the photon packet tries 
 *      to survive a roulette.
 *****************************************************************/
void photon::Roulette(cube ***&pca, input &rin)
{
  	double weight=its_weight;
	
	if( weight == 0.0)
	{      
    		Set_dead(1);
		DeleteCubeArray(pca, rin);
	}
  	else 
		if(RandomNum() < CHANCE) /* survived the roulette.*/
    			its_weight /= CHANCE;
  		else
		{ 
    			its_dead=1;
			DeleteCubeArray(pca, rin);
		}
}




/***********************************************************
*	Returns 1 if the photon has escaped the phantom
*	and 0 otherwise.
*
*	Kills the photon if its scapes the boundaries edges
*	in the x and y directions.
************************************************************/	
Boolean photon::OutOfBoundaries( input &rin, cube ***&pca, output &rout)
{
	double ix=its_ix;
	double ixmax=its_ixmax;
	double ixmin=its_ixmin;
	double iy=its_iy;
	double iymax=its_iymax;
	double iymin=its_iymin;
	
	if ( (ix<iymin) || (ix>ixmax ) || (iy<iymin) || (iy>iymax )  )
	{	
		nrerror("\nThe photon has escaped and it shouldn't have.");
		return 1;
	}
	else
		return 0;
		
		
	
	
	
}




/***********************************************************
*	Do one simulation run.
*************************************************************/
void photon::DoOneRun(input &rin, output &rout, cube ***&pca)
{
	rout.InitOutputData(rin);
	long i_photon=rin.Get_nphotons();
	
	do
	{	
		LaunchPhoton(rin, pca, rout);
		if ( its_dead==0 )
		{
			StepMoveHopeCheckDropSpin(rin, rout, pca);
			//cout<<i_photon<<endl;
		}
		
	}while( --i_photon);
}
	




/*************************************************************
*	Random number generator from Numerical Recipes in C.
**************************************************************/
float photon::ran3(int *idum)	
{
	static int inext,inextp;
  	static long ma[56];
  	static int iff=0;
  	long mj,mk;
  	int i,ii,k;
  
  	if (*idum < 0 || iff == 0) 
	{
    		iff=1;
    		mj=MSEED-(*idum < 0 ? -*idum : *idum);
    		mj %= MBIG;
    		ma[55]=mj;
    		mk=1;
    		for (i=1;i<=54;i++) 
		{
      			ii=(21*i) % 55;
      			ma[ii]=mk;
      			mk=mj-mk;
      			if (mk < MZ) 
				mk += MBIG;
      			mj=ma[ii];
    		}
    		for (k=1;k<=4;k++)
      			for (i=1;i<=55;i++) 
			{
        			ma[i] -= ma[1+(i+30) % 55];
        			if (ma[i] < MZ)
					 ma[i] += MBIG;
      			}
    		inext=0;
    		inextp=31;
    		*idum=1;
  	}
  	if (++inext == 56)
		inext=1;
  	if (++inextp == 56) 
		inextp=1;
  	mj=ma[inext]-ma[inextp];
  	if (mj < MZ) 
		mj += MBIG;
  	ma[inext]=mj;
  return mj*FAC;
}



double photon::RandomNum(void)
{
	static Boolean first_time=1;
  	static int idum;      /* seed for ran3. */
  
  	if(first_time) 
	{
		#if STANDARDTEST /* Use fixed seed to test the program. */
    		idum = - 1;
		#else
    		idum = -(int)time(NULL)%(1<<15);
          		/* use 16-bit integer as the seed. */
		#endif
    		ran3(&idum);
    		first_time = 0;
    		idum = 1;
  	}
  
  	return( (double)ran3(&idum) );
}


/******************************************************
* Function use to check values during run time. 
* It can be place at any part of the program.
******************************************************/
void photon::ListPhValues()
{
	cout<<"\nits_x:\t"<<its_x<<endl;
	cout<<"\nits_y:\t"<<its_y<<endl;
	cout<<"\nits_z:\t"<<its_z<<endl;
	cout<<"\nits_ix:\t"<<its_ix<<endl;
	cout<<"\nits_iy:\t"<<its_iy<<endl;
	cout<<"\nits_iz:\t"<<its_iz<<endl;
	cout<<"\nits_ux:\t"<<its_ux<<endl;
	cout<<"\nits_uy:\t"<<its_uy<<endl;
	cout<<"\nits_uz:\t"<<its_uz<<endl;
	cout<<"\nits_weight:\t"<<its_weight<<endl;
	cout<<"\nits_layer:\t"<<its_layer<<endl;
	
	
}
	
	
	
	
	
	
	
	
	
	
	
	

	
	
