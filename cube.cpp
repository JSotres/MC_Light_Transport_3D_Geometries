/*********************************************************
*	This file includes the functions of the
*	cube class.
**********************************************************/
#include "transport.h"



/***********************************************************
*      Report error message to stderr, then exit the program
*      with signal 1.
****/
void cube::nrerror(char error_text[])
     
{
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}




void cube::Set_properties(input &rin )
{
	its_dx=rin.Get_dx(rin.Get_layer());
	its_dy=rin.Get_dy(rin.Get_layer());
	its_dz=rin.Get_dz(rin.Get_layer());
	its_n=rin.Get_n(rin.Get_layer());
	its_g=rin.Get_g(rin.Get_layer());
	its_mua=rin.Get_mua(rin.Get_layer());
	its_mus=rin.Get_mus(rin.Get_layer());
	
	
}





