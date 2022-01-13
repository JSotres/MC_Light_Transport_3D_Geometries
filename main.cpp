/************************************************
* main function
*************************************************/
#include "transport.h"







int main()
{
	input in;
	
	photon ph;
	
	cube ***pcarray;
	
	output out;
	
	in.ReadParam();
	
	out.SetFile();
	
	ph.DoOneRun(in, out, pcarray);
	
	out.SumScaleResult(in);		
	
	out.FreeMemory(in);
	
	in.DeleteInputPointers();
	
	cout<<"\nnph: "<<in.Get_nphotons()<<endl;
	
	return 0;
}





