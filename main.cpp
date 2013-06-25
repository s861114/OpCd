#include <stdio.h>
#include <stdlib.h>
#include "wksp.h"
#include "function.h"

int main(void)
{
	int N,N2;
	int Nlayer[8]={1,2,3,3,4,4,5,5};
	
	for(int indc=1;indc<2;indc++){
		WKSP wksp;	
		N=Nlayer[indc];
		N2=N*2;
		wksp.N=N;
		wksp.N2=N2;
		wksp.initial_malloc(N,N2);
		wksp.slfcssnt_Ef(N,N2);
		double Ef=wksp.frmlvl_skr(N,N2,2e+16);
		printf("Ef=%f",Ef);
		//	wksp.opdc();	
		
	}

	return 0;
}
