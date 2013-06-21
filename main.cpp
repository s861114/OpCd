#include <stdio.h>
#include <stdlib.h>
#include "wksp.h"
#include "function.h"

int main(void)
{
	WKSP wksp;
	
//	int indc=1;
//	int	indc_nlayer_arr[8]={1,2,3,3,4,4,5,5};
//	wksp.N=indc_nlayer_arr[indc];
//	wksp.N2=wksp.N*2;
//	printf("%d %d \n",wksp.N,wksp.N2);

	//Here wksp.N2=4, wksp.N=2;

	wksp.initial_malloc();
	wksp.set_H_AB();
//	wksp.diag_slfcssnt(2,0,0,0,0);
	
	return 0;
}
