#include <stdio.h>
#include <stdlib.h>
#include "wksp.h"
#include "function.h"

int main(void)
{
	WKSP wksp;
	wksp.set_H_AB();
	wksp.band_cal();	
	double Ef=wksp.frmlvl_skr(2e+16);
	printf("Ef=%f",Ef);
//	wksp.opdc();	
	

	return 0;
}
