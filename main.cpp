#include <stdio.h>
#include <stdlib.h>
#include "wksp.h"
#include "function.h"

int main(void)
{
	WKSP wksp;
	wksp.set_H_AB();
	wksp.band_cal();
	wksp.frmlvl_skr();
//	wksp.opdc();	
	

	return 0;
}
