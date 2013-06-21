#include <stdio.h>
#include <stdlib.h>
#include "wksp.h"
#include "function.h"

int main(void)
{
	WKSP wksp;
	//wksp.set_H_A();
	wksp.diag_slfcssnt(3,0,0,0,0);
	
	return 0;
}
