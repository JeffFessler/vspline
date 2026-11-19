/***************************** mexdep.c **********************************
* SUBROUTINES WHICH DEPEND ON WHETHER OR NOT MEX FILES ARE USED.
* USES fprintf(), calloc(), free()
* Jeff Fessler
* COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
* LELAND STANFORD JUNIOR UNIVERSITY. ALL RIGHTS RESERVED.
*************************************************************************/
#include "mexdep.h"


/*
* ALLOCATION WITH CHECKING
*/
char *alloc_check(cint n1, cint n2)
{
	char *p;
	if ( !(p = mycalloc(n1, n2)) )
		Exit("Can't allocate memory")
	return p;
}
