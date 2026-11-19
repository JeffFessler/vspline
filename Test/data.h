/***************************** data.h ************************************
* DECLARATIONS FOR DATA THAT IS READ FROM A FILE
* Jeff Fessler
* COPYRIGHT 1989, THE BOARD OF TRUSTEES OF THE
* LELAND STANFORD JUNIOR UNIVERSITY.  ALL RIGHTS RESERVED.
*************************************************************************/

typedef struct
{
	ALPHA_SET alpha;
	VECT ti;
	CMAT yi;
	COVAR *covar;
	int M, N;
	truf flag_iid, flag_equ;
} FILEDATA;
