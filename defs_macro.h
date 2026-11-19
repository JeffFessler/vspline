// defs_macro.h
// favorite macro definitions
// Copyright 1992, Jeff Fessler, University of Michigan

/*
* Some include files vary between systems
*/

#ifndef DEFS_MACRO
#define DEFS_MACRO

#define BZERO(p,n)	bzero((char *) (p), (int) ((n) * sizeof(*(p))))
#define BCOPY(s,d,n)	bcopy((char *) (s), (char *) (d), (int) ((n) * sizeof(*(s))))
#define	CALLOC(n,s)	calloc((unsigned) (n), (unsigned) (s))
#define	MALLOC(n,s)	malloc((unsigned) ((n)*(s)))
#define	ALLOC(n,s)	CALLOC(n,s)

#define	SWAB(s,d,n)	swab((char *)(s), (char *)(d), n)

#define	NONULL(str)	((str) ? (str) : "null")

#define Argflt(n,def)	argc > n ? atof(argv[n]) : def
#define Argint(n,def)	argc > n ? atoi(argv[n]) : def
#define Argstr(n,def)	argc > n ? argv[n] : def
#define Argstn(n)	(argc > n && argv[n][0] != '-') ? argv[n] : NULL

#ifndef ERROR_MSG
#	define ERROR_MSG(msg)	{ fputs(msg, stderr); fputs("\n", stderr); }
#endif


#ifndef NOTE
#	define NOTE(msg)	{ char Zstr[1024];	\
	sprintf(Zstr, "NOTE in %s line %d:", __FILE__, __LINE__);	\
	ERROR_MSG(Zstr) ERROR_MSG(msg) }
#endif

#ifndef NOTE1
#	define NOTE1(msg, arg)	{ char Zstr[1024];	\
	sprintf(Zstr, "NOTE in %s line %d:", __FILE__, __LINE__);	\
	ERROR_MSG(Zstr)		\
	sprintf(Zstr, msg, arg);	\
	ERROR_MSG(Zstr) }
#endif

#ifndef NOTE2
#	define NOTE2(msg, arg1, arg2)	{ char Zstr[1024];	\
	sprintf(Zstr, "NOTE in %s line %d:", __FILE__, __LINE__);	\
	ERROR_MSG(Zstr)		\
	sprintf(Zstr, msg, arg1, arg2);	\
	ERROR_MSG(Zstr) }
#endif


#ifndef WARN
#	define WARN(msg)	{ char Zstr[1024];	\
	sprintf(Zstr, "WARNING in %s line %d:", __FILE__, __LINE__);	\
	ERROR_MSG(Zstr) ERROR_MSG(msg) }
#endif

#ifndef WARN1
#	define WARN1(msg, arg)	{ char Zstr[1024];	\
	sprintf(Zstr, "WARNING in %s line %d:", __FILE__, __LINE__);	\
	ERROR_MSG(Zstr)		\
	sprintf(Zstr, msg, arg);	\
	ERROR_MSG(Zstr) }
#endif

#ifndef WARN2
#	define WARN2(msg, arg1, arg2)	{ char Zstr[1024];	\
	sprintf(Zstr, "WARNING in %s line %d:", __FILE__, __LINE__);	\
	ERROR_MSG(Zstr)		\
	sprintf(Zstr, msg, arg1, arg2);	\
	ERROR_MSG(Zstr) }
#endif

#ifndef FAIL
#	define FAIL(msg)	{ char Zstr[1024]; \
	sprintf(Zstr, "ERROR in %s line %d:", __FILE__, __LINE__); \
	ERROR_MSG(Zstr) ERROR_MSG(msg) return FAILURE;}
#endif

#ifndef FAIL1
#	define FAIL1(msg, arg)	{ char Zstr[1024]; \
	sprintf(Zstr, "ERROR in %s line %d:", __FILE__, __LINE__); \
	ERROR_MSG(Zstr)	\
	sprintf(Zstr, msg, arg);	\
	ERROR_MSG(Zstr) return FAILURE; }
#endif

#ifndef FAIL2
#	define FAIL2(msg, arg1, arg2)	{ char Zstr[1024]; \
	sprintf(Zstr, "ERROR in %s line %d:", __FILE__, __LINE__); \
	ERROR_MSG(Zstr)	\
	sprintf(Zstr, msg, arg1, arg2);	\
	ERROR_MSG(Zstr) return FAILURE; }
#endif

#ifndef EXIT
#	define EXIT(msg)	{ char Zstr[1024]; \
	sprintf(Zstr, "ERROR in %s line %d", __FILE__, __LINE__); \
	ERROR_MSG(Zstr) ERROR_MSG(msg) exit(-1);}
#endif


/*
* And those pesky IO functions....
*/
#if defined(titan)
#define Fread(p,s,n,f)		fread ((char *)(p), (SIZE_T)(s), (int)(n), f)
#define Fwrite(p,s,n,f)		fwrite((char *)(p), (SIZE_T)(s), (int)(n), f)
#else
# if defined(sun)
# define Fread(p,s,n,f)		fread ((char *)(p), (int)(s), (int)(n), f)
# define Fwrite(p,s,n,f)	fwrite((char *)(p), (int)(s), (int)(n), f)
# else
# define Fread(p,s,n,f)		fread ((void *)(p), (SIZE_T)(s), (SIZE_T)(n), f)
# define Fwrite(p,s,n,f)	fwrite((void *)(p), (SIZE_T)(s), (SIZE_T)(n), f)
# define MEMSET(p,c,n)		memset((void *)(p), (int)(c), (SIZE_T)(n))
# define MEMCPY(t,f,n)		memcpy((void *)(t), (void *)(f), (SIZE_T)(n))
# define STRNCAT(s1,s2,n)	strncat(s1, s2, (int) (n))
# endif
#endif

#if defined(titan) | defined(sun)
#define MEMSET(p,c,n)		memset((char *)(p), (int)(c), (int)(n))
#define MEMCPY(t,f,n)		memcpy((char *)(t), (char *)(f), (int)(n))
#define STRNCAT(s1,s2,n)	strncat(s1, s2, (SIZE_T) (n))
#endif

/*
* The suffix "0" means return 0 on error.
*/
#define	FOPEN_0(fp,n,t)		{if ( !(fp = fopen(n, t)) )	FAIL1("fopen '%s'", n)}
#define	FFLUSH_0(fp)		{if (fflush(fp))		FAIL("fflush")}
#define	FSKIP_0(fp,n)		{if (fseek(fp, (long) (n), 1))	FAIL("fseek1")}
#define	FCLOSE_0(fp)		{if (fclose(fp))		FAIL("fclose")}
#define	FREAD_0(p,s,n,fp)	{int ZF = Fread(p, s, n, fp); if (ZF != (n)) \
				FAIL1("fread only found %d bytes", ZF)}
#define	FWRITE_0(p,s,n,fp)	{if ((n) != Fwrite(p, s, n, fp)) FAIL("fwrite")}

/*
* The suffix "X" means exit() on error.
*/
#define	FOPEN_X(fp,n,t)		{if ( !(fp = fopen(n, t)) )	EXIT("fopen")}
#define	FFLUSH_X(fp)		{if (fflush(fp))		EXIT("fflush")}
#define	FCLOSE_X(fp)		{if (fclose(fp))		EXIT("fclose")}
#define	FREAD_X(p,s,n,fp)	{if ((n) != Fread(p, s, n, fp))	EXIT("fread")}
#define	FWRITE_X(p,s,n,fp)	{if ((n) != Fwrite(p, s, n, fp)) EXIT("fwrite")}


/*
* The remaining macro definitions should be environment independent
*/

#define MAX_LINE_LENGTH 256	// MAX LENGTH FOR LINES OF ASCII FILES

#define SQR(x)		((x) * (x))
#define ODD(n)		((n)/2 != ((n)+1)/2)
#define EVEN(n)		((n)/2 == ((n)+1)/2)
#ifndef MAX
#	define MAX(x,y)	(((x) > (y)) ? (x) : (y))
#endif
#ifndef MIN
#	define MIN(x,y)	(((x) < (y)) ? (x) : (y))
#endif
#define ROUND(x)	floor(0.5 + (x))
#ifdef _HPUX_SOURCE
#	define RAND01		drand48()
#	define SEEDER(s)	srand48(s)
#else
#	define RAND01		(random() / 2147483647.0)
#	define SEEDER(s)	srandom(s)
#endif

#define FILE_READ	"r"
#define FILE_WRITE	"w"
#define FILE_APPEND	"a"

#ifndef TRUE
#	define	TRUE	1
#endif
#ifndef FALSE
#	define	FALSE	0
#endif

#define	SUCCESS	1
#define	FAILURE	0

#define	ALLOC_0(ptr,type,n,s)	\
	if ( !((ptr) = (type *) ALLOC(n,s)) )	FAIL2("alloc %d %d", n,s)

// uint will come from sys/types.h
typedef int sof; // SUCCESS OR FAILURE
typedef short		word;	// two byte integers
typedef unsigned short	uword;	// two byte integers, unsigned
typedef unsigned char	byte;
typedef int tristate;		// 0, 1, or 2
typedef double TYPE_CALC;	// fastest floating point type

#endif // DEFS_MACRO
