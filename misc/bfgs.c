/* bfgs.f -- translated by f2c (version 19991025).
*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "bfgs.h"

void Yanggradient (int n,  double x[],  int need[], double f0, 
		   double g[], double (*fun)( double x[]), 
		   double space[], int central,
		    double* lowbound,  double* upbound)
{

//f0=fun(x) is given for Central=0
  
  int i,j;
  double *x0=space, *x1=space+n, eh0, eh01, eh;
  eh0=eh01=1.e-8;

  if (central) {
    for (i=0;i<n;i++)  {
      if (need[i]) {
	for (j=0;j<n;j++)  x0[j]=x1[j]=x[j];
	eh=pow(eh01*(fabs(x[i])+1), 0.67);
	x0[i]-=eh; x1[i]+=eh;
	if (x0[i]<lowbound[i])
	  {x1[i]+=eh; g[i] = ((*fun)(x1)-f0)/(eh*2.0);}
	else if (x1[i]>upbound[i])
	  {x0[i]-=eh; g[i] = (f0-(*fun)(x0))/(eh*2.0);}
	else
	  g[i] = ((*fun)(x1) - (*fun)(x0))/(eh*2.0);
      }
    }
  }
  else {
    for (i=0;i<n;i++)  {
      if (need[i]) {
	for (j=0;j<n;j++)  x1[j]=x[j];
	//eh=eh0*(fabs(x[i])+1);
	eh=2.0*pow(eh0*(fabs(x[i])+1), 0.67);
	if (x1[i]+eh>upbound[i])
	  {
	    x1[i]-=eh;
	    g[i] = (f0-(*fun)(x1))/eh;
	  }
	else
	  {
	    x1[i]+=eh;
	    g[i] = ((*fun)(x1)-f0)/eh;
	  }
      }
    }
  }
}

static double oldf0;
static int CENTRALMODE=1;

void getgradient(int npar,  double invec[],  int need_gradient[], 
		 double outvec[], double(*func)( double []), 
		  double* lowbound,  double* upbound)
     
{
  double f0;
  int i;
  double *space =malloc((2*npar+10)*sizeof(double));

  CENTRALMODE=1;
  //  printf("getgradient npar=%i\n", npar);
  for (i=0; i<npar; i++)
    if (need_gradient[i]) break;
  if (i==npar) goto checkbounds;
  
  f0 = func(invec);
  
  if (CENTRALMODE > 0)
    Yanggradient(npar, invec, need_gradient, f0, outvec, func, 
		 space, 1, lowbound, upbound);
  else 
    Yanggradient(npar, invec, need_gradient, f0, outvec, func, 
		 space, 0, lowbound, upbound);
  if (CENTRALMODE == 0 && fabs(oldf0-f0) < 0.1)
    CENTRALMODE = 1;
  else if (CENTRALMODE > 0 && fabs(oldf0-f0) > 1.0)
    CENTRALMODE = 0;
  oldf0 = f0;

 checkbounds:
  for (i=0; i<npar; i++) {
    if (invec[i] <= lowbound[i] && outvec[i] > 0.0)
      outvec[i] = 0.0;
    if (invec[i] >=upbound[i] && outvec[i] < 0.0)
      outvec[i] = 0.0;
  }
  free(space);
}

static int numpar0=0;
static double (*func0)( double x[])=NULL;
static double *lowbound0=NULL, *upbound0=NULL;
static int *need0=NULL;

void numeric_likelihood( double *invec, double *outvec) {
  getgradient(numpar0, invec, need0, outvec, func0, lowbound0, upbound0);
}


typedef int logical;


//nbd is a vector of integers of dimension numpars.
//nbd[i]=0 if there are no bounds for parameter i, 
//      =1 if there are only lower bounds
//      =2 if there are both lower/upper bounds
//      =3 if there is only upper bound                
//or send nbd=NULL is equivalent to nbd=2 for all parameters
//
//noisy=0 => no output, noisy=1 => one line of output, noisy < 99 some output,
//noisy>=100 probably too much
//
//dfun is derivative function or send NULL to use numerical derivative
//(getgradient function above)
double findmax_bfgs(int numpars, double *invec, double (*fun)( double x[]),
		    void (*dfun)( double x[], double y[]),
		    double *lowbound, double *upbound, int *nbd, int noisy) {
  int i, m, isave[44], *nbd0, *iwa;
  double *grad, *wa, factr, pgtol, dsave[29], like;
  char task[60], csave[60];
  void (*dfun0)( double x[], double y[]);
  int setulb_();
  void my_strcpy(char *s1, char *s2);
  logical lsave[4];

  m=MVAL;
  factr = FACTR;
  pgtol = PGTOL;

  if (dfun==NULL) {
    dfun0=numeric_likelihood;
    func0 = fun;
    numpar0 = numpars;
    lowbound0 = lowbound;
    upbound0 = upbound;
    need0 = malloc(numpars*sizeof(int));
    for (i=0; i<numpars; i++) need0[i]=1;
  }
  else dfun0 = dfun;
  if (nbd==NULL) {
    nbd0 = malloc(numpars*sizeof(int));
    for (i=0; i<numpars; i++) nbd0[i]=2;
  }
  else nbd0 = nbd;
  
  grad = malloc(numpars*sizeof(double));
  wa = malloc(((2*m+4)*numpars + 12*m*m + 12*m)*sizeof(double));
  iwa = malloc(3*numpars*sizeof(int));
  my_strcpy(task, "START");
  for (i=5; i<60; i++) task[i]=' ';
  like = (*fun)(invec);
  (*dfun0)(invec, grad);
  while (1) {
    setulb_(&numpars, &m, invec, lowbound, upbound, nbd0, &like, 
	    grad, &factr, &pgtol, wa, iwa, task, &noisy, csave, lsave, 
	    isave, dsave);
    if (task[0]=='F' && task[1]=='G') {
      //      printf("\t");
      like = (*fun)(invec);
      (*dfun0)(invec, grad);
      if (noisy > 1) printf("like=%f\n", like);
      continue;
    }
    else if (strncmp(task, "NEW_X", 5)==0) {
      //      if (noisy > 1) printf("\n");
      continue;
    }
    else break;
  }
  if (noisy > 1) printf("\n");
  free(wa);
  free(iwa);
  free(grad);
  if (dfun==NULL) {
    lowbound0=NULL;
    upbound0=NULL;
    func0=NULL;
    free(need0);

    need0=NULL;
    numpar0 = 0;
  }
  if (nbd==NULL) free(nbd0);
  return -like;
}


void my_strcpy(char *s1, char *s2) {
  int i, len=strlen(s2);
  for (i=0; i<len; i++)
    s1[i]=s2[i];
}


typedef double doublereal;
typedef int integer;
typedef int ftnlen;


#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)
#define bit_test(a,b)	((a) >> (b) & 1)
#define bit_clear(a,b)	((a) & ~((uinteger)1 << (b)))
#define bit_set(a,b)	((a) |  ((uinteger)1 << (b)))


/* Table of ant values */

static doublereal c_b9 = 0.;
static integer c__1 = 1;
//static integer c__9 = 9;
static integer c__11 = 11;
//static integer c__3 = 3;
static doublereal c_b275 = .001;
static doublereal c_b276 = .9;
static doublereal c_b277 = .1;
//static integer c__5 = 5;

#define FALSE_ (0)
#define TRUE_ (1)

/* Subroutine */ int setulb_(int *n, int *m, double *x, double *l, double *u, 
			     int *nbd, double *f, double *g, double *factr, 
			     double *pgtol, double *wa, int *iwa, 
			     char *task, int *iprint, char *csave, int *lsave,
			     int *isave, double *dsave)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    //    integer s_cmp();

    /* Local variables */
    static integer lsnd, lsgo, lygo, l1, l2, l3, ld, lr, lt;
    extern /* Subroutine */ int mainlb_();
    static integer lz, lwa, lsg, lyg, lwn, lss, lws, lwt, lsy, lwy, lyy;

/*     ************ */

/*     Subroutine setulb */

/*     This subroutine partitions the working arrays wa and iwa, and */
/*       then uses the limited memory BFGS method to solve the bound */
/*       rained optimization problem by calling mainlb. */
/*       (The direct method will be used in the subspace minimization.) */

/*     n is an integer variable. */
/*       On entry n is the dimension of the problem. */
/*       On exit n is unchanged. */

/*     m is an integer variable. */
/*       On entry m is the maximum number of variable metric corrections */
/*         used to define the limited memory matrix. */
/*       On exit m is unchanged. */

/*     x is a double precision array of dimension n. */
/*       On entry x is an approximation to the solution. */
/*       On exit x is the current approximation. */

/*     l is a double precision array of dimension n. */
/*       On entry l is the lower bound on x. */
/*       On exit l is unchanged. */

/*     u is a double precision array of dimension n. */
/*       On entry u is the upper bound on x. */
/*       On exit u is unchanged. */

/*     nbd is an integer array of dimension n. */
/*       On entry nbd represents the type of bounds imposed on the */
/*         variables, and must be specified as follows: */
/*         nbd(i)=0 if x(i) is unbounded, */
/*                1 if x(i) has only a lower bound, */
/*                2 if x(i) has both lower and upper bounds, and */
/*                3 if x(i) has only an upper bound. */
/*       On exit nbd is unchanged. */

/*     f is a double precision variable. */
/*       On first entry f is unspecified. */
/*       On final exit f is the value of the function at x. */

/*     g is a double precision array of dimension n. */
/*       On first entry g is unspecified. */
/*       On final exit g is the value of the gradient at x. */

/*     factr is a double precision variable. */
/*       On entry factr >= 0 is specified by the user.  The iteration */
/*         will stop when */

/*         (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch */

/*         where epsmch is the machine precision, which is automatically */
/*         generated by the code. Typical values for factr: 1.d+12 for */
/*         low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely */
/*         high accuracy. */
/*       On exit factr is unchanged. */

/*     pgtol is a double precision variable. */
/*       On entry pgtol >= 0 is specified by the user.  The iteration */
/*         will stop when */

/*                 max{|proj g_i | i = 1, ..., n} <= pgtol */

/*         where pg_i is the ith component of the projected gradient. */
/*       On exit pgtol is unchanged. */

/*     wa is a double precision working array of length */
/*       (2mmax + 4)nmax + 12mmax^2 + 12mmax. */

/*     iwa is an integer working array of length 3nmax. */

/*     task is a working string of characters of length 60 indicating */
/*       the current job when entering and quitting this subroutine. */

/*     iprint is an integer variable that must be set by the user. */
/*       It controls the frequency and type of output generated: */
/*        iprint<0    no output is generated; */
/*        iprint=0    print only one line at the last iteration; */
/*        0<iprint<99 print also f and |proj g| every iprint iterations; */
/*        iprint=99   print details of every iteration except n-vectors; */
/*        iprint=100  print also the changes of active set and final x; */
/*        iprint>100  print details of every iteration including x and g; */
/*       When iprint > 0, the file iterate.dat will be created to */
/*                        summarize the iteration. */

/*     csave is a working string of characters of length 60. */

/*     lsave is a logical working array of dimension 4. */
/*       On exit with 'task' = NEW_X, the following information is */
/*                                                             available: */
/*         If lsave(1) = .true.  then  the initial X has been replaced by */
/*                                     its projection in the feasible set; */
/*         If lsave(2) = .true.  then  the problem is rained; */
/*         If lsave(3) = .true.  then  each variable has upper and lower */
/*                                     bounds; */

/*     isave is an integer working array of dimension 44. */
/*       On exit with 'task' = NEW_X, the following information is */
/*                                                             available: */
/*         isave(22) = the total number of intervals explored in the */
/*                         search of Cauchy points; */
/*         isave(26) = the total number of skipped BFGS updates before */
/*                         the current iteration; */
/*         isave(30) = the number of current iteration; */
/*         isave(31) = the total number of BFGS updates prior the current */
/*                         iteration; */
/*         isave(33) = the number of intervals explored in the search of */
/*                         Cauchy point in the current iteration; */
/*         isave(34) = the total number of function and gradient */
/*                         evaluations; */
/*         isave(36) = the number of function value or gradient */
/*                                  evaluations in the current iteration; */
/*         if isave(37) = 0  then the subspace argmin is within the box; */
/*         if isave(37) = 1  then the subspace argmin is beyond the box; */
/*         isave(38) = the number of free variables in the current */
/*                         iteration; */
/*         isave(39) = the number of active raints in the current */
/*                         iteration; */
/*         n + 1 - isave(40) = the number of variables leaving the set of */
/*                           active raints in the current iteration; */
/*         isave(41) = the number of variables entering the set of active */
/*                         raints in the current iteration. */

/*     dsave is a double precision working array of dimension 29. */
/*       On exit with 'task' = NEW_X, the following information is */
/*                                                             available: */
/*         dsave(1) = current 'theta' in the BFGS matrix; */
/*         dsave(2) = f(x) in the previous iteration; */
/*         dsave(3) = factr*epsmch; */
/*         dsave(4) = 2-norm of the line search direction vector; */
/*         dsave(5) = the machine precision epsmch generated by the code; */
/*         dsave(7) = the accumulated time spent on searching for */
/*                                                         Cauchy points; */
/*         dsave(8) = the accumulated time spent on */
/*                                                 subspace minimization; */
/*         dsave(9) = the accumulated time spent on line search; */
/*         dsave(11) = the slope of the line search function at */
/*                                  the current point of line search; */
/*         dsave(12) = the maximum relative step length imposed in */
/*                                                           line search; */
/*         dsave(13) = the infinity norm of the projected gradient; */
/*         dsave(14) = the relative step length in the line search; */
/*         dsave(15) = the slope of the line search function at */
/*                                 the starting point of the line search; */
/*         dsave(16) = the square of the 2-norm of the line search */
/*                                                      direction vector. */

/*     Subprograms called: */

/*       L-BFGS-B Library ... mainlb. */


/*     References: */

/*       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited */
/*       memory algorithm for bound rained optimization'', */
/*       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208. */

/*       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a */
/*       limited memory FORTRAN code for solving bound rained */
/*       optimization problems'', Tech. Report, NAM-11, EECS Department, */
/*       Northwestern University, 1994. */

/*       (Postscript files of these papers are available via anonymous */
/*        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.) */

/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --iwa;
    --g;
    --nbd;
    --u;
    --l;
    --x;
    --wa;
    --lsave;
    --isave;
    --dsave;

    /* Function Body */
    if (strncmp(task, "START", 5) == 0) {
	isave[1] = *m * *n;
/* Computing 2nd power */
	i__1 = *m;
	isave[2] = i__1 * i__1;
/* Computing 2nd power */
	i__1 = *m;
	isave[3] = i__1 * i__1 << 2;
	isave[4] = 1;
	isave[5] = isave[4] + isave[1];
	isave[6] = isave[5] + isave[1];
	isave[7] = isave[6] + isave[2];
	isave[8] = isave[7] + isave[2];
	isave[9] = isave[8] + isave[2];
	isave[10] = isave[9] + isave[2];
	isave[11] = isave[10] + isave[3];
	isave[12] = isave[11] + isave[3];
	isave[13] = isave[12] + *n;
	isave[14] = isave[13] + *n;
	isave[15] = isave[14] + *n;
	isave[16] = isave[15] + *n;
	isave[17] = isave[16] + (*m << 3);
	isave[18] = isave[17] + *m;
	isave[19] = isave[18] + *m;
	isave[20] = isave[19] + *m;
    }
    l1 = isave[1];
    l2 = isave[2];
    l3 = isave[3];
    lws = isave[4];
    lwy = isave[5];
    lsy = isave[6];
    lss = isave[7];
    lyy = isave[8];
    lwt = isave[9];
    lwn = isave[10];
    lsnd = isave[11];
    lz = isave[12];
    lr = isave[13];
    ld = isave[14];
    lt = isave[15];
    lwa = isave[16];
    lsg = isave[17];
    lsgo = isave[18];
    lyg = isave[19];
    lygo = isave[20];
    mainlb_(n, m, &x[1], &l[1], &u[1], &nbd[1], f, &g[1], factr, pgtol, &wa[
	    lws], &wa[lwy], &wa[lsy], &wa[lss], &wa[lyy], &wa[lwt], &wa[lwn], 
	    &wa[lsnd], &wa[lz], &wa[lr], &wa[ld], &wa[lt], &wa[lwa], &wa[lsg],
	     &wa[lsgo], &wa[lyg], &wa[lygo], &iwa[1], &iwa[*n + 1], &iwa[(*n 
	    << 1) + 1], task, iprint, csave, &lsave[1], &isave[22], &dsave[1],
	     (ftnlen)60, (ftnlen)60);
    return 0;
} /* setulb_ */

/* ======================= The end of setulb ============================= */
/* Subroutine */ int mainlb_(n, m, x, l, u, nbd, f, g, factr, pgtol, ws, wy, 
	sy, ss, yy, wt, wn, snd, z__, r__, d__, t, wa, sg, sgo, yg, ygo, 
	index, iwhere, indx2, task, iprint, csave, lsave, isave, dsave, 
	task_len, csave_len)
integer *n, *m;
doublereal *x, *l, *u;
integer *nbd;
doublereal *f, *g, *factr, *pgtol, *ws, *wy, *sy, *ss, *yy, *wt, *wn, *snd, *
	z__, *r__, *d__, *t, *wa, *sg, *sgo, *yg, *ygo;
integer *index, *iwhere, *indx2;
char *task;
integer *iprint;
char *csave;
logical *lsave;
integer *isave;
doublereal *dsave;
ftnlen task_len;
ftnlen csave_len;
{
    /* Format strings */
  /*    static char fmt_1002[] = "(/,\002At iterate\002,i5,4x,\002f= \002,1p,d12\
.5,4x,\002|proj g|= \002,1p,d12.5)";
    static char fmt_1003[] = "(2(1x,i4),5x,\002-\002,5x,\002-\002,3x,\002\
-\002,5x,\002-\002,5x,\002-\002,8x,\002-\002,3x,1p,2(1x,d10.3))";
    static char fmt_1001[] = "(//,\002ITERATION \002,i5)";
    static char fmt_1005[] = "(/,\002 Singular triangular system detected\
;\002,/,\002   refresh the lbfgs memory and restart the iteration.\002)";
    static char fmt_1006[] = "(/,\002 Nonpositive definiteness in Cholesky f\
actorization in formk;\002,/,\002   refresh the lbfgs memory and restart the\
 iteration.\002)";
    static char fmt_1008[] = "(/,\002 Bad direction in the line search;\002,\
/,\002   refresh the lbfgs memory and restart the iteration.\002)";
    static char fmt_1004[] = "(\002  ys=\002,1p,e10.3,\002  -gs=\002,1p,e10.\
3,\002 BFGS update SKIPPED\002)";
    static char fmt_1007[] = "(/,\002 Nonpositive definiteness in Cholesky f\
actorization in formt;\002,/,\002   refresh the lbfgs memory and restart the\
iteration.\002)";*/

    /* System generated locals */
    integer ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset, 
	    ss_dim1, ss_offset, yy_dim1, yy_offset, wt_dim1, wt_offset, 
	    wn_dim1, wn_offset, snd_dim1, snd_offset, i__1;
    doublereal d__1, d__2;
    int timer_();

    /* Builtin functions */
    //    integer s_cmp();
    //    /* Subroutine */ int s_copy();
    //    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static integer head;
    static doublereal fold;
    static integer nact;
    static doublereal ddum;
    extern doublereal ddot_();
    static integer info;
    static doublereal time;
    static integer nfgv, ifun, iter, nint;
    static char word[3];
    static doublereal time1, time2;
    static integer i__, iback, k;
    extern /* Subroutine */ int dscal_();
    static doublereal gdold;
    static integer nfree;
    static logical boxed;
    static integer itail;
    static doublereal theta;
    extern /* Subroutine */ int freev_(), dcopy_();
    static doublereal dnorm;
    extern /* Subroutine */ int formk_();
    static integer nskip, iword;
    extern /* Subroutine */ int formt_(), subsm_();
    static doublereal xstep, stpmx;
    extern /* Subroutine */ int prn1lb_(), prn2lb_(), prn3lb_();
    static doublereal gd, dr, rr;
    static integer ileave;
    extern /* Subroutine */ int errclb_();
    static integer itfile;
    static doublereal cachyt, epsmch;
    static logical updatd;
    static doublereal sbtime;
    extern /* Subroutine */ int active_();
    static logical prjctd;
    static integer iupdat;
    extern doublereal dpmeps_();
    static logical cnstnd;
    static doublereal sbgnrm;
    static integer nenter;
    static doublereal lnscht;
    extern /* Subroutine */ int cauchy_(), cmprlb_(), lnsrlb_(), matupd_();
    static integer nintol;
    extern /* Subroutine */ int projgr_();
    static doublereal dtd;
    static integer col;
    static doublereal tol;
    static logical wrk;
    static doublereal stp, cpu1, cpu2;

    /* Fortran I/O blocks */
    /*    static cilist io___62 = { 0, 6, 0, fmt_1002, 0 };
    static cilist io___63 = { 0, 0, 0, fmt_1003, 0 };
    static cilist io___64 = { 0, 6, 0, fmt_1001, 0 };
    static cilist io___66 = { 0, 6, 0, fmt_1005, 0 };
    static cilist io___68 = { 0, 6, 0, fmt_1006, 0 };
    static cilist io___69 = { 0, 6, 0, fmt_1005, 0 };
    static cilist io___71 = { 0, 6, 0, fmt_1008, 0 };
    static cilist io___75 = { 0, 6, 0, fmt_1004, 0 };
    static cilist io___76 = { 0, 6, 0, fmt_1007, 0 };*/

/*     ************ */

/*     Subroutine mainlb */

/*     This subroutine solves bound rained optimization problems by */
/*       using the compact formula of the limited memory BFGS updates. */

/*     n is an integer variable. */
/*       On entry n is the number of variables. */
/*       On exit n is unchanged. */

/*     m is an integer variable. */
/*       On entry m is the maximum number of variable metric */
/*          corrections allowed in the limited memory matrix. */
/*       On exit m is unchanged. */

/*     x is a double precision array of dimension n. */
/*       On entry x is an approximation to the solution. */
/*       On exit x is the current approximation. */

/*     l is a double precision array of dimension n. */
/*       On entry l is the lower bound of x. */
/*       On exit l is unchanged. */

/*     u is a double precision array of dimension n. */
/*       On entry u is the upper bound of x. */
/*       On exit u is unchanged. */

/*     nbd is an integer array of dimension n. */
/*       On entry nbd represents the type of bounds imposed on the */
/*         variables, and must be specified as follows: */
/*         nbd(i)=0 if x(i) is unbounded, */
/*                1 if x(i) has only a lower bound, */
/*                2 if x(i) has both lower and upper bounds, */
/*                3 if x(i) has only an upper bound. */
/*       On exit nbd is unchanged. */

/*     f is a double precision variable. */
/*       On first entry f is unspecified. */
/*       On final exit f is the value of the function at x. */

/*     g is a double precision array of dimension n. */
/*       On first entry g is unspecified. */
/*       On final exit g is the value of the gradient at x. */

/*     factr is a double precision variable. */
/*       On entry factr >= 0 is specified by the user.  The iteration */
/*         will stop when */

/*         (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch */

/*         where epsmch is the machine precision, which is automatically */
/*         generated by the code. */
/*       On exit factr is unchanged. */

/*     pgtol is a double precision variable. */
/*       On entry pgtol >= 0 is specified by the user.  The iteration */
/*         will stop when */

/*                 max{|proj g_i | i = 1, ..., n} <= pgtol */

/*         where pg_i is the ith component of the projected gradient. */
/*       On exit pgtol is unchanged. */

/*     ws, wy, sy, and wt are double precision working arrays used to */
/*       store the following information defining the limited memory */
/*          BFGS matrix: */
/*          ws, of dimension n x m, stores S, the matrix of s-vectors; */
/*          wy, of dimension n x m, stores Y, the matrix of y-vectors; */
/*          sy, of dimension m x m, stores S'Y; */
/*          ss, of dimension m x m, stores S'S; */
/* 	   yy, of dimension m x m, stores Y'Y; */
/*          wt, of dimension m x m, stores the Cholesky factorization */
/*                                  of (theta*S'S+LD^(-1)L'); see eq. */
/*                                  (2.26) in [3]. */

/*     wn is a double precision working array of dimension 2m x 2m */
/*       used to store the LEL^T factorization of the indefinite matrix */
/*                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ] */
/*                     [L_a -R_z           theta*S'AA'S ] */

/*       where     E = [-I  0] */
/*                     [ 0  I] */

/*     snd is a double precision working array of dimension 2m x 2m */
/*       used to store the lower triangular part of */
/*                 N = [Y' ZZ'Y   L_a'+R_z'] */
/*                     [L_a +R_z  S'AA'S   ] */

/*     z(n),r(n),d(n),t(n),wa(8*m) are double precision working arrays. */
/*       z is used at different times to store the Cauchy point and */
/*       the Newton point. */

/*     sg(m),sgo(m),yg(m),ygo(m) are double precision working arrays. */

/*     index is an integer working array of dimension n. */
/*       In subroutine freev, index is used to store the free and fixed */
/*          variables at the Generalized Cauchy Point (GCP). */

/*     iwhere is an integer working array of dimension n used to record */
/*       the status of the vector x for GCP computation. */
/*       iwhere(i)=0 or -3 if x(i) is free and has bounds, */
/*                 1       if x(i) is fixed at l(i), and l(i) .ne. u(i) */
/*                 2       if x(i) is fixed at u(i), and u(i) .ne. l(i) */
/*                 3       if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i) */
/*                -1       if x(i) is always free, i.e., no bounds on it. */

/*     indx2 is an integer working array of dimension n. */
/*       Within subroutine cauchy, indx2 corresponds to the array iorder. */
/*       In subroutine freev, a list of variables entering and leaving */
/*       the free set is stored in indx2, and it is passed on to */
/*       subroutine formk with this information. */

/*     task is a working string of characters of length 60 indicating */
/*       the current job when entering and leaving this subroutine. */

/*     iprint is an INTEGER variable that must be set by the user. */
/*       It controls the frequency and type of output generated: */
/*        iprint<0    no output is generated; */
/*        iprint=0    print only one line at the last iteration; */
/*        0<iprint<99 print also f and |proj g| every iprint iterations; */
/*        iprint=99   print details of every iteration except n-vectors; */
/*        iprint=100  print also the changes of active set and final x; */
/*        iprint>100  print details of every iteration including x and g; */
/*       When iprint > 0, the file iterate.dat will be created to */
/*                        summarize the iteration. */

/*     csave is a working string of characters of length 60. */

/*     lsave is a logical working array of dimension 4. */

/*     isave is an integer working array of dimension 23. */

/*     dsave is a double precision working array of dimension 29. */


/*     Subprograms called */

/*       L-BFGS-B Library ... cauchy, subsm, lnsrlb, formk, */

/*        errclb, prn1lb, prn2lb, prn3lb, active, projgr, */

/*        freev, cmprlb, matupd, formt. */

/*       Minpack2 Library ... timer, dpmeps. */

/*       Linpack Library ... dcopy, ddot. */


/*     References: */

/*       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited */
/*       memory algorithm for bound rained optimization'', */
/*       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208. */

/*       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN */
/*       Subroutines for Large Scale Bound rained Optimization'' */
/*       Tech. Report, NAM-11, EECS Department, Northwestern University, */
/*       1994. */

/*       [3] R. Byrd, J. Nocedal and R. Schnabel "Representations of */
/*       Quasi-Newton Matrices and their use in Limited Memory Methods'', */
/*       Mathematical Programming 63 (1994), no. 4, pp. 129-156. */

/*       (Postscript files of these papers are available via anonymous */
/*        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.) */

/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --indx2;
    --iwhere;
    --index;
    --t;
    --d__;
    --r__;
    --z__;
    --g;
    --nbd;
    --u;
    --l;
    --x;
    --ygo;
    --yg;
    --sgo;
    --sg;
    --wa;
    snd_dim1 = 2 * *m;
    snd_offset = 1 + snd_dim1 * 1;
    snd -= snd_offset;
    wn_dim1 = 2 * *m;
    wn_offset = 1 + wn_dim1 * 1;
    wn -= wn_offset;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1 * 1;
    wt -= wt_offset;
    yy_dim1 = *m;
    yy_offset = 1 + yy_dim1 * 1;
    yy -= yy_offset;
    ss_dim1 = *m;
    ss_offset = 1 + ss_dim1 * 1;
    ss -= ss_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1 * 1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1 * 1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1 * 1;
    ws -= ws_offset;
    --lsave;
    --isave;
    --dsave;

    /* Function Body */
    if (strncmp(task, "START", 5) == 0) {
      timer_(&time1);
/*        Generate the current machine precision. */
      epsmch = dpmeps_();
/*        Initialize counters and scalars when task='START'. */
/*           for the limited memory BFGS matrices: */
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = FALSE_;
/*           for operation counts: */
	iter = 0;
	nfgv = 0;
	nint = 0;
	nintol = 0;
	nskip = 0;
	nfree = *n;
/*           for stopping tolerance: */
	tol = *factr * epsmch;
/*           for measuring running time: */
	cachyt = 0.;
	sbtime = 0.;
	lnscht = 0.;
/*           'word' records the status of subspace solutions. */
	my_strcpy(word, "---");
	//	s_copy(word, "---", (ftnlen)3, (ftnlen)3);
/*           'info' records the termination information. */
	info = 0;
/*        Check the input arguments for errors. */
	errclb_(n, m, factr, &l[1], &u[1], &nbd[1], task, &info, &k, (ftnlen)
		60);
	if (strncmp(task, "ERROR", 5)==0) {
	    prn3lb_(n, &x[1], f, task, iprint, &info, &itfile, &iter, &nfgv, &
		    nintol, &nskip, &nact, &sbgnrm, &c_b9, &nint, word, &
		    iback, &stp, &xstep, &k, &cachyt, &sbtime, &lnscht, (
		    ftnlen)60, (ftnlen)3);
	    return 0;
	}
	prn1lb_(n, m, &l[1], &u[1], &x[1], iprint, &itfile, &epsmch);
/*        Initialize iwhere & project x onto the feasible set. */
	active_(n, &l[1], &u[1], &nbd[1], &x[1], &iwhere[1], iprint, &prjctd, 
		&cnstnd, &boxed);
/*        The end of the initialization. */
    } else {
/*          restore local variables. */
	prjctd = lsave[1];
	cnstnd = lsave[2];
	boxed = lsave[3];
	updatd = lsave[4];
	nintol = isave[1];
	itfile = isave[3];
	iback = isave[4];
	nskip = isave[5];
	head = isave[6];
	col = isave[7];
	itail = isave[8];
	iter = isave[9];
	iupdat = isave[10];
	nint = isave[12];
	nfgv = isave[13];
	info = isave[14];
	ifun = isave[15];
	iword = isave[16];
	nfree = isave[17];
	nact = isave[18];
	ileave = isave[19];
	nenter = isave[20];
	theta = dsave[1];
	fold = dsave[2];
	tol = dsave[3];
	dnorm = dsave[4];
	epsmch = dsave[5];
	cpu1 = dsave[6];
	cachyt = dsave[7];
	sbtime = dsave[8];
	lnscht = dsave[9];
	time1 = dsave[10];
	gd = dsave[11];
	stpmx = dsave[12];
	sbgnrm = dsave[13];
	stp = dsave[14];
	gdold = dsave[15];
	dtd = dsave[16];
/*        After returning from the driver go to the point where execution */
/*        is to resume. */
	if (strncmp(task, "FG_LN", 5)==0) {
	    goto L666;
	}
	if (strncmp(task, "NEW_X", 5)==0) {
	    goto L777;
	}
	if (strncmp(task, "FG_ST", 5)==0) {
	    goto L111;
	}
	if (strncmp(task, "STOP", 4)==0) {
	  if (strstr(task, "CPU")!=NULL) {
/*                                          restore the previous iterate. */
		dcopy_(n, &t[1], &c__1, &x[1], &c__1);
		dcopy_(n, &r__[1], &c__1, &g[1], &c__1);
		*f = fold;
	    }
	    goto L999;
	}
    }
/*     Compute f0 and g0. */
    my_strcpy(task, "FG_START");
/*          return to the driver to calculate f and g; reenter at 111. */
    goto L1000;
L111:
    nfgv = 1;
/*     Compute the infinity norm of the (-) projected gradient. */
    projgr_(n, &l[1], &u[1], &nbd[1], &x[1], &g[1], &sbgnrm);
    if (*iprint >= 1) {
      printf("At iterate %4i f=%f\t|proj g| = %e\n", iter, *f, sbgnrm);
      /*	s_wsfe(&io___62);
	do_fio(&c__1, (char *)&iter, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*f), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&sbgnrm, (ftnlen)sizeof(doublereal));
	e_wsfe();
	io___63.ciunit = itfile;
      
	s_wsfe(&io___63);
	do_fio(&c__1, (char *)&iter, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nfgv, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&sbgnrm, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*f), (ftnlen)sizeof(doublereal));
	e_wsfe();*/
    }
    if (sbgnrm <= *pgtol) {
/*                                terminate the algorithm. */
      my_strcpy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL");
      //	s_copy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL", (
      //		ftnlen)60, (ftnlen)48);
	goto L999;
    }
/* ----------------- the beginning of the loop -------------------------- */
L222:
    if (*iprint >= 99) {
      i__1 = iter + 1;
      printf("ITERATION %i\n", i__1);
      /*	s_wsfe(&io___64);
	i__1 = iter + 1;
	do_fio(&c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	e_wsfe();*/
    }
    iword = -1;

    if (! cnstnd && col > 0) {
/*                                            skip the search for GCP. */
	dcopy_(n, &x[1], &c__1, &z__[1], &c__1);
	wrk = updatd;
	nint = 0;
	goto L333;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*     Compute the Generalized Cauchy Point (GCP). */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    timer_(&cpu1);
    cauchy_(n, &x[1], &l[1], &u[1], &nbd[1], &g[1], &indx2[1], &iwhere[1], &t[
	    1], &d__[1], &z__[1], m, &wy[wy_offset], &ws[ws_offset], &sy[
	    sy_offset], &wt[wt_offset], &theta, &col, &head, &wa[1], &wa[(*m 
	    << 1) + 1], &wa[(*m << 2) + 1], &wa[*m * 6 + 1], &nint, &sg[1], &
	    yg[1], iprint, &sbgnrm, &info);
    if (info != 0) {
/*         singular triangular system detected; refresh the lbfgs memory. */
	if (*iprint >= 1) {
	  printf("Singular triangular system detected; refresh the lbfgs memory and restart the iteration.\n");
	  //	    s_wsfe(&io___66);
	  //	    e_wsfe();
	}
	info = 0;
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = FALSE_;
	timer_(&cpu2);
	cachyt = cachyt + cpu2 - cpu1;
	goto L222;
    }
    timer_(&cpu2);
    cachyt = cachyt + cpu2 - cpu1;
    nintol += nint;
/*     Count the entering and leaving variables for iter > 0; */
/*     find the index set of free and active variables at the GCP. */
    freev_(n, &nfree, &index[1], &nenter, &ileave, &indx2[1], &iwhere[1], &
	    wrk, &updatd, &cnstnd, iprint, &iter);
    nact = *n - nfree;
L333:
/*     If there are no free variables or B=I, then */
/*                                        skip the subspace minimization. */
    if (nfree == 0 || col == 0) {
	goto L555;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*     Subspace minimization. */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    timer_(&cpu1);
/*     Form  the LEL^T factorization of the indefinite */
/*       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ] */
/*                     [L_a -R_z           theta*S'AA'S ] */
/*       where     E = [-I  0] */
/*                     [ 0  I] */
    if (wrk) {
	formk_(n, &nfree, &index[1], &nenter, &ileave, &indx2[1], &iupdat, &
		updatd, &wn[wn_offset], &snd[snd_offset], m, &ws[ws_offset], &
		wy[wy_offset], &sy[sy_offset], &theta, &col, &head, &info);
    }
    if (info != 0) {
/*          nonpositive definiteness in Cholesky factorization; */
/*          refresh the lbfgs memory and restart the iteration. */
	if (*iprint >= 1) {
	  printf("Nonpositive defiteness in Cholesky factorization; refresh the lbfgs memory ane restart the iteration\n");
	  //	    s_wsfe(&io___68);
	  //	    e_wsfe();
	}
	info = 0;
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = FALSE_;
	timer_(&cpu2);
	sbtime = sbtime + cpu2 - cpu1;
	goto L222;
    }
/*        compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x) */
/*                                                   from 'cauchy'). */
    cmprlb_(n, m, &x[1], &g[1], &ws[ws_offset], &wy[wy_offset], &sy[sy_offset]
	    , &wt[wt_offset], &z__[1], &r__[1], &wa[1], &index[1], &theta, &
	    col, &head, &nfree, &cnstnd, &info);
    if (info != 0) {
	goto L444;
    }
/*       call the direct method. */
    subsm_(n, m, &nfree, &index[1], &l[1], &u[1], &nbd[1], &z__[1], &r__[1], &
	    ws[ws_offset], &wy[wy_offset], &theta, &col, &head, &iword, &wa[1]
	    , &wn[wn_offset], iprint, &info);
L444:
    if (info != 0) {
/*          singular triangular system detected; */
/*          refresh the lbfgs memory and restart the iteration. */
	if (*iprint >= 1) {
	  printf("Singular triangular system detected; refresh the lbfgs memory and restart the iteration\n");
	  //	    s_wsfe(&io___69);
	  //	    e_wsfe();
	}
	info = 0;
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = FALSE_;
	timer_(&cpu2);
	sbtime = sbtime + cpu2 - cpu1;
	goto L222;
    }
    timer_(&cpu2);
    sbtime = sbtime + cpu2 - cpu1;
L555:
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*     Line search and optimality tests. */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     Generate the search direction d:=z-x. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__] = z__[i__] - x[i__];
/* L40: */
    }
    timer_(&cpu1);
L666:
    lnsrlb_(n, &l[1], &u[1], &nbd[1], &x[1], f, &fold, &gd, &gdold, &g[1], &
	    d__[1], &r__[1], &t[1], &z__[1], &stp, &dnorm, &dtd, &xstep, &
	    stpmx, &iter, &ifun, &iback, &nfgv, &info, task, &boxed, &cnstnd, 
	    csave, &isave[22], &dsave[17], (ftnlen)60, (ftnlen)60);
    if (info != 0 || iback >= 20) {
/*          restore the previous iterate. */
	dcopy_(n, &t[1], &c__1, &x[1], &c__1);
	dcopy_(n, &r__[1], &c__1, &g[1], &c__1);
	*f = fold;
	if (col == 0) {
/*             abnormal termination. */
	    if (info == 0) {
		info = -9;
/*                restore the actual number of f and g evaluations etc. */
		--nfgv;
		--ifun;
		--iback;
	    }
	    my_strcpy(task, "ABNORMAL_TERMINATION_IN_LNSRCH");
	    //	    s_copy(task, "ABNORMAL_TERMINATION_IN_LNSRCH", (ftnlen)60, (
	    //		    ftnlen)30);
	    ++iter;
	    goto L999;
	} else {
/*             refresh the lbfgs memory and restart the iteration. */
	  if (*iprint >= 1) {
	    printf("refresh the lbfgs memory and restart the iteration\n");
	    //		s_wsfe(&io___71);
	    //		e_wsfe();
	    }
	    if (info == 0) {
		--nfgv;
	    }
	    info = 0;
	    col = 0;
	    head = 1;
	    theta = 1.;
	    iupdat = 0;
	    updatd = FALSE_;
	    my_strcpy(task, "RESTART_FROM_LNSRCH");
	    //	    s_copy(task, "RESTART_FROM_LNSRCH", (ftnlen)60, (ftnlen)19);
	    timer_(&cpu2);
	    lnscht = lnscht + cpu2 - cpu1;
	    goto L222;
	}
    } else if (strncmp(task, "FG_LN", 5)==0) {
/*          return to the driver for calculating f and g; reenter at 666. */
	goto L1000;
    } else {
/*          calculate and print out the quantities related to the new X. */
	timer_(&cpu2);
	lnscht = lnscht + cpu2 - cpu1;
	++iter;
/*        Compute the infinity norm of the projected (-)gradient. */
	projgr_(n, &l[1], &u[1], &nbd[1], &x[1], &g[1], &sbgnrm);
/*        Print iteration information. */
	prn2lb_(n, &x[1], f, &g[1], iprint, &itfile, &iter, &nfgv, &nact, &
		sbgnrm, &nint, word, &iword, &iback, &stp, &xstep, (ftnlen)3);
	goto L1000;
    }
L777:
/*     Test for termination. */
    if (sbgnrm <= *pgtol) {
/*                                terminate the algorithm. */
      my_strcpy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL");
	goto L999;
    }
/* Computing MAX */
    d__1 = abs(fold), d__2 = abs(*f), d__1 = max(d__1,d__2);
    ddum = max(d__1,1.);
    if (fold - *f <= tol * ddum) {
/*                                        terminate the algorithm. */
      my_strcpy(task, "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH");
	if (iback >= 10) {
	    info = -5;
	}
/*           i.e., to issue a warning if iback>10 in the line search. */
	goto L999;
    }
/*     Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__[i__] = g[i__] - r__[i__];
/* L42: */
    }
    rr = ddot_(n, &r__[1], &c__1, &r__[1], &c__1);
    if (stp == 1.) {
	dr = gd - gdold;
	ddum = -gdold;
    } else {
	dr = (gd - gdold) * stp;
	dscal_(n, &stp, &d__[1], &c__1);
	ddum = -gdold * stp;
    }
    if (dr <= epsmch * ddum) {
/*                            skip the L-BFGS update. */
	++nskip;
	updatd = FALSE_;
	if (*iprint >= 1) {
	  printf("ys=%f, gs=%f, BFGS update SKIPPED\n", dr, ddum);
	  /*s_wsfe(&io___75);
	    do_fio(&c__1, (char *)&dr, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&ddum, (ftnlen)sizeof(doublereal));
	    e_wsfe();*/
	}
	goto L888;
    }
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */

/*     Update the L-BFGS matrix. */

/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
    updatd = TRUE_;
    ++iupdat;
/*     Update matrices WS and WY and form the middle matrix in B. */
    matupd_(n, m, &ws[ws_offset], &wy[wy_offset], &sy[sy_offset], &ss[
	    ss_offset], &d__[1], &r__[1], &itail, &iupdat, &col, &head, &
	    theta, &rr, &dr, &stp, &dtd);
/*     Form the upper half of the pds T = theta*SS + L*D^(-1)*L'; */
/*        Store T in the upper triangular of the array wt; */
/*        Cholesky factorize T to J*J' with */
/*           J' stored in the upper triangular of wt. */
    formt_(m, &wt[wt_offset], &sy[sy_offset], &ss[ss_offset], &col, &theta, &
	    info);
    if (info != 0) {
/*          nonpositive definiteness in Cholesky factorization; */
/*          refresh the lbfgs memory and restart the iteration. */
	if (*iprint >= 1) {
	  printf("Nonpositive definiteness in Cholesky factorization; refresh the lbfgs memory and restart the iteration.\n");
	  //	    s_wsfe(&io___76);
	  //	    e_wsfe();
	}
	info = 0;
	col = 0;
	head = 1;
	theta = 1.;
	iupdat = 0;
	updatd = FALSE_;
	goto L222;
    }
/*     Now the inverse of the middle matrix in B is */
/*       [  D^(1/2)      O ] [ -D^(1/2)  D^(-1/2)*L' ] */
/*       [ -L*D^(-1/2)   J ] [  0        J'          ] */
L888:
/* -------------------- the end of the loop ----------------------------- */
    goto L222;
L999:
    timer_(&time2);
    time = time2 - time1;
    prn3lb_(n, &x[1], f, task, iprint, &info, &itfile, &iter, &nfgv, &nintol, 
	    &nskip, &nact, &sbgnrm, &time, &nint, word, &iback, &stp, &xstep, 
	    &k, &cachyt, &sbtime, &lnscht, (ftnlen)60, (ftnlen)3);
L1000:
/*     Save local variables. */
    lsave[1] = prjctd;
    lsave[2] = cnstnd;
    lsave[3] = boxed;
    lsave[4] = updatd;
    isave[1] = nintol;
    isave[3] = itfile;
    isave[4] = iback;
    isave[5] = nskip;
    isave[6] = head;
    isave[7] = col;
    isave[8] = itail;
    isave[9] = iter;
    isave[10] = iupdat;
    isave[12] = nint;
    isave[13] = nfgv;
    isave[14] = info;
    isave[15] = ifun;
    isave[16] = iword;
    isave[17] = nfree;
    isave[18] = nact;
    isave[19] = ileave;
    isave[20] = nenter;
    dsave[1] = theta;
    dsave[2] = fold;
    dsave[3] = tol;
    dsave[4] = dnorm;
    dsave[5] = epsmch;
    dsave[6] = cpu1;
    dsave[7] = cachyt;
    dsave[8] = sbtime;
    dsave[9] = lnscht;
    dsave[10] = time1;
    dsave[11] = gd;
    dsave[12] = stpmx;
    dsave[13] = sbgnrm;
    dsave[14] = stp;
    dsave[15] = gdold;
    dsave[16] = dtd;
    return 0;
} /* mainlb_ */

/* ======================= The end of mainlb ============================= */
/* Subroutine */ int active_(n, l, u, nbd, x, iwhere, iprint, prjctd, cnstnd, 
	boxed)
integer *n;
doublereal *l, *u;
integer *nbd;
doublereal *x;
integer *iwhere, *iprint;
logical *prjctd, *cnstnd, *boxed;
{
    /* Format strings */
  //    static char fmt_1001[] = "(/,\002At X0 \002,i9,\002 variables are exactly at the bounds\002)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle(), s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static integer nbdd, i__;

    /* Fortran I/O blocks */
    /*    static cilist io___81 = { 0, 6, 0, 0, 0 };
    static cilist io___82 = { 0, 6, 0, 0, 0 };
    static cilist io___83 = { 0, 6, 0, fmt_1001, 0 };*/


/*     ************ */

/*     Subroutine active */

/*     This subroutine initializes iwhere and projects the initial x to */
/*       the feasible set if necessary. */

/*     iwhere is an integer array of dimension n. */
/*       On entry iwhere is unspecified. */
/*       On exit iwhere(i)=-1  if x(i) has no bounds */
/*                         3   if l(i)=u(i) */
/*                         0   otherwise. */
/*       In cauchy, iwhere is given finer gradations. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
/*     Initialize nbdd, prjctd, cnstnd and boxed. */
    /* Parameter adjustments */
    --iwhere;
    --x;
    --nbd;
    --u;
    --l;

    /* Function Body */
    nbdd = 0;
    *prjctd = FALSE_;
    *cnstnd = FALSE_;
    *boxed = TRUE_;
/*     Project the initial x to the easible set if necessary. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (nbd[i__] > 0) {
	    if (nbd[i__] <= 2 && x[i__] <= l[i__]) {
		if (x[i__] < l[i__]) {
		    *prjctd = TRUE_;
		    x[i__] = l[i__];
		}
		++nbdd;
	    } else if (nbd[i__] >= 2 && x[i__] >= u[i__]) {
		if (x[i__] > u[i__]) {
		    *prjctd = TRUE_;
		    x[i__] = u[i__];
		}
		++nbdd;
	    }
	}
/* L10: */
    }
/*     Initialize iwhere and assign values to cnstnd and boxed. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (nbd[i__] != 2) {
	    *boxed = FALSE_;
	}
	if (nbd[i__] == 0) {
/*                                this variable is always free */
	    iwhere[i__] = -1;
/*           otherwise set x(i)=mid(x(i), u(i), l(i)). */
	} else {
	    *cnstnd = TRUE_;
	    if (nbd[i__] == 2 && u[i__] - l[i__] <= 0.) {
/*                   this variable is always fixed */
		iwhere[i__] = 3;
	    } else {
		iwhere[i__] = 0;
	    }
	}
/* L20: */
    }
    if (*iprint >= 0) {
	if (*prjctd) {
	  printf("The initial X is infeasible. Restart with its projection\n");
	  /*	    s_wsle(&io___81);
	    do_lio(&c__9, &c__1, "The initial X is infeasible.  Restart with\
 its projection.", (ftnlen)58);
 e_wsle();*/
	}
	if (! (*cnstnd)) {
	  printf("The problem is unrained\n");
	  /*	    s_wsle(&io___82);
	    do_lio(&c__9, &c__1, "This problem is unrained.", (ftnlen)30)
		    ;
		    e_wsle();*/
	}
    }
    if (*iprint > 0) {
      printf("ITERATION %i\n", nbdd);
      /*	s_wsfe(&io___83);
	do_fio(&c__1, (char *)&nbdd, (ftnlen)sizeof(integer));
	e_wsfe();*/
    }
    return 0;
} /* active_ */

/* ======================= The end of active ============================= */
/* Subroutine */ int bmv_(m, sy, wt, col, v, p, info)
integer *m;
doublereal *sy, *wt;
integer *col;
doublereal *v, *p;
integer *info;
{
    /* System generated locals */
    integer sy_dim1, sy_offset, wt_dim1, wt_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer i__, k;
    extern /* Subroutine */ int dtrsl_();
    static integer i2;
    static doublereal sum;

/*     ************ */

/*     Subroutine bmv */

/*     This subroutine computes the product of the 2m x 2m middle matrix */
/* 	in the compact L-BFGS formula of B and a 2m vector v; */
/* 	it returns the product in p. */

/*     m is an integer variable. */
/*       On entry m is the maximum number of variable metric corrections */
/*         used to define the limited memory matrix. */
/*       On exit m is unchanged. */

/*     sy is a double precision array of dimension m x m. */
/*       On entry sy specifies the matrix S'Y. */
/*       On exit sy is unchanged. */

/*     wt is a double precision array of dimension m x m. */
/*       On entry wt specifies the upper triangular matrix J' which is */
/*         the Cholesky factor of (thetaS'S+LD^(-1)L'). */
/*       On exit wt is unchanged. */

/*     col is an integer variable. */
/*       On entry col specifies the number of s-vectors (or y-vectors) */
/*         stored in the compact L-BFGS formula. */
/*       On exit col is unchanged. */

/*     v is a double precision array of dimension 2col. */
/*       On entry v specifies vector v. */
/*       On exit v is unchanged. */

/*     p is a double precision array of dimension 2col. */
/*       On entry p is unspecified. */
/*       On exit p is the product Mv. */

/*     info is an integer variable. */
/*       On entry info is unspecified. */
/*       On exit info = 0       for normal return, */
/*                    = nonzero for abnormal return when the system */
/*                                to be solved by dtrsl is singular. */

/*     Subprograms called: */

/*       Linpack ... dtrsl. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1 * 1;
    wt -= wt_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1 * 1;
    sy -= sy_offset;
    --p;
    --v;

    /* Function Body */
    if (*col == 0) {
	return 0;
    }
/*     PART I: solve [  D^(1/2)      O ] [ p1 ] = [ v1 ] */
/*                   [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ]. */
/* 	solve Jp2=v2+LD^(-1)v1. */
    p[*col + 1] = v[*col + 1];
    i__1 = *col;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i2 = *col + i__;
	sum = 0.;
	i__2 = i__ - 1;
	for (k = 1; k <= i__2; ++k) {
	    sum += sy[i__ + k * sy_dim1] * v[k] / sy[k + k * sy_dim1];
/* L10: */
	}
	p[i2] = v[i2] + sum;
/* L20: */
    }
/*     Solve the triangular system */
    dtrsl_(&wt[wt_offset], m, col, &p[*col + 1], &c__11, info);
    if (*info != 0) {
	return 0;
    }
/*     	solve D^(1/2)p1=v1. */
    i__1 = *col;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = v[i__] / sqrt(sy[i__ + i__ * sy_dim1]);
/* L30: */
    }
/*     PART II: solve [ -D^(1/2)   D^(-1/2)*L'  ] [ p1 ] = [ p1 ] */
/*                    [  0         J'           ] [ p2 ]   [ p2 ]. */
/*       solve J^Tp2=p2. */
    dtrsl_(&wt[wt_offset], m, col, &p[*col + 1], &c__1, info);
    if (*info != 0) {
	return 0;
    }
/*       compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2) */
/*                 =-D^(-1/2)p1+D^(-1)L'p2. */
    i__1 = *col;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = -p[i__] / sqrt(sy[i__ + i__ * sy_dim1]);
/* L40: */
    }
    i__1 = *col;
    for (i__ = 1; i__ <= i__1; ++i__) {
	sum = 0.;
	i__2 = *col;
	for (k = i__ + 1; k <= i__2; ++k) {
	    sum += sy[k + i__ * sy_dim1] * p[*col + k] / sy[i__ + i__ * 
		    sy_dim1];
/* L50: */
	}
	p[i__] += sum;
/* L60: */
    }
    return 0;
} /* bmv_ */

/* ======================== The end of bmv =============================== */
/* Subroutine */ int cauchy_(n, x, l, u, nbd, g, iorder, iwhere, t, d__, xcp, 
	m, wy, ws, sy, wt, theta, col, head, p, c__, wbp, v, nint, sg, yg, 
	iprint, sbgnrm, info)
integer *n;
doublereal *x, *l, *u;
integer *nbd;
doublereal *g;
integer *iorder, *iwhere;
doublereal *t, *d__, *xcp;
integer *m;
doublereal *wy, *ws, *sy, *wt, *theta;
integer *col, *head;
doublereal *p, *c__, *wbp, *v;
integer *nint;
doublereal *sg, *yg;
integer *iprint;
doublereal *sbgnrm;
integer *info;
{
    /* Format strings */
  /*    static char fmt_3010[] = "(/,\002---------------- CAUCHY entered--------\
-----------\002)";
    static char fmt_1010[] = "(\002Cauchy X =  \002,/,(4x,1p,6(1x,d11.4)))";
    static char fmt_4011[] = "(/,\002Piece    \002,i3,\002 --f1, f2 at start\
 point \002,1p,2(1x,d11.4))";
    static char fmt_5010[] = "(\002Distance to the next break point =  \002,\
1p,d11.4)";
    static char fmt_6010[] = "(\002Distance to the stationary point =  \002,\
1p,d11.4)";
    static char fmt_4010[] = "(\002Piece    \002,i3,\002 --f1, f2 at start p\
oint \002,1p,2(1x,d11.4))";
    static char fmt_2010[] = "(/,\002---------------- exit CAUCHY-----------\
    -----------\002,/)";*/

    /* System generated locals */
    integer wy_dim1, wy_offset, ws_dim1, ws_offset, sy_dim1, sy_offset, 
	    wt_dim1, wt_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle(), s_wsfe(), e_wsfe(), do_fio();

    /* Local variables */
    static doublereal dibp;
    extern doublereal ddot_();
    static integer iter;
    static doublereal zibp, tsum, dibp2;
    static integer i__, j;
    static logical bnded;
    extern /* Subroutine */ int dscal_();
    static doublereal neggi;
    static integer nfree;
    static doublereal bkmin;
    static integer nleft;
    extern /* Subroutine */ int dcopy_(), daxpy_();
    static doublereal f1, f2, dt, tj, tl;
    static integer nbreak, ibkmin;
    static doublereal tu;
    extern /* Subroutine */ int hpsolb_();
    static integer pointr;
    static doublereal tj0;
    static logical xlower, xupper;
    static integer ibp;
    static doublereal dtm;
    extern /* Subroutine */ int bmv_();
    static doublereal wmc, wmp, wmw;
    static integer col2;

    /* Fortran I/O blocks */
    /*    static cilist io___88 = { 0, 6, 0, 0, 0 };
    static cilist io___96 = { 0, 6, 0, fmt_3010, 0 };
    static cilist io___105 = { 0, 6, 0, fmt_1010, 0 };
    static cilist io___109 = { 0, 6, 0, 0, 0 };
    static cilist io___116 = { 0, 6, 0, fmt_4011, 0 };
    static cilist io___117 = { 0, 6, 0, fmt_5010, 0 };
    static cilist io___118 = { 0, 6, 0, fmt_6010, 0 };
    static cilist io___121 = { 0, 6, 0, 0, 0 };
    static cilist io___126 = { 0, 6, 0, 0, 0 };
    static cilist io___127 = { 0, 6, 0, 0, 0 };
    static cilist io___128 = { 0, 6, 0, fmt_4010, 0 };
    static cilist io___129 = { 0, 6, 0, fmt_6010, 0 };
    static cilist io___130 = { 0, 6, 0, fmt_1010, 0 };
    static cilist io___131 = { 0, 6, 0, fmt_2010, 0 };*/


/*     ************ */

/*     Subroutine cauchy */

/*     For given x, l, u, g (with sbgnrm > 0), and a limited memory */
/*       BFGS matrix B defined in terms of matrices WY, WS, WT, and */
/*       scalars head, col, and theta, this subroutine computes the */
/*       generalized Cauchy point (GCP), defined as the first local */
/*       minimizer of the quadratic */

/*                  Q(x + s) = g's + 1/2 s'Bs */

/*       along the projected gradient direction P(x-tg,l,u). */
/*       The routine returns the GCP in xcp. */

/*     n is an integer variable. */
/*       On entry n is the dimension of the problem. */
/*       On exit n is unchanged. */

/*     x is a double precision array of dimension n. */
/*       On entry x is the starting point for the GCP computation. */
/*       On exit x is unchanged. */

/*     l is a double precision array of dimension n. */
/*       On entry l is the lower bound of x. */
/*       On exit l is unchanged. */

/*     u is a double precision array of dimension n. */
/*       On entry u is the upper bound of x. */
/*       On exit u is unchanged. */

/*     nbd is an integer array of dimension n. */
/*       On entry nbd represents the type of bounds imposed on the */
/*         variables, and must be specified as follows: */
/*         nbd(i)=0 if x(i) is unbounded, */
/*                1 if x(i) has only a lower bound, */
/*                2 if x(i) has both lower and upper bounds, and */
/*                3 if x(i) has only an upper bound. */
/*       On exit nbd is unchanged. */

/*     g is a double precision array of dimension n. */
/*       On entry g is the gradient of f(x).  g must be a nonzero vector. */
/*       On exit g is unchanged. */

/*     iorder is an integer working array of dimension n. */
/*       iorder will be used to store the breakpoints in the piecewise */
/*       linear path and free variables encountered. On exit, */
/*         iorder(1),...,iorder(nleft) are indices of breakpoints */
/*                                which have not been encountered; */
/*         iorder(nleft+1),...,iorder(nbreak) are indices of */
/*                                     encountered breakpoints; and */
/*         iorder(nfree),...,iorder(n) are indices of variables which */
/*                 have no bound raits along the search direction. */

/*     iwhere is an integer array of dimension n. */
/*       On entry iwhere indicates only the permanently fixed (iwhere=3) */
/*       or free (iwhere= -1) components of x. */
/*       On exit iwhere records the status of the current x variables. */
/*       iwhere(i)=-3  if x(i) is free and has bounds, but is not moved */
/*                 0   if x(i) is free and has bounds, and is moved */
/*                 1   if x(i) is fixed at l(i), and l(i) .ne. u(i) */
/*                 2   if x(i) is fixed at u(i), and u(i) .ne. l(i) */
/*                 3   if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i) */
/*                 -1  if x(i) is always free, i.e., it has no bounds. */

/*     t is a double precision working array of dimension n. */
/*       t will be used to store the break points. */

/*     d is a double precision array of dimension n used to store */
/*       the Cauchy direction P(x-tg)-x. */

/*     xcp is a double precision array of dimension n used to return the */
/*       GCP on exit. */

/*     m is an integer variable. */
/*       On entry m is the maximum number of variable metric corrections */
/*         used to define the limited memory matrix. */
/*       On exit m is unchanged. */

/*     ws, wy, sy, and wt are double precision arrays. */
/*       On entry they store information that defines the */
/*                             limited memory BFGS matrix: */
/*         ws(n,m) stores S, a set of s-vectors; */
/*         wy(n,m) stores Y, a set of y-vectors; */
/*         sy(m,m) stores S'Y; */
/*         wt(m,m) stores the */
/*                 Cholesky factorization of (theta*S'S+LD^(-1)L'). */
/*       On exit these arrays are unchanged. */

/*     theta is a double precision variable. */
/*       On entry theta is the scaling factor specifying B_0 = theta I. */
/*       On exit theta is unchanged. */

/*     col is an integer variable. */
/*       On entry col is the actual number of variable metric */
/*         corrections stored so far. */
/*       On exit col is unchanged. */

/*     head is an integer variable. */
/*       On entry head is the location of the first s-vector (or y-vector) */
/*         in S (or Y). */
/*       On exit col is unchanged. */

/*     p is a double precision working array of dimension 2m. */
/*       p will be used to store the vector p = W^(T)d. */

/*     c is a double precision working array of dimension 2m. */
/*       c will be used to store the vector c = W^(T)(xcp-x). */

/*     wbp is a double precision working array of dimension 2m. */
/*       wbp will be used to store the row of W corresponding */
/*         to a breakpoint. */

/*     v is a double precision working array of dimension 2m. */

/*     nint is an integer variable. */
/*       On exit nint records the number of quadratic segments explored */
/*         in searching for the GCP. */

/*     sg and yg are double precision arrays of dimension m. */
/*       On entry sg  and yg store S'g and Y'g correspondingly. */
/*       On exit they are unchanged. */

/*     iprint is an INTEGER variable that must be set by the user. */
/*       It controls the frequency and type of output generated: */
/*        iprint<0    no output is generated; */
/*        iprint=0    print only one line at the last iteration; */
/*        0<iprint<99 print also f and |proj g| every iprint iterations; */
/*        iprint=99   print details of every iteration except n-vectors; */
/*        iprint=100  print also the changes of active set and final x; */
/*        iprint>100  print details of every iteration including x and g; */
/*       When iprint > 0, the file iterate.dat will be created to */
/*                        summarize the iteration. */

/*     sbgnrm is a double precision variable. */
/*       On entry sbgnrm is the norm of the projected gradient at x. */
/*       On exit sbgnrm is unchanged. */

/*     info is an integer variable. */
/*       On entry info is 0. */
/*       On exit info = 0       for normal return, */
/*                    = nonzero for abnormal return when the the system */
/*                              used in routine bmv is singular. */

/*     Subprograms called: */

/*       L-BFGS-B Library ... hpsolb, bmv. */

/*       Linpack ... dscal dcopy, daxpy. */


/*     References: */

/*       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited */
/*       memory algorithm for bound rained optimization'', */
/*       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208. */

/*       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN */
/*       Subroutines for Large Scale Bound rained Optimization'' */
/*       Tech. Report, NAM-11, EECS Department, Northwestern University, */
/*       1994. */

/*       (Postscript files of these papers are available via anonymous */
/*        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.) */

/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
/*     Check the status of the variables, reset iwhere(i) if necessary; */
/*       compute the Cauchy direction d and the breakpoints t; initialize */
/*       the derivative f1 and the vector p = W'd (for theta = 1). */
    /* Parameter adjustments */
    --xcp;
    --d__;
    --t;
    --iwhere;
    --iorder;
    --g;
    --nbd;
    --u;
    --l;
    --x;
    --yg;
    --sg;
    --v;
    --wbp;
    --c__;
    --p;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1 * 1;
    wt -= wt_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1 * 1;
    sy -= sy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1 * 1;
    ws -= ws_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1 * 1;
    wy -= wy_offset;

    /* Function Body */
    if (*sbgnrm <= 0.) {
	if (*iprint >= 0) {
	  printf("Subgnorm = 0. GCP = X.\n");
	  /*	    s_wsle(&io___88);
	    do_lio(&c__9, &c__1, "Subgnorm = 0.  GCP = X.", (ftnlen)23);
	    e_wsle();*/
	}
	dcopy_(n, &x[1], &c__1, &xcp[1], &c__1);
	return 0;
    }
    bnded = TRUE_;
    nfree = *n + 1;
    nbreak = 0;
    ibkmin = 0;
    bkmin = 0.;
    col2 = *col << 1;
    f1 = 0.;
    if (*iprint >= 99) {
      printf("---------------- CAUCHY entered-------------------\n");
      /*	s_wsfe(&io___96);
		e_wsfe();*/
    }
/*     We set p to zero and build it up as we determine d. */
    i__1 = col2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	p[i__] = 0.;
/* L20: */
    }
/*     In the following loop we determine for each variable its bound */
/*        status and its breakpoint, and update p accordingly. */
/*        Smallest breakpoint is identified. */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	neggi = -g[i__];
	if (iwhere[i__] != 3 && iwhere[i__] != -1) {
/*             if x(i) is not a ant and has bounds, */
/*             compute the difference between x(i) and its bounds. */
	    if (nbd[i__] <= 2) {
		tl = x[i__] - l[i__];
	    }
	    if (nbd[i__] >= 2) {
		tu = u[i__] - x[i__];
	    }
/*           If a variable is close enough to a bound */
/*             we treat it as at bound. */
	    xlower = nbd[i__] <= 2 && tl <= 0.;
	    xupper = nbd[i__] >= 2 && tu <= 0.;
/*              reset iwhere(i). */
	    iwhere[i__] = 0;
	    if (xlower) {
		if (neggi <= 0.) {
		    iwhere[i__] = 1;
		}
	    } else if (xupper) {
		if (neggi >= 0.) {
		    iwhere[i__] = 2;
		}
	    } else {
		if (abs(neggi) <= 0.) {
		    iwhere[i__] = -3;
		}
	    }
	}
	pointr = *head;
	if (iwhere[i__] != 0 && iwhere[i__] != -1) {
	    d__[i__] = 0.;
	} else {
	    d__[i__] = neggi;
	    f1 -= neggi * neggi;
/*             calculate p := p - W'e_i* (g_i). */
	    i__2 = *col;
	    for (j = 1; j <= i__2; ++j) {
		p[j] += wy[i__ + pointr * wy_dim1] * neggi;
		p[*col + j] += ws[i__ + pointr * ws_dim1] * neggi;
		pointr = pointr % *m + 1;
/* L40: */
	    }
	    if (nbd[i__] <= 2 && nbd[i__] != 0 && neggi < 0.) {
/*                                 x(i) + d(i) is bounded; compute t(i). */
		++nbreak;
		iorder[nbreak] = i__;
		t[nbreak] = tl / (-neggi);
		if (nbreak == 1 || t[nbreak] < bkmin) {
		    bkmin = t[nbreak];
		    ibkmin = nbreak;
		}
	    } else if (nbd[i__] >= 2 && neggi > 0.) {
/*                                 x(i) + d(i) is bounded; compute t(i). */
		++nbreak;
		iorder[nbreak] = i__;
		t[nbreak] = tu / neggi;
		if (nbreak == 1 || t[nbreak] < bkmin) {
		    bkmin = t[nbreak];
		    ibkmin = nbreak;
		}
	    } else {
/*                x(i) + d(i) is not bounded. */
		--nfree;
		iorder[nfree] = i__;
		if (abs(neggi) > 0.) {
		    bnded = FALSE_;
		}
	    }
	}
/* L50: */
    }
/*     The indices of the nonzero components of d are now stored */
/*       in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n). */
/*       The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin. */
    if (*theta != 1.) {
/*                   complete the initialization of p for theta not= one. */
	dscal_(col, theta, &p[*col + 1], &c__1);
    }
/*     Initialize GCP xcp = x. */
    dcopy_(n, &x[1], &c__1, &xcp[1], &c__1);
    if (nbreak == 0 && nfree == *n + 1) {
/*                  is a zero vector, return with the initial xcp as GCP. */
	if (*iprint > 100) {
	  printf("Cauchy X = \n");
	  i__1 = *n;
	  for (i__ = 1; i__ <= i__1; ++i__) 
	    printf("\t%e", xcp[i__]);
	  printf("\n");
	  /*	    s_wsfe(&io___105);
		    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&xcp[i__], (ftnlen)sizeof(doublereal));
	    }
	    e_wsfe();*/
	}
	return 0;
    }
/*     Initialize c = W'(xcp - x) = 0. */
    i__1 = col2;
    for (j = 1; j <= i__1; ++j) {
	c__[j] = 0.;
/* L60: */
    }
/*     Initialize derivative f2. */
    f2 = -(*theta) * f1;
    if (*col > 0) {
	bmv_(m, &sy[sy_offset], &wt[wt_offset], col, &p[1], &v[1], info);
	if (*info != 0) {
	    return 0;
	}
	f2 -= ddot_(&col2, &v[1], &c__1, &p[1], &c__1);
    }
    dtm = -f1 / f2;
    tsum = 0.;
    *nint = 1;
    if (*iprint >= 99) {
      printf("There are %i breakpoints\n", nbreak);
      /*s_wsle(&io___109);
	do_lio(&c__9, &c__1, "There are ", (ftnlen)10);
	do_lio(&c__3, &c__1, (char *)&nbreak, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, "  breakpoints ", (ftnlen)14);
	e_wsle();*/
    }
/*     If there are no breakpoints, locate the GCP and return. */
    if (nbreak == 0) {
	goto L888;
    }
    nleft = nbreak;
    iter = 1;
    tj = 0.;
/* ------------------- the beginning of the loop ------------------------- */
L777:
/*     Find the next smallest breakpoint; */
/*       compute dt = t(nleft) - t(nleft + 1). */
    tj0 = tj;
    if (iter == 1) {
/*         Since we already have the smallest breakpoint we need not do */
/*         heapsort yet. Often only one breakpoint is used and the */
/*         cost of heapsort is avoided. */
	tj = bkmin;
	ibp = iorder[ibkmin];
    } else {
	if (iter == 2) {
/*             Replace the already used smallest breakpoint with the */
/*             breakpoint numbered nbreak > nlast, before heapsort call. */
	    if (ibkmin != nbreak) {
		t[ibkmin] = t[nbreak];
		iorder[ibkmin] = iorder[nbreak];
	    }
/*        Update heap structure of breakpoints */
/*           (if iter=2, initialize heap). */
	}
	i__1 = iter - 2;
	hpsolb_(&nleft, &t[1], &iorder[1], &i__1);
	tj = t[nleft];
	ibp = iorder[nleft];
    }
    dt = tj - tj0;
    if (dt != 0. && *iprint >= 100) {
      printf("Piece %i --f1, f2 at start point %e %e\n",
	     *nint, f1, f2);
      /*	s_wsfe(&io___116);
	do_fio(&c__1, (char *)&(*nint), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&f1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&f2, (ftnlen)sizeof(doublereal));
	e_wsfe();*/

      printf("Distance to the next break point = %e\n", dt);
      printf("Distance to the stationary point = %e\n", dtm);
      /*	s_wsfe(&io___117);
	do_fio(&c__1, (char *)&dt, (ftnlen)sizeof(doublereal));
	e_wsfe();
	s_wsfe(&io___118);
	do_fio(&c__1, (char *)&dtm, (ftnlen)sizeof(doublereal));
	e_wsfe();*/
    }
/*     If a minimizer is within this interval, locate the GCP and return. */
    if (dtm < dt) {
	goto L888;
    }
/*     Otherwise fix one variable and */
/*       reset the corresponding component of d to zero. */
    tsum += dt;
    --nleft;
    ++iter;
    dibp = d__[ibp];
    d__[ibp] = 0.;
    if (dibp > 0.) {
	zibp = u[ibp] - x[ibp];
	xcp[ibp] = u[ibp];
	iwhere[ibp] = 2;
    } else {
	zibp = l[ibp] - x[ibp];
	xcp[ibp] = l[ibp];
	iwhere[ibp] = 1;
    }
    if (*iprint >= 100) {
      printf("Variable %i is fixed\n", ibp);
      /*	s_wsle(&io___121);
	do_lio(&c__9, &c__1, "Variable  ", (ftnlen)10);
	do_lio(&c__3, &c__1, (char *)&ibp, (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, "  is fixed.", (ftnlen)11);
	e_wsle();*/
    }
    if (nleft == 0 && nbreak == *n) {
/*                                             all n variables are fixed, */
/*                                                return with xcp as GCP. */
	dtm = dt;
	goto L999;
    }
/*     Update the derivative information. */
    ++(*nint);
/* Computing 2nd power */
    d__1 = dibp;
    dibp2 = d__1 * d__1;
/*     Update f1 and f2. */
/*        temporarily set f1 and f2 for col=0. */
    f1 = f1 + dt * f2 + dibp2 - *theta * dibp * zibp;
    f2 -= *theta * dibp2;
    if (*col > 0) {
/*                          update c = c + dt*p. */
	daxpy_(&col2, &dt, &p[1], &c__1, &c__[1], &c__1);
/*           choose wbp, */
/*           the row of W corresponding to the breakpoint encountered. */
	pointr = *head;
	i__1 = *col;
	for (j = 1; j <= i__1; ++j) {
	    wbp[j] = wy[ibp + pointr * wy_dim1];
	    wbp[*col + j] = *theta * ws[ibp + pointr * ws_dim1];
	    pointr = pointr % *m + 1;
/* L70: */
	}
/*           compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'. */
	bmv_(m, &sy[sy_offset], &wt[wt_offset], col, &wbp[1], &v[1], info);
	if (*info != 0) {
	    return 0;
	}
	wmc = ddot_(&col2, &c__[1], &c__1, &v[1], &c__1);
	wmp = ddot_(&col2, &p[1], &c__1, &v[1], &c__1);
	wmw = ddot_(&col2, &wbp[1], &c__1, &v[1], &c__1);
/*           update p = p - dibp*wbp. */
	d__1 = -dibp;
	daxpy_(&col2, &d__1, &wbp[1], &c__1, &p[1], &c__1);
/*           complete updating f1 and f2 while col > 0. */
	f1 += dibp * wmc;
	f2 = f2 + dibp * 2. * wmp - dibp2 * wmw;
    }
    if (nleft > 0) {
	dtm = -f1 / f2;
	goto L777;
/*                 to repeat the loop for unsearched intervals. */
    } else if (bnded) {
	f1 = 0.;
	f2 = 0.;
	dtm = 0.;
    } else {
	dtm = -f1 / f2;
    }
/* ------------------- the end of the loop ------------------------------- */
L888:
    if (*iprint >= 99) {
      printf("GCP found in this segment\n");
      /*	s_wsle(&io___126);
	e_wsle();
	s_wsle(&io___127);
	do_lio(&c__9, &c__1, "GCP found in this segment", (ftnlen)25);
	e_wsle();*/
      printf("Piece %i --f1, f2 at start point %e %e\n",
	     *nint, f1, f2);

      /*	s_wsfe(&io___128);
	do_fio(&c__1, (char *)&(*nint), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&f1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&f2, (ftnlen)sizeof(doublereal));
	e_wsfe();*/
      printf("Distance to the stationary point = %e\n", dtm);
      /*	s_wsfe(&io___129);
	do_fio(&c__1, (char *)&dtm, (ftnlen)sizeof(doublereal));
	e_wsfe();*/
    }
    if (dtm <= 0.) {
	dtm = 0.;
    }
    tsum += dtm;
/*     Move free variables (i.e., the ones w/o breakpoints) and */
/*       the variables whose breakpoints haven't been reached. */
    daxpy_(n, &tsum, &d__[1], &c__1, &xcp[1], &c__1);
L999:
/*     Update c = c + dtm*p = W'(x^c - x) */
/*       which will be used in computing r = Z'(B(x^c - x) + g). */
    if (*col > 0) {
	daxpy_(&col2, &dtm, &p[1], &c__1, &c__[1], &c__1);
    }
    if (*iprint > 100) {
      printf("Cauchy X =\n");
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__) 
	printf("\t%e", xcp[i__]);
      printf("\n");
      /*	s_wsfe(&io___130);
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, (char *)&xcp[i__], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();*/
    }
    if (*iprint >= 99) {
      printf("---------------- exit CAUCHY----------------------\n");
      /*
	s_wsfe(&io___131);
	e_wsfe();*/
    }
    return 0;
} /* cauchy_ */

/* ====================== The end of cauchy ============================== */
/* Subroutine */ int cmprlb_(n, m, x, g, ws, wy, sy, wt, z__, r__, wa, index, 
	theta, col, head, nfree, cnstnd, info)
integer *n, *m;
doublereal *x, *g, *ws, *wy, *sy, *wt, *z__, *r__, *wa;
integer *index;
doublereal *theta;
integer *col, *head, *nfree;
logical *cnstnd;
integer *info;
{
    /* System generated locals */
    integer ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset, 
	    wt_dim1, wt_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k;
    static doublereal a1, a2;
    static integer pointr;
    extern /* Subroutine */ int bmv_();

/*     ************ */

/*     Subroutine cmprlb */

/*       This subroutine computes r=-Z'B(xcp-xk)-Z'g by using */
/*         wa(2m+1)=W'(xcp-x) from subroutine cauchy. */

/*     Subprograms called: */

/*       L-BFGS-B Library ... bmv. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --index;
    --r__;
    --z__;
    --g;
    --x;
    --wa;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1 * 1;
    wt -= wt_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1 * 1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1 * 1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1 * 1;
    ws -= ws_offset;

    /* Function Body */
    if (! (*cnstnd) && *col > 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    r__[i__] = -g[i__];
/* L26: */
	}
    } else {
	i__1 = *nfree;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = index[i__];
	    r__[i__] = -(*theta) * (z__[k] - x[k]) - g[k];
/* L30: */
	}
	bmv_(m, &sy[sy_offset], &wt[wt_offset], col, &wa[(*m << 1) + 1], &wa[
		1], info);
	if (*info != 0) {
	    *info = -8;
	    return 0;
	}
	pointr = *head;
	i__1 = *col;
	for (j = 1; j <= i__1; ++j) {
	    a1 = wa[j];
	    a2 = *theta * wa[*col + j];
	    i__2 = *nfree;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		k = index[i__];
		r__[i__] = r__[i__] + wy[k + pointr * wy_dim1] * a1 + ws[k + 
			pointr * ws_dim1] * a2;
/* L32: */
	    }
	    pointr = pointr % *m + 1;
/* L34: */
	}
    }
    return 0;
} /* cmprlb_ */

/* ======================= The end of cmprlb ============================= */
/* Subroutine */ int errclb_(n, m, factr, l, u, nbd, task, info, k, task_len)
integer *n, *m;
doublereal *factr, *l, *u;
integer *nbd;
char *task;
integer *info, *k;
ftnlen task_len;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    //    /* Subroutine */ int s_copy();

    /* Local variables */
    static integer i__;

/*     ************ */

/*     Subroutine errclb */

/*     This subroutine checks the validity of the input data. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
/*     Check the input arguments for errors. */
    /* Parameter adjustments */
    --nbd;
    --u;
    --l;

    /* Function Body */
    if (*n <= 0) {
      my_strcpy(task, "ERROR: N .LE. 0");
      //	s_copy(task, "ERROR: N .LE. 0", (ftnlen)60, (ftnlen)15);
    }
    if (*m <= 0) {
      my_strcpy(task, "ERROR: M .LE. 0");
      //	s_copy(task, "ERROR: M .LE. 0", (ftnlen)60, (ftnlen)15);
    }
    if (*factr < 0.) {
      my_strcpy(task, "ERROR: FACTR .LT. 0");
      //	s_copy(task, "ERROR: FACTR .LT. 0", (ftnlen)60, (ftnlen)19);
    }
/*     Check the validity of the arrays nbd(i), u(i), and l(i). */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (nbd[i__] < 0 || nbd[i__] > 3) {
/*                                                   return */
	  my_strcpy(task, "ERROR: INVALID NBD");
	  //	    s_copy(task, "ERROR: INVALID NBD", (ftnlen)60, (ftnlen)18);
	    *info = -6;
	    *k = i__;
	}
	if (nbd[i__] == 2) {
	    if (l[i__] > u[i__]) {
/*                                    return */
	      my_strcpy(task, "ERROR: NO FEASIBLE SOLUTION");
	      //		s_copy(task, "ERROR: NO FEASIBLE SOLUTION", (ftnlen)60, (
	      //			ftnlen)27);
		*info = -7;
		*k = i__;
	    }
	}
/* L10: */
    }
    return 0;
} /* errclb_ */

/* ======================= The end of errclb ============================= */
/* Subroutine */ int formk_(n, nsub, ind, nenter, ileave, indx2, iupdat, 
	updatd, wn, wn1, m, ws, wy, sy, theta, col, head, info)
integer *n, *nsub, *ind, *nenter, *ileave, *indx2, *iupdat;
logical *updatd;
doublereal *wn, *wn1;
integer *m;
doublereal *ws, *wy, *sy, *theta;
integer *col, *head, *info;
{
    /* System generated locals */
    integer wn_dim1, wn_offset, wn1_dim1, wn1_offset, ws_dim1, ws_offset, 
	    wy_dim1, wy_offset, sy_dim1, sy_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer dend, pend;
    extern doublereal ddot_();
    static integer upcl;
    static doublereal temp1, temp2, temp3, temp4;
    static integer i__, k;
    extern /* Subroutine */ int dpofa_(), dcopy_(), dtrsl_();
    static integer ipntr, jpntr, k1, m2, dbegin, is, js, iy, jy, pbegin, is1, 
	    js1, col2;

/*     ************ */

/*     Subroutine formk */

/*     This subroutine forms  the LEL^T factorization of the indefinite */

/*       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ] */
/*                     [L_a -R_z           theta*S'AA'S ] */
/*                                                    where E = [-I  0] */
/*                                                              [ 0  I] */
/*     The matrix K can be shown to be equal to the matrix M^[-1]N */
/*       occurring in section 5.1 of [1], as well as to the matrix */
/*       Mbar^[-1] Nbar in section 5.3. */

/*     n is an integer variable. */
/*       On entry n is the dimension of the problem. */
/*       On exit n is unchanged. */

/*     nsub is an integer variable */
/*       On entry nsub is the number of subspace variables in free set. */
/*       On exit nsub is not changed. */

/*     ind is an integer array of dimension nsub. */
/*       On entry ind specifies the indices of subspace variables. */
/*       On exit ind is unchanged. */

/*     nenter is an integer variable. */
/*       On entry nenter is the number of variables entering the */
/*         free set. */
/*       On exit nenter is unchanged. */

/*     ileave is an integer variable. */
/*       On entry indx2(ileave),...,indx2(n) are the variables leaving */
/*         the free set. */
/*       On exit ileave is unchanged. */

/*     indx2 is an integer array of dimension n. */
/*       On entry indx2(1),...,indx2(nenter) are the variables entering */
/*         the free set, while indx2(ileave),...,indx2(n) are the */
/*         variables leaving the free set. */
/*       On exit indx2 is unchanged. */

/*     iupdat is an integer variable. */
/*       On entry iupdat is the total number of BFGS updates made so far. */
/*       On exit iupdat is unchanged. */

/*     updatd is a logical variable. */
/*       On entry 'updatd' is true if the L-BFGS matrix is updatd. */
/*       On exit 'updatd' is unchanged. */

/*     wn is a double precision array of dimension 2m x 2m. */
/*       On entry wn is unspecified. */
/*       On exit the upper triangle of wn stores the LEL^T factorization */
/*         of the 2*col x 2*col indefinite matrix */
/*                     [-D -Y'ZZ'Y/theta     L_a'-R_z'  ] */
/*                     [L_a -R_z           theta*S'AA'S ] */

/*     wn1 is a double precision array of dimension 2m x 2m. */
/*       On entry wn1 stores the lower triangular part of */
/*                     [Y' ZZ'Y   L_a'+R_z'] */
/*                     [L_a+R_z   S'AA'S   ] */
/*         in the previous iteration. */
/*       On exit wn1 stores the corresponding updated matrices. */
/*       The purpose of wn1 is just to store these inner products */
/*       so they can be easily updated and inserted into wn. */

/*     m is an integer variable. */
/*       On entry m is the maximum number of variable metric corrections */
/*         used to define the limited memory matrix. */
/*       On exit m is unchanged. */

/*     ws, wy, sy, and wtyy are double precision arrays; */
/*     theta is a double precision variable; */
/*     col is an integer variable; */
/*     head is an integer variable. */
/*       On entry they store the information defining the */
/*                                          limited memory BFGS matrix: */
/*         ws(n,m) stores S, a set of s-vectors; */
/*         wy(n,m) stores Y, a set of y-vectors; */
/*         sy(m,m) stores S'Y; */
/*         wtyy(m,m) stores the Cholesky factorization */
/*                                   of (theta*S'S+LD^(-1)L') */
/*         theta is the scaling factor specifying B_0 = theta I; */
/*         col is the number of variable metric corrections stored; */
/*         head is the location of the 1st s- (or y-) vector in S (or Y). */
/*       On exit they are unchanged. */

/*     info is an integer variable. */
/*       On entry info is unspecified. */
/*       On exit info =  0 for normal return; */
/*                    = -1 when the 1st Cholesky factorization failed; */
/*                    = -2 when the 2st Cholesky factorization failed. */

/*     Subprograms called: */

/*       Linpack ... dcopy, dpofa, dtrsl. */


/*     References: */
/*       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited */
/*       memory algorithm for bound rained optimization'', */
/*       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208. */

/*       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a */
/*       limited memory FORTRAN code for solving bound rained */
/*       optimization problems'', Tech. Report, NAM-11, EECS Department, */
/*       Northwestern University, 1994. */

/*       (Postscript files of these papers are available via anonymous */
/*        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.) */

/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
/*     Form the lower triangular part of */
/*               WN1 = [Y' ZZ'Y   L_a'+R_z'] */
/*                     [L_a+R_z   S'AA'S   ] */
/*        where L_a is the strictly lower triangular part of S'AA'Y */
/*              R_z is the upper triangular part of S'ZZ'Y. */
    /* Parameter adjustments */
    --indx2;
    --ind;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1 * 1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1 * 1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1 * 1;
    ws -= ws_offset;
    wn1_dim1 = 2 * *m;
    wn1_offset = 1 + wn1_dim1 * 1;
    wn1 -= wn1_offset;
    wn_dim1 = 2 * *m;
    wn_offset = 1 + wn_dim1 * 1;
    wn -= wn_offset;

    /* Function Body */
    if (*updatd) {
	if (*iupdat > *m) {
/*                                 shift old part of WN1. */
	    i__1 = *m - 1;
	    for (jy = 1; jy <= i__1; ++jy) {
		js = *m + jy;
		i__2 = *m - jy;
		dcopy_(&i__2, &wn1[jy + 1 + (jy + 1) * wn1_dim1], &c__1, &wn1[
			jy + jy * wn1_dim1], &c__1);
		i__2 = *m - jy;
		dcopy_(&i__2, &wn1[js + 1 + (js + 1) * wn1_dim1], &c__1, &wn1[
			js + js * wn1_dim1], &c__1);
		i__2 = *m - 1;
		dcopy_(&i__2, &wn1[*m + 2 + (jy + 1) * wn1_dim1], &c__1, &wn1[
			*m + 1 + jy * wn1_dim1], &c__1);
/* L10: */
	    }
	}
/*          put new rows in blocks (1,1), (2,1) and (2,2). */
	pbegin = 1;
	pend = *nsub;
	dbegin = *nsub + 1;
	dend = *n;
	iy = *col;
	is = *m + *col;
	ipntr = *head + *col - 1;
	if (ipntr > *m) {
	    ipntr -= *m;
	}
	jpntr = *head;
	i__1 = *col;
	for (jy = 1; jy <= i__1; ++jy) {
	    js = *m + jy;
	    temp1 = 0.;
	    temp2 = 0.;
	    temp3 = 0.;
/*             compute element jy of row 'col' of Y'ZZ'Y */
	    i__2 = pend;
	    for (k = pbegin; k <= i__2; ++k) {
		k1 = ind[k];
		temp1 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
/* L15: */
	    }
/*             compute elements jy of row 'col' of L_a and S'AA'S */
	    i__2 = dend;
	    for (k = dbegin; k <= i__2; ++k) {
		k1 = ind[k];
		temp2 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
		temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
/* L16: */
	    }
	    wn1[iy + jy * wn1_dim1] = temp1;
	    wn1[is + js * wn1_dim1] = temp2;
	    wn1[is + jy * wn1_dim1] = temp3;
	    jpntr = jpntr % *m + 1;
/* L20: */
	}
/*          put new column in block (2,1). */
	jy = *col;
	jpntr = *head + *col - 1;
	if (jpntr > *m) {
	    jpntr -= *m;
	}
	ipntr = *head;
	i__1 = *col;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    is = *m + i__;
	    temp3 = 0.;
/*             compute element i of column 'col' of R_z */
	    i__2 = pend;
	    for (k = pbegin; k <= i__2; ++k) {
		k1 = ind[k];
		temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
/* L25: */
	    }
	    ipntr = ipntr % *m + 1;
	    wn1[is + jy * wn1_dim1] = temp3;
/* L30: */
	}
	upcl = *col - 1;
    } else {
	upcl = *col;
    }
/*       modify the old parts in blocks (1,1) and (2,2) due to changes */
/*       in the set of free variables. */
    ipntr = *head;
    i__1 = upcl;
    for (iy = 1; iy <= i__1; ++iy) {
	is = *m + iy;
	jpntr = *head;
	i__2 = iy;
	for (jy = 1; jy <= i__2; ++jy) {
	    js = *m + jy;
	    temp1 = 0.;
	    temp2 = 0.;
	    temp3 = 0.;
	    temp4 = 0.;
	    i__3 = *nenter;
	    for (k = 1; k <= i__3; ++k) {
		k1 = indx2[k];
		temp1 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
		temp2 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
/* L35: */
	    }
	    i__3 = *n;
	    for (k = *ileave; k <= i__3; ++k) {
		k1 = indx2[k];
		temp3 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
		temp4 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
/* L36: */
	    }
	    wn1[iy + jy * wn1_dim1] = wn1[iy + jy * wn1_dim1] + temp1 - temp3;
	    wn1[is + js * wn1_dim1] = wn1[is + js * wn1_dim1] - temp2 + temp4;
	    jpntr = jpntr % *m + 1;
/* L40: */
	}
	ipntr = ipntr % *m + 1;
/* L45: */
    }
/*       modify the old parts in block (2,1). */
    ipntr = *head;
    i__1 = *m + upcl;
    for (is = *m + 1; is <= i__1; ++is) {
	jpntr = *head;
	i__2 = upcl;
	for (jy = 1; jy <= i__2; ++jy) {
	    temp1 = 0.;
	    temp3 = 0.;
	    i__3 = *nenter;
	    for (k = 1; k <= i__3; ++k) {
		k1 = indx2[k];
		temp1 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
/* L50: */
	    }
	    i__3 = *n;
	    for (k = *ileave; k <= i__3; ++k) {
		k1 = indx2[k];
		temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
/* L51: */
	    }
	    if (is <= jy + *m) {
		wn1[is + jy * wn1_dim1] = wn1[is + jy * wn1_dim1] + temp1 - 
			temp3;
	    } else {
		wn1[is + jy * wn1_dim1] = wn1[is + jy * wn1_dim1] - temp1 + 
			temp3;
	    }
	    jpntr = jpntr % *m + 1;
/* L55: */
	}
	ipntr = ipntr % *m + 1;
/* L60: */
    }
/*     Form the upper triangle of WN = [D+Y' ZZ'Y/theta   -L_a'+R_z' ] */
/*                                     [-L_a +R_z        S'AA'S*theta] */
    m2 = *m << 1;
    i__1 = *col;
    for (iy = 1; iy <= i__1; ++iy) {
	is = *col + iy;
	is1 = *m + iy;
	i__2 = iy;
	for (jy = 1; jy <= i__2; ++jy) {
	    js = *col + jy;
	    js1 = *m + jy;
	    wn[jy + iy * wn_dim1] = wn1[iy + jy * wn1_dim1] / *theta;
	    wn[js + is * wn_dim1] = wn1[is1 + js1 * wn1_dim1] * *theta;
/* L65: */
	}
	i__2 = iy - 1;
	for (jy = 1; jy <= i__2; ++jy) {
	    wn[jy + is * wn_dim1] = -wn1[is1 + jy * wn1_dim1];
/* L66: */
	}
	i__2 = *col;
	for (jy = iy; jy <= i__2; ++jy) {
	    wn[jy + is * wn_dim1] = wn1[is1 + jy * wn1_dim1];
/* L67: */
	}
	wn[iy + iy * wn_dim1] += sy[iy + iy * sy_dim1];
/* L70: */
    }
/*     Form the upper triangle of WN= [  LL'            L^-1(-L_a'+R_z')] */
/*                                    [(-L_a +R_z)L'^-1   S'AA'S*theta  ] */
/*        first Cholesky factor (1,1) block of wn to get LL' */
/*                          with L' stored in the upper triangle of wn. */
    dpofa_(&wn[wn_offset], &m2, col, info);
    if (*info != 0) {
	*info = -1;
	return 0;
    }
/*        then form L^-1(-L_a'+R_z') in the (1,2) block. */
    col2 = *col << 1;
    i__1 = col2;
    for (js = *col + 1; js <= i__1; ++js) {
	dtrsl_(&wn[wn_offset], &m2, col, &wn[js * wn_dim1 + 1], &c__11, info);
/* L71: */
    }
/*     Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the */
/*        upper triangle of (2,2) block of wn. */
    i__1 = col2;
    for (is = *col + 1; is <= i__1; ++is) {
	i__2 = col2;
	for (js = is; js <= i__2; ++js) {
	    wn[is + js * wn_dim1] += ddot_(col, &wn[is * wn_dim1 + 1], &c__1, 
		    &wn[js * wn_dim1 + 1], &c__1);
/* L74: */
	}
/* L72: */
    }
/*     Cholesky factorization of (2,2) block of wn. */
    dpofa_(&wn[*col + 1 + (*col + 1) * wn_dim1], &m2, col, info);
    if (*info != 0) {
	*info = -2;
	return 0;
    }
    return 0;
} /* formk_ */

/* ======================= The end of formk ============================== */
/* Subroutine */ int formt_(m, wt, sy, ss, col, theta, info)
integer *m;
doublereal *wt, *sy, *ss;
integer *col;
doublereal *theta;
integer *info;
{
    /* System generated locals */
    integer wt_dim1, wt_offset, sy_dim1, sy_offset, ss_dim1, ss_offset, i__1, 
	    i__2, i__3;

    /* Local variables */
    static doublereal ddum;
    static integer i__, j, k;
    extern /* Subroutine */ int dpofa_();
    static integer k1;

/*     ************ */

/*     Subroutine formt */

/*       This subroutine forms the upper half of the pos. def. and symm. */
/*         T = theta*SS + L*D^(-1)*L', stores T in the upper triangle */
/*         of the array wt, and performs the Cholesky factorization of T */
/*         to produce J*J', with J' stored in the upper triangle of wt. */

/*     Subprograms called: */

/*       Linpack ... dpofa. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
/*     Form the upper half of  T = theta*SS + L*D^(-1)*L', */
/*        store T in the upper triangle of the array wt. */
    /* Parameter adjustments */
    ss_dim1 = *m;
    ss_offset = 1 + ss_dim1 * 1;
    ss -= ss_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1 * 1;
    sy -= sy_offset;
    wt_dim1 = *m;
    wt_offset = 1 + wt_dim1 * 1;
    wt -= wt_offset;

    /* Function Body */
    i__1 = *col;
    for (j = 1; j <= i__1; ++j) {
	wt[j * wt_dim1 + 1] = *theta * ss[j * ss_dim1 + 1];
/* L52: */
    }
    i__1 = *col;
    for (i__ = 2; i__ <= i__1; ++i__) {
	i__2 = *col;
	for (j = i__; j <= i__2; ++j) {
	    k1 = min(i__,j) - 1;
	    ddum = 0.;
	    i__3 = k1;
	    for (k = 1; k <= i__3; ++k) {
		ddum += sy[i__ + k * sy_dim1] * sy[j + k * sy_dim1] / sy[k + 
			k * sy_dim1];
/* L53: */
	    }
	    wt[i__ + j * wt_dim1] = ddum + *theta * ss[i__ + j * ss_dim1];
/* L54: */
	}
/* L55: */
    }
/*     Cholesky factorize T to J*J' with */
/*        J' stored in the upper triangle of wt. */
    dpofa_(&wt[wt_offset], m, col, info);
    if (*info != 0) {
	*info = -3;
    }
    return 0;
} /* formt_ */

/* ======================= The end of formt ============================== */
/* Subroutine */ int freev_(n, nfree, index, nenter, ileave, indx2, iwhere, 
	wrk, updatd, cnstnd, iprint, iter)
integer *n, *nfree, *index, *nenter, *ileave, *indx2, *iwhere;
logical *wrk, *updatd, *cnstnd;
integer *iprint, *iter;
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    static integer iact, i__, k;

    /* Fortran I/O blocks */
    /*    static cilist io___168 = { 0, 6, 0, 0, 0 };
    static cilist io___169 = { 0, 6, 0, 0, 0 };
    static cilist io___170 = { 0, 6, 0, 0, 0 };
    static cilist io___172 = { 0, 6, 0, 0, 0 };*/


/*     ************ */

/*     Subroutine freev */

/*     This subroutine counts the entering and leaving variables when */
/*       iter > 0, and finds the index set of free and active variables */
/*       at the GCP. */

/*     cnstnd is a logical variable indicating whether bounds are present */

/*     index is an integer array of dimension n */
/*       for i=1,...,nfree, index(i) are the indices of free variables */
/*       for i=nfree+1,...,n, index(i) are the indices of bound variables */
/*       On entry after the first iteration, index gives */
/*         the free variables at the previous iteration. */
/*       On exit it gives the free variables based on the determination */
/*         in cauchy using the array iwhere. */

/*     indx2 is an integer array of dimension n */
/*       On entry indx2 is unspecified. */
/*       On exit with iter>0, indx2 indicates which variables */
/*          have changed status since the previous iteration. */
/*       For i= 1,...,nenter, indx2(i) have changed from bound to free. */
/*       For i= ileave+1,...,n, indx2(i) have changed from free to bound. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --iwhere;
    --indx2;
    --index;

    /* Function Body */
    *nenter = 0;
    *ileave = *n + 1;
    if (*iter > 0 && *cnstnd) {
/*                           count the entering and leaving variables. */
	i__1 = *nfree;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    k = index[i__];
	    if (iwhere[k] > 0) {
		--(*ileave);
		indx2[*ileave] = k;
		if (*iprint >= 100) {
		  printf("Variable %i leaves the set of free variables\n",
			 k);
		  /*
		    s_wsle(&io___168);
		    do_lio(&c__9, &c__1, "Variable ", (ftnlen)9);
		    do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
		    do_lio(&c__9, &c__1, " leaves the set of free variables", 
			    (ftnlen)33);
			    e_wsle();*/
		}
	    }
/* L20: */
	}
	i__1 = *n;
	for (i__ = *nfree + 1; i__ <= i__1; ++i__) {
	    k = index[i__];
	    if (iwhere[k] <= 0) {
		++(*nenter);
		indx2[*nenter] = k;
		if (*iprint >= 100) {
		  printf("Variable %i enters the set of free variables\n",k);
		  /*		    s_wsle(&io___169);
		    do_lio(&c__9, &c__1, "Variable ", (ftnlen)9);
		    do_lio(&c__3, &c__1, (char *)&k, (ftnlen)sizeof(integer));
		    do_lio(&c__9, &c__1, " enters the set of free variables", 
			    (ftnlen)33);
			    e_wsle();*/
		}
	    }
/* L22: */
	}
	if (*iprint >= 99) {
	  i__1 = *n + 1 - *ileave;
	  printf("%i variables leave; %i variables enter\n",
		 i__1, *nenter);
	  /*	    s_wsle(&io___170);
	    i__1 = *n + 1 - *ileave;
	    do_lio(&c__3, &c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " variables leave; ", (ftnlen)18);
	    do_lio(&c__3, &c__1, (char *)&(*nenter), (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, " variables enter", (ftnlen)16);
	    e_wsle();*/
	}
    }
    *wrk = *ileave < *n + 1 || *nenter > 0 || *updatd;
/*     Find the index set of free and active variables at the GCP. */
    *nfree = 0;
    iact = *n + 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (iwhere[i__] <= 0) {
	    ++(*nfree);
	    index[*nfree] = i__;
	} else {
	    --iact;
	    index[iact] = i__;
	}
/* L24: */
    }
    if (*iprint >= 99) {
      i__1 = *iter + 1;
      printf(" %i variables are free at GCP %i\n", *nfree, i__1);
      /*	s_wsle(&io___172);
	do_lio(&c__3, &c__1, (char *)&(*nfree), (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " variables are free at GCP ", (ftnlen)27);
	i__1 = *iter + 1;
	do_lio(&c__3, &c__1, (char *)&i__1, (ftnlen)sizeof(integer));
	e_wsle();*/
    }
    return 0;
} /* freev_ */

/* ======================= The end of freev ============================== */
/* Subroutine */ int hpsolb_(n, t, iorder, iheap)
integer *n;
doublereal *t;
integer *iorder, *iheap;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static doublereal ddum;
    static integer i__, j, k, indxin, indxou;
    static doublereal out;

/*     ************ */

/*     Subroutine hpsolb */

/*     This subroutine sorts out the least element of t, and puts the */
/*       remaining elements of t in a heap. */

/*     n is an integer variable. */
/*       On entry n is the dimension of the arrays t and iorder. */
/*       On exit n is unchanged. */

/*     t is a double precision array of dimension n. */
/*       On entry t stores the elements to be sorted, */
/*       On exit t(n) stores the least elements of t, and t(1) to t(n-1) */
/*         stores the remaining elements in the form of a heap. */

/*     iorder is an integer array of dimension n. */
/*       On entry iorder(i) is the index of t(i). */
/*       On exit iorder(i) is still the index of t(i), but iorder may be */
/*         permuted in accordance with t. */

/*     iheap is an integer variable specifying the task. */
/*       On entry iheap should be set as follows: */
/*         iheap .eq. 0 if t(1) to t(n) is not in the form of a heap, */
/*         iheap .ne. 0 if otherwise. */
/*       On exit iheap is unchanged. */


/*     References: */
/*       Algorithm 232 of CACM (J. W. J. Williams): HEAPSORT. */

/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */

/*     ************ */
    /* Parameter adjustments */
    --iorder;
    --t;

    /* Function Body */
    if (*iheap == 0) {
/*        Rearrange the elements t(1) to t(n) to form a heap. */
	i__1 = *n;
	for (k = 2; k <= i__1; ++k) {
	    ddum = t[k];
	    indxin = iorder[k];
/*           Add ddum to the heap. */
	    i__ = k;
L10:
	    if (i__ > 1) {
		j = i__ / 2;
		if (ddum < t[j]) {
		    t[i__] = t[j];
		    iorder[i__] = iorder[j];
		    i__ = j;
		    goto L10;
		}
	    }
	    t[i__] = ddum;
	    iorder[i__] = indxin;
/* L20: */
	}
    }
/*     Assign to 'out' the value of t(1), the least member of the heap, */
/*        and rearrange the remaining members to form a heap as */
/*        elements 1 to n-1 of t. */
    if (*n > 1) {
	i__ = 1;
	out = t[1];
	indxou = iorder[1];
	ddum = t[*n];
	indxin = iorder[*n];
/*        Restore the heap */
L30:
	j = i__ + i__;
	if (j <= *n - 1) {
	    if (t[j + 1] < t[j]) {
		++j;
	    }
	    if (t[j] < ddum) {
		t[i__] = t[j];
		iorder[i__] = iorder[j];
		i__ = j;
		goto L30;
	    }
	}
	t[i__] = ddum;
	iorder[i__] = indxin;
/*     Put the least member in t(n). */
	t[*n] = out;
	iorder[*n] = indxou;
    }
    return 0;
} /* hpsolb_ */

/* ====================== The end of hpsolb ============================== */
/* Subroutine */ int lnsrlb_(n, l, u, nbd, x, f, fold, gd, gdold, g, d__, r__,
	 t, z__, stp, dnorm, dtd, xstep, stpmx, iter, ifun, iback, nfgv, info,
	 task, boxed, cnstnd, csave, isave, dsave, task_len, csave_len)
integer *n;
doublereal *l, *u;
integer *nbd;
doublereal *x, *f, *fold, *gd, *gdold, *g, *d__, *r__, *t, *z__, *stp, *dnorm,
	 *dtd, *xstep, *stpmx;
integer *iter, *ifun, *iback, *nfgv, *info;
char *task;
logical *boxed, *cnstnd;
char *csave;
integer *isave;
doublereal *dsave;
ftnlen task_len;
ftnlen csave_len;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    //    integer s_cmp();
    double sqrt();
    //    /* Subroutine */ int s_copy();

    /* Local variables */
    extern doublereal ddot_();
    static integer i__;
    extern /* Subroutine */ int dcopy_();
    static doublereal a1, a2;
    extern /* Subroutine */ int dcsrch_();

/*     ********** */

/*     Subroutine lnsrlb */

/*     This subroutine calls subroutine dcsrch from the Minpack2 library */
/*       to perform the line search.  Subroutine dscrch is safeguarded so */
/*       that all trial points lie within the feasible region. */

/*     Subprograms called: */

/*       Minpack2 Library ... dcsrch. */

/*       Linpack ... dtrsl, ddot. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ********** */
    /* Parameter adjustments */
    --z__;
    --t;
    --r__;
    --d__;
    --g;
    --x;
    --nbd;
    --u;
    --l;
    --isave;
    --dsave;

    /* Function Body */
    if (strncmp(task, "FG_LN", 5)==0) {
	goto L556;
    }
    *dtd = ddot_(n, &d__[1], &c__1, &d__[1], &c__1);
    *dnorm = sqrt(*dtd);
/*     Determine the maximum step length. */
    *stpmx = 1e10;
    if (*cnstnd) {
	if (*iter == 0) {
	    *stpmx = 1.;
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		a1 = d__[i__];
		if (nbd[i__] != 0) {
		    if (a1 < 0. && nbd[i__] <= 2) {
			a2 = l[i__] - x[i__];
			if (a2 >= 0.) {
			    *stpmx = 0.;
			} else if (a1 * *stpmx < a2) {
			    *stpmx = a2 / a1;
			}
		    } else if (a1 > 0. && nbd[i__] >= 2) {
			a2 = u[i__] - x[i__];
			if (a2 <= 0.) {
			    *stpmx = 0.;
			} else if (a1 * *stpmx > a2) {
			    *stpmx = a2 / a1;
			}
		    }
		}
/* L43: */
	    }
	}
    }
    if (*iter == 0 && ! (*boxed)) {
/* Computing MIN */
	d__1 = 1. / *dnorm;
	*stp = min(d__1,*stpmx);
    } else {
	*stp = 1.;
    }
    dcopy_(n, &x[1], &c__1, &t[1], &c__1);
    dcopy_(n, &g[1], &c__1, &r__[1], &c__1);
    *fold = *f;
    *ifun = 0;
    *iback = 0;
    my_strcpy(csave, "START");
    //    s_copy(csave, "START", (ftnlen)60, (ftnlen)5);
L556:
    *gd = ddot_(n, &g[1], &c__1, &d__[1], &c__1);
    if (*ifun == 0) {
	*gdold = *gd;
	if (*gd >= 0.) {
/*                               the directional derivative >=0. */
/*                               Line search is impossible. */
	    *info = -4;
	    return 0;
	}
    }
    dcsrch_(f, gd, stp, &c_b275, &c_b276, &c_b277, &c_b9, stpmx, csave, &
	    isave[1], &dsave[1], (ftnlen)60);
    *xstep = *stp * *dnorm;
    if (strncmp(csave, "CONV", 4) != 0 && strncmp(csave, "WARN", 4)!=0) {
	my_strcpy(task, "FG_LNSRCH");
	++(*ifun);
	++(*nfgv);
	*iback = *ifun - 1;
	if (*stp == 1.) {
	    dcopy_(n, &z__[1], &c__1, &x[1], &c__1);
	} else {
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		x[i__] = *stp * d__[i__] + t[i__];
/* L41: */
	    }
	}
    } else {
	my_strcpy(task, "NEW_X");
    }
    return 0;
} /* lnsrlb_ */

/* ======================= The end of lnsrlb ============================= */
/* Subroutine */ int matupd_(n, m, ws, wy, sy, ss, d__, r__, itail, iupdat, 
	col, head, theta, rr, dr, stp, dtd)
integer *n, *m;
doublereal *ws, *wy, *sy, *ss, *d__, *r__;
integer *itail, *iupdat, *col, *head;
doublereal *theta, *rr, *dr, *stp, *dtd;
{
    /* System generated locals */
    integer ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset, 
	    ss_dim1, ss_offset, i__1, i__2;

    /* Local variables */
    extern doublereal ddot_();
    static integer j;
    extern /* Subroutine */ int dcopy_();
    static integer pointr;

/*     ************ */

/*     Subroutine matupd */

/*       This subroutine updates matrices WS and WY, and forms the */
/*         middle matrix in B. */

/*     Subprograms called: */

/*       Linpack ... dcopy, ddot. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
/*     Set pointers for matrices WS and WY. */
    /* Parameter adjustments */
    --r__;
    --d__;
    ss_dim1 = *m;
    ss_offset = 1 + ss_dim1 * 1;
    ss -= ss_offset;
    sy_dim1 = *m;
    sy_offset = 1 + sy_dim1 * 1;
    sy -= sy_offset;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1 * 1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1 * 1;
    ws -= ws_offset;

    /* Function Body */
    if (*iupdat <= *m) {
	*col = *iupdat;
	*itail = (*head + *iupdat - 2) % *m + 1;
    } else {
	*itail = *itail % *m + 1;
	*head = *head % *m + 1;
    }
/*     Update matrices WS and WY. */
    dcopy_(n, &d__[1], &c__1, &ws[*itail * ws_dim1 + 1], &c__1);
    dcopy_(n, &r__[1], &c__1, &wy[*itail * wy_dim1 + 1], &c__1);
/*     Set theta=yy/ys. */
    *theta = *rr / *dr;
/*     Form the middle matrix in B. */
/*        update the upper triangle of SS, */
/*                                         and the lower triangle of SY: */
    if (*iupdat > *m) {
/*                              move old information */
	i__1 = *col - 1;
	for (j = 1; j <= i__1; ++j) {
	    dcopy_(&j, &ss[(j + 1) * ss_dim1 + 2], &c__1, &ss[j * ss_dim1 + 1]
		    , &c__1);
	    i__2 = *col - j;
	    dcopy_(&i__2, &sy[j + 1 + (j + 1) * sy_dim1], &c__1, &sy[j + j * 
		    sy_dim1], &c__1);
/* L50: */
	}
    }
/*        add new information: the last row of SY */
/*                                             and the last column of SS: */
    pointr = *head;
    i__1 = *col - 1;
    for (j = 1; j <= i__1; ++j) {
	sy[*col + j * sy_dim1] = ddot_(n, &d__[1], &c__1, &wy[pointr * 
		wy_dim1 + 1], &c__1);
	ss[j + *col * ss_dim1] = ddot_(n, &ws[pointr * ws_dim1 + 1], &c__1, &
		d__[1], &c__1);
	pointr = pointr % *m + 1;
/* L51: */
    }
    if (*stp == 1.) {
	ss[*col + *col * ss_dim1] = *dtd;
    } else {
	ss[*col + *col * ss_dim1] = *stp * *stp * *dtd;
    }
    sy[*col + *col * sy_dim1] = *dr;
    return 0;
} /* matupd_ */

/* ======================= The end of matupd ============================= */
/* Subroutine */ int prn1lb_(n, m, l, u, x, iprint, itfile, epsmch)
integer *n, *m;
doublereal *l, *u, *x;
integer *iprint, *itfile;
doublereal *epsmch;
{
    /* Format strings */
  /*    static char fmt_7001[] = "(\002RUNNING THE L-BFGS-B CODE\002,/,/,\002   \
        * * *\002,/,/,\002Machine precision =\002,1p,d10.3)";
    static char fmt_2001[] = "(\002RUNNING THE L-BFGS-B CODE\002,/,/,\002it \
   = iteration number\002,/,\002nf    = number of function evaluations\002,/,\
\002nint  = number of segments explored during the Cauchy search\002,/,\002n\
act  = number of active bounds at the generalized Cauchy point\002,/,\002sub\
   = manner in which the subspace minimization terminated:\002,/,\002       \
 con = converged, bnd = a bound was reached\002,/,\002itls  = number of iter\
ations performed in the line search\002,/,\002stepl = step length used\002,/,\
\002tstep = norm of the displacement (total step)\002,/,\002projg = norm of \
the projected gradient\002,/,\002f     = function value\002,/,/,\002        \
   * * *\002,/,/,\002Machine precision =\002,1p,d10.3)";
    static char fmt_9001[] = "(/,3x,\002it\002,3x,\002nf\002,2x,\002nint\002\
,2x,\002nact\002,2x,\002sub\002,2x,\002itls\002,2x,\002stepl\002,4x,\002tstep\
\002,5x,\002projg\002,8x,\002f\002)";
static char fmt_1004[] = "(/,a4,1p,6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))";*/

    /* System generated locals */
  //    integer i__1;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe(), s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    //    static integer i__;

    /* Fortran I/O blocks */
    /*    static cilist io___185 = { 0, 6, 0, fmt_7001, 0 };
    static cilist io___186 = { 0, 6, 0, 0, 0 };
    static cilist io___187 = { 0, 0, 0, fmt_2001, 0 };
    static cilist io___188 = { 0, 0, 0, 0, 0 };
    static cilist io___189 = { 0, 0, 0, fmt_9001, 0 };
    static cilist io___190 = { 0, 6, 0, fmt_1004, 0 };
    static cilist io___192 = { 0, 6, 0, fmt_1004, 0 };
    static cilist io___193 = { 0, 6, 0, fmt_1004, 0 };*/


/*     ************ */

/*     Subroutine prn1lb */

/*     This subroutine prints the input data, initial point, upper and */
/*       lower bounds of each variable, machine precision, as well as */
/*       the headings of the output. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --x;
    --u;
    --l;

    /* Function Body */
    if (*iprint >= 0) {
      printf("RUNNING THE L-BFGS-B CODE\n\n\t* * *\n\nMachine precision=%e\n",
	     *epsmch);
      /*
	s_wsfe(&io___185);
	do_fio(&c__1, (char *)&(*epsmch), (ftnlen)sizeof(doublereal));
	e_wsfe();*/
      printf("N = %i    M = %i\n", *n, *m);
      /*
	s_wsle(&io___186);
	do_lio(&c__9, &c__1, "N = ", (ftnlen)4);
	do_lio(&c__3, &c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, "    M = ", (ftnlen)8);
	do_lio(&c__3, &c__1, (char *)&(*m), (ftnlen)sizeof(integer));
	e_wsle();*/
      /*	if (*iprint >= 1) {

	    io___187.ciunit = *itfile;
	    s_wsfe(&io___187);
	    do_fio(&c__1, (char *)&(*epsmch), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    io___188.ciunit = *itfile;
	    s_wsle(&io___188);
	    do_lio(&c__9, &c__1, "N = ", (ftnlen)4);
	    do_lio(&c__3, &c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	    do_lio(&c__9, &c__1, "    M = ", (ftnlen)8);
	    do_lio(&c__3, &c__1, (char *)&(*m), (ftnlen)sizeof(integer));
	    e_wsle();
	    io___189.ciunit = *itfile;
	    s_wsfe(&io___189);
	    e_wsfe();
	    if (*iprint > 100) {
		s_wsfe(&io___190);
		do_fio(&c__1, "L =", (ftnlen)3);
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    do_fio(&c__1, (char *)&l[i__], (ftnlen)sizeof(doublereal))
			    ;
		}
		e_wsfe();
		s_wsfe(&io___192);
		do_fio(&c__1, "X0 =", (ftnlen)4);
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal))
			    ;
		}
		e_wsfe();
		s_wsfe(&io___193);
		do_fio(&c__1, "U =", (ftnlen)3);
		i__1 = *n;
		for (i__ = 1; i__ <= i__1; ++i__) {
		    do_fio(&c__1, (char *)&u[i__], (ftnlen)sizeof(doublereal))
			    ;
		}
		e_wsfe();
	    }
	    }*/
    }
    return 0;
} /* prn1lb_ */

/* ======================= The end of prn1lb ============================= */
/* Subroutine */ int prn2lb_(n, x, f, g, iprint, itfile, iter, nfgv, nact, 
	sbgnrm, nint, word, iword, iback, stp, xstep, word_len)
integer *n;
doublereal *x, *f, *g;
integer *iprint, *itfile, *iter, *nfgv, *nact;
doublereal *sbgnrm;
integer *nint;
char *word;
integer *iword, *iback;
doublereal *stp, *xstep;
ftnlen word_len;
{
    /* Format strings */
  /*    static char fmt_2001[] = "(/,\002At iterate\002,i5,4x,\002f= \002,1p,d12\
.5,4x,\002|proj g|= \002,1p,d12.5)";
    static char fmt_1004[] = "(/,a4,1p,6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))";
    static char fmt_3001[] = "(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),1\
    p,2(1x,d10.3))";*/

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    //    /* Subroutine */ int s_copy();
    integer s_wsle(), do_lio(), e_wsle(), s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static integer imod, i__;

    /* Fortran I/O blocks */
    /*    static cilist io___194 = { 0, 6, 0, 0, 0 };
    static cilist io___195 = { 0, 6, 0, fmt_2001, 0 };
    static cilist io___196 = { 0, 6, 0, fmt_1004, 0 };
    static cilist io___198 = { 0, 6, 0, fmt_1004, 0 };
    static cilist io___200 = { 0, 6, 0, fmt_2001, 0 };
    static cilist io___201 = { 0, 0, 0, fmt_3001, 0 };*/


/*     ************ */

/*     Subroutine prn2lb */

/*     This subroutine prints out new information after a successful */
/*       line search. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
/*           'word' records the status of subspace solutions. */
    /* Parameter adjustments */
    --g;
    --x;

    /* Function Body */
    if (*iword == 0) {
/*                            the subspace minimization converged. */
      my_strcpy(word, "con");
      //	s_copy(word, "con", (ftnlen)3, (ftnlen)3);
    } else if (*iword == 1) {
      /*                          the subspace minimization stopped at a bound. */
      my_strcpy(word, "bnd");
      //	s_copy(word, "bnd", (ftnlen)3, (ftnlen)3);
    } else if (*iword == 5) {
/*                             the truncated Newton step has been used. */
      my_strcpy(word, "TNT");
      //	s_copy(word, "TNT", (ftnlen)3, (ftnlen)3);
    } else {
      my_strcpy(word, "---");
      //	s_copy(word, "---", (ftnlen)3, (ftnlen)3);
    }
    if (*iprint >= 99) {
      printf("LINE SEARCH %i times; norm of step = %e\n",
	     *iback, *xstep);
      /*	s_wsle(&io___194);
	do_lio(&c__9, &c__1, "LINE SEARCH", (ftnlen)11);
	do_lio(&c__3, &c__1, (char *)&(*iback), (ftnlen)sizeof(integer));
	do_lio(&c__9, &c__1, " times; norm of step = ", (ftnlen)23);
	do_lio(&c__5, &c__1, (char *)&(*xstep), (ftnlen)sizeof(doublereal));
	e_wsle();*/
      printf("At iterate %i f=%e |proj g| = %e\n",
	     *iter, *f, *sbgnrm);
      /*	s_wsfe(&io___195);
	do_fio(&c__1, (char *)&(*iter), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*f), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*sbgnrm), (ftnlen)sizeof(doublereal));
	e_wsfe();*/
	if (*iprint > 100) {
	  printf("\nX=");
	  i__1 = *n;
	  for (i__ = 1; i__ <= i__1; ++i__) 
	    printf("\t%e", x[i__]);
	  printf("\n");
	  /*
	    s_wsfe(&io___196);
	    do_fio(&c__1, "X =", (ftnlen)3);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	    }
	    e_wsfe();*/
	  printf("\nG=");
	  i__1 = *n;
	  for (i__1 = 1; i__ <= i__1; ++i__) 
	    printf("\t%e", g[i__]);
	  printf("\n");
	  /*
	    s_wsfe(&io___198);
	    do_fio(&c__1, "G =", (ftnlen)3);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&g[i__], (ftnlen)sizeof(doublereal));
	    }
	    e_wsfe();*/
	}
    } else if (*iprint > 0) {
	imod = *iter % *iprint;
	if (imod == 0) {
	  /*	    s_wsfe(&io___200);
	    do_fio(&c__1, (char *)&(*iter), (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&(*f), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*sbgnrm), (ftnlen)sizeof(doublereal));
	    e_wsfe();*/
	}
    }
    if (*iprint >= 1) {
      /*printf("%i\t%i\t%i\t%i\t%s\t%i\t%e\t%e\t%e\t%e\n",
	     iter, nfgv, nint, nact, word, back, stp, xstep, sbgnrm, f);
	     io___201.ciunit = *itfile;
	s_wsfe(&io___201);
	do_fio(&c__1, (char *)&(*iter), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*nfgv), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*nint), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*nact), (ftnlen)sizeof(integer));
	do_fio(&c__1, word, (ftnlen)3);
	do_fio(&c__1, (char *)&(*iback), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*stp), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*xstep), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*sbgnrm), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*f), (ftnlen)sizeof(doublereal));
	e_wsfe();*/
    }
    return 0;
} /* prn2lb_ */

/* ======================= The end of prn2lb ============================= */
/* Subroutine */ int prn3lb_(n, x, f, task, iprint, info, itfile, iter, nfgv, 
	nintol, nskip, nact, sbgnrm, time, nint, word, iback, stp, xstep, k, 
	cachyt, sbtime, lnscht, task_len, word_len)
integer *n;
doublereal *x, *f;
char *task;
integer *iprint, *info, *itfile, *iter, *nfgv, *nintol, *nskip, *nact;
doublereal *sbgnrm, *time;
integer *nint;
char *word;
integer *iback;
doublereal *stp, *xstep;
integer *k;
doublereal *cachyt, *sbtime, *lnscht;
ftnlen task_len;
ftnlen word_len;
{
    /* Format strings */
  /*    static char fmt_3003[] = "(/,\002           * * *\002,/,/,\002Tit   = to\
tal number of iterations\002,/,\002Tnf   = total number of function evaluati\
ons\002,/,\002Tnint = total number of segments explored during\002,\002 Cauc\
hy searches\002,/,\002Skip  = number of BFGS updates skipped\002,/,\002Nact \
 = number of active bounds at final generalized\002,\002 Cauchy point\002,/\
,\002Projg = norm of the final projected gradient\002,/,\002F     = final fu\
nction value\002,/,/,\002           * * *\002)";
    static char fmt_3004[] = "(/,3x,\002N\002,3x,\002Tit\002,2x,\002Tnf\002,\
2x,\002Tnint\002,2x,\002Skip\002,2x,\002Nact\002,5x,\002Projg\002,8x,\002\
F\002)";
    static char fmt_3005[] = "(i5,2(1x,i4),(1x,i6),(2x,i4),(1x,i5),1p,2(2x,d\
10.3))";
    static char fmt_1004[] = "(/,a4,1p,6(1x,d11.4),/,(4x,1p,6(1x,d11.4)))";
    static char fmt_3009[] = "(/,a60)";
    static char fmt_9011[] = "(/,\002 Matrix in 1st Cholesky factorization i\
n formk is not Pos. Def.\002)";
    static char fmt_9012[] = "(/,\002 Matrix in 2st Cholesky factorization i\
n formk is not Pos. Def.\002)";
    static char fmt_9013[] = "(/,\002 Matrix in the Cholesky factorization i\
n formt is not Pos. Def.\002)";
    static char fmt_9014[] = "(/,\002 Derivative >= 0, backtracking line sea\
rch impossible.\002,/,\002   Previous x, f and g restored.\002,/,\002 Possib\
le causes: 1 error in function or gradient evaluation;\002,/,\002           \
       2 rounding errors dominate computation.\002)";
    static char fmt_9015[] = "(/,\002 Warning:  more than 10 function and gr\
adient\002,/,\002   evaluations in the last line search.  Termination\002,/\
,\002   may possibly be caused by a bad search direction.\002)";
    static char fmt_9018[] = "(/,\002 The triangular system is singular.\002)"
	    ;
    static char fmt_9019[] = "(/,\002 Line search cannot locate an adequate \
point after 20 function\002,/,\002  and gradient evaluations.  Previous x, f\
 and g restored.\002,/,\002 Possible causes: 1 error in function or gradient\
 evaluation;\002,/,\002                  2 rounding error dominate computati\
on.\002)";
    static char fmt_3007[] = "(/,\002 Cauchy                time\002,1p,e10.\
3,\002 seconds.\002,/\002 Subspace minimization time\002,1p,e10.3,\002 secon\
ds.\002,/\002 Line search           time\002,1p,e10.3,\002 seconds.\002)";
    static char fmt_3008[] = "(/,\002 Total User time\002,1p,e10.3,\002 seco\
nds.\002,/)";
    static char fmt_3002[] = "(2(1x,i4),2(1x,i5),2x,a3,1x,i4,1p,2(2x,d7.1),6\
    x,\002-\002,10x,\002-\002)";*/

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    //    integer s_cmp(), s_wsfe(), e_wsfe(), do_fio(), s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    static integer i__;

    /* Fortran I/O blocks */
    /*    static cilist io___202 = { 0, 6, 0, fmt_3003, 0 };
    static cilist io___203 = { 0, 6, 0, fmt_3004, 0 };
    static cilist io___204 = { 0, 6, 0, fmt_3005, 0 };
    static cilist io___205 = { 0, 6, 0, fmt_1004, 0 };
    static cilist io___207 = { 0, 6, 0, 0, 0 };
    static cilist io___208 = { 0, 6, 0, fmt_3009, 0 };
    static cilist io___209 = { 0, 6, 0, fmt_9011, 0 };
    static cilist io___210 = { 0, 6, 0, fmt_9012, 0 };
    static cilist io___211 = { 0, 6, 0, fmt_9013, 0 };
    static cilist io___212 = { 0, 6, 0, fmt_9014, 0 };
    static cilist io___213 = { 0, 6, 0, fmt_9015, 0 };
    static cilist io___214 = { 0, 6, 0, 0, 0 };
    static cilist io___215 = { 0, 6, 0, 0, 0 };
    static cilist io___216 = { 0, 6, 0, fmt_9018, 0 };
    static cilist io___217 = { 0, 6, 0, fmt_9019, 0 };
    static cilist io___218 = { 0, 6, 0, fmt_3007, 0 };
    static cilist io___219 = { 0, 6, 0, fmt_3008, 0 };
    static cilist io___220 = { 0, 0, 0, fmt_3002, 0 };
    static cilist io___221 = { 0, 0, 0, fmt_3009, 0 };
    static cilist io___222 = { 0, 0, 0, fmt_9011, 0 };
    static cilist io___223 = { 0, 0, 0, fmt_9012, 0 };
    static cilist io___224 = { 0, 0, 0, fmt_9013, 0 };
    static cilist io___225 = { 0, 0, 0, fmt_9014, 0 };
    static cilist io___226 = { 0, 0, 0, fmt_9015, 0 };
    static cilist io___227 = { 0, 0, 0, fmt_9018, 0 };
    static cilist io___228 = { 0, 0, 0, fmt_9019, 0 };
    static cilist io___229 = { 0, 0, 0, fmt_3008, 0 };*/


/*     ************ */

/*     Subroutine prn3lb */

/*     This subroutine prints out information when either a built-in */
/*       convergence test is satisfied or when an error message is */
/*       generated. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    if (strncmp(task, "ERROR", 5) == 0) {
	goto L999;
    }
    if (*iprint >= 0) {
      printf("\n\t* * *\n\nTit   = total number of iterations\nTnf   = total number of function evaluations\nTnint = total number of segments explored during Cauchy searches\nSkip  = number of BFGS updates skipped\nNact  = number of active bounds at final generalized Cuachy point\nProjg = norm of the final projected gradient\nF     = final function value\n\n\t* * *\n");
      /*	s_wsfe(&io___202);
		e_wsfe();*/
      printf("N\tTit\tTnf\tTnint\tSkip\tNact\tProjg\tF\n");
      printf("%i\t%i\t%i\t%i\t%i\t%i\t%e\t%e\n",
	     *n, *iter, *nfgv, *nintol, *nskip, *nact, *sbgnrm, *f);
      /*	s_wsfe(&io___203);
	e_wsfe();
	s_wsfe(&io___204);
	do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*iter), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*nfgv), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*nintol), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*nskip), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*nact), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*sbgnrm), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*f), (ftnlen)sizeof(doublereal));
	e_wsfe();*/
	if (*iprint >= 100) {
	  printf("X = ");
	  i__1 = *n;
	  for (i__ = 1; i__ <= i__1; ++i__) 
	    printf("%e", x[i__]);
	  printf("\n");
	  /*	    s_wsfe(&io___205);
	    do_fio(&c__1, "X =", (ftnlen)3);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	    }
	    e_wsfe();*/
	}
	if (*iprint >= 1) {
	  printf(" F = %f\n", *f);
	  /*	    s_wsle(&io___207);
	    do_lio(&c__9, &c__1, " F =", (ftnlen)4);
	    do_lio(&c__5, &c__1, (char *)&(*f), (ftnlen)sizeof(doublereal));
	    e_wsle();*/
	}
    }
L999:
    if (*iprint >= 0) {
      printf("%s\n", task);
      /*	s_wsfe(&io___208);
	do_fio(&c__1, task, (ftnlen)60);
	e_wsfe();*/
	if (*info != 0) {
	  if (*info == -1) {
	    printf("Matrix in 1st Cholesky factorization in formk is not Pos. Def.\n");
	    /*	    
		s_wsfe(&io___209);
		e_wsfe();*/
	    }
	    if (*info == -2) {
	      printf(" Matrix in 2st Cholesky factorization in formk is not Pos. Def.\n");
	      /*		s_wsfe(&io___210);
				e_wsfe();*/
	    }
	    if (*info == -3) {
	      printf("Matrix in the Cholesky factorization in formt is not Pos. Def.\n");
	      /*		s_wsfe(&io___211);
				e_wsfe();*/
	    }
	    if (*info == -4) {
	      printf("Derivative >= 0, backtracking line search impossible.\n, Previous x, f and g restored.\nPossible causes: 1 error in function or gradient evaluation;\n2 rounding errors dominate computation.\n");
	      /*		s_wsfe(&io___212);
				e_wsfe();*/
	    }
	    if (*info == -5) {
	      printf("Warning:  more than 10 function and gradient evaluations in the last line search.  Termination may possibly be caused by a bad search direction.\n");
	      /*		s_wsfe(&io___213);
				e_wsfe();*/
	    }
	    if (*info == -6) {
	      printf(" Input nbd(%i) is invalid\n", *k);
	      /*		s_wsle(&io___214);
		do_lio(&c__9, &c__1, " Input nbd(", (ftnlen)11);
		do_lio(&c__3, &c__1, (char *)&(*k), (ftnlen)sizeof(integer));
		do_lio(&c__9, &c__1, ") is invalid.", (ftnlen)13);
		e_wsle();*/
	    }
	    if (*info == -7) {
	      printf(" l(%i) > u(%i). No feasible solutions.\n",
		     *k, *k);
	      /*
		s_wsle(&io___215);
		do_lio(&c__9, &c__1, " l(", (ftnlen)3);
		do_lio(&c__3, &c__1, (char *)&(*k), (ftnlen)sizeof(integer));
		do_lio(&c__9, &c__1, ") > u(", (ftnlen)6);
		do_lio(&c__3, &c__1, (char *)&(*k), (ftnlen)sizeof(integer));
		do_lio(&c__9, &c__1, ").  No feasible solution.", (ftnlen)25);
		e_wsle();*/
	    }
	    if (*info == -8) {
	      printf("The triangular system is singular.\n");
	      /*		s_wsfe(&io___216);
				e_wsfe();*/
	    }
	    if (*info == -9) {
	      printf("Line search cannot locate an adequate point after 20 function and gradient evaluations.  Previous x, f and g restored.  Possible causes: 1 error in function or gradient evaluation; 2 rounding error dominate computation.\n");
	      /*		s_wsfe(&io___217);
				e_wsfe();*/
	    }
	}
	if (*iprint >= 1) {
	  printf("Cauchy time %.0f seconds.\nSubspace minimization time %.0f seconds.\nLine search time %.0f seconds.\n", *cachyt, *sbtime, *lnscht);
	  /*	    s_wsfe(&io___218);
	    do_fio(&c__1, (char *)&(*cachyt), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*sbtime), (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&(*lnscht), (ftnlen)sizeof(doublereal));
	    e_wsfe();*/
	}
	printf("Total User time %.0f seconds.\n", *time);
	/*	s_wsfe(&io___219);
	do_fio(&c__1, (char *)&(*time), (ftnlen)sizeof(doublereal));
	e_wsfe();*/
	/*	if (*iprint >= 1) {
		if (*info == -4 || *info == -9) {
		io___220.ciunit = *itfile;
		s_wsfe(&io___220);
		do_fio(&c__1, (char *)&(*iter), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*nfgv), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*nint), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*nact), (ftnlen)sizeof(integer));
		do_fio(&c__1, word, (ftnlen)3);
		do_fio(&c__1, (char *)&(*iback), (ftnlen)sizeof(integer));
		do_fio(&c__1, (char *)&(*stp), (ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&(*xstep), (ftnlen)sizeof(doublereal));
		e_wsfe();
		}
		io___221.ciunit = *itfile;
	    s_wsfe(&io___221);
	    do_fio(&c__1, task, (ftnlen)60);
	    e_wsfe();
	    if (*info != 0) {
		if (*info == -1) {
		    io___222.ciunit = *itfile;
		    s_wsfe(&io___222);
		    e_wsfe();
		}
		if (*info == -2) {
		    io___223.ciunit = *itfile;
		    s_wsfe(&io___223);
		    e_wsfe();
		}
		if (*info == -3) {
		    io___224.ciunit = *itfile;
		    s_wsfe(&io___224);
		    e_wsfe();
		}
		if (*info == -4) {
		    io___225.ciunit = *itfile;
		    s_wsfe(&io___225);
		    e_wsfe();
		}
		if (*info == -5) {
		    io___226.ciunit = *itfile;
		    s_wsfe(&io___226);
		    e_wsfe();
		}
		if (*info == -8) {
		    io___227.ciunit = *itfile;
		    s_wsfe(&io___227);
		    e_wsfe();
		}
		if (*info == -9) {
		    io___228.ciunit = *itfile;
		    s_wsfe(&io___228);
		    e_wsfe();
		}
		}
	    io___229.ciunit = *itfile;
	    s_wsfe(&io___229);
	    do_fio(&c__1, (char *)&(*time), (ftnlen)sizeof(doublereal));
	    e_wsfe();
	    }*/
    }
/* L3006: */
    return 0;
} /* prn3lb_ */

/* ======================= The end of prn3lb ============================= */
/* Subroutine */ int projgr_(n, l, u, nbd, x, g, sbgnrm)
integer *n;
doublereal *l, *u;
integer *nbd;
doublereal *x, *g, *sbgnrm;
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__;
    static doublereal gi;

/*     ************ */

/*     Subroutine projgr */

/*     This subroutine computes the infinity norm of the projected */
/*       gradient. */


/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --g;
    --x;
    --nbd;
    --u;
    --l;

    /* Function Body */
    *sbgnrm = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gi = g[i__];
	if (nbd[i__] != 0) {
	    if (gi < 0.) {
		if (nbd[i__] >= 2) {
/* Computing MAX */
		    d__1 = x[i__] - u[i__];
		    gi = max(d__1,gi);
		}
	    } else {
		if (nbd[i__] <= 2) {
/* Computing MIN */
		    d__1 = x[i__] - l[i__];
		    gi = min(d__1,gi);
		}
	    }
	}
/* Computing MAX */
	d__1 = *sbgnrm, d__2 = abs(gi);
	*sbgnrm = max(d__1,d__2);
/* L15: */
    }
    return 0;
} /* projgr_ */

/* ======================= The end of projgr ============================= */
/* Subroutine */ int subsm_(n, m, nsub, ind, l, u, nbd, x, d__, ws, wy, theta,
	 col, head, iword, wv, wn, iprint, info)
integer *n, *m, *nsub, *ind;
doublereal *l, *u;
integer *nbd;
doublereal *x, *d__, *ws, *wy, *theta;
integer *col, *head, *iword;
doublereal *wv, *wn;
integer *iprint, *info;
{
    /* Format strings */
  /*    static char fmt_1001[] = "(/,\002----------------SUBSM entered----------\
-------\002,/)";
    static char fmt_1002[] = "(\002ALPHA = \002,f7.5,\002 backtrack to the B\
OX\002)";
    static char fmt_1003[] = "(\002Subspace solution X =  \002,/,(4x,1p,6(1x\
,d11.4)))";
    static char fmt_1004[] = "(/,\002----------------exit SUBSM ------------\
    --------\002,/)";*/

    /* System generated locals */
    integer ws_dim1, ws_offset, wy_dim1, wy_offset, wn_dim1, wn_offset, i__1, 
	    i__2;

    /* Builtin functions */
    integer s_wsfe(), e_wsfe(), do_fio(), s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    static doublereal temp1, temp2;
    static integer i__, j, k;
    static doublereal alpha;
    extern /* Subroutine */ int dtrsl_();
    static integer m2;
    static doublereal dk;
    static integer js, jy, pointr, ibd, col2;

    /* Fortran I/O blocks */
    /*    static cilist io___232 = { 0, 6, 0, fmt_1001, 0 };
    static cilist io___246 = { 0, 6, 0, fmt_1002, 0 };
    static cilist io___247 = { 0, 6, 0, 0, 0 };
    static cilist io___248 = { 0, 6, 0, fmt_1003, 0 };
    static cilist io___249 = { 0, 6, 0, fmt_1004, 0 };*/


/*     ************ */

/*     Subroutine subsm */

/*     Given xcp, l, u, r, an index set that specifies */
/* 	the active set at xcp, and an l-BFGS matrix B */
/* 	(in terms of WY, WS, SY, WT, head, col, and theta), */
/* 	this subroutine computes an approximate solution */
/* 	of the subspace problem */

/*     	(P)   min Q(x) = r'(x-xcp) + 1/2 (x-xcp)' B (x-xcp) */

/*             subject to l<=x<=u */
/* 	  	        x_i=xcp_i for all i in A(xcp) */

/* 	along the subspace unrained Newton direction */

/* 	   d = -(Z'BZ)^(-1) r. */

/*       The formula for the Newton direction, given the L-BFGS matrix */
/*       and the Sherman-Morrison formula, is */

/* 	   d = (1/theta)r + (1/theta*2) Z'WK^(-1)W'Z r. */

/*       where */
/*                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ] */
/*                     [L_a -R_z           theta*S'AA'S ] */

/*     Note that this procedure for computing d differs */
/*     from that described in [1]. One can show that the matrix K is */
/*     equal to the matrix M^[-1]N in that paper. */

/*     n is an integer variable. */
/*       On entry n is the dimension of the problem. */
/*       On exit n is unchanged. */

/*     m is an integer variable. */
/*       On entry m is the maximum number of variable metric corrections */
/*         used to define the limited memory matrix. */
/*       On exit m is unchanged. */

/*     nsub is an integer variable. */
/*       On entry nsub is the number of free variables. */
/*       On exit nsub is unchanged. */

/*     ind is an integer array of dimension nsub. */
/*       On entry ind specifies the coordinate indices of free variables. */
/*       On exit ind is unchanged. */

/*     l is a double precision array of dimension n. */
/*       On entry l is the lower bound of x. */
/*       On exit l is unchanged. */

/*     u is a double precision array of dimension n. */
/*       On entry u is the upper bound of x. */
/*       On exit u is unchanged. */

/*     nbd is a integer array of dimension n. */
/*       On entry nbd represents the type of bounds imposed on the */
/*         variables, and must be specified as follows: */
/*         nbd(i)=0 if x(i) is unbounded, */
/*                1 if x(i) has only a lower bound, */
/*                2 if x(i) has both lower and upper bounds, and */
/*                3 if x(i) has only an upper bound. */
/*       On exit nbd is unchanged. */

/*     x is a double precision array of dimension n. */
/*       On entry x specifies the Cauchy point xcp. */
/*       On exit x(i) is the minimizer of Q over the subspace of */
/*                                                        free variables. */

/*     d is a double precision array of dimension n. */
/*       On entry d is the reduced gradient of Q at xcp. */
/*       On exit d is the Newton direction of Q. */

/*     ws and wy are double precision arrays; */
/*     theta is a double precision variable; */
/*     col is an integer variable; */
/*     head is an integer variable. */
/*       On entry they store the information defining the */
/*                                          limited memory BFGS matrix: */
/*         ws(n,m) stores S, a set of s-vectors; */
/*         wy(n,m) stores Y, a set of y-vectors; */
/*         theta is the scaling factor specifying B_0 = theta I; */
/*         col is the number of variable metric corrections stored; */
/*         head is the location of the 1st s- (or y-) vector in S (or Y). */
/*       On exit they are unchanged. */

/*     iword is an integer variable. */
/*       On entry iword is unspecified. */
/*       On exit iword specifies the status of the subspace solution. */
/*         iword = 0 if the solution is in the box, */
/*                 1 if some bound is encountered. */

/*     wv is a double precision working array of dimension 2m. */

/*     wn is a double precision array of dimension 2m x 2m. */
/*       On entry the upper triangle of wn stores the LEL^T factorization */
/*         of the indefinite matrix */

/*              K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ] */
/*                  [L_a -R_z           theta*S'AA'S ] */
/*                                                    where E = [-I  0] */
/*                                                              [ 0  I] */
/*       On exit wn is unchanged. */

/*     iprint is an INTEGER variable that must be set by the user. */
/*       It controls the frequency and type of output generated: */
/*        iprint<0    no output is generated; */
/*        iprint=0    print only one line at the last iteration; */
/*        0<iprint<99 print also f and |proj g| every iprint iterations; */
/*        iprint=99   print details of every iteration except n-vectors; */
/*        iprint=100  print also the changes of active set and final x; */
/*        iprint>100  print details of every iteration including x and g; */
/*       When iprint > 0, the file iterate.dat will be created to */
/*                        summarize the iteration. */

/*     info is an integer variable. */
/*       On entry info is unspecified. */
/*       On exit info = 0       for normal return, */
/*                    = nonzero for abnormal return */
/*                                  when the matrix K is ill-conditioned. */

/*     Subprograms called: */

/*       Linpack dtrsl. */


/*     References: */

/*       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited */
/*       memory algorithm for bound rained optimization'', */
/*       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208. */



/*                           *  *  * */

/*     NEOS, November 1994. (Latest revision June 1996.) */
/*     Optimization Technology Center. */
/*     Argonne National Laboratory and Northwestern University. */
/*     Written by */
/*                        Ciyou Zhu */
/*     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal. */


/*     ************ */
    /* Parameter adjustments */
    --d__;
    --x;
    --nbd;
    --u;
    --l;
    wn_dim1 = 2 * *m;
    wn_offset = 1 + wn_dim1 * 1;
    wn -= wn_offset;
    --wv;
    wy_dim1 = *n;
    wy_offset = 1 + wy_dim1 * 1;
    wy -= wy_offset;
    ws_dim1 = *n;
    ws_offset = 1 + ws_dim1 * 1;
    ws -= ws_offset;
    --ind;

    /* Function Body */
    if (*nsub <= 0) {
	return 0;
    }
    if (*iprint >= 99) {
      printf("----------------SUBSM entered-----------------\n");
      /*	s_wsfe(&io___232);
		e_wsfe();*/
    }
/*     Compute wv = W'Zd. */
    pointr = *head;
    i__1 = *col;
    for (i__ = 1; i__ <= i__1; ++i__) {
	temp1 = 0.;
	temp2 = 0.;
	i__2 = *nsub;
	for (j = 1; j <= i__2; ++j) {
	    k = ind[j];
	    temp1 += wy[k + pointr * wy_dim1] * d__[j];
	    temp2 += ws[k + pointr * ws_dim1] * d__[j];
/* L10: */
	}
	wv[i__] = temp1;
	wv[*col + i__] = *theta * temp2;
	pointr = pointr % *m + 1;
/* L20: */
    }
/*     Compute wv:=K^(-1)wv. */
    m2 = *m << 1;
    col2 = *col << 1;
    dtrsl_(&wn[wn_offset], &m2, &col2, &wv[1], &c__11, info);
    if (*info != 0) {
	return 0;
    }
    i__1 = *col;
    for (i__ = 1; i__ <= i__1; ++i__) {
	wv[i__] = -wv[i__];
/* L25: */
    }
    dtrsl_(&wn[wn_offset], &m2, &col2, &wv[1], &c__1, info);
    if (*info != 0) {
	return 0;
    }
/*     Compute d = (1/theta)d + (1/theta**2)Z'W wv. */
    pointr = *head;
    i__1 = *col;
    for (jy = 1; jy <= i__1; ++jy) {
	js = *col + jy;
	i__2 = *nsub;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    k = ind[i__];
	    d__[i__] = d__[i__] + wy[k + pointr * wy_dim1] * wv[jy] / *theta 
		    + ws[k + pointr * ws_dim1] * wv[js];
/* L30: */
	}
	pointr = pointr % *m + 1;
/* L40: */
    }
    i__1 = *nsub;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__[i__] /= *theta;
/* L50: */
    }
/*     Backtrack to the feasible region. */
    alpha = 1.;
    temp1 = alpha;
    i__1 = *nsub;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = ind[i__];
	dk = d__[i__];
	if (nbd[k] != 0) {
	    if (dk < 0. && nbd[k] <= 2) {
		temp2 = l[k] - x[k];
		if (temp2 >= 0.) {
		    temp1 = 0.;
		} else if (dk * alpha < temp2) {
		    temp1 = temp2 / dk;
		}
	    } else if (dk > 0. && nbd[k] >= 2) {
		temp2 = u[k] - x[k];
		if (temp2 <= 0.) {
		    temp1 = 0.;
		} else if (dk * alpha > temp2) {
		    temp1 = temp2 / dk;
		}
	    }
	    if (temp1 < alpha) {
		alpha = temp1;
		ibd = i__;
	    }
	}
/* L60: */
    }
    if (alpha < 1.) {
	dk = d__[ibd];
	k = ind[ibd];
	if (dk > 0.) {
	    x[k] = u[k];
	    d__[ibd] = 0.;
	} else if (dk < 0.) {
	    x[k] = l[k];
	    d__[ibd] = 0.;
	}
    }
    i__1 = *nsub;
    for (i__ = 1; i__ <= i__1; ++i__) {
	k = ind[i__];
	x[k] += alpha * d__[i__];
/* L70: */
    }
    if (*iprint >= 99) {
      if (alpha < 1.) {
	printf("ALPHA = %e backtrack to the BOX\n", alpha);
	    /*
	    s_wsfe(&io___246);
	    do_fio(&c__1, (char *)&alpha, (ftnlen)sizeof(doublereal));
	    e_wsfe();*/
	} else {
	  printf("SM solution inside the box\n");
	  /*	    s_wsle(&io___247);
	    do_lio(&c__9, &c__1, "SM solution inside the box", (ftnlen)26);
	    e_wsle();*/
	}
	if (*iprint > 100) {
	  i__1 = *n;
	  for (i__ = 1; i__ <= i__1; ++i__) 
	    printf("\t%e", x[i__]);
	  printf("\n");
	  /*
	    s_wsfe(&io___248);
	    i__1 = *n;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		do_fio(&c__1, (char *)&x[i__], (ftnlen)sizeof(doublereal));
	    }
	    e_wsfe();*/
	}
    }
    if (alpha < 1.) {
	*iword = 1;
    } else {
	*iword = 0;
    }
    if (*iprint >= 99) {
      printf("----------------exit SUBSM --------------------\n");
      /*	s_wsfe(&io___249);
		e_wsfe();*/
    }
    return 0;
} /* subsm_ */

/* ====================== The end of subsm =============================== */
/* Subroutine */ int dcsrch_(f, g, stp, ftol, gtol, xtol, stpmin, stpmax, 
	task, isave, dsave, task_len)
doublereal *f, *g, *stp, *ftol, *gtol, *xtol, *stpmin, *stpmax;
char *task;
integer *isave;
doublereal *dsave;
ftnlen task_len;
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    //    integer s_cmp();
    //    /* Subroutine */ int s_copy();

    /* Local variables */
    static integer stage;
    static doublereal finit, ginit, width, ftest, gtest, stmin, stmax, width1,
	     fm, gm, fx, fy, gx, gy;
    static logical brackt;
    extern /* Subroutine */ int dcstep_();
    static doublereal fxm, fym, gxm, gym, stx, sty;

/*     ********** */

/*     Subroutine dcsrch */

/*     This subroutine finds a step that satisfies a sufficient */
/*     decrease condition and a curvature condition. */

/*     Each call of the subroutine updates an interval with */
/*     endpoints stx and sty. The interval is initially chosen */
/*     so that it contains a minimizer of the modified function */

/*           psi(stp) = f(stp) - f(0) - ftol*stp*f'(0). */

/*     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the */
/*     interval is chosen so that it contains a minimizer of f. */

/*     The algorithm is designed to find a step that satisfies */
/*     the sufficient decrease condition */

/*           f(stp) <= f(0) + ftol*stp*f'(0), */

/*     and the curvature condition */

/*           abs(f'(stp)) <= gtol*abs(f'(0)). */

/*     If ftol is less than gtol and if, for example, the function */
/*     is bounded below, then there is always a step which satisfies */
/*     both conditions. */

/*     If no step can be found that satisfies both conditions, then */
/*     the algorithm stops with a warning. In this case stp only */
/*     satisfies the sufficient decrease condition. */

/*     A typical invocation of dcsrch has the following outline: */

/*     task = 'START' */
/*  10 continue */
/*        call dcsrch( ... ) */
/*        if (task .eq. 'FG') then */
/*           Evaluate the function and the gradient at stp */
/*           goto 10 */
/*           end if */

/*     NOTE: The user must no alter work arrays between calls. */

/*     The subroutine statement is */

/*        subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax, */
/*                          task,isave,dsave) */
/*     where */

/*       f is a double precision variable. */
/*         On initial entry f is the value of the function at 0. */
/*            On subsequent entries f is the value of the */
/*            function at stp. */
/*         On exit f is the value of the function at stp. */

/* 	g is a double precision variable. */
/*         On initial entry g is the derivative of the function at 0. */
/*            On subsequent entries g is the derivative of the */
/*            function at stp. */
/*         On exit g is the derivative of the function at stp. */

/* 	stp is a double precision variable. */
/*         On entry stp is the current estimate of a satisfactory */
/*            step. On initial entry, a positive initial estimate */
/*            must be provided. */
/*         On exit stp is the current estimate of a satisfactory step */
/*            if task = 'FG'. If task = 'CONV' then stp satisfies */
/*            the sufficient decrease and curvature condition. */

/*       ftol is a double precision variable. */
/*         On entry ftol specifies a nonnegative tolerance for the */
/*            sufficient decrease condition. */
/*         On exit ftol is unchanged. */

/*       gtol is a double precision variable. */
/*         On entry gtol specifies a nonnegative tolerance for the */
/*            curvature condition. */
/*         On exit gtol is unchanged. */

/* 	xtol is a double precision variable. */
/*         On entry xtol specifies a nonnegative relative tolerance */
/*            for an acceptable step. The subroutine exits with a */
/*            warning if the relative difference between sty and stx */
/*            is less than xtol. */
/*         On exit xtol is unchanged. */

/* 	stpmin is a double precision variable. */
/*         On entry stpmin is a nonnegative lower bound for the step. */
/*         On exit stpmin is unchanged. */

/* 	stpmax is a double precision variable. */
/*         On entry stpmax is a nonnegative upper bound for the step. */
/*         On exit stpmax is unchanged. */

/*       task is a character variable of length at least 60. */
/*         On initial entry task must be set to 'START'. */
/*         On exit task indicates the required action: */

/*            If task(1:2) = 'FG' then evaluate the function and */
/*            derivative at stp and call dcsrch again. */

/*            If task(1:4) = 'CONV' then the search is successful. */

/*            If task(1:4) = 'WARN' then the subroutine is not able */
/*            to satisfy the convergence conditions. The exit value of */
/*            stp contains the best point found during the search. */

/*            If task(1:5) = 'ERROR' then there is an error in the */
/*            input arguments. */

/*         On exit with convergence, a warning or an error, the */
/*            variable task contains additional information. */

/*       isave is an integer work array of dimension 2. */

/*       dsave is a double precision work array of dimension 13. */

/*     Subprograms called */

/* 	MINPACK-2 ... dcstep */

/*     MINPACK-1 Project. June 1983. */
/*     Argonne National Laboratory. */
/*     Jorge J. More' and David J. Thuente. */

/*     MINPACK-2 Project. October 1993. */
/*     Argonne National Laboratory and University of Minnesota. */
/*     Brett M. Averick, Richard G. Carter, and Jorge J. More'. */

/*     ********** */
/*     Initialization block. */
    /* Parameter adjustments */
    --dsave;
    --isave;

    /* Function Body */
    if (strncmp(task, "START", 5) == 0) {
/*        Check the input arguments for errors. */
	if (*stp < *stpmin) {
	    my_strcpy(task, "ERROR: STP .LT. STPMIN");
	}
	if (*stp > *stpmax) {
	    my_strcpy(task, "ERROR: STP .GT. STPMAX");
	}
	if (*g >= 0.) {
	    my_strcpy(task, "ERROR: INITIAL G .GE. ZERO");
	}
	if (*ftol < 0.) {
	    my_strcpy(task, "ERROR: FTOL .LT. ZERO");
	}
	if (*gtol < 0.) {
	    my_strcpy(task, "ERROR: GTOL .LT. ZERO");
	}
	if (*xtol < 0.) {
	  my_strcpy(task, "ERROR: XTOL .LT. ZERO");
	}
	if (*stpmin < 0.) {
	  my_strcpy(task, "ERROR: STPMIN .LT. ZERO");
	}
	if (*stpmax < *stpmin) {
	  my_strcpy(task, "ERROR: STPMAX .LT. STPMIN");
	}
	/*        Exit if there are errors on input. */
	if (strncmp(task, "ERROR", 5) == 0) {
	  return 0;
	}
	/*        Initialize local variables. */
	brackt = FALSE_;
	stage = 1;
	finit = *f;
	ginit = *g;
	gtest = *ftol * ginit;
	width = *stpmax - *stpmin;
	width1 = width / .5;
/*        The variables stx, fx, gx contain the values of the step, */
/*        function, and derivative at the best step. */
/*        The variables sty, fy, gy contain the value of the step, */
/*        function, and derivative at sty. */
/*        The variables stp, f, g contain the values of the step, */
/*        function, and derivative at stp. */
	stx = 0.;
	fx = finit;
	gx = ginit;
	sty = 0.;
	fy = finit;
	gy = ginit;
	stmin = 0.;
	stmax = *stp + *stp * 4.;
	task[0]='F';
	task[1]='G';
	goto L1000;
    } else {
/*        Restore local variables. */
	if (isave[1] == 1) {
	    brackt = TRUE_;
	} else {
	    brackt = FALSE_;
	}
	stage = isave[2];
	ginit = dsave[1];
	gtest = dsave[2];
	gx = dsave[3];
	gy = dsave[4];
	finit = dsave[5];
	fx = dsave[6];
	fy = dsave[7];
	stx = dsave[8];
	sty = dsave[9];
	stmin = dsave[10];
	stmax = dsave[11];
	width = dsave[12];
	width1 = dsave[13];
    }
/*     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the */
/*     algorithm enters the second stage. */
    ftest = finit + *stp * gtest;
    if (stage == 1 && *f <= ftest && *g >= 0.) {
	stage = 2;
    }
/*     Test for warnings. */
    if (brackt && (*stp <= stmin || *stp >= stmax)) {
      my_strcpy(task, "WARNING: ROUNDING ERRORS PREVENT PROGRESS");
    }
    if (brackt && stmax - stmin <= *xtol * stmax) {
	my_strcpy(task, "WARNING: XTOL TEST SATISFIED");
    }
    if (*stp == *stpmax && *f <= ftest && *g <= gtest) {
	my_strcpy(task, "WARNING: STP = STPMAX");
    }
    if (*stp == *stpmin && (*f > ftest || *g >= gtest)) {
	my_strcpy(task, "WARNING: STP = STPMIN");
    }
/*     Test for convergence. */
    if (*f <= ftest && abs(*g) <= *gtol * (-ginit)) {
	my_strcpy(task, "CONVERGENCE");
    }
/*     Test for termination. */
    if (strncmp(task, "WARN", 4) == 0 || strncmp(task, "CONV", 4)==0) {
	goto L1000;
    }
/*     A modified function is used to predict the step during the */
/*     first stage if a lower function value has been obtained but */
/*     the decrease is not sufficient. */
    if (stage == 1 && *f <= fx && *f > ftest) {
/*        Define the modified function and derivative values. */
	fm = *f - *stp * gtest;
	fxm = fx - stx * gtest;
	fym = fy - sty * gtest;
	gm = *g - gtest;
	gxm = gx - gtest;
	gym = gy - gtest;
/*        Call dcstep to update stx, sty, and to compute the new step. */
	dcstep_(&stx, &fxm, &gxm, &sty, &fym, &gym, stp, &fm, &gm, &brackt, &
		stmin, &stmax);
/*        Reset the function and derivative values for f. */
	fx = fxm + stx * gtest;
	fy = fym + sty * gtest;
	gx = gxm + gtest;
	gy = gym + gtest;
    } else {
/*       Call dcstep to update stx, sty, and to compute the new step. */
	dcstep_(&stx, &fx, &gx, &sty, &fy, &gy, stp, f, g, &brackt, &stmin, &
		stmax);
    }
/*     Decide if a bisection step is needed. */
    if (brackt) {
	if ((d__1 = sty - stx, abs(d__1)) >= width1 * .66) {
	    *stp = stx + (sty - stx) * .5;
	}
	width1 = width;
	width = (d__1 = sty - stx, abs(d__1));
    }
/*     Set the minimum and maximum steps allowed for stp. */
    if (brackt) {
	stmin = min(stx,sty);
	stmax = max(stx,sty);
    } else {
	stmin = *stp + (*stp - stx) * 1.1;
	stmax = *stp + (*stp - stx) * 4.;
    }
/*     Force the step to be within the bounds stpmax and stpmin. */
    *stp = max(*stp,*stpmin);
    *stp = min(*stp,*stpmax);
/*     If further progress is not possible, let stp be the best */
/*     point obtained during the search. */
    if ((brackt && (*stp <= stmin || *stp >= stmax)) || (brackt && stmax - stmin 
	    <= *xtol * stmax)) {
	*stp = stx;
    }
/*     Obtain another function and derivative. */
    task[0]='F';
    task[1]='G';
L1000:
/*     Save local variables. */
    if (brackt) {
	isave[1] = 1;
    } else {
	isave[1] = 0;
    }
    isave[2] = stage;
    dsave[1] = ginit;
    dsave[2] = gtest;
    dsave[3] = gx;
    dsave[4] = gy;
    dsave[5] = finit;
    dsave[6] = fx;
    dsave[7] = fy;
    dsave[8] = stx;
    dsave[9] = sty;
    dsave[10] = stmin;
    dsave[11] = stmax;
    dsave[12] = width;
    dsave[13] = width1;
    return 0;
} /* dcsrch_ */

/* ====================== The end of dcsrch ============================== */
/* Subroutine */ int dcstep_(stx, fx, dx, sty, fy, dy, stp, fp, dp, brackt, 
	stpmin, stpmax)
doublereal *stx, *fx, *dx, *sty, *fy, *dy, *stp, *fp, *dp;
logical *brackt;
doublereal *stpmin, *stpmax;
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static doublereal sgnd, stpc, stpf, stpq, p, q, gamma, r__, s, theta;

/*     ********** */

/*     Subroutine dcstep */

/*     This subroutine computes a safeguarded step for a search */
/*     procedure and updates an interval that contains a step that */
/*     satisfies a sufficient decrease and a curvature condition. */

/*     The parameter stx contains the step with the least function */
/*     value. If brackt is set to .true. then a minimizer has */
/*     been bracketed in an interval with endpoints stx and sty. */
/*     The parameter stp contains the current step. */
/*     The subroutine assumes that if brackt is set to .true. then */

/*           min(stx,sty) < stp < max(stx,sty), */

/*     and that the derivative at stx is negative in the direction */
/*     of the step. */

/*     The subroutine statement is */

/*       subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt, */
/*                         stpmin,stpmax) */

/*     where */

/*       stx is a double precision variable. */
/*         On entry stx is the best step obtained so far and is an */
/*            endpoint of the interval that contains the minimizer. */
/*         On exit stx is the updated best step. */

/*       fx is a double precision variable. */
/*         On entry fx is the function at stx. */
/*         On exit fx is the function at stx. */

/*       dx is a double precision variable. */
/*         On entry dx is the derivative of the function at */
/*            stx. The derivative must be negative in the direction of */
/*            the step, that is, dx and stp - stx must have opposite */
/*            signs. */
/*         On exit dx is the derivative of the function at stx. */

/*       sty is a double precision variable. */
/*         On entry sty is the second endpoint of the interval that */
/*            contains the minimizer. */
/*         On exit sty is the updated endpoint of the interval that */
/*            contains the minimizer. */

/*       fy is a double precision variable. */
/*         On entry fy is the function at sty. */
/*         On exit fy is the function at sty. */

/*       dy is a double precision variable. */
/*         On entry dy is the derivative of the function at sty. */
/*         On exit dy is the derivative of the function at the exit sty. */

/*       stp is a double precision variable. */
/*         On entry stp is the current step. If brackt is set to .true. */
/*            then on input stp must be between stx and sty. */
/*         On exit stp is a new trial step. */

/*       fp is a double precision variable. */
/*         On entry fp is the function at stp */
/*         On exit fp is unchanged. */

/*       dp is a double precision variable. */
/*         On entry dp is the the derivative of the function at stp. */
/*         On exit dp is unchanged. */

/*       brackt is an logical variable. */
/*         On entry brackt specifies if a minimizer has been bracketed. */
/*            Initially brackt must be set to .false. */
/*         On exit brackt specifies if a minimizer has been bracketed. */
/*            When a minimizer is bracketed brackt is set to .true. */

/*       stpmin is a double precision variable. */
/*         On entry stpmin is a lower bound for the step. */
/*         On exit stpmin is unchanged. */

/*       stpmax is a double precision variable. */
/*         On entry stpmax is an upper bound for the step. */
/*         On exit stpmax is unchanged. */

/*     MINPACK-1 Project. June 1983 */
/*     Argonne National Laboratory. */
/*     Jorge J. More' and David J. Thuente. */

/*     MINPACK-2 Project. October 1993. */
/*     Argonne National Laboratory and University of Minnesota. */
/*     Brett M. Averick and Jorge J. More'. */

/*     ********** */
    sgnd = *dp * (*dx / abs(*dx));
/*     First case: A higher function value. The minimum is bracketed. */
/*     If the cubic step is closer to stx than the quadratic step, the */
/*     cubic step is taken, otherwise the average of the cubic and */
/*     quadratic steps is taken. */
    if (*fp > *fx) {
	theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
/* Computing MAX */
	d__1 = abs(theta), d__2 = abs(*dx), d__1 = max(d__1,d__2), d__2 = abs(
		*dp);
	s = max(d__1,d__2);
/* Computing 2nd power */
	d__1 = theta / s;
	gamma = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
	if (*stp < *stx) {
	    gamma = -gamma;
	}
	p = gamma - *dx + theta;
	q = gamma - *dx + gamma + *dp;
	r__ = p / q;
	stpc = *stx + r__ * (*stp - *stx);
	stpq = *stx + *dx / ((*fx - *fp) / (*stp - *stx) + *dx) / 2. * (*stp 
		- *stx);
	if ((d__1 = stpc - *stx, abs(d__1)) < (d__2 = stpq - *stx, abs(d__2)))
		 {
	    stpf = stpc;
	} else {
	    stpf = stpc + (stpq - stpc) / 2.;
	}
	*brackt = TRUE_;
/*     Second case: A lower function value and derivatives of opposite */
/*     sign. The minimum is bracketed. If the cubic step is farther from */
/*     stp than the secant step, the cubic step is taken, otherwise the */
/*     secant step is taken. */
    } else if (sgnd < 0.) {
	theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
/* Computing MAX */
	d__1 = abs(theta), d__2 = abs(*dx), d__1 = max(d__1,d__2), d__2 = abs(
		*dp);
	s = max(d__1,d__2);
/* Computing 2nd power */
	d__1 = theta / s;
	gamma = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
	if (*stp > *stx) {
	    gamma = -gamma;
	}
	p = gamma - *dp + theta;
	q = gamma - *dp + gamma + *dx;
	r__ = p / q;
	stpc = *stp + r__ * (*stx - *stp);
	stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
	if ((d__1 = stpc - *stp, abs(d__1)) > (d__2 = stpq - *stp, abs(d__2)))
		 {
	    stpf = stpc;
	} else {
	    stpf = stpq;
	}
	*brackt = TRUE_;
/*     Third case: A lower function value, derivatives of the same sign, */
/*     and the magnitude of the derivative decreases. */
    } else if (abs(*dp) < abs(*dx)) {
/*        The cubic step is computed only if the cubic tends to infinity */
/*        in the direction of the step or if the minimum of the cubic */
/*        is beyond stp. Otherwise the cubic step is defined to be the */
/*        secant step. */
	theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
/* Computing MAX */
	d__1 = abs(theta), d__2 = abs(*dx), d__1 = max(d__1,d__2), d__2 = abs(
		*dp);
	s = max(d__1,d__2);
/*        The case gamma = 0 only arises if the cubic does not tend */
/*        to infinity in the direction of the step. */
/* Computing MAX */
/* Computing 2nd power */
	d__3 = theta / s;
	d__1 = 0., d__2 = d__3 * d__3 - *dx / s * (*dp / s);
	gamma = s * sqrt((max(d__1,d__2)));
	if (*stp > *stx) {
	    gamma = -gamma;
	}
	p = gamma - *dp + theta;
	q = gamma + (*dx - *dp) + gamma;
	r__ = p / q;
	if (r__ < 0. && gamma != 0.) {
	    stpc = *stp + r__ * (*stx - *stp);
	} else if (*stp > *stx) {
	    stpc = *stpmax;
	} else {
	    stpc = *stpmin;
	}
	stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
	if (*brackt) {
/*           A minimizer has been bracketed. If the cubic step is */
/*           closer to stp than the secant step, the cubic step is */
/*           taken, otherwise the secant step is taken. */
	    if ((d__1 = stpc - *stp, abs(d__1)) < (d__2 = stpq - *stp, abs(
		    d__2))) {
		stpf = stpc;
	    } else {
		stpf = stpq;
	    }
	    if (*stp > *stx) {
/* Computing MIN */
		d__1 = *stp + (*sty - *stp) * .66;
		stpf = min(d__1,stpf);
	    } else {
/* Computing MAX */
		d__1 = *stp + (*sty - *stp) * .66;
		stpf = max(d__1,stpf);
	    }
	} else {
/*           A minimizer has not been bracketed. If the cubic step is */
/*           farther from stp than the secant step, the cubic step is */
/*           taken, otherwise the secant step is taken. */
	    if ((d__1 = stpc - *stp, abs(d__1)) > (d__2 = stpq - *stp, abs(
		    d__2))) {
		stpf = stpc;
	    } else {
		stpf = stpq;
	    }
	    stpf = min(*stpmax,stpf);
	    stpf = max(*stpmin,stpf);
	}
/*     Fourth case: A lower function value, derivatives of the same sign, */
/*     and the magnitude of the derivative does not decrease. If the */
/*     minimum is not bracketed, the step is either stpmin or stpmax, */
/*     otherwise the cubic step is taken. */
    } else {
	if (*brackt) {
	    theta = (*fp - *fy) * 3. / (*sty - *stp) + *dy + *dp;
/* Computing MAX */
	    d__1 = abs(theta), d__2 = abs(*dy), d__1 = max(d__1,d__2), d__2 = 
		    abs(*dp);
	    s = max(d__1,d__2);
/* Computing 2nd power */
	    d__1 = theta / s;
	    gamma = s * sqrt(d__1 * d__1 - *dy / s * (*dp / s));
	    if (*stp > *sty) {
		gamma = -gamma;
	    }
	    p = gamma - *dp + theta;
	    q = gamma - *dp + gamma + *dy;
	    r__ = p / q;
	    stpc = *stp + r__ * (*sty - *stp);
	    stpf = stpc;
	} else if (*stp > *stx) {
	    stpf = *stpmax;
	} else {
	    stpf = *stpmin;
	}
    }
/*     Update the interval which contains a minimizer. */
    if (*fp > *fx) {
	*sty = *stp;
	*fy = *fp;
	*dy = *dp;
    } else {
	if (sgnd < 0.) {
	    *sty = *stx;
	    *fy = *fx;
	    *dy = *dx;
	}
	*stx = *stp;
	*fx = *fp;
	*dx = *dp;
    }
/*     Compute the new step. */
    *stp = stpf;
    return 0;
} /* dcstep_ */


typedef double real;

/* ====================== The end of dcstep ============================== */
/* Subroutine */ int timer_(ttime)
doublereal *ttime;
{
  //    static real temp;
    extern doublereal etime_();
    //    static real tarray[2];
    time_t currtime;
/*     ********* */

/*     Subroutine timer */

/*     This subroutine is used to determine user time. In a typical */
/*     application, the user time for a code segment requires calls */
/*     to subroutine timer to determine the initial and final time. */

/*     The subroutine statement is */

/*       subroutine timer(ttime) */

/*     where */

/*       ttime is an output variable which specifies the user time. */

/*     Argonne National Laboratory and University of Minnesota. */
/*     MINPACK-2 Project. */

/*     Modified October 1990 by Brett M. Averick. */

/*     ********** */
/*     The first element of the array tarray specifies user time */
//    temp = time_(tarray);
//    *ttime = (doublereal) tarray[0];
    time(&currtime);
    *ttime = (double)currtime;
    return 0;
} /* timer_ */

/* ====================== The end of timer =============================== */
doublereal dnrm2_(n, x, incx)
integer *n;
doublereal *x;
integer *incx;
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer i__;
    static doublereal scale;

/*     ********** */

/*     Function dnrm2 */

/*     Given a vector x of length n, this function calculates the */
/*     Euclidean norm of x with stride incx. */

/*     The function statement is */

/*       double precision function dnrm2(n,x,incx) */

/*     where */

/*       n is a positive integer input variable. */

/*       x is an input array of length n. */

/*       incx is a positive integer variable that specifies the */
/*         stride of the vector. */

/*     Subprograms called */

/*       FORTRAN-supplied ... abs, max, sqrt */

/*     MINPACK-2 Project. February 1991. */
/*     Argonne National Laboratory. */
/*     Brett M. Averick. */

/*     ********** */
    /* Parameter adjustments */
    --x;

    /* Function Body */
    ret_val = 0.;
    scale = 0.;
    i__1 = *n;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MAX */
	d__2 = scale, d__3 = (d__1 = x[i__], abs(d__1));
	scale = max(d__2,d__3);
/* L10: */
    }
    if (scale == 0.) {
	return ret_val;
    }
    i__2 = *n;
    i__1 = *incx;
    for (i__ = 1; i__1 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__1) {
/* Computing 2nd power */
	d__1 = x[i__] / scale;
	ret_val += d__1 * d__1;
/* L20: */
    }
    ret_val = scale * sqrt(ret_val);
    return ret_val;
} /* dnrm2_ */

/* ====================== The end of dnrm2 =============================== */
doublereal dpmeps_()
{
    /* Initialized data */

    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal two = 2.;

    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static doublereal beta;
    static integer irnd;
    static doublereal temp, temp1, a, b;
    static integer i__;
    static doublereal betah;
    static integer ibeta, negep;
    static doublereal tempa;
    static integer itemp, it;
    static doublereal betain;

/*     ********** */

/*     Subroutine dpeps */

/*     This subroutine computes the machine precision parameter */
/*     dpmeps as the smallest floating point number such that */
/*     1 + dpmeps differs from 1. */

/*     This subroutine is based on the subroutine machar described in */

/*     W. J. Cody, */
/*     MACHAR: A subroutine to dynamically determine machine parameters, */
/*     ACM Transactions on Mathematical Software, 14, 1988, pages 303-311. */

/*     The subroutine statement is: */

/*       subroutine dpeps(dpmeps) */

/*     where */

/*       dpmeps is a double precision variable. */
/*         On entry dpmeps need not be specified. */
/*         On exit dpmeps is the machine precision. */

/*     MINPACK-2 Project. February 1991. */
/*     Argonne National Laboratory and University of Minnesota. */
/*     Brett M. Averick. */

/*     ******* */
/*     determine ibeta, beta ala malcolm. */
    a = one;
    b = one;
L10:
    a += a;
    temp = a + one;
    temp1 = temp - a;
    if (temp1 - one == zero) {
	goto L10;
    }
L20:
    b += b;
    temp = a + b;
    itemp = (integer) (temp - a);
    if (itemp == 0) {
	goto L20;
    }
    ibeta = itemp;
    beta = (doublereal) ibeta;
/*     determine it, irnd. */
    it = 0;
    b = one;
L30:
    ++it;
    b *= beta;
    temp = b + one;
    temp1 = temp - b;
    if (temp1 - one == zero) {
	goto L30;
    }
    irnd = 0;
    betah = beta / two;
    temp = a + betah;
    if (temp - a != zero) {
	irnd = 1;
    }
    tempa = a + beta;
    temp = tempa + betah;
    if (irnd == 0 && temp - tempa != zero) {
	irnd = 2;
    }
/*     determine dpmeps. */
    negep = it + 3;
    betain = one / beta;
    a = one;
    i__1 = negep;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a *= betain;
/* L40: */
    }
L50:
    temp = one + a;
    if (temp - one != zero) {
	goto L60;
    }
    a *= beta;
    goto L50;
L60:
    ret_val = a;
    if (ibeta == 2 || irnd == 0) {
	goto L70;
    }
    a = a * (one + a) / two;
    temp = one + a;
    if (temp - one != zero) {
	ret_val = a;
    }
L70:
    return ret_val;
} /* dpmeps_ */

/* ====================== The end of dpmeps ============================== */
/* Subroutine */ int daxpy_(n, da, dx, incx, dy, incy)
integer *n;
doublereal *da, *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;


/*     ant times a vector plus a vector. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*da == 0.) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[iy] += *da * dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[i__] += *da * dx[i__];
/* L30: */
    }
    if (*n < 4) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 4) {
	dy[i__] += *da * dx[i__];
	dy[i__ + 1] += *da * dx[i__ + 1];
	dy[i__ + 2] += *da * dx[i__ + 2];
	dy[i__ + 3] += *da * dx[i__ + 3];
/* L50: */
    }
    return 0;
} /* daxpy_ */

/* ====================== The end of daxpy =============================== */
/* Subroutine */ int dcopy_(n, dx, incx, dy, incy)
integer *n;
doublereal *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, m, ix, iy, mp1;


/*     copies a vector, x, to a vector, y. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[iy] = dx[ix];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return 0;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 7;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[i__] = dx[i__];
/* L30: */
    }
    if (*n < 7) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 7) {
	dy[i__] = dx[i__];
	dy[i__ + 1] = dx[i__ + 1];
	dy[i__ + 2] = dx[i__ + 2];
	dy[i__ + 3] = dx[i__ + 3];
	dy[i__ + 4] = dx[i__ + 4];
	dy[i__ + 5] = dx[i__ + 5];
	dy[i__ + 6] = dx[i__ + 6];
/* L50: */
    }
    return 0;
} /* dcopy_ */

/* ====================== The end of dcopy =============================== */
doublereal ddot_(n, dx, incx, dy, incy)
integer *n;
doublereal *dx;
integer *incx;
doublereal *dy;
integer *incy;
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val;

    /* Local variables */
    static integer i__, m;
    static doublereal dtemp;
    static integer ix, iy, mp1;


/*     forms the dot product of two vectors. */
/*     uses unrolled loops for increments equal to one. */
/*     jack dongarra, linpack, 3/11/78. */


    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    ret_val = 0.;
    dtemp = 0.;
    if (*n <= 0) {
	return ret_val;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments */
/*          not equal to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[ix] * dy[iy];
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    ret_val = dtemp;
    return ret_val;

/*        code for both increments equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dtemp += dx[i__] * dy[i__];
/* L30: */
    }
    if (*n < 5) {
	goto L60;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i__ = mp1; i__ <= i__1; i__ += 5) {
	dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] + dx[
		i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + dx[i__ + 
		4] * dy[i__ + 4];
/* L50: */
    }
L60:
    ret_val = dtemp;
    return ret_val;
} /* ddot_ */

/* ====================== The end of ddot ================================ */
/* Subroutine */ int dpofa_(a, lda, n, info)
doublereal *a;
integer *lda, *n, *info;
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    extern doublereal ddot_();
    static integer j, k;
    static doublereal s, t;
    static integer jm1;


/*     dpofa factors a double precision symmetric positive definite */
/*     matrix. */

/*     dpofa is usually called by dpoco, but it can be called */
/*     directly with a saving in time if  rcond  is not needed. */
/*     (time for dpoco) = (1 + 18/n)*(time for dpofa) . */

/*     on entry */

/*        a       double precision(lda, n) */
/*                the symmetric matrix to be factored.  only the */
/*                diagonal and upper triangle are used. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*     on return */

/*        a       an upper triangular matrix  r  so that  a = trans(r)*r */
/*                where  trans(r)  is the transpose. */
/*                the strict lower triangle is unaltered. */
/*                if  info .ne. 0 , the factorization is not complete. */

/*        info    integer */
/*                = 0  for normal return. */
/*                = k  signals an error condition.  the leading minor */
/*                     of order  k  is not positive definite. */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas ddot */
/*     fortran sqrt */

/*     internal variables */

/*     begin block with ...exits to 40 */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1 * 1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = k - 1;
	    t = a[k + j * a_dim1] - ddot_(&i__3, &a[k * a_dim1 + 1], &c__1, &
		    a[j * a_dim1 + 1], &c__1);
	    t /= a[k + k * a_dim1];
	    a[k + j * a_dim1] = t;
	    s += t * t;
/* L10: */
	}
L20:
	s = a[j + j * a_dim1] - s;
/*     ......exit */
	if (s <= 0.) {
	    goto L40;
	}
	a[j + j * a_dim1] = sqrt(s);
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* dpofa_ */

/* ====================== The end of dpofa =============================== */
/* Subroutine */ int dscal_(n, da, dx, incx)
integer *n;
doublereal *da, *dx;
integer *incx;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, m, nincx, mp1;


/*     scales a vector by a ant. */
/*     uses unrolled loops for increment equal to one. */
/*     jack dongarra, linpack, 3/11/78. */
/*     modified 3/93 to return if incx .le. 0. */


    /* Parameter adjustments */
    --dx;

    /* Function Body */
    if (*n <= 0 || *incx <= 0) {
	return 0;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    nincx = *n * *incx;
    i__1 = nincx;
    i__2 = *incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
	dx[i__] = *da * dx[i__];
/* L10: */
    }
    return 0;

/*        code for increment equal to 1 */


/*        clean-up loop */

L20:
    m = *n % 5;
    if (m == 0) {
	goto L40;
    }
    i__2 = m;
    for (i__ = 1; i__ <= i__2; ++i__) {
	dx[i__] = *da * dx[i__];
/* L30: */
    }
    if (*n < 5) {
	return 0;
    }
L40:
    mp1 = m + 1;
    i__2 = *n;
    for (i__ = mp1; i__ <= i__2; i__ += 5) {
	dx[i__] = *da * dx[i__];
	dx[i__ + 1] = *da * dx[i__ + 1];
	dx[i__ + 2] = *da * dx[i__ + 2];
	dx[i__ + 3] = *da * dx[i__ + 3];
	dx[i__ + 4] = *da * dx[i__ + 4];
/* L50: */
    }
    return 0;
} /* dscal_ */

/* ====================== The end of dscal =============================== */
/* Subroutine */ int dtrsl_(t, ldt, n, b, job, info)
doublereal *t;
integer *ldt, *n;
doublereal *b;
integer *job, *info;
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2;

    /* Local variables */
    static integer case__;
    extern doublereal ddot_();
    static doublereal temp;
    static integer j;
    extern /* Subroutine */ int daxpy_();
    static integer jj;



/*     dtrsl solves systems of the form */

/*                   t * x = b */
/*     or */
/*                   trans(t) * x = b */

/*     where t is a triangular matrix of order n. here trans(t) */
/*     denotes the transpose of the matrix t. */

/*     on entry */

/*         t         double precision(ldt,n) */
/*                   t contains the matrix of the system. the zero */
/*                   elements of the matrix are not referenced, and */
/*                   the corresponding elements of the array can be */
/*                   used to store other information. */

/*         ldt       integer */
/*                   ldt is the leading dimension of the array t. */

/*         n         integer */
/*                   n is the order of the system. */

/*         b         double precision(n). */
/*                   b contains the right hand side of the system. */

/*         job       integer */
/*                   job specifies what kind of system is to be solved. */
/*                   if job is */

/*                        00   solve t*x=b, t lower triangular, */
/*                        01   solve t*x=b, t upper triangular, */
/*                        10   solve trans(t)*x=b, t lower triangular, */
/*                        11   solve trans(t)*x=b, t upper triangular. */

/*     on return */

/*         b         b contains the solution, if info .eq. 0. */
/*                   otherwise b is unaltered. */

/*         info      integer */
/*                   info contains zero if the system is nonsingular. */
/*                   otherwise info contains the index of */
/*                   the first zero diagonal element of t. */

/*     linpack. this version dated 08/14/78 . */
/*     g. w. stewart, university of maryland, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,ddot */
/*     fortran mod */

/*     internal variables */


/*     begin block permitting ...exits to 150 */

/*        check for zero diagonal elements. */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1 * 1;
    t -= t_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (*info = 1; *info <= i__1; ++(*info)) {
/*     ......exit */
	if (t[*info + *info * t_dim1] == 0.) {
	    goto L150;
	}
/* L10: */
    }
    *info = 0;

/*        determine the task and go to it. */

    case__ = 1;
    if (*job % 10 != 0) {
	case__ = 2;
    }
    if (*job % 100 / 10 != 0) {
	case__ += 2;
    }
    switch ((int)case__) {
	case 1:  goto L20;
	case 2:  goto L50;
	case 3:  goto L80;
	case 4:  goto L110;
    }

/*        solve t*x=b for t lower triangular */

L20:
    b[1] /= t[t_dim1 + 1];
    if (*n < 2) {
	goto L40;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	temp = -b[j - 1];
	i__2 = *n - j + 1;
	daxpy_(&i__2, &temp, &t[j + (j - 1) * t_dim1], &c__1, &b[j], &c__1);
	b[j] /= t[j + j * t_dim1];
/* L30: */
    }
L40:
    goto L140;

/*        solve t*x=b for t upper triangular. */

L50:
    b[*n] /= t[*n + *n * t_dim1];
    if (*n < 2) {
	goto L70;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
	j = *n - jj + 1;
	temp = -b[j + 1];
	daxpy_(&j, &temp, &t[(j + 1) * t_dim1 + 1], &c__1, &b[1], &c__1);
	b[j] /= t[j + j * t_dim1];
/* L60: */
    }
L70:
    goto L140;

/*        solve trans(t)*x=b for t lower triangular. */

L80:
    b[*n] /= t[*n + *n * t_dim1];
    if (*n < 2) {
	goto L100;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
	j = *n - jj + 1;
	i__2 = jj - 1;
	b[j] -= ddot_(&i__2, &t[j + 1 + j * t_dim1], &c__1, &b[j + 1], &c__1);
	b[j] /= t[j + j * t_dim1];
/* L90: */
    }
L100:
    goto L140;

/*        solve trans(t)*x=b for t upper triangular. */

L110:
    b[1] /= t[t_dim1 + 1];
    if (*n < 2) {
	goto L130;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j - 1;
	b[j] -= ddot_(&i__2, &t[j * t_dim1 + 1], &c__1, &b[1], &c__1);
	b[j] /= t[j + j * t_dim1];
/* L120: */
    }
L130:
L140:
L150:
    return 0;
} /* dtrsl_ */


