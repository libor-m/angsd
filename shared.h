
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <map>
#include <zlib.h>
#include "bambi_interface.h"
#include "bams.h" // <- only used for getting the header, 
#include "argStruct.h"

#ifndef _types_t
#define _types_t
//typedef short unsigned int suint;
typedef unsigned int suint;

#define LENS 10000
#define GZOPT "w6h"
// struct for covar class - classy




extern int refToInt[256];
extern char intToRef[5];



typedef struct{
  double *pml; 
  double *pmlun;
  double *pEM; 
  double *pEMun; 
  double *pmlSNP; 
  double *pmlunSNP;
  double *pEMSNP; 
  double *pEMunSNP;
}funkyPrintFreq;

typedef struct{

  double **stat;
  int **keepInd;
  double *freq;
  double *lrt_snp;
  int *highHe;
  int *highHo;
}funkyPrintAsso;


typedef struct {
  funkyPrintFreq *freq;
  funkyPrintAsso *asso;
  
}funkyPrint;


typedef struct {
  //dynamic
  //loci *sites;//<- this will become obsolete
  int refId;//<- this will be the future 4486
  int *posi;//<- this will be the future 4486

  suint **counts;//[nsites][5xiNind] !!!! nope! only 4
  double **likes;
  double **post;
  double *phat;
  
  char *major;
  char *minor;
  char *ref;
  char *anc;

  //primitives
  int numSites;
  int nInd;
  int **depth;

  //for print
  funkyPrint *results;

  //sites removed in analysis and print
  int *keepSites;
  int chunkNumber;
  
  //stuff needed for bamreader
  fcb *for_callback;
  chunkyT *chk;
  int killSig;
  
  //extra stuff associated with each analysis module
  void **extras;
  
}funkyPars;

funkyPars *allocFunkyPars();
void deallocFunkyPars(funkyPars *p);


void init(argStruct *arguments);//intialize all needed objects
void destroy();//destroy all initialized objects
void selector(funkyPars *p);

size_t fsize(const char* fname);
int fexists(const char* str);//{///@param str Filename given as a string.
FILE *openFile(const char* a,const char* b);
gzFile openFileGz(const char* a,const char* b,const char *mode);
FILE *getFILE(const char*fname,const char* mode);
gzFile getGz(const char*fname,const char* mode);

#define LOGMAX 20000   // pre-computed logfactorial 
#endif
