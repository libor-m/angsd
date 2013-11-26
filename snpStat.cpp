
#include <cmath>
#include <ctype.h>
#include "analysisFunction.h"
#include "shared.h"
#include "fet.c"



class snpStat:public general{
public:
  int doSnpStat;
  
  snpStat(const char *outfiles,argStruct *arguments,int inputtype);
  ~snpStat();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  
  
};
void snpStat::printArg(FILE *fp){
  fprintf(fp,"doSnpStat=%d\n",doSnpStat);
  
}
void snpStat::run(funkyPars *pars){
  if(!doSnpStat)
    return;
  chunkyT *chk = pars->chk;
  
  if(doSnpStat==1){
    //loop over sites;
      for(int s=0;s<pars->numSites;s++){
	if(pars->keepSites[s]==0)
	  continue;
	//loop over samples
	int cnts[4]={0,0,0,0};
	for(int i=0;i<pars->nInd;i++){
	  tNode nd = chk->nd[s][i];
	  for(int l=0;l<nd.l;l++){
	    int obB = refToInt[nd.seq[l]];
	    int strand = (isupper(nd.seq[l])==0)<<1;
	    if(obB==4)
	      continue;
	    if((obB!=pars->major[s] && obB!=pars->minor[s]) )
	      continue;
	    if(obB==pars->major[s])
	      strand +=1;
	    cnts[strand]++;
	  }
	}
	double left,right,twotail,prob;
	prob = kt_fisher_exact(cnts[0], cnts[1], cnts[2], cnts[3], &left, &right, &twotail);
	fprintf(stderr,"res\t posi=%d\t%d %d %d %d %f\n",pars->posi[s], cnts[0],cnts[1],cnts[2],cnts[3],twotail);
      }
  }
  
}

void snpStat::clean(funkyPars *fp){
  if(!doSnpStat)
    return;

  
}

void snpStat::print(funkyPars *fp){
  if(!doSnpStat)
    return;
      
}


void snpStat::getOptions(argStruct *arguments){
  //default
  doSnpStat=0;

  //from command line
  doSnpStat=angsd::getArg("-doSnpStat",doSnpStat,arguments);
  if(doSnpStat==-999){
    doSnpStat=0;
    printArg(stderr);
    exit(0);
  }
  if(doSnpStat==0)
    return;
  int domajmin=0;
  //from command line
  domajmin=angsd::getArg("-doSNP",domajmin,arguments);
  if(!domajmin){
    fprintf(stderr,"-doSnpStat require -doSNP\n");
    exit(0);
  }
    
  printArg(arguments->argumentFile);

}


snpStat::snpStat(const char *outfiles,argStruct *arguments,int inputtype){
  getOptions(arguments);
  if(doSnpStat)
    fprintf(stderr,"running doSnpStat=%d\n",doSnpStat);
}

snpStat::~snpStat(){


}
