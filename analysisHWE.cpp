/*
  thorfinn thorfinn@binf.ku.dk 19dec 2012

  Authors of this file
  Filipe,anders

  part of angsd

 */
#include <cmath>
#include "shared.h"

#include "analysisFunction.h"
#include "general.h"
#include "analysisHWE.h"

void hwe::printArg(FILE *argFile){
  fprintf(argFile,"-------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doHWE\t%d\n",doHWE);
  fprintf(argFile,"\n");
}

void hwe::getOptions(argStruct *arguments){
  int tmpDoMaf =0;
  int GL =0;
  doHWE=angsd::getArg("-doHWE",doHWE,arguments);
  if(doHWE==0)
    return;
  GL=angsd::getArg("-GL",GL,arguments);
  tmpDoMaf=angsd::getArg("-doMaf",tmpDoMaf,arguments);
  if(doHWE&&(tmpDoMaf==0)){
    fprintf(stderr,"You supplied -doHWE, you should also choose -doMaf\n");
    exit(0);
  }
  if(arguments->inputtype==5){
    fprintf(stderr,"Error: you cannot estimate HWE based on posterior probabilities (beagle)\n");
    exit(0);
  }
  if(arguments->inputtype==1&&GL==0){
    fprintf(stderr,"Error: you need to estimate genotype likelihoods (GL) to estimate HWE \n");
    exit(0);
  }

}

hwe::hwe(const char *outfiles,argStruct *arguments,int inputtype){

  doHWE=0;
  
  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-doHWE")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }




  getOptions(arguments);
  printArg(arguments->argumentFile);
  if(doHWE==0)
    return;

  //make output files
  const char* postfix;
  postfix=".hwe.gz";
  if(doHWE>0){
    outfileZ = openFileGz(outfiles,postfix,GZOPT);
    //print header
    gzprintf(outfileZ,"Chromo\tPosition\tMajor\tMinor\tFreq\thweFreq\tF\tLRT\n");
  }
}


hwe::~hwe(){

  if(doHWE==0)
    return;
  if(doHWE>0)
    gzclose(outfileZ);
}


void hwe::clean(funkyPars *pars){
  if(doHWE==0)
    return;

  funkyHWE *hweStruct =(funkyHWE *) pars->extras[index];
  delete[] hweStruct->freq;
  delete[] hweStruct->F;
  delete[] hweStruct->like0;
  delete[] hweStruct->likeF;
  delete hweStruct;
  
}

void hwe::print(funkyPars *pars){
  if(doHWE<=0)
    return;

  funkyHWE *hweStruct = (funkyHWE *) pars->extras[index];//new

  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0) 
      continue;
    gzprintf(outfileZ,"%s\t%d\t%c\t%c\t%f\t%f\t%f\t%f\n",header->name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],pars->results->asso->freq[s],hweStruct->freq[s],hweStruct->F[s],2*hweStruct->like0[s]-2*hweStruct->likeF[s]);

  }

}


void hwe::run(funkyPars *pars){
 
  if(doHWE==0)
    return;

  //  pars->hweStruct = new funkyHWE;//old
  funkyHWE *hweStruct = new funkyHWE;//new

  double *freq = new double[pars->numSites];
  double *F = new double[pars->numSites];
  double *like0 = new double[pars->numSites];
  double *likeF = new double[pars->numSites];

  double **loglike3;
  loglike3=angsd::get3likes(pars);

  for(int s=0;s<pars->numSites;s++){
    if(pars->keepSites[s]==0) 
      continue;
    
    //start parameter
    double x[2];
    x[0]=0.05;
    x[1]=0.05;
    estHWE(x,loglike3[s],pars->nInd);
    freq[s]=x[0];
    F[s]=x[1];
    likeF[s] = HWE_like(x,loglike3[s],pars->nInd);
    x[1]=0.0;
    x[0]=pars->results->asso->freq[s];
    like0[s] = HWE_like(x,loglike3[s],pars->nInd);
    //    fprintf(stderr,"%f\t%f\n",x[0],x[1]);
    //fprintf(stderr,"%f\t%f\t%f\n",loglike3[s][0],loglike3[s][1],loglike3[s][2]);

  }
  //old
  //  pars->hweStruct->freq=freq;
  //  pars->hweStruct->F=F;


  hweStruct->freq=freq;
  hweStruct->F=F;
  hweStruct->like0=like0;
  hweStruct->likeF=likeF;
  pars->extras[index] = hweStruct;


  for(int s=0;s<pars->numSites;s++)
    delete[] loglike3[s];
  delete[] loglike3;

}



void hwe::HWE_EM(double *x,double *loglike,int nInd){
  double freq=x[0];
  double F=x[1];
  double p0=(pow(1-freq,2)+freq*(1-freq)*F);
  double p1=(2*freq*(1-freq)-2*freq*(1-freq)*F);
  double p2=(pow(freq,2)+freq*(1-freq)*F);

  double freq2=0;
  double F2=0;
  double norm;
  double oHO=0;
  double eHO=nInd*(pow(freq,2)+pow(1-freq,2));

  for(int i=0;i<nInd;i++){
    norm=angsd::addProtect3(log(p0)+loglike[i*3+0],log(p1)+loglike[i*3+1],log(p2)+loglike[i*3+2]);
    freq2+=exp(log(p1)+loglike[i*3+1]-norm)+exp(log(2)+log(p2)+loglike[i*3+2]-norm);
    oHO+=exp(log(p0)+loglike[i*3+0]-norm)+exp(log(p2)+loglike[i*3+2]-norm);
  //oHO+=(p0*loglike[i*3+0]+p2*loglike[i*3+2])/(2*norm);
    //fprintf(stderr,"1 %f 2 %f 3 %f\n",loglike[i*3+0],loglike[i*3+1],loglike[i*3+2]);
  }
  F2=(oHO-eHO)/(nInd-eHO);
  if(F2<0)
    F2=0;
  if(F2>1)
    F2=1;
  
  freq2=freq2/(2*nInd);
  //fprintf(stderr,"norm %f\tx[0] %f\tx[1] %f\t1 %f 2 %f 3 %f oHO %f eHO %f freq %f\n",norm,x[0],x[1],p0,p1,p2,oHO,eHO,freq2);

  if(freq2<0.0001)
    F2=0;
  //    F2=0;
  x[0]=freq2;
  x[1]=F2;
  if(freq2>1.0000001){
    fprintf(stderr,"something is wrong i HWE\t freq %f\n",freq2);
    fflush(stderr);
    exit(0);

  }
}




double hwe::HWE_like(double *x,double *loglike,int nInd){
  double freq=x[0];
  double F=x[1];
  double p0=(pow(1-freq,2)+freq*(1-freq)*F);
  double p1=(2*freq*(1-freq)-2*freq*(1-freq)*F);
  double p2=(pow(freq,2)+freq*(1-freq)*F);
  double totalLogLike=0;
  for(int i=0;i<nInd;i++)
    totalLogLike+=angsd::addProtect3(log(p0)+loglike[i*3+0],log(p1) + loglike[i*3+1],log(p2) + loglike[i*3+2]);
  return totalLogLike;
}

void hwe::estHWE(double *x,double *loglike,int nInd){
  double l=HWE_like(x,loglike,nInd);
  int iter=50;
  double d=l;
  int printer=0;
  for(int i=0;i<iter;i++){
    HWE_EM(x,loglike,nInd);
    l=HWE_like(x,loglike,nInd);
    if(d>l+0.01){
      //   fprintf(stderr,"d %f\tl %f\n",d,l);
      printer=1;
    }
    d=l;
  }
  if(printer & 0){
    x[0]=0.05;
    x[1]=0.05;
    l=HWE_like(x,loglike,nInd);
    for(int i=0;i<iter;i++){
      fprintf(stderr,"like %d\t%f\tf %f\tF %f\t%d\n",i,l,x[0],x[1],nInd);
      HWE_EM(x,loglike,nInd);
      l=HWE_like(x,loglike,nInd);
    }

    for(int ind=0;ind<nInd;ind++){
      fprintf(stdout,"likes %f, %f, %f\n",loglike[ind*3+0],loglike[ind*3+1],loglike[ind*3+2]);
    }
    fflush(stderr);
    exit(0);

  }

  //print here
  //  fprintf(stderr,"like %d\t%f\tf %f\tF %f\t%d\n",iter,l,x[0],x[1],nInd);
}
/*
x<-c(0.05,0.05)
for(tal in 1:200){
l<-getLike(x,exp(ll))
cat(x[1]," ",x[2]," ",l,"\n")
x<-em(x,exp(ll))
}

 */
