/*
  thorfinn thorfinn@binf.ku.dk dec17 2012
 
    
  anders albrecht@binf.ku.dk made this.

  part of angsd
*/


#include "shared.h"
#include <cmath>
#include "analysisFunction.h"


class ansErr:public general{
public:
  int doAnsError;
  int nInd;
  int minQ;
  int sample;
  //none optional stuff
  FILE *outfile;
  FILE *outfile2;
  int currentChr;
  ansErr(const char *outfiles,argStruct *arguments,int inputtype);
  ~ansErr();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void print(funkyPars *pars);
  void clean(funkyPars *pars);
  void printArg(FILE *argFile);

  size_t **alleleCounts; //[ind][125]; 
  size_t **alleleCountsChr; //[ind][125]; 
};

void ansErr::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doAnsError\t%d\n",doAnsError);
  fprintf(argFile,"\t1: EM v1\n");
  fprintf(argFile,"\t2: EM v2\n");
  fprintf(argFile,"\t-minQ\t\t%d\t(remove bases with qscore<minQ)\n",minQ);
  fprintf(argFile,"\t-sample\t\t%d\tSample only 1 read\n",sample);
  fprintf(argFile,"\n");
}

void ansErr::getOptions(argStruct *arguments){

  //from command line
  doAnsError=angsd::getArg("-doAnsError",doAnsError,arguments);
  minQ=angsd::getArg("-minQ",minQ,arguments);
  sample=angsd::getArg("-sample",sample,arguments);
  nInd=arguments->nInd;
  
  //inputtype 0) soap 1) samglf 2) samglfClean 3) tglf 4)sim type 5) beagle
  if(doAnsError){
    if(arguments->inputtype!=0&&arguments->inputtype!=7){
      fprintf(stderr,"Error: bam or soap input needed for -doAnsError \n");
      exit(0);
    }
  }

}

ansErr::ansErr(const char *outfiles,argStruct *arguments,int inputtype){
  doAnsError=0;
  minQ=MINQ;//<-general.h
  sample=1;
  currentChr=-1;
  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-doAnsError")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
  getOptions(arguments);
  printArg(arguments->argumentFile);

  if(doAnsError==0)
    return;

 

  //make output files
  const char* postfix;
  postfix=".ansError";
  outfile = openFile(outfiles,postfix);
  const char* postfix2;
  postfix2=".ansErrorChr";
  outfile2 = openFile(outfiles,postfix2);

  //allocate allele counts
  alleleCounts = new size_t *[nInd];
  for(int i=0;i<nInd;i++)
    alleleCounts[i] = new size_t [256];
  for(int i=0;i<nInd;i++)
    for(int j=0;j<256;j++)
      alleleCounts[i][j]=0;

  alleleCountsChr = new size_t *[nInd];
  for(int i=0;i<nInd;i++)
    alleleCountsChr[i] = new size_t [256];
  for(int i=0;i<nInd;i++)
    for(int j=0;j<256;j++)
      alleleCountsChr[i][j]=0;

}


ansErr::~ansErr(){

  if(doAnsError==0)
    return;

  if(doAnsError==1){
    for(int i=0;i<nInd;i++){
      for(int j=0;j<125;j++)
	fprintf(outfile,"%lu\t",alleleCounts[i][j]);
      fprintf(outfile,"\n");
    }
  }


  for(int i=0;i<nInd;i++)
    delete[]  alleleCounts[i];
  delete [] alleleCounts; 

  for(int i=0;i<nInd;i++)
    delete[]  alleleCountsChr[i];
  delete [] alleleCountsChr; 

  fclose(outfile);
  fclose(outfile2);
}


void ansErr::clean(funkyPars *pars){

}



void ansErr::print(funkyPars *pars){


  if(doAnsError==0)
    return;

  if(doAnsError==2){
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
      fprintf(outfile,"%c",intToRef[pars->anc[s]]);
      fprintf(outfile,"\t%c",intToRef[pars->ref[s]]);

     for(int i=0;i<pars->nInd;i++){
	
	//      tNode &nd = pars->chk->nd[s][i];
	int allele=4;
	for(int j=0;j<pars->chk->nd[s][i].l;j++){
	  if(pars->chk->nd[s][i].qs[j]>=minQ){
	     allele=refToInt[pars->chk->nd[s][i].seq[j]];
	    break;
	  }
	}

	fprintf(outfile,"\t%c",intToRef[allele]);
      }
      fprintf(outfile,"\n");
    }
  }
  else if(doAnsError==1){
   
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;

      if(sample){
	for(int i=0;i<pars->nInd;i++){
	  for(int j=0;j<pars->chk->nd[s][i].l;j++){
	    if(pars->chk->nd[s][i].qs[j]>=minQ){
	      alleleCounts[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i].seq[j]]]++;
	      alleleCountsChr[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i].seq[j]]]++;
	      break;
	    }
	  }
	}
      }
      else{

	for(int i=0;i<pars->nInd;i++){
	  for(int j=0;j<pars->chk->nd[s][i].l;j++){
	    if(pars->chk->nd[s][i].qs[j]>=minQ){
	      alleleCounts[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i].seq[j]]]++;
	      alleleCountsChr[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i].seq[j]]]++;
	    }
	  }
	}
      }
    }
  }


  if(doAnsError==1){
    if(currentChr==-1)
      currentChr=pars->refId;
    if(currentChr!=pars->refId){
      fprintf(outfile2,"Chr: \t %s\n",header->name[currentChr]);
      for(int i=0;i<nInd;i++){
	for(int j=0;j<125;j++)
	  fprintf(outfile2,"%lu\t",alleleCountsChr[i][j]);
	fprintf(outfile2,"\n");
      }
      for(int i=0;i<nInd;i++)
	for(int j=0;j<256;j++)
	  alleleCountsChr[i][j]=0;
      currentChr=pars->refId;
    }

  }
}


void ansErr::run(funkyPars *pars){

  if(doAnsError==0)
    return;

}

