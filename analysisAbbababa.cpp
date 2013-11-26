/*
  thorfinn thorfinn@binf.ku.dk 
  anders albrecht@binf.ku.dk made this.
  part of angsd
*/


#include "shared.h"
#include <cmath>
#include "analysisFunction.h"


class abbababa:public general{
public:
  int doAbbababa;
  int nInd;
  int minQ;
  int sample;
  //none optional stuff
  gzFile outfileZ;
  FILE *outfile;
  int currentChr;
  abbababa(const char *outfiles,argStruct *arguments,int inputtype);
  ~abbababa();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void print(funkyPars *pars);
  void clean(funkyPars *pars);
  void printArg(FILE *argFile);

  size_t **alleleCounts; //[ind][125]; 
  size_t **alleleCountsChr; //[ind][125]; 
};

void abbababa::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doAbbababa\t%d\n",doAbbababa);
  fprintf(argFile,"\t1: EM v1\n");
  fprintf(argFile,"\t2: EM v2\n");
  fprintf(argFile,"\t-minQ\t\t%d\t(remove bases with qscore<minQ)\n",minQ);
  fprintf(argFile,"\t-sample\t\t%d\tSample only 1 read\n",sample);
  fprintf(argFile,"\n");
}

void abbababa::getOptions(argStruct *arguments){

  //from command line
  doAbbababa=angsd::getArg("-doAbbababa",doAbbababa,arguments);
  minQ=angsd::getArg("-minQ",minQ,arguments);
  sample=angsd::getArg("-sample",sample,arguments);
  nInd=arguments->nInd;
  
  //inputtype 0) soap 1) samglf 2) samglfClean 3) tglf 4)sim type 5) beagle
  if(doAbbababa){
    if(arguments->inputtype!=0&&arguments->inputtype!=7){
      fprintf(stderr,"Error: bam or soap input needed for -doAbbababa \n");
      exit(0);
    }
  }
  if(doAbbababa & nInd >1){
    fprintf(stderr,"Error: Only use a single individual\n");
  }
}

abbababa::abbababa(const char *outfiles,argStruct *arguments,int inputtype){
  doAbbababa=0;
  minQ=MINQ;//<-general.h
  sample=1;
  currentChr=-1;
  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-doAbbababa")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
  getOptions(arguments);
  printArg(arguments->argumentFile);

  if(doAbbababa==0)
    return;

 

  //make output files
  const char* postfix;
  postfix=".abbababa.gz";
  outfileZ = openFileGz(outfiles,postfix,GZOPT);
  const char* postfix2;
  postfix2=".abbababa2";
  outfile = openFile(outfiles,postfix2);

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


abbababa::~abbababa(){

  if(doAbbababa==0)
    return;

  if(doAbbababa==1){
    for(int i=0;i<nInd;i++){
      for(int j=0;j<125;j++)
	gzprintf(outfileZ,"%lu\t",alleleCounts[i][j]);
      gzprintf(outfileZ,"\n");
    }
  }


  for(int i=0;i<nInd;i++)
    delete[]  alleleCounts[i];
  delete [] alleleCounts; 

  for(int i=0;i<nInd;i++)
    delete[]  alleleCountsChr[i];
  delete [] alleleCountsChr; 

  gzclose(outfileZ);
  fclose(outfile);
}


void abbababa::clean(funkyPars *pars){

}



void abbababa::print(funkyPars *pars){


  if(doAbbababa==0)
    return;

  /*
  for(int i=0;i<pars->nInd;i++){
    for(int j=0;j<pars->chk->nd[s][i].l;j++){
      if(pars->chk->nd[s][i].qs[j]>=minQ){
	alleleCounts[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i].seq[j]]]++;
	alleleCountsChr[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i].seq[j]]]++;
	break;
      }
    }
  }
  */  


  if(doAbbababa==1){
    if(currentChr==-1)
      currentChr=pars->refId;
    if(currentChr!=pars->refId){
      fprintf(outfile,"Chr: \t %s\n",header->name[currentChr]);
      for(int i=0;i<nInd;i++){
	for(int j=0;j<125;j++)
	  fprintf(outfile,"%lu\t",alleleCountsChr[i][j]);
	fprintf(outfile,"\n");
      }
      for(int i=0;i<nInd;i++)
	for(int j=0;j<256;j++)
	  alleleCountsChr[i][j]=0;
      currentChr=pars->refId;
    }

  }
}


void abbababa::run(funkyPars *pars){

  if(doAbbababa==0)
    return;

}

