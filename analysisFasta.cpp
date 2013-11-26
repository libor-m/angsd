/*
  make a fasta file from a bam file. 
  thorfinn thorfinn@binf.ku.dk dec17 2012   
  anders albrecht@binf.ku.dk made this.
  part of angsd
*/


#include "shared.h"
#include <cmath>
#include "analysisFunction.h"
#include <cstdlib>
#include "kstring.h"//<-used for buffered output
#include "general.h"
typedef struct {
  char *seq;
  int start;
  int stop;
}funkyFasta;


class fasta:public general{
private:
  kstring_t bufstr;
  size_t currentPos;
  int currentChr;
public:
  int doFasta;
  int minQ;
  //none optional stuff
  FILE *outfile;
  FILE *outfile2;

  fasta(const char *outfiles,argStruct *arguments,int inputtype);
  ~fasta();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void print(funkyPars *pars);
  void clean(funkyPars *pars);
  void printArg(FILE *argFile);

};

void fasta::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doFasta\t%d\n",doFasta);
  fprintf(argFile,"\t1: use the most common base\n");
  fprintf(argFile,"\t2: use a random base\n");
  fprintf(argFile,"\t-minQ\t\t%d\t(remove bases with qscore<minQ)\n",minQ);
  fprintf(argFile,"\n");
}

void fasta::getOptions(argStruct *arguments){

  //from command line
  doFasta=angsd::getArg("-doFasta",doFasta,arguments);
  minQ=angsd::getArg("-minQ",minQ,arguments);
    
  //inputtype 0) soap 1) samglf 2) samglfClean 3) tglf 4)sim type 5) beagle
  if(doFasta){
    if(arguments->inputtype!=0&&arguments->inputtype!=7){
      fprintf(stderr,"Error: bam or soap input needed for -doFasta \n");
      exit(0);
    }
  }

}

fasta::fasta(const char *outfiles,argStruct *arguments,int inputtype){
  doFasta=0;
  minQ=MINQ;//<-general.h
  currentChr=-1;
  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-doFasta")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  
  getOptions(arguments);
  printArg(arguments->argumentFile);

  if(doFasta==0)
    return;

 

  //make output files
  const char* postfix;
  postfix=".fa";
  outfile = openFile(outfiles,postfix);
  const char* postfix2;
  postfix2=".fastaChr";
  outfile2 = openFile(outfiles,postfix2);

  bufstr.s=NULL;
  bufstr.m=0;
  bufstr.l=0;

  currentPos=0;

}


fasta::~fasta(){

  if(doFasta==0)
    return;

  /*
 if(doFasta==1){
    for(int i=0;i<nInd;i++){
      for(int j=0;j<125;j++)
	//	fprintf(outfile,"%lu\t",alleleCounts[i][j]);
	//      fprintf(outfile,"\n");
    }
  }
  */


  while(currentChr<header->n_ref){
    fprintf(stderr,"currentpos %lu %d\n",currentPos,header->l_ref[currentChr]); 
    for(int p=currentPos;p<header->l_ref[currentChr];p++){
      //if()
      fprintf(outfile,"N");
    }
    
    currentPos=0;
    currentChr++;
    if(currentChr<header->n_ref)
      fprintf(outfile,"\n>%s\n",header->name[currentChr]);
 
  }

  fclose(outfile);
  fclose(outfile2);
  if(bufstr.s!=NULL)
    free(bufstr.s);

}


void fasta::clean(funkyPars *pars){

  if(doFasta==0)
    return;

  funkyFasta *fastaStruct =(funkyFasta *) pars->extras[index];
  delete[] fastaStruct->seq;
  delete fastaStruct;

}



void fasta::print(funkyPars *pars){

  if(doFasta==0)
    return;

 funkyFasta *fastaStruct = (funkyFasta *) pars->extras[index];//new


 if(currentChr==-1){
   currentPos=0;
   currentChr=0;
   fprintf(outfile,">%s\n",header->name[currentChr]);
 }

 while(currentChr!=pars->refId){
   //   fprintf(stderr,"currentpos %lu %d\n",currentPos,header->l_ref[currentChr]);
  
   for(int p=currentPos;p<header->l_ref[currentChr];p++)
     fprintf(outfile,"N");

   currentPos=0;
   currentChr++;
   fprintf(outfile,"\n>%s\n",header->name[currentChr]);
 }
 for(int p=currentPos;p<fastaStruct->start;p++){
   fprintf(outfile,"N");
   currentPos++;
 }
 int c=0;
 for(int p=currentPos;p<fastaStruct->stop;p++){
   if(p=pars->posi[c]){
     fprintf(outfile,"%c",fastaStruct->seq[c]);
     //     fprintf(outfile,"%d%c\n",currentPos,fastaStruct->seq[c]);
     c++;
   }
   else
     fprintf(outfile,"N");
   currentPos++;
 }

  /*
  if(doFasta==0){//use all reads
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
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
  if(doFasta==1){//random number read
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
      for(int i=0;i<pars->nInd;i++){
	if(pars->chk->nd[s][i].l==0)
	  continue;
	int j = std::rand() % pars->chk->nd[s][i].l;
	if(pars->chk->nd[s][i].qs[j]>=minQ){
	  alleleCounts[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i].seq[j]]]++;
	  alleleCountsChr[i][pars->anc[s]*25+pars->ref[s]*5+refToInt[pars->chk->nd[s][i].seq[j]]]++;
	}
      }
    }
  }
  if(doFasta==2){//take first read
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
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
  }
  
   

  if(doFasta==1){
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

  */
}


void fasta::run(funkyPars *pars){

  if(doFasta==0)
    return;

  funkyFasta *fastaStruct = new funkyFasta;//new
  char *seq;
  int start;
  int stop;
  start = pars->posi[0];
  stop = pars->posi[pars->numSites-1];

  seq = new char[stop-start+1];

  for(int m=0;m<stop-start+1;m++)
    seq[m]='N';



  if(doFasta==2){//random number read
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
      if(pars->chk->nd[s][0].l==0)
	continue;
      int j = std::rand() % pars->chk->nd[s][0].l;
      if(pars->chk->nd[s][0].qs[j]>=minQ){
	seq[pars->posi[s]-start] = pars->chk->nd[s][0].seq[j];
      }
    }
  }

  fastaStruct->seq=seq;
  fastaStruct->start=start;
  fastaStruct->stop=stop;
  pars->extras[index] = fastaStruct;

}

