
#include "shared.h"
#include "general.h"
#include "analysisFunction.h" 


//these are the major builtin analysis that angsd can perform
#include "analysisKeepList.h"
#include "analysisMajorMinor.cpp"
#include "analysisMaf.h"
#include "analysisEstError.h"
#include "analysisEstLikes.cpp"
#include "analysisAsso.cpp"
#include "analysisHWE.h"
#include "analysisAnsError.cpp"
#include "analysisCallGenotypes.cpp"
#include "getFasta.h"//for reading fasta; ancestral and refernce
#include "analysisCount.cpp" //generate counts from reads
//extrastuff
#include "angsd_realSFS.cpp"
#include "analysisCovar.cpp" 

#include "thorfinn.h"
#include "snpStat.h"

#include "snptools.cpp" //<-implemenation of some stuff from snptools. 
#include "hetplas.cpp" //<-implementation of hetero plasmic
int general::tot_index =0;
aHead *general::header = NULL;
std::map<char *,int,ltstr> *general::revMap = NULL;
general **extra(int &nItem,const char *outfiles,int inputtype,argStruct *arguments){
  int nit=0;
  //  printHd(hd,stderr);
  //change the number of method when adding a new one
  general **tskStuff =new general*[17];
  tskStuff[nit++] = new filter(arguments);
  tskStuff[nit++] = new getFasta(arguments);
  tskStuff[nit++] = new countCls(outfiles,arguments,inputtype);
  tskStuff[nit++] = new error(outfiles,arguments,inputtype);
  tskStuff[nit++] = new likeClass(outfiles,arguments,inputtype);
  tskStuff[nit++] = new majorminor(outfiles,arguments,inputtype);
  tskStuff[nit++] = new frequency(outfiles,arguments,inputtype);
  tskStuff[nit++] = new asso(outfiles,arguments,inputtype);
  tskStuff[nit++] = new hwe(outfiles,arguments,inputtype);
  tskStuff[nit++] = new ansErr(outfiles,arguments,inputtype);
  tskStuff[nit++] = new callGenotypes(outfiles,arguments,inputtype);
  tskStuff[nit++] = new realSFS(outfiles,arguments,inputtype);
  tskStuff[nit++] = new covar(outfiles);
  tskStuff[nit++] = new thorfinn(outfiles,arguments,inputtype);
  tskStuff[nit++] = new snpStat(outfiles,arguments,inputtype);
  tskStuff[nit++] = new snptools(outfiles,arguments,inputtype);
  tskStuff[nit++] = new hetplas(outfiles,arguments,inputtype);
  //add yours here:


  //don't touch below
  nItem = nit;
  return tskStuff;
}
