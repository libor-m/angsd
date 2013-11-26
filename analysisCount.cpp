/*
  thorfinn thorfinn@binf.ku.dk dec17 2012
  
  part of angsd

  Class that works with counts

  This function


*/


class countCls : public general{
private:
  const char* postfix1; //.pos
  const char* postfix2;//.counts;
  const char* postfix3;//.qs
  const char* postfix4;//.depthSample
  const char* postfix5;//.depthGlobal
  int dumpCounts;
  int doCounts;
  int doQsDist;//0=nothing 1=overall 
  //

  FILE *oFileCountsBin;
  FILE *oFileCountsPos;
  FILE *oFileQs;

  FILE *oFileGlobDepth;
  FILE *oFileSamplDepth;

  size_t *qsDist;
  int nInd;
  int minInd;
  int minQ;
   int trim;
  int setMaxDepth;
  //from depth class
  int doDepth;
  int maxDepth;
  int minDepth;
  size_t **depthCount;
  size_t *globCount;

public:

  
  //none optional stuff
  
  countCls(const char *outfiles,argStruct *arguments,int inputtype);
  ~countCls();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);

};
void countCls::printArg(FILE *argFile){
  fprintf(argFile,"---------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doCounts\t%d\t(Count the number A,C,G,T. All sites, All samples)\n",doCounts);
  
  
  fprintf(argFile,"\t-minQ\t\t%d\t(remove bases with qscore<minQ)\n",minQ);
  fprintf(argFile,"\t-setMaxDepth\t%d\t(-1 indicates no filtering)\n",setMaxDepth);
  fprintf(argFile,"\t-trim\t\t%d\t(trim ends of reads)\n",trim);
  fprintf(argFile,"\t-minDepth\t%d\t(bin together high depths)\n",minDepth);
  fprintf(argFile,"\t-minInd\t\t%d\t(0 indicates no filtering)\n",minInd);
  fprintf(argFile,"\nFiledumping:\n");
  fprintf(argFile,"\t-doDepth\t%d\t(dump distribution of seqdepth)\t%s,%s\n",doDepth,postfix4,postfix5);
  fprintf(argFile,"\t  -maxDepth\t%d\t(bin together high depths)\n",maxDepth);
  
  fprintf(argFile,"\t-doQsDist\t%d\t(dump distribution of qscores)\t%s\n",doQsDist,postfix3);  
  fprintf(argFile,"\t-dumpCounts\t%d\n",dumpCounts);
  fprintf(argFile,"\t  1: total seqdepth for site\t%s\n",postfix1);
  fprintf(argFile,"\t  2: seqdepth persample\t\t%s,%s\n",postfix1,postfix2);
  fprintf(argFile,"\t  3: A,C,G,T overall all samples\t%s,%s\n",postfix1,postfix2);
  fprintf(argFile,"\t  4: A,C,G,T for all samples\t%s,%s\n",postfix1,postfix2);

  fprintf(argFile,"\n");
  fprintf(argFile,"\n");

 
}

int calcSum(suint *counts,int len){
  int tmp=0;
  for(int i=0;i<len;i++)
    tmp += counts[i];
  return tmp;
}

void printCounts(char *chr,int *posi,suint **counts,int nSites,size_t nInd,FILE *oFileCountsPos,FILE *oFileCountsBin,int dumpType,int *keepSites){
  for(int s=0;s<nSites;s++){
    if(keepSites[s]==0)
      continue;
    fprintf(oFileCountsPos,"%s\t%d\t%d\n",chr,posi[s]+1,calcSum(counts[s],4*nInd));
    
    //if we need per sample info
    if(oFileCountsBin) {
      if(dumpType==4)//count A,C,G,T
	for(int i=0;i<4*nInd;i++)
	  fprintf(oFileCountsBin,"%u\t",counts[s][i]);
      else if(dumpType==2){//print A+C+G+T
	for(int n=0;n<nInd;n++)
	  fprintf(oFileCountsBin,"%u\t",counts[s][n*4]+counts[s][n*4+1]+counts[s][n*4+2]+counts[s][n*4+3]);
      }else{//overall sum of A,C,G,T
	size_t tsum[4]={0,0,0,0};
	for(int i=0;i<4*nInd;i++)
	  tsum[i%4] +=counts[s][i];
	fprintf(oFileCountsBin,"%zu\t%zu\t%zu\t%zu",tsum[0],tsum[1],tsum[2],tsum[3]);
      }
      fprintf(oFileCountsBin,"\n");	
    }
  }
}


void countCls::getOptions(argStruct *arguments){

  //from command line
  doCounts=angsd::getArg("-doCounts",doCounts,arguments);
  dumpCounts=angsd::getArg("-dumpCounts",dumpCounts,arguments);
  minQ=angsd::getArg("-minQ",minQ,arguments);
  trim=angsd::getArg("-trim",trim,arguments);
  doQsDist=angsd::getArg("-doQsDist",doQsDist,arguments);
  minInd = angsd::getArg("-minInd",minInd,arguments);
  setMaxDepth = angsd::getArg("-setMaxDepth",setMaxDepth,arguments);
  doDepth=angsd::getArg("-doDepth",doDepth,arguments);
  maxDepth=angsd::getArg("-maxDepth",maxDepth,arguments);
  minDepth=angsd::getArg("-minDepth",minDepth,arguments);

  if(dumpCounts&&doCounts==0){
    fprintf(stderr,"You must supply -doCounts if you want to dumpcounts\n");
    exit(0);
  }

  if(doDepth!=0&&doCounts==0){
    fprintf(stderr,"Must supply -doCounts 1 if you want depth distribution");
    exit(0);
  }

 if(doQsDist!=0&&doCounts==0){
    fprintf(stderr,"Must supply -doCounts 1 if you want qscore distribution");
    exit(0);
  }



}
//constructor
countCls::countCls(const char *outfiles,argStruct *arguments,int inputtype){
  nInd=arguments->nInd;
  minInd = 0;
  minDepth =-1;
  trim =0;
  dumpCounts =0;
  doCounts = 0;
  doQsDist = 0;
  minQ =13;
  doDepth = 0;
  maxDepth = 100;
  setMaxDepth = -1;
  
  //make output files
  postfix1=".pos";
  postfix2=".counts";
  postfix3=".qs";
  postfix4=".depthSample";
  postfix5=".depthGlobal";


  //from command line
  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-doCounts")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  getOptions(arguments);
  printArg(arguments->argumentFile);
  
  oFileCountsPos = oFileCountsBin = oFileQs = NULL;

  if(dumpCounts){
    oFileCountsPos=openFile(outfiles,postfix1);
    fprintf(oFileCountsPos,"chr\tpos\ttotDepth\n");
    if(dumpCounts>1)
      oFileCountsBin = openFile(outfiles,postfix2);
    if(dumpCounts==2)
      for(int i=0;i<arguments->nInd;i++)
	fprintf(oFileCountsBin,"ind%dTotDepth\t",i);
    if(dumpCounts==3)
      fprintf(oFileCountsBin,"totA\ttotC\ttotG\ttotT");
    if(dumpCounts==4)
      for(int i=0;i<arguments->nInd;i++)
	fprintf(oFileCountsBin,"ind%d_A\tind%d_C\tind%d_G\tind%d_T\t",i,i,i,i);
    if(dumpCounts>1)
      fprintf(oFileCountsBin,"\n");
    
  }

  if(doQsDist){
    //datastructures needed
    qsDist = new size_t[256];
    memset(qsDist,0,256*sizeof(size_t));
    //prepare outputfile
    oFileQs = openFile(outfiles,postfix3);
    fprintf(oFileQs,"qscore\tcounts\n");
  }
  if(doDepth){
    depthCount=new size_t *[arguments->nInd];
    for(int i=0;i<nInd;i++)
      depthCount[i]=new size_t[maxDepth+1];
    for(int i=0;i<nInd;i++)
      for(int j=0;j<maxDepth+1;j++)
	depthCount[i][j]=0;
    
    globCount = new size_t[maxDepth+1];
    memset(globCount,0,sizeof(size_t)*(maxDepth+1));

    oFileSamplDepth = openFile(outfiles,postfix4);
    oFileGlobDepth = openFile(outfiles,postfix5);
  }
  
}


void printQs(FILE *fp,size_t *ary){
  int firstidx=0;
  for(int i=0;i<256;i++)
     if(ary[i]!=0){
      firstidx=i;
      break;
    }
  int lastidx=255;
  for(int i=255;i>=0;i--)
    if(ary[i]!=0){
      lastidx=i;
      break;
    }
  for(int i=firstidx;i<=lastidx;i++)
    fprintf(fp,"%d\t%lu\n",i,ary[i]);
  
}


countCls::~countCls(){
  if(oFileCountsBin)
    fclose(oFileCountsBin);
  if(oFileCountsPos)
    fclose(oFileCountsPos);
  if(doQsDist){
    printQs(oFileQs,qsDist);
    fclose(oFileQs);
    delete[] qsDist;
  }

  if(doDepth){
    for(int i=0;i<nInd;i++){
      for(int j=0;j<maxDepth+1;j++){
	fprintf(oFileSamplDepth,"%lu\t",depthCount[i][j]);
      }
      fprintf(oFileSamplDepth,"\n");
    }
    //thorfinn
    for(int j=0;j<maxDepth+1;j++)
      fprintf(oFileGlobDepth,"%lu\t",globCount[j]);
    fprintf(oFileGlobDepth,"\n");
  

    //clean depthCount
    for(int i=0;i<nInd;i++)
      delete[]  depthCount[i];
    delete[] depthCount; 
    
    fclose(oFileSamplDepth);
    fclose(oFileGlobDepth);
  }


}

void countQs(const chunkyT *chk,size_t *ret,int minQ,int trim,int *keepSites){
  
  suint **cnts = new suint*[chk->nSites];
  for(int s=0;s<chk->nSites;s++){
    if(keepSites[s]==0)
      continue;
    //loop over sites
    for(int n=0;n<chk->nSamples;n++){
      //loop over samples
      for(int l=0;l<chk->nd[s][n].l;l++){
	//loop over persample reads for this position/sample
	if(chk->nd[s][n].qs[l] <minQ||chk->nd[s][n].posi[l]<trim||chk->nd[s][n].isop[l]<trim)
	  continue;
	ret[chk->nd[s][n].qs[l]]++;
	
      }
    }
  }
}



void countCls::print(funkyPars *pars){

  if(dumpCounts)
    printCounts(header->name[pars->refId],pars->posi,pars->counts,pars->numSites,pars->nInd,oFileCountsPos,oFileCountsBin,dumpCounts,pars->keepSites);

  if(doQsDist)
    countQs(pars->chk,qsDist,minQ,trim,pars->keepSites);
  
  if(doDepth!=0){
    assert(pars->counts!=NULL);
    for(int s=0;s<pars->numSites;s++) {
      if(pars->keepSites[s]==0)
	continue; 
      for(int i=0;i<pars->nInd;i++){
	int sum=0;
	for(int a=0;a<4;a++)
	  sum+=pars->counts[s][i*4+a];
	if(sum>maxDepth){
	  sum=maxDepth;
	}
	depthCount[i][sum]++;	
      }
    }
    //thorfinn below
    
    for(int s=0;s<pars->numSites;s++) {
      if(pars->keepSites[s]==0)
	continue; 
      int sum=0;
      for(int i=0;i<4*pars->nInd;i++){
	sum+=pars->counts[s][i];
	if(sum>maxDepth){
	  sum=maxDepth;
	}
      }
      globCount[sum]++;
    } 
  }

}


void countCls::clean(funkyPars *pars){
  if(doCounts||dumpCounts){
    for(int i=0;i<pars->numSites;i++)
      delete [] pars->counts[i];
    delete [] pars->counts;
  }
}


//dragon update with keeplist so we only count necessary sites
suint **countNucs(const chunkyT *chk,int minQ,int trim,int *keepSites){
  suint **cnts = new suint*[chk->nSites];
  for(int s=0;s<chk->nSites;s++){
    cnts[s] = new suint[4*chk->nSamples];
    if(keepSites[s]==0)
      continue;
    memset(cnts[s],0,4*chk->nSamples*sizeof(suint));
    //loop over samples
    for(int n=0;n<chk->nSamples;n++){
      //loop over persample reads
      for(int l=0;l<chk->nd[s][n].l;l++){
	int allele = refToInt[chk->nd[s][n].seq[l]];
	if(chk->nd[s][n].qs[l] <minQ||chk->nd[s][n].posi[l]<trim||chk->nd[s][n].isop[l]<trim||allele==4){
	  //fprintf(stderr,"ind %d,allele %d\tminQ %d\tQ %d\tMapQ %d\n",n,allele,minQ,chk->nd[s][n].mapQ[l]);
	  continue;
	}
	cnts[s][4*n+allele]++;
      }
    }
  }
  return cnts;
}




void countCls::run(funkyPars *pars){
  if(doCounts==0)
    return;
  assert(pars->chk!=NULL&&pars->counts==NULL);
  pars->counts = countNucs(pars->chk,minQ,trim,pars->keepSites);
  // fprintf(stderr,"%d\n",pars->keepSites[0]);
  //modify keepsites;
  if(minInd!=0) {
    for(int i=0;i<pars->numSites;i++){
      if(pars->keepSites[i]==0)
	continue;
      //THIS IS REALLY STUPID but lets count number of samples wiht info
      int nDep =0;
      for(int s=0;s<pars->nInd;s++){
	int dep=0;
	for(int j=0;j<4;j++)
	  dep += pars->counts[i][s*4+j];
	if(dep)
	  nDep++;
      }
      //nDep is now the number of sapmles wiht info
      if(nDep<minInd)
	pars->keepSites[i] = 0;
      else
	pars->keepSites[i] = nDep;
      
    }
  }
  if(setMaxDepth!=-1){
    for(int s=0;s<pars->numSites;s++){
      size_t totSum = calcSum(pars->counts[s],4*nInd);
      if(totSum>setMaxDepth)
	pars->keepSites[s]=0;
      else{
	suint *ps = pars->counts[s];
	for(int i=0;i<pars->nInd;i++){
	  int iSum = ps[i*4]+ps[i*4+1]+ps[i*4+2]+ps[i*4+2]+ps[i*4+3];
	  if(totSum>iSum){
	    pars->keepSites[s]=0;
	    break;
	  }
	}
      }
    }
  }
  //fprintf(stderr,"minDepth=%d nInd=%d pars->numsites=%d\n",minDepth,nInd,pars->numSites);
  if(minDepth!=-1){
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
      size_t totSum = calcSum(pars->counts[s],4*pars->nInd);
      if(totSum<minDepth)
	pars->keepSites[s]=0;
      else{
	int nDep =0;
	for(int i=0;i<pars->nInd;i++){
	  int dep=0;
	  for(int j=0;j<4;j++)
	    dep += pars->counts[s][i*4+j];
	  if(dep)
	    nDep++;
	}
	pars->keepSites[s]= nDep;

      }
    }
  }
  
}
