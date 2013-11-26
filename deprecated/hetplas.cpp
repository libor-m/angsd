#include "analysisFunction.h"
#include "shared.h"
#include <ctype.h>
#define MAX_QS 500


static size_t fullDist[MAX_QS][MAX_QS][2];//<-strand
static size_t majDist[MAX_QS][MAX_QS][2];//<-strand
static size_t minDist[MAX_QS][MAX_QS][2];//<-strand

//struct which contains the results for a single chunk
typedef struct{
  double **freq;//contains the estimate of the freq of the 4 alleles and the unknown minor freq;
  double **llh;//contains the LRT and the llhNull llhAlt;
  int **nItr;
  int **maxbases;//contains the allele of the maxbase
  double **diff;
  size_t **qsMajor;
  size_t **qsMinor;
  size_t **majPosi;
  size_t **minPosi;
}resStruct;


double **gen_probs(){
  double **ret =new double*[MAX_QS];
  for(int i=0;i<MAX_QS;i++)
    ret[i] = new double[16];

  

  for(int i=0;i<MAX_QS;i++){
    double p1 = log(1-pow(10,-i/10.0));
    double p3 = log((1-exp(p1))/3.0);
    double val= std::max(p1,p3);
    val=std::min(p1-val,p3-val);
    for(int ii=0;ii<16;ii++){
      ret[i][ii] = exp(val);
    }
    ret[i][0]=ret[i][5]=ret[i][10]=ret[i][15] =exp(0);
  }
  
  return ret;
}


class hetplas:public general{
private:
  const char *outfile_c;
  double **probs;
  FILE *fp;
  FILE *fptest;
  FILE *fptest2;
  FILE *fptestpo;
  FILE *fptestop;
  FILE *fpwhichmax;
  int minQ;
public:
  void calcQsDists(int,size_t*,size_t*,size_t *,size_t*,tNode *);
  int doHetPlas;
  int makellhs(tNode*,double**,int*);
  void doNew(funkyPars *pars);
  hetplas(const char *outfiles,argStruct *arguments,int inputtype);
  ~hetplas();
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  
  
};
void hetplas::printArg(FILE *fp){
  fprintf(fp,"--------------\n%s:\n",__FILE__);
  fprintf(fp,"\t-doHetPlas=%d\n",doHetPlas);
  fprintf(fp,"\t-minQ=%d\n",minQ);
}

int hetplas::makellhs(tNode *tn,double **liks,int *oldsize){
  if(*oldsize<tn->l){
    if(*oldsize>0){
      for(int i=0;i<4;i++)
	delete [] liks[i];
    }
    for(int i=0;i<4;i++)
      liks[i] = new double[tn->l];
    *oldsize = tn->l;
  }
  //  fprintf(stderr,"tn->l=%d\n",tn->l);
  int depth =0;
  for(int l=0;l<tn->l;l++){
    //fprintf(stderr,"seq=%c int=%d qs=%d\n",tn->seq[l],refToInt[tn->seq[l]],tn->qs[l]);
    if(refToInt[tn->seq[l]]!=4&&tn->qs[l]>=minQ){
      for(int ll=0;ll<4;ll++){
     	liks[ll][depth] = probs[tn->qs[l]][refToInt[tn->seq[l]]*4+ll];
      }
      depth++;
    }//else
     // fprintf(stderr,"skipping\n");
  }
  return depth;
}

//program will calculate the likelihood of the x' given the liks, seqdepth is the dimension
double like(double *x,double **liks,int seqdepth){
  double totLik=0;
  for(int s=0;s<seqdepth;s++){
    double tmp=0;
    for(int i=0;i<4;i++)
      tmp+=x[i]*liks[i][s];
    totLik += log(tmp);
  }
  return totLik;
}

// function will perform em estimation of the frequencines
//function will modfy x,itr,and diff

double em3(double *x,double **liks,int seqdepth,int maxIter,double tole,int &itr,double &diff){
  itr=0;
  double likev=-1e12;
  
  while(itr<maxIter){
    double newx[4]={0,0,0,0};
    memcpy(newx,x,sizeof(double)*4);
    for(int s=0;s<seqdepth;s++){
      double tsum=0;
      for(int i=0;i<4;i++)
	tsum += x[i]*liks[i][s];
      for(int i=0;i<4;i++)
	newx[i] += x[i]*liks[i][s]/tsum;
    }
    double tsum =0;
    for(int i=0;i<4;i++)
      tsum +=newx[i];
    for(int i=0;i<4;i++)
      newx[i] /=tsum;
	
    double newlike=like(newx,liks,seqdepth);
    diff=fabs(newlike-likev);
    if(diff<tole){
      likev=newlike;
      memcpy(x,newx,sizeof(double)*4);
      //      fprintf(stderr,"breaking\n");
      break;
    }

    if(1&&(newlike<likev)){
      fprintf(stderr,"Problems: newlike=%f oldlike=%f diff=%e itr=%d\n",newlike,likev,diff,itr);
      break;
      exit(0);
    }
    assert(newlike>=likev);
    //fprintf(stderr,"lik=%f newlike=%f\n",likev,newlike);

    likev=newlike;
    memcpy(x,newx,sizeof(double)*4);
    itr++;
  }
  //fprintf(stderr,"itr=%d like=%f par=(%f,%f,%f,%f)\n",itr,likev,x[0],x[1],x[2],x[3]);

}

void hetplas::calcQsDists(int major,size_t *qsMajor,size_t *qsMinor,size_t *majPosi,size_t*minPosi,tNode *tn){

  for(int l=0;l<tn->l;l++){
    int oBase = refToInt[tn->seq[l]];
    if(refToInt[tn->seq[l]]!=4&&tn->qs[l]>=minQ){
      if(oBase==major){
	qsMajor[tn->qs[l]]++;
	//fprintf(stderr,"%d\n",tn->posi[l]);
	majPosi[tn->posi[l]]++;
      }else{
	qsMinor[tn->qs[l]]++;
	minPosi[tn->posi[l]]++;
      }
    }
  }
}



void hetplas::doNew(funkyPars *pars){
  //fprintf(stderr,"pars->numSites=%d\n",pars->numSites);
  //  fflush(stderr);
  resStruct *rs = new resStruct;
  rs->freq=  new double*[pars->numSites]; 
  rs->llh= new double*[pars->numSites]; 
  rs->nItr = new int*[pars->numSites];
  rs->diff = new double *[pars->numSites];
  rs->qsMajor = new size_t *[pars->numSites];
  rs->qsMinor = new size_t *[pars->numSites];
  rs->majPosi = new size_t *[pars->numSites];
  rs->minPosi = new size_t *[pars->numSites];
  rs->maxbases=new int*[pars->numSites];

  double **liks = new double*[4];
  int oldsize =0;
  for(int s=0;s<pars->numSites;s++){
    //    fprintf(stderr,"AAAA %d %d\n",pars->posi[s]+1,pars->keepSites[s]);
    if(pars->keepSites[s]) {
      rs->freq[s] = new double[5*pars->nInd];  
      rs->llh[s] = new double [3*pars->nInd];  
      rs->diff[s] = new double [pars->nInd];  
      rs->nItr[s] = new int [pars->nInd];  
      rs->qsMajor[s] = new size_t[MAX_QS];
      memset(rs->qsMajor[s],0,MAX_QS*sizeof(size_t));
      rs->qsMinor[s] = new size_t[MAX_QS];
      memset(rs->qsMinor[s],0,MAX_QS*sizeof(size_t));
      rs->majPosi[s]=new size_t[MAX_QS];
      memset(rs->majPosi[s],0,MAX_QS*sizeof(size_t));
      rs->minPosi[s]=new size_t[MAX_QS];
      memset(rs->minPosi[s],0,MAX_QS*sizeof(size_t));
      rs->maxbases[s] = new int[pars->nInd];
      memset(rs->maxbases[s],0,sizeof(int)*pars->nInd);      
      for(int i=0;i<pars->nInd;i++){
	int seqdepth=makellhs(&pars->chk->nd[s][i],liks,&oldsize);
	if(seqdepth==0){
	  pars->keepSites[s]=0;
	  continue;//DRAGON should only skip ind, is ok with 1 sample
	}
	double par[4]={0.25,0.25,0.25,0.25};
	//¯	fprintf(stderr,"lik=%f\n",like(par,liks,seqdepth));
	em3(par,liks,seqdepth,100,1e-6,rs->nItr[s][i],rs->diff[s][i]);
	int maxBase=angsd::whichMax(par,4);
	rs->maxbases[s][i]=maxBase;
	//assert(maxBase!=-1);//<- if vals are equal...
	if(maxBase==-1){
	  fprintf(stderr,"CHECK SITE: %d  par=(%f,%f,%f,%f) seqdepth=%d\n",pars->posi[s]+1,par[0],par[1],par[2],par[3],seqdepth);
	  pars->keepSites[s] =0;
	  continue;
	}
	double unknownMinorSum =0;
	for(int um=0;um<4;um++)
	  if(um!=maxBase)
	    unknownMinorSum += par[um];
	memcpy(rs->freq[s]+4*i,par,sizeof(double)*4);
	rs->freq[s][4*i+4]=unknownMinorSum;
	double llhAlt = like(par,liks,seqdepth);
	par[0]=par[1]=par[2]=par[3]=0.0;
	par[maxBase]=1.0;
	double llhNeu=like(par,liks,seqdepth);
	rs->llh[s][i*3]=2*llhNeu-2*llhAlt;
	rs->llh[s][i*3+1]=llhNeu;
	rs->llh[s][i*3+2]=llhAlt;
	calcQsDists(maxBase,rs->qsMajor[s],rs->qsMinor[s],rs->majPosi[s],rs->minPosi[s],&pars->chk->nd[s][i]);
      }
      
    }
  }
  for(int i=0;i<4;i++)
    delete [] liks[i];
  delete [] liks;
  pars->extras[index] = rs;
}

void hetplas::run(funkyPars *pars){
  //  fprintf(stderr,"pars->numSites=%d\n",pars->numSites);
  //exit(0);
  if(doHetPlas==0)
    return ;
  if(!doHetPlas||pars->numSites==0){
    fprintf(stderr,"skipping due to nos sites\n");
    return;
  }
  assert(pars->numSites>0);
  assert(pars->chk!=NULL);
  if(doHetPlas){
    //exit(0);
    doNew(pars);
  }
  
   
}

void hetplas::clean(funkyPars *pars){
  if(doHetPlas==0)
    return;
  if(doHetPlas){
    resStruct *rs = (resStruct *)pars->extras[index];
    for(int i=0;i<pars->numSites;i++)
      if(pars->keepSites[i]){
	delete [] rs->freq[i];
	delete [] rs->llh[i];
	delete [] rs->diff[i];
	delete [] rs->nItr[i];
      }
    delete [] rs->freq;
    delete [] rs->llh;
    delete [] rs->nItr;
    delete [] rs->diff;
    delete rs;
  }
    

}

void hetplas::print(funkyPars *pars){
  if(doHetPlas==0)
    return;
  if(doHetPlas){
    resStruct *rs =(resStruct *) pars->extras[index];
    for(int s=0;s<pars->numSites;s++)
      if(pars->keepSites[s]){
	fprintf(fp,"%s\t%d",header->name[pars->refId],pars->posi[s]+1);
	for(int i=0;i<pars->nInd;i++){
	  //write freq
	  for(int j=0;j<5;j++)
	    fprintf(fp,"\t%f",rs->freq[s][i*4+j]);
	  for(int j=0;j<3;j++)
	    fprintf(fp,"\t%f",rs->llh[s][i*3+j]);
	  fprintf(fp,"\t%d\t%e",rs->nItr[s][i],rs->diff[s][i]);
	}
	fprintf(fp,"\n");
      }
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]){
	for(int q=0;q<MAX_QS;q++)
	  fprintf(fptest,"%zu\t",rs->qsMajor[s][q]);
	fprintf(fptest,"\n");
      }
    }
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]){
	for(int q=0;q<MAX_QS;q++)
	  fprintf(fptest2,"%zu\t",rs->qsMinor[s][q]);
	fprintf(fptest2,"\n");
      }
    }
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]){
	for(int q=0;q<MAX_QS;q++)
	  fprintf(fptestpo,"%zu\t",rs->majPosi[s][q]);
	fprintf(fptestpo,"\n");
      }
    }
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]){
	for(int q=0;q<MAX_QS;q++)
	  fprintf(fptestop,"%zu\t",rs->minPosi[s][q]);
	fprintf(fptestop,"\n");
      }
    }
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
      for(int i=0;i<pars->nInd;i++) {
	tNode *tn = &pars->chk->nd[s][i];
	for(int l=0;l<tn->l;l++){
	  if(tn->qs[l]<minQ||refToInt[tn->seq[l]]==4)
	    continue;
	  int b1 = refToInt[tn->seq[l]];
	  int st = !isupper(tn->seq[l]);
	  fullDist[tn->posi[l]][tn->qs[l]][st]++;
	  if(b1==rs->maxbases[s][0])
	    majDist[tn->posi[l]][tn->qs[l]][st]++;
	  else
	    minDist[tn->posi[l]][tn->qs[l]][st]++;
	}

      }
    }
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0)
	continue;
      fprintf(fpwhichmax,"%d\t%d\n",pars->posi[s]+1,rs->maxbases[s][0]);
    }
  }
}

void hetplas::getOptions(argStruct *arguments){
  //from command line
  

  doHetPlas=angsd::getArg("-doHetPlas",doHetPlas,arguments);
 
  if(doHetPlas==0)
    return;

  minQ = angsd::getArg("-minQ",minQ,arguments);
  
  printArg(arguments->argumentFile);

}


hetplas::hetplas(const char *outfiles,argStruct *arguments,int inputtype){
  doHetPlas =0;
  minQ = 13;
  probs=NULL;
  outfile_c = outfiles;
  
  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-doHetPlas")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }


  getOptions(arguments);
  if(doHetPlas)
    fprintf(stderr,"running doHetPlas=%d\n",doHetPlas);
  if(doHetPlas){
    probs=gen_probs();
    fp=openFile(outfiles,".hetGL");
    fptest=openFile(outfiles,".hetTest");
    fptest2=openFile(outfiles,".hetTest2");
    fptestpo=openFile(outfiles,".posi");
    fptestop=openFile(outfiles,".isop");
    fpwhichmax=openFile(outfiles,".whichmax");
  } 
}


void printStuff(FILE *fp,size_t mat[MAX_QS][MAX_QS][2]){
  for(int st=0;st<2;st++){
    for(int i=0;i<MAX_QS;i++){
      for(int j=0;j<MAX_QS;j++)
	fprintf(fp,"%zu\t",mat[i][j][st]);
      fprintf(fp,"%d\n",st);
    }
  }

}


hetplas::~hetplas(){
  if(doHetPlas){
    fclose(fp);
    fclose(fptest);
    fclose(fptest2);
    fclose(fptestpo);
    fclose(fptestop);
    FILE *myfp=openFile(outfile_c,".fullDist");
    printStuff(myfp,fullDist);
    fclose(myfp);
    myfp=openFile(outfile_c,".majDist");
    printStuff(myfp,majDist);
    fclose(myfp);
    myfp=openFile(outfile_c,".minDist");
    printStuff(myfp,minDist);
    fclose(myfp);
  }
  
}


