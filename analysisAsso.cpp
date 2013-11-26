/*
  thorfinn thorfinn@binf.ku.dk dec17 2012
  
  has modified so no sites[].chromo etc are used
    
  anders albrecht@binf.ku.dk made this.

  part of angsd

  DRAGON positions are not offset correctly
*/
#include <cmath>
#include <zlib.h>
#include "kstring.h"
#include "shared.h"
#include "analysisFunction.h"


class asso:public general{
private:
  kstring_t bufstr;
public:
  //none optional stuff
  gzFile *MultiOutfile;
  int doPrint;
  int minCov; //not for users
  int doMaf;
  int dynCov;//not for users
  int doAsso;
  int doPost;
  int GL;
  int sitePerm;  //not for users
  int isBinary;
  int minHigh;
  int minCount;
  int adjust;  //not for users
  int model;
  void run(funkyPars  *pars);
  void print(funkyPars *pars);  
  void clean(funkyPars *pars);  
  void getOptions(argStruct *arguments);
  void printArg(FILE *argFile);

  asso(const char *outfiles,argStruct *arguments,int inputtype);
  ~asso();
  //other stuff
  char *covfile;
  char *yfile;
  
  angsd::Matrix<double> ymat;
  angsd::Matrix<double> covmat;
  void scoreAsso(funkyPars  *pars);
  void frequencyAsso(funkyPars  *pars);
  double doAssociation(funkyPars *pars,double *post,double *y,int keepInd,int *keepList,double freq,int s);
  void getFit(double *res,double *Y,double *covMatrix,int nInd,int nEnv);
  double normScoreEnv(double *post,int numInds, double *y, double *ytilde,double *cov,int nEnv,double freq,funkyPars *pars,int s);
  double binomScoreEnv(double *post,int numInds, double *y, double *ytilde,double *cov,int nEnv,double freq,funkyPars *pars,int s);
  void printDoAsso(funkyPars *pars);
};

void asso::printArg(FILE *argFile){
  fprintf(argFile,"-------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-doAsso\t%d\n",doAsso);
  fprintf(argFile,"\t1: Frequency Test (Known Major and Minor)\n");
  fprintf(argFile,"\t2: Score Test\n");
  fprintf(argFile,"\t3: Frequency Test (Unknown Minor)\t\n");
  fprintf(argFile,"  Frequency Test Options:\n");
  fprintf(argFile,"\t-yBin\t\t%s\t(File containing disease status)\t\n\n",yfile);
  fprintf(argFile,"  Score Test Options:\n");
  fprintf(argFile,"\t-yBin\t\t%s\t(File containing disease status)\n",yfile);
  fprintf(argFile,"\t-yQuant\t\t%s\t(File containing phenotypes)\n",yfile);
  fprintf(argFile,"\t-minHigh\t%d\t(Require atleast minHigh number of high credible genotypes)\n",minHigh);
  fprintf(argFile,"\t-minCount\t%d\t(Require this number of minor alleles, estimated from MAF)\n",minCount);
  fprintf(argFile,"\t-cov\t\t%s\t(File containing additional covariates)\n",covfile);
  fprintf(argFile,"\t-model\t%d\n",model);
  fprintf(argFile,"\t1: Additive/Log-Additive (Default)\n");
  fprintf(argFile,"\t2: Dominant\n");
  fprintf(argFile,"\t3: Recessive\n\n");
  fprintf(argFile,"Examples:\n\tPerform Frequency Test\n\t  \'./angsd -yBin pheno.ybin -doAsso 1 -GL 1 -out out -doMajorMinor 1 -minLRT 24 -doMaf 2 -doSNP 1 -bam bam.filelist'\n");
  fprintf(argFile,"\tPerform Score Test\n\t  \'./angsd -yBin pheno.ybin -doAsso 2 -GL 1 -doPost 1 -out out -doMajorMinor 1 -minLRT 24 -doMaf 2 -doSNP 1 -bam bam.filelist'\n");
  fprintf(argFile,"\n");
}


void asso::getOptions(argStruct *arguments){


  doAsso=angsd::getArg("-doAsso",doAsso,arguments);

  doMaf=angsd::getArg("-doMaf",doMaf,arguments);

  adjust=angsd::getArg("-adjust",adjust,arguments);
  model=angsd::getArg("-model",model,arguments);
  minCov=angsd::getArg("-minCov",minCov,arguments);
  dynCov=angsd::getArg("-dynCov",dynCov,arguments);
  minHigh=angsd::getArg("-minHigh",minHigh,arguments);
  doPrint=angsd::getArg("-doPrint",doPrint,arguments);
  minCount=angsd::getArg("-minCount",minCount,arguments);
  sitePerm=angsd::getArg("-sitePerm",sitePerm,arguments);
  GL=angsd::getArg("-GL",GL,arguments);
  covfile=angsd::getArg("-cov",covfile,arguments);
  doPost=angsd::getArg("-doPost",doPost,arguments);
  yfile=angsd::getArg("-yBin",yfile,arguments);
  if(yfile!=NULL)
    isBinary=1;
  yfile=angsd::getArg("-yQuant",yfile,arguments);

  if(doPrint)
    fprintf(stderr,"finished [%s]\t[%s]\n",__FILE__,__FUNCTION__);


  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  if(doAsso && doMaf==0){
    fprintf(stderr,"Error: you must estimate the maf (-doMaf) in order to perform association \n");
    exit(0);

  }

  if(doAsso && yfile==NULL){
    fprintf(stderr,"Error: you must provide a phenotype file (-yBin or -yQuant) to perform association \n");
    exit(0);
   }
 if(doAsso && arguments->inputtype==5&&doAsso==1){
    fprintf(stderr,"Error: Only doAsso=2 can be performed on posterior input\n");
    exit(0);
  }
  if(doAsso && arguments->inputtype!=5&&(doAsso==2)&&doPost==0){
    fprintf(stderr,"Error: For doAsso=2 you must estimate the posterior probabilites for the genotypes (doPost!=0) \n");
    exit(0);
  }  

}
 
asso::asso(const char *outfiles,argStruct *arguments,int inputtype){
  bufstr.s=NULL;bufstr.l=bufstr.m=0;
  //default
  model=1;
  doPrint=0;
  doAsso=0;
  GL=0;
  doPost=0;
  isBinary=0;
  sitePerm=0;//not for users
  covfile=NULL;
  yfile=NULL;
  minHigh=10;
  minCount=10;
  dynCov=0;//not for users
  minCov=5;//not for users
  adjust=1;//not for users
  doMaf=0;
  //from command line


  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-doAsso")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }
  

  getOptions(arguments);
  printArg(arguments->argumentFile);

  if(doAsso==0)
    return;


  //read phenotype
  if(isBinary)
    ymat = angsd::getMatrix(yfile,1,100000);
  else
    ymat = angsd::getMatrix(yfile,0,100000);

  //read covariates 
  if(covfile!=NULL)
    covmat = angsd::getMatrix(covfile,0,100000);
  else{
    covmat.x=0;
    covmat.y=0;
    covmat.matrix=NULL;
  }
  if(covfile!=NULL&&(covmat.x!=ymat.x)){
    fprintf(stderr,"The number of covariates (%d) does not match the number of phenotypes (%d)\n",covmat.x,ymat.x);
    exit(0);
  }

  //inputtype 0) soap 1) samglf 2) samglfClean 3) tglf 4)sim type 5) beagle
 
  if(!isBinary&&doAsso==1){
    fprintf(stderr,"Error: Only doAsso=2 can be performed on quantitative traits\n");
    exit(0);
  }

  //make output files
  MultiOutfile = new gzFile[ymat.y];
  const char* postfix;
  postfix=".lrt";

  if(covfile!=NULL&&minCov>0){
    int keepList[ymat.x];
    for(int i=0 ; i < ymat.x;i++) {
      keepList[i]=1;
      for(int yi=0;yi<ymat.y;yi++) {
	if(ymat.matrix[i][yi]==-999)
	  keepList[i]=0;
      }
      for(int ci=0;ci<covmat.y;ci++) {
	if(covmat.matrix[i][ci]==-999)
	  keepList[i]=0;
      }
    }
    int nCov=0;
    int count[covmat.y];
    for(int ci=0;ci<covmat.y;ci++) {
      count[ci]=0;
      for(int i=0 ; i < ymat.x;i++) {
	if(keepList[i]==0)
	  continue;
	if(covmat.matrix[i][ci]!=0){
	  count[ci]++;
	}
      }
  
      if(count[ci]<minCov){
	fprintf(stderr,"Error: Cov #%d only has %d non zero entries\n",ci,count[ci]);
      }
      else
	nCov++;

    }
    if(!dynCov&&covmat.y!=nCov){
      fprintf(stderr,"Error: Creating new covariant matrix with %d columns\n",nCov);
      exit(0);

    }
    else if(covmat.y!=nCov){
      //      angsd::printMatrix(covmat,stderr);
      fprintf(stderr,"Error: Creating new covariant matrix with %d columns\n",nCov);
      angsd::Matrix<double> newmat;
      newmat.x=covmat.x;
      newmat.y=nCov;
      newmat.matrix=new double*[covmat.x];
      for(int xi=0;xi<covmat.x;xi++){
	newmat.matrix[xi] = new double[nCov];
	int tempCount=0;
	for(int ci=0;ci<covmat.y;ci++){
	  if(count[ci]>minCov){
	    newmat.matrix[xi][tempCount]=covmat.matrix[xi][ci];
	    tempCount++;
	  }
	}
      }
      angsd::deleteMatrix(covmat);
      covmat=newmat;
      
    }
  }

  //open outfiles
  for(int i=0;i<ymat.y;i++){
    char ary[5000];
    snprintf(ary,5000,"%s%d.gz",postfix,i);
    MultiOutfile[i] = openFileGz(outfiles,ary,GZOPT);
  }

  //print header
  for(int yi=0;yi<ymat.y;yi++){
    if(doAsso==2)
      gzprintf(MultiOutfile[yi],"Chromosome\tPosition\tMajor\tMinor\tFrequency\tN\tLRT\thighHe\thighHo\n");
    else
      gzprintf(MultiOutfile[yi],"Chromosome\tPosition\tMajor\tMinor\tFrequency\tLRT\n");
  }
}


asso::~asso(){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);


  if(doAsso==0)
    return;
  for(int i=0;i<ymat.y;i++)
    gzclose(MultiOutfile[i]);
  delete [] MultiOutfile;

  if(covfile!=NULL)
    angsd::deleteMatrix(covmat);
  angsd::deleteMatrix(ymat);

}


void asso::clean(funkyPars *pars){
  if(doAsso==0)
    return;


  for(int yi=0;yi<ymat.y;yi++)
     delete[] pars->results->asso->stat[yi];

  delete[] pars->results->asso->stat;
  

  if(doAsso==2){
    delete[]  pars->results->asso->highHe;
    delete[]  pars->results->asso->highHo;
    
    for( int yi =0;yi<ymat.y;yi++)
      delete[] pars->results->asso->keepInd[yi];
    delete[] pars->results->asso->keepInd;
    
  }

}



void asso::print(funkyPars *pars){
  if(doAsso==0)
    return;

  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  printDoAsso(pars);
}


void asso::run(funkyPars *pars){


  if(doAsso==0)
    return;
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

 
  if(doAsso==1||doAsso==3){
    frequencyAsso(pars);
  }
  else if(doAsso==2){
    pars->results->asso->highHe=new int[pars->numSites];
    pars->results->asso->highHo=new int[pars->numSites];

    scoreAsso(pars);
  }

}

void asso::frequencyAsso(funkyPars  *pars){

  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  if(pars->nInd!=ymat.x){
    fprintf(stderr,"The number of sequenced individuals (%d) does not match the number of phenotypes (%d)\n",pars->nInd,ymat.x);

    fflush(stderr);
    exit(0);
  }

  if(ymat.y!=1){
    fprintf(stderr,"Only one phenotype allowed for doAsso 1 or 3 \n");
    fflush(stderr);
    exit(0);
  }

  double **stat = new double*[ymat.y];
  for(int yi=0;yi<ymat.y;yi++)
    stat[yi] = new double[pars->numSites];
 
 
  
  int y[pars->nInd];
  for(int i=0;i<pars->nInd;i++)
    y[i]=ymat.matrix[i][0];

  double **like0;//genotype likelihood for controls
  double **like1;//genotype likelihood for cases
  double **likeAll;//genotype likelihood for cases and controls
  int Ncases=0; //number of cases
  int Ncontrols=0; //number of cases
  int Nall=0; //number of cases and controls
  int cases[pars->nInd];
  int controls[pars->nInd];
  int all[pars->nInd];

  for(int i=0;i<pars->nInd;i++){
    cases[i]=0;
    controls[i]=0;
    all[i]=0;
    if((int)y[i]==1){
      Ncases++;
      cases[i]=1;
      Nall++;
      all[i]=1;
    }
    if((int)y[i]==0){
      Ncontrols++;
      controls[i]=1;
      Nall++;
      all[i]=1;
    }
    if(doPrint)
      fprintf(stderr,"all, case, control: %d %d %d\n",all[i],cases[i],controls[i]);
  }
  if(doPrint)
    fprintf(stderr,"count complete [%s]\t[%s]\n",__FILE__,__FUNCTION__);


  if(doAsso==1){
    like0=angsd::get3likes(pars,controls);
    like1=angsd::get3likes(pars,cases);
    likeAll=angsd::get3likes(pars,all);
  }
  if(doAsso==3){//use all 10 genotype likes
    like0=angsd::getlikes(pars,controls);
    like1=angsd::getlikes(pars,cases);
    likeAll=angsd::getlikes(pars,all);
  }
 if(doPrint)
    fprintf(stderr,"like complete [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  for(int s=0;s<pars->numSites;s++){//loop overs sites
    stat[0][s]=-999;
    if(pars->keepSites[s]==0)
      continue;
    if(doAsso==1){

      if(doPrint)
	fprintf(stderr,"do freq [%s]\t[%s]\n",__FILE__,__FUNCTION__);


      double score0=frequency::likeFixedMinor(frequency::likeFixedMinor_bfgs(like0[s],Ncontrols),like0[s],Ncontrols);
      //likelihood for the cases
      double score1=frequency::likeFixedMinor(frequency::likeFixedMinor_bfgs(like1[s],Ncases),like1[s],Ncases);
      //likelihood for all individuals
      double scoreNull=frequency::likeFixedMinor(frequency::likeFixedMinor_bfgs(likeAll[s],Nall),likeAll[s],Nall);
      //likelhood ratio statistics \sim chi^2
      double LRT=-2*(score0+score1-scoreNull);
      stat[0][s]=LRT;
    }
    if(doAsso==3){
      double score0=frequency::likeNoFixedMinor(frequency::likeNoFixedMinor_bfgs(like0[s],Ncontrols,pars->major[s]),like0[s],Ncontrols,pars->major[s]);
      //likelihood for the cases
      double score1=frequency::likeNoFixedMinor(frequency::likeNoFixedMinor_bfgs(like1[s],Ncases,pars->major[s]),like1[s],Ncases,pars->major[s]);
      //likelihood for all individuals
      double scoreNull=frequency::likeNoFixedMinor(frequency::likeNoFixedMinor_bfgs(likeAll[s],Nall,pars->major[s]),likeAll[s],Nall,pars->major[s]);
      //likelhood ratio statistics \sim chi^2
      double LRT=-2*(score0+score1-scoreNull);
      stat[0][s]=LRT;
    }
   

  }
  pars->results->asso->stat=stat;
  for(int s=0;s<pars->numSites;s++){
    delete[] like0[s];
    delete[] like1[s];
    delete[] likeAll[s];
  }
   delete[] like0;
   delete[] like1;
   delete[] likeAll;


 if(doPrint)
    fprintf(stderr,"finish [%s]\t[%s]\n",__FILE__,__FUNCTION__);

}

void asso::scoreAsso(funkyPars  *pars){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  if(pars->nInd!=ymat.x){
    fprintf(stderr,"The number of sequenced individuals (%d) does not match the number of phenotypes (%d)\n",pars->nInd,ymat.x);
    exit(0);
  }


  int **keepInd  = new int*[ymat.y];
  double **stat = new double*[ymat.y];
  for(int yi=0;yi<ymat.y;yi++){
    stat[yi] = new double[pars->numSites];
    keepInd[yi]= new int[pars->numSites];
  }
  
  for(int s=0;s<pars->numSites;s++){//loop overs sites
    if(pars->keepSites[s]==0)
      continue;
    
    
    int *keepListAll = new int[pars->nInd];
    for(int i=0 ; i<pars->nInd ;i++){
      keepListAll[i]=1;

    }

    for(int yi=0;yi<ymat.y;yi++) { //loop over phenotypes
      int *keepList = new int[pars->nInd];
      keepInd[yi][s]=0;
      for(int i=0 ; i<pars->nInd ;i++) {
	keepList[i]=1;
	if(keepListAll[i]==0||ymat.matrix[i][yi]==-999)
	  keepList[i]=0;
	if(covfile!=NULL)
	  for(int ci=0;ci<covmat.y;ci++) {
	    if(covmat.matrix[i][ci]==-999)
	      keepList[i]=0;
	  }


	if(keepList[i]==1)
	  keepInd[yi][s]++;
      }  
      double *y = new double[pars->nInd];
      for(int i=0 ; i<pars->nInd ;i++)
	y[i]=ymat.matrix[i][yi]; 
 

      stat[yi][s]=doAssociation(pars,pars->post[s],y,keepInd[yi][s],keepList,pars->results->asso->freq[s],s);
      
      //cleanup
       delete [] y;
      delete [] keepList;

    } //phenotypes end
 
    delete [] keepListAll;
  } // sites end

  pars->results->asso->stat=stat;
  pars->results->asso->keepInd=keepInd;
}


void asso::getFit(double *res,double *Y,double *covMatrix,int nInd,int nEnv){

  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  double Xt_y[nEnv];
  double invXtX_Xt_y[nEnv];
  for(int i=0;i<nEnv;i++)
    Xt_y[i]=0;
  for(int i=0;i<nEnv;i++)
    invXtX_Xt_y[i]=0;


 
  //get t(X)%*%y
  for(int x=0;x<nEnv;x++)//col X
    for(int i=0;i<nInd;i++)
      Xt_y[x]+=covMatrix[x*nInd+i]*Y[i];

  //get inv(t(X)%*%X)
  double XtX[nEnv*nEnv];
  for(int i=0;i<nEnv*nEnv;i++)
    XtX[i]=0;

  for(int x=0;x<nEnv;x++)//col X
    for(int y=0;y<nEnv;y++)//row Xt
      for(int i=0;i<nInd;i++)
	XtX[x*nEnv+y]+=covMatrix[y*nInd+i]*covMatrix[x*nInd+i];

  double workspace[2*nEnv];
  angsd::matinv(XtX, nEnv, nEnv, workspace);


  //get (inv(t(X)%*%X))%*%(t(X)%*%y)
 for(int x=0;x<nEnv;x++)//col X
   for(int y=0;y<nEnv;y++)//row Xt
      invXtX_Xt_y[x]+=XtX[y*nEnv+x]*Xt_y[y];

  //get X%*%(inv(t(X)%*%X))%*%(t(X)%*%y)
   for(int j=0;j<nInd;j++)//row Xt
     res[j]=0;
 
   for(int j=0;j<nInd;j++)//row Xt
     for(int x=0;x<nEnv;x++)
       res[j]+=covMatrix[x*nInd+j]*invXtX_Xt_y[x];
 
}



double asso::doAssociation(funkyPars *pars,double *postOrg,double *yOrg,int keepInd,int *keepList,double freq,int s){
  if(doPrint)
    fprintf(stderr,"Staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  double covMatrix[(covmat.y+1)*keepInd];
  double y[keepInd];
  double post[keepInd*3];
  int count=0;
  for(int i=0;i<pars->nInd;i++){
    if(keepList[i]){
      y[count]=yOrg[i];
      for(int g=0;g<3;g++)
	post[count*3+g]=postOrg[i*3+g];
      count++;
    }

  }

  if(count!=keepInd){
    fprintf(stderr,"[%s] wrong number of non missing\n",__FUNCTION__);
    fflush(stderr);
    exit(0);
  }

  int nEnv;
  if(adjust==1)
    nEnv=(covmat.y+1);
  else
    nEnv=1;



  int num;
  for(int j=0;j<keepInd;j++)
    covMatrix[j]=1;
  if(covmat.matrix!=NULL){
    num=0;
    for(int j=0;j<covmat.x;j++){
      if(keepList[j]==0)
	continue;
      for(int i=1;i<nEnv;i++)
	covMatrix[i*keepInd+num]=covmat.matrix[j][i-1];   
      num++;
    }
  }

  //fprintf(stderr,"number of inds %d\n",keepInd);

  // permutation
  if(sitePerm){
    if((covmat.y+1)==1){
      for(int i=0 ; i<keepInd ;i++){	
	int j = rand() % (keepInd);
	angsd::swapDouble(y[j],y[i]); 
      }
    }
    else{
      int col0=0; 
      for(int i=0 ; i<covmat.x ;i++) {
	if(keepList[i]==0)
	  continue;
	if(covmat.matrix[i][0]<0.5)
	  col0++;
      }
      for(int i=0 ; i<col0 ;i++) {
	int j = rand() % (col0);
	angsd::swapDouble(y[j],y[i]);
      }
      for(int i=0 ; i<keepInd-col0 ;i++) {
	int j = rand() % (keepInd-col0);
	angsd::swapDouble(y[j+col0],y[i+col0]);
      }
      if(col0<500||keepInd-col0<500){
	fprintf(stderr,"colTrouble %d %d\n",col0,keepInd-col0);
      }
    }

  }

  //


 
  double *yfit = new double[keepInd];
  if(nEnv==1){
    double mean=0;
    for(int i=0;i<keepInd;i++)
      mean+=y[i];
    mean=mean/keepInd;
    for(int i=0;i<keepInd;i++)
      yfit[i]=mean;
  }
  else
    getFit(yfit,y,covMatrix,keepInd,nEnv);
  
  //for(int i=0;i<keepInd;i++)
  //   fprintf(stdout,"%f\t",y[i]);
  // exit(0);


  if(model==2){
    for(int i=0 ; i<keepInd ;i++) {
      post[i*3+1]+=post[i*3+2];
      post[i*3+2]=0;
    }
  }
  if(model==3){
    for(int i=0 ; i<keepInd ;i++) {
      post[i*3+0]+=post[i*3+1];
      post[i*3+1]=post[i*3+2];
      post[i*3+2]=0;
    }
  }


  double stat;
  if(isBinary)
    stat = binomScoreEnv(post,keepInd,y,yfit,covMatrix,nEnv,freq,pars,s);
  else
    stat = normScoreEnv(post,keepInd,y,yfit,covMatrix,nEnv,freq,pars,s);

  delete[] yfit;
  return stat;

}




double asso::normScoreEnv(double *post,int numInds, double *y, double *ytilde,double *cov,int nEnv,double freq,funkyPars *pars,int s){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  int rankProb=0;
  double sum=0;
  double *Ex = angsd::allocArray<double>(numInds,0);
  double *Ex2 = angsd::allocArray<double>(numInds,0);
  double U=0;
  int highHE=0;
  int highHO=0;
  
  double sumEx=0;
  
  for(int i=0;i<numInds;i++)
    sum+=pow(y[i]-ytilde[i],2);
  double var=sum/(numInds-nEnv);
  
  double Vaa[nEnv*nEnv];
  for(int x=0;x<nEnv*nEnv;x++)
    Vaa[x]=0;
  double Vab[nEnv];
  for(int x=0;x<nEnv;x++)
    Vab[x]=0;
  
  
  for(int i=0 ; i<numInds ;i++) {
    Ex[i]=post[i*3+1]+2*post[i*3+2];
    Ex2[i]=post[i*3+1]+4*post[i*3+2];
    U+=Ex[i]*(y[i]-ytilde[i])/var;
    
    //Vaa<-Vaa+1/var*Xe[tal,]%*%t(Xe[tal,])
    for(int Nx=0;Nx<nEnv;Nx++)
      for(int Ny=0;Ny<nEnv;Ny++){
	Vaa[Nx*nEnv+Ny]+= (1/var)*cov[Nx*numInds+i]*cov[Ny*numInds+i];
      }
    
    for(int x=0;x<nEnv;x++)
      Vab[x]+= (1/var)*Ex[i]*cov[x*numInds+i];
    
    //Vab<-Vab+1/var*Ex[tal]*cbind(Xe[tal,])
    
    if(post[i*3+1]>0.90)
      highHE++;
    if(post[i*3+0]>0.90||post[i*3+2]>0.90)
      highHO++;
  }//recursion done
  
  pars->results->asso->highHe[s] = highHE;
  pars->results->asso->highHo[s] = highHO;
  
  for(int i =0; i<numInds;i++)
    sumEx+=Ex[i];
  
  
  
  
  //  double Vab=sumEx/var;
  double Vbb=0;
  for(int i =0; i<numInds;i++)
    Vbb+=(1/var-pow(y[i]-ytilde[i],2)/pow(var,2))*Ex2[i]+pow(y[i]-ytilde[i],2)/pow(var,2)*pow(Ex[i],2);
  
  //I<-Vbb-t(Vab)%*%MASS::ginv(Vaa)%*%Vab
  double workspace[2*nEnv];
  rankProb=angsd::matinv(Vaa, nEnv, nEnv, workspace);
  
  double I =0;
  
  //inv(Vaa)%*%Vab
  double invVaa_Vab[nEnv];
  for(int x=0;x<nEnv;x++)
    invVaa_Vab[x]=0;
  
  //NB! Vaa is now the inverse 
  for(int Nx=0;Nx<nEnv;Nx++)
    for(int Ny=0;Ny<nEnv;Ny++)
      invVaa_Vab[Nx]+=Vaa[Nx*nEnv+Ny]*Vab[Ny];
  
  //I<-t(Vab)%*%MASS::ginv(Vaa)%*%Vab
  for(int x=0;x<nEnv;x++)
    I+=Vab[x]*invVaa_Vab[x];
  //I<-Vbb-t(Vab)%*%MASS::ginv(Vaa)%*%Vab
  //s  fprintf(stderr,"tVab_invVaa_Vab: %f\n",I);  
  I=Vbb-I;


  //the observed varians of the dispersion 
  double Vbs=0;
  for(int i =0; i<numInds;i++)
    Vbs+=Ex[i]*(y[i]-ytilde[i])/pow(var,2);
  
  double Vss=0;
  for(int i =0; i<numInds;i++)
    Vss+=pow(y[i],2)+2*y[i]*ytilde[i];
  Vss=Vss*pow(var,-3)-numInds/(4*M_PI*pow(var,2));

  //fprintf(stderr,"Vbs %f Vss %f\n",Vbs,Vss);
  I=I-pow(Vbs,2)/Vss;

  
  double lrt =pow(U,2)/I;

  if((highHE < minHigh) || (highHO<minHigh) )
    lrt=-999;//set_snan(lrt);
  if(freq*numInds*2 < minCount || (1-freq)*numInds*2 < minCount)
    lrt=-999;//set_snan(lrt);      set_snan(lrt);
  if(rankProb!=0)
    lrt=-99;
  //      fprintf(lrtfile,"\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",numInds,ytilde,var,U,Vaa,Vab,Vbb,I,lrt);
  //fprintf(lrtfile,"\t%d\t%f",numInds,lrt);
  /*
  if(verbose>0)
    fprintf(fp,"\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",numInds,ytilde[0],var,U,Vaa[0],Vab[0],Vbb,I,lrt);
    else
    fprintf(fp,"\t%d\t%f",numInds,lrt);
  */

  if((0||lrt>1000||I<-0.01||lrt<0)&&lrt!=-999&&lrt!=-99){//!std::isnan(lrt)){
    for(int i=0 ; i<numInds ;i++) {
      fprintf(stderr,"y: %f\t  yfit: %f \t post %f %f %f\tEx %f %f\tU %f cov: ",y[i],ytilde[i],post[i*3+0],post[i*3+1],post[i*3+2],Ex[i],Ex2[i],U);
      for(int j=0;j<nEnv;j++)
	fprintf(stderr,"%f\t",cov[j*numInds+i]);
      fprintf(stderr,"\n");
 
    }
    for(int j=0;j<pow(nEnv,2);j++)
      fprintf(stderr,"Vaa: %f\t",Vaa[j]); 
    fprintf(stderr,"rank %d\tlrt: %f\t",rankProb,lrt); 
    fprintf(stderr,"\n");                                                                                                            
 
  
    for(int j=0;j<nEnv;j++)
      fprintf(stderr,"Vab: %f\t",Vab[j]);
    fprintf(stderr,"\n");

    fflush(stderr);
    fflush(stdout);
    exit(0);
  }
  
  delete []  Ex;
  delete []  Ex2;
  return lrt;

}
  


double asso::binomScoreEnv(double *post,int numInds, double *y, double *ytilde,double *cov,int nEnv,double freq,funkyPars *pars,int s){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  double *Ex = angsd::allocArray<double>(numInds,0);
  double *Ex2 = angsd::allocArray<double>(numInds,0);
  int rankProb=0;

  double U=0;
  int highHE=0;
  int highHO=0;
  double sumEx=0;
  double Vaa[nEnv*nEnv];
  for(int x=0;x<nEnv*nEnv;x++)
    Vaa[x]=0;
  double Vab[nEnv];
  for(int x=0;x<nEnv;x++)
    Vab[x]=0;


  for(int i=0 ; i<numInds ;i++) {
    Ex[i]=post[i*3+1]+2*post[i*3+2];
    Ex2[i]=post[i*3+1]+4*post[i*3+2];
    U+=Ex[i]*(y[i]-ytilde[i]);   

  // Vaa<-Vaa+yTilde[i]*(1-yTilde[i])*A[i,]%*%t(A[i,])
    for(int Nx=0;Nx<nEnv;Nx++)
      for(int Ny=0;Ny<nEnv;Ny++){
	Vaa[Nx*nEnv+Ny]+= ytilde[i]*(1-ytilde[i])*cov[Nx*numInds+i]*cov[Ny*numInds+i];
      }

    for(int x=0;x<nEnv;x++)
      Vab[x]+= ytilde[i]*(1-ytilde[i])*Ex[i]*cov[x*numInds+i];

    //Vba<-Vba+yTilde[i]*(1-yTilde[i])*A[i,]*Ex[i]

    if(post[i*3+1]>0.9)
      highHE++;
    if(post[i*3+0]>0.9||post[i*3+2]>0.9)
      highHO++;
  }//recursion done
  pars->results->asso->highHe[s] = highHE;
  pars->results->asso->highHo[s] = highHO;
 

    for(int i =0; i<numInds;i++)
      sumEx+=Ex[i];
    // double Vaa=ytilde[0]*(1-ytilde[0])*numInds;
    //    double Vab=ytilde[0]*(1-ytilde[0])*sumEx;
    double Vbb=0;
    for(int i =0; i<numInds;i++)
      Vbb+=(ytilde[i]*(1-ytilde[i])-pow(y[i]-ytilde[i],2))*Ex2[i]+pow(y[i]-ytilde[i],2)*pow(Ex[i],2);

    double workspace[2*nEnv];
    rankProb=angsd::matinv(Vaa, nEnv, nEnv, workspace);
 

    double I =0;

    double invVaa_Vab[nEnv];
    for(int x=0;x<nEnv;x++)
      invVaa_Vab[x]=0;

    for(int Nx=0;Nx<nEnv;Nx++)
      for(int Ny=0;Ny<nEnv;Ny++)
	invVaa_Vab[Nx]+=Vaa[Nx*nEnv+Ny]*Vab[Ny];

    for(int x=0;x<nEnv;x++)
      I+=Vab[x]*invVaa_Vab[x];
    I=Vbb-I;

    double lrt =pow(U,2)/I;

    if(highHE <minHigh||highHO <minHigh)
      lrt=-999;
    //freq*numInds*2 is the expected number of minor alleles
    if(freq*numInds*2 < minCount || (1-freq)*numInds*2 < minCount)
      lrt=-999;
    if(rankProb!=0)
      lrt=-99;
    //dispersion matrix has zero corners

    if((lrt>1000||I<-0.01||lrt<0)&&lrt!=-999&&lrt!=-99){//!std::isnan(lrt)){
      for(int i=0 ; i<numInds ;i++) {
	fprintf(stderr,"y: %f\t  post %f %f %f\tEx %f %f\tU %f\n",y[i],post[i*3+0],post[i*3+1],post[i*3+2],Ex[i],Ex2[i],U);
      }
      exit(0);
    }

    delete []  Ex;
    delete []  Ex2;
   
    return lrt;
}





void asso::printDoAsso(funkyPars *pars){
  if(doPrint)
    fprintf(stderr,"staring [%s]\t[%s]\n",__FILE__,__FUNCTION__);

  for(int yi=0;yi<ymat.y;yi++){
    bufstr.l=0;
    for(int s=0;s<pars->numSites;s++){
      if(pars->keepSites[s]==0){//will skip sites that have been removed      
	continue;
     } 
      if(doAsso==2){
	ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%d\t%f\t%d\t%d\n",header->name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],pars->results->asso->freq[s],pars->results->asso->keepInd[yi][s],pars->results->asso->stat[yi][s],pars->results->asso->highHe[s],pars->results->asso->highHo[s]);
	//	fprintf(MultiOutfile[yi],"%s\t%d\t%c\t%c\t%f\t%d\t%f\t%d\t%d\n",header->name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],pars->results->asso->freq[s],pars->results->asso->keepInd[yi][s],pars->results->asso->stat[yi][s],pars->results->asso->highHe[s],pars->results->asso->highHo[s]);
      }else{
	ksprintf(&bufstr,"%s\t%d\t%c\t%c\t%f\t%f\n",header->name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],pars->results->asso->freq[s],pars->results->asso->stat[yi][s]);
	//	fprintf(MultiOutfile[yi],"%s\t%d\t%c\t%c\t%f\t%f\n",header->name[pars->refId],pars->posi[s]+1,intToRef[pars->major[s]],intToRef[pars->minor[s]],pars->results->asso->freq[s],pars->results->asso->stat[yi][s]);
      }
    }
    gzwrite(MultiOutfile[yi],bufstr.s,bufstr.l);
  }
}
