/*
  small class to filter positions and assign major minor.
  This class can be improved very much.

  Thorfinn 7march 2013

  On the most basic level then this class will remove those sites with an effective sample size<minInd

  It can also be used for filtering away those sites not included in the -filter file.keep

  Or you can supply a bimfile, and then the major/minor from the bim file will be used as major/minor throughout the program


 */
#include <cassert>
#include "pthread.h"
#include "shared.h"
#include "general.h"
#include "analysisFunction.h"
#include "analysisKeepList.h"

#define BIN ".bin"
#define IDX ".idx"


//return one+two
char *append(const char *one,const char*two){
  char *ret = new char[strlen(one)+strlen(two)+1];
  strcpy(ret,one);
  strcpy(ret+strlen(ret),two);
  return ret;
}


filt *filt_read(const char *fname){
  fprintf(stderr,"\t->[%s] Reading binary representation:%s\n",__FILE__,fname);
  filt *ret = new filt;
  ret->bg =NULL;
  ret->fp =NULL;
  ret->keeps = ret->major = ret->minor=NULL;
  ret->nCols =0;
  ret->curLen =0;
  char *bin_name=append(fname,BIN);
  char *idx_name=append(fname,IDX);

  ret->fp= fopen(idx_name,"r");
  ret->bg=bgzf_open(bin_name,"r");


  while(1){
    int chrId;
    if(0==fread(&chrId,sizeof(int),1,ret->fp))
      break;
    asdf_dats tmp;
    if(1!=fread(&tmp.offs,sizeof(int64_t),1,ret->fp)){
      fprintf(stderr,"Problem reading chunk from binary file:\n");
      exit(0);
    }
    if(1!=fread(&tmp.len,sizeof(int),1,ret->fp)){
      fprintf(stderr,"Problem reading chunk from binary file:\n");
      exit(0);
    }
    if(1!=fread(&ret->nCols,sizeof(int),1,ret->fp)){
      fprintf(stderr,"Problem reading chunk from binary file:\n");
      exit(0);
    }
    ret->offs[chrId]=tmp;
  }
  
  std::map<int,asdf_dats>::const_iterator it;
  //  fprintf(stderr,"nCols: %d\n",ret->nCols);
  for(it=ret->offs.begin();0&&it!=ret->offs.end();++it)
    fprintf(stderr,"id:%d offs:%zu len:%d\n",it->first,it->second.offs,it->second.len);
  fprintf(stderr,"\t-> [%s] nChr: %zu loaded from binary filter file\n",__FILE__,ret->offs.size());
  if(ret->nCols==4)
    fprintf(stderr,"Filterfile contains major/minor information\n");

  return ret;
}





void filt_gen(const char *fname,std::map<char*,int,ltstr>* revMap,aHead *hd){
  fprintf(stderr,"Filterfile: %s supplied will generate binary representations...\n",fname);
  std::map<char*,int,ltstr>::const_iterator it;
  for(it= revMap->begin();0&&it!=revMap->end();++it)
    fprintf(stderr,"%s->%d->%d\n",it->first,it->second,hd->l_ref[it->second]);
  
  gzFile gz = Z_NULL;
  gz = gzopen(fname,"r");
  if(gz==Z_NULL){
    fprintf(stderr,"Problem opening file:%s\n",fname);
    exit(0);
  }

  char* outnames_bin = append(fname,BIN);
  char* outnames_idx = append(fname,IDX);
    
  const char *delims = "\t \n";
  BGZF *cfpD = bgzf_open(outnames_bin,"w9");
  FILE *fp =fopen(outnames_idx,"w");
  
  





  std::map <int,char> mm;//simple structure to check that input has been sorted by chr/contig
  char *ary = NULL;
  char *major = NULL;
  char *minor = NULL;
  int last=-1;
  int nCols = -1;
  char buf[LENS];
  //  gzread(gz,buf,LENS);

  ////  char *rd=gzgets(gz,buf,LENS);
  //fprintf(stderr,"rd=%s LENS=%d\n",buf,LENS);
  while(gzgets(gz,buf,LENS)){
    char chr[LENS] ;int id=-1;
    int posi=-1;
    char maj='N';
    char min='N';
    
    int nRead=sscanf(buf,"%s\t%d\t%c\t%c\n",chr,&posi,&maj,&min);
    
    posi--;
    assert(nRead!=0);
    if(nCols==-1)
      nCols=nRead;
#if 0
    fprintf(stderr,"nRead=%d: %s %d %c %c\n",nRead,chr,posi,maj,min);
    exit(0);
#endif    
    assert(nRead==nRead);
    it = revMap->find(chr);
    if(it==revMap->end()){
      fprintf(stderr,"chr: %s from filterfile: %s doesn't exist in index, will exit()\n",chr,fname);
      exit(0);
    }else
      id=it->second;
    
    /*
      1. if we have observed a chromo change then dump data
      2. realloc
      3. plug in values
      (need to check that we have alloced a data <=> last!=-1)
     */
    if(last!=-1 &&last!=id){
      //      fprintf(stderr,"writing index: last=%d id=%d\n",last,id);
      assert(ary!=NULL);
      //write data and index stuff
      int64_t retVal =bgzf_tell(cfpD);
      fwrite(&last,1,sizeof(int),fp);
      fwrite(&retVal,1,sizeof(int64_t),fp);
      fwrite(&hd->l_ref[last],1,sizeof(int),fp);
      fwrite(&nCols,1,sizeof(int),fp);
      bgzf_write(cfpD,ary,hd->l_ref[last]);//write len of chr
      if(nCols==4){
	bgzf_write(cfpD,major,hd->l_ref[last]);//write len of chr
	bgzf_write(cfpD,minor,hd->l_ref[last]);//write len of chr
      }
    }
    if(last!=id){
      fprintf(stderr,"Allocing chr:%s\n",chr);
      std::map<int,char>::iterator it=mm.find(id);
      if(it!=mm.end()){
	fprintf(stderr,"filter file, doesn't look sorted by chr, will exit()");
	exit(0);
      }else
	mm[id]=1;
      last=id;
      ary=new char[hd->l_ref[last]];
      memset(ary,0,hd->l_ref[last]);
      //      fprintf(stderr,"chr: %s id=%d len=%d\n",chr,last,hd->l_ref[last]);
      if(nCols==4){
	major=new char[hd->l_ref[last]];
	memset(major,0,hd->l_ref[last]);
	minor=new char[hd->l_ref[last]];
	memset(minor,0,hd->l_ref[last]);
      }
      
    }else
      if(posi > hd->l_ref[id]){
	fprintf(stderr,"Position in filter file:%s is after end of chromosome? Will exit\n",fname);
	exit(0);
      }
      ary[posi] = 1;
      if(nCols==4){
	major[posi] = refToInt[maj];
	minor[posi] = refToInt[min];
      }
  }
  if(last!=-1){
    //fprintf(stderr,"writing index\n");
    assert(ary!=NULL);
    //write data and index stuff
    int64_t retVal =bgzf_tell(cfpD);
    fwrite(&last,1,sizeof(int),fp);
    fwrite(&retVal,1,sizeof(int64_t),fp);
    fwrite(&hd->l_ref[last],1,sizeof(int),fp);
    fwrite(&nCols,1,sizeof(int),fp);
    bgzf_write(cfpD,ary,hd->l_ref[last]);//write len of chr
    if(nCols==4){
      bgzf_write(cfpD,major,hd->l_ref[last]);//write len of chr
      bgzf_write(cfpD,minor,hd->l_ref[last]);//write len of chr
    }
  }

  fprintf(stderr,"Filtering complete: Observed: %zu different chromosomes from file:%s\n",mm.size(),fname);

  mm.clear();
  gzclose(gz);fclose(fp);bgzf_close(cfpD);

}




filt *filt_init(const char *fname,std::map<char*,int,ltstr>* revMap,aHead *hd){
  char *bin_name=append(fname,BIN);
  char *idx_name=append(fname,IDX);
  if(!fexists(bin_name)||!fexists(idx_name))
    filt_gen(fname,revMap,hd);
  return  filt_read(fname);
}



filter::~filter(){
  if(fl!=NULL){
    if(fl->keeps)
      free(fl->keeps);
    if(fl->major)
      free(fl->major);
    if(fl->minor)
      free(fl->minor);
    fclose(fl->fp);
    bgzf_close(fl->bg);
    delete fl;
  }
}


filter::filter(argStruct *arguments){
  //below if shared for all analysis classes
  header = arguments->hd;
  revMap = arguments->revMap;
  fl = NULL;
  //his is used by this class
  keepsChr = NULL;
  curChr = -1;
  fp = NULL;
  minInd = 0;
  fname = NULL;
  doMajorMinor =0;
  
  if(arguments->argc==2){
    if(!strcmp(arguments->argv[1],"-filter")){
      printArg(stdout);
      exit(0);
    }else
      return;
  }

  //get options and print them
  getOptions(arguments);
  printArg(arguments->argumentFile);

}


void filter::printArg(FILE *argFile){
  fprintf(argFile,"--------------\n%s:\n",__FILE__);
  fprintf(argFile,"\t-filter\t\t%s \n",fname);
  fprintf(argFile,"\t-doMajorMinor\t%d\t\n",doMajorMinor);
  fprintf(argFile,"\t1: Infer major and minor from GL\n");
  fprintf(argFile,"\t2: Infer major and minor from allele counts\n");
  fprintf(argFile,"\t3: use major and minor from bim file (requires -filter afile.bim)\n");
  fprintf(argFile,"\t4: Use reference allele as major (requires -ref)\n");
  fprintf(argFile,"\t5: Use ancestral allele as major (requires -anc)\n");
  fprintf(argFile,"\t-minInd\t\t%d\tOnly use site if atleast minInd of samples has data\n",minInd);  fprintf(argFile,"\n");
}

void filter::getOptions(argStruct *arguments){
  fname=angsd::getArg("-filter",fname,arguments);

  if(fname!=NULL)  
    fl = filt_init(fname,revMap,header);
  if(fl!=NULL)
    fprintf(stderr,"\t-> [%s] -filter is still beta, use at own risk...\n",__FILE__);


  //1=bim 2=keep
  doMajorMinor = angsd::getArg("-doMajorMinor",doMajorMinor,arguments);
  if(doMajorMinor==3 && fl!=NULL&& fl->nCols!=4){
    fprintf(stderr,"Must supply -filter with a file containing major and minor if -doMajorMinor 3\n");
  }
  if(doMajorMinor!=3 && fl!=NULL&& fl->nCols==4){
    fprintf(stderr,"Filter file contains major/minor information to use these in analysis supper \'-doMajorMinor 3\'\n");
  }
  
  minInd = angsd::getArg("-minInd",minInd,arguments);
}


void filter::run(funkyPars *p){
  //  fprintf(stderr,"nsites=%d\n",p->numSites);
  p->keepSites=new int[p->numSites];
  
  for(int s=0;s<p->numSites;s++){
    p->keepSites[s]=p->nInd;
    //    p->results->freq->keepInd[s]=nInd;  
  }


  if(fl!=NULL && fl->nCols==4){
    //    fprintf(stderr,"aloocating for major and minor\n");
    p->major = new char [p->numSites];
    p->minor = new char [p->numSites];
    for(int i=0;i<p->numSites;i++){
      p->major[i] = 4;
      p->minor[i] = 4;
    }
  }

  if(fl!=NULL) {
    for(int s=0;s<p->numSites;s++){

      if(fl->keeps[p->posi[s]]==0){
	//	fprintf(stderr,"Plugging inf vals std\n");
	p->keepSites[s] =0;
      }
      if(p->keepSites[s] && fl->nCols==4){
	//fprintf(stderr,"Plugging inf vals std majorminor\n");
	p->major[s] = fl->major[p->posi[s]];
	p->minor[s] = fl->minor[p->posi[s]];
      }
    }
  }
  //how set the keepsites according the effective sample size persite
  if(0!=minInd){
    if(p->chk!=NULL){
      //loop over sites;
      for(int s=0;s<p->numSites;s++){
	if(p->keepSites[s]==0)
	  continue;
	int nInfo =0;
	tNode *tn = p->chk->nd[s];
	//loop over samples;
	for(int i=0;i<p->nInd;i++)
	  if(tn[i].l!=0)
	    nInfo++;
	//fprintf(stdout,"%d %d\n",nInfo,minInd);
	if(minInd<=nInfo)
	  p->keepSites[s] =nInfo;
	else
	  p->keepSites[s] =0;
      }
    
    }
  }
}
void filter::print(funkyPars *p){
}

void filter::clean(funkyPars *p){
  
}


void filter::readSites(int refId) {
  //  fprintf(stderr,"[%s].%s():%d -> refId:%d\n",__FILE__,__FUNCTION__,__LINE__,refId);
  if(fl==NULL)
    return;
  std::map<int,asdf_dats> ::iterator it = fl->offs.find(refId);
  if(it==fl->offs.end()){
    fprintf(stderr,"[%s.%s():%d] Problem finding chrId: %d in indexed filtering file, consider running with '-r chr:' \n",__FILE__,__FUNCTION__,__LINE__,refId);
    exit(0);
  }
  bgzf_seek(fl->bg,it->second.offs,SEEK_SET);
  if(it->second.len>fl->curLen) 
    fl->keeps=(char*) realloc(fl->keeps,it->second.len);
  bgzf_read(fl->bg,fl->keeps,it->second.len);

  if(fl->nCols==4){
    if(it->second.len>fl->curLen) {
      fl->major = (char*) realloc(fl->major,it->second.len);
      fl->minor = (char*) realloc(fl->minor,it->second.len);
    }
    bgzf_read(fl->bg,fl->major,it->second.len);
    bgzf_read(fl->bg,fl->minor,it->second.len);
  }
 
}
