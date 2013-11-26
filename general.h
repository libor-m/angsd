#pragma once
#include "shared.h"
#include "bams.h"

class general{
public:
  static aHead *header;//contains the header of a single bam;
  static std::map<char *,int,ltstr> *revMap;
  int index;
  static int tot_index;
  //  virtuel general()
  virtual void run(funkyPars *f)=0;
  virtual void print( funkyPars *f)=0;
  //  virtual void printArg(const char *fname)=0; <-maybe include
  virtual void clean(funkyPars *f)=0;
  general(){index=tot_index++;};
  virtual ~general(){};
};


general **extra(int &nItem,const char *outfiles,int inputtype,argStruct *arguments);
