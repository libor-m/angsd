#include "generals.h"
#include "shared.h"

class fsb : generals{

  void run(funkyPars *p);
  void print(funkyPars *p);
  void clean();

};

void fsb::clean(){

}

void fsb::print(){
 
}

char *calcStuff(funkyPars *p){
  new ret = new char [p->numSites];
  for(int i=0;i<p->nInd;i++){
    for(int s=0;s<p->numSites;s++){
      int tmp[2]={0,0};
      nodeT tsk =chk->nd[i][s];      
      for (int c=0;c<tsk.l;c++){
	if(isUpper(tsk.seq[c]))
	  tmp[0]++;
	else
	  tmp[1]++;
	
      }
      if(tmp[0]!=0 && tmp[1]!=0)
	ret[s] = 0;

    }

      
}


void fsb::run(funkyPars *p){

  char *strIsBad = calcStuff(p);


  for(int i=0;i<p->numSites;p++){
    if(strIsBad[i])
      p->keepSites[i] =0;
    else
      p->keepSites[i] =1;
  }


}
