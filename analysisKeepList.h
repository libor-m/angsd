#pragma once

typedef struct{
  int major;
  int minor;
}mm;


//this below is uberaddhoc. Under the assumption of the .bim file being ordered by chromosome. 
//We can avoid the first key. And simply compare by single ints eg the position, this should be done DRAGON

struct cmp_ints {
  bool operator()(const mm& en,const  mm& to) {
    if(en.major!=to.major)
      return en.major<to.major;
    else
      return en.minor<to.minor;
  }
};

typedef std::map<mm,mm,cmp_ints > fMap;





class filter : public general{
  int doMajorMinor;
  int doFilter;
  char *fname;
  fMap fm;
  FILE *fp;

  char *keepsChr;//<- a char array of length=reflength, used for .keeps
  int minInd;
  //  int nInd;
  int curChr;
public:
  void readSites();
  //none optional stuff
  FILE *outfile;
  filter(argStruct *arguments);
  void getOptions(argStruct *arguments);
  void run(funkyPars  *pars);
  void clean(funkyPars *pars);
  void print(funkyPars *pars);
  void printArg(FILE *argFile);
  ~filter();
};
