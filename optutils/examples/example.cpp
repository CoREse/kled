#include "../OptHelper.h"
#include <stdio.h>

int main(int argc, const char ** argv)
{
    OptHelper OH=OptHelper("example nothing");
    int arg=0;
    const char * args="";
    bool argb=true;
    double argd;
    OH.addOpt('c', "cc", 1, 0, "1234556",'i',&arg);
    OH.addOpt('d', "double", 1, 0, "double",'F',&argd);
    OH.addOpt('b', 0, 2, "bool", "bool arg",'b',&argb);
    OH.addOpt('s', "string", 2, 0, "string arg",'s',&args);
    for (int i=0;i<argc;++i) printf("%s ",argv[i]);
    printf("\n");
    OH.getOpts(argc,argv);
    for (int i=0;i<argc;++i) printf("%s ",argv[i]);
    printf("\n");
    printf("number: %d, %lf, %s, %d\n", arg,argd, args,argb);
    for (int i=0;i<OH.Args.size();++i) printf("%s\n",OH.Args[i]);
    OH.showhelp();
    return 0;
}