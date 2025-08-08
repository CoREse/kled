#ifndef OPT_HELPER_H
#define OPT_HELPER_H

#include <stdio.h>
#include <vector>
#include <string>

struct OptEntry
{
    char ShortOpt;
    const char * LongOpt;
    int NeedArg;
    const char * ArgName;
    const char * HelpDescription;
    char DataType;//'i' for int, 'F' for double, 's' for c string, 'S' for string, 'b' for bool;
    void * Data;
    void * DefaultData;//if NULL then don't init the Data;
    bool Multi;
    std::string DefaultValue;
    void * (*funAfter) (void *);
    void * funArgs;
    OptEntry(char short_opt, const char * long_opt, int need_arg, const char * arg_name, const char * help_description, char data_type, void * data, bool multi=false, void *(*funAfter)(void *)=NULL, void* funArgs=NULL);//data should be allocated pointer with correct data type
    std::string getOptName() const;
};

class OptHelper
{
    std::vector<OptEntry> Opts;
    const char * Usage;
    void initData();
    public:
    std::vector<const char *> Args;
    std::vector<std::string> OptPairs;
    OptHelper(const char * usage=NULL);
    /*
    data and defaultdata should be allocated pointer with correct data type, recommend to initialize data outside then pass it to both data and default data
    if multi is true, then data should be a pointer to std::vector of the right type, so OptHelper can store options occurred multiple times. if not set multi, will overide the option time by time
    */
    void addOpt(char short_opt, const char * long_opt, int need_arg, const char * arg_name, const char * help_description, char data_type, void * data, bool multi=false, void *(*funAfter)(void *)=NULL, void* funArgs=NULL);
    int getOpts(int argc, const char ** argv);
    void setUsage(const char * Usage);
    void showhelp(FILE* File=stderr);
};

#endif