#ifndef OPT_HELPER_H
#define OPT_HELPER_H

#include <stdio.h>
#include <vector>

struct OptEntry
{
    char ShortOpt;
    const char * LongOpt;
    int NeedArg;
    const char * ArgName;
    const char * HelpDescription;
    char DataType;//'i' for int, 'F' for double, 's' for string, 'b' for bool;
    void * Data;
    void * DefaultData;//if NULL then don't init the Data;
    bool Multi;
    OptEntry(char short_opt, const char * long_opt, int need_arg, const char * arg_name, const char * help_description, char data_type, void * data, bool multi=false);//data should be allocated pointer with correct data type
};

class OptHelper
{
    std::vector<OptEntry> Opts;
    const char * Usage;
    void initData();
    public:
    std::vector<const char *> Args;
    OptHelper(const char * usage=NULL);
    /*
    data and defaultdata should be allocated pointer with correct data type, recommend to initialize data outside then pass it to both data and default data
    if multi is true, then data should be a pointer to std::vector of the right type, so OptHelper can store options occurred multiple times. if not set multi, will overide the option time by time
    */
    void addOpt(char short_opt, const char * long_opt, int need_arg, const char * arg_name, const char * help_description, char data_type, void * data, bool multi=false);
    int getOpts(int argc, const char ** argv);
    void setUsage(const char * Usage);
    void showhelp(FILE* File=stderr);
};

#endif