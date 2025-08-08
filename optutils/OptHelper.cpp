#include "OptHelper.h"
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <string>
using namespace std;

OptEntry::OptEntry(char short_opt, const char * long_opt, int need_arg, const char * arg_name, const char * help_description, char data_type, void * data, bool multi, void*(*funAfter)(void *), void * funArgs)//data should be allocated pointer with correct data type
    :ShortOpt(short_opt),LongOpt(long_opt),NeedArg(need_arg),ArgName(arg_name),HelpDescription(help_description),DataType(data_type),Data(data),DefaultData(Data), Multi(multi), funAfter(funAfter), funArgs(funArgs){
        DefaultValue="";
        if (data==NULL) return;
        if (data_type=='s')
        {
            if ((*(const char **)data)==NULL) return;
            DefaultValue+="\"";
            DefaultValue+=(*(const char **)data);
            DefaultValue+="\"";
        }
        else if (data_type=='S')
        {
            DefaultValue+="\"";
            DefaultValue+=(*(string*)data);
            DefaultValue+="\"";
        }
        else if (data_type=='i')
        {
            DefaultValue=to_string(*(int *)data);
        }
        else if (data_type=='F')
        {
            DefaultValue=to_string(*(double *)data);
        }
        else if (data_type=='b')
        {
            if (*((bool *)data)==true) DefaultValue="true";
            else DefaultValue="false";
        }
    }

string OptEntry::getOptName() const
{
    if (ShortOpt!=0) return string("-")+ShortOpt;
    return string("--")+LongOpt;
}

void OptHelper::addOpt(char short_opt, const char * long_opt, int need_arg, const char * arg_name, const char * help_description, char data_type, void * data, bool multi, void*(*funAfter)(void *), void * funArgs)//data should be allocated pointer with correct data type
{
    Opts.push_back(OptEntry(short_opt, long_opt, need_arg, arg_name, help_description, data_type, data, multi, funAfter, funArgs));
}

char * get_short_opts(const vector<OptEntry> &options)
{
    char* result=(char*) malloc(1000);//>256*3

    unsigned j=0;
    for (unsigned i = 0; i<options.size();++i)
    {
        if (options[i].ShortOpt != 0)
        {
            result[j++]=options[i].ShortOpt;
            switch (options[i].NeedArg)
            {
            case no_argument:
                break;
            case optional_argument:
                result[j++] = ':';
            case required_argument:
                result[j++] = ':';
                break;
            }
        }
    }
    result[j]='\0';
    return result;
}

struct option *get_long_opts(const vector<OptEntry> &options)
{
    unsigned count = 0;
    unsigned i = 0;
    struct option *result = (struct option *)calloc((options.size() + 1) , sizeof(struct option));
    unsigned j = 0;
    for (i = 0; i<options.size(); ++i)
    {
        if (options[i].LongOpt != NULL)
        {
            result[j].name = options[i].LongOpt;
            result[j].has_arg = options[i].NeedArg;
            result[j].flag = NULL;
            if (options[i].ShortOpt!=0)
                result[j++].val = options[i].ShortOpt;
            else
                result[j++].val = i+1000;
        }
    }
    result[j].name = NULL;
    return result;
}

int OptHelper::getOpts(int argc, const char ** argv)
{
    char * Shorts=get_short_opts(Opts);
    struct option * Longs=get_long_opts(Opts);

    int longind;
    int shortopt;
    int ShortIndex[256];
    for (int i=0;i<256;++i)
    {
        ShortIndex[i]=-1;
    }
    for (int i=0;i<Opts.size();++i)
    {
        if (Opts[i].ShortOpt!=0)
            ShortIndex[Opts[i].ShortOpt]=i;
    }
    int Index;
    while ((shortopt=getopt_long(argc, (char * const *)argv, Shorts, Longs, &longind))>=0)
    {
        if (shortopt==':' || shortopt=='?')
        {
            free(Shorts);
            free(Longs);
            return 0;
        }
        if (shortopt>=1000) Index=shortopt-1000;
        else Index=ShortIndex[shortopt];
        assert(Index>=0);
        if (optarg==0)//no arg
        {
            if (Opts[Index].DataType=='b' && Opts[Index].Multi==false)
            {
                *(bool*)Opts[Index].Data=!(*(bool*)Opts[Index].Data);
            }
            else
            {
                char * OptName=(char *)malloc(Opts[Index].LongOpt!=NULL?strlen(Opts[Index].LongOpt)+10:10);
                if (Opts[Index].ShortOpt!=0)
                {
                    OptName[0]='-';
                    OptName[1]=Opts[Index].ShortOpt;
                    OptName[2]='\0';
                }
                if (Opts[Index].LongOpt!=NULL)
                {
                    strcat(OptName," --");
                    strcat(OptName,Opts[Index].LongOpt);
                }
                fprintf(stderr, "[WARN] Parsing %s with no arg!\n",OptName);
                free(OptName);
            }
            OptPairs.push_back(Opts[Index].getOptName());
        }
        else
        {
            if (Opts[Index].Multi)
            {
                switch (Opts[Index].DataType)
                {
                    case 'i':
                    ((vector<int> *)(Opts[Index].Data))->push_back(atoi(optarg));
                    break;
                    case 'F':
                    ((vector<double> *)(Opts[Index].Data))->push_back(atof(optarg));
                    break;
                    case 'b':
                    ((vector<bool>*)(Opts[Index].Data))->push_back(atoi(optarg));
                    break;
                    case 's':
                    ((vector<const char *>*)(Opts[Index].Data))->push_back(optarg);
                    //*(const char**)Opts[Index].Data=optarg;
                    break;
                    case 'S':
                    ((vector<string>*)(Opts[Index].Data))->push_back(optarg);
                    break;
                }
            }
            else
            {
                switch (Opts[Index].DataType)
                {
                    case 'i':
                    *(int*)(Opts[Index].Data)=atoi(optarg);
                    break;
                    case 'F':
                    *(double*)(Opts[Index].Data)=atof(optarg);
                    break;
                    case 'b':
                    *(bool*)(Opts[Index].Data)=atoi(optarg);
                    break;
                    case 's':
                    *(const char**)(Opts[Index].Data)=optarg;
                    break;
                    case 'S':
                    *(string*)(Opts[Index].Data)=optarg;
                    break;
                }
            }
            OptPairs.push_back(Opts[Index].getOptName()+" "+string(optarg));
        }
        if (Opts[Index].funAfter!=NULL) Opts[Index].funAfter(Opts[Index].funArgs);
    }
    if (shortopt<0)
    {
        for (int i=optind;i<argc;++i) Args.push_back(argv[i]);
    }

    free(Shorts);
    free(Longs);
    return 1;
}

OptHelper::OptHelper(const char * usage)
:Usage(usage){}

void OptHelper::setUsage(const char * usage)
{
    Usage=usage;
}

void OptHelper::showhelp(FILE* file)
{
    if (Usage!=NULL) fprintf(file,"%s\n",Usage);
    fprintf(file,"Options:\n");
    for (unsigned i=0;i<Opts.size(); ++i)
    {
        fprintf(file,"  ");
        unsigned j = 0;
        if (Opts[i].ShortOpt != 0)
        {
            fprintf(file, "-%c", Opts[i].ShortOpt);
            j = 2;
            if (Opts[i].NeedArg != no_argument)
            {
                if (Opts[i].NeedArg == optional_argument)
                {
                    fprintf(file, "[");
                    ++j;
                }
                else
                {
                    fprintf(file, " ");
                    ++j;
                }
                if (Opts[i].ArgName == NULL)
                {
                    fprintf(file, "ARG");
                    j += 3;
                }
                else
                {
                    unsigned k = 0;
                    while (Opts[i].ArgName[k] != '\0')
                    {
                        fprintf(file, "%c", Opts[i].ArgName[k++]);
                        if ((j + k) % 28 == 0)
                            fprintf(file, "\n  ");
                    }
                    j += k;
                }
                if (Opts[i].NeedArg == optional_argument)
                {
                    fprintf(file, "]");
                    ++j;
                    if (j % 28 == 0)
                        fprintf(file, "\n  ");
                }
            }
        }

        if (Opts[i].LongOpt != NULL)
        {
            if (Opts[i].ShortOpt != 0)
            {
                fprintf(file, ",");
                ++j;
                if (j % 28 == 0)
                    fprintf(file, "\n  ");
                fprintf(file, " ");
                ++j;
                if (j % 28 == 0)
                    fprintf(file, "\n  ");
            }
            fprintf(file, "-");
            ++j;
            if (j % 28 == 0)
                fprintf(file, "\n  ");
            fprintf(file, "-");
            ++j;
            if (j % 28 == 0)
                fprintf(file, "\n  ");
            unsigned k = 0;
            while (Opts[i].LongOpt[k] != '\0')
            {
                fprintf(file, "%c", Opts[i].LongOpt[k++]);
                if ((j + k) % 28 == 0)
                    fprintf(file, "\n  ");
            }
            j += k;
            if (Opts[i].NeedArg != no_argument)
            {
                if (Opts[i].NeedArg == optional_argument)
                {
                    fprintf(file, "[");
                    ++j;
                    if (j % 28 == 0)
                        fprintf(file, "\n  ");
                }

                fprintf(file, "=");
                ++j;
                if (j % 28 == 0)
                    fprintf(file, "\n  ");
                char * arg_name=(char*) malloc(1024);/*
                char * argdefault=(char*) malloc(64);
                if (Opts[i].DataType=='i'|| Opts[i].DataType=='b')
                {
                    snprintf(argdefault,64,"%d",*(int*)Opts[i].Data);
                }
                else if(Opts[i].DataType=='F')
                {
                    snprintf(argdefault,64,"%lf",*(double*)Opts[i].Data);
                }
                else
                {
                    snprintf(argdefault,64,"%s",*(const char **)Opts[i].Data);
                }
                argdefault[63]=0;*/
                strcpy(arg_name,Opts[i].ArgName == NULL ? "ARG" : Opts[i].ArgName);
                //strcat(arg_name,"(");
                //strcat(arg_name,argdefault);
                //strcat(arg_name,")");
                unsigned k = 0;
                while (arg_name[k] != '\0')
                {
                    fprintf(file, "%c", arg_name[k++]);
                    if ((j + k) % 28 == 0)
                        fprintf(file, "\n  ");
                }
                j += k;
                free(arg_name);
                //free(argdefault);
                if (Opts[i].NeedArg == optional_argument)
                {
                    fprintf(file, "]");
                    ++j;
                }
            }
        }
        for (unsigned k = 0; k < 32 - (j % 28); ++k)
            fprintf(file, " ");
        j = 32;
        if (Opts[i].HelpDescription != NULL)
            fprintf(file, Opts[i].HelpDescription);
            if (Opts[i].DefaultValue!="") fprintf(file, " (%s).",Opts[i].DefaultValue.c_str());
        fprintf(file, "\n");
    }
}

bool assignWithType(void * a, const void * b, char Type)
{
    if (Type=='i')
    {
        *((int*)a)=*((const int*)b);
        return true;
    }
    if(Type=='F')
    {
        *((double*)a)=*((const double*)b);
        return true;
    }
    if(Type=='s')
    {
        *((const char **)a)=*((const char **)b);
        return true;
    }
    if(Type=='b')
    {
        *((bool*)a)=*((const bool*)b);
        return true;
    }
    return false;
}

void OptHelper::initData()
{
    for (int i=0;i<Opts.size();++i)
    {
        if (Opts[i].Data!=NULL && Opts[i].DefaultData!=NULL) assignWithType(Opts[i].Data,Opts[i].DefaultData,Opts[i].DataType);
    }
}