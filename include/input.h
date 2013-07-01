#ifndef INPUT
#define INPUT

#include <map>
#include <string>
#include <vector>

// forward declaration to avoid circular dependency
class cmpi;

class cinput
{
  public:
    cinput(cmpi *);
    ~cinput();
    int readinifile ();
    int readproffile();
    int getItem(int *        , std::string, std::string, std::string);
    int getItem(int *        , std::string, std::string, std::string, int);
    int getItem(double *     , std::string, std::string, std::string);
    int getItem(double *     , std::string, std::string, std::string, double);
    int getItem(bool *       , std::string, std::string, std::string);
    int getItem(bool *       , std::string, std::string, std::string, bool);
    int getItem(std::string *, std::string, std::string, std::string);
    int getItem(std::string *, std::string, std::string, std::string, std::string);

    // list retrieval
    int getList(std::vector<int> *        , std::string, std::string, std::string);
    int getList(std::vector<std::string> *, std::string, std::string, std::string);

    int getProf(double *     , std::string, int size);
    int clear();
    int printUnused();

  private:
    cmpi *mpi;

    template <class valuetype>
    int parseItem(valuetype *, std::string, std::string, std::string, bool, valuetype);

    template <class valuetype>
    int parseList(std::vector<valuetype> *, std::string, std::string, std::string);

    int checkItemExists(std::string, std::string, std::string el="default");
    int checkItem(int *        , std::string, std::string, std::string el="default");
    int checkItem(double *     , std::string, std::string, std::string el="default");
    int checkItem(bool *       , std::string, std::string, std::string el="default");
    int checkItem(std::string *, std::string, std::string, std::string el="default");

    // list retrieval
    int checkList(std::vector<int> *        , std::string, std::string, std::string el="default");
    int checkList(std::vector<std::string> *, std::string, std::string, std::string el="default");

    struct inputtype
    {
      std::string data;
      bool isused;
    };
    typedef std::map<std::string, inputtype> inputmap1d;
    typedef std::map<std::string, inputmap1d> inputmap2d;
    typedef std::map<std::string, inputmap2d> inputmap;
    inputmap inputlist;
    typedef std::map<std::string, std::vector<double> > profmap;
    profmap proflist;
    std::vector<std::string> varnames;
    std::vector<double> varvalues;
    std::string isused;
};
#endif
