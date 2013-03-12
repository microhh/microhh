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
    int readinifile (std::string);
    int readproffile(std::string);
    int getItem(int *        , std::string, std::string, std::string el ="default");
    int getItem(double *     , std::string, std::string, std::string el ="default");
    int getItem(bool *       , std::string, std::string, std::string el ="default");
    int getItem(std::string *, std::string, std::string, std::string el ="default");
    int getItem(int *        , std::string, std::string, int        , std::string el ="default");
    int getItem(double *     , std::string, std::string, double     , std::string el ="default");
    int getItem(bool *       , std::string, std::string, bool       , std::string el ="default");
    int getItem(std::string *, std::string, std::string, std::string, std::string el ="default");
    int getProf(double *     , std::string, int size);
    int clear();

  private:
    cmpi *mpi;
    int checkItemExists(std::string, std::string, std::string el="default");
    int checkItem(int *   , std::string, std::string, std::string el="default");
    int checkItem(double *, std::string, std::string, std::string el="default");
    int checkItem(bool *  , std::string, std::string, std::string el="default");
    int checkItem(std::string *  , std::string, std::string, std::string el="default");
    typedef std::map<std::string, std::string> inputmap1D;
    typedef std::map<std::string, inputmap1D> inputmap2D;
    typedef std::map<std::string, inputmap2D> inputmap;
    inputmap inputlist;
    typedef std::map<std::string, std::vector<double> > profmap;
    profmap proflist;
    std::vector<std::string> varnames;
    std::vector<double> varvalues;
};
#endif
