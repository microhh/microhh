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
    int getItem(int *   , std::string, std::string);
    int getItem(double *, std::string, std::string);
    int getItem(bool *  , std::string, std::string);
    int getItem(int *   , std::string, std::string, int);
    int getItem(double *, std::string, std::string, double);
    int getItem(bool *  , std::string, std::string, bool);
    int getProf(double *, std::string, int size);
    int clear();

  private:
    cmpi *mpi;
    int checkItemExists(std::string, std::string);
    int checkItem(int *   , std::string, std::string);
    int checkItem(double *, std::string, std::string);
    int checkItem(bool *  , std::string, std::string);
    typedef std::map<std::string, std::string> inputmapinner;
    typedef std::map<std::string, inputmapinner> inputmap;
    inputmap inputlist;
    typedef std::map<std::string, std::vector<double> > profmap;
    profmap proflist;
};
#endif
