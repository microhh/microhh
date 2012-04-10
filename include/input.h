#ifndef INPUT
#define INPUT

#include <map>
#include <string>

class cinput
{
  public:
    cinput();
    ~cinput();
    int readinifile();
    int getItem(int *, std::string, std::string, bool);
    int getItem(double *, std::string, std::string, bool);
    int getItem(bool *, std::string, std::string, bool);

  private:
    int checkItemExists(std::string, std::string, bool);
    typedef std::map<std::string, std::string> inputmapinner;
    typedef std::map<std::string, inputmapinner> inputmap;
    inputmap inputlist;
};
#endif
