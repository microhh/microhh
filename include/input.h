#ifndef INPUT
#define INPUT

#include <map>
#include <string>

class cinput
{
  public:
    cinput();
    ~cinput();
    int getItem(int *, std::string, std::string);
    int getItem(double *, std::string, std::string);
    int getItem(bool *, std::string, std::string);

  private:
    std::map<std::string, std::map<std::string, std::string> > inputlist;
};
#endif
