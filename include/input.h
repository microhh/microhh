#ifndef INPUT
#define INPUT

#include <map>
#include <vector>

class Input
{
    public:
        Input(const std::string&);
        template<typename T> T get_item(const std::string&, const std::string&, const std::string&);
        template<typename T> T get_item(const std::string&, const std::string&, const std::string&, const T);
        template<typename T> std::vector<T> get_list(const std::string&, const std::string&, const std::string&);
        void print_itemlist();

        typedef std::map<std::string, std::map< std::string, std::map<std::string, std::string>>> Itemlist;

    private:
        Itemlist itemlist;
};
#endif
