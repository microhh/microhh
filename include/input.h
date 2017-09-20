#ifndef INPUT
#define INPUT

#include <map>
#include <vector>

class Master;

class Input
{
    public:
        Input(Master&, const std::string&);
        template<typename T> T get_item(const std::string&, const std::string&, const std::string&);
        template<typename T> T get_item(const std::string&, const std::string&, const std::string&, const T);
        template<typename T> std::vector<T> get_list(const std::string&, const std::string&, const std::string&);
        template<typename T> std::vector<T> get_list(const std::string&, const std::string&, const std::string&, const std::vector<T>);
        // void print_itemlist();
        void print_unused_items();
        void flag_as_used(const std::string&, const std::string&, const std::string&);

        typedef std::map<std::string, std::map< std::string, std::map<std::string, std::pair<std::string, bool>>>> Itemlist;

    private:
        Master& master;
        Itemlist itemlist;
};
#endif
