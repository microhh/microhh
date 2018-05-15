#ifndef DATA_BLOCK
#include <map>
#include <vector>
class Master;

class Data_block
{
    public:
        Data_block(Master& master, const std::string&);
        template<typename T> void get_vector(std::vector<T>&,
                                             const std::string&,
                                             const size_t,
                                             const size_t,
                                             const size_t);
         std::vector<std::string> get_headers();
         int get_vector_length(const std::string&);
    private:
        Master& master;
        std::map<std::string, std::vector<std::string>> data_series;
        std::vector<std::string> header_items;

};
#endif
