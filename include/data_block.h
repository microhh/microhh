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
    private:
        Master& master;
        std::map<std::string, std::vector<std::string>> data_series;
};
#endif

