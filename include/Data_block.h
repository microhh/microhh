#ifndef DATA_BLOCK
// #include <map>

class Data_block
{
    public:
        Data_block(const std::string&);
        template<typename T> void get_vector(std::vector<T>&,
                                             const std::string&,
                                             const size_t,
                                             const size_t,
                                             const size_t);
    private:
        std::map<std::string, std::vector<std::string>> data_series;
};
#endif

