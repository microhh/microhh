#ifndef CONVERT
namespace Convert
{
    template<typename T>
    inline void check_item(const T& t) {}
    
    template<>
    inline void check_item(const std::string& s)
    {
        // Check whether string is empty.
        if (s.empty())
            throw std::runtime_error("Illegal string");
    
        // Return string if all characters are alphanumeric.
        if (find_if(s.begin(), s.end(), [](const char c) { return !std::isalnum(c); }) == s.end())
            return;
        else
            throw std::runtime_error("Illegal string: " + s);
    }

    template<typename T>
    inline T get_item_from_stream(std::istringstream& ss)
    {
        // Read the item from the stringstream, operator >> trims automatically.
        T item;
        if (!(ss >> item))
            throw std::runtime_error("Item does not match type");
    
        // Check whether stringstream is empty, if not type is incorrect.
        std::string dummy;
        if (ss >> dummy)
            throw std::runtime_error("Item does not match type");
    
        return item;
    }
}
#endif
