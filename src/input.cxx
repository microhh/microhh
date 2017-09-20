#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include "input.h"
#include "convert.h"
#include "master.h"

namespace
{
    using namespace Convert;

    std::string get_item_string(const Input::Itemlist& itemlist,
                                const std::string& blockname,
                                const std::string& itemname,
                                const std::string& subitemname)
    {
        auto itblock = itemlist.find(blockname);
        if (itblock == itemlist.end())
        {
            std::string error_message = "Block \"" + blockname + "\" does not exist";
            throw std::runtime_error(error_message);
        }

        auto ititem = itblock->second.find(itemname);
        if (ititem == itblock->second.end())
        {
            std::string error_message = "Item \"" + itemname + "\" in group \"" + blockname + "\" does not exist";
            throw std::runtime_error(error_message);
        }

        // If the subitem is empty, it finds the default option "".
        auto itsubitem = ititem->second.find(subitemname);
        if (itsubitem == ititem->second.end())
        {
            itsubitem = ititem->second.find(std::string(""));
            if (itsubitem == ititem->second.end())
            {
                std::string error_message = "Subitem \"" + subitemname + "\" in item \"" + itemname + "\" does not exist";
                throw std::runtime_error(error_message);
            }
        }

        return itsubitem->second;
    }
}

Input::Input(Master& master, const std::string& file_name) : master(master)
{
    std::string blockname;

    // Read file and throw exception on error.
    std::ifstream infile;
    infile.open(file_name);
    if (!infile)
        throw std::runtime_error("Illegal file name");

    std::string line;

    while (std::getline(infile, line))
    {
        // Strip of the comments.
        std::vector<std::string> strings;
        boost::split(strings, line, boost::is_any_of("#"));

        // Keep part that is not comment.
        if (strings.size() >= 2)
            line = strings[0];

        // Strip of all the whitespace.
        boost::trim(line);

        // Split string on = character.
        strings.clear();
        boost::split(strings, line, boost::is_any_of("="));

        // In case of one element, check whether we have a block header.
        if (strings.size() == 1)
        {
            std::string header = strings[0];
            boost::trim(header);

            // If only an empty line remains, jump to the next line.
            if (header.empty())
                continue;

            // Store the block name and check for validity of string.
            // TODO: use .front() and back() C++11 once Intel supports it.
            if (header.at(0) == '[' && header.at(header.size()-1) == ']')
            {
                blockname = header.substr(1, header.size()-2);
                check_item(blockname);
            }
            else
            {
                std::string error_message = "Line \"" + line + "\" is illegal";
                throw std::runtime_error(error_message);
            }
        }
        // Read item.
        else if (strings.size() == 2)
        {
            if (blockname.empty())
                throw std::runtime_error("No block name found");

            std::string left = strings[0];
            std::string right = strings[1];
            boost::trim(left);
            boost::trim(right);

            // Check if suboptions are defined.
            std::string leftsub;
            std::size_t openpos = left.find_first_of("[");
            std::size_t closepos = left.find_last_of("]");
            if (openpos != std::string::npos && closepos != std::string::npos)
            {
                if (openpos < closepos)
                {
                    // Save the suboption and check string.
                    leftsub = left.substr(openpos+1, closepos-openpos-1);
                    check_item(leftsub);

                    // Strip of suboption.
                    left = left.substr(0, openpos);
                }
            }

            check_item(left);

            // Leave the checking of the right string for later
            // when the type is known.
            itemlist[blockname][left][leftsub] = right;
        }
        // Throw an error.
        else
        {
            std::string error_message = "Line \"" + line + "\" is illegal";
            throw std::runtime_error(error_message);
        }
    }
}

void Input::print_itemlist()
{
    // Print the list as a test.
    for (auto& b : itemlist)
        for (auto& i : b.second)
            for (auto& is : i.second)
                std::cout << b.first << "," << i.first << "," << is.first << "," << is.second << ";" << std::endl;
}

namespace
{
    template<typename T>
    T convert_value_to_item(const std::string& value)
    {
        std::istringstream ss(value);

        T item = get_item_from_stream<T>(ss);
        check_item<T>(item);

        return item;
    }
}

template<typename T>
T Input::get_item(const std::string& blockname,
                  const std::string& itemname,
                  const std::string& subitemname)
{
    std::string value = get_item_string(itemlist, blockname, itemname, subitemname);
    T item = convert_value_to_item<T>(value); 

    std::string itemout = "[" + blockname + "][" + itemname + "]";
    if (!subitemname.empty())
        itemout += "[" + subitemname + "]";

    std::ostringstream ss;
    ss << std::left << std::setw(30) << itemout << "= " 
       << std::right << std::setw(11) << std::setprecision(5) << std::boolalpha << item 
       << std::endl;

    master.print_message(ss);

    return item;
}

template<typename T>
T Input::get_item(const std::string& blockname,
                  const std::string& itemname,
                  const std::string& subitemname,
                  const T default_value)
{
    T item;
    std::string itemqualifier;

    try
    {
        std::string value = get_item_string(itemlist, blockname, itemname, subitemname);
        item = convert_value_to_item<T>(value); 
    }
    catch (std::runtime_error& e)
    {
        item = default_value;
        itemqualifier = "(default)";
    }

    std::string itemout = "[" + blockname + "][" + itemname + "]";
    if (!subitemname.empty())
        itemout += "[" + subitemname + "]";

    std::ostringstream ss;
    ss << std::left << std::setw(30) << itemout << "= " 
       << std::right << std::setw(11) << std::setprecision(5) << std::boolalpha << item 
       << "   " << itemqualifier << std::endl;

    master.print_message(ss);

    return item;
}

template<typename T>
std::vector<T> Input::get_list(const std::string& blockname,
                               const std::string& itemname,
                               const std::string& subitemname)
{
    std::string value = get_item_string(itemlist, blockname, itemname, subitemname);

    std::vector<std::string> listitems;
    boost::split(listitems, value, boost::is_any_of(","));

    std::vector<T> list;
    for (std::string itemstring : listitems)
    {
        std::istringstream ss(itemstring);
        T item = get_item_from_stream<T>(ss);
        check_item(item);
        list.push_back(item);
    }

    return list;
}

template<typename T>
std::vector<T> Input::get_list(const std::string& blockname,
                               const std::string& itemname,
                               const std::string& subitemname,
                               const std::vector<T> default_value)
{
    bool set_default_value = false;
    std::string value;

    try
    {
        value = get_item_string(itemlist, blockname, itemname, subitemname);
    }
    catch (std::runtime_error& e)
    {
        set_default_value = true;
        // itemqualifier = "(default)";
    }

    if (!set_default_value)
    {
        std::vector<std::string> listitems;
        boost::split(listitems, value, boost::is_any_of(","));

        std::vector<T> list;
        for (std::string itemstring : listitems)
        {
            std::istringstream ss(itemstring);
            T item = get_item_from_stream<T>(ss);
            check_item(item);
            list.push_back(item);
        }

        return list;
    }
    else
        return default_value;
}

// Explicitly instantiate templates.
template bool Input::get_item<bool>(const std::string&, const std::string&, const std::string&);
template int Input::get_item<int>(const std::string&, const std::string&, const std::string&);
template double Input::get_item<double>(const std::string&, const std::string&, const std::string&);
template float Input::get_item<float>(const std::string&, const std::string&, const std::string&);
template std::string Input::get_item<std::string>(const std::string&, const std::string&, const std::string&);

template bool Input::get_item<bool>(const std::string&, const std::string&, const std::string&, const bool);
template int Input::get_item<int>(const std::string&, const std::string&, const std::string&, const int);
template double Input::get_item<double>(const std::string&, const std::string&, const std::string&, const double);
template float Input::get_item<float>(const std::string&, const std::string&, const std::string&, const float);
template std::string Input::get_item<std::string>(const std::string&, const std::string&, const std::string&, const std::string);

template std::vector<int> Input::get_list<int>(const std::string&, const std::string&, const std::string&);
template std::vector<double> Input::get_list<double>(const std::string&, const std::string&, const std::string&);
template std::vector<float> Input::get_list<float>(const std::string&, const std::string&, const std::string&);
template std::vector<std::string> Input::get_list<std::string>(const std::string&, const std::string&, const std::string&);

template std::vector<int> Input::get_list<int>(const std::string&, const std::string&, const std::string&, const std::vector<int>);
template std::vector<double> Input::get_list<double>(const std::string&, const std::string&, const std::string&, const std::vector<double>);
template std::vector<float> Input::get_list<float>(const std::string&, const std::string&, const std::string&, const std::vector<float>);
template std::vector<std::string> Input::get_list<std::string>(const std::string&, const std::string&, const std::string&, const std::vector<std::string>);
