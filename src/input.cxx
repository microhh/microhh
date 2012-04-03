#include <cstdio>
#include <map>
#include <sstream>
#include "input.h"

cinput::cinput()
{
  std::printf("Creating instance of object input\n");

  // construct input list before implementing file reading
  // setup Moser 180 case
  inputlist["grid"]["itot"] = "64";
  inputlist["grid"]["jtot"] = "64";
  inputlist["grid"]["ktot"] = "64";

  inputlist["grid"]["xsize"] = "6.28";
  inputlist["grid"]["ysize"] = "3.14";
  inputlist["grid"]["zsize"] = "2.";

  // set the time properties
  inputlist["time"]["runtime"] = "10000.";
  inputlist["time"]["cflmax" ] = "0.8";
  inputlist["time"]["adaptivestep" ] = "true";
  // end setup Moser case
}

cinput::~cinput()
{
  std::printf("Destroying instance of object input\n");
}

// overloaded return functions
int cinput::getItem(int *value, std::string cat, std::string item)
{
  std::stringstream ss(std::string(inputlist[cat][item]));
  ss >> *value;

  return 0;
}

int cinput::getItem(double *value, std::string cat, std::string item)
{
  std::stringstream ss(std::string(inputlist[cat][item]));
  ss >> *value;

  return 0;
}

int cinput::getItem(bool *value, std::string cat, std::string item)
{
  std::string itemvalue = inputlist[cat][item];
  std::stringstream ss(itemvalue);
  if(itemvalue == "true" || itemvalue == "false")
    ss >> std::boolalpha >> *value;
  else if(itemvalue == "1" || itemvalue == "0")
    ss >> std::noboolalpha >> *value;
  else
    return 1;

  return 0;
}

