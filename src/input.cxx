#include <cstdio>
#include <cstring>
#include <map>
#include <string>
#include <sstream>
#include <fstream>
#include "input.h"

cinput::cinput()
{
  std::printf("Creating instance of object input\n");
}

cinput::~cinput()
{
  std::printf("Destroying instance of object input\n");
}

int cinput::readinifile()
{
  // construct input list before implementing file reading
  /*// setup Moser 180 case
  inputlist["grid"]["itot"] = "64";
  inputlist["grid"]["jtot"] = "64";
  inputlist["grid"]["ktot"] = "64";

  inputlist["grid"]["xsize"] = "6.28";
  inputlist["grid"]["ysize"] = "3.14";
  inputlist["grid"]["zsize"] = "2.";

  // set the time properties
  inputlist["time"]["runtime"] = "10000.";
  inputlist["time"]["cflmax" ] = "0.8";
  inputlist["time"]["adaptivestep" ] = "true";*/
  // end setup Moser case
 
  char inputline[256], temp1[256], block[256], item[256], lhs[256], rhs[256];

  // read the input file
  FILE *inputfile;
  inputfile = fopen("microhh.ini", "r");
  int n;
  bool blockset = false;
  int  nerrors  = 0;
 
  std::printf("Processing ini file\n");

  // check the cases: comments, empty line, block, value, rubbish
  while(std::fgets(inputline, 256, inputfile) != NULL)
  {
    // check for empty line
    n = std::sscanf(inputline, " %s ", temp1);
    if(n == 0) 
      continue; 

    // check for comments
    n = std::sscanf(inputline, " #%[^\n]", temp1);
    if(n > 0)
      continue;

    n = std::sscanf(inputline, " [%[^]]] ", temp1);
    if(n == 1)
    {
      n = std::sscanf(temp1, "%s %s", block);
      if(n == 1)
      {
        std::printf("Found block [%s]\n", block);
        blockset = true;
      }
      else
      {
        std::printf("ERROR in block specification [%s]\n", temp1);
        nerrors++;
      }
      continue;
    }
    // read items
    n = std::sscanf(inputline, "%[^=] = %[^\n]", temp1, rhs);
    if(n == 2)
    {
      n = std::sscanf(temp1, " %[a-zA-Z0-9_] %s", lhs);
      if(n == 1)
      {
        if(!blockset)
        {
          std::printf("ERROR item [?][%s] = \"%s\"\n", lhs, rhs);
          nerrors++;
          continue;
        }

        std::printf("Found item  [%s][%s] = \"%s\"\n", block, lhs, rhs);
        std::string blockstring(block);
        std::string itemstring(lhs);
        std::string valuestring(rhs);
        inputlist[blockstring][itemstring] = valuestring;
      }
      else
      {
        n = std::sscanf(inputline, "%[^=]", temp1);
        std::printf("ERROR item  [%s][%s]\n", block, temp1);
        nerrors++;
      }
    }

    // throw exception
    else
    {
      n = std::sscanf(inputline, "%[^\n]", temp1);
      if(n > 0)
      {
        std::printf("ERROR \"%s\" is illegal input\n", temp1);
        nerrors++;
      }
    }
  }
  fclose(inputfile);

  std::printf("Ini file has been processed with %d errors\n", nerrors);
  return nerrors;
}

// overloaded return functions
int cinput::getItem(int *value, std::string cat, std::string item)
{
  std::string inputstring = inputlist[cat][item];

  std::istringstream ss(inputstring);
  ss >> *value;

  return 0;
}

int cinput::getItem(double *value, std::string cat, std::string item)
{
  std::string inputstring = inputlist[cat][item];

  std::istringstream ss(inputstring);
  ss >> *value;

  return 0;
}

int cinput::getItem(bool *value, std::string cat, std::string item)
{
  std::string itemvalue = inputlist[cat][item];
  std::istringstream ss(itemvalue);
  if(itemvalue == "true" || itemvalue == "false")
    ss >> std::boolalpha >> *value;
  else if(itemvalue == "1" || itemvalue == "0")
    ss >> std::noboolalpha >> *value;

  return 0;
}

