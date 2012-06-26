#include <cstdio>
#include <cstring>
#include <map>
#include "input.h"

cinput::cinput()
{
  std::printf("Creating instance of object input\n");
}

cinput::~cinput()
{
  std::printf("Destroying instance of object input\n");
}

int cinput::readinifile(std::string inputfilename)
{
  char inputline[256], temp1[256], block[256], lhs[256], rhs[256], dummy[256];

  // read the input file
  FILE *inputfile;
  inputfilename += ".ini";
  inputfile = fopen(inputfilename.c_str(), "r");

  if(inputfile == NULL)
  {
    std::printf("ERROR \"%s\" does not exist\n", inputfilename.c_str());
    return 1;
  }

  int n;
  bool blockset = false;
  int  nerrors  = 0;
 
  std::printf("Processing inifile \"%s\"\n", inputfilename.c_str());

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
      n = std::sscanf(temp1, "%s %s", block, dummy);
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
      n = std::sscanf(temp1, " %[a-zA-Z0-9_()] %s", lhs, dummy);
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

  std::printf("Inifile has been processed with %d errors\n", nerrors);
  return nerrors;
}

int cinput::readproffile(std::string inputfilename)
{
  char inputline[256], temp1[256], block[256], lhs[256], rhs[256], dummy[256];
  char *substring;
  int n;

  // read the input file
  FILE *inputfile;
  inputfilename += ".prof";
  inputfile = fopen(inputfilename.c_str(), "r");

  if(inputfile == NULL)
  {
    std::printf("ERROR \"%s\" does not exist\n", inputfilename.c_str());
    return 1;
  }

  std::printf("Processing proffile \"%s\"\n", inputfilename.c_str());
  int nvar = 0;

  // first find the header
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

    // read the header
    // read the first substring
    substring = std::strtok(inputline, " ,;\t\n");
    while(substring != NULL)
    {
      nvar++;

      if(!std::isalpha(substring[0]))
      {
        std::printf("ERROR \"%s\" is not a variable name\n", substring);
        return 1;
      }
        
      std::printf("Found variable \"%s\"\n", substring);

      // store the variable name

      // read the next substring
      substring = std::strtok(NULL, " ,;\t\n");
    }

    if(nvar == 0)
    {
      std::printf("ERROR no variable names in header\n");
      return 1;
    }

    // step out of the fgets loop
    break;
  }
  
  // second read the data
  // continue reading
  int ncols;
  double datavalue;

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

    // read the data
    ncols = 0;
    // read the first substring
    substring = std::strtok(inputline, " ,;\t\n");
    while(substring != NULL)
    {
      ncols++;

      // scan the line, while checking that the whole string has been read
      n = std::sscanf(inputline, " %lf %[^\n]", &datavalue, temp1);

      if(n != 1)
      {
        std::printf("ERROR \"%s\" is not a correct data value\n", substring);
        return 1;
      }

      // store the data
        
      // read the next substring
      substring = std::strtok(NULL, " ,;\t\n");
    }

    if(nvar == 0)
    {
      std::printf("ERROR no variable names in header\n");
      return 1;
    }
    else if(ncols != nvar)
    {
      std::printf("ERROR %d data columns, but %d defined variables\n", ncols, nvar);
      return 1;
    }
  }

  return 0;
}


int cinput::checkItemExists(std::string cat, std::string item)
{  
  inputmap::const_iterator it1 = inputlist.find(cat);

  bool readerror = false;

  if(it1 != inputlist.end())
  {
    inputmapinner::const_iterator it2 = it1->second.find(item);

    if(it2 == it1->second.end())
      readerror = true;
  }
  else
    readerror = true;

  if(readerror)
    return 1;

  return 0;
}

// overloaded return functions
int cinput::getItem(int *value, std::string cat, std::string item)
{
  if(checkItemExists(cat, item))
  {
    std::printf("ERROR [%s][%s] does not exist\n", cat.c_str(), item.c_str());
    return 1;
  }
  else
    if(checkItem(value, cat, item))
      return 1;
  
  return 0;
}

int cinput::getItem(int *value, std::string cat, std::string item, int def)
{
  if(checkItemExists(cat, item))
  {
    std::printf("WARNING [%s][%s] does not exist, default value of %d used\n", cat.c_str(), item.c_str(), def);
    *value = def;
    return 0;
  }
  else
    if(checkItem(value, cat, item))
      return 1;
  
  return 0;
}

int cinput::checkItem(int *value, std::string cat, std::string item)
{
  char inputstring[256], temp[256];
  std::strcpy(inputstring, inputlist[cat][item].c_str());

  int inputint;
  int n = std::sscanf(inputstring, " %d %[^\n] ", &inputint, temp);

  if(n == 1)
    *value = inputint;
  else
  {
    if(std::strcmp(inputstring,""))
    {
      std::printf("ERROR [%s][%s] = \"%s\" is not of type INT\n", cat.c_str(), item.c_str(), inputstring);
      return 1;
    }
  }

  return 0;
}

int cinput::getItem(double *value, std::string cat, std::string item)
{
  if(checkItemExists(cat, item))
  {
    std::printf("ERROR [%s][%s] does not exist\n", cat.c_str(), item.c_str());
    return 1;
  }
  else
    if(checkItem(value, cat, item))
      return 1;
  
  return 0;
}

int cinput::getItem(double *value, std::string cat, std::string item, double def)
{
  if(checkItemExists(cat, item))
  {
    std::printf("WARNING [%s][%s] does not exist, default value of %f used\n", cat.c_str(), item.c_str(), def);
    *value = def;
    return 0;
  }
  else
    if(checkItem(value, cat, item))
      return 1;
  
  return 0;
}

int cinput::checkItem(double *value, std::string cat, std::string item)
{
  char inputstring[256], temp[256];
  std::strcpy(inputstring, inputlist[cat][item].c_str());

  double inputdouble;
  int n = std::sscanf(inputstring, " %lf %[^\n] ", &inputdouble, temp);

  // catch the situation where a double is closed with a ".", which is not read by sscanf's %f
  if(n == 1 || (n == 2 && !std::strcmp(".", temp)))
    *value = inputdouble;
  else
  {
    if(std::strcmp(inputstring,""))
    {
      std::printf("ERROR [%s][%s] = \"%s\" is not of type DOUBLE\n", cat.c_str(), item.c_str(), inputstring);
      return 1;
    }
  }

  return 0;
}

int cinput::getItem(bool *value, std::string cat, std::string item)
{
  if(checkItemExists(cat, item))
  {
    std::printf("ERROR [%s][%s] does not exist\n", cat.c_str(), item.c_str());
    return 1;
  }
  else
    if(checkItem(value, cat, item))
      return 1;
  
  return 0;
}

int cinput::getItem(bool *value, std::string cat, std::string item, bool def)
{
  if(checkItemExists(cat, item))
  {
    std::printf("WARNING [%s][%s] does not exist, default value of %d used\n", cat.c_str(), item.c_str(), def);
    *value = def;
    return 0;
  }
  else
    if(checkItem(value, cat, item))
      return 1;
  
  return 0;
}

int cinput::checkItem(bool *value, std::string cat, std::string item)
{
  char inputstring[256], inputbool[256], temp[256];
  std::strcpy(inputstring, inputlist[cat][item].c_str());

  int n = std::sscanf(inputstring, " %s %[^\n] ", inputbool, temp);

  bool boolerror = false;

  if(n == 1)
  {
    if(std::strcmp("true", inputbool) == 0 || 
       std::strcmp("True", inputbool) == 0 || 
       std::strcmp("TRUE", inputbool) == 0 || 
       std::strcmp("1"   , inputbool) == 0 )
      *value = true;
    else if(std::strcmp("false", inputbool) == 0 ||
            std::strcmp("False", inputbool) == 0 || 
            std::strcmp("FALSE", inputbool) == 0 || 
            std::strcmp("0    ", inputbool) == 0 )
      *value = false;
    else 
      boolerror = true;
  }
  
  if(n != 1 || boolerror)
  {
    if(std::strcmp(inputstring,""))
    {
      std::printf("ERROR [%s][%s] = \"%s\" is not of type BOOL\n", cat.c_str(), item.c_str(), inputstring);
      return 1;
    }
  }

  return 0;
}

