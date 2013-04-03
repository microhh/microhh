#include <cstdio>
#include <cstring>
#include <cctype>
#include <map>
#include "input.h"
#include "mpiinterface.h"
#include <algorithm>
#include <string>

cinput::cinput(cmpi *mpiin)
{
  mpi = mpiin;
}

cinput::~cinput()
{
}

int cinput::clear()
{
  inputlist.clear();
  proflist.clear();
  return 0;
}

int cinput::readinifile()
{
  char inputline[256], temp1[256], block[256], lhs[256], rhs[256], dummy[256], element[256];

  // read the input file
  FILE *inputfile;
  std::string inputfilename = mpi->simname + ".ini";

  if(mpi->mpiid == 0)
  {
    inputfile = fopen(inputfilename.c_str(), "r");
    if(inputfile == NULL)
    {
      std::printf("ERROR \"%s\" does not exist\n", inputfilename.c_str());
      return 1;
    }
  }

  int n;
  bool blockset = false;
  int nerrors = 0;
  int nlines  = 0;
  int nline;

  if(mpi->mpiid == 0)
  {
    std::printf("Processing inifile \"%s\"\n", inputfilename.c_str());
    while(std::fgets(inputline, 256, inputfile) != NULL)
      nlines++;
    std::printf("Inifile contains %d lines\n", nlines);
    rewind(inputfile);
  }
  mpi->broadcast(&nlines, 1);

  // check the cases: comments, empty line, block, value, rubbish
  for(int nn=0; nn<nlines; nn++)
  {
    nline = nn+1;
    if(mpi->mpiid == 0)
    {
      // fetch a line and broadcast it
      std::fgets(inputline, 256, inputfile);
    }
    mpi->broadcast(inputline, 256);

    for(int i = 0; inputline[i] != '\0'; i++)
    {
      inputline[i] = tolower(inputline[i]);
    }

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
        blockset = true;
      }
      else
      {
        if(mpi->mpiid == 0) std::printf("ERROR line %d: illegal block specification [%s]\n", nline, temp1);
        return 1;
      }
      continue;
    }
    // read items
    n = std::sscanf(inputline, "%[^=] = %[^\n]", temp1, rhs);
    if(n == 2)
    {
      
      n = std::sscanf(temp1, " %[a-zA-Z0-9_()][%[^]]] %s", lhs, element,dummy);
      if(n <= 2)
      {
        if(!blockset)
        {
          if(mpi->mpiid == 0) std::printf("ERROR line %d: illegal item [?][%s] = \"%s\"\n", nline, lhs, rhs);
          nerrors++;
          return 1;
        }

        if(n ==1)
        {
          std::strcpy(element,"default");
        }
        std::string blockstring(block);
        std::string itemstring(lhs);
        std::string elementstring(element);
        std::string valuestring(rhs);
        if(checkItemExists(blockstring,itemstring,elementstring))
        {
          inputlist[blockstring][itemstring][elementstring].data   = valuestring;
          inputlist[blockstring][itemstring][elementstring].isused = false;
        }
        else
        {
          if(mpi->mpiid == 0) std::printf("ERROR line %d: Item [%s][%s][%s] defined for the second time\n", nline, block, lhs, element);
          return 1;
        }
      }
      else
      {
        n = std::sscanf(inputline, "%[^=]", temp1);
        if(mpi->mpiid == 0) std::printf("ERROR line %d: illegal item  [%s][%s]\n", nline, block, temp1);
        nerrors++;
      }
    }

    // throw exception
    else
    {
      n = std::sscanf(inputline, "%[^\n]", temp1);
      if(n > 0)
      {
        if(mpi->mpiid == 0) std::printf("ERROR line %d: \"%s\" is illegal input\n", nline, temp1);
        nerrors++;
      }
    }
  }

  if(mpi->mpiid == 0)
  {
    std::printf("Inifile has been processed with %d errors\n", nerrors);
    fclose(inputfile);
  }

  return nerrors;
}

int cinput::readproffile()
{
  char inputline[256], temp1[256];
  char *substring;
  int n;

  // read the input file
  FILE *inputfile;
  std::string inputfilename = mpi->simname + ".prof";

  if(mpi->mpiid == 0)
  {
    inputfile = fopen(inputfilename.c_str(), "r");
    if(inputfile == NULL)
    {
      std::printf("ERROR \"%s\" does not exist\n", inputfilename.c_str());
      return 1;
    }
  }

  int nlines = 0;
  int nline;
  int nvar   = 0;
  // std::vector<std::string> varnames;

  if(mpi->mpiid == 0)
  {
    std::printf("Processing proffile \"%s\"\n", inputfilename.c_str());
    while(std::fgets(inputline, 256, inputfile) != NULL)
      nlines++;
    std::printf("Inifile contains %d lines\n", nlines);
    rewind(inputfile);
  }
  mpi->broadcast(&nlines, 1);

  int nn;

  // first find the header
  for(nn=0; nn<nlines; nn++)
  {
    nline = nn+1;
    if(mpi->mpiid == 0)
    {
      // fetch a line and broadcast it
      std::fgets(inputline, 256, inputfile);
    }
    mpi->broadcast(inputline, 256);

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
        if(mpi->mpiid == 0)
        {
          std::printf("ERROR at line %d: \"%s\" is not a variable name\n", nline, substring);
          fclose(inputfile);
        }
        return 1;
      }
        
      if(mpi->mpiid == 0) std::printf("Found variable \"%s\"\n", substring);

      // temporarily store the variable name
      varnames.push_back(std::string(substring));

      // read the next substring
      substring = std::strtok(NULL, " ,;\t\n");
    }

    if(nvar == 0)
    {
      if(mpi->mpiid == 0)
      {
        std::printf("ERROR no variable names in header\n");
        fclose(inputfile);
      }
      return 1;
    }

    // step out of the fgets loop
    break;
  }
  
  // second read the data
  // continue reading
  int ncols;
  double datavalue;

  // std::vector<double> varvalues;

  // continue the loop from the exit value of nn
  for(nn++; nn<nlines; nn++)
  {
    nline = nn+1;
    if(mpi->mpiid == 0)
    {
      // fetch a line and broadcast it
      std::fgets(inputline, 256, inputfile);
    }
    mpi->broadcast(inputline, 256);

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
    varvalues.clear();
    // read the first substring
    substring = std::strtok(inputline, " ,;\t\n");
    while(substring != NULL)
    {
      ncols++;

      // scan the line, while checking that the whole string has been read
      n = std::sscanf(substring, " %lf %[^\n]", &datavalue, temp1);

      if(n != 1)
      {
        if(mpi->mpiid == 0)
        {
          std::printf("ERROR line %d: \"%s\" is not a correct data value\n", nline, substring);
          fclose(inputfile);
        }
        return 1;
      }

      // temporarily store the data
      varvalues.push_back(datavalue);
        
      // read the next substring
      substring = std::strtok(NULL, " ,;\t\n");
    }

    if(ncols != nvar)
    {
      if(mpi->mpiid == 0)
      {
        std::printf("ERROR line %d: %d data columns, but %d defined variables\n", nline, ncols, nvar);
        fclose(inputfile);
      }
      return 1;
    }

    // store the data
    for(n=0; n<nvar; n++)
      proflist[varnames[n]].push_back(varvalues[n]);
  }

  if(mpi->mpiid == 0)
    fclose(inputfile);

  return 0;
}


int cinput::checkItemExists(std::string cat, std::string item, std::string el)
{  
  inputmap::const_iterator it1 = inputlist.find(cat);

  bool readerror = false;

  if(it1 != inputlist.end())
  {
    inputmap2d::const_iterator it2 = it1->second.find(item);

    if(it2 != it1->second.end())
    {
      inputmap1d::const_iterator it3 = it2->second.find(el);
      if(it3 == it2->second.end())
        readerror = true;
    }
    else
      readerror = true;
  }
  else
    readerror = true;

  if(readerror)
    return 1;

  return 0;
}

// overloaded return functions
int cinput::getItem(int *value, std::string cat, std::string item, std::string el)
{
  char cwhy[256]="";
  char cwho[256]="";

  strcat(cwho, "[");
  strcat(cwho,cat.c_str());
  strcat(cwho,"][");
  strcat(cwho,item.c_str());
  strcat(cwho,"]");

  if(!el.empty())
  {
    strcat(cwho,"[");
    strcat(cwho,el.c_str());
    strcat(cwho,"]");
    if(!checkItemExists(cat, item, el))
    {
      if(checkItem(value, cat, item, el))
        return 1;
      strcpy(cwhy, "element specific value");
    }
  }
  if(!strcmp(cwhy,""))
  {
    if(checkItemExists(cat, item))
    {
      if(mpi->mpiid == 0) std::printf("ERROR [%s][%s] does not exist\n", cat.c_str(), item.c_str());
      return 1;
    }
    else
    {
      if(checkItem(value, cat, item))
        return 1;
      strcpy(cwhy, "global value");
    }
  }
  strncat(cwho,"                          ",30-strlen(cwho));
  if(mpi->mpiid == 0) std::printf("%s= %9d   (%s)\n", cwho, *value, cwhy);
  return 0;
}

int cinput::getItem(int *value, std::string cat, std::string item, std::string el, int def)
{
  char cwhy[256]="";
  char cwho[256]="";

  strcat(cwho, "[");
  strcat(cwho,cat.c_str());
  strcat(cwho,"][");
  strcat(cwho,item.c_str());
  strcat(cwho,"]");

  if(!el.empty())
  {
    strcat(cwho,"[");
    strcat(cwho,el.c_str());
    strcat(cwho,"]");
    if(!checkItemExists(cat, item, el))
    {
      if(checkItem(value, cat, item, el))
        return 1;
      strcpy(cwhy, "element specific value");
    }
  }
  if(!strcmp(cwhy,""))
  {
    if(checkItemExists(cat, item))
    {
      *value = def;
      strcpy(cwhy, "default value");
    }
    else
    {
      if(checkItem(value, cat, item))
        return 1;
      strcpy(cwhy, "global value");
    }
  }
  strncat(cwho,"                          ",30-strlen(cwho));
  if(mpi->mpiid == 0) std::printf("%s= %9d   (%s)\n", cwho, *value, cwhy);
  return 0;
}

int cinput::checkItem(int *value, std::string cat, std::string item, std::string el)
{
  char inputstring[256], temp[256];
  std::strcpy(inputstring, inputlist[cat][item][el].data.c_str());

  int inputint;
  int n = std::sscanf(inputstring, " %d %[^\n] ", &inputint, temp);

  if(n == 1)
    *value = inputint;
  else
  {
    if(std::strcmp(inputstring,""))
    {
      if (el == "default")
      {
        if(mpi->mpiid == 0) std::printf("ERROR [%s][%s] = \"%s\" is not of type INT\n", cat.c_str(), item.c_str(), inputstring);
      }
      else
      {
        if(mpi->mpiid == 0) std::printf("ERROR [%s][%s][%s] = \"%s\" is not of type INT\n", cat.c_str(), item.c_str(), el.c_str(), inputstring);
      }
      return 1;
    }
  }
  inputlist[cat][item][el].isused = true;

  return 0;
}

int cinput::getItem(double *value, std::string cat, std::string item, std::string el)
{
  char cwhy[256]="";
  char cwho[256]="";

  strcat(cwho, "[");
  strcat(cwho,cat.c_str());
  strcat(cwho,"][");
  strcat(cwho,item.c_str());
  strcat(cwho,"]");

  if(!el.empty())
  {
    strcat(cwho,"[");
    strcat(cwho,el.c_str());
    strcat(cwho,"]");
    if(!checkItemExists(cat, item, el))
    {
      if(checkItem(value, cat, item, el))
        return 1;
      strcpy(cwhy, "element specific value");
    }
  }
  if(!strcmp(cwhy,""))
  {
    if(checkItemExists(cat, item))
    {
      if(mpi->mpiid == 0) std::printf("ERROR [%s][%s] does not exist\n", cat.c_str(), item.c_str());
      return 1;
    }
    else
    {
      if(checkItem(value, cat, item))
        return 1;
      strcpy(cwhy, "global value");
    }
  }
  strncat(cwho,"                          ",30-strlen(cwho));
  if(mpi->mpiid == 0) std::printf("%s= %9.4G   (%s)\n", cwho, *value, cwhy);
  return 0;
}

int cinput::getItem(double *value, std::string cat, std::string item, std::string el, double def)
{
  char cwhy[256]="";
  char cwho[256]="";

  strcat(cwho, "[");
  strcat(cwho,cat.c_str());
  strcat(cwho,"][");
  strcat(cwho,item.c_str());
  strcat(cwho,"]");

  if(!el.empty())
  {
    strcat(cwho,"[");
    strcat(cwho,el.c_str());
    strcat(cwho,"]");
    if(!checkItemExists(cat, item, el))
    {
      if(checkItem(value, cat, item, el))
        return 1;
      strcpy(cwhy, "element specific value");
    }
  }
  if(!strcmp(cwhy,""))
  {
    if(checkItemExists(cat, item))
    {
      *value = def;
      strcpy(cwhy, "default value");
    }
    else
    {
      if(checkItem(value, cat, item))
        return 1;
      strcpy(cwhy, "global value");
    }
  }
  strncat(cwho,"                          ",30-strlen(cwho));
  if(mpi->mpiid == 0) std::printf("%s= %9.4G   (%s)\n", cwho,*value, cwhy);
  return 0;
}

int cinput::checkItem(double *value, std::string cat, std::string item, std::string el)
{
  char inputstring[256], temp[256];
  std::strcpy(inputstring, inputlist[cat][item][el].data.c_str());

  double inputdouble;
  int n = std::sscanf(inputstring, " %lf %[^\n] ", &inputdouble, temp);
  // catch the situation where a double is closed with a ".", which is not read by sscanf's %f
  if(n == 1 || (n == 2 && !std::strcmp(".", temp)))
    *value = inputdouble;
  else
  {
    if(std::strcmp(inputstring,""))
    {
      if(el == "default")
      {
        if(mpi->mpiid == 0) std::printf("ERROR [%s][%s] = \"%s\" is not of type DOUBLE\n", cat.c_str(), item.c_str(), inputstring);
      }
      else
      {
        if(mpi->mpiid == 0) std::printf("ERROR [%s][%s][%s] = \"%s\" is not of type DOUBLE\n", cat.c_str(), item.c_str(), el.c_str(), inputstring);
      }
      return 1;
    }
  }
  inputlist[cat][item][el].isused = true;

  return 0;
}

int cinput::getItem(bool *value, std::string cat, std::string item, std::string el)
{
  char cwhy[256]="";
  char cwho[256]="";

  strcat(cwho, "[");
  strcat(cwho,cat.c_str());
  strcat(cwho,"][");
  strcat(cwho,item.c_str());
  strcat(cwho,"]");

  if(!el.empty())
  {
    strcat(cwho,"[");
    strcat(cwho,el.c_str());
    strcat(cwho,"]");
    if(!checkItemExists(cat, item, el))
    {
      if(checkItem(value, cat, item, el))
        return 1;
      strcpy(cwhy, "element specific value");
    }
  }
  if(!strcmp(cwhy,""))
  {
    if(checkItemExists(cat, item))
    {
      if(mpi->mpiid == 0) std::printf("ERROR [%s][%s] does not exist\n", cat.c_str(), item.c_str());
      return 1;
    }
    else
    {
      if(checkItem(value, cat, item))
        return 1;
      strcpy(cwhy, "global value");
    }
  }
  strncat(cwho,"                          ",30-strlen(cwho));
  if(mpi->mpiid == 0) std::printf("%s= %s   (%s)\n", cwho,(*value) ? "     true" : "    false", cwhy);
  return 0;
}

int cinput::getItem(bool *value, std::string cat, std::string item, std::string el, bool def)
{
  char cwhy[256]="";
  char cwho[256]="";

  strcat(cwho, "[");
  strcat(cwho,cat.c_str());
  strcat(cwho,"][");
  strcat(cwho,item.c_str());
  strcat(cwho,"]");

  if(!el.empty())
  {
    strcat(cwho,"[");
    strcat(cwho,el.c_str());
    strcat(cwho,"]");
    if(!checkItemExists(cat, item, el))
    {
      if(checkItem(value, cat, item, el))
        return 1;
      strcpy(cwhy, "element specific value");
    }
  }
  if(!strcmp(cwhy,""))
  {
    if(checkItemExists(cat, item))
    {
      *value = def;
      strcpy(cwhy, "default value");
    }
    else
    {
      if(checkItem(value, cat, item))
        return 1;
      strcpy(cwhy, "global value");
    }
  }
  strncat(cwho,"                          ",30-strlen(cwho));
  if(mpi->mpiid == 0) std::printf("%s= %s   (%s)\n", cwho,(*value) ? "     true" : "    false", cwhy);
  return 0;
}

int cinput::checkItem(bool *value, std::string cat, std::string item, std::string el)
{
  char inputstring[256], inputbool[256], temp[256];
  std::strcpy(inputstring, inputlist[cat][item][el].data.c_str());

  int n = std::sscanf(inputstring, " %s %[^\n] ", inputbool, temp);

  bool boolerror = false;

  if(n == 1)
  {
    if(std::strcmp("true", inputbool) == 0 ||
       std::strcmp("1"   , inputbool) == 0 )
      *value = true;
    else if(std::strcmp("false", inputbool) == 0 ||
            std::strcmp("0    ", inputbool) == 0 )
      *value = false;
    else 
      boolerror = true;
  }
  
  if(n != 1 || boolerror)
  {
    if(std::strcmp(inputstring,""))
    {
      if (el == "default")
      {
        if(mpi->mpiid == 0) std::printf("ERROR [%s][%s] = \"%s\" is not of type BOOL\n", cat.c_str(), item.c_str(), inputstring);
      }
      else
      {
        if(mpi->mpiid == 0) std::printf("ERROR [%s][%s][%s] = \"%s\" is not of type BOOL\n", cat.c_str(), item.c_str(), el.c_str(), inputstring);
      }
      return 1;
    }
  }
  inputlist[cat][item][el].isused = true;

  return 0;
}

int cinput::getItem(std::string *value, std::string cat, std::string item, std::string el)
{
  char cwhy[256]="";
  char cwho[256]="";

  strcat(cwho, "[");
  strcat(cwho,cat.c_str());
  strcat(cwho,"][");
  strcat(cwho,item.c_str());
  strcat(cwho,"]");

  if(!el.empty())
  {
    strcat(cwho,"[");
    strcat(cwho,el.c_str());
    strcat(cwho,"]");
    if(!checkItemExists(cat, item, el))
    {
      if(checkItem(value, cat, item, el))
        return 1;
      strcpy(cwhy, "element specific value");
    }
  }
  if(!strcmp(cwhy,""))
  {
    if(checkItemExists(cat, item))
    {
      if(mpi->mpiid == 0) std::printf("ERROR [%s][%s] does not exist\n", cat.c_str(), item.c_str());
      return 1;
    }
    else
    {
      if(checkItem(value, cat, item))
        return 1;
      strcpy(cwhy, "global value");
    }
  }
  strncat(cwho,"                          ",30-strlen(cwho));
  std::printf("%s= %9s   (%s)\n", cwho,value->c_str(), cwhy);
  return 0;
}

int cinput::getItem(std::string *value, std::string cat, std::string item, std::string el, std::string def)
{
  char cwhy[256]="";
  char cwho[256]="";

  strcat(cwho, "[");
  strcat(cwho,cat.c_str());
  strcat(cwho,"][");
  strcat(cwho,item.c_str());
  strcat(cwho,"]");

  if(!el.empty())
  {
    strcat(cwho,"[");
    strcat(cwho,el.c_str());
    strcat(cwho,"]");
    if(!checkItemExists(cat, item, el))
    {
      if(checkItem(value, cat, item, el))
        return 1;
      strcpy(cwhy, "element specific value");
    }
  }
  if(!strcmp(cwhy,""))
  {
    if(checkItemExists(cat, item))
    {
      *value = def;
      strcpy(cwhy, "default value");
    }
    else
    {
      if(checkItem(value, cat, item))
        return 1;
      strcpy(cwhy, "global value");
    }
  }
  strncat(cwho,"                          ",30-strlen(cwho));
  std::printf("%s= %9s   (%s)\n", cwho, value->c_str(), cwhy);
  return 0;
}

int cinput::checkItem(std::string *value, std::string cat, std::string item, std::string el)
{
  std::string inputstring, dummy;
  inputstring = inputlist[cat][item][el].data;

  int n = std::sscanf(inputstring.c_str(), "%s %s",value->c_str(), dummy.c_str());

  if(n == 1)
    *value = inputstring;
  else
  {
    if(!inputstring.empty())
    {
      if (el == "default")
      {
        if(mpi->mpiid == 0) std::printf("ERROR [%s][%s] = \"%s\" is not of type STRING\n", cat.c_str(), item.c_str(), inputstring.c_str());
      }
      else
      {
        if(mpi->mpiid == 0) std::printf("ERROR [%s][%s][%s] = \"%s\" is not of type STRING\n", cat.c_str(), item.c_str(), el.c_str(), inputstring.c_str());
      }
      return 1;
    }
  }
  inputlist[cat][item][el].isused = true;

  return 0;
}

// list retrieval function
int cinput::getItem(std::vector<std::string> *value, std::string cat, std::string item, std::string def)
{
  char cwho[256]="";

  strcat(cwho, "[");
  strcat(cwho,cat.c_str());
  strcat(cwho,"][");
  strcat(cwho,item.c_str());
  strcat(cwho,"]");
  strncat(cwho,"                          ",30-strlen(cwho));
  if(checkItemExists(cat, item))
  {
    std::printf("%s   NOT FOUND\n",cwho);
  }
  else
  {
    if(checkItem(value, cat, item))
      return 1;
    std::printf("%s= ",cwho);
    for(std::vector<std::string>::iterator it = value->begin(); it !=value->end()-1; ++it)
    {
      std::printf("%s, ",it->c_str());
    }
    std::printf("%s\n",(value->end()-1)->c_str());
  }

  return 0;
}

int cinput::checkItem(std::vector<std::string> *value, std::string cat, std::string item, std::string el)
{
  std::string inputstring, dummy;
  inputstring = inputlist[cat][item][el].data;

  char temp1[256];
  char * temp2;
  std::strcpy(temp1, inputstring.c_str());

  // first, split string on the delimiter
  temp2 = std::strtok(temp1, ",");

  while(temp2 != NULL)
  {
    // read in the string part in temp1
    int n = std::sscanf(temp2, "%s %s", temp1, dummy.c_str());

    // store the contents in the vector, or throw exception
    if(n == 1)
      value->push_back(temp1);
    else
    {
      if(!inputstring.empty())
      {
        if (el == "default")
          if(mpi->mpiid == 0) std::printf("ERROR [%s][%s] = \"%s\" is not a list of type STRING\n", cat.c_str(), item.c_str(), inputstring.c_str());
        else
          if(mpi->mpiid == 0) std::printf("ERROR [%s][%s][%s] = \"%s\" is not a list of type STRING\n", cat.c_str(), item.c_str(), el.c_str(), inputstring.c_str());
        // empty the vector
        value->clear();
        return 1;
      }
    }

    // retrieve the next raw substring
    temp2 = std::strtok(NULL, ",");
  }
  inputlist[cat][item][el].isused = true;

  return 0;
}

int cinput::printUnused()
{
  for(inputmap::iterator it1=inputlist.begin(); it1!=inputlist.end(); ++it1)
  {
    for(inputmap2d::iterator it2=it1->second.begin(); it2!=it1->second.end(); ++it2)
    {
      for(inputmap1d::iterator it3=it2->second.begin(); it3!=it2->second.end(); ++it3)
      {
        if(!it3->second.isused)
        {
          if(it3->first == "default")
          {
            if(mpi->mpiid == 0) std::printf("WARNING [%s][%s] = \"%s\" is not used\n", it1->first.c_str(), it2->first.c_str(), it3->second.data.c_str());
          }
          else
          {
            if(mpi->mpiid == 0) std::printf("WARNING [%s][%s][%s] = \"%s\" is not used\n", it1->first.c_str(), it2->first.c_str(), it3->first.c_str(), it3->second.data.c_str());
          }
        }
      }
    }
  }
  return 0;
}

int cinput::getProf(double *data, std::string varname, int kmaxin)
{
  profmap::const_iterator it = proflist.find(varname);

  if(it != proflist.end())
  {
    int profsize = proflist[varname].size();
    if(profsize < kmaxin)
    {
      if(mpi->mpiid == 0) std::printf("ERROR only %d of %d levels can be read for variable \"%s\"\n", profsize, kmaxin, varname.c_str());
      return 1;
    }
    if(profsize > kmaxin)
      if(mpi->mpiid == 0) std::printf("WARNING %d is larger than the number of grid points %d for variable \"%s\"\n", profsize, kmaxin, varname.c_str());

    for(int k=0; k<kmaxin; k++)
      data[k] = proflist[varname][k];

    if(mpi->mpiid == 0) std::printf("Variable \"%s\" has been read from the input\n", varname.c_str());
  }
  else
  { 
    if(mpi->mpiid == 0) std::printf("WARNING no profile data for variable \"%s\", values set to zero\n", varname.c_str());
    for(int k=0; k<kmaxin; k++)
      data[k] = 0.;
  }

  return 0;
}

