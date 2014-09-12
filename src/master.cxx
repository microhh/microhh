#include <cstdarg>
#include <cstdio>
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
#include "master.h"

void cmaster::printMessage(const char *format, ...)
{
  if(mpiid == 0)
  {
    va_list args;
    va_start(args, format);
    std::vfprintf(stdout, format, args);
    va_end(args);
  }
}

void cmaster::printWarning(const char *format, ...)
{
  std::string warningstr("WARNING: ");
  warningstr += std::string(format);

  const char *warningformat = warningstr.c_str();

  if(mpiid == 0)
  {
    va_list args;
    va_start(args, format);
    std::vfprintf(stdout, warningformat, args);
    va_end(args);
  }
}

void cmaster::printError(const char *format, ...)
{
  std::string errorstr("ERROR: ");
  errorstr += std::string(format);

  const char *errorformat = errorstr.c_str();

  if(mpiid == 0)
  {
    va_list args;
    va_start(args, format);
    std::vfprintf(stdout, errorformat, args);
    va_end(args);
  }
}
