#include <cstdarg>
#include <cstdio>
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
