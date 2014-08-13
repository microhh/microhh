#include "fields.h"
#include "grid.h"

int cfields::prepareGPU()
{
  const int nmemsize = grid->ncells*sizeof(double);
  return 0;
}
