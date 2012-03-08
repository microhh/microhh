#include <iostream>
#include "grid.h"
#include "fields.h"

fields::fields(grid *dnsgrid)
{
  std::cout << "Creating fields" << std::endl;
  // allocate memory for 3d arrays
  std::cout << "Allocating 6 x " << dnsgrid->ncells << std::endl;
  flow = new double[dnsgrid->ncells*6];

  // set pointers to correct location
  u  = &flow[dnsgrid->ncells*0];
  v  = &flow[dnsgrid->ncells*1];
  w  = &flow[dnsgrid->ncells*2];
  ut = &flow[dnsgrid->ncells*3];
  vt = &flow[dnsgrid->ncells*4];
  wt = &flow[dnsgrid->ncells*5];

  // set all values to 0
  for(int i=0; i<dnsgrid->ncells; i++)
    flow[i] = 0.;
}

fields::~fields()
{
  std::cout << "Deleting fields" << std::endl;
  delete[] flow;
}

