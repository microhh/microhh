#include <cstdio>
#include "input.h"
#include "grid.h"
#include "fields.h"

int main()
{
  // create the class objects
  cinput  input;
  cgrid   grid;
  cfields fields(&grid);
  
  // read the input data and terminate on error
  if(input.readinifile() != 0)
    return 1;
  if(grid.readinifile(&input) != 0)
    return 1;
  
  // initialize the objects, allocate the required memory
  grid.initgrid();
  fields.initfields();

  // fill the fields with data
  grid.creategrid();
  fields.createfields();

  // store the data on disk
  fields.boundary();
  grid.save();
  fields.save(0);

  return 0;
}

