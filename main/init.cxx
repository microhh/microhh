#include <cstdio>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "pres.h"

int main()
{
  // create the class objects
  cinput  input;
  cgrid   grid;
  cfields fields(&grid);
  cpres   pres  (&grid, &fields);
  
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
  if(grid.save())
    return 1;
  if(fields.save(0))
    return 1;
  if(pres.save())
    return 1;

  return 0;
}

