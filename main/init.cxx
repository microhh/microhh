#include <cstdio>
#include "input.h"
#include "grid.h"
#include "fields.h"

int main()
{
  // INIT
  // read the input data
  cinput input;
  if(input.readinifile())
    return 1;

  // initialize the MPI interface
  // create the objects, read the inputdata
  cgrid grid;
  if(grid.readinifile(&input))
    return 1;
  
  cfields fields(&grid);

  // initialize the objects, allocate the required memory
  grid.initgrid();
  fields.initfields();

  // create the objects, fill the fields with data
  grid.creategrid();
  fields.createfields();

  // store the data
  fields.boundary();
  grid.save();
  fields.save(0);
  fields.save(0);
  fields.save(0);
  fields.save(0);
  fields.save(0);
  // END INIT
  return 0;
}

