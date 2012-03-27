#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "timeloop.h"
#include "advec.h"
#include "diff.h"
#include "force.h"
#include "pres.h"
#include "timeint.h"

int main()
{
  // INIT
  // initialize the MPI interface
  // create the objects, read the inputdata
  cgrid grid;
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
  fields.u->save(0);
  fields.v->save(0);
  fields.w->save(0);
  fields.p->save(0);
  fields.s->save(0);
  // END INIT
  return 0;
}

