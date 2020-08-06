#ifndef DEFINE_BOOL_H
#define DEFINE_BOOL_H

#ifdef USE_CBOOL
#define BOOL_TYPE signed char
#else
#define BOOL_TYPE int
#endif

#endif
