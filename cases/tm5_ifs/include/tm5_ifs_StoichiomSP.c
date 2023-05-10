/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Sparse Stoichiometric Data Structures File                       */
/*                                                                  */
/* Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor  */
/*       (http://www.cs.vt.edu/~asandu/Software/KPP)                */
/* KPP is distributed under GPL, the general public licence         */
/*       (http://www.gnu.org/copyleft/gpl.html)                     */
/* (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa           */
/* (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech            */
/*     With important contributions from:                           */
/*        M. Damian, Villanova University, USA                      */
/*        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany */
/*                                                                  */
/* File                 : tm5_ifs_StoichiomSP.c                     */
/* Time                 : Tue Aug 31 13:57:36 2021                  */
/* Working directory    : /home/WUR/krol005/kpp/examples            */
/* Equation file        : tm5_ifs.kpp                               */
/* Output root filename : tm5_ifs                                   */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tm5_ifs_Parameters.h"
#include "tm5_ifs_Global.h"
#include "tm5_ifs_Sparse.h"



/* Row-compressed sparse data for the Jacobian of reaction products JVRP */
 /* Beginning of rows in JVRP */

  int  CROW_JVRP[] = {
       0,  2,  4,  6,  7,  9, 10, 12, 14, 16, 17, 19,
      21, 23, 25, 27, 29, 31, 33, 35, 36, 37, 38, 39,
      40, 41, 43, 44, 45, 47, 49, 51, 52, 54, 55, 57,
      57, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66 }; 

 /* Column indices in JVRP */

  int  ICOL_JVRP[] = {
      11, 15, 10, 11, 10, 15, 10,  2, 15, 15, 11, 12,
      10, 12, 14, 15, 15,  9, 12,  9, 10,  9, 10,  7,
      15,  7, 15,  3, 15,  6, 15,  8, 11,  8, 15, 11,
      14,  6,  6,  7,  2,  5, 14,  4,  4, 11, 14, 12,
      13, 13, 14,  1, 10, 13,  1,  8, 13, 11, 12, 14,
       5,  0,  7,  2,  6,  4 }; 

 /* Row indices in JVRP */

  int  IROW_JVRP[] = {
       0,  0,  1,  1,  2,  2,  3,  4,  4,  5,  6,  6,
       7,  7,  8,  8,  9, 10, 10, 11, 11, 12, 12, 13,
      13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19,
      20, 21, 22, 23, 24, 25, 25, 26, 27, 28, 28, 29,
      29, 30, 30, 31, 32, 32, 33, 34, 34, 37, 38, 39,
      40, 41, 42, 43, 44, 45 }; 




/*  Stoichiometric Matrix in Compressed Column Sparse Format        */

 /* Beginning of columns in STOICM */

  int  CCOL_STOICM[] = {
       0,  3,  6,  8, 10, 13, 15, 18, 22, 25, 27, 32,
      35, 38, 42, 44, 47, 51, 56, 60, 62, 65, 68, 70,
      74, 76, 79, 82, 85, 88, 91, 94, 97,100,102,107,
     108,109,110,111,112,113,114,115,116,117,118 }; 

 /* Row indices in STOICM */

  int  IROW_STOICM[] = {
      10, 11, 15, 10, 11, 15, 10, 15,  2, 10,  2, 10,
      15, 10, 15, 11, 12, 14, 10, 12, 14, 15,  0, 14,
      15,  9, 15,  6,  9, 10, 12, 14,  7,  9, 10,  6,
       9, 10,  6,  7,  9, 15,  7, 15,  3, 10, 15,  3,
       6, 10, 15,  5,  8,  9, 11, 15,  5,  8,  9, 15,
      11, 15, 11, 12, 14,  3,  6, 10,  3,  6,  6,  7,
      10, 15,  2, 15,  4,  5, 14,  4,  5, 14,  4,  5,
      14, 11, 13, 14, 12, 13, 14,  1, 13, 14,  1, 13,
      14,  0, 10, 13,  0,  1,  5,  8,  9, 13, 14,  8,
      12, 11, 12, 14,  5,  0,  7,  2,  6,  4 }; 

 /* Column indices in STOICM */

  int  ICOL_STOICM[] = {
       0,  0,  0,  1,  1,  1,  2,  2,  3,  3,  4,  4,
       4,  5,  5,  6,  6,  6,  7,  7,  7,  7,  8,  8,
       8,  9,  9, 10, 10, 10, 10, 10, 11, 11, 11, 12,
      12, 12, 13, 13, 13, 13, 14, 14, 15, 15, 15, 16,
      16, 16, 16, 17, 17, 17, 17, 17, 18, 18, 18, 18,
      19, 19, 20, 20, 20, 21, 21, 21, 22, 22, 23, 23,
      23, 23, 24, 24, 25, 25, 25, 26, 26, 26, 27, 27,
      27, 28, 28, 28, 29, 29, 29, 30, 30, 30, 31, 31,
      31, 32, 32, 32, 33, 33, 34, 34, 34, 34, 34, 35,
      36, 37, 38, 39, 40, 41, 42, 43, 44, 45 }; 

 /* Stoichiometric Matrix in compressed column format */

  double  STOICM[] = {
         1,   -1,   -1,   -1,   -1,    1,   -1,   -1,
         1,   -2,   -1,    1,   -1,    1,   -1,   -1,
        -1,    1,   -1,   -1,    1,    1,    1,   -1,
        -1,    1,   -1,    1,   -1,    1,   -1,    1,
         1,   -1,   -1,    1,   -1,   -1,  0.4,   -1,
       0.6, -0.6,   -1,-0.77,   -1,    1,   -1,    1,
        -1,    1,   -1,  0.1,   -1, 0.75,   -1,  0.8,
       0.1,   -1,    1,   -1,   -1,    2,    1,    1,
        -1,    1,   -1,    2,    1,   -1,    1,   -1,
         1,    1,   -1,    2,    1,   -1,   -1,   -1,
         1,    1,   -1,    1,    1,   -1,    1,   -1,
        -1,   -1,    2,    1,   -1,   -1,   -1,    1,
         1,    1,   -1,   -1,    2,   -1,  0.1,   -1,
         1,   -1,    1,    1,    1,   -1,   -1,   -1,
        -1,   -1,   -1,   -1,   -1,   -1 }; 

