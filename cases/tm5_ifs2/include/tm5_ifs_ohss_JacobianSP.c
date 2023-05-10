/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Sparse Jacobian Data Structures File                             */
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
/* File                 : tm5_ifs_ohss_JacobianSP.c                 */
/* Time                 : Sat Apr  3 22:32:57 2021                  */
/* Working directory    : /home/WUR/krol005/kpp/examples            */
/* Equation file        : tm5_ifs_ohss.kpp                          */
/* Output root filename : tm5_ifs_ohss                              */
/*                                                                  */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "tm5_ifs_ohss_Parameters.h"
#include "tm5_ifs_ohss_Global.h"
#include "tm5_ifs_ohss_Sparse.h"



/* Sparse Jacobian Data                                             */

 /* Row indexes of the LU Jacobian of variables */

  int  LU_IROW[] = {
       0,  0,  0,  1,  1,  2,  2,  2,  2,  3,  3,  3,
       4,  4,  4,  5,  5,  5,  6,  6,  6,  6,  6,  7,
       7,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,  8,
       8,  9,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10,
      10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11 }; 

 /* Column indexes of the LU Jacobian of variables */

  int  LU_ICOL[] = {
       0,  3,  7,  1,  2,  2,  4,  8, 10,  3,  6,  7,
       4, 10, 11,  5,  6,  9,  3,  5,  6,  7,  9,  3,
       5,  6,  7,  8,  9, 10, 11,  6,  7,  8,  9, 10,
      11,  5,  6,  7,  8,  9, 10, 11,  4,  5,  6,  7,
       8,  9, 10, 11,  1,  2,  4,  8,  9, 10, 11 }; 

 /* Compressed row indexes of the LU Jacobian of variables */

  int  LU_CROW[] = {
       0,  3,  5,  9, 12, 15, 18, 23, 31, 37, 44, 52,
      59 }; 

 /* Diagonal indexes of the LU Jacobian of variables */

  int  LU_DIAG[] = {
       0,  3,  5,  9, 12, 15, 20, 26, 33, 41, 50, 58,
      59 }; 

