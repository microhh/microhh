/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/*                                                                  */
/* Sparse Hessian Data Structures File                              */
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
/* File                 : tm5_ifs_HessianSP.c                       */
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



/* Hessian Sparse Data                                              */
/*                                                                  */
 /* Index i of Hessian element d^2 f_i/dv_j.dv_k */

  int  IHESS_I[] = {
       0,  0,  1,  2,  2,  3,  3,  4,  5,  5,  5,  5,
       6,  6,  6,  6,  7,  7,  8,  8,  8,  9,  9,  9,
       9,  9,  9, 10, 10, 10, 10, 10, 10, 10, 10, 10,
      10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 13,
      13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14,
      14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15
      }; 

 /* Index j of Hessian element d^2 f_i/dv_j.dv_k */

  int  IHESS_J[] = {
      10, 14, 13,  2, 10,  3,  6,  5,  5,  8,  8,  8,
       6,  7,  9,  9,  7,  9,  8,  8,  8,  7,  8,  8,
       8,  9,  9,  2,  3,  6,  9,  9, 10, 10, 10, 10,
      10, 11,  8, 10, 11, 11, 11,  9, 10, 11, 12,  8,
      10, 11, 12, 13,  5,  8,  9, 10, 11, 11, 12, 13,
      14,  2,  3,  6,  7,  8,  8, 10, 10, 10, 11, 14
      }; 

 /* Index k of Hessian element d^2 f_i/dv_j.dv_k */

  int  IHESS_K[] = {
      13, 15, 14, 15, 10, 15, 15, 14, 14, 11, 13, 15,
      15, 15, 10, 12, 15, 10, 11, 13, 15, 15, 11, 13,
      15, 10, 12, 15, 15, 15, 10, 12, 10, 11, 12, 13,
      15, 15, 11, 11, 12, 14, 15, 12, 12, 12, 13, 13,
      13, 14, 13, 14, 14, 13, 12, 12, 12, 14, 13, 14,
      15, 15, 15, 15, 15, 11, 15, 11, 12, 15, 15, 15
      }; 

