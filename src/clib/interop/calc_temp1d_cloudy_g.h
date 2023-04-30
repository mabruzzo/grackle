#ifndef FORTRAN_INTERFACE
/***********************************************************************
/
/ Declare function to compute temperature of a 1D slice using a cloudy table
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this
/ software.
************************************************************************/

#ifndef __CALC_TEMP1D_CLOUDY_G_H
#define __CALC_TEMP1D_CLOUDY_G_H

#include <stdint.h> // int32_t

// in the future, we will define GR_RESTRICT to be the restrict keyword
// introduced in C99 (or compiler-specific alternatives). For now, haven't done
// this to avoid headaches with the C++ unit-tests (restrict is not recognized
// by C++, but compiler-specific alternatives are)
//
// I don't actually think it would be problematic to use restrict within source
// files and skip it in header files (as long as assumptions about no-pointer
// aliasing) are enforced 
#ifndef GR_RESTRICT
#define GR_RESTRICT /* ... */
#endif /* GR_RESTRICT */

/// Calculate temperature and mean molecular weight for tabulated cooling.
/// @details This performs a calculation for a 1D slice from a 3D field.
///
/// @author Britton Smith
/// @author Matthew Abruzzo
/// @date May, 2015
/// @remark This was originally written by Britton Smith in May 2015 in
///     Fortran. It was transcribed to C by Matthew Abruzzo in April 2023.
///
/// @param[in]  d  3D density field
/// @param[in]  metal  3D metal density field
/// @param[in]  e  3D specific internal energy field
/// @param[in]  rhoH  (precomputed) total H mass density. This only holds
///     values for the 1D slice
/// @param[in]  in,jn,kn  dimensions of 3D fields (1D array hold ``in`` items)
/// @param[in]  is,ie  start and (inclusive) end indices of active region
///     (zero-based)
/// @param[in]  j,k  indices along other dimensions (one-based)
/// @param[out] tgas  1D array to store output temperature values
/// @param[out] mmw  1D array to store output mean molecular weight values
/// @param[in]  dom  unit conversion to proper number density in code units
/// @param[in]  zr  current redshift
/// @param[in]  temstart, temend  start and end of temperature range for rate
///     table
/// @param[in]  gamma  adiabatic index
/// @param[in]  utem  temperature units
/// @param[in]  imetal  flag if metal field is active (0 = no, 1 = yes)
/// @param[in]  clGridRank  rank of cloudy cooling data grid
/// @param[in]  clGridDim  array containing dimensions of cloudy data
/// @param[in]  clPar1, clPar2, clPar3  arrays containing cloudy grid parameter
///     values.
/// @param[in]  clDataSize  total size of flattened mmw data array
/// @param[in]  clMMW  cloudy mmw data
/// @param[in]  itmask  iteration mask
///
/// @note
/// None of the pointers are allowed to alias. Other than ``clPar2`` (when
/// ``clGridRank < 2``) or ``clPar3`` (when ``clGridRank < 3``), no parameters
/// should be passed ``NULL``
void calc_temp1d_cloudy_g(
        const gr_float* d, const gr_float* metal, const gr_float* e,
        const double* rhoH,
        int in, int jn, int kn, int is, int ie, int j, int k,
        GR_RESTRICT double* tgas, GR_RESTRICT double* mmw,
        double dom, double zr, double temstart, double temend, double gamma,
        double utem, int imetal,
        long long clGridRank,
        const long long* clGridDim,
        const double* clPar1,
        const double* clPar2,
        const double* clPar3,
        long long clDataSize,
        const double* clMMW,
        const int32_t* itmask);

#endif /*__CALC_TEMP1D_CLOUDY_G_H */

#else

! we explicitly avoid the header guards inclusion guards here
! - it is not clear how nicely those would play with the including
!   this file in multiple fortran functions/subroutines in a single
!   file
!
! before including this header, be sure that:
!   1. you have invoked ``USE ISO_C_BINDING``
!   2. grackle_fortran_types.def has already been included

      interface
        subroutine calc_temp1d_cloudy_g(d, metal, e, rhoH,
     &                in, jn, kn, is, ie, j, k,
     &                tgas, mmw, dom, zr, 
     &                temstart, temend,
     &                gamma, utem, imetal,
     &                clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3,
     &                clDataSize, clMMW,
     &                itmask) bind(C)
          IMPORT
          R_PREC, INTENT(IN) :: d(*)
          R_PREC, INTENT(IN) :: metal(*)
          R_PREC, INTENT(IN) :: e(*)
          REAL(C_DOUBLE), INTENT(IN) :: rhoH(*)
          INTEGER(C_INT), VALUE, INTENT(IN) :: in
          INTEGER(C_INT), VALUE, INTENT(IN) :: jn
          INTEGER(C_INT), VALUE, INTENT(IN) :: kn
          INTEGER(C_INT), VALUE, INTENT(IN) :: is
          INTEGER(C_INT), VALUE, INTENT(IN) :: ie
          INTEGER(C_INT), VALUE, INTENT(IN) :: j
          INTEGER(C_INT), VALUE, INTENT(IN) :: k
          REAL(C_DOUBLE), INTENT(OUT) :: tgas(*)
          REAL(C_DOUBLE), INTENT(OUT) :: mmw(*)
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: dom
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: zr
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: temstart
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: temend
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: gamma
          REAL(C_DOUBLE), VALUE, INTENT(IN) :: utem
          INTEGER(C_INT), VALUE, INTENT(IN) :: imetal
          INTEGER(C_LONG_LONG), VALUE, INTENT(IN) :: clGridRank
          INTEGER(C_LONG_LONG), INTENT(IN) :: clGridDim(*)
          REAL(C_DOUBLE), INTENT(IN) :: clPar1(*)
          REAL(C_DOUBLE), INTENT(IN) :: clPar2(*)
          REAL(C_DOUBLE), INTENT(IN) :: clPar3(*)
          INTEGER(C_LONG_LONG), VALUE, INTENT(IN) :: clDataSize
          REAL(C_DOUBLE), INTENT(IN) :: clMMW(*)
          INTEGER(C_INT32_T), INTENT(IN) :: itmask(*)
        end subroutine calc_temp1d_cloudy_g
      end interface

#endif
