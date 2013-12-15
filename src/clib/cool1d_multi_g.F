#include "fortran.def"
#include "phys_const.def"

!=======================================================================
!/////////////////////  SUBROUTINE COOL1D_MULTI_G  \\\\\\\\\\\\\\\\\\\\\

      subroutine cool1d_multi_g(
     &                d, e, u, v, w, de, HI, HII, HeI, HeII, HeIII,
     &                in, jn, kn, nratec, 
     &                iexpand, ispecies, imetal, imcool, idust, idim, 
     &                is, ie, j, k, ih2co, ipiht, iter, igammah,
     &                aye, temstart, temend, z_solar,
     &                utem, uxyz, uaye, urho, utim,
     &                gamma,
     &                ceHIa, ceHeIa, ceHeIIa, ciHIa, ciHeIa, 
     &                ciHeISa, ciHeIIa, reHIIa, reHeII1a, 
     &                reHeII2a, reHeIIIa, brema, compa, gammaha,
     &                comp_xraya, comp_temp,
     &                piHI, piHeI, piHeII, comp1, comp2,
     &                HM, H2I, H2II, DI, DII, HDI, metal,
     &                hyd01ka, h2k01a, vibha, rotha, rotla,
     &                hyd01k, h2k01, vibh, roth, rotl,
     &                gpldla, gphdla, gpldl, gphdl,
     &                hdltea, hdlowa, hdlte, hdlow,
     &                gaHIa, gaH2a, gaHea, gaHpa, gaela,
     &                gasgra, 
     &                ceHI, ceHeI, ceHeII, ciHI, ciHeI, ciHeIS, ciHeII,
     &                reHII, reHeII1, reHeII2, reHeIII, brem,
     &                indixe, t1, t2, logtem, tdef, edot,
     &                tgas, tgasold, p2d, tdust, metallicity, rhoH, 
     &                ih2optical, iciecool, ciecoa, cieco,
     &                icmbTfloor, iClHeat,
     &                clEleFra, clGridRank, clGridDim,
     &                clPar1, clPar2, clPar3, clPar4, clPar5,
     &                clDataSize, clCooling, clHeating,
     &                itmask)

!  SOLVE RADIATIVE COOLING/HEATING EQUATIONS
!
!  written by: Yu Zhang, Peter Anninos and Tom Abel
!  date:       
!  modified1: January, 1996 by Greg Bryan; adapted to KRONOS
!  modified2: October, 1996 by GB; moved to AMR
!  modified3: February, 2003 by Robert Harkness; iteration mask
!  modified6: September, 2009 by BDS to include cloudy cooling
!
!  PURPOSE:
!    Solve the energy cooling equations.
!
!  INPUTS:
!    is,ie   - start and end indicies of active region (zero-based!)
!
!  PARAMETERS:
!

!-----------------------------------------------------------------------

      implicit NONE

!  Arguments

      integer in, jn, kn, is, ie, j, k, nratec, idim,
     &        iexpand, ih2co, ipiht, ispecies, imcool, idust, 
     &        imetal, igammah, ih2optical, iciecool
      real    aye, temstart, temend, z_solar,
     &        utem, uxyz, uaye, urho, utim,
     &        gamma
      real    d(in,jn,kn),    e(in,jn,kn),
     &        u(in,jn,kn),    v(in,jn,kn),     w(in,jn,kn),
     &       de(in,jn,kn),   HI(in,jn,kn),   HII(in,jn,kn),
     &      HeI(in,jn,kn), HeII(in,jn,kn), HeIII(in,jn,kn)
      real    HM(in,jn,kn),  H2I(in,jn,kn), H2II(in,jn,kn),
     &        DI(in,jn,kn),  DII(in,jn,kn), HDI(in,jn,kn),
     &        metal(in,jn,kn)
      real    hyd01ka(nratec), h2k01a(nratec), vibha(nratec), 
     &        rotha(nratec), rotla(nratec), gpldla(nratec),
     &        gphdla(nratec), hdltea(nratec), hdlowa(nratec)
      real    gaHIa(nratec), gaH2a(nratec), gaHea(nratec),
     &        gaHpa(nratec), gaela(nratec), gasgra(nratec), 
     &        ciecoa(nratec)
      real    ceHIa(nratec), ceHeIa(nratec), ceHeIIa(nratec),
     &        ciHIa(nratec), ciHeIa(nratec), ciHeISa(nratec), 
     &        ciHeIIa(nratec), reHIIa(nratec), reHeII1a(nratec), 
     &        reHeII2a(nratec), reHeIIIa(nratec), brema(nratec)
      real    compa, piHI, piHeI, piHeII, comp_xraya, comp_temp
      real    gammaha

!  Cloudy cooling data

      integer icmbTfloor, iClHeat, clGridRank, clDataSize
      integer clGridDim(clGridRank)
      real clEleFra
      real clPar1(clGridDim(1)), clPar2(clGridDim(2)),
     &     clPar3(clGridDim(3)), clPar4(clGridDim(4)),
     &     clPar5(clGridDim(5))
      real clCooling(clDataSize), clHeating(clDataSize)

!  Parameters

      double precision mh
      real    ZSOLAR, mu_metal
      parameter (mh = mass_h)      !DPC
      parameter (mu_metal = 16.d0)    ! approx. mean molecular weight of metals

!  Locals

      integer i, j1, iter
      real dom, qq, vibl, logtem0, logtem9, dlogtem, energy, zr
      real dt2, ttmin, comp1, comp2, scoef, acoef, HIdot,
     &     hdlte1, hdlow1, gamma2, x, fudge, fH2,
     &     gphdl1, factor, dom_inv, tau, ciefudge
      double precision coolunit, dbase1, tbase1, xbase1, rtunits,
     &     nH2, nother

!  Slice locals
 
      integer indixe(in)
      real t1(in), t2(in), logtem(in), tdef(in), p2d(in),
     &     tgas(in), tgasold(in), tdust(in), rhoH(in), nh(in),
     &     metallicity(in)
      double precision edot(in)

!  Cooling/heating slice locals

      double precision ceHI(in), ceHeI(in), ceHeII(in),
     &     ciHI(in), ciHeI(in), ciHeIS(in), ciHeII(in),
     &     reHII(in), reHeII1(in), reHeII2(in), reHeIII(in),
     &     brem(in), cieco(in)
      real hyd01k(in), h2k01(in), vibh(in), roth(in), rotl(in),
     &     gpldl(in), gphdl(in), hdlte(in), hdlow(in), metalc(in)
      real gaHI(in), gaH2(in), gaHe(in), gaHp(in), gael(in),
     &     galdl(in), gasgr(in), gasgr_tdust(in)

!  Metal cooling locals

      integer tcmb_idx
      real logxe0, logxe1, dlogxea, xe, log_xe, xe1, dlogT, dlogxe
      real logTcmb, dlogTcmb, metalc_cmb, metalfree

!  Iteration mask

      logical itmask(in)

!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\/////////////////////////////////
!=======================================================================

!     Set log values of start and end of lookup tables

      logtem0 = log(temstart)
      logtem9 = log(temend)
      dlogtem= (log(temend) - log(temstart))/real(nratec-1)

!     Set units

      dom      = urho*(aye**3)/mh
      dom_inv  = 1.d0/dom
      tbase1   = utim
      xbase1   = uxyz/(aye*uaye)    ! uxyz is [x]*a      = [x]*[a]*a'        '
      dbase1   = urho*(aye*uaye)**3 ! urho is [dens]/a^3 = [dens]/([a]*a')^3 '
      coolunit = (uaye**5 * xbase1**2 * mh**2) / (tbase1**3 * dbase1)
      zr       = 1.d0/(aye*uaye) - 1.d0
      fudge    = 1.d0

!     Set compton cooling coefficients (and temperature)

      comp1 = compa  * (1.d0 + zr)**4
      comp2 = 2.73d0 * (1.d0 + zr)

!     Compute Pressure

      do i = is+1, ie+1
         if ( itmask(i) ) then
            p2d(i) = (gamma - 1.d0)*d(i,j,k)*e(i,j,k)
         end if
      enddo

!     Compute temperature and H number density

      do i = is+1, ie+1
         if ( itmask(i) ) then
            tgas(i) = 
     &           (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.d0 +
     &           HI(i,j,k) + HII(i,j,k) + de(i,j,k)
            rhoH(i) = HI(i,j,k) + HII(i,j,k)
         end if
      enddo

!          (include molecular hydrogen, but ignore deuterium)

      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            if ( itmask(i) ) then
               tgas(i) = tgas(i) +
     &              HM(i,j,k) + (H2I(i,j,k) + H2II(i,j,k))/2.d0
               rhoH(i) = rhoH(i) + H2I(i,j,k) + H2II(i,j,k)
            end if
         enddo
      endif

!     Include metal species and calculate metallicity
      
      if (imetal .eq. 1) then
         do i = is+1, ie+1
            tgas(i) = tgas(i) + metal(i,j,k)/mu_metal
            metallicity(i) = metal(i,j,k) / d(i,j,k) / z_solar
         enddo
      endif

      do i = is+1, ie+1
         if ( itmask(i) ) then
            tgas(i) = max(p2d(i)*utem/tgas(i), temstart)
            nh(i) = rhoH(i) * dom
         end if
      enddo

!     Correct temperature for gamma from H2

      if (ispecies .gt. 1) then
         do i = is+1, ie+1
            if ( itmask(i) ) then
            nH2 = 0.5d0*(H2I(i,j,k) + H2II(i,j,k))
            nother = (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.d0 +
     &           HI(i,j,k) + HII(i,j,k) + de(i,j,k)
            if (nH2/nother .gt. 1.0d-3) then
               x = 6100.d0/tgas(i) ! not quite self-consistent
               if (x .gt. 10.d0) then
                  gamma2 = 0.5d0*5.d0
               else
                  gamma2 = 0.5d0*(5.d0 + 2.d0*x**2 * 
     &                 exp(x)/(exp(x)-1)**2)
               endif
            else
               gamma2 = 2.5d0
            endif
            gamma2 = 1.d0 + (nH2 + nother)/
     &                      (nH2*gamma2 + nother/(gamma-1.d0))
            tgas(i) = tgas(i) * (gamma2 - 1.d0)/(gamma - 1.d0)
            end if
         enddo
      endif

!     If this is the first time through, just set tgasold to tgas

      if (iter .eq. 1) then
         do i = is+1, ie+1
            if ( itmask(i) ) then
            tgasold(i) = tgas(i)
            end if
         enddo
      endif

!     --- 6 species cooling ---

      do i = is+1, ie+1
         if ( itmask(i) ) then

!        Compute log temperature and truncate if above/below table max/min

         logtem(i) = log(0.5d0*(tgas(i)+tgasold(i)))
         logtem(i) = max(logtem(i), logtem0)
         logtem(i) = min(logtem(i), logtem9)

!        Compute index into the table and precompute parts of linear interp

         indixe(i) = min(nratec-1,
     &                   max(1,int((logtem(i)-logtem0)/dlogtem)+1))
         t1(i) = (logtem0 + (indixe(i) - 1)*dlogtem)
         t2(i) = (logtem0 + (indixe(i)    )*dlogtem)
         tdef(i) = (logtem(i) - t1(i)) / (t2(i) - t1(i))

!        Lookup cooling values and do a linear temperature in log(T)

         ceHI(i) = ceHIa(indixe(i)) + tdef(i)
     &         *(ceHIa(indixe(i)+1) -ceHIa(indixe(i)))
         ceHeI(i) = ceHeIa(indixe(i)) + tdef(i)
     &         *(ceHeIa(indixe(i)+1) -ceHeIa(indixe(i)))
         ceHeII(i) = ceHeIIa(indixe(i)) + tdef(i)
     &         *(ceHeIIa(indixe(i)+1) -ceHeIIa(indixe(i)))
         ciHI(i) = ciHIa(indixe(i)) + tdef(i)
     &         *(ciHIa(indixe(i)+1) -ciHIa(indixe(i)))
         ciHeI(i) = ciHeIa(indixe(i)) + tdef(i)
     &         *(ciHeIa(indixe(i)+1) -ciHeIa(indixe(i)))
         ciHeIS(i) = ciHeISa(indixe(i)) + tdef(i)
     &         *(ciHeISa(indixe(i)+1) -ciHeISa(indixe(i)))
         ciHeII(i) = ciHeIIa(indixe(i)) + tdef(i)
     &         *(ciHeIIa(indixe(i)+1) -ciHeIIa(indixe(i)))
         reHII(i) = reHIIa(indixe(i)) + tdef(i)
     &         *(reHIIa(indixe(i)+1) -reHIIa(indixe(i)))
         reHeII1(i)=reHeII1a(indixe(i)) + tdef(i)
     &        *(reHeII1a(indixe(i)+1)-reHeII1a(indixe(i)))
         reHeII2(i)=reHeII2a(indixe(i)) + tdef(i)
     &        *(reHeII2a(indixe(i)+1)-reHeII2a(indixe(i)))
         reHeIII(i)=reHeIIIa(indixe(i)) + tdef(i)
     &        *(reHeIIIa(indixe(i)+1)-reHeIIIa(indixe(i)))
         brem(i) = brema(indixe(i)) + tdef(i)
     &         *(brema(indixe(i)+1) -brema(indixe(i)))

         end if
      enddo

!     Compute the cooling function

      do i = is+1, ie+1
         if ( itmask(i) ) then
         edot(i) = (

!                    Collisional excitations

     &             - ceHI  (i)*HI  (i,j,k)*de(i,j,k)              ! ce of HI
     &             - ceHeI (i)*HeII(i,j,k)*de(i,j,k)**2*dom/4.d0  ! ce of HeI
     &             - ceHeII(i)*HeII(i,j,k)*de(i,j,k)/4.d0         ! ce of HeII

!                    Collisional ionizations

     &             - ciHI  (i)*HI  (i,j,k)*de(i,j,k)              ! ci of HI
     &             - ciHeI (i)*HeI (i,j,k)*de(i,j,k)/4.d0         ! ci of HeI
     &             - ciHeII(i)*HeII(i,j,k)*de(i,j,k)/4.d0         ! ci of HeII
     &             - ciHeIS(i)*HeII(i,j,k)*de(i,j,k)**2*dom/4.d0  ! ci of HeIS

!                    Recombinations

     &             - reHII  (i)*HII  (i,j,k)*de(i,j,k)           ! re of HII
     &             - reHeII1(i)*HeII (i,j,k)*de(i,j,k)/4.d0      ! re of HeII
     &             - reHeII2(i)*HeII (i,j,k)*de(i,j,k)/4.d0      ! re of HeII
     &             - reHeIII(i)*HeIII(i,j,k)*de(i,j,k)/4.d0      ! re of HeIII

!                    Compton cooling or heating

     &             - comp1*(tgas(i)-comp2)*de(i,j,k)*dom_inv

!                    X-ray compton heating

     &             - comp_xraya*(tgas(i)-comp_temp)*de(i,j,k)*dom_inv

!                    Bremsstrahlung

     &             - brem(i)*(HII(i,j,k)+HeII(i,j,k)/4.d0+HeIII(i,j,k))
     &                        *de(i,j,k)

!                    Photoelectric heating by UV-irradiated dust

     &             + float(igammah)*gammaha*(HI(i,j,k)+HII(i,j,k))
     &             *dom_inv)

         if (edot(i) .ne. edot(i)) then
            write(6,*) 'NaN in edot[1]: ', i, j, k, edot(i), 
     &           HI(i,j,k), HII(i,j,k), HeI(i,j,k), HeII(i,j,k), 
     &           HeIII(i,j,k), de(i,j,k), d(i,j,k), 
     &           tgas(i), p2d(i)
         endif
         
         end if
      enddo
     
!     --- H2 cooling ---

      if (ispecies .gt. 1) then

#define USE_GLOVER_ABEL2008
#ifdef USE_GLOVER_ABEL2008
         do i = is+1, ie+1
            if ( itmask(i) ) then
            gaHI(i) = gaHIa(indixe(i)) + tdef(i)
     &         *(gaHIa(indixe(i)+1) - gaHIa(indixe(i)))
            gaH2(i) = gaH2a(indixe(i)) + tdef(i)
     &         *(gaH2a(indixe(i)+1) - gaH2a(indixe(i)))
            gaHe(i) = gaHea(indixe(i)) + tdef(i)
     &         *(gaHea(indixe(i)+1) - gaHea(indixe(i)))
            gaHp(i) = gaHpa(indixe(i)) + tdef(i)
     &         *(gaHpa(indixe(i)+1) - gaHpa(indixe(i)))
            gael(i) = gaela(indixe(i)) + tdef(i)
     &         *(gaela(indixe(i)+1) - gaela(indixe(i)))
            gphdl(i) = gphdla(indixe(i)) + tdef(i)
     &         *(gphdla(indixe(i)+1) - gphdla(indixe(i)))
            cieco(i) = ciecoa(indixe(i)) + tdef(i)
     &         *(ciecoa(indixe(i)+1) - ciecoa(indixe(i)))
            end if
         enddo

         do i = is+1, ie+1
            if ( itmask(i) ) then
#ifdef OPTICAL_DEPTH_FUDGE
            nH2 = 0.5d0*H2I(i,j,k)
            nother = (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.d0 +
     &                HI(i,j,k) + HII(i,j,k) + de(i,j,k)
            fH2 = nH2/(nH2 + nother)
            fudge = sqrt((40.d0 * 10.d0**(4.8d0 * 
     &               sqrt(max(log10(tgas(i)),2.d0)-2.d0)) / fH2**2)/
     &               ((nH2 + nother)*dom) )
            fudge = min(fudge, 1.d0)
#endif /* OPTICAL_DEPTH_FUDGE */
            ! Note that this optical depth approximation comes from
            ! RA04.
            if (ih2optical.eq.1) then
                fudge = (0.76d0*d(i,j,k)*dom/8.d9)**(-0.45d0)
                fudge = min(fudge, 1.d0)
            else
                fudge = 1.d0
            endif
            galdl(i) = gaHI(i) * HI(i,j,k)  + gaH2(i) * H2I(i,j,k)
     &               + gaHe(i) * HeI(i,j,k) + gaHp(i) * HII(i,j,k)
     &               + gael(i) * de(i,j,k)
            gphdl1 = gphdl(i)/dom
            edot(i) = edot(i) - float(ih2co)*fudge*H2I(i,j,k)*
     &           gphdl(i)/(1.d0 + gphdl1/galdl(i)) / (2.d0*dom)

            end if
         enddo
#else

#define USE_GALLI_PALLA1999
#define NO_OPTICAL_DEPTH_FUDGE

!        Use the Galli and Palla (1999) cooling rates for molecular H.

#ifdef USE_GALLI_PALLA1999

         do i = is+1, ie+1
            if ( itmask(i) ) then
            gpldl(i) = gpldla(indixe(i)) + tdef(i)
     &         *(gpldla(indixe(i)+1) - gpldla(indixe(i)))
            gphdl(i) = gphdla(indixe(i)) + tdef(i)
     &         *(gphdla(indixe(i)+1) - gphdla(indixe(i)))
            cieco(i) = ciecoa(indixe(i)) + tdef(i)
     &         *(ciecoa(indixe(i)+1) - ciecoa(indixe(i)))
            end if
         enddo

         do i = is+1, ie+1
            if ( itmask(i) ) then

#ifdef OPTICAL_DEPTH_FUDGE
            nH2 = 0.5d0*H2I(i,j,k)
            nother = (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.d0 +
     &                HI(i,j,k) + HII(i,j,k) + de(i,j,k)
            fH2 = nH2/(nH2 + nother)
            fudge = sqrt((40.d0 * 10.d0**(4.8d0 * 
     &               sqrt(max(log10(tgas(i)),2.d0)-2.d0)) / fH2**2)/
     &               ((nH2 + nother)*dom) )
            fudge = min(fudge, 1.d0)
#endif /* OPTICAL_DEPTH_FUDGE */
            ! Note that this optical depth approximation comes from
            ! RA04.
            if (ih2optical.eq.1) then
                fudge = (0.76d0*d(i,j,k)*dom/8.d9)**(-0.45d0)
                fudge = min(fudge, 1.d0)
            else
                fudge = 1.d0
            endif
            gphdl1 = gphdl(i)/(HI(i,j,k)*dom)
            edot(i) = edot(i) - float(ih2co)*fudge*H2I(i,j,k)*
     &           gphdl(i)/(1.d0 + gphdl1/gpldl(i)) / (2.d0*dom)

            end if
         enddo

#else /* USE_GALLI_PALLA1999 */

         do i = is+1, ie+1
            if ( itmask(i) ) then
            hyd01k(i) = hyd01ka(indixe(i)) + tdef(i)
     &         *(hyd01ka(indixe(i)+1)-hyd01ka(indixe(i)))
            h2k01(i) = h2k01a(indixe(i)) + tdef(i)
     &         *(h2k01a(indixe(i)+1) - h2k01a(indixe(i)))
            vibh(i) = vibha(indixe(i)) + tdef(i)
     &         *(vibha(indixe(i)+1) - vibha(indixe(i)))
            roth(i) = rotha(indixe(i)) + tdef(i)
     &         *(rotha(indixe(i)+1) - rotha(indixe(i)))
            rotl(i) = rotla(indixe(i)) + tdef(i)
     &         *(rotla(indixe(i)+1) - rotla(indixe(i)))
            cieco(i) = ciecoa(indixe(i)) + tdef(i)
     &         *(ciecoa(indixe(i)+1) - ciecoa(indixe(i)))
            end if
         enddo

         do i = is+1, ie+1
            if ( itmask(i) ) then
            qq   = 1.2d0*(HI(i,j,k)*dom)**0.77d0 + 
     &                (H2I(i,j,k)*dom/2.d0)**0.77d0
            vibl = (HI(i,j,k)*hyd01k(i) + 
     &             H2I(i,j,k)/2.d0*h2k01(i))
     &             *dom*8.18d-13

#ifdef OPTICAL_DEPTH_FUDGE
            nH2 = 0.5d0*H2I(i,j,k)
            nother = (HeI(i,j,k) + HeII(i,j,k) + HeIII(i,j,k))/4.d0 +
     &                HI(i,j,k) + HII(i,j,k) + de(i,j,k)
            fH2 = nH2/(nH2 + nother)
            fudge = sqrt((40.d0 * 10.d0**(4.8d0 * 
     &               sqrt(max(log10(tgas(i)),2.d0)-2.d0)) / fH2**2)/
     &               ((nH2 + nother)*dom) )
            fudge = min(fudge, 1.d0)
#endif /* OPTICAL_DEPTH_FUDGE */

            edot(i) = edot(i) - float(ih2co)*fudge*H2I(i,j,k)*(
     &          vibh(i)/(1.d0+vibh(i)/max(   vibl     ,tiny)) +
     &          roth(i)/(1.d0+roth(i)/max(qq*rotl(i),tiny))     
     &                                                      )/2.d0/dom
            end if
         enddo

#endif /* USE_GALLI_PALLA1999 */
#endif /* USE_GLOVER_ABEL2008 */

c     CIE
c     cooling from H2-H2 and He-H2 collisional induced emission comes
C     with its own radiative transfer correction as discussed in
C     Ripamonti & Abel 2003
         if (iciecool.eq.1) then
            do i = is+1, ie+1
            if (itmask(i)) then
c     Only calculate if H2I(i) is a substantial fraction
              if (d(i,j,k)*dom.gt.1e10) then
                ciefudge = 1.d0
                tau = ((d(i,j,k)/2d16)*dom)**2.8d0  ! 2e16 is in units of cm^-3
                tau = max(tau, 1.d-5)
                ciefudge = min((1.d0-exp(-tau))/tau,1.d0)
c               Matt's attempt at a second exponentialier cutoff
                tau = ((d(i,j,k)/2.d18)*dom)**8.d0  ! 2e18 is in units of cm^-3
                tau = max(tau, 1.d-5)
                ciefudge = ciefudge*min((1.-exp(-tau))/tau,1.d0)
c               ciefudge, which is applied to the continuum, is applied to edot
                edot(i) = ciefudge*(edot(i) - 
     &                  H2I(i,j,k)*(d(i,j,k)*cieco(i)))
              endif
            endif
            enddo
         endif

      endif

!     --- Cooling from HD ---

      if (ispecies .gt. 2) then
         do i = is+1, ie+1
            if ( itmask(i) ) then
c CMB cooling floor
               if (tgas(i) .gt. comp2) then
                  hdlte(i) = hdltea(indixe(i)) + tdef(i)
     &            *(hdltea(indixe(i)+1) - hdltea(indixe(i)))
                  hdlow(i) = hdlowa(indixe(i)) + tdef(i)
     &            *(hdlowa(indixe(i)+1) - hdlowa(indixe(i)))
               else
                  hdlte(i) = tiny
                  hdlow(i) = tiny
               endif
            end if
         enddo

         do i = is+1, ie+1
            if ( itmask(i) ) then
c  old (incorrect) way:
c               hdlte1 = hdlte(i)/(HDI(i,j,k)*dom/2.0)
c               hdlow1 = max(hdlow(i), tiny)
c               edot(i) = edot(i) - HDI(i,j,k)*
c     .                     (hdlte1/(1.0 + hdlte1/hdlow1)/(2.0*dom))
c  new (correct) way: (april 4, 2007)
               hdlte1 = hdlte(i)/(HI(i,j,k)*dom)
               hdlow1 = max(hdlow(i), tiny)
               edot(i) = edot(i) - HDI(i,j,k)*
     .                     (hdlte(i)/(1.d0 + hdlte1/hdlow1)) /(3.d0*dom)
            end if
         enddo
      endif

!     --- Gas to grain heat transfer ---

      if (idust .gt. 0) then

!     Look up gas/grain heat transfer rates

         do i = is+1, ie+1
            if ( itmask(i) ) then
               gasgr(i) = gasgra(indixe(i)) + tdef(i)
     &              *(gasgra(indixe(i)+1) -gasgra(indixe(i)))
               gasgr_tdust(i) = gasgr(i) * coolunit / mh
            endif
         enddo

!     Compute dust temperature

         call calc_tdust_1d_g(tdust, tgas, nh, gasgr_tdust, 
     &           itmask, comp2, in, is, ie, j, k)

!     Calculate cooling rate

         do i = is+1, ie+1
            if ( itmask(i) ) then
               edot(i) = edot(i) - 
     &              gasgr(i) * (tgas(i) - tdust(i)) * 
     &              metallicity(i) * rhoH(i) * rhoH(i)
            endif
         enddo

      endif

!     --- Compute (external) radiative heating terms ---
!     Photoionization heating

      do i = is+1, ie+1
         if ( itmask(i) ) then
            edot(i) = edot(i) + float(ipiht)*(
     &           + piHI  *HI  (i,j,k) ! pi of HI
     &           + piHeI *HeI (i,j,k)*0.25d0 ! pi of HeI
     &           + piHeII*HeII(i,j,k)*0.25d0 ! pi of HeII
     &           )/dom
         end if
      enddo

!     --- Cloudy metal cooling and heating ---

      if (imcool .eq. 1) then

         call cool1d_cloudy_g(d, de, rhoH, metallicity,
     &        in, jn, kn, is, ie, j, k,
     &        logtem, edot, comp2, ispecies, dom, zr,
     &        icmbTfloor, iClHeat, 
     &        clEleFra, clGridRank, clGridDim,
     &        clPar1, clPar2, clPar3, clPar4, clPar5,
     &        clDataSize, clCooling, clHeating,
     &        itmask)

      endif

!     Set tgasold

      do i=is+1, ie+1
         if ( itmask(i) ) then
         tgasold(i) = tgas(i)
         end if
      enddo

      return
      end