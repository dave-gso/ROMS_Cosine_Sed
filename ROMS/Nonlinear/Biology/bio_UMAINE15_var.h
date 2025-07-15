/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2012 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Assigns metadata indices to the CoSiNE ecosystem model            **
**  variables that are used in input and output NetCDF files.         **
**  The metadata nformation is read from "varinfo.dat".               **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
!  By: PENG XIU 12/2013                                                !
************************************************************************
*/

!
!  Model state biological tracers.
!

              CASE ('idTvar(iNO3_)')
                idTvar(iNO3_)=varid
              CASE ('idTvar(iNH4_)')
                idTvar(iNH4_)=varid
              CASE ('idTvar(iSiOH)')
                idTvar(iSiOH)=varid
              CASE ('idTvar(iSphy)')
                idTvar(iSphy)=varid
              CASE ('idTvar(iLphy)')
                idTvar(iLphy)=varid
              CASE ('idTvar(iSzoo)')
                idTvar(iSzoo)=varid
              CASE ('idTvar(iLzoo)')
                idTvar(iLzoo)=varid
              CASE ('idTvar(iSDet)')
                idTvar(iSDet)=varid
              CASE ('idTvar(iopal)')
                idTvar(iopal)=varid
              CASE ('idTvar(iPO4_)')
                idTvar(iPO4_)=varid
#ifdef OXYGEN
              CASE ('idTvar(iOxyg)')
                idTvar(iOxyg)=varid
#endif
#ifdef CARBON
              CASE ('idTvar(iTIC_)')
                idTvar(iTIC_)=varid
              CASE ('idTvar(iTAlk)')
                idTvar(iTAlk)=varid
! AKB 6/22/23
! Added in CO2 downward surface flux when CARBON is on
! Unsure if order matters here so maybe this will throw an error

              CASE ('idfgCO2')
                idfgCO2=varid
              CASE ('idpCO2s')
                idpCO2s=varid
#endif
              CASE ('idTvar(iChl1)')
                idTvar(iChl1)=varid
              CASE ('idTvar(iChl2)')
                idTvar(iChl2)=varid
#ifdef OPTICS_OP1
              CASE ('idkdPAR')
                idkdPAR=varid
#endif
#ifdef HAB
              CASE ('idTvar(iHphy)')
                idTvar(iHphy)=varid
              CASE ('idTvar(iChl3)')
                idTvar(iChl3)=varid
#endif

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI                                                                                   
              CASE ('idTads(iNO3_)')                                                                  
                idTads(iNO3_)=varid                                                                   
              CASE ('idTads(iNH4_)')                                                                  
                idTads(iNH4_)=varid                                                                   
              CASE ('idTads(iSiOH)')                                                                  
                idTads(iSiOH)=varid                                                                   
              CASE ('idTads(iSphy)')                                                                  
                idTads(iSphy)=varid                                                                   
              CASE ('idTads(iLphy)')                                                                  
                idTads(iLphy)=varid                                                                   
              CASE ('idTads(iSzoo)')                                                                  
                idTads(iSzoo)=varid                                                                   
              CASE ('idTads(iLzoo)')                                                                  
                idTads(iLzoo)=varid
              CASE ('idTads(iSDet)')
                idTads(iSDet)=varid
              CASE ('idTads(iopal)')
                idTads(iopal)=varid
              CASE ('idTads(iPO4_)')
                idTads(iPO4_)=varid
#ifdef OXYGEN
               CASE ('idTads(iOxyg)')
                idTads(iOxyg)=varid   
#endif                                                                
# ifdef CARBON
                CASE ('idTads(iTIC_)')
                idTads(iTIC_)=varid
                CASE ('idTads(iTAlk)')
                idTads(iTAlk)=varid
# endif
                CASE ('idTads(iChl1)')
                idTads(iChl1)=varid
                CASE ('idTads(iChl2)')
                idTads(iChl2)=varid
#endif

! Do it again for climatologies
              CASE ('idTclm(iNO3_)')
                idTclm(iNO3_)=varid
              CASE ('idTclm(iNH4_)')
                idTclm(iNH4_)=varid
              CASE ('idTclm(iSiOH)')
                idTclm(iSiOH)=varid
              CASE ('idTclm(iSphy)')
                idTclm(iSphy)=varid
              CASE ('idTclm(iLphy)')
                idTclm(iLphy)=varid
              CASE ('idTclm(iSzoo)')
                idTclm(iSzoo)=varid
              CASE ('idTclm(iLzoo)')
                idTclm(iLzoo)=varid
              CASE ('idTclm(iSDet)')
                idTclm(iSDet)=varid
              CASE ('idTclm(iopal)')
                idTclm(iopal)=varid
              CASE ('idTclm(iPO4_)')
                idTclm(iPO4_)=varid
#ifdef OXYGEN
              CASE ('idTclm(iOxyg)')
                idTclm(iOxyg)=varid
#endif
#ifdef CARBON
              CASE ('idTclm(iTIC_)')
                idTclm(iTIC_)=varid
              CASE ('idTclm(iTAlk)')
                idTclm(iTAlk)=varid
#endif
              CASE ('idTclm(iChl1)')
                idTclm(iChl1)=varid
              CASE ('idTclm(iChl2)')
                idTclm(iChl2)=varid
                
!
!  Biological tracers open boundary conditions.
!

              CASE ('idTbry(iwest,iNO3_)')
                idTbry(iwest,iNO3_)=varid
              CASE ('idTbry(ieast,iNO3_)')
                idTbry(ieast,iNO3_)=varid
              CASE ('idTbry(isouth,iNO3_)')
                idTbry(isouth,iNO3_)=varid
              CASE ('idTbry(inorth,iNO3_)')
                idTbry(inorth,iNO3_)=varid
              CASE ('idTbry(iwest,iNH4_)')
                idTbry(iwest,iNH4_)=varid
              CASE ('idTbry(ieast,iNH4_)')
                idTbry(ieast,iNH4_)=varid
              CASE ('idTbry(isouth,iNH4_)')
                idTbry(isouth,iNH4_)=varid
              CASE ('idTbry(inorth,iNH4_)')
                idTbry(inorth,iNH4_)=varid
              CASE ('idTbry(iwest,iSiOH)')
                idTbry(iwest,iSiOH)=varid
              CASE ('idTbry(ieast,iSiOH)')
                idTbry(ieast,iSiOH)=varid
              CASE ('idTbry(isouth,iSiOH)')
                idTbry(isouth,iSiOH)=varid
              CASE ('idTbry(inorth,iSiOH)')
                idTbry(inorth,iSiOH)=varid
              CASE ('idTbry(iwest,iSphy)')
                idTbry(iwest,iSphy)=varid
              CASE ('idTbry(ieast,iSphy)')
                idTbry(ieast,iSphy)=varid
              CASE ('idTbry(isouth,iSphy)')
                idTbry(isouth,iSphy)=varid
              CASE ('idTbry(inorth,iSphy)')
                idTbry(inorth,iSphy)=varid
              CASE ('idTbry(iwest,iLphy)')
                idTbry(iwest,iLphy)=varid
              CASE ('idTbry(ieast,iLphy)')
                idTbry(ieast,iLphy)=varid
              CASE ('idTbry(isouth,iLphy)')
                idTbry(isouth,iLphy)=varid
              CASE ('idTbry(inorth,iLphy)')
                idTbry(inorth,iLphy)=varid
              CASE ('idTbry(iwest,iSzoo)')
                idTbry(iwest,iSzoo)=varid
              CASE ('idTbry(ieast,iSzoo)')
                idTbry(ieast,iSzoo)=varid
              CASE ('idTbry(isouth,iSzoo)')
                idTbry(isouth,iSzoo)=varid
              CASE ('idTbry(inorth,iSzoo)')
                idTbry(inorth,iSzoo)=varid
              CASE ('idTbry(iwest,iLzoo)')
                idTbry(iwest,iLzoo)=varid
              CASE ('idTbry(ieast,iLzoo)')
                idTbry(ieast,iLzoo)=varid
              CASE ('idTbry(isouth,iLzoo)')
                idTbry(isouth,iLzoo)=varid
              CASE ('idTbry(inorth,iLzoo)')
                idTbry(inorth,iLzoo)=varid
              CASE ('idTbry(iwest,iSDet)')
                idTbry(iwest,iSDet)=varid
              CASE ('idTbry(ieast,iSDet)')
                idTbry(ieast,iSDet)=varid
              CASE ('idTbry(isouth,iSDet)')
                idTbry(isouth,iSDet)=varid
              CASE ('idTbry(inorth,iSDet)')
                idTbry(inorth,iSDet)=varid
              CASE ('idTbry(iwest,iopal)')
                idTbry(iwest,iopal)=varid
              CASE ('idTbry(ieast,iopal)')
                idTbry(ieast,iopal)=varid
              CASE ('idTbry(isouth,iopal)')
                idTbry(isouth,iopal)=varid
              CASE ('idTbry(inorth,iopal)')
                idTbry(inorth,iopal)=varid
              CASE ('idTbry(iwest,iPO4_)')
                idTbry(iwest,iPO4_)=varid
              CASE ('idTbry(ieast,iPO4_)')
                idTbry(ieast,iPO4_)=varid
              CASE ('idTbry(isouth,iPO4_)')
                idTbry(isouth,iPO4_)=varid
              CASE ('idTbry(inorth,iPO4_)')
                idTbry(inorth,iPO4_)=varid
#ifdef OXYGEN
              CASE ('idTbry(iwest,iOxyg)')
                idTbry(iwest,iOxyg)=varid
              CASE ('idTbry(ieast,iOxyg)')
                idTbry(ieast,iOxyg)=varid
              CASE ('idTbry(isouth,iOxyg)')
                idTbry(isouth,iOxyg)=varid
              CASE ('idTbry(inorth,iOxyg)')
                idTbry(inorth,iOxyg)=varid
#endif
#ifdef CARBON
              CASE ('idTbry(iwest,iTIC_)')
                idTbry(iwest,iTIC_)=varid
              CASE ('idTbry(ieast,iTIC_)')
                idTbry(ieast,iTIC_)=varid
              CASE ('idTbry(isouth,iTIC_)')
                idTbry(isouth,iTIC_)=varid
              CASE ('idTbry(inorth,iTIC_)')
                idTbry(inorth,iTIC_)=varid
              CASE ('idTbry(iwest,iTAlk)')
                idTbry(iwest,iTAlk)=varid
              CASE ('idTbry(ieast,iTAlk)')
                idTbry(ieast,iTAlk)=varid
              CASE ('idTbry(isouth,iTAlk)')
                idTbry(isouth,iTAlk)=varid
              CASE ('idTbry(inorth,iTAlk)')
                idTbry(inorth,iTAlk)=varid
#endif
              CASE ('idTbry(iwest,iChl1)')
                idTbry(iwest,iChl1)=varid
              CASE ('idTbry(ieast,iChl1)')
                idTbry(ieast,iChl1)=varid
              CASE ('idTbry(isouth,iChl1)')
                idTbry(isouth,iChl1)=varid
              CASE ('idTbry(inorth,iChl1)')
                idTbry(inorth,iChl1)=varid
              CASE ('idTbry(iwest,iChl2)')
                idTbry(iwest,iChl2)=varid
              CASE ('idTbry(ieast,iChl2)')
                idTbry(ieast,iChl2)=varid
              CASE ('idTbry(isouth,iChl2)')
                idTbry(isouth,iChl2)=varid
              CASE ('idTbry(inorth,iChl2)')
                idTbry(inorth,iChl2)=varid
#ifdef HAB
              CASE ('idTbry(iwest,iChl3)')
                idTbry(iwest,iChl3)=varid
              CASE ('idTbry(ieast,iChl3)')
                idTbry(ieast,iChl3)=varid
              CASE ('idTbry(isouth,iChl3)')
                idTbry(isouth,iChl3)=varid
              CASE ('idTbry(inorth,iChl3)')
                idTbry(inorth,iChl3)=varid
              CASE ('idTbry(iwest,iHphy)')
                idTbry(iwest,iHphy)=varid
              CASE ('idTbry(ieast,iHphy)')
                idTbry(ieast,iHphy)=varid
              CASE ('idTbry(isouth,iHphy)')
                idTbry(isouth,iHphy)=varid
              CASE ('idTbry(inorth,iHphy)')
                idTbry(inorth,iHphy)=varid
#endif
!
!  Biological tracers point Source/Sinks (river runoff).
!

              CASE ('idRtrc(iNO3_)')
                idRtrc(iNO3_)=varid
              CASE ('idRtrc(iNH4_)')
                idRtrc(iNH4_)=varid
              CASE ('idRtrc(iSiOH)')
                idRtrc(iSiOH)=varid
              CASE ('idRtrc(iSphy)')
                idRtrc(iSphy)=varid
              CASE ('idRtrc(iLphy)')
                idRtrc(iLphy)=varid
              CASE ('idRtrc(iSzoo)')
                idRtrc(iSzoo)=varid
              CASE ('idRtrc(iLzoo)')
                idRtrc(iLzoo)=varid
              CASE ('idRtrc(iSDet)')
                idRtrc(iSDet)=varid
              CASE ('idRtrc(iopal)')
                idRtrc(iopal)=varid
              CASE ('idRtrc(iPO4_)')
                idRtrc(iPO4_)=varid
# ifdef OXYGEN
              CASE ('idRtrc(iOxyg)')
                idRtrc(iOxyg)=varid
# endif
# ifdef CARBON
              CASE ('idRtrc(iTIC_)')
                idRtrc(iTIC_)=varid
              CASE ('idRtrc(iTAlk)')
                idRtrc(iTAlk)=varid
# endif
              CASE ('idRtrc(iChl1)')
                idRtrc(iChl1)=varid
              CASE ('idRtrc(iChl2)')
                idRtrc(iChl2)=varid
# ifdef HAB
              CASE ('idRtrc(iHphy)')
                idRtrc(iHphy)=varid
              CASE ('idRtrc(iChl3)')
                idRtrc(iChl3)=varid
# endif
                
#ifdef DIAGNOSTICS_BIO
# ifdef SEDBIO
              CASE ('iDbio2(ibNO3fx)')
                iDbio2(ibNO3fx)=varid
              CASE ('iDbio2(ibNH4fx)')
                iDbio2(ibNH4fx)=varid
              CASE ('iDbio2(ibPO4fx)')
                iDbio2(ibPO4fx)=varid
              CASE ('iDbio2(ibSiO2fx)')
                iDbio2(ibSiO2fx)=varid
              CASE ('iDbio2(ibPONfx)')
                iDbio2(ibPONfx)=varid
              CASE ('iDbio2(ibPSifx)')
                iDbio2(ibPSifx)=varid
#  ifdef CARBON
              CASE ('iDbio2(ibTICfx)')
                iDbio2(ibTICfx)=varid
              CASE ('iDbio2(ibAlkfx)')
                iDbio2(ibAlkfx)=varid
#  endif
#  ifdef OXYGEN
              CASE ('iDbio2(ibO2fx)')
                iDbio2(ibO2fx)=varid
#  endif
# endif

# ifdef CARBON
              CASE ('iDbio2(iCOfx)')
                iDbio2(iCOfx)=varid
              CASE ('iDbio2(ipCO2)')
                iDbio2(ipCO2)=varid
# endif
# ifdef OXYGEN
              CASE ('iDbio2(iO2fx)')
                iDbio2(iO2fx)=varid
# endif       
              CASE ('iDbio3(iPPro1)')
                iDbio3(iPPro1)=varid
              CASE ('iDbio3(iPPro2)')
                iDbio3(iPPro2)=varid
              CASE ('iDbio3(iNO3u)')
                iDbio3(iNO3u)=varid
# ifdef OXYGEN
              CASE ('iDbio3(iO2pr)')
                iDbio3(iO2pr)=varid
              CASE ('iDbio3(initri)')
                iDbio3(initri)=varid
              CASE ('iDbio3(iremin)')
                iDbio3(iremin)=varid
              CASE ('iDbio3(izoopl)')
                iDbio3(izoopl)=varid
#  ifdef PHYTO_RESP
		  CASE ('iDbio3(iphyres)')
		    iDbio3(iphyres)=varid
#  endif
# endif
#endif

#ifdef SEDBIO
        ! Particulate matter variables
            CASE ('idbPM(ibPOC)')
                idbPM(ibPOC)=varid
            CASE ('idbPM(ibPON)')
                idbPM(ibPON)=varid
            CASE ('idbPM(ibPOP)')
                idbPM(ibPOP)=varid
            CASE ('idbPM(ibPSi)')
                idbPM(ibPSi)=varid
        ! Pore water concentration variables
            CASE ('idbPW(ibNO3)')
                idbPW(ibNO3)=varid
            CASE ('idbPW(ibNH4)')
                idbPW(ibNH4)=varid
            CASE ('idbPW(ibPO4)')
                idbPW(ibPO4)=varid
            CASE ('idbPW(ibSi)')
                idbPW(ibSi)=varid
            CASE ('idbPW(ibTIC)')
                idbPW(ibTIC)=varid
            CASE ('idbPW(ibAlk)')
                idbPW(ibAlk)=varid
        ! Benthic particulate matter decay rates
            CASE ('idbDR(ibUC)')
                idbDR(ibUC)=varid
            CASE ('idbDR(ibUN)')
                idbDR(ibUN)=varid
            CASE ('idbDR(ibUP)')
                idbDR(ibUP)=varid
            CASE ('idbDR(ibUS)')
                idbDR(ibUS)=varid
        ! Benthic fluxes
            CASE ('idbSF(ibJNO3)')
                idbSF(ibJNO3)=varid
            CASE ('idbSF(ibJNH4)')
                idbSF(ibJNH4)=varid
            CASE ('idbSF(ibJPO4)')
                idbSF(ibJPO4)=varid
            CASE ('idbSF(ibJSi)')
                idbSF(ibJSi)=varid
            CASE ('idbSF(ibSOD)')
                idbSF(ibSOD)=varid
            CASE ('idbSF(ibJTIC)')
                idbSF(ibJTIC)=varid
            CASE ('idbSF(ibJAlk)')
                idbSF(ibJAlk)=varid
#endif
