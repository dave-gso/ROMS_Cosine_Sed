/*
** svn $Id: bio_umaine_def.h 645 2013-01-22 23:21:54Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2013 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Defines UMAINE CoSiNE ecosystem model input parameters in         **
**  output NetCDF files. It is included in routine "def_info.F".      **
**                                                                    **
!  By: PENG XIU 12/2013                                                !
************************************************************************
*/

!
!  Define UMaine CoSiNE ecosystem parameters.
!
      Vinfo( 1)='BioIter'
      Vinfo( 2)='number of iterations to achieve convergence'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='reg1'
      Vinfo( 2)='microzooplankton excretion rate to ammonium'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='reg2'
      Vinfo( 2)='mesozooplankton excretion rate to ammonium'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='gmaxs1'
      Vinfo( 2)='maximum specific growth rate of small phytoplankton'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='gmaxs2'
      Vinfo( 2)='maximum specific growth rate of diatom'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

#ifdef PHYTO_RESP
	Vinfo( 1)='rrb1'
	Vinfo( 2)='basal respiration rate of small phyto.'
	Vinfo( 3)='day-1'
	status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
	IF (exit_flag.ne.NoError) RETURN

	Vinfo( 1)='rrb2'
	Vinfo( 2)='basal respiration rate of diatom'
	Vinfo( 3)='day-1'
	status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
	IF (exit_flag.ne.NoError) RETURN

	Vinfo( 1)='rrg1'
	Vinfo( 2)='growth-dependent respiration rate of small phyto.'
	Vinfo( 3)='unitless'
	status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
	IF (exit_flag.ne.NoError) RETURN

	Vinfo( 1)='rrg2'
	Vinfo( 2)='growth-dependent respiration rate of diatom'
	Vinfo( 3)='unitless'
	status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
	IF (exit_flag.ne.NoError) RETURN
#endif
      Vinfo( 1)='beta1'
      Vinfo( 2)='microzooplankton maximum grazing rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='beta2'
      Vinfo( 2)='mesozooplankton maximum grazing rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='akz1'
      Vinfo( 2)='half saturation for microzooplankton grazing'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='akz2'
      Vinfo( 2)='half saturation for mesozooplankton grazing'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='PARfrac'
      Vinfo( 2)='photosynthetically available radiation fraction'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='amaxs1'
      Vinfo( 2)='initial slope of P-I curve of small phytoplankton'
      Vinfo( 3)='meter2 watts-1 day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='amaxs2'
      Vinfo( 2)='initial slope of P-I curve of diatom'
      Vinfo( 3)='meter2 watts-1 day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN


      Vinfo( 1)='parsats1'
      Vinfo( 2)='PAR saturation onset parameter for iSphy'
      Vinfo( 3)='watts m-2'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='parsats2'
      Vinfo( 2)='PAR saturation onset parameter for iLphy'
      Vinfo( 3)='watts m-2'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='pis1'
      Vinfo( 2)='ammonium inhibition parameter for iSphy'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='pis2'
      Vinfo( 2)='ammonium inhibition parameter for iLphy'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='akno3s1'
      Vinfo( 2)='half saturation of nitrate uptake by iSphy'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='akno3s2'
      Vinfo( 2)='half saturation of nitrate uptake by iLphy'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='aknh4s1'
      Vinfo( 2)='half saturation of ammonium uptake by iSphy'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='aknh4s2'
      Vinfo( 2)='half saturation of ammonium uptake by iLphy'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='akpo4s1'
      Vinfo( 2)='half saturation of phosphate uptake by iSphy'
      Vinfo( 3)='millimole_P meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='akpo4s2'
      Vinfo( 2)='half saturation of phosphate uptake by iLphy'
      Vinfo( 3)='millimole_P meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='akco2s1'
      Vinfo( 2)='half saturation of co2 uptake by iSphy'
      Vinfo( 3)='millimole_C meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='akco2s2'
      Vinfo( 2)='half saturation of co2 uptake by iLphy'
      Vinfo( 3)='millimole_C meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='aksio4s2'
      Vinfo( 2)='half saturation of silicate uptake by iLphy'
      Vinfo( 3)='millimole_Si meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='akox'
      Vinfo( 2)='half saturation of oxidation'
      Vinfo( 3)='millimole_O meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ak1'
      Vinfo( 2)='light attenuation coefficient of water'
      Vinfo( 3)='meter-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ak2'
      Vinfo( 2)='phytoplankton specific light attenuation rate'
      Vinfo( 3)='meter2 millimole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='bgamma0'
      Vinfo( 2)='mesozooplankton specific mortality rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='bgamma1'
      Vinfo( 2)='grazing efficiency of microzooplankton'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='bgamma2'
      Vinfo( 2)='grazing efficiency of mesozooplankton'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='bgamma3'
      Vinfo( 2)='death rate of small phytoplankton'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='bgamma4'
      Vinfo( 2)='death rate of large phytoplankton'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='bgamma5'
      Vinfo( 2)='decay rate of nitrogenous detritus'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='bgamma5s'
      Vinfo( 2)='decay rate of silicious detritus'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='bgamma6'
      Vinfo( 2)='aggregation rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='bgamma7'
      Vinfo( 2)='nitrification rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='wsd'
      Vinfo( 2)='sinking velocity of detritus'
      Vinfo( 3)='meter day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='wsdsi'
      Vinfo( 2)='sinking velocity of detritus silicate'
      Vinfo( 3)='meter day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='wsp'
      Vinfo( 2)='sinking velocity of large phytoplankton'
      Vinfo( 3)='meter day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='pco2a'
      Vinfo( 2)='air CO2 partial pressure'
      Vinfo( 3)='parts per million by volume'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='si2n'
      Vinfo( 2)='silicate to nitrogen ratio'
      Vinfo( 3)='mole_Si mole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='p2n'
      Vinfo( 2)='Phosphorus to nitrogen ratio'
      Vinfo( 3)='mole_P mole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='o2no'
      Vinfo( 2)='oxygen to nitrate ratio'
      Vinfo( 3)='mole_O2 mole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='o2nh'
      Vinfo( 2)='oxygen to ammonium ratio'
      Vinfo( 3)='mole_O2 mole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='c2n'
      Vinfo( 2)='carbon to nitrogen ratio'
      Vinfo( 3)='mole_C mole_N-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ro5'
      Vinfo( 2)='grazing preference for diatom'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ro6'
      Vinfo( 2)='grazing preference for mesozooplankton'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ro7'
      Vinfo( 2)='grazing preference for detritus'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return
      
      Vinfo( 1)='Chl2cs1_m'
      Vinfo( 2)='Maximum chlorophyll to carbon ratio for sphy'
      Vinfo( 3)='mg_Chl/mg_C'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return
      
      Vinfo( 1)='Chl2cs2_m'
      Vinfo( 2)='Maximum chlorophyll to carbon ratio for diatom'
      Vinfo( 3)='mg_Chl/mg_C'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

#ifdef HAB
      Vinfo( 1)='gmaxs3'
      Vinfo( 2)='maximum specific growth rate of HAB phyto'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN
# ifdef PHYTO_RESP
	Vinfo( 1)='rrb3'
	Vinfo( 2)='basal respiration rate of HAB phyto'
	Vinfo( 3)='day-1'
	status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
	IF (exit_flag.ne.NoError) RETURN

	Vinfo( 1)='rrg3'
	Vinfo( 2)='growth-dependent respiration rate of HAB phyto.'
	Vinfo( 3)='unitless'
	status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
	IF (exit_flag.ne.NoError) RETURN
# endif
      Vinfo( 1)='amaxs3'
      Vinfo( 2)='initial slope of P-I curve of HAB phytoplankton'
      Vinfo( 3)='meter2 watts-1 day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='parsats3'
      Vinfo( 2)='PAR saturation onset parameter for iHphy'
      Vinfo( 3)='watts m-2'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='pis3'
      Vinfo( 2)='ammonium inhibition parameter for iHphy'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='akno3s3'
      Vinfo( 2)='half saturation of nitrate uptake by iHphy'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='aknh4s3'
      Vinfo( 2)='half saturation of ammonium uptake by iHphy'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='akpo4s3'
      Vinfo( 2)='half saturation of phosphate uptake by iHphy'
      Vinfo( 3)='millimole_P meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='akco2s3'
      Vinfo( 2)='half saturation of co2 uptake by iHphy'
      Vinfo( 3)='millimole_C meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='aksio4s3'
      Vinfo( 2)='half saturation of silicate uptake by iHphy'
      Vinfo( 3)='millimole_Si meter-3'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='bgamma4s3'
      Vinfo( 2)='death rate of HAB phytoplankton'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='wsp3'
      Vinfo( 2)='sinking velocity of HAB phytoplankton'
      Vinfo( 3)='meter day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) RETURN

      Vinfo( 1)='ro5H'
      Vinfo( 2)='grazing preference for HAB phyto'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='Chl2cs3_m'
      Vinfo( 2)='Maximum chlorophyll to carbon ratio for HAB phyto'
      Vinfo( 3)='mg_Chl/mg_C'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return
#endif

#ifdef SEDBIO
      Vinfo( 1)='bUmax'
      Vinfo( 2)='Maximum POM decay rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/sedbiodim/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='bUmaxSi'
      Vinfo( 2)='Maximum particulate silicate decay rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/sedbiodim/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='bdep'
      Vinfo( 2)='Thickness of benthic sediment layer'
      Vinfo( 3)='meter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='balpha'
      Vinfo( 2)='Ratio of bacterial growth rate to maximum decay rate'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='bw'
      Vinfo( 2)='Sedimentation rate'
      Vinfo( 3)='centimeter/year'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='btheta_diag'
      Vinfo( 2)='Base for temperature dependence of diagenesis'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='bnit'
      Vinfo( 2)='Nitrification rate at reference temp. and unlimited O2'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='btheta_nit'
      Vinfo( 2)='Base for temperature dependence of nitrification'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='bdo_nit'
      Vinfo( 2)='half saturation O2 concentration for nitrification'
      Vinfo( 3)='mg/liter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='bdenit'
      Vinfo( 2)='Denitrification rate at reference temp. and anoxia'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='btheta_denit'
      Vinfo( 2)='Base for temperature dependence of denitrification'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='bpsi_n'
      Vinfo( 2)='Constant for oxygen dependence of denitrification'
      Vinfo( 3)='mg/liter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='bdo_c'
      Vinfo( 2)='half saturation constant for carbon oxidation'
      Vinfo( 3)='mg/liter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='bao2'
      Vinfo( 2)='stoichiometric coeff. for oxygen consumed by nitrif.'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='bpi'
      Vinfo( 2)='fraction of PO4 dissolved in pore water under anoxia'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='bpsi_p'
      Vinfo( 2)='Constant for oxygen dependence of dissolved PO4'
      Vinfo( 3)='mg/liter'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='bfc'
      Vinfo( 2)='Particulate carbon partitioning coefficient'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/sedbiodim/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='bfn'
      Vinfo( 2)='Particulate nitrogen partitioning coefficient'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/sedbiodim/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='bfp'
      Vinfo( 2)='Particulate phosphorus partitioning coefficient'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/sedbiodim/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

      Vinfo( 1)='bfs'
      Vinfo( 2)='Particulate silicon partitioning coefficient'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/sedbiodim/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return

#endif
#ifdef OPTICS_OP1
      Vinfo( 1)='optic_upd_fac'
      Vinfo( 2)='compute kdPAR every optic_upd_fac timesteps'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (exit_flag.ne.NoError) return
#endif

	Vinfo( 1)='q10_phyto_prod'
	Vinfo( 2)='Q10 for phytoplankton production'
	Vinfo( 3)='nondimensional'
	status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
	IF (exit_flag.ne.NoError) return

	Vinfo( 1)='q10_phyto_resp'
	Vinfo( 2)='Q10 for phytoplankton respiration'
	Vinfo( 3)='nondimensional'
	status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
	IF (exit_flag.ne.NoError) return

	Vinfo( 1)='q10_bact_remin'
	Vinfo( 2)='Q10 for bacterial remineralization'
	Vinfo( 3)='nondimensional'
	status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
	IF (exit_flag.ne.NoError) return

	Vinfo( 1)='q10_opal_dissol'
	Vinfo( 2)='Q10 for opal dissolution'
	Vinfo( 3)='nondimensional'
	status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
	IF (exit_flag.ne.NoError) return

	Vinfo( 1)='q10_zoopl_resp'
	Vinfo( 2)='Q10 for zooplankton respiration'
	Vinfo( 3)='nondimensional'
	status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
	IF (exit_flag.ne.NoError) return

	Vinfo( 1)='q10_zoopl_graz'
	Vinfo( 2)='Q10 for zooplankton grazing'
	Vinfo( 3)='nondimensional'
	status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
	IF (exit_flag.ne.NoError) return

	Vinfo( 1)='q10_nitrif'
	Vinfo( 2)='Q10 for nitrification'
	Vinfo( 3)='nondimensional'
	status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
	IF (exit_flag.ne.NoError) return
