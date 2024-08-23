!
! Copyright (c) 2024 Graziano Giuliani from UNESCO ICTP
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
module mod_rad_common

  use mod_constants
  use iso_fortran_env

  implicit none

  private

  real(rkx) , dimension(nspi) , parameter :: wavmin = &
    [ 0.200_rkx , 0.245_rkx , 0.265_rkx , 0.275_rkx , 0.285_rkx ,&
      0.295_rkx , 0.305_rkx , 0.350_rkx , 0.640_rkx , 0.700_rkx ,&
      0.701_rkx , 0.701_rkx , 0.701_rkx , 0.701_rkx , 0.702_rkx ,&
      0.702_rkx , 2.630_rkx , 4.160_rkx , 4.160_rkx ]

  real(rkx) , dimension(nspi) , parameter :: wavmax = &
    [ 0.245_rkx , 0.265_rkx , 0.275_rkx , 0.285_rkx , 0.295_rkx ,&
      0.305_rkx , 0.350_rkx , 0.640_rkx , 0.700_rkx , 5.000_rkx ,&
      5.000_rkx , 5.000_rkx , 5.000_rkx , 5.000_rkx , 5.000_rkx ,&
      5.000_rkx , 2.860_rkx , 4.550_rkx , 4.550_rkx ]

  real(rkx) , dimension(n_greenhouse_gases) , parameter :: cgunit = &
    [ 1.0e-6, 1.0e-9, 1.0e-9, 1.0e-12, 1.0e-12 ]

  type ghg_mf
    character(len=8) :: sname
    integer(ik4) :: year
    integer(ik4) :: month
    real(rkx) , dimension(:,:) , pointer :: gmf
  end type ghg_mf

  real(rkx) , dimension(:) , pointer :: heppatsi => null()
  real(rkx) , parameter :: tsifac = 0.9965_rkx

  character(len=maxpath) :: ghgpath
  type(ghg_mf) ::ghgc

  public :: wavmin , wavmax
  public :: initghg , ghgval
  public :: read_solarforcing , solar_irradiance

  contains

  subroutine initghg(path,scenario)
    implicit none
    character(len=*) , intent(in) :: path
    character(len=*) , intent(in) :: scenario
    ghgpath = path
    call load_scenario(scenario,1950,1,ghgc)
  end subroutine initghg

  real(rkx) function ghgval(igas,year,month,lat)
    implicit none
    integer(ik4) , intent(in) :: igas , year , month
    real(rkx) , intent(in) :: lat
    integer(ik4) :: ilat
    if ( ghgc%year /= year .or. ghgc%month /= month ) then
      call load_scenario(ghgc%sname,year,month,ghgc)
    end if
    ilat = int((lat+89.75_rkx)/0.5_rkx) + 1
    ghgval = ghgc%gmf(ilat,igas) * cgunit(igas)
  end function ghgval

  subroutine load_scenario(sname,year,month,ghgmf)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: sname
    integer(ik4) , intent(in) :: year , month
    type(ghg_mf) , intent(inout) :: ghgmf
    integer(ik4) , parameter :: imax = 5
    integer(ik4) , parameter :: smax = 9

    character(len=1024) :: filename
    integer(ik4) :: imod , ires , itim , ityp
    integer(ik4) :: i , ierr
    integer(ik4) :: ncid , varid , dimid , irec
    integer(ik4) , dimension(2) :: istart , icount
    integer(ik4) , save :: nlat

    character(len=*) , parameter , dimension(1) :: resolution = &
     ['0p5x360deg_']
    character(len=*) , parameter , dimension(1) :: inptype = &
     ['_input4MIPs_GHGConcentrations_']
    character(len=*) , parameter , dimension(2) :: timeperiod = &
     ['000001-201412','201501-250012']
    character(len=*) , dimension(imax) , parameter :: var_name = &
     ['mole_fraction_of_carbon_dioxide_in_air', &
      'mole_fraction_of_nitrous_oxide_in_air ', &
      'mole_fraction_of_methane_in_air       ', &
      'mole_fraction_of_cfc11_in_air         ', &
      'mole_fraction_of_cfc12_in_air         ' ]
    character(len=*) , dimension(imax) , parameter :: varname = &
     ['mole-fraction-of-carbon-dioxide-in-air', &
      'mole-fraction-of-nitrous-oxide-in-air ', &
      'mole-fraction-of-methane-in-air       ', &
      'mole-fraction-of-cfc11-in-air         ', &
      'mole-fraction-of-cfc12-in-air         ' ]
    character(len=*) , dimension(smax) , parameter :: modelname = &
     ['CMIP_UoM-CMIP-1-2-0_gr-                            ', &
      'ScenarioMIP_UoM-IMAGE-ssp119-1-2-1_gr-             ', &
      'ScenarioMIP_UoM-IMAGE-ssp126-1-2-1_gr-             ', &
      'ScenarioMIP_UoM-MESSAGE-GLOBIOM-ssp245-1-2-1_gr-   ', &
      'ScenarioMIP_UoM-AIM-ssp370-1-2-1_gr-               ', &
      'ScenarioMIP_UoM-GCAM4-ssp434-1-2-1_gr-             ', &
      'ScenarioMIP_UoM-GCAM4-ssp460-1-2-1_gr-             ', &
      'ScenarioMIP_UoM-REMIND-MAGPIE-ssp534-over-1-2-1_gr-', &
      'ScenarioMIP_UoM-REMIND-MAGPIE-ssp585-1-2-1_gr-     ' ]

    ires = 1
    ityp = 1
    if ( year < 2015 ) then
      imod = 1
      itim = 1
      irec = max(year*12 + month,1)
    else
      itim = 2
      select case (sname)
        case ('SSP119', 'ssp119')
          imod = 2
        case ('SSP126', 'ssp126')
          imod = 3
        case ('SSP245', 'ssp245')
          imod = 4
        case ('SSP370', 'ssp370')
          imod = 5
        case ('SSP434', 'ssp434')
          imod = 6
        case ('SSP460', 'ssp460')
          imod = 7
        case ('SSP534', 'ssp534')
          imod = 8
        case ('SSP585', 'ssp585')
          imod = 9
        case default
          imod = 0
          write (error_unit,*) 'Unsupported emission scenario: ', sname
          write (error_unit,*) 'Use one in SSP supported values:'
          stop __FILE__//': UNSUPPORTED EMISSION SCENARIO'
      end select
      irec = min((year-2015)*12 + month,5832)
    end if
    ghgmf%sname = sname
    ghgmf%year = year
    ghgmf%month = month
    do i = 1 , imax
      filename = trim(varname(i)) // trim(inptype(ityp)) // &
        trim(modelname(imod)) // trim(resolution(ires)) // &
        trim(timeperiod(itim)) // '.nc'
      ierr = nf90_open(filename, nf90_nowrite, ncid)
      if ( ierr /= nf90_noerr ) then
        filename = trim(ghgpath) // '/' // 'GHG' // '/' // filename
        ierr = nf90_open(filename, nf90_nowrite, ncid)
        if ( ierr /= nf90_noerr ) then
          write (error_unit, *) nf90_strerror(ierr) , trim(filename)
          stop __FILE__//': CANNOT OPEN FILE '//trim(filename)
        end if
      end if
      if ( .not. associated(ghgmf%gmf) ) then
        ierr = nf90_inq_dimid(ncid,'lat',dimid)
        if ( ierr /= nf90_noerr ) then
          write (error_unit, *) nf90_strerror(ierr) , trim(filename)
          stop __FILE__//': DIMENSION LAT SEARCH ERROR'
        end if
        ierr = nf90_inquire_dimension(ncid,dimid,len=nlat)
        if ( ierr /= nf90_noerr ) then
          write (error_unit, *) nf90_strerror(ierr) , trim(filename)
          stop __FILE__//': DIMENSION LAT READ ERROR'
        end if
        allocate(ghgmf%gmf(nlat,imax))
      end if
      ierr = nf90_inq_varid(ncid,var_name(i),varid)
      if ( ierr /= nf90_noerr ) then
        write (error_unit, *) nf90_strerror(ierr) , trim(filename)
        stop __FILE__//': VARIABLE '//var_name(i)//' FIND ERROR'
      end if
      istart(1) = 1
      istart(2) = irec
      icount(1) = nlat
      icount(2) = 1
      ierr = nf90_get_var(ncid,varid,ghgmf%gmf(:,i),istart,icount)
      if ( ierr /= nf90_noerr ) then
        write (error_unit, *) nf90_strerror(ierr) , trim(filename)
        stop __FILE__//': VARIABLE '//var_name(i)//' READ ERROR'
      end if
      ierr = nf90_close(ncid)
      if ( ierr /= nf90_noerr ) then
        write (error_unit, *) nf90_strerror(ierr) , trim(filename)
        stop __FILE__//': CANNOT CLOSE FILE'
      end if
    end do
    write(output_unit,'(a,i0.4,i0.2,a,a)') &
          ' Load GHG for ',year,month,' ',trim(sname)
  end subroutine load_scenario

  subroutine read_solarforcing(path)
    use netcdf
    implicit none
    character(len=*) :: path
    character(len=maxpath) :: heppafile
    integer(ik4) :: iret , ncid , idimid , ntime , ivar

    write(output_unit,*) 'Read solar forcing total irradiance data...'
    heppafile = path//'/'//'SOLAR'//'/'// &
      'solarforcing-ref-mon_input4MIPs_solar_CMIP_SOLARIS-HEPPA'// &
      '-3-2_gn_185001-229912.nc'
    iret = nf90_open(heppafile,nf90_nowrite,ncid)
    if ( iret /= nf90_noerr ) then
      write (error_unit, *) trim(heppafile) , ': ', nf90_strerror(iret)
      stop __FILE__//': CANNOT OPEN SOLARFORCING FILE'
    end if
    iret = nf90_inq_dimid(ncid,'time',idimid)
    if ( iret /= nf90_noerr ) then
      write (error_unit, *) trim(heppafile) , ': ', nf90_strerror(iret)
      stop __FILE__//': CANNOT FIND TIME DIMENSION'
    end if
    iret = nf90_inquire_dimension(ncid,idimid,len=ntime)
    if ( iret /= nf90_noerr ) then
      write (error_unit, *) trim(heppafile) , ': ', nf90_strerror(iret)
      stop __FILE__//': CANNOT READ TIME DIMENSION'
    end if
    iret = nf90_inq_varid(ncid,'tsi',ivar)
    if ( iret /= nf90_noerr ) then
      write (error_unit, *) trim(heppafile) , ': ', nf90_strerror(iret)
      stop __FILE__//': VARIABLE TSI NOT FOUND'
    end if
    allocate(heppatsi(ntime))
    iret = nf90_get_var(ncid,ivar,heppatsi)
    if ( iret /= nf90_noerr ) then
      write (error_unit, *) trim(heppafile) , ': ', nf90_strerror(iret)
      stop __FILE__//': CANNOT READ VARIABLE TSI'
    end if
    iret = nf90_close(ncid)
    if ( iret /= nf90_noerr ) then
      write (error_unit, *) trim(heppafile) , ': ', nf90_strerror(iret)
      stop __FILE__//': CANNOT CLOSE SOLARFORCING FILE'
    end if
  end subroutine read_solarforcing

  real(rkx) function solar_irradiance(year,month)
    implicit none
    integer(ik4) , intent(in) :: year , month
    integer(ik4) :: idate
    idate = (year-1850)*12 + month
    if ( idate < 1 ) then
      idate = mod(idate,132)+1
    else if ( idate > 5400 ) then
      idate = 5400 - 132 + mod(idate,132)
    end if
    solar_irradiance = tsifac*heppatsi(idate)
  end function solar_irradiance

end module mod_rad_common
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
