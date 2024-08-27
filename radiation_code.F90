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
program radiation_code
  use mod_constants
  use mod_dimensions
  use mod_rad_common
  use mod_rad_radiation
  use mod_rad_aerosol
  use mod_sunorbit
  use iso_fortran_env

  implicit none

  integer :: year_start = 1950
  integer :: year_end   = 2024

  character(len=*) , parameter :: e5_base = 'data'
  character(len=*) , parameter :: cmip6_base = 'data'
  character(len=*) , parameter :: macv2_base = 'data/MACV2'
  character(len=*) , parameter :: e5_mask = e5_base//'/fixed/landmask.nc'
  character(len=*) , parameter :: e5_geop = e5_base//'/fixed/geopot.nc'
  character(len=maxpath) :: e5_ps , e5_t , e5_qv , e5_qc , e5_qi
  character(len=maxpath) :: e5_cld , e5_o3
  character(len=maxpath) :: e5_adirsw , e5_adifsw , e5_adirlw , e5_adiflw
  character(len=maxpath) :: e5_ts
  character(len=16) :: timepart
  type(radtype) :: rt
  integer(ik4) :: rgrid , reduced_points , nlev
  integer(ik4) :: n , k , it , ith , ntimes
  integer(ik4) :: n1 , n2
  integer(ik4) :: iyear , imonth , mday
  integer(ik8) :: istep
  real(rkx) , dimension(:) , pointer :: am , bm , ai , bi
  real(rkx) , dimension(:,:) , pointer :: dum
  real(rkx) :: lwc , iwc , tcels , tempc , fsr , desr , aiwc , biwc
  real(rkx) :: kabsi , kabsl , kabs , arg , cldemis
  real(rk8) :: eccen , obliq , mvelp , obliqr , lambm0 , mvelpp , calday
  real(rk8) :: declin , eccf
  integer(ik8) :: tstart , tstop , trate

  call system_clock(count_rate=trate)
  call initghg(cmip6_base,'SSP370')
  call read_solarforcing(cmip6_base)

  call read_dimensions(e5_geop,rgrid,nlev,reduced_points)

  write(output_unit,*) 'Number of reduced grid points: ', rgrid
  write(output_unit,*) 'Number of latitudes          : ', reduced_points
  write(output_unit,*) 'Number of vertical levels    : ', nlev

  call init_dimensions(rgrid,nlev)

  write(output_unit,*) 'Allocating vertical grids....'
  allocate(am(nlev),bm(nlev))
  allocate(ai(nlev+1),bi(nlev+1))

  write(output_unit,*) 'Allocating radiation type....'

  ! This must be for MPI par....
  n1 = 271000
  n2 = 272000

  call allocate_radtype(rt,n1,n2)
  call allocate_aerosol(n1,n2)
  allocate(dum(nlev,n1:n2))
  call read_geolocation(e5_mask,rt%dlat,rt%dlon)

  print *, 'Latitude extrema are  :', rt%dlat(n1), rt%dlat(n2)
  print *, 'Longitude extrema are :', rt%dlon(n1), rt%dlon(n2)

  call read_static_2di(e5_mask,'lsm',rt%ioro)
  call read_vertical_grid(e5_geop,am,bm,ai,bi)
  call read_static_2d(e5_geop,'z',rt%ht)

  do concurrent ( n = n1:n2 )
    rt%xptrop(n) = 25.0e3_rkx - 15.0e3_rkx*cos(rt%dlat(n)*degrad)**2
    rt%ht(n) = rt%ht(n)/egrav
  end do

  call init_aerosol(n1,n2,rt%ht,rt%dlat,rt%dlon,macv2_base)

  ! Call radiation code

  iyear = year_start
  istep = 0

  do iyear = year_start , year_end

    call orb_params(iyear,eccen,obliq,mvelp,obliqr,lambm0,mvelpp)

    calday = 0.0_rk8

    do imonth = 1 , 12

      rt%scon = solar_irradiance(iyear,imonth) * 1000.0_rkx ! cgs
      print *, 'Solar constant in CGS: ',rt%scon

      write(timepart,'(i0.4,a,i0.2)') iyear, '_', imonth
      e5_ps = e5_base//'/mlevel/lnps_'//trim(timepart)//'.nc'
      e5_t = e5_base//'/mlevel/t_'//trim(timepart)//'.nc'
      e5_qv = e5_base//'/mlevel/q_'//trim(timepart)//'.nc'
      e5_qc = e5_base//'/mlevel/clwc_'//trim(timepart)//'.nc'
      e5_qi = e5_base//'/mlevel/ciwc_'//trim(timepart)//'.nc'
      e5_cld = e5_base//'/mlevel/cc_'//trim(timepart)//'.nc'
      e5_o3 = e5_base//'/mlevel/o3_'//trim(timepart)//'.nc'
      e5_adirsw = e5_base//'/surface/aluvp_'//trim(timepart)//'.nc'
      e5_adifsw = e5_base//'/surface/aluvd_'//trim(timepart)//'.nc'
      e5_adirlw = e5_base//'/surface/alnip_'//trim(timepart)//'.nc'
      e5_adiflw = e5_base//'/surface/alnid_'//trim(timepart)//'.nc'
      e5_ts = e5_base//'/surface/skt_'//trim(timepart)//'.nc'

      ntimes = ntimes_in_file(e5_ps)

      do it = 1 , ntimes
        !
        !   Reset all arrays
        !
        rt%alb(:) = 0.0_rkx
        rt%albc(:) = 0.0_rkx
        rt%flns(:) = 0.0_rkx
        rt%flnsc(:) = 0.0_rkx
        rt%flnt(:) = 0.0_rkx
        rt%lwout(:) = 0.0_rkx
        rt%lwin(:) = 0.0_rkx
        rt%flntc(:) = 0.0_rkx
        rt%flwds(:) = 0.0_rkx
        rt%fsds(:) = 0.0_rkx
        rt%fsnirt(:) = 0.0_rkx
        rt%fsnirtsq(:) = 0.0_rkx
        rt%fsnrtc(:) = 0.0_rkx
        rt%fsns(:) = 0.0_rkx
        rt%fsnsc(:) = 0.0_rkx
        rt%fsnt(:) = 0.0_rkx
        rt%fsntc(:) = 0.0_rkx
        rt%solin(:) = 0.0_rkx
        rt%solout(:) = 0.0_rkx
        rt%soll(:) = 0.0_rkx
        rt%solld(:) = 0.0_rkx
        rt%sols(:) = 0.0_rkx
        rt%solsd(:) = 0.0_rkx
        rt%qrl(:,:) = 0.0_rkx
        rt%qrs(:,:) = 0.0_rkx

        ! No LW aerosol effect

        rt%aertrlw(:,:,:) = 1.0_rkx

        mday = (it-1)/4
        call orb_decl(calday,eccen,mvelpp,lambm0, &
                  obliqr,declin,eccf)
        rt%iyear = iyear
        rt%eccf = eccf
        rt%calday = calday

        ! RegCM has a defualt of 18 hours here...
        rt%labsem = (mod(istep,3) == 0)
        if ( rt%labsem ) print *, 'Computing Absorbtion emission...'

        print *, '########################################################'
        print *, iyear , imonth, mday+1, calday

        ! Time in hourly file consistent with 6h in mlev
        ith = (it-1) * 6 + 1

        call read_field_2d(e5_ps,'lnsp',rt%ps,1,it)
        call read_field_3d(e5_t,'t',rt%t,it)
        call read_field_3d(e5_qv,'q',rt%q,it)
        call read_field_3d(e5_qc,'clwc',rt%ql,it)
        call read_field_3d(e5_qi,'ciwc',rt%qi,it)
        call read_field_3d(e5_cld,'cc',dum,it)
        call read_field_3d(e5_o3,'o3',rt%o3vmr,it)

        do concurrent ( n = n1:n2 )
          rt%ps(n) = exp(rt%ps(n))
          do k = 1 , nlev
            rt%pmid(k,n) = am(k)+ bm(k)*rt%ps(n)
            rt%o3vmr(k,n) = rt%o3vmr(k,n) * (amd/amo3)
            rt%q(k,n) = rt%q(k,n) / (1.0_rkx-rt%q(k,n))
            rt%qi(k,n) = rt%qi(k,n) / (1.0_rkx-rt%qi(k,n))
            rt%ql(k,n) = rt%ql(k,n) / (1.0_rkx-rt%ql(k,n))
            if ( rt%qi(k,n) > 1.0e-11_rkx ) then
              if ( rt%ql(k,n) > 1.0e-11_rkx ) then
                rt%fice(k,n) = rt%qi(k,n) / (rt%ql(k,n)+rt%qi(k,n))
              else
                rt%fice(k,n) = 1.0_rkx
              end if
            else
              rt%fice(k,n) = 0.0_rkx
            end if
            rt%rho(k,n) = (rt%pmid(k,n))/(rt%t(k,n)*rgas)
            rt%rh(k,n) = min( 1.0_rkx,max( 0.0_rkx, &
              (rt%q(k,n)/pfwsat(rt%t(k,n),rt%pmid(k,n))) ) )
            rt%pmln(k,n) = log(rt%pmid(k,n))
          end do
          rt%pint(1,n) = 1.0e-10_rkx
          rt%piln(1,n) = log(rt%pint(1,n))
          do k = 2 , nlev+1
            rt%pint(k,n) = ai(k)+ bi(k)*rt%ps(n)
            rt%piln(k,n) = log(rt%pint(k,n))
          end do
          do k = 2 , nlev
            rt%cld(k,n) = dum(k-1,n) + dum(k,n) - (dum(k-1,n)*dum(k,n))
          end do
          rt%zq(nlev+1,n) = rt%ht(n)
          do k = nlev, 1 , -1
            rt%zq(k,n) = rt%zq(k+1,n) + rovg*rt%t(k,n) * &
                  log(rt%pint(k+1,n)/rt%pint(k,n))
          end do
          do k = 1 , nlev
            rt%dz(k,n) = rt%zq(k,n) - rt%zq(k+1,n)
            rt%za(k,n) = 0.5_rkx * (rt%zq(k,n)+rt%zq(k+1,n))
            ! If some cloud present
            if ( dum(k,n) > 0.001_rkx ) then
              ! In cloud here , clwp must be in g/m2 (mm)
              rt%clwp(k,n) = 1000.0_rkx * rt%rho(k,n) * &
                (rt%qi(k,n)+rt%ql(k,n)) / dum(k,n) * rt%dz(k,n)
            else
              ! Safe low
              rt%clwp(k,n) = 1.0e-10_rkx
            end if
          end do
        end do
        do n = n1 , n2
          do k = 1 , nlev
            ! If some clouds here
            if ( dum(k,n) > 0.001_rkx ) then
              lwc = (rt%ql(k,n) * rt%rho(k,n) * 1000.0_rkx)/dum(k,n)
              if ( rt%ioro(n) == 1 ) then
                rt%rel(k,n) = min(4.0_rkx + 7.0_rkx*lwc,15.0_rkx)
              else
                rt%rel(k,n) = min(5.5_rkx + 9.5_rkx*lwc,18.0_rkx)
              end if
              iwc = (rt%qi(k,n) * rt%rho(k,n) * 1000.0_rkx)/dum(k,n)
              tempc = rt%t(k,n) - 83.15_rkx
              tcels = tempc - tzero
              fsr = 1.2351_rkx + 0.0105_rkx * tcels
              aiwc = 45.8966_rkx * iwc**0.2214_rkx
              biwc = 0.7957_rkx * iwc**0.2535_rkx
              desr = fsr*(aiwc+biwc*tempc)
              desr = max(30.0_rkx,min(155.0_rkx,desr))
              rt%rei(k,n) = 0.64952_rkx*desr
              kabsi = 0.005_rkx + (1.0_rkx/rt%rei(k,n))
              kabsl = 0.31_rkx * exp(-0.08_rkx*rt%rel(k,n))
              kabs = kabsl*(1.0_rkx-rt%fice(k,n)) + (kabsi*rt%fice(k,n))
              arg = min(kabs*rt%clwp(k,n),25.0_rkx)
              cldemis = max(1.0_rkx - exp(-arg),0.00_rkx)
              rt%effcld(k,n) = rt%cld(k,n)*cldemis
            else
              rt%effcld(k,n) = 0.0_rkx
              rt%rel(k,n) = 4.0_rkx
              rt%rei(k,n) = 0.5_rkx
            end if
          end do
        end do
        call read_field_2d(e5_adirsw,'aluvp',rt%adirsw,0,ith)
        call read_field_2d(e5_adifsw,'aluvd',rt%adifsw,0,ith)
        call read_field_2d(e5_adirlw,'alnip',rt%adirlw,0,ith)
        call read_field_2d(e5_adiflw,'alnid',rt%adiflw,0,ith)
        call read_field_2d(e5_ts,'skt',rt%ts,0,ith)
        do concurrent ( n = n1:n2 )
          rt%emiss(n) = 1.0_rkx
          rt%asw(n) = max(rt%adirsw(n),rt%adifsw(n))
          rt%alw(n) = max(rt%adirlw(n),rt%adiflw(n))
          rt%czen(n) = max(min(orb_cosz(calday,rt%dlat(n)*degrad, &
                           rt%dlon(n)*degrad,declin),1.0_rkx),0.0_rkx)
        end do
        do concurrent ( n = n1:n2 )
          if ( rt%czen(n) > 0.0_rkx ) then
            rt%czengt0(n) = .true.
          else
            rt%czengt0(n) = .false.
          end if
        end do

        call system_clock(tstart)
        call radctl(rt,iyear,imonth)
        call system_clock(tstop)

        print *, '########################################################'
        print *, 'Elapsed time : ',real(tstop-tstart,rk8)/real(trate,rk8), &
          ' s'
        print *, '########################################################'
        print *, 'TOA SW incoming flux (W/m2) : ', rt%solin(n1)
        print *, 'TOA LW incoming flux (W/m2) : ', rt%lwin(n1)
        print *, 'TOA SW outgoing flux (W/m2) : ', rt%solout(n1)
        print *, 'TOA LW outgoing flux (W/m2) : ', rt%lwout(n1)
        print *, '########################################################'
        print *, 'TOA LW net radiation        : ', rt%flnt(n1)
        print *, 'TOA SW net radiation        : ', rt%fsnt(n1)
        print *, 'SRF LW net radiation        : ', rt%flns(n1)
        print *, 'SRF SW net radiation        : ', rt%fsns(n1)
        print *, 'SRF LW downwelling          : ', rt%flwds(n1)
        print *, 'SRF SW downwelling          : ', rt%fsds(n1)
        print *, 'TOA LW absorbed flux        : ', rt%fsnirt(n1)
        print *, 'TOA LW absorbed flux > 0.7u : ', rt%fsnirtsq(n1)
        print *, '########################################################'
        print *, 'TOA LW Clearsky absorbed    : ', rt%fsnrtc(n1)
        print *, 'TOA SW Clearsky absorbed    : ', rt%fsntc(n1)
        print *, 'TOA LW Clearsky net         : ', rt%flntc(n1)
        print *, 'SRF LW Clearsky net         : ', rt%flnsc(n1)
        print *, 'SRF SW Clearsky net         : ', rt%fsnsc(n1)
        print *, '########################################################'
        print *, 'SRF LW direct into soil     : ', rt%soll(n1)
        print *, 'SRF LW diffuse into soil    : ', rt%solld(n1)
        print *, 'SRF SW direct into soil     : ', rt%sols(n1)
        print *, 'SRF SW diffuse into soil    : ', rt%solsd(n1)
        print *, '########################################################'
        print *, 'TOA SW Aerosol effect       : ', rt%aeradfo(n1)
        print *, 'SRF SW Aerosol effect       : ', rt%aeradfos(n1)
        print *, 'TOA LW Aerosol effect       : ', rt%aerlwfo(n1)
        print *, 'SRF LW Aerosol effect       : ', rt%aerlwfos(n1)
        print *, '########################################################'
        print *, 'Radiation H/C in LW (max)   : ',  maxval(rt%qrl(:,n1))
        print *, 'Radiation H/C in LW (min)   : ',  minval(rt%qrl(:,n1))
        print *, 'Radiation H/C in SW (max)   : ',  maxval(rt%qrs(:,n1))
        print *, 'Radiation H/C in SW (min)   : ',  minval(rt%qrs(:,n1))
        print *, '########################################################'

        calday = calday + 1.0_rkx/4.0_rkx
        istep = istep + 1
      end do
    end do
  end do

  call deallocate_radtype(rt)
  call deallocate_aerosol( )
  deallocate(dum)
  deallocate(am,bm)
  deallocate(ai,bi)

  contains

  integer function ntimes_in_file(fname) result(nt)
    use netcdf
    character(len=*) , intent(in) :: fname
    integer :: ncid , ncstat , idimid
    ncstat = nf90_open(fname,nf90_nowrite,ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot open file '//trim(fname)
    ncstat = nf90_inq_dimid(ncid, 'time', idimid)
    if ( ncstat /= nf90_noerr ) stop 'No time dimension in '//trim(fname)
    ncstat = nf90_inquire_dimension(ncid, idimid, len=nt)
    if ( ncstat /= nf90_noerr ) stop 'rgrid dim time error in '//trim(fname)
    ncstat = nf90_close(ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot close file '//trim(fname)
  end function ntimes_in_file

  subroutine read_geolocation(fname,lat,lon)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: fname
    real(rkx) , dimension(n1:n2) , intent(out) :: lat , lon
    real(rk8) , dimension(reduced_points) :: rlat
    real(rk8) , dimension(nr) :: flat , flon
    integer , dimension(reduced_points) :: nlons
    real(rkx) , dimension(reduced_points) :: dlons
    integer :: ncid , ncvar , ncstat
    integer :: i , j , n , istart

    ncstat = nf90_open(fname,nf90_nowrite,ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot open file '//trim(fname)
    ncstat = nf90_inq_varid(ncid,'reduced_points',ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable reduced_points in file '//trim(fname)
    ncstat = nf90_get_var(ncid,ncvar,nlons)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot read variable reduced_points in file '//trim(fname)
    ncstat = nf90_inq_varid(ncid,'lat',ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable lat in file '//trim(fname)
    ncstat = nf90_get_var(ncid,ncvar,rlat)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot read variable lat in file '//trim(fname)
    ncstat = nf90_close(ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot close file '//trim(fname)
    do i = 1 , reduced_points
      dlons(i) = 360_rkx / real(nlons(i),rkx)
    end do
    istart = 1
    do i = 1 , reduced_points
      do j = istart , istart+nlons(i)-1
        flat(j) = rlat(i)
        !flon(j) = 0.5_rkx*dlons(i) + (j-istart) * dlons(i)
        ! According to EMOS, 0 is always present (WHY, THOUGH?)
        flon(j) = (j-istart) * dlons(i)
      end do
      istart = istart + nlons(i)
    end do
    do n = n1 , n2
      lat(n) = flat(n)
      lon(n) = flon(n)
    end do
  end subroutine read_geolocation

  subroutine read_vertical_grid(fname,am,bm,ai,bi)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: fname
    real(rkx) , dimension(nlev) :: am , bm
    real(rkx) , dimension(nlev+1) :: ai , bi
    integer :: ncid , ncvar , ncstat
    ncstat = nf90_open(fname,nf90_nowrite,ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot open file '//trim(fname)
    ncstat = nf90_inq_varid(ncid,'hyam',ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable hyam in file '//trim(fname)
    ncstat = nf90_get_var(ncid,ncvar,am)
    if ( ncstat /= nf90_noerr ) stop 'error read var hyam in '//trim(fname)
    ncstat = nf90_inq_varid(ncid,'hybm',ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable hybm in file '//trim(fname)
    ncstat = nf90_get_var(ncid,ncvar,bm)
    if ( ncstat /= nf90_noerr ) stop 'error read var hybm in '//trim(fname)
    ncstat = nf90_inq_varid(ncid,'hyai',ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable hyai in file '//trim(fname)
    ncstat = nf90_get_var(ncid,ncvar,ai)
    if ( ncstat /= nf90_noerr ) stop 'error read var hyai in '//trim(fname)
    ncstat = nf90_inq_varid(ncid,'hybi',ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable hybi in file '//trim(fname)
    ncstat = nf90_get_var(ncid,ncvar,bi)
    if ( ncstat /= nf90_noerr ) stop 'error read var hybi in '//trim(fname)
    ncstat = nf90_close(ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot close file '//trim(fname)
  end subroutine read_vertical_grid

  subroutine read_static_2di(fname,vname,f)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: fname , vname
    integer(ik4) , dimension(n1:n2) , intent(out) :: f
    integer :: ncid , ncvar , ncstat
    integer , dimension(2) :: istart , icount

    ncstat = nf90_open(fname,nf90_nowrite,ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot open file '//trim(fname)
    ncstat = nf90_inq_varid(ncid,vname,ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable '//trim(vname)//' in file '//trim(fname)
    istart(1) = n1
    istart(2) = 1
    icount(1) = n2-n1+1
    icount(2) = 1
    ncstat = nf90_get_var(ncid,ncvar,f,istart,icount)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot read variable '//trim(vname)//' in file '//trim(fname)
    ncstat = nf90_close(ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot close file '//trim(fname)
  end subroutine read_static_2di

  subroutine read_static_2d(fname,vname,f)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: fname , vname
    real(rkx) , dimension(n1:n2) , intent(out) :: f
    integer :: ncid , ncvar , ncstat
    integer , dimension(2) :: istart , icount

    ncstat = nf90_open(fname,nf90_nowrite,ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot open file '//trim(fname)
    ncstat = nf90_inq_varid(ncid,vname,ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable '//trim(vname)//' in file '//trim(fname)
    istart(1) = n1
    istart(2) = 1
    icount(1) = n2-n1+1
    icount(2) = 1
    ncstat = nf90_get_var(ncid,ncvar,f,istart,icount)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot read variable '//trim(vname)//' in file '//trim(fname)
    ncstat = nf90_close(ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot close file '//trim(fname)
  end subroutine read_static_2d

  subroutine read_field_2d(fname,vname,f,il,it)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: fname , vname
    real(rkx) , dimension(n1:n2) , intent(out) :: f
    integer(ik4) , intent(in) :: il, it
    integer :: ncid , ncvar , ncstat

    ncstat = nf90_open(fname,nf90_nowrite,ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot open file '//trim(fname)
    ncstat = nf90_inq_varid(ncid,vname,ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable '//trim(vname)//' in file '//trim(fname)
    if ( il == 0 ) then
      notime : block
        integer , dimension(2) :: istart , icount
        istart(1) = n1
        istart(2) = it
        icount(1) = n2-n1+1
        icount(2) = 1
        ncstat = nf90_get_var(ncid,ncvar,f,istart,icount)
        if ( ncstat /= nf90_noerr ) &
          stop 'cannot read variable '//trim(vname)//' in file '//trim(fname)
      end block notime
    else
      hastime : block
        integer , dimension(3) :: istart , icount
        istart(1) = n1
        istart(2) = il
        istart(3) = it
        icount(1) = n2-n1+1
        icount(2) = 1
        icount(3) = 1
        ncstat = nf90_get_var(ncid,ncvar,f,istart,icount)
        if ( ncstat /= nf90_noerr ) &
          stop 'cannot read variable '//trim(vname)//' in file '//trim(fname)
      end block hastime
    end if
    ncstat = nf90_close(ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot close file '//trim(fname)
  end subroutine read_field_2d

  subroutine read_field_3d(fname,vname,f,it)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: fname , vname
    real(rkx) , dimension(nlev,n1:n2) , intent(out) :: f
    integer(ik4) , intent(in) :: it
    real(rkx) , dimension(n1:n2,nlev) :: rf
    integer :: ncid , ncvar , ncstat
    integer :: i , k

    ncstat = nf90_open(fname,nf90_nowrite,ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot open file '//trim(fname)
    ncstat = nf90_inq_varid(ncid,vname,ncvar)
    if ( ncstat /= nf90_noerr ) &
      stop 'cannot find variable '//trim(vname)//' in file '//trim(fname)
    if ( it == 0 ) then
      notime : block
        integer , dimension(3) :: istart , icount
        istart(1) = n1
        istart(2) = 1
        istart(3) = 1
        icount(1) = n2-n1+1
        icount(2) = nlev
        icount(3) = 1
        ncstat = nf90_get_var(ncid,ncvar,rf,istart,icount)
        if ( ncstat /= nf90_noerr ) &
          stop 'cannot read variable '//trim(vname)//' in file '//trim(fname)
      end block notime
    else
      hastime : block
        integer , dimension(3) :: istart , icount
        istart(1) = n1
        istart(2) = 1
        istart(3) = it
        icount(1) = n2-n1+1
        icount(2) = nlev
        icount(3) = 1
        ncstat = nf90_get_var(ncid,ncvar,rf,istart,icount)
        if ( ncstat /= nf90_noerr ) &
          stop 'cannot read variable '//trim(vname)//' in file '//trim(fname)
      end block hastime
    end if
    ncstat = nf90_close(ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot close file '//trim(fname)
    do concurrent ( i = n1:n2, k = 1:nlev )
      f(k,i) = rf(i,k)
    end do
  end subroutine read_field_3d

  subroutine read_dimensions(fname,rgrid,nlev,reduced_points)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: fname
    integer(ik4) , intent(out) :: rgrid , nlev , reduced_points
    integer :: ncid , ncstat , idimid

    ncstat = nf90_open(fname,nf90_nowrite,ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot open file '//trim(fname)
    ncstat = nf90_inq_dimid(ncid, 'rgrid', idimid)
    if ( ncstat /= nf90_noerr ) stop 'No rgrid dimension in '//trim(fname)
    ncstat = nf90_inquire_dimension(ncid, idimid, len=rgrid)
    if ( ncstat /= nf90_noerr ) stop 'rgrid dim read error in '//trim(fname)
    ncstat = nf90_inq_dimid(ncid, 'nhym', idimid)
    if ( ncstat /= nf90_noerr ) stop 'No nhym dimension in '//trim(fname)
    ncstat = nf90_inquire_dimension(ncid, idimid, len=nlev)
    if ( ncstat /= nf90_noerr ) stop 'nhym dim read error in '//trim(fname)
    ncstat = nf90_inq_dimid(ncid, 'reduced_points', idimid)
    if ( ncstat /= nf90_noerr ) stop 'No lat dimension in '//trim(fname)
    ncstat = nf90_inquire_dimension(ncid, idimid, len=reduced_points)
    if ( ncstat /= nf90_noerr ) stop 'lat dim read error in '//trim(fname)
    ncstat = nf90_close(ncid)
    if ( ncstat /= nf90_noerr ) stop 'cannot close file '//trim(fname)
  end subroutine read_dimensions

  ! Computes saturation pressurre
  ! Reference:  Polynomial approximations from:
  !             Piotr J. Flatau, et al.,1992:  Polynomial fits to saturation
  !             vapor pressure.  Journal of Applied Meteorology, 31, 1507-1513.
  !
!DIR$ ATTRIBUTES FORCEINLINE :: pfesat
  pure elemental real(rkx) function pfesat(t,p) result(es)
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: t , p ! Temperature (K) , Pressure (Pa)

    real(rk8) :: td , t_limit , esat
    !
    ! For water vapor (temperature range 0C-100C)
    !
    real(rk8) , parameter :: a0 =  0.611213476e+03_rk8
    real(rk8) , parameter :: a1 =  0.444007856e+02_rk8
    real(rk8) , parameter :: a2 =  0.143064234e+01_rk8
    real(rk8) , parameter :: a3 =  0.264461437e-01_rk8
    real(rk8) , parameter :: a4 =  0.305903558e-03_rk8
    real(rk8) , parameter :: a5 =  0.196237241e-05_rk8
    real(rk8) , parameter :: a6 =  0.892344772e-08_rk8
    real(rk8) , parameter :: a7 = -0.373208410e-10_rk8
    real(rk8) , parameter :: a8 =  0.209339997e-13_rk8
    !
    ! For ice (temperature range -75C-0C)
    !
    real(rk8) , parameter :: c0 =  0.611123516e+03_rk8
    real(rk8) , parameter :: c1 =  0.503109514e+02_rk8
    real(rk8) , parameter :: c2 =  0.188369801e+01_rk8
    real(rk8) , parameter :: c3 =  0.420547422e-01_rk8
    real(rk8) , parameter :: c4 =  0.614396778e-03_rk8
    real(rk8) , parameter :: c5 =  0.602780717e-05_rk8
    real(rk8) , parameter :: c6 =  0.387940929e-07_rk8
    real(rk8) , parameter :: c7 =  0.149436277e-09_rk8
    real(rk8) , parameter :: c8 =  0.262655803e-12_rk8

    t_limit = t - tzero
    if ( t_limit > 100.0_rk8 ) t_limit = 100.0_rk8
    if ( t_limit < -75.0_rk8 ) t_limit = -75.0_rk8
    td = t_limit
    if ( td >= 0.0_rk8 ) then
      esat = a0 + td*(a1 + td*(a2 + td*(a3 + td*(a4 &
         + td*(a5 + td*(a6 + td*(a7 + td*a8)))))))
    else
      esat = c0 + td*(c1 + td*(c2 + td*(c3 + td*(c4 &
         + td*(c5 + td*(c6 + td*(c7 + td*c8)))))))
    end if
    es = real(min(esat,0.15_rk8*p),rkx) ! pa
  end function pfesat

!DIR$ ATTRIBUTES FORCEINLINE :: pfwsat
  pure elemental real(rkx) function pfwsat(t,p,e) result(ws) ! In kg/kg
!$acc routine seq
    implicit none
    real(rkx) , intent(in) :: t             ! Temperature (K)
    real(rkx) , intent(in) :: p             ! Pressure (Pa)
    real(rkx) , intent(in) , optional :: e  ! Saturated vapor pressure (Pa)
    real(rkx) :: es ! , qs
    if ( present(e) ) then
      es = e
    else
      es = pfesat(t,p)
    end if
    ws = ep2 * (es / (p - es))
    ! Bolton 1980
    ! qs = ep2 * es / (p - (0.378_rkx * es))
    ! ws = qs * ( d_one - qs)
  end function pfwsat

end program radiation_code

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
