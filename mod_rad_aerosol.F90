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
module mod_rad_aerosol

  use mod_constants
  use mod_rad_common
  use mod_dimensions
  use mod_sunorbit
  use mod_simple_plumes

  implicit none

  private

  public :: allocate_aerosol , deallocate_aerosol , init_aerosol , aeroppt

  real(rk8) , dimension(:) , pointer :: dnovrnr4 , latr4 , lonr4 , altr4
  real(rk8) , dimension(:) , pointer :: lambdaw
  real(rk8) , dimension(:,:) , pointer :: zr4 , dzr4
  real(rk8) , dimension(:,:,:) , pointer :: extprofr4
  real(rk8) , dimension(:,:,:) , pointer :: asyprofr4
  real(rk8) , dimension(:,:,:) , pointer :: ssaprofr4

  character(len=maxpath) :: macv2sp_hist , macv2sp_scen

  contains

  subroutine allocate_aerosol(n1,n2)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2
    allocate(dnovrnr4(n1:n2))
    allocate(latr4(n1:n2))
    allocate(lonr4(n1:n2))
    allocate(altr4(n1:n2))
    allocate(zr4(n1:n2,kz))
    allocate(dzr4(n1:n2,kz))
    allocate(lambdaw(nspi))
    allocate(extprofr4(n1:n2,kz,nspi))
    allocate(asyprofr4(n1:n2,kz,nspi))
    allocate(ssaprofr4(n1:n2,kz,nspi))
  end subroutine allocate_aerosol

  subroutine deallocate_aerosol
    implicit none
    deallocate(dnovrnr4)
    deallocate(latr4)
    deallocate(lonr4)
    deallocate(altr4)
    deallocate(lambdaw)
    deallocate(extprofr4)
    deallocate(asyprofr4)
    deallocate(ssaprofr4)
  end subroutine deallocate_aerosol

  subroutine init_aerosol(n1,n2,ht,lat,lon,macv2_path)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2
    real(rkx) , intent(in) , dimension(n1:n2) :: ht , lat , lon
    character(len=*) , intent(in) :: macv2_path
    integer(ik4) :: n
    do concurrent ( n = n1:n2 )
      latr4(n) = lat(n)
      lonr4(n) = lon(n)
      altr4(n) = ht(n)
    end do
    do n = 1 , nspi
      lambdaw(n) = (wavmin(n)+wavmax(n))*d_half*d_1000
    end do
    macv2sp_hist = trim(macv2_path)//'/MACv2.0-SP_v1.nc'
    macv2sp_scen = trim(macv2_path)//'/MACv2.0-SP_AIM-SSP3-Ref-SPA0.nc'
  end subroutine init_aerosol

  subroutine aeroppt(n1,n2,iyear,cday,z,dz,ftota3d,gtota3d,tauasc3d,tauxar3d)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2 , iyear
    real(rkx) , intent(in) :: cday
    real(rkx) , intent(in) , dimension(kz,n1:n2) :: z
    real(rkx) , intent(in) , dimension(kz,n1:n2) :: dz
    real(rkx) , intent(out) , dimension(kz,n1:n2,nspi) :: ftota3d
    real(rkx) , intent(out) , dimension(kz,n1:n2,nspi) :: gtota3d
    real(rkx) , intent(out) , dimension(kz,n1:n2,nspi) :: tauasc3d
    real(rkx) , intent(out) , dimension(kz,n1:n2,nspi) :: tauxar3d

    integer(ik4) :: m , n , k
    real(rk8) :: yf

    yf = iyear + cday/dayspy

    do concurrent ( n = n1:n2 , k = 1:kz )
      zr4(n,k) = z(k,n)
      dzr4(n,k) = dz(k,n)
    end do

    do m = 1 , nspi
      ! We have it for wl < 700.0
      if ( lambdaw(m) < 700.0_rkx ) then
        call sp_aop_profile(macv2sp_hist,macv2sp_scen,      &
                            kz,(n2-n1+1),lambdaw(m),        &
                            altr4,lonr4,latr4,yf,zr4,dzr4,  &
                            dnovrnr4,extprofr4(:,:,m),      &
                            ssaprofr4(:,:,m),asyprofr4(:,:,m))
      else
        extprofr4(:,:,m) = 0.0_rkx
        ssaprofr4(:,:,m) = 0.0_rkx
        asyprofr4(:,:,m) = 0.0_rkx
      end if
    end do
    ! Prepare for radiation code
    do concurrent ( m = 1:nspi, n = n1:n2 )
      do k = 1 , kz
        tauxar3d(k,n,m) = extprofr4(n,k,m)
        tauasc3d(k,n,m) = ssaprofr4(n,k,m)*extprofr4(n,k,m)
        gtota3d(k,n,m)  = asyprofr4(n,k,m)*ssaprofr4(n,k,m)*extprofr4(n,k,m)
        ftota3d(k,n,m)  = asyprofr4(n,k,m)*asyprofr4(n,k,m) * &
                          ssaprofr4(n,k,m)*extprofr4(n,k,m)
      end do
    end do
  end subroutine aeroppt

end module mod_rad_aerosol

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
