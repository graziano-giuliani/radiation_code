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
module mod_dimensions

  use mod_constants

  implicit none

  public

  integer :: kz
  integer :: kzp1
  integer :: kzp2
  integer :: kzp3
  integer :: kzp4
  integer :: kzm1
  integer :: kzm2
  integer :: nr

  contains

  subroutine init_dimensions(nrmax,kmax)
    implicit none
    integer , intent(in) :: nrmax , kmax
    kz = kmax
    kzp1 = kz+1
    kzp2 = kz+2
    kzp3 = kz+3
    kzp4 = kz+4
    kzm1 = kz-1
    kzm2 = kz-2
    nr = nrmax
  end subroutine init_dimensions

end module mod_dimensions

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
