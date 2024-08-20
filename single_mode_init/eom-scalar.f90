#include "macros.h"
#include "fldind.h"
#define FIELD_TYPE Field_Model
module eom
  use constants
  use fftw3
  implicit none

  real(dl), dimension(1:nVar), target :: yvec
  type(transformPair1D) :: tPair
  type(transformPair1D) :: vPair

contains
  
  subroutine derivs(yc, yp)
    real(dl), dimension(:), intent(in) :: yc
    real(dl), dimension(:), intent(out) :: yp
    integer :: x

    yp(TIND) = 1._dl  ! Uncomment to track time as a variable
    yp(FLD) = yc(DFLD)

    if (interaction == 0.) then !! free field: m^2 phi(x) term in eom
      yp(DFLD) = - m2Bare * yc(FLD)

   else if (interaction == 1.) then !! omega(x) m^2 phi(x) term in the eom 
      do x = 1, nLat
          if (x < nLat/2) then
               yp(nLat+x) = - m2Bare * Omega1(x - nLat/2) * yc(x)
          else
               yp(nLat+x) = - m2Bare * Omega2(x - nLat/2) * yc(x)
          endif
      enddo
   endif

    tPair%realSpace(:) = yc(FLD)
    call laplacian_1d_wtype(tPair,dk)
    yp(DFLD) = yp(DFLD) + tPair%realSpace(:)
  end subroutine derivs

  real(dl) function Omega1(x)
    integer, intent(in) :: x
    Omega1 = 1._dl / (1._dl + exp(-2._dl * lambda * (x + nLat/frac/2._dl)))
  end function Omega1
  
  real(dl) function Omega2(x)
    integer, intent(in) :: x
    Omega2 = 1._dl / (1._dl + exp(2._dl * lambda * (x - nLat/frac/2._dl)))
  end function Omega2

end module eom
