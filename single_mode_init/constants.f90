module constants
    use, intrinsic :: iso_c_binding
    implicit none

! General
    integer, parameter :: dl = kind(1.d0)
    real(dl), parameter :: twopi = 6.2831853071795864769252867665590
    complex(C_DOUBLE_COMPLEX), parameter :: iImag = (0._dl,1._dl)

! Lattice parameters
    integer, parameter :: nLat = 512
    real(dl), parameter :: len = 1000._dl
    real(dl), parameter :: dx = len/dble(nLat)
    real(dl), parameter :: dk = twopi/len
    integer, parameter :: nyq = nLat/2+1
    integer, parameter :: nTime = nLat

! Simulation parameters
    integer, parameter :: nFld = 1
    integer, parameter :: nVar = 2*nFld*nLat+1
    integer, parameter :: minSim = 0, nSims = nyq-1

! Field definition
    real(dl), parameter :: phi0 = 1._dl
!    integer, parameter :: kcmin = 2, kcmax = nyq
    integer, parameter :: kcmin = 1, kcmax = 1
    integer, parameter :: randomQ = 0 ! 1 multiplies mode amplitude by random deviate; 0 means no randomness

! Potential Parameters
    integer, parameter :: frac = 2
    real(dl), parameter :: lambda = 0.05_dl
    real(dl), parameter :: m2Bare = 1._dl
    real(dl), parameter :: interaction = 1._dl
!    real(dl), parameter :: BEClambda = 10

end module constants
