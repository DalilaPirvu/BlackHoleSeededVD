module constants
    use, intrinsic :: iso_c_binding
    implicit none

! Constants, general
    integer, parameter :: dl = kind(1.d0)
    real(dl), parameter :: twopi = 6.2831853071795864769252867665590
    complex(C_DOUBLE_COMPLEX), parameter :: iImag = (0._dl,1._dl)

! Lattice
    integer, parameter :: nLat = 512
    real(dl), parameter :: len = 100._dl
    real(dl), parameter :: dx = len/dble(nLat)
    real(dl), parameter :: dk = twopi/len
    integer, parameter :: nyq = nLat/2+1
    integer, parameter :: nTime = nLat

! Simulation
    integer, parameter :: nFld = 1
    integer, parameter :: nVar = 2*nFld*nLat+1
    integer, parameter :: minSim = 0
    integer, parameter :: nSims = 10

! Field definition
    real(dl), parameter :: phi0 = twopi/10._dl
    real(dl), parameter :: mean_phi = twopi*0.5_dl
    integer, parameter :: kcmin = 2, kcmax = nyq
    ! 1 multiplies mode amplitude by random deviate; 0 means no randomness
    integer, parameter :: randomQ = 1

! Potential Parameters
    integer, parameter :: frac = 2
    real(dl), parameter :: lambda = 0.028_dl
    real(dl), parameter :: m2Bare = 1._dl
    real(dl), parameter :: V0 = 1._dl
    real(dl), parameter :: delta = (m2Bare * phi0**2. / V0 + 1._dl)**0.5_dl
    ! interaction = 0 means constant mass term interaction
    ! interaction = 1 means x-dependent mass term interaction
    ! interaction = 2 means BEC potential with x-dependence
    real(dl), parameter :: interaction = 2._dl

end module constants
