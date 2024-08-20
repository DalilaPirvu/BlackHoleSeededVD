#include "macros.h"
#define SMOOTH 1

program Gross_Pitaevskii_1d
  ! the following modules are incorporated:
  use, intrinsic :: iso_c_binding
  use gaussianRandomField
  use integrator
  use constants
  use eom
  implicit none

  real(dl), pointer :: time
  integer :: alph = 64
  integer :: sim
  integer, parameter :: inFile = 70, cpFile = 71
  real(dl), dimension(:,:), pointer :: fld

  fld(1:nLat,1:2) => yvec(1:nVar-1) ! store the field in yvec
  time => yvec(nVar) ! last position stores the time?
  call setup(nVar)

  call initialize_rand(93286123,12)
  do sim = minSim, nSims-1 ! run nSims simulations for the parameters, each with its output files
       call initialise_fields(fld)
       if (sim >= minSim) then
            call time_evolve(sim)
            print*, "Simulation ", sim+1, " in ", nSims, " done!"
       endif
  enddo
  print*, "Done! m2bare", m2Bare, "delta =", delta

contains

  subroutine initialise_fields(fld)
    real(dl), dimension(:,:), intent(inout) :: fld

    fld(:,1) = mean_phi
    fld(:,2) = 0._dl
    yvec(2*nLat+1) = 0._dl ! Add a tcur pointer here
    call initialize_fluctuations(fld)
  end subroutine initialise_fields


  subroutine time_evolve(sim)
    real(dl) :: dt, dtout
    integer :: i, j, outsize, sim

    dt = dx/alph
    if (dt > dx) print*, "Warning, violating Courant condition" !i.e. alph > 1
    outsize = alph
    dtout = dt*outsize

    do i = 1, nTime
       call output_fields(fld, dt, dtout, sim)
       do j = 1, outsize
          call gl10(yvec, dt)
       enddo
   enddo
  end subroutine time_evolve

  subroutine initialize_fluctuations(fld)
    real(dl), dimension(:,:), intent(inout) :: fld
    real(dl), dimension(1:size(fld(:,1)/2+1)) :: spec1, spec2
    real(dl), dimension(1:size(fld(:,1)/2+1)) :: fldspec, momspec
    real(dl), dimension(1:size(fld(:,1))) :: dfFld, dfMom
    real(dl) :: norm
    real(dl), allocatable, dimension(:) :: w2
    real(dl), allocatable, dimension(:,:) :: transfMatrix

    norm = 1._dl / phi0 / sqrt(4._dl * len)
    spec1 = 0._dl
    spec2 = 0._dl
    fldspec = 0._dl
    momspec = 0._dl

    call get_basis(w2, transfMatrix)

    spec1(kcmin:kcmax) = norm / w2(kcmin:kcmax)**0.25 / (exp(twopi * w2(kcmin:kcmax)**0.5 / lambda) - 1._dl)**0.5
    spec2(kcmin:kcmax) = spec1(kcmin:kcmax) * w2(kcmin:kcmax)**0.5

    ! change basis
    spec1 = matmul(transpose(transfMatrix), spec1(:nyq))
    spec2 = matmul(transpose(transfMatrix), spec2(:nyq))

    ! set mean field to zero; and all frequencies beyond nyquist, or k_max
    fldspec(:nyq) = spec1(:nyq)
    momspec(:nyq) = spec2(:nyq)

    call generate_1dGRF(dfFld, fldspec(:nyq))
    call generate_1dGRF(dfMom, momspec(:nyq))
    
    fld(:,1) = fld(:,1) + dfFld(:)
    fld(:,2) = fld(:,2) + dfMom(:)
  end subroutine initialize_fluctuations



  subroutine output_fields(fld, dt, dtout, sim)
    real(dl), dimension(1:nLat, 1:2) :: fld
    real(dl) :: dt, dtout
    logical :: o
    integer :: ii, sim
    integer, parameter :: oFile = 98

    inquire(file='/gpfs/dpirvu/dilatonBH/BEC_sim'//trim(str(sim))//'_N'//trim(str(nLat))//'_L'//trim(str(int(len)))//'_lambda'//trim(real_str(lambda))//'_fields.dat', opened=o)
    if (.not.o) then
       open(unit=oFile,file='/gpfs/dpirvu/dilatonBH/BEC_sim'//trim(str(sim))//'_N'//trim(str(nLat))//'_L'//trim(str(int(len)))//'_lambda'//trim(real_str(lambda))//'_fields.dat')
       write(oFile,*) "# Lattice Parameters nLat = ", nLat, "dx = ", dx, "lenLat = ", len, "nyquist = ", nyq
       write(oFile,*) "# Time Stepping parameters dt = ", dt, "dt_out = ", dtout
       write(oFile,*) "# Other Parameters m2Bare = ", m2Bare, "dk = ", dk, "lambda = ", lambda
    endif

    do ii = 1, nLat
       write(oFile,*) fld(ii,:)
    enddo
  end subroutine output_fields


  subroutine setup(nVar)
    integer, intent(in) :: nVar
    call init_integrator(nVar)
    call initialize_transform_1d(tPair,nLat)
    call initialize_transform_1d(vPair,nLat)
  end subroutine setup

  character(20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  end function str

  character(20) function real_str(k)
    real(dl), intent(in) :: k
    write (real_str, '(f12.3)') k
    real_str = adjustl(real_str)
  end function real_str

  subroutine get_basis(w2, transfMatrix)
  !!!! Note: these have accuracy less than machine precision of 15 digits; 
     real(dl), allocatable, dimension(:), intent(inout) :: w2
     real(dl), allocatable, dimension(:,:), intent(inout) :: transfMatrix
     integer :: ii, jj, kk, nn, mm, ll
          open(unit=99,file='../data/frequencies_N'//trim(str(nLat))//'_L'//trim(str(int(len)))//'_lambda'//trim(real_str(lambda))//'.txt', status='old', action='read')
     read(99, *) nn
     if (nn == nyq) then
         allocate(w2(nn))
     else
         print*, 'Data file not compatible.'
     endif
     do ii = 1, nn
         read(99,*) w2(ii)
     enddo
     close(99)
    open(unit=88,file='../data/transfMatrix_N'//trim(str(nLat))//'_L'//trim(str(int(len)))//'_lambda'//trim(real_str(lambda))//'.txt', status='old', action='read')
     read(88, *) mm
     read(88, *) ll
     allocate(transfMatrix(mm,ll))
     do jj = 1, mm
        do kk = 1, ll
           read(88, *) transfMatrix(jj,kk)
        enddo
     enddo
     close(88)
  end subroutine get_basis


end program Gross_Pitaevskii_1d
