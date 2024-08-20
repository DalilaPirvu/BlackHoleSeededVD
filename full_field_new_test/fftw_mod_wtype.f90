!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for implementing 1-dimensional derivatives and antiderivatives using pseudospectral methods
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> @author
!> Jonathan Braden, University College London
!>
! DESCRIPTION
!> @brief
!> An implementation of pseudospectral based one-dimensional derivatives and antiderivatives
!>
!> This module is used to compute derivatives and antiderivatives of one-dimensional fields
!> using pseudospectral approximations.
!> The required FFT's are implemented using the open source package FFTW.
!>
!> A typical usage is as follows
!> @code{.f90}
!>  type(transformPair1D) :: tPair
!>  real(C_DOUBLE), parameter :: dk = 1._dl
!>  call initialize_transform_1d(tPair, 1024)
!>  tPair%realSpace(:) = initialField(:)
!>  call derivative_n_1d(tPair,dk)
!> @endcode
!> @warning When calling any of the derivative subroutines, the field to be differentiated must already be stored in tPair%realSpace.
!> @warning The differentiated function is stored in tPair%realSpace after the subroutine call.
!>
!> The FFT convention used in FFTW is (in 1D with N lattice sites)
!> \f[ g_k = \sum_l e^{i2\pi k l / N} f(x_l) \f]
!> so that the nth derivative multiplies the Fourier amplitudes by \f[A^{(n)} = \left(i dk \right)^n \f].
!> As well, FFTW computes an unnormalised transform,
!> \f[ f(x_l) = \sum_k e^{-i2\pi kl /N} g_k \f]
!> so that the final answer must me multiplied by \f$ N_{\mathrm{lat}}^{-1} \f$
!> where \f$ N_{lat} \f$ is the total number of lattice sites.
!> In the actual implementation of the derivative routines below, these are implemented as an normalisation factor \f[ \mathrm{norm} = N_{lat}^{-1}\left( i dk \right)^n \f]
!> multiplying each Fourier amplitude.
!>
!> Additionally, for inverse derivatives (which constitute elliptic equations which require boundary conditions),
!> the boundary conditions are automatically assumed to be periodic since this module is based on Fourier pseudospectral methods.
!> In particular, any mean value for a field we attempt to compute and inverse derivative of will be ignored, since it will not posess continuous solutions
!> within the space of continuous functions.
!> @warning The implicit periodic boundary conditions must be kept in mind when using this module to solve elliptic equations.
!>
! TO DO
!> @todo
!> @arg Include sample use cases in overall module documentation.
!> @arg Implement higher dimensional calculations
!> @arg Include required OpenMP directives
!> @arg Implement MPI communication for parallelisation
!> @arg Use separate flags for FFTW OpenMP parallellism, and parallellism in the loops
!> @arg Allow user to specify level of tuning of FFTW through a flag the user passes in during initialisation
!> @arg Link documentation of derivative functions to a single function where I'll show Fourier conventions, etc.
!>      Or else document this in the start of the module?  Improve the description of the calculation above
!> @arg Since higher dimensions allow for vector derivatives (ie. the gradient), implement a type that allows for this to be stored easily instead of the ugly rigging that is being used now.  This will introduce additional tuning parameters based on convenience, speed and memory
!> @arg Finish code example in module documentation
!> @arg Either include in this module, or make a separate convenience module, where a transform pair type is defined along with some subroutines that call this module
!>      Then have simple calls like derivative(f, df, dk) without referencing the transform, with f the input field, df the differentiated field and dk the fourier spacing
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module fftw3
  use, intrinsic :: iso_c_binding
#ifdef USEOMP
  use omp_lib
#endif
  use constants
  implicit none
  include 'fftw3.f03'

  type transformPair1D
     integer :: nx, nnx
     real(C_DOUBLE), pointer :: realSpace(:)
     complex(C_DOUBLE_COMPLEX), pointer :: specSpace(:)
     type(C_PTR) :: planf, planb
     type(C_PTR), private :: rPtr, sPtr
  end type transformPair1D
  
contains
   
!****************************!
! Setup Transformation Pairs !
!****************************!

  !something to initialise openmpi
  subroutine boot_openmp(nThread)
  integer, optional, intent(in) :: nThread
  integer :: errorOMP
#ifdef USEOMP
    errorOMP = fftw_init_threads()
    if (errorOMP == 0) then
       print*,"Error initializing OpenMP threading for FFTW"
       stop
    endif
    if (present(nThread)) then
       errorOMP = nThread
    else
       print*,"Defaulting to OMP_NUM_THREADS environment variable"
       errorOMP = omp_get_max_threads()
    endif
    call fftw_plan_with_nthreads(errorOMP)
    print*,"FFTW booted using ",errorOMP," threads"
#endif
  end subroutine boot_openmp

#define RSPACE1D create_transform_1d%realSpace
#define SSPACE1D create_transform_1d%specSPACE

  function create_transform_1d(n)
    type(transformPair1D) :: create_transform_1d
    integer, intent(in) :: n
    call allocate_1d_array(n, RSPACE1D, SSPACE1D, create_transform_1d%rPtr, create_transform_1d%sPtr)
    create_transform_1d%nx = n
    create_transform_1d%planf = fftw_plan_dft_r2c_1d(n, RSPACE1D, SSPACE1D, FFTW_MEASURE)
    create_transform_1d%planb = fftw_plan_dft_c2r_1d(n, SSPACE1D, RSPACE1D, FFTW_MEASURE)
  end function create_transform_1d

  subroutine initialize_transform_1d(this, n) !called with (tPair, nLat)
  ! Creates a tPair which is a transformPair1D type of object
  ! A transform pair is for doing 1D transforms using FFTW stuff
    type(transformPair1D), intent(out) :: this
    integer, intent(in) :: n
  
  !the % takes out that one argument of the object
    this%nx = n
    this%nnx = n/2 + 1
    call allocate_1d_array(n, this%realSpace, this%specSpace, this%rPtr , this%sPtr)
    
    this%planf = fftw_plan_dft_r2c_1d(n, this%realSpace, this%specSpace, FFTW_MEASURE)
    this%planb = fftw_plan_dft_c2r_1d(n, this%specSpace, this%realSpace, FFTW_MEASURE)
  end subroutine initialize_transform_1d

  subroutine destroy_transform_1d(this)
  ! Destroys the tPair it is called with 
    type(transformPair1D), intent(inout) :: this
    call fftw_destroy_plan(this%planf)
    call fftw_destroy_plan(this%planb)
    call fftw_free(this%rPtr)
    call fftw_free(this%sPtr)
  end subroutine destroy_transform_1d
  
  subroutine allocate_1d_array(L, arr, Fk, fptr, fkptr)
    integer, intent(in) :: L
    real(C_DOUBLE), pointer :: arr(:)
    complex(C_DOUBLE_COMPLEX), pointer :: Fk(:)
    type(C_PTR) :: fptr, fkptr
    integer :: LL
    
    LL = L/2+1
    fptr = fftw_alloc_real(int(L, C_SIZE_T))
    call c_f_pointer(fptr, arr, [L])
    
    fkptr = fftw_alloc_complex(int(L, C_SIZE_T))
    call c_f_pointer(fkptr, Fk, [LL])
  end subroutine allocate_1d_array

  subroutine destroy_arrays_1d(planf,planb,fptr,fkptr)
    type(C_PTR) :: planf, planb
    type(C_PTR) :: fptr, fkptr
    call fftw_destroy_plan(planf)
    call fftw_destroy_plan(planb)
    call fftw_free(fptr)
    call fftw_free(fkptr)
  end subroutine destroy_arrays_1d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine laplacian_1d_wtype(tPair, dk)
    type(transformPair1D), intent(inout) :: tPair
    real(dl), intent(in) :: dk
    integer :: i
    complex(C_DOUBLE_COMPLEX) :: norm

    norm = -dk**2 / dble(tPair%nx)
    call fftw_execute_dft_r2c(tPair%planf, tPair%realSpace, tPair%specSpace)
    do i=1,tPair%nnx
       tPair%specSpace(i) = (i-1)**2 * norm * tPair%specSpace(i) !the (i-1) is squared compared to above
    enddo
    call fftw_execute_dft_c2r(tPair%planb, tPair%specSpace, tPair%realSpace)
  end subroutine laplacian_1d_wtype

end module fftw3