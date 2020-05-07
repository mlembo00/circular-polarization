!! MLembo 19-06-17 !!
!! Power Spectra evaluation - Faraday conversion !!



!!! --------------------- START modules --------------------- !!! 

module READER

contains
  
  subroutine FILEREADER_camb(filename, x, y)
    implicit none
    
    character(LEN=1024), allocatable :: filename
    real(4) :: temp, z
    real(8), dimension(:) :: x, y
    integer :: l, i, nel, myunit

    nel=size(x)
    
    open (newunit=myunit, file=trim(filename), status='old', action='read')
    
    do i=1,nel
       read(myunit,*) l, temp, x(i), y(i), z
    enddo

    close(myunit)
  end subroutine FILEREADER_camb


!   subroutine FILEREADER_Ftable(filename, F)
!     implicit none
    
!     character(LEN=1024), allocatable :: filename
!     real(8), dimension(:,:) :: F
!     integer :: dim, i, myunit

!     dim=size(F,1)
    
!     open (newunit=myunit, file=trim(filename), status='old', action='read')
        
!     do i=1,dim
!       read(myunit,*) F(i,:)
!     enddo

!     close(myunit)
!   end subroutine FILEREADER_Ftable
  
! end module READER


  subroutine FILEREADER_EBtable(filename, diagonal, diagonal_p1, diagonal_m1)
    implicit none
    
    character(LEN=1024), allocatable :: filename
    real(8), dimension(:) :: diagonal, diagonal_p1, diagonal_m1
    integer :: i, nel, myunit

    nel=size(diagonal)
    !write(*,*) nel
    
    open (newunit=myunit, file=trim(filename), status='old', action='read')
    
    do i=1,nel
      read(myunit,*) diagonal(i), diagonal_p1(i), diagonal_m1(i)
      !write(*,*) diagonal(i), diagonal_p1(i), diagonal_m1(i)
    enddo

    close(myunit)
  end subroutine FILEREADER_EBtable

  subroutine FILEREADER_Vtable(filename, diagonal, diagonal_p1, diagonal_m1, diagonal_p2, diagonal_m2)
    implicit none
    
    character(LEN=1024), allocatable :: filename
    real(8), dimension(:) :: diagonal, diagonal_p1, diagonal_m1, diagonal_p2, diagonal_m2
    integer :: i, nel, myunit

    nel=size(diagonal)
    !write(*,*) nel
    
    open (newunit=myunit, file=trim(filename), status='old', action='read')
    
    do i=1,nel
      read(myunit,*) diagonal(i), diagonal_p1(i), diagonal_m1(i), diagonal_p2(i), diagonal_m2(i)
      !write(*,*) diagonal(i), diagonal_p1(i), diagonal_m1(i), diagonal_p2(i), diagonal_m2(i)
    enddo

    close(myunit)
  end subroutine FILEREADER_Vtable
  
end module READER



!!! ---------------------- END modules ---------------------- !!!


!!! --------------------- START program --------------------- !!!

program main
  use READER
#ifdef USEMPI  
  use mpi
#endif
  implicit none
  include "omp_lib.h"

!!! -------------- Variable Declaration - START ------------- !!!  
  
  !integer, parameter  :: dp = selected_real_kind(6, 9)
  complex, parameter :: ii=(0d0,1.d0)
  real(8), parameter :: pi=3.14159
  integer, parameter :: dp = SELECTED_REAL_KIND(14)
  double precision, parameter :: tol=1.0d-8

  integer :: feedback
  real :: timing
  integer :: t0, t1, clock_max, clock_rate  !Timing
  integer*4 :: argc
  integer(dp) :: l, i, l1max, lmax
  real(dp) :: f1_chi, f2_chi, fourpi, twopi, prefac, prefac_p1, prefac_m1, prefac_p2, prefac_m2
  real(dp), dimension(:), allocatable :: tildecl_EE, tildecl_BB, cl_EE, cl_BB, cl_VV
  real(dp), dimension(:), allocatable :: K11,K22p1,K22m1,K33p1,K33m1,K44,K44p2,K44m2
  character(len=1024), allocatable :: namefile1, namefile2, namefile3
  character(len=1024) :: parfile, infile_camb, outfile
  integer :: myunit

  namelist /input/ f1_chi, f2_chi, lmax, infile_camb, outfile, feedback

  fourpi = 4*pi
  twopi = 2*pi

!!! -------------- Variable Declaration - END --------------- !!! 


!!! ---------------- Initial setting - START ---------------- !!!
  
  argc = command_argument_count()     
  if (argc /= 1) then
    write(*,*) 'Run ./clV <param.ini>'
  else
    CALL get_command_argument(1,parfile)
  endif

  open(newunit=myunit,file=trim(parfile))
  read(myunit,nml=input)
  close(myunit)

  namefile1 = infile_camb
  namefile2 = 'fw3j_ee_bb.txt'
  namefile3 = 'fw3j_vv.txt'
   
  if (allocated(tildecl_EE)) deallocate(tildecl_EE)
  if (allocated(tildecl_BB)) deallocate(tildecl_BB)
  allocate(tildecl_EE(lmax+2),tildecl_BB(lmax+2))
  call FILEREADER_camb(namefile1,tildecl_EE,tildecl_BB)


  if (allocated(K11)) deallocate(K11)
  allocate(K11(lmax+2))
  if (allocated(K22p1)) deallocate(K22p1)
  allocate(K22p1(lmax+2))
  if (allocated(K22m1)) deallocate(K22m1)
  allocate(K22m1(lmax+2))
  call FILEREADER_EBtable(namefile2,K11,K22p1,K22m1) 

  if (allocated(K33p1)) deallocate(K33p1)
  allocate(K33p1(lmax+2))
  if (allocated(K33m1)) deallocate(K33m1)
  allocate(K33m1(lmax+2))
  if (allocated(K44)) deallocate(K44)
  allocate(K44(lmax+2))
  if (allocated(K44p2)) deallocate(K44p2)
  allocate(K44p2(lmax+2))
  if (allocated(K44m2)) deallocate(K44m2)
  allocate(K44m2(lmax+2))
  call FILEREADER_Vtable(namefile3,K44,K33p1,K33m1,K44p2,K44m2)


 write(*,*) 'checks ---->'
 write(*,*) 'fw3j_ee_bb:', K11(1), K22p1(1), K22m1(1)
 write(*,*) 'fw3j_vv:', K33p1(1), K33m1(1), K44(1), K44p2(1), K44m2(1)

  if (allocated(cl_EE)) deallocate(cl_EE)
  allocate(cl_EE(lmax))
  cl_EE = 0
  if (allocated(cl_BB)) deallocate(cl_BB)
  allocate(cl_BB(lmax))
  cl_BB = 0
  if (allocated(cl_VV)) deallocate(cl_VV)
  allocate(cl_VV(lmax))
  cl_VV = 0

  t0=0
  call system_clock(t0, clock_rate, clock_max)

!!! ----------------- Initial setting - END ----------------- !!!
  

if (feedback .ge. 1) then
  write(*,*) 'You are computing the power spectra with the following chi_ij parameters:'
  write(*,*) 'f1_chi, function of chi_ij with i .neq, j --> ',f1_chi
  write(*,*) 'f2_chi, function of chi_ij with i .eq, j --> ',f2_chi

endif


!!! ----------------- other parameters ----------------- !!!

  do i = 1,lmax-1, 1
    l = i+1
    prefac = twopi/(l*(l+1))
    prefac_p1 = twopi/((l+2)*(l+1))
    prefac_m1 = twopi/(l*(l-1))
    prefac_p2 = twopi/((l+2)*(l+3))
    prefac_m2 = twopi/((l-2)*(l-1))

    cl_EE(i) = prefac*tildecl_EE(i) + 1/fourpi * f1_chi * &
          &(K11(i)*prefac*tildecl_EE(i)+K22p1(i+1)*prefac_p1*tildecl_BB(i+1)+K22m1(i)*prefac_m1*tildecl_BB(i-1)) 
    
    cl_BB(i) = prefac*tildecl_BB(i) + 1/fourpi * f1_chi * &
          &(K11(i)*prefac*tildecl_BB(i)+K22p1(i+1)*prefac_p1*tildecl_EE(i+1)+K22m1(i)*prefac_m1*tildecl_EE(i-1)) 

    cl_VV(i) = 1/pi * f2_chi * (K44(i)*prefac*tildecl_BB(i)+K44p2(i+2)*prefac_p2*tildecl_BB(i+2)+K44m2(i)*prefac_m2*tildecl_BB(i-2)+&
                                 & K33p1(i+1)*prefac_p1*tildecl_EE(i+1)+K33m1(i)*prefac_m1*tildecl_EE(i-1))
  enddo !closing l loop

  open (newunit=myunit, file=trim(outfile))
  write(myunit, *) '## ell - EE - BB - VV ##'
    do l = 1, lmax-1, 1
      if (mod(l,10) .eq. 0 .and. (feedback .ge. 3)) write(*,*) 'Writing multipole: ', l+1
      write(myunit,'(1I6, 3E15.7)') l+1, cl_EE(l)*(l+1)*(l+2)/twopi, cl_BB(l)*(l+1)*(l+2)/twopi, cl_VV(l)*(l+1)*(l+2)/twopi     
    enddo
  close(myunit)


  deallocate(cl_EE, cl_BB, cl_VV, tildecl_EE, tildecl_BB)
    

  !Timing
  t1=0
  call system_clock(t1, clock_rate, clock_max)
  timing = real(t1 - t0)/real(clock_rate)
  if (feedback .ge. 1) write (*,*) 'Total time = ', timing, 'seconds'

  if (feedback .ge. 1) write(*,*) 'Finish ----> Cl from lmin =', 2, 'to lmax =', lmax

end program main

!!! ---------------------- END program ---------------------- !!!

