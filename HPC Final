Program b_mpi
  use mpi
  implicit none
  integer,parameter :: dp=selected_real_kind(15,300)
  integer,parameter :: Nr=100, Nz=200, Phi_0=1000
  integer, parameter :: r_small = 10, z_small = 50
  integer :: r,z,k,ierror
  real(dp),allocatable,dimension(:,:) :: Phi, Final_Phi
  real(dp) :: start,finish, r_in=10._dp, R_out= 100._dp, x,U

  integer :: my_rank, num_procs, up_rank,down_rank
  integer :: lbound,ubound,divide
  integer,dimension(1:MPI_STATUS_SIZE) :: send_status,recv_status
  integer :: downsend_request,upsend_request,downrecv_request,uprecv_request

  call MPI_init(ierror)
  if (ierror/=MPI_success) stop 'Error in MPI_init'

  call MPI_comm_size(MPI_comm_world,num_procs,ierror)
  if (ierror/=MPI_success) stop 'Error in MPI_comm_size'

  call MPI_comm_rank(MPI_comm_world,my_rank,ierror)
  if (ierror/=MPI_success) stop 'Error in MPI_comm_rank'

  divide = Nz/num_procs
  if (my_rank==0) then
   lbound = (my_rank * divide) 
   ubound = lbound + divide 
  end if
 
  if (my_rank>0) then
   lbound = (my_rank * divide) + 1
   ubound = lbound + divide-1
  end if

  if (mod(Nz,num_procs)/=0) then
   stop 'Grid cannot be divided into this number of processors'
  end if

  if (my_rank==0) then
   allocate(Phi(0:Nr,lbound:ubound+1), stat = ierror)
   if (ierror/=0) stop 'Error in allocating Phi'

  else if (0<my_rank.and.my_rank<num_procs-1) then
   allocate(Phi(0:Nr,lbound-1:ubound+1), stat=ierror)
   if (ierror/=0) stop 'Error in allocating Phi'
  
  else
   allocate(Phi(0:Nr,lbound-1:ubound+1), stat = ierror)
   if (ierror/=0) stop 'Error in allocating Phi'
  end if
  
  allocate(Final_Phi(0:Nr,0:Nz), stat = ierror)
  if (ierror/=0) stop 'Error in allocating Final_Phi'

  Phi = 0._dp

  up_rank = my_rank + 1
  down_rank = my_rank - 1
 
  if (lbound<=z_small) then
    if (ubound<=z_small) then
    Phi(0:r_small,lbound:ubound) = Phi_0
    else
    Phi(0:r_small,lbound:z_small) = Phi_0
    end if
  end if
  
  

  start = mpi_wtime()
 
  z_plane: do r=r_small+1,Nr
   x=r
   Phi(r,0) = Phi_0*((log(R_out)-log(x))/(log(R_out)-log(r_in)))
  end do z_plane

  iterations: do k=1,20000

  !Upward Comms
  if (my_rank<num_procs-1) then
   call MPI_ISSEND(Phi(:,ubound),Nr+1,MPI_Double_precision,up_rank,1,MPI_comm_world,upsend_request,ierror)
   if (ierror/=MPI_success) stop 'Error in MPI_ISSEND upwards'
   
   call MPI_irecv(Phi(:,ubound+1),Nr+1,MPI_double_precision,up_rank,2,MPI_comm_world,uprecv_request,ierror)
   if (ierror/=MPI_success) stop 'Error in MPI_Irecv upwards'
  end if

  !Downward Comms
  if (my_rank>0) then
   call MPI_irecv(Phi(:,lbound-1),Nr+1,MPI_double_precision,down_rank,1,MPI_comm_world,downrecv_request,ierror)
   if (ierror/=MPI_success) stop 'Error in MPI_irecv downwards'
 
   call MPI_ISSEND(Phi(:,lbound),Nr+1,MPI_double_precision,down_rank,2,MPI_comm_world,downsend_request,ierror)
   if (ierror/=MPI_success) stop 'Error in MPI_ISSEND downwards'
  end if

  r_zero: do z = lbound+1,ubound-1
   U = (2._dp/3._dp)*Phi(1,z) + (1._dp/6._dp)*(Phi(0,z+1) + Phi(0,z-1))
   Phi(0,z) = (Phi(0,z) + 1.2_dp*(U-Phi(0,z)))
  end do r_zero

  Main_grid: do r = 1,Nr
   do z = lbound+1,ubound-1
    U = 0.25_dp*(Phi(r+1,z) +Phi(r-1,z) + Phi(r,z+1) + Phi(r,z-1)) &
    & + (1/(8*r)) * (Phi(r+1,z) - Phi(r-1,z))
    Phi(r,z) = Phi(r,z) + 1.2_dp*(U-Phi(r,z))
   end do 
  end do Main_grid

  if (my_rank<num_procs-1) then
   call MPI_wait(upsend_request,send_status,ierror)
   if (ierror/=MPI_success) stop 'Error in MPI_wait upsend'

   call	MPI_wait(uprecv_request,recv_status,ierror)
   if (ierror/=MPI_success) stop 'Error in MPI_wait uprecv'
  end if

  if (my_rank>0) then
   call	MPI_wait(downrecv_request,recv_status,ierror)
   if (ierror/=MPI_success) stop 'Error in MPI_wait downrecv'

   call	MPI_wait(downsend_request,send_status,ierror)
   if (ierror/=MPI_success) stop 'Error in MPI_wait downsend'
  end if
 
  r_zero: do z = lbound,ubound,(ubound-lbound)
   U = (2._dp/3._dp)*Phi(1,z) + (1._dp/6._dp)*(Phi(0,z+1) + Phi(0,z-1))
   Phi(0,z) = (Phi(0,z) + 1.2_dp*(U-Phi(0,z)))
  end do

  if (my_rank>0) then
   do r = 1,Nr
    U = 0.25_dp*(Phi(r+1,lbound) + Phi(r-1,lbound) + Phi(r,lbound+1) + Phi(r,lbound-1)) &
    & + (1/(8*r))*(Phi(r+1,lbound)+Phi(r-1,lbound))
    Phi(r,1) = Phi(r,1) + 1.2_dp*(U-Phi(r,1))
   end do
  end if

  if (my_rank<num_procs-1) then
   do r = 1,Nr
    U = 0.25_dp*(Phi(r+1,ubound) + Phi(r-1,ubound) + Phi(r,ubound-1) + Phi(r,ubound+1)) &
    & + (1/(8*r))*(Phi(r+1,ubound) + Phi(r-1,ubound))
    Phi(r,ubound) = Phi(r,ubound) + 1.2_dp*(U-Phi(r,ubound))
   end do
  end if

  if (lbound<=z_small) then
   if (ubound<=z_small) then
    Phi(0:r_small,lbound:ubound) = Phi_0
   else
    Phi(0:r_small,lbound:z_small) = Phi_0
   end if
  end if

   end do iterations
  
    call MPI_gather(Phi(0:Nr,lbound:ubound),(Nr+1)*((ubound-lbound)+1),MPI_double_precision, &
    & Final_Phi(0:Nr,0:Nz),(Nr+1)*((ubound-lbound)+1),MPI_double_precision,0,MPI_comm_world,ierror)
    if (ierror/=MPI_success) stop 'Error in MPI_gather'
    
    finish = MPI_wtime()

    if (my_rank==0) then
     print *, 'A=', my_rank, Final_Phi(0,75), 'B=', my_rank, Final_Phi(60,50), 'C=', my_rank, Final_Phi(50,125)
     print *, 'Time taken=', finish-start,' s'
     call pgm(Final_Phi)
    end if

    

    deallocate(Phi,stat=ierror)
    if (ierror/=0) stop 'Error in deallocating Phi array'

    deallocate(Final_Phi,stat=ierror)
    if (ierror/=0) stop 'Error in deallocating Final_Phi array'

    call MPI_finalize(ierror)
    if (ierror/=MPI_success) stop 'Error in MPI_finalize'


contains
subroutine pgm(grid)
  implicit none
  integer, parameter :: dp = selected_real_kind(15,300)
  !set up the pgmfile
  integer,parameter :: Nr=100,Nz=200
  integer :: i,j,k,out_unit,Phi_max,Phi_min,max_greys,ierr
  integer,dimension(:,:),allocatable :: pixels
  real(dp), dimension(:,:) :: grid

  allocate(pixels(0:Nr,0:Nz),stat=ierr)
  if (ierr/=0) stop 'Error in allocating pixels'

  max_greys=255
  Phi_max=maxval(grid)
  Phi_min=minval(grid)

  do j=1,Nz
     do i=1,Nr
        pixels(i,j)=int((grid(i,j)-Phi_min)*max_greys/(Phi_max-Phi_min)) !Tmin<T<Tmax
     end do
  end do

  out_unit=10
  open(file='b_mpi.pgm',unit=out_unit)
  write (out_unit,11) 'P2'                 !pgm magic number
  write (out_unit,12) Nr,Nz                !width, height
  write (out_unit,13) max_greys            !max gray value
  do j=1,Nz
     do i=1,Nr
        write (out_unit,14) (pixels(i,j))  !each line < 70 chars
     end do
     write (out_unit,14) (pixels(k,j),k=i,Nr)
  end do
  close (unit=out_unit)

11 format(a2)
12 format(i10,1x,i10)
13 format(i10)
14 format (15(1x,i3))

end subroutine

end program 
