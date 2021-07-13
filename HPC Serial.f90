program serial
    implicit none
    integer, parameter :: dp=selected_real_kind(15,300)
    integer, parameter :: Nr = 100, Nz = 200, Phi_0 = 1000
    integer, parameter :: r_small = 10, z_small = 50
    real(dp) :: r_in = 10._dp, R_out = 100._dp, x
    integer, parameter :: delta = 1.0e-3_dp
    integer :: i,j,k,istat
    real(dp) :: U
    real(dp), allocatable, dimension(:,:) :: Phi_new
    
    !allocate two required Phi arrays
    allocate(Phi_new(0:Nr,0:Nz), stat=istat)
    if (istat/=0) stop 'Error in allocating Phi array'

    !allocate(Phi_new(0:Nr,0:Nz), stat=istat)
    !if (istat/=0) stop 'Error in allocating Phi array'
    
    !Set boundary conditions
    Phi_new(:,:) = 0._dp
    !Set Plane conditiions for z=0
    
    !Start Big loop
    iterations: do j = 1,20000

    z_plane: do i = r_small+1,Nr
     x = i
     Phi_new(i,0) = Phi_0*((log(R_out)-log(x))/(log(R_out)-log(r_in)))
    end do z_plane
    
    Phi_new(0:r_small,0:z_small) = Phi_0

    Main_grid: do i = 1,Nr
      do k = 1,Nz
        U = 0.25_dp*(Phi_new(i+1,k)+Phi_new(i-1,k)+Phi_new(i,k+1)&
        & +Phi_new(i,k-1)) + (1/8*k)*(Phi_new(i,k+1) - Phi_new(i,k-1))
        Phi_new(i,k) = Phi_new(i,k) + 1.5_dp*(U - Phi_new(i,k))
      end do
    end do Main_grid

    end do iterations

  
    print *,'A=', Phi_new(0,75), 'B=', Phi_new(60,50), 'C=', Phi_new(50,125)
    !call pgm()

    deallocate(Phi_new, stat=istat)
    if(istat/=0) stop 'Error in deallocating Phi_new'
 
    !deallocate(Phi_new, stat=istat)
    !if(istat/=0) stop 'Error in deallocating Phi_new'


end program
