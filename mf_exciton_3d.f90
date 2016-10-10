!---------------!
program exciton3d
!---------------!
implicit none


! main variables
real(8) :: amp, t_e, t_h, eps_0_e, eps_0_h, u_0, alfa, beta, tolerance
integer :: counter, lat_size, a, lat_const, statistical_sample
! work variables
integer :: i,j,k,l,m,n,ii,num,conv_iteration, NN
integer :: hole_index, el_index
real(8) :: error_e, error_h, Eg_0, Eg_coulomb, en_addit
real(8) :: xx, norm_e, norm_h, E_binding, e_bind_avg, e_bind_std 
real(8) :: coulomb
real :: t1, t2
real(8), allocatable :: onsit_e(:,:),onsit_h(:,:),hop_e(:,:),hop_h(:,:)
real(8), allocatable :: n_e(:,:,:), n_h(:,:,:), pd_e(:,:,:), pd_h(:,:,:), &
           H_e(:,:), H_h(:,:), n_e_prev(:,:,:), n_h_prev(:,:,:), fluc(:,:,:), &
           eig_val_e(:), eig_val_h(:), eig_vec_e(:,:), eig_vec_h(:,:), &
           e_bind_array(:)
!real(8), allocatable :: m0(:,:),m1(:,:),m2(:,:),eig(:) 
real :: bin(300)
logical :: output_display, random_seed_from_clock
character(40) :: tag
character(40) :: tag2


amp = 0.                         ! Fluctuation amplitude
lat_size = 12                    ! Lattice size is NxN  
t_e = 0.20                       ! hopping e, every other calculation is done on this scale
t_h = 0.14                       ! hopping h 
eps_0_e = 2.0                    ! mid band energy (onsite) of conduction band
eps_0_h = -1.7                   ! mid band energy (onsite) of valence band 
u_0 = 0.8                        ! coulomb parameter
alfa = 0.75                      ! inverse of Bohr radius of exciton
a = 1                            ! nearest neighbour distance
lat_const = 1                    ! lattice constant     
beta = .24                       ! scf memory parameter
tolerance = 1.e-4                ! convergence charge density difference tolerance
statistical_sample = 1           ! number of samples
output_display = .true.
random_seed_from_clock = .true.

call CPU_TIME( t1 ) 
write(tag,*) amp
if (output_display) then  !Doubt
    write(*,*) "# fluc_ampl = ", amp, " lat_size = ", lat_size, " u_0 = ", u_0
    write(*,*) "# beta = ", beta, " tolerance = ", tolerance
endif

open(33,file='logfile.dat',status='replace')
write(33,*) "Logfile"
write(33,*) "t_e = ", t_e, "t_h = ", t_h
write(33,*) "eps_0_e = ", eps_0_e, "eps_0_h = ", eps_0_h
write(33,*) "fluctuation amplitude = ", amp, " lat_size = ", lat_size, " u_0 =", u_0
write(33,*) "beta = ", beta, " tolerance = ", tolerance
write(33,*) " "
write(33,*) " #*#*#*#*#*#*#*#*#*#*#*#*#*#*#"

NN = lat_size
hole_index = NN**3
el_index = 1
allocate(e_bind_array(statistical_sample))
counter = 1  
do while (counter <= statistical_sample) ! statistical sample loop
    !write(*,*) "ciao", 3**2
    allocate(onsit_e(lat_size**3,lat_size**3)); allocate(onsit_h(lat_size**3,lat_size**3)) ! matrix construction for 3d cubic lattice
    onsit_e = 0.d0 ; onsit_h = 0.d0
    allocate(hop_e(lat_size**3,lat_size**3)); allocate(hop_h(lat_size**3,lat_size**3))
    hop_e = 0.d0 ; hop_h = 0.d0
    allocate(n_e_prev(lat_size,lat_size,lat_size)); allocate(n_h_prev(lat_size,lat_size,lat_size))
    n_e_prev = 0.d0 ; n_h_prev = 0.d0
    allocate(n_e(lat_size,lat_size,lat_size)); allocate(n_h(lat_size,lat_size,lat_size))
    allocate(pd_e(lat_size,lat_size,lat_size)); allocate(pd_h(lat_size,lat_size,lat_size))
    allocate(H_e(lat_size**3,lat_size**3)); allocate(H_h(lat_size**3,lat_size**3))
    H_e = 0.d0 ; H_h = 0.d0
    allocate(fluc(lat_size,lat_size,lat_size))
    en_addit = 0.d0

    if (random_seed_from_clock) then
        call init_random_seed()
    endif
    do i = 1, NN 
        do j = 1, NN
            do k = 1, NN
                CALL RANDOM_NUMBER(xx) ! at this stage xx is a random number from [0,1)
                xx = -amp + 2*amp*xx   ! at this stage xx is a random number from [-amp,amp)
                fluc(i,j,k) = xx
            enddo
        enddo    
    enddo

    !row interactions e/h
    do k = 1, NN    ! the third dimension
        do i = 1, NN ! x,y dimension
            do j = 1, NN-1
                hop_e( (k-1)*NN**2+(i-1)*NN+j+1 , (k-1)*NN**2+(i-1)*NN+j   ) = -t_e
                hop_e( (k-1)*NN**2+(i-1)*NN+j   , (k-1)*NN**2+(i-1)*NN+j+1 ) = -t_e
                hop_h( (k-1)*NN**2+(i-1)*NN+j+1 , (k-1)*NN**2+(i-1)*NN+j   ) = -t_h
                hop_h( (k-1)*NN**2+(i-1)*NN+j   , (k-1)*NN**2+(i-1)*NN+j+1 ) = -t_h
            enddo
        enddo
    enddo    
    !column interactions e/h
    do k = 1, NN    ! the third dimension
        do j = 1, NN
            do i = 1, NN-1
                hop_e( (k-1)*NN**2+(i-1)*NN+j , ((NN**2)*(k-1))+((i)*NN)+j ) = -t_e
                hop_e( ((NN**2)*(k-1))+((i)*NN)+j , (k-1)*NN**2+(i-1)*NN+j ) = -t_e
                hop_h( (k-1)*NN**2+(i-1)*NN+j , ((NN**2)*(k-1))+((i)*NN)+j ) = -t_h
                hop_h( ((NN**2)*(k-1))+((i)*NN)+j , (k-1)*NN**2+(i-1)*NN+j ) = -t_h
            enddo
        enddo
    enddo
    !vertical hopping interaction e/h
    do i = 1, NN    
        do j = 1, NN
            do k = 1, NN-1
                hop_e( ((NN**2)*(k-1))+(i)*NN+j , ((NN**2)*(k))+((i)*NN)+j ) = -t_e
                hop_e( ((NN**2)*(k))+((i)*NN)+j , ((NN**2)*(k-1))+(i)*NN+j ) = -t_e
                hop_h( ((NN**2)*(k-1))+(i)*NN+j , ((NN**2)*(k))+((i)*NN)+j ) = -t_h
                hop_h( ((NN**2)*(k))+((i)*NN)+j , ((NN**2)*(k-1))+(i)*NN+j ) = -t_h
            enddo
        enddo
    enddo
    ! onsite interaction e/h
    do i = 1, NN
        do j = 1, NN
            do k = 1, NN
                onsit_e( (k-1)*NN**2+(i-1)*NN+j, (k-1)*NN**2+(i-1)*NN+j ) = eps_0_e + fluc(i,j,k)
                onsit_h( (k-1)*NN**2+(i-1)*NN+j, (k-1)*NN**2+(i-1)*NN+j ) = eps_0_h + fluc(i,j,k)
            enddo
        enddo
    enddo
    
    H_e = onsit_e + hop_e
    H_h = onsit_h + hop_h

    allocate( eig_val_e(NN**3) ); allocate( eig_val_h(NN**3) )
    allocate( eig_vec_e(NN**3,NN**3) ); allocate( eig_vec_h(NN**3,NN**3) )
    eig_vec_e = H_e
    eig_vec_h = H_h
    
    ii=NN**3   ! Doubt
    call diasym(eig_vec_e,eig_val_e,ii)
    call diasym(eig_vec_h,eig_val_h,ii)
    Eg_0 = eig_val_e(el_index)-eig_val_h(hole_index)

    if (counter==1) then
        write(33,*) "Energy Gap Eg_0 = ", Eg_0
    endif
    !do i=1,NN**3
    !    write(26,*) eig_val_e(i),eig_val_h(i)
    !enddo
    call get_dos(NN,eig_val_e,eig_val_h,"0")

!    write(10,20)i,m1(:,i)
    !xx = 0
    do i = 1, NN
        do j = 1, NN
            do k = 1, NN
                n_e(i,j,k)=eig_vec_e( (k-1)*NN**2+(i-1)*NN+j , el_index )**2
                !n_e(i,j,k)=eig_vec_e( el_index , (k-1)*NN**2+(i-1)*NN+j )**2
                n_h(i,j,k)=eig_vec_h( (k-1)*NN**2+(i-1)*NN+j , hole_index )**2
                !xx=xx+n_e(i,j,k)
            enddo
        enddo
    enddo
    !write (*,*) "Eg_0 = ",Eg_0," norm = ",xx

    call print_density(n_e,n_h,NN,tag)
    conv_iteration = 0

    error_e = tolerance+1 ; error_h = tolerance+1
    do while (error_e > tolerance)  ! scf loop
        norm_e =0.d0; norm_h =0.d0
        do i = 1, NN
            do j = 1, NN
                do k = 1, NN
                    norm_e = norm_e + beta*n_e(i,j,k)+(1-beta)*n_e_prev(i,j,k)
                    norm_h = norm_h + beta*n_h(i,j,k)+(1-beta)*n_h_prev(i,j,k)
                enddo
            enddo
        enddo
        do i = 1, NN
            do j = 1, NN
                do k = 1, NN
                    pd_e(i,j,k) = (beta*n_e(i,j,k)+(1-beta)*n_e_prev(i,j,k))/norm_e
                    pd_h(i,j,k) = (beta*n_h(i,j,k)+(1-beta)*n_h_prev(i,j,k))/norm_h
                enddo
            enddo
        enddo
        do i = 1, NN
            do j = 1, NN
                do k = 1, NN
                    onsit_e( (k-1)*NN**2+(i-1)*NN+j, (k-1)*NN**2+(i-1)*NN+j ) = eps_0_e + fluc(i,j,k)
                    onsit_h( (k-1)*NN**2+(i-1)*NN+j, (k-1)*NN**2+(i-1)*NN+j ) = eps_0_h + fluc(i,j,k)
                enddo
            enddo
        enddo
        do i = 1, NN 
            do j = 1, NN
                do k = 1, NN
                    do l = 1, NN
                        do m = 1, NN
                            do n = 1, NN
                                onsit_e( (k-1)*NN**2+(i-1)*NN+j, (k-1)*NN**2+(i-1)*NN+j ) = & 
                                  onsit_e( (k-1)*NN**2+(i-1)*NN+j, (k-1)*NN**2+(i-1)*NN+j ) - &
                                  pd_h(l, m, n) * coulomb(i,j,k,l,m,n,u_0,alfa)
                                onsit_h( (k-1)*NN**2+(i-1)*NN+j, (k-1)*NN**2+(i-1)*NN+j ) = & 
                                  onsit_h( (k-1)*NN**2+(i-1)*NN+j, (k-1)*NN**2+(i-1)*NN+j ) + &
                                  pd_e(l, m, n) * coulomb(l,m,n,i,j,k,u_0,alfa)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
        H_e = onsit_e + hop_e
        H_h = onsit_h + hop_h

        eig_vec_e = H_e
        eig_vec_h = H_h

        ii=NN**3
        call diasym(eig_vec_e,eig_val_e,ii)
        call diasym(eig_vec_h,eig_val_h,ii)
        Eg_coulomb = eig_val_e(el_index)-eig_val_h(hole_index)

        xx = 0.d0
        n_e_prev(:,:,:) = n_e
        n_h_prev(:,:,:) = n_h
        do i = 1, NN
            do j = 1, NN
                do k = 1, NN
                    n_e(i,j,k)=eig_vec_e( (k-1)*NN**2+(i-1)*NN+j , el_index )**2
                    n_h(i,j,k)=eig_vec_h( (k-1)*NN**2+(i-1)*NN+j , hole_index )**2
                enddo
            enddo
        enddo
        ! calculate the scf error
        error_e = 0.d0 ; error_h = 0.d0
        do i = 1, NN
            do j = 1, NN
                do k = 1, NN
                    error_e = error_e + abs(n_e(i,j,k)-n_e_prev(i,j,k))
                    error_h = error_h + abs(n_h(i,j,k)-n_h_prev(i,j,k))
                enddo
            enddo
        enddo

        if (output_display) then
            !write(*,"(I3,3x,f6.4,1x,f7.5,1x,f7.5)") conv_iteration, Eg_coulomb, error_e, error_h
            write(*,*) conv_iteration, Eg_coulomb, error_e, error_h
        endif

        conv_iteration = conv_iteration + 1
        call get_dos(NN,eig_val_e,eig_val_h,"U")

    end do ! scf loop

    ! debug
    call print_density(n_e,n_h,NN,tag)

    ! calculate the additional energy term
    do i = 1, NN
        do j = 1, NN
            do k = 1, NN
                do l = 1, NN
                    do m = 1, NN
                        do n = 1, NN
                            en_addit = en_addit + &
                              coulomb(l,m,n,i,j,k,u_0,alfa)*n_e(i,j,k)*n_h(l,m,n)
                        enddo
                    enddo
                enddo
            enddo
        enddo    
    enddo        
    E_binding = Eg_0 - Eg_coulomb - en_addit
    e_bind_array(counter) = E_binding

    write(tag2,*) counter
    open(34,file='worklog'// trim(adjustl(tag2)) //'.dat')
    write(34,*) "#Iteration      E_binding"
    write(34,*) counter, E_binding
    close(34)

    deallocate(n_e); deallocate(n_h)
    deallocate(n_e_prev); deallocate(n_h_prev)
    deallocate(fluc)
    deallocate(pd_e); deallocate(pd_h)
    deallocate(H_e); deallocate(H_h)
    deallocate(onsit_e); deallocate(onsit_h)
    deallocate(hop_e); deallocate(hop_h)   
    deallocate(eig_val_e); deallocate(eig_val_h)   
    deallocate(eig_vec_e); deallocate(eig_vec_h)   
    counter = counter +1

end do ! statistical sample loop

do i = 1, 300
    bin(i) = -15+((0.1)*i)
enddo
!call Distribute(eig_val_e,eig_val_h,NN**3,bin,300)

call CPU_TIME( t2 )
e_bind_avg = 0.d0
e_bind_std = 0.d0
do i = 1, statistical_sample
    e_bind_avg = e_bind_avg + e_bind_array(i)
enddo
e_bind_avg = e_bind_avg/statistical_sample
do i = 1, statistical_sample
    e_bind_std = e_bind_std + (e_bind_array(i)-e_bind_avg)**2/statistical_sample
enddo
e_bind_std = (e_bind_std)**.5

write(33,*) "Binding_energy    ", e_bind_avg, &
            "Standard_deviation   ", e_bind_std
write(33,*) "Binding energies of the samples"
do i = 1, statistical_sample
    write(33,*) i, e_bind_array(i)
enddo
write(33,*) "Cputime = ", t2-t1
close(33)
deallocate(e_bind_array)


end program exciton3d
!-------------------!

!---------------------------------------------------------!
!Calls the LAPACK diagonalization subroutine DSYEV        !
!input:  a(n,n) = real symmetric matrix to be diagonalized!
!            n  = size of a                               !
!output: a(n,n) = orthonormal eigenvectors of a           !
!        eig(n) = eigenvalues of a in ascending order     !
!---------------------------------------------------------!
 subroutine diasym(a,eig,num)
 implicit none

 integer num,l,inf
 real*8  a(num,num),eig(num),work(num*(3+num/2))

 l=num*(3+num/2)
 call dsyev('V','U',num,a,num,eig,work,l,inf)

 end subroutine diasym
!---------------------!

subroutine print_density(n_e,n_h,N,tag)
    implicit none
    integer :: i, j, k, N
    !real(8), allocatable :: n_e(:,:,:), n_h(:,:,:)
    real(8) :: n_e(N,N,N), n_h(N,N,N)
    character(40) :: tag

    !allocate(n_e(N,N,N)); allocate(n_h(N,N,N))
    open(10,file='density.dat',status='replace')
    !open(10,file='density'// trim(adjustl(tag)) //'.dat',status='replace')
    do i = 1, N ! Doubt
        do j = 1, N
            do k = 1, N
                write(10,*) i,j,k,n_e(i,j,k),n_h(i,j,k)
            enddo
            write(10,*) ""
        enddo
    enddo
    close(10)
    end subroutine print_density

real(8) function coulomb(i,j,k,l,m,n,u_0,alfa)
    implicit none
    integer, intent(in) :: i,j,k,l,m,n
    real(8), intent(in) :: u_0,alfa
    
    !coulomb = (alfa*u_0)/( dsqrt( (i-k)**2 + (j-l)**2) )
    if (i==l .and. j==m .and. k==n) then
        coulomb = u_0
    else
        coulomb = (alfa*u_0)/( ( (i-l)**2 + (j-m)**2 + (k-n)**2) )**.5
    endif
end function coulomb

subroutine get_dos(NN,eig_val_e,eig_val_h,tag)
    implicit none
    integer, intent(in) :: NN
    real(8), intent(in) :: eig_val_e(NN**3),eig_val_h(NN**3)
    character(1), intent(in) :: tag
    integer :: i,i2,npt,neigenv,lor
    real(8) :: de,sigma,gaussian,lorentzian,constant,res,emin,emax
    real(8), allocatable :: en(:),dos(:),test(:)
    emin = -15 !eig_val_h(1)
    emax = 15  !eig_val_e(NN**3)
    neigenv = NN**3 
    npt = 1000 ! hard-coded: number of division of the energy (x) axis
    !sigma = 0.6*abs(eig_val_h(NN**3)-eig_val_h(NN**3-4))
    sigma = 0.25           ! hard-coded: broadening
    allocate(en(npt))
    allocate(dos(npt))
    de=(emax-emin)/(npt-1)
    do i=1,npt
        en(i)=emin+de*(i-1)
    enddo
    dos(:)=0.d0
    do i=1,neigenv
        do i2=1,npt
            dos(i2)=dos(i2)+lorentzian(en(i2),eig_val_h(i),sigma)
            dos(i2)=dos(i2)+lorentzian(en(i2),eig_val_e(i),sigma)
        enddo
    enddo
    !open(57,file="dos.dat")
    open(57,file='dos'//  trim(adjustl(tag)) //'.dat')
    write(57,*) "## emin,emax=", emin,emax
    write(57,*) "## npt,sigma=", npt,sigma
    do i2=1,npt-1
        write(57,*)  en(i2),dos(i2)
    enddo
    close(57)
    deallocate(en,dos)
end subroutine get_dos

function gaussian(x,x0,sigma)
 IMPLICIT NONE
 real (8), intent(in) :: x, x0, sigma
 real (8) :: k,gaussian
 gaussian=exp(-(x-x0)**2/(2*sigma**2))/(sigma*sqrt(2.d0*3.1415927))
end function gaussian

function lorentzian(x,x0,sigma)
 IMPLICIT NONE
 real (8), intent(in) :: x, x0, sigma
 real (8) :: k,lorentzian
 lorentzian=sigma/2/3.1415927/((x-x0)**2+(sigma/2)**2)
end function lorentzian

!!!
!!!
!!!  Pankaj subroutines
!!!
!!!
subroutine Distribute(X1,X2,states,rang,M)
    implicit none
    integer, intent(in) :: states                     ! # of scores
    real(8), intent(in) :: X1(states), X2(states)     ! input eigenvalues
    integer, intent(in) :: M                          ! # of ranges
    real, intent(in)    :: rang(M)                    ! range array
    integer             :: i, j  
    real                :: Bucket(M)                  ! counting bucket
    real                :: lorentzian(3,M), interpolation(2*M+1)
    real                :: gamma_lorenz
      
    gamma_lorenz = 0.1
      
    do i = 1, M                                       ! clear buckets
        Bucket(i) = 0
    enddo
    do i = 1, M                   
        do j = 1, 3
            lorentzian(j,i) = 0
        enddo
    enddo
    do i = 1, (2*M)+1
        interpolation(i) = 0
    enddo
    do i = 1, states                                  ! for each input score
        do j = 1, M                                   ! determine the bucket
            if (X1(i) < rang(j)) then
                Bucket(j) = Bucket(j) + 1
                exit
            END IF               
        enddo                       
    enddo
    write(*,*) 'DEBUG'
    do i = 1, states                                  ! for each input score
        do j = 1, M                                   ! determine the bucket
            if (X2(i) < rang(j)) then
                Bucket(j) = Bucket(j) + 1
                exit
            END IF               
        enddo                       
    enddo
    !!!!!!Lorentzian:
    do j = 1, M
        do i = 1, 3
            lorentzian(i,j) = (Bucket(j)*((gamma_lorenz/2)**2))/((((i-2)*(gamma_lorenz/2))**2)+((gamma_lorenz/2)**2))
        enddo
    enddo    
    interpolation(1) = lorentzian(1,1)
    do i = 1, M-1
        interpolation((2*i)+1) = (lorentzian(1,i+1)+lorentzian(3,i-1))/2
    enddo
    interpolation(2*M+1) = lorentzian(3,M)
    do i = 1, M
        interpolation(2*i) = (interpolation((2*i)-1)+interpolation((2*i)+1))/2
    enddo
    call Dos_file(interpolation, (2*M)+1)           ! print a histogram
end subroutine Distribute

!END MODULE  Density_of_state

subroutine Dos_file(coun, div)
    implicit none
    integer :: i, div
    real :: coun(div)
    open(10,file='Dos_file.dat',status='replace')
    do i = 1, div 
        write(10,*) -15+(0.05*(i-1)), "    ", coun(i)
    enddo
    close(10)
end subroutine Dos_file 

SUBROUTINE init_random_seed()
INTEGER :: i, n, clock
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
          
CALL SYSTEM_CLOCK(COUNT=clock)
          
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)
          
DEALLOCATE(seed)
END SUBROUTINE

