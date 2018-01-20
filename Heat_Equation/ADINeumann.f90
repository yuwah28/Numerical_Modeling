     Module parameters
       implicit none
       INTEGER, PARAMETER :: Nx = 100, Ny = 100 , timemax = 200 , t_steps = 1
       REAL, PARAMETER :: xa = 0.d0 ,  xb = 6.d0 , ya = 0.d0 , yb = 6.d0 !'domain'
       REAL, PARAMETER :: dx = (xb-xa)/Nx , dy = (yb-ya)/Ny , dt = 0.01 , h = dx**2.d0
       REAL, PARAMETER :: rat = dt/(2.d0*dx**2.d0)   !'physical paramters'
       REAL, PARAMETER :: pi = 4.d0*atan(1.d0)
        !==For computing exact solution
       INTEGER, PARAMETER :: Mi = 30 , Nj = 30
       INTEGER, PARAMETER :: x0 = xb/2.d0 , y0 = yb/2.d0
       REAL, PARAMETER :: Amp = 2.d0
     end module parameters



     Program main
       use parameters
       implicit none

       INTEGER i , j , t , k , dim , frame  
       INTEGER, PARAMETER :: dimen = Nx+Ny+2
       REAL(8) :: time0 , time , ssum
       REAL(8), ALLOCATABLE ::  x(:) , y(:)  
       REAL(8) xmesh(Nx+1,Ny+1) , ymesh(Nx+1,Ny+1)      
       REAL(8), ALLOCATABLE :: Csym(:,:) , BC1(:,:) , BC2(:,:)
       REAL(8), ALLOCATABLE :: uold(:,:) , uhat(:,:), udel(:,:) , Bn_prime(:,:) , Bn(:,:) , u(:,:)
       REAL(8), ALLOCATABLE :: a(:) , b(:) , c(:) , r(:) , v(:) !'tridag entries'
       !==For computing exact solution
       INTEGER m , n
       REAL(8), ALLOCATABLE :: uint(:,:) !, QSumm(:) , QSumn(:)
       REAL(8) :: QSumx , QSumy
       REAL(8), ALLOCATABLE :: ufun(:,:) , uexact(:,:)
       REAL(8) :: error , errorfun ,  maxnorm , energy

       allocate(x(1:Nx+1),y(1:Ny+1))
       allocate(Csym(1:Ny+1,1:Ny+1),BC1(1:Nx+1,1:Ny+1),BC2(1:Ny+1,1:Nx+1))
       allocate(uold(1:Nx+1,1:Ny+1),uint(1:Nx+1,1:Ny+1),uhat(1:Nx+1,1:Ny+1),udel(1:Nx+1,1:Ny+1))
       allocate(Bn_prime(1:Nx+1,1:Ny+1),Bn(1:Nx+1,1:Ny+1),u(1:Nx+1,1:Ny+1))
       allocate(a(dimen) , b(dimen) , c(dimen) , r(dimen) , v(dimen))
       !==For computing exact solution
       !allocate(QSumm(1:Mi),QSumn(1:Nj))
       allocate(ufun(1:Nx+1,1:Ny+1),uexact(1:Nx+1,1:Ny+1)) !'initial vlaue for exact soln'

       !meshing
       do j = 1 , Ny+1
          do i = 1 , Nx+1
             x(i) = xa + real(i-1)*dx
             y(j) = ya + real(j-1)*dy
             xmesh(i,j) = x(i)
             ymesh(i,j) = y(j)
          enddo
       enddo
   
       
       !==Computing Fourier coefficients for exact soln==
       !call simpson(x,y,QSumm,QSumn)


       !Initial Condition
       !==Open File==
       open(unit=100, file='initialN.dat')
       ! Header for contour.dat (fileformat tecplot)
       write(100,*) 'VARIABLES = "X", "Y", "U"'
       write(100,*) 'Zone F=point, I =', Nx+1, 'J =',Ny+1
       !==========================================================
       do i = 1 , Nx+1
          do j = 1 , Ny+1
             uold(i,j) = Amp*exp(-0.5*((x(i)-x0)**2.d0 + (y(j)-y0)**2.d0))
             uint(i,j) = uold(i,j)
             uexact(i,j) = uold(i,j)
             u(i,j) = uold(i,j)
          enddo
       enddo

       do j = 1 , Ny+1
          do i = 1 , Nx+1
             write(100,*) xmesh(i,j),ymesh(i,j),uold(i,j)
          enddo
       enddo
       !==End I.C==

       !==Euler Step==
       ! do i = 2 , Nx
       !     do j = 2 , Ny
       !         u(i,j) = uold(i,j) + (dt/h) * ( uold(i-1,j) -2.d0*uold(i,j) +uold(i+1,j) + &
       !                 uold(i,j-1) -2.d0*uold(i,j) + uold(i,j+1) )
       !     enddo
       ! enddo
        !==End Euler==
       


        !==Open File==
        !open(unit=200, file='heatout.dat')
        !! Header for contour.dat (fileformat tecplot)
        !write(200,*) 'VARIABLES = "X", "Y", "U", "W"'
        !write(200,*) 'Zone F=point, I =', Nx+1, 'J =',Ny+1

        !==Error Norm==
        open(unit=300, file='normN.dat')

        !==Simulation/Movie:
        time0 = 0.d0
        time = 0.d0
        frame = 0
        call output_grid(frame,time,x,y,u,uexact)
        print "(a,i5,a,i5,a,e16.8)","Writing frame ",frame," during step tn=",0," t=",time


        !Time-evolution
        !==========================================================
        do t = 1 ,timemax


           time = time0 + real(t)*dt

           !Boundary Conditions
           !==========================================================
           !'Define Ghost Cells for Neumann'
            ! do j = 1 , Ny+1
            !    uold(1,j) = uold(2,j)
            !    uold(Nx+1,j) = uold(Nx,j)
            !enddo
            !do i = 1 , Nx+1
            !    uold(i,1) = uold(i,2)
            !    uold(i,Ny+1) = uold(i,Ny)
            !enddo
            !==End Ghost Cells==


            !==Start ADI 1st step==
            !for some fixed j
            do j = 2, Ny

                Bn_prime(1,j) = (1-4.d0*rat+4.d0*rat**2.d0)*uold(1,j) + &
                    (rat-2.d0*rat**2.d0) * uold(2,j) + (rat-2.d0*rat**2.d0) * uold(2,j) + &
                    (rat-2.d0*rat**2.d0) * uold(1,j-1) + (rat-2.d0*rat**2.d0) * uold(1,j+1) + &
                    (rat**2.d0)*(uold(2,j-1) +uold(2,j-1) +uold(2,j+1) +uold(2,j+1) )

                Bn_prime(Nx+1,j) = (1-4.d0*rat+4.d0*rat**2.d0)*uold(Nx+1,j) + &
                    (rat-2.d0*rat**2.d0) * uold(Nx,j) + (rat-2.d0*rat**2.d0) * uold(Nx,j) + &
                    (rat-2.d0*rat**2.d0) * uold(Nx+1,j-1) + (rat-2.d0*rat**2.d0) * uold(Nx+1,j+1) + &
                    (rat**2.d0)*(uold(Nx,j-1) +uold(Nx,j-1) +uold(Nx,j+1) +uold(Nx,j+1) )

                Bn_prime(2,j) = (1-4.d0*rat+4.d0*rat**2.d0)*uold(2,j) + &
                    (rat-2.d0*rat**2.d0) * uold(1,j) + (rat-2.d0*rat**2.d0) * uold(3,j) + &
                    (rat-2.d0*rat**2.d0) * uold(2,j-1) + (rat-2.d0*rat**2.d0) * uold(2,j+1) + &
                    (rat**2.d0)*(uold(1,j-1) +uold(3,j-1) +uold(1,j+1) +uold(3,j+1) )

                Bn_prime(Nx,j) = (1-4.d0*rat+4.d0*rat**2.d0)*uold(Nx,j) + &
                    (rat-2.d0*rat**2.d0) * uold(Nx-1,j) + (rat-2.d0*rat**2.d0) * uold(Nx+1,j) + &
                    (rat-2.d0*rat**2.d0) * uold(Nx,j-1) + (rat-2.d0*rat**2.d0) * uold(Nx,j+1) + &
                    (rat**2.d0)*(uold(Nx-1,j-1) +uold(Nx+1,j-1) +uold(Nx-1,j+1) +uold(Nx+1,j+1) )

                !'Bn(i,j) given j for i = 2,Nx'

                do i  =  3, Nx-1
                    Bn(i,j) = (1-4.d0*rat+4.d0*rat**2.d0)*uold(i,j) + &
                            (rat-2.d0*rat**2.d0) * uold(i-1,j) + (rat-2.d0*rat**2.d0) * uold(i+1,j) + &
                            (rat-2.d0*rat**2.d0) * uold(i,j-1) + (rat-2.d0*rat**2.d0) * uold(i,j+1) + &
                            (rat**2.d0)*(uold(i-1,j-1) +uold(i+1,j-1) +uold(i-1,j+1) +uold(i+1,j+1) )
                enddo
                    Bn(2,j) = Bn_prime(2,j) + (rat/(1+2.d0*rat))*Bn_prime(1,j)

                    Bn(Nx,j) = Bn_prime(Nx,j) + (rat/(1+2.d0*rat))*Bn_prime(Nx+1,j)


                !'Compute for Sym A (Nx-1 x Nx -1) for Direchlet'
                !'Compute for Sym A (Nx+1 x Nx +1) for Neumann'
                dim = Nx-1

                do i = 1 , Nx-1
                    a(i) = -rat
                    b(i) = (1.d0+2.d0*rat)
                    c(i) = -rat
                enddo
                b(1) = (-rat**2.d0)/(1+2.d0*rat) + (1+2.d0*rat)
                b(Nx-1) = (-rat**2.d0)/(1+2.d0*rat) + (1+2.d0*rat)

                do i = 1 , Nx-1
                    r(i) = Bn(i+1,j)
                enddo

                call tridag(a,b,c,r,v,dim)

                do i = 1 , Nx-1
                    uhat(i+1,j) = v(i)
                    !print*, 'uhat', t,i,j,uhat(i,j)
                enddo

                 !print*, 'uhat', i,j,uhat(i,j)
          enddo
          !==End 1st step==


          !==Start ADI 2nd step==
          !for some fixed i
          do i = 2, Nx

             !'Compute for Sym C (Ny-1 x Ny -1)'
             dim = Ny-1
             !=='Dirac'==
             !do j = 1 , Ny-1
             !    a(j) = -rat
             !    b(j) = (1.d0+2.d0*rat)
             !    c(j) = -rat
             !enddo
             !==end dirac==

             !=='Neumann'==
             do j = 1 , Ny-1
                a(j) = -rat
                b(j) = (1.d0+2.d0*rat)
                c(j) = -rat
             enddo
             b(1) = 1.d0 + rat
             b(Ny-1) = 1.d0 + rat
             !==end neumann==

             do j = 1 , Ny-1
                r(j) = uhat(i,j+1)
             enddo

             call tridag(a,b,c,r,v,dim)

             do j = 1 , Ny-1
                u(i,j+1) = v(j)
             enddo

         enddo

         !==Augmented Neumann==
         do i = 1 , Nx+1
            u(i,1) = u(i,2)
            u(i,Ny+1) = u(i,Ny)
         enddo

         do j = 1 ,Ny+1
            u(1,j) = u(2,j)
            u(Nx+1,j) = u(Nx,j)
         enddo
         !==end augmentation
         !==End 2st step==

         !==Updata u(i,j)==
         uold = u

         !==Construct Exact equation==
         ssum = 0.d0

         do i = 1 , Nx+1
             do j = 1 , Ny+1

                uexact(i,j) = QSumx(x(i),time)*QSumy(y(j),time)*Amp/(4.d0*pi*time)

             enddo
         enddo

         !==end construction==


         !do i = 1 , Nx+1
         !   do j = 1 , Ny+1
         !       print*, i,j,x(i),y(j),u(i,j),uexact(i,j),abs(u(i,j) - uexact(i,j))
         !   enddo
         !enddo
         !==Error Norm==
         maxnorm = 0.d0
         do j = 1 , Ny+1
            do i = 1 , Nx+1
                error = abs(u(i,j) - uexact(i,j))
                if(error.gt.maxnorm) then
                maxnorm = error
                endif
                   !print*, i , j, maxnorm
            enddo
         enddo
     
         !==Check Energy==
         call energy_conserv(uexact,energy)
         !==End energy check==


         write(300,*) t , maxnorm , energy


         !==Write out u every t_steps==
         if (mod(t,t_steps) == 0) then
            frame = frame + 1
            call output_grid(frame,time,x,y,u,uexact)
            print "(a,i5,a,i5,a,e16.8)","Writing frame ",frame," during step tn=",t," t=",time
         endif


         !==Record Output Data at Final Time==
         !if(t.eq.20) then
         !    do j = 1 , Ny+1
         !        do i = 1 , Nx+1
         !            write(200,*) xmesh(i,j),ymesh(i,j),u(i,j),uexact(i,j)
         !        enddo
         !    enddo
         !endif


      enddo
      !==End time-evolution==

        !close(200)

    End Program




    !===Subroutine-dance===
    SUBROUTINE tridag(a,b,c,r,v,dim)
      INTEGER dim, NMAX
      REAL(8), intent(in) :: a(dim) , b(dim) , c(dim) , r(dim)
      REAL(8), intent(out) ::  v(dim)
      PARAMETER (NMAX=5000)
      !Solves for a vector solution v(1:dim) of length n the tridiagonal linear set given by eqn(2.4.1)
      !a(1:dim) , b(1:dim) , c(1:dim) , and r(1:dim) are input vectors and are not modified.
      !Parameter: NMAX is the maximum expected vallues of n.
      !======================================================
      ! [ b c 0 0 0][v1]   [r1]
      ! [ a b c 0 0][v2]   [r2]
      ! [ 0 a b c 0][v3] = [r3]   eqn(2.4.1)
      ! [ 0 0 a b c][v4]   [r4]
      ! [ 0 0 0 a b][v5]   [r5]
      !======================================================
      INTEGER j
      REAL bet , gam(NMAX)

      if(b(1).eq.0.) stop   'tridag:rewrite equation'
      ! if this happens then you should rewrite your equations as a set of order N-1 with u2 trivially eliminated
      bet = b(1)
      v(1) = r(1)/bet

      do j = 2 , dim
         gam(j) = c(j-1)/bet
         bet = b(j) - a(j)*gam(j)
         if(bet.eq.0.)stop  'tridag failed'
         v(j) = (r(j)-a(j)*v(j-1))/bet
      enddo

      do j = dim-1,1,-1
         v(j) = v(j) - gam(j+1)*v(j+1)
      enddo

      return
    End SUBROUTINE tridag




    !===Function for Exact Soln===
    FUNCTION Qsumx(xin,tau)
      use parameters
      implicit none

      INTEGER i
      REAL(8) :: Qsumx
      REAL(8), intent(in) :: xin , tau
      REAL(8) :: sm2, sm4 
      REAL(8) :: w(1:Nx+1)

      do i = 1 , Nx+1
         w(i) = xa + real(i-1)*dx
      enddo

      !===Sum_x phi_m===!
      !--compute s2
         sm2 = 0.d0

         do i=1 , Nx/2-1
            !~: s2 = s2 + fun(2*i)
            !~: Here fun(i) := (exp( -((xin - w(i))**2.d0)/(4.d0*tau)) + exp( -((xin + w(i))**2.d0)/(4.d0*tau)))*exp(-0.5*(w(i)-x0)**2.d0)
            sm2 = sm2 + (exp( -((xin - w(2*i))**2.d0)/(4.d0*tau)) + &
                    exp( -((xin + w(2*i))**2.d0)/(4.d0*tau)))*exp(-0.5*(w(2*i)-x0)**2.d0)
         enddo

         !--compute s4
         sm4 = 0.d0

         do i=1, Nx/2
            !~: s4 = s4 + fun(2*i-1)
            sm4 = sm4 + (exp( -((xin - w(2*i-1))**2.d0)/(4.d0*tau)) + &
                    exp( -((xin + w(2*i-1))**2.d0)/(4.d0*tau)))*exp(-0.5*(w(2*i-1)-x0)**2.d0)
         enddo

         !--compute total sum
         !~: dx/3.d0 * !(fun(1) + 2.d0*s2 + 4.d0*s4 + fun(Nx+1))
         QSumx = (dx/3.d0) * ((exp( -((xin - w(1))**2.d0)/(4.d0*tau)) + &
                    exp( -((xin + w(1))**2.d0)/(4.d0*tau)))*exp(-0.5*(w(1)-x0)**2.d0) + &
                    2.d0*sm2 + 4.d0*sm4 + &
                    (exp( -((xin - w(Nx+1))**2.d0)/(4.d0*tau)) + &
                    exp( -((xin + w(Nx+1))**2.d0)/(4.d0*tau)))*exp(-0.5*(w(Nx+1)-x0)**2.d0))

    End FUNCTION Qsumx

    !===Function for Exact Soln===
    FUNCTION QSumy(yin,tau)
        use parameters
        implicit none

        INTEGER i
        REAL(8) :: Qsumy
        REAL(8), intent(in) :: yin , tau
        REAL(8) :: sn2 , sn4
        REAL(8) :: z(1:Ny+1)

        do i = 1 , Ny+1
           z(i) = ya + real(i-1)*dy
        enddo

        !===Sum_x phi_m===!
        !--compute s2
        sn2 = 0.d0

        do i=1 , Ny/2-1
            !~: s2 = s2 + fun(2*i)
            !~: Here fun(i) := (exp( -((yin - z(i))**2.d0)/(4.d0*tau)) + exp( -((yin + z(i))**2.d0)/(4.d0*tau)))*exp(-0.5*(z(i)-y0)**2.d0)
            sn2 = sn2 + (exp( -((yin - z(2*i))**2.d0)/(4.d0*tau)) + &
                    exp( -((yin + z(2*i))**2.d0)/(4.d0*tau)))*exp(-0.5*(z(2*i)-y0)**2.d0)
        enddo

        !--compute s4
        sn4 = 0.d0

        do i=1, Ny/2
            !~: s4 = s4 + fun(2*i-1)
            sn4 = sn4 + (exp( -((yin - z(2*i-1))**2.d0)/(4.d0*tau)) + &
                    exp( -((yin + z(2*i-1))**2.d0)/(4.d0*tau)))*exp(-0.5*(z(2*i-1)-y0)**2.d0)
        enddo

        !--compute total sum
        !~: dx/3.d0 * !(fun(1) + 2.d0*s2 + 4.d0*s4 + fun(Nx+1))
        QSumy = (dy/3.d0) * ((exp( -((yin - z(1))**2.d0)/(4.d0*tau)) + &
                exp( -((yin + z(1))**2.d0)/(4.d0*tau)))*exp(-0.5*(z(1)-y0)**2.d0) + &
                2.d0*sn2 + 4.d0*sn4 + &
                (exp( -((yin - z(Ny+1))**2.d0)/(4.d0*tau)) + &
                exp( -((yin + z(Ny+1))**2.d0)/(4.d0*tau)))*exp(-0.5*(z(Ny+1)-y0)**2.d0))

    End FUNCTION QSumy


    !===Energy Conservation==
    SUBROUTINE energy_conserv(u,energy)
        use parameters
        implicit none

        INTEGER i , j
        REAL(8), intent(in) :: u(1:Nx+1,1:Ny+1)
        REAL(8) :: em2 , em4
        REAL(8) :: energy_x(1:Ny+1)
        REAL(8) :: en2 , en4
        REAL(8), intent(out) ::  energy


        !===Sum_x phi_m===!
        !--compute s2
        do j = 1 , Ny+1

            em2 = 0.d0

            do i=1 , Nx/2-1
                !~: s2 = s2 + fun(2*i)
                !~: Here fun(i) := u(i,j)
                em2 = em2 + u(2*i,j)
            enddo

            !--compute s4
            em4 = 0.d0

            do i=1, Nx/2
                !~: s4 = s4 + fun(2*i-1)
                em4 = em4 + u(2*i-1,j)
            enddo

            !--compute total sum
            !~: dx/3.d0 * !(fun(1) + 2.d0*s2 + 4.d0*s4 + fun(Nx+1))
            energy_x(j) = (dx/3.d0) * ( u(1,j) + 2.d0*em2 + 4.d0*em4 + u(Nx+1,j) )

        enddo


        en2 = 0.d0

        do i=1 , Ny/2-1
            !~: s2 = s2 + fun(2*i)
            !~: Here fun(i) := energy_x(i)
            en2 = en2 + energy_x(2*j)
        enddo

        !--compute s4
        en4 = 0.d0

        do i=1, Ny/2
            !~: s4 = s4 + fun(2*i-1)
            en4 = en4 + energy_x(2*j-1)
        enddo

        !--compute total sum
        !~: dx/3.d0 * !(fun(1) + 2.d0*s2 + 4.d0*s4 + fun(Nx+1))
        energy = (dy/3.d0) * ( energy_x(1) + 2.d0*en2 + 4.d0*en4 + energy_x(Ny+1) )

    end subroutine energy_conserv
    !===End energy==


    !===Subroutine-dance===
    SUBROUTINE output_grid(frame,time,x,y,u,uexact)
      use parameters
      implicit none

      ! Arguments
      integer, intent(in) :: frame
      double precision, intent(in) :: time
      double precision :: x(1:Nx+1) , y(1:Ny+1)
      double precision, dimension(1:Nx+1,1:Ny+1) :: u , uexact

      ! Locals
      integer :: i,j


       ! Open output file and write file header:

       if (time==0) then
          open(unit=70,file='outputN.dat',access='sequential',status='unknown')
          write(70,*) ' VARIABLES= "x", "y", "u"'
          write(70,101) time,Nx+1,Ny+1
       endif

101    FORMAT('ZONE T="t = ',e26.16,'"',' F=POINT, I=',I5,' J=', I5)

       if (time>0) then
          write(70,101) time,Nx+1,Ny+1
       endif

       ! Write out data
       do j=1,Ny+1
          do i=1,Nx+1

             ! Reset to zero if exponent is too large
             if (abs(u(i,j)) < 1d-99) u(i,j) = 0.d0
             if (abs(uexact(i,j)) < 1d-99) uexact(i,j) = 0.d0

             write(70,"(5e26.16)") x(i),y(j),u(i,j)

          enddo
       enddo

     end subroutine output_grid
