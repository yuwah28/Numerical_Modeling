     Module parameters
       implicit none
       INTEGER, PARAMETER :: Nx = 150, Ny = 150 , timemax =2000 , t_steps = 1
       REAL, PARAMETER :: xa = 0.d0 ,  xb = 6.d0 , ya = 0.d0 , yb = 6.d0 !'domain'
       REAL, PARAMETER :: dx = (xb-xa)/Nx , dy = (yb-ya)/Ny , dt = 0.0001
       REAL, PARAMETER :: rat = dt/(2.d0*dx**2.d0)   !'physical paramters'
       REAL, PARAMETER :: pi = 4.d0*atan(1.d0)
        !==For computing exact solution
       INTEGER, PARAMETER :: Mi = 20 , Nj = 20
       INTEGER, PARAMETER :: x0 = xb/2.d0 , y0 = yb/2.d0
       REAL, PARAMETER :: Amp = 2.d0 ,sigma = 0.01
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
       REAL(8), ALLOCATABLE :: uold(:,:) , uhat(:,:), udel(:,:) , Bn(:,:) , u(:,:)
       REAL(8), ALLOCATABLE :: a(:) , b(:) , c(:) , r(:) , v(:) !'tridag entries'
       !==For computing exact solution
       INTEGER m , n
       REAL(8), ALLOCATABLE :: uint(:,:) , QSumm(:) , QSumn(:)
       REAL(8), ALLOCATABLE :: ufun(:,:) , uexact(:,:) 
       REAL(8) :: error , errorfun ,  maxnorm , maxnormfun , norm , norm2

       allocate(x(1:Nx+1),y(1:Ny+1))
       allocate(Csym(1:Ny+1,1:Ny+1),BC1(1:Nx+1,1:Ny+1),BC2(1:Ny+1,1:Nx+1))
       allocate(uold(1:Nx+1,1:Ny+1),uint(1:Nx+1,1:Ny+1),uhat(1:Nx+1,1:Ny+1),udel(1:Nx+1,1:Ny+1))
       allocate(Bn(1:Nx+1,1:Ny+1),u(1:Nx+1,1:Ny+1))
       allocate(a(dimen) , b(dimen) , c(dimen) , r(dimen) , v(dimen))
       !==For computing exact solution
       allocate(QSumm(1:Mi),QSumn(1:Nj))
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
       call simpson(x,y,QSumm,QSumn)


 
       !Initial Condition
       !==Open File==
       open(unit=100, file='initialY.dat')
       ! Header for contour.dat (fileformat tecplot)
       write(100,*) 'VARIABLES = "X", "Y", "U"'
       write(100,*) 'Zone F=point, I =', Nx+1, 'J =',Ny+1
       !==========================================================
       do i = 1 , Nx+1
          do j = 1 , Ny+1
             uold(i,j) = Amp*exp(-sigma*((x(i)-x0)**2.d0 + (y(j)-y0)**2.d0))
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



       !Boundary Conditions
       !==========================================================
       !'actual solution'
       do i = 1 , Nx+1
          u(i,1) = 0.d0
          u(i,Ny+1) = 0.d0
       enddo

       do j = 1 ,Ny+1
          u(1,j) = 0.d0
          u(Nx+1,j) = 0.d0
       enddo

       !B.C for uhat'
       !'Matrix C'
       Csym = 0.d0
       do i = 1 , Ny+1
          Csym(i,i) = 1.d0+2.d0*rat
       enddo
       do i = 1 , Ny
          Csym(i,i+1) = -rat
       enddo
       do i = 2, Ny+1
          Csym(i,i-1) = -rat
       enddo
       
       !'left'
       do j = 1 , Ny+1
          uhat(1,j) = 0.d0
       enddo

       do j = 1 , Ny+1
          !'Martix Multiplication for'
          do k = 1 , Ny+1
             uhat(1,j) = uhat(1,j) + Csym(j,k)*u(1,k)
          enddo
       enddo

       !'right'
       do j = 1 , Ny+1
          uhat(Nx+1,j) = 0.d0
       enddo

       do j = 1 , Ny+1
          do k = 1 , Ny+1
             uhat(Nx+1,j) = uhat(Nx+1,j) + Csym(j,k)*u(Nx+1,k)
          enddo
       enddo


       !'BC1'
       do i = 3 , Nx-1
          do j = 1 , Ny+1
             BC1(i,j) = 0.d0
          enddo
       enddo
       do j = 1 , Ny+1
          BC1(2,j) = rat*uhat(1,j)
          BC1(Nx,j) = rat*uhat(Nx+1,j)
       enddo

       !'BC2'
       do j = 3 , Ny-1
          do i = 1 , Nx+1
             BC2(i,j) = 0.d0
          enddo
       enddo
       do i = 1 , Nx+1
          BC2(i,2) = rat*u(i,1)
          BC2(i,Ny) = rat*u(i,Ny+1)
       enddo
       !==End B.Cs==



       !==Open File==
       !open(unit=200, file='heatout.dat')
       !! Header for contour.dat (fileformat tecplot)
       !write(200,*) 'VARIABLES = "X", "Y", "U", "W"'
       !write(200,*) 'Zone F=point, I =', Nx+1, 'J =',Ny+1

      

       !==Simulation/Movie:
       time0 = 0.d0
       time = 0.d0
       frame = 0
       call output_grid(frame,time,x,y,u,uexact)
       print "(a,i5,a,i5,a,e16.8)","Writing frame ",frame," during step tn=",0," t=",time

       !==Error Norm==
       open(unit=300, file='normY.dat')

       if (time .eq. 0.d0) then

            !==Construct Exact equation==
            ssum = 0.d0

            do i = 1 , Nx+1
                do j = 1 , Ny+1

                    ssum = 0.d0
                    do m = 1 , Mi
                        do n = 1 , Nj
                            ssum = ssum +  QSumm(m)*QSumn(n) *sin(real(m)*pi*x(i)/xb)*sin(real(n)*pi*y(j)/yb)* &
                                exp(-time*(pi**2.d0)*((real(m)/xb)**2.d0+(real(n)/yb)**2.d0))
                        enddo
                    enddo
                    uexact(i,j) = Amp*ssum * 4.d0/(xb*yb)
                enddo
            enddo
           !== write initial error
            norm = 0.d0

           !==(Error) Max-Norm==
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

            !==2norm
            !do j = 1 , Ny+1
            !    do i = 1 , Nx+1
            !        norm = norm + abs(uint(i,j) - uexact(i,j))**2.d0
            !    enddo
            !enddo

            !norm2 = sqrt(norm)


            write(300,*) t , maxnorm !, maxnormfun

       endif


       !Time-evolution
       !==========================================================
       do t = 1 ,timemax

          time = time0 + real(t)*dt

          !==Start ADI 1st step==
          !for some fixed j
          do j = 2, Ny

             !'Bn(i,j) given j for i = 2,Nx'
             do i  = 2 , Nx
                Bn(i,j) = (1-4.d0*rat+4.d0*rat**2.d0)*uold(i,j) + &
                            (rat-2.d0*rat**2.d0) * uold(i-1,j) + (rat-2.d0*rat**2.d0) * uold(i+1,j) + &
                            (rat-2.d0*rat**2.d0) * uold(i,j-1) + (rat-2.d0*rat**2.d0) * uold(i,j+1) + &
                            (rat**2.d0)*(uold(i-1,j-1) +uold(i+1,j-1) +uold(i-1,j+1) +uold(i+1,j+1) )
             enddo



             !'Compute for Sym A (Nx-1 x Nx -1)'

             dim = Nx-1

             do i = 1 , Nx-1
                a(i) = -rat
                b(i) = (1.d0+2.d0*rat)
                c(i) = -rat
             enddo

             do i = 1 , Nx-1
                r(i) = (Bn(i+1,j) + BC1(i+1,j))
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
             do j = 1 , Ny-1
                a(j) = -rat
                b(j) = (1.d0+2.d0*rat)
                c(j) = -rat
            enddo

            do j = 1 , Ny-1
                r(j) = (uhat(i,j+1) + BC2(i,j+1))
            enddo


            call tridag(a,b,c,r,v,dim)

            do j = 1 , Ny-1
               u(i,j+1) = v(j)
            enddo

         enddo
         !==End 2st step==


         !==Updata u(i,j)==
         uold = u


         !==Construct Exact equation==
         ssum = 0.d0

         do i = 1 , Nx+1
             do j = 1 , Ny+1

             ssum = 0.d0
                 do m = 1 , Mi
                     do n = 1 , Nj
                         ssum = ssum +  QSumm(m)*QSumn(n) *sin(real(m)*pi*x(i)/xb)*sin(real(n)*pi*y(j)/yb)* &
                             exp(-time*(pi**2.d0)*((real(m)/xb)**2.d0+(real(n)/yb)**2.d0))
                     enddo
                 enddo

             uexact(i,j) = Amp*ssum * 4.d0/(xb*yb)
             enddo
         enddo

         !==end construction==

         !==Construct fundamental solution==
         do i = 1 , Nx+1
            do j = 1 , Ny+1
                ufun(i,j) = ( (1.d0/sqrt(2.d0*pi*time))**2 )*exp(-(x(i)-x0)**2/(2.d0*time))*exp(-(y(j)-y0)**2/(2.d0*time))
            enddo
         enddo
         !==end==

         !do i = 1 , Nx+1
         !   do j = 1 , Ny+1
         !       print*, i,j,x(i),y(j),u(i,j),uexact(i,j),abs(u(i,j) - uexact(i,j))
         !   enddo
         !enddo
         !==(Error) Max-Norm==
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
     
        !==(Error) 2-Norm==
         !norm = 0.d0
         !do j = 1 , Ny+1
         !   do i = 1 , Nx+1
         !       norm = norm + abs(u(i,j) - uexact(i,j))**2.d0
         !   enddo
         !enddo

         !norm2 = sqrt(norm)


         write(300,*) t , maxnorm !, maxnormfun

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




    !===Subroutine-dance===
    SUBROUTINE simpson(x,y,QSumm,QSumn)
      use parameters
      implicit none

      INTEGER i , j , m , n
      REAL(8), intent(in) :: x(1:Nx+1) , y(1:Ny+1)
      REAL(8) :: sm2, sm4 , sn2 , sn4
      REAL(8), intent(out) ::  QSumm(1:Mi) , QSumn(1:Nj)

    
      !===Sum_x phi_m===!
      !--compute s2
      do m = 1 , Mi

         sm2 = 0.d0

         do i=1 , Nx/2-1
            !~: s2 = s2 + fun(2*i)
            !~: Here fun(i) := exp( -sigma*(x(i)-xb)**2)*sin(m*pi*x(i)/xb)
            sm2 = sm2 + exp( -sigma*(x(2*i)-x0)**2.d0) * sin(real(m)*pi*x(2*i)/xb)
         enddo

         !--compute s4
         sm4 = 0.d0

         do i=1, Nx/2
            !~: s4 = s4 + fun(2*i-1)
            sm4 = sm4 + exp( -sigma*(x(2*i-1)-x0)**2.d0) * sin(real(m)*pi*x(2*i-1)/xb)
         enddo

         !--compute total sum
         !~: dx/3.d0 * !(fun(1) + 2.d0*s2 + 4.d0*s4 + fun(Nx+1))
         QSumm(m) = (dx/3.d0) * (exp( -sigma*(x(1)-x0)**2.d0) * sin(real(m)*pi*x(1)/xb) + 2.d0*sm2 + 4.d0*sm4 + &
              exp( -sigma*(x(Nx+1)-x0)**2.d0) * sin(real(m)*pi*x(Nx+1)/xb) )

      enddo



      !===Sum_y phi_n===!
      !--compute s2
      do n = 1 , Nj

         sn2 = 0.d0

         do j=1 , Ny/2-1
            sn2 = sn2 + exp( -sigma*(y(2*j)-y0)**2.d0) * sin(real(n)*pi*y(2*j)/yb)
         enddo

         !--compute s4
         sn4 = 0.d0

         do j=1, Ny/2
            sn4 = sn4 + exp( -sigma*(y(2*j-1)-y0)**2.d0) * sin(real(n)*pi*y(2*j-1)/yb)
         enddo

         !--compute total sum
         !~: dy/3.d0 * !(fun(1) + 2.d0*s2 + 4.d0*s4 + fun(Ny+1))
         QSumn(n) = (dy/3.d0) * (exp( -sigma*(y(1)-y0)**2.d0) * sin(real(n)*pi*y(1)/yb) + 2.d0*sn2 + 4.d0*sn4 + &
              exp( -sigma*(y(Ny+1)-y0)**2.d0) * sin(real(n)*pi*y(Ny+1)/yb) )

      enddo

    end subroutine simpson




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
          open(unit=70,file='outputY.dat',access='sequential',status='unknown')
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
