      module libs
      use iso_c_binding

      abstract interface
        function getfunval_proc(r1) result(funval) bind(c)
          import :: c_double 
          real(c_double), intent(in), value :: r1 
          real(c_double) :: funval 
        end function  
      end interface 

      contains

        subroutine integrate(getfunvalfromc, x0, x1, n, outputreal) 
     \    bind(c, name='Integrate')
          type(c_funptr), intent(in), value  :: getfunvalfromc
          real(c_double), intent(in)         :: x0, x1  
          integer(c_int), intent(in)         :: n 
          procedure(getfunval_proc), pointer :: getfunval 
          real(c_double), intent(out)        :: outputreal 
          real*8                             :: s, h, xi 
          integer*4                          :: m 
          m = n  
          call c_f_procpointer(cptr=getfunvalfromc, fptr=getfunval)
          if((m/2)*2 .ne. m) m = m + 1 
          s = .0d0 
          h = (x1-x0)/dble(m)
          do i=2, m-2, 2
            xi = x0 + dble(i)*h 
            s = s + 2.d0 * getfunval(xi) + 4.d0 * getfunval(xi+h)
          enddo 
          outputreal = ( s+getfunval(x0)+getfunval(x1)+
     \                   4.d0*getfunval(x0+h) ) * h/3.d0
          nullify(getfunval)
        end subroutine 

        subroutine get_determinant(n, a, det) 
     \    bind(c, name='GetDeterminant')
          integer(c_int), intent(in) :: n 
          real(c_double), intent(in) :: a(n,n)
          real(c_double), intent(out) :: det 
          det = deter(n, a)
        end subroutine 

        recursive subroutine fft(n, x) 
     \    bind(c, name='FastFourierTransform')
          integer(c_int), intent(in)                             :: n  
          complex(c_double_complex), intent(inout), dimension(n) :: x
          complex*16, allocatable, dimension(:) :: even, odd 
          complex*16 t 
          real*8 pi
          if(n.le.1) return 
          allocate(odd((n+1)/2))
          allocate(even(n/2))
          pi = dacos(-1.d0)
          odd = x(1:n:2)
          even = x(2:n:2)
          call fft((n+1)/2, odd)
          call fft(n/2, even)
          do i=1,n/2
            t = exp(dcmplx(.0d0, -2.d0*pi*dble(i-1)/dble(n)))*even(i)
            x(i)     = odd(i) + t 
            x(i+n/2) = odd(i) - t  
          enddo 
          deallocate(odd)
          deallocate(even)
        end subroutine 

        subroutine least_square_fit(n, x, y, m, q, r2, sem, seq) 
     \           bind(c, name='LeastSquareFit')
          real(c_double), intent(out), dimension(n) :: x, y
          real(c_double), intent(out)               :: m, q, r2, sem, 
     \  seq
          real*8 ssxx, ssyy, ssxy, s, xm, ym
          xm = sum(x(:)) / dble(n)
          ym = sum(y(:)) / dble(n)
          ssxx = sum((x(:) - xm)**2)
          ssyy = sum((y(:) - ym)**2)
          ssxy = sum((x(:) - xm)*(y(:) - ym))
          m = ssxy / ssxx 
          q = ym - m*xm 
          r2 = ssxy**2 / (ssxx*ssyy)
          s = dsqrt((ssyy - m*ssxy)/dble(n-2))
          sem = s / dsqrt(ssxx)
          seq = s * dsqrt(1.d0/dble(n) + xm**2/ssxx)   
        end subroutine

        subroutine jacobi_diag(a, q, w) 
     \           bind(c, name='JacobiDiagonalization')
          real(c_double), intent(inout), dimension(3,3) :: a
          real(c_double), intent(out), dimension(3,3)   :: q
          real(c_double), intent(out), dimension(3)     :: w
          integer, parameter :: n = 3
          real*8 sd, so
          real*8 s, c, t
          real*8 g, h, z, theta
          real*8 thresh
          integer i, x, y, r
          do x=1,n
            q(x,x) = 1.d0
            do y=1,x-1
              q(x,y) = .0d0
              q(y,x) = .0d0 
            enddo
          enddo
          do x=1,n
            w(x) = a(x,x)
          enddo
          sd = .0d0 
          do x=1,n 
            sd = sd + abs(w(x))
          enddo 
          sd = sd**2
          do i=1,50
            so = .0d0  
            do x=1,n 
              do y=x+1,n 
                so = so + abs(a(x,y))
              enddo 
            enddo 
            if (so .eq. .0d0) then
              call eigsrt(n, w, q) 
              return
            endif
            if (i .lt. 4) then
              thresh = .2d0*so/n**2 
            else 
              thresh = .0d0 
            endif 
            do x=1,n 
              do y=x+1,n 
                g = 100.d0 * abs(a(x,y))
                if (i .gt. 4 .and. abs(w(x))+g .eq. abs(w(x))
     \                     .and. abs(w(y))+g .eq. abs(w(y))) then 
                  a(x,y) = .0
                else if (abs(a(x,y)) .gt. thresh) then 
                  h = w(y) - w(x)
                  if (abs(h)+g .eq. abs(h)) then 
                    t = a(x,y) / h 
                  else 
                    theta = .5d0 * h / a(x,y)
                    if (theta .lt. .0d0) then 
                      t = -1.d0 / (dsqrt(1.d0 + theta**2) - theta)
                    else 
                      t = +1.d0 / (dsqrt(1.d0 + theta**2) + theta)
                    endif 
                  endif 
                  c = 1.d0 / dsqrt(1.d0 + t**2)
                  s = t*c 
                  z = t*a(x,y)
                  a(x,y) = .0d0 
                  w(x) = w(x) - z 
                  w(y) = w(y) + z 
                  do r=1,x-1 
                    t = a(r,x)
                    a(r,x) = c*t - s*a(r,y)
                    a(r,y) = s*t + c*a(r,y)
                  enddo 
                  do r=x+1,y-1 
                    t = a(x,r)
                    a(x,r) = c*t - s*a(r,y)
                    a(r,y) = s*t + c*a(r,y)
                  enddo 
                  do r=y+1,n
                    t = a(x,r)
                    a(x,r) = c*t - s*a(y,r)
                    a(y,r) = s*t + c*a(y,r)
                  enddo 
                  do r=1,n 
                    t = q(r,x)
                    q(r,x) = c*t - s*q(r,y)
                    q(r,y) = s*t + c*q(r,y)
                  enddo 
                endif 
              enddo 
            enddo 
          enddo 
          print *, 'Jacobi Diagonalization: No convergence'
        end subroutine

        subroutine gaussj(n, m, a, b) bind(c, name='GaussJordan')
          !!! Pag. 1014 Numerical Recipes in Fortran 90
          integer(c_int), intent(in)                    :: n, m 
          real(c_double), intent(inout), dimension(n,n) :: a 
          real(c_double), intent(inout), dimension(n,m) :: b 
          integer, dimension(n) :: ipiv, indxr, indxc 
          logical, dimension(n) :: lpiv 
          real*8 pivinv 
          real*8, dimension(n) :: dumc
          integer, target, dimension(2) :: irc
          integer, pointer :: irow, icol 
          irow => irc(1)
          icol => irc(2)
          ipiv = 0 
          do i=1,n
            lpiv = (ipiv == 0)
            irc = maxloc(array=abs(a), 
     \                 mask=outerand(lpiv, lpiv))
            ipiv(icol) = ipiv(icol) + 1
            if(ipiv(icol) .gt. 1) then
              print '(a)', 'nrerror: gaussj: singular matrix (1)' 
            endif  
            if(irow .ne. icol) then
              call swap(a(irow,:), a(icol,:))
              call swap(b(irow,:), b(icol,:))
            endif 
            indxr(i) = irow 
            indxc(i) = icol 
            if(a(icol,icol) .eq. .0d0) then 
              print '(a)', 'nrerror: gaussj: singular matrix (2)'
            endif 
            pivinv = 1.d0 / a(icol, icol)
            a(icol, icol) = 1.d0 
            a(icol, :) = a(icol, :) * pivinv 
            b(icol, :) = b(icol, :) * pivinv 
            dumc = a(:, icol)
            a(:, icol) = .0d0 
            a(icol, icol) = pivinv 
            a(:icol-1, :) = a(:icol-1, :) - outerprod(dumc(:icol-1), 
     \                                              a(icol,:))
            b(:icol-1, :) = b(:icol-1, :) - outerprod(dumc(:icol-1), 
     \                                              b(icol,:))
            a(icol+1:, :) = a(icol+1:, :) - outerprod(dumc(icol+1:), 
     \                                              a(icol,:))
            b(icol+1:, :) = b(icol+1:, :) - outerprod(dumc(icol+1:), 
     \                                              b(icol,:))
          enddo   
          do l=n,1,-1
            call swap(a(:,indxr(l)), a(:,indxc(l))) 
          enddo 
        end subroutine 

        function outerprod(a, b)
          real*8, dimension(:), intent(in)    :: a, b 
          real*8, dimension(size(a), size(b)) :: outerprod
          outerprod = spread(a,dim=2,ncopies=size(b)) *
     \                spread(b,dim=1,ncopies=size(a))
        end function 

        function outerand(a, b)
          logical, dimension(:), intent(in)    :: a, b 
          logical, dimension(size(a), size(b)) :: outerand 
          outerand = spread(a,dim=2,ncopies=size(b)) .and. 
     \               spread(b,dim=1,ncopies=size(a))
        end function 

        subroutine swap(a, b)
          real*8, dimension(:), intent(inout) :: a, b 
          real*8, dimension(size(a))          :: dum 
          dum = a 
          a = b 
          b = dum 
        end subroutine 

        subroutine swap0(a, b)
          real*8, intent(inout) :: a, b 
          real*8 dum 
          dum = a 
          a = b 
          b = dum 
        end subroutine 

        subroutine eigsrt(n, d, v)
          integer, intent(in)                   :: n
          real*8, dimension(n), intent(inout)   :: d 
          real*8, dimension(n,n), intent(inout) :: v 
          do i=1,n-1 
            j = imaxloc(d(i:n)) + i - 1
            if(j.ne.i) then
              call swap0(d(i), d(j))
              call swap(v(:,i), v(:,j)) 
            endif 
          enddo
        end subroutine 

        function imaxloc(arr)
          real*8, intent(in), dimension(:) :: arr 
          integer imaxloc 
          integer, dimension(1) :: imax 
          imax = maxloc(arr)
          imaxloc = imax(1)
        end function 

        subroutine assert(l1, string)
          logical, intent(in) :: l1 
          character(len=*), intent(in) :: string 
          if(.not.l1) then
            print *, 'nrerror: an asseration failed with this tag:', 
     \                string
            stop  
          endif 
        end subroutine 

        recursive function deter(n, mat) result( accum )
          integer                      :: n 
          real*8, dimension(n,n)      :: mat
          real*8                      :: accum   
          real*8, dimension(n-1, n-1) :: submat
          integer sgn 
          if(n.eq.1) then 
            accum = mat(1,1)
          else 
            accum = .0d0 
            sgn = 1
            do i=1,n
              submat(:, :i-1) = mat(2:,   :i-1)
              submat(:,i:   ) = mat(2:,i+1:   )
              accum = accum + sgn*mat(1, i)*deter(n-1, submat) 
              sgn = -sgn
            enddo  
          endif 
        end function 

      end module
