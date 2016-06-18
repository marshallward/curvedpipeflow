      program pipe

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c pipe.f
c Start Date: 31.1.2002
c Author: Marshall Ward
c
c For Part III Maths Essay, "Fast Flow in Curved Pipes: a Conundrum?" set by
c Stephen Cowley, at the University of Cambridge, UK.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Variables:
c
c  integer - algorithm information
c  -------------------------------
c  k     - Maximum order of calculation. See introductory notes
c  m     - number of points on the grid (including endpoints) - MUST BE ODD
c   (dgbsv variables)
c  kl    - number of subdiagonal
c  ku    - number of superdiagonals
c  nrhs  - number of diagonals for "B" (bb)
c  ldab  - leading dimension of ab ( = 2*kl+ku+1)
c  ipiv  - n-array of pivot indices
c  ldb   - leading dimension of bb ( = n)
c  info  - a status variable
c
c  doublereal - algorithm information
c  ----------------------------------
c  ab    - matrix A in linear eqn to be solved for X, A*X=B
c  b     - vector B in linear eqn to be solved for X, A*X=B
c
c
c  doublereal - grid information
c  ------------------------------
c  (lower bound is always y=0)
c  inf   - approximate value for upper bound y->infinity (about y=35)
c  h     - step size = (inf - 0)/(m-1)
c
c  doublereal - velocity components
c  --------------------------------
c  u     - tangential parallel component of BL velocity. See intro notes.
c  v     - (negative) radial component of BL velocity. See intro notes.
c  w     - cross-sectional component of BL velocity. See intro notes.
c  wc    - cross-sectional component of core velocity. See intro notes.c  wc0   - lowest order value of wc in x-expansion
c  q     - integral of u (to find wc0)
c
c
c Warnings:
c  Numerical integrator requires an even number of intervals
c   -> m must be odd
c  I use 11 nodes/10 intervals in the cubic splines
c   -> I should use a multiple of k*10+1 grid points
c  I do a inf/(m-1)... I don't really know what effect this has
c
c Performance:
c  inf = 35 seems to give machine double precision
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer order, m, n, iterate
      parameter (order = 0, m = 101, n = 3*m, iterate = 20)

      integer kl, ku, nrhs, ldab, ipiv(n), ldb, info
      parameter (kl = 5, ku = 5, nrhs = 1, ldab = 2*kl+ku+1, ldb = n)

      double precision inf, h, wc0, fac
      parameter (inf = 2.0D1, wc0 = 1.250D0, h = inf/dble(m-1))

      double precision u(0:order, 0:m-1), v(0:order, 0:m-1),
     & w(0:order, 0:m-1), wc(0:order), q(0:order)
      double precision ab(ldab, n), b(ldb, nrhs)

c Local variables
      integer i, j

c Initialize parameters
      call initVelocity(u, v, w, wc, wc0, order, m, h)

c DumbInit
c      do 5 i = 0, m-1
c         u(i) = 0
c         v(i) = (line)
c         w(i) = (line)
c 5    continue

c Solve lowest order - nonlinear problem, iterative
      do 10 j = 1, iterate
         call setupZeroOrder(ab, b, h, kl, ku, u, v, w, wc, order, m)

         call dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
         do 20 i = 1, m
           u(0,i-1) = u(0,i-1) + b(3*i-2,1)
           v(0,i-1) = v(0,i-1) + b(3*i-1,1)
           w(0,i-1) = w(0,i-1) + b(3*i-0,1)
 20      continue
 10   continue

c solve ith order problem - linear, fast!
      do 30 i = 1, order
c   Get wc(k)
         call integrateU(u, q, h, i-1, order, m)
         call getWc(q, wc, i, order)

c   Solve u(k,i) et. al.
         call setupHighOrder(ab, b, h, kl, ku, u, v, w, wc, order, m,
     &                         i)

         call dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)

         do 40 j = 0, m-1
            u(i,j) = b(3*j+1,1)
            v(i,j) = b(3*j+2,1)
            w(i,j) = b(3*j+3,1)
 40      continue

 30   continue

c Output results
c      do 50 i = 0, order
c         write(*,*) wc(i)
c         write(*,*) dble(0)
c 50   continue


c      do 60 i = 0, int(inf)
c         write(*,*) u(order,i*(m-1)/int(inf)),
c     &              v(order,i*(m-1)/int(inf)),
c     &              w(order,i*(m-1)/int(inf))
c 60   continue

      do 70 i = 0, (m-1)
          write(*,*) v(0,i)
 70   continue

      stop
      end

c End pipe

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c initVelocity - Uses a cubic spline interpolation between 11 points
c  from x=0 to x=10 (stored on file 'nodes').
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine initVelocity(u, v, w, wc, wc0, k, m, h)
      integer m
      double precision u(0:k, 0:m-1), v(0:k, 0:m-1), w(0:k, 0:m-1),
     & wc(0:k), wc0, h
c Local variables
      integer i, j, r, file
      double precision node(3,0:10), d(11), e(10),
     & fu(11,1), fv(11,1), fw(11,1),
     & a(3,0:9), b(3,0:9), c(3,0:9)
      integer n, nrhs, ldb, info
      parameter(n = 11, nrhs = 1, ldb = 11, file = 99)
      double precision uamp, wamp, vmax
      parameter(wamp = 7.7971D0, vmax = 1.483238407D0)

      uamp = (wc0 / vmax)*wamp

c Fix boundary conditions and make initial guesses

c r is the number of points in each subinterval, x=j to x=j+1
      r = 1/h

c u's nodes
c learn to read files someday, retard
      open(file, FILE='nodes', STATUS='OLD')

      do 100 i = 0, 10
         read(file,500) node(1,i), node(2,i), node(3,i)
 100  continue
 500  format(3(E24.4))

      close(file)

c Setup matrix for spline coefficients
c BAD GUESSES: A and B for v and w
      fu(1,1) = 6*((node(1,1)-node(1,0)) - 1.327501)
      fv(1,1) = 6*((node(2,1)-node(2,0)) + 0.4214)
      fw(1,1) = 6*((node(3,1)-node(3,0)) - 0.9085)
      do 10 i = 2, n-1
         fu(i,1) = 6*(node(1,i-2)-2*node(1,i-1)+node(1,i))
         fv(i,1) = 6*(node(2,i-2)-2*node(2,i-1)+node(2,i))
         fw(i,1) = 6*(node(3,i-2)-2*node(3,i-1)+node(3,i))
 10   continue
      fu(n,1) = 6*(-4.060965e-05 - (node(1,n-1)-node(1,n-2)))
      fv(n,1) = 6*(-3.1986e-05 - (node(2,n-1)-node(2,n-2)))
      fw(n,1) = 6*(9.7376e-06 - (node(3,n-1)-node(3,n-2)))

c Solve for curvatures
      d(1) = 2
      e(1) = 1
      do 20 i = 2, n-1
         d(i) = 4
         e(i) = 1
 20   continue
      d(n) = 2
      call dptsv(n, nrhs, d, e, fu, ldb, info)

      d(1) = 2
      e(1) = 1
      do 30 i = 2, n-1
         d(i) = 4
         e(i) = 1
 30   continue
      d(n) = 2
      call dptsv(n, nrhs, d, e, fv, ldb, info)

      d(1) = 2
      e(1) = 1
      do 40 i = 2, n-1
         d(i) = 4
         e(i) = 1
 40   continue
      d(n) = 2
      call dptsv(n, nrhs, d, e, fw, ldb, info)

c Calculate coefficients
      do 50 i = 0, n-2
         a(1,i) = (fu(i+2,1)-fu(i+1,1))/dble(6)
         b(1,i) = fu(i+1,1)/dble(2)
         c(1,i) = (node(1,i+1)-node(1,i))
     &               - (2*fu(i+1,1)+fu(i+2,1))/dble(6)

         a(2,i) = (fv(i+2,1)-fv(i+1,1))/dble(6)
         b(2,i) = fv(i+1,1)/dble(2)
         c(2,i) = (node(2,i+1)-node(2,i))
     &               - (2*fv(i+1,1)+fv(i+2,1))/dble(6)

         a(3,i) = (fw(i+2,1)-fw(i+1,1))/dble(6)
         b(3,i) = fw(i+1,1)/dble(2)
         c(3,i) = (node(3,i+1)-node(3,i))
     &               - (2*fw(i+1,1)+fw(i+2,1))/dble(6)
 50   continue

c Determine values for (u,v,w)
      do 60 i = 0, n-2
         do 70 j = 0, r
            u(0,i*r+j) = a(1,i)*(h*j)**3 + b(1,i)*(h*j)**2
     &                    + c(1,i)*(h*j) + node(1,i)
            v(0,i*r+j) = a(2,i)*(h*j)**3 + b(2,i)*(h*j)**2
     &                    + c(2,i)*(h*j) + node(2,i)
            w(0,i*r+j) = a(3,i)*(h*j)**3 + b(3,i)*(h*j)**2
     &                    + c(3,i)*(h*j) + node(3,i)
 70      continue
 60   continue

      do 80 i = 10*r, m-1
c         u(0,i) = uamp*(i*h)*exp(-vmax*i*h)
c         v(0,i) = -vmax + uamp/vmax*i*h*exp(-vmax*i*h)
c         w(0,i) = wc0 - wamp*exp(-vmax*i*h)
         u(0,i) = 0
         v(0,i) = -vmax
         w(0,i) = wc0
 80   continue

      wc(0) = wc0

      return
      end

c End initVelocity

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c setupZeroOrder - fill AB and B with the coefficients from the finite
c  difference equations. This should be called at the beginning of
c  every iteration.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine setupZeroOrder(ab, b, h, kl, ku, u, v, w, wc, k, m)
      integer kl, ku, k, m
      double precision ab(2*kl+ku+1, 3*m), b(3*m, 1), h,
     & u(0:k, 0:m-1), v(0:k, 0:m-1), w(0:k, 0:m-1), wc(0:k)
c Local variables
      integer i, j

c Some unassigned elements need to be zero. This step could use a lot\
c of improvement, but it's probably not worth the effort
c
c In fact, if my suspicion is correct and I can reduce the bandwidth from
c  6 to 4, and then this is completely unnecessary. Keep note of that.
      do 10 j = 1, 3*m
         do 20 i = kl+1, 2*kl+ku+1
            ab(i,j) = 0
 20      continue
 10   continue

c Assign values to AB
c Recall that AB(kl+ku+1+i-j,j) = A(i,j)
c  This compresses data by stacking diagonals on top of each other, with some
c  column offset

c Setup the lower boundary of AB
      ab(kl+ku+1,1) = 1
      ab(kl+ku+1,2) = 1
      ab(kl+ku+1,3) = 1

c Setup the middle grid points of AB
c DONE
      do 30 i = 2, m-1

         ab(kl+ku+1-2+5,3*i-5) = 1 + h/2*v(0,i-1)
         ab(kl+ku+1-2+4,3*i-4) = 0
         ab(kl+ku+1-2+3,3*i-3) = 0
         ab(kl+ku+1-2+2,3*i-2) = -2 - 2*(h**2)*u(0,i-1)
         ab(kl+ku+1-2+1,3*i-1) = -h/2*(u(0,i-0)-u(0,i-2))
         ab(kl+ku+1-2-0,3*i-0) = -(h**2)*w(0,i-1)
         ab(kl+ku+1-2-1,3*i+1) = 1 - h/2*v(0,i-1)
         ab(kl+ku+1-2-2,3*i+2) = 0
         ab(kl+ku+1-2-3,3*i+3) = 0

         ab(kl+ku+1-1+5,3*i-5) = h/2
         ab(kl+ku+1-1+4,3*i-4) = -1
         ab(kl+ku+1-1+3,3*i-3) = 0
         ab(kl+ku+1-1+2,3*i-2) = h/2
         ab(kl+ku+1-1+1,3*i-1) = 1
         ab(kl+ku+1-1-0,3*i-0) = 0
         ab(kl+ku+1-1-1,3*i+1) = 0
         ab(kl+ku+1-1-2,3*i+2) = 0
         ab(kl+ku+1-1-3,3*i+3) = 0

         ab(kl+ku+1-0+5,3*i-5) = 0
         ab(kl+ku+1-0+4,3*i-4) = 0
         ab(kl+ku+1-0+3,3*i-3) = 1 + h/2*v(0,i-1)
         ab(kl+ku+1-0+2,3*i-2) = 0
         ab(kl+ku+1-0+1,3*i-1) = -h/2*(w(0,i-0)-w(0,i-2))
         ab(kl+ku+1-0-0,3*i-0) = -2
         ab(kl+ku+1-0-1,3*i+1) = 0
         ab(kl+ku+1-0-2,3*i+2) = 0
         ab(kl+ku+1-0-3,3*i+3) = 1 - h/2*v(0,i-1)

 30   continue

c Setup the upper boundary of AB
c Might need to use a backward difference on v...
      ab(kl+ku+1,3*m-2) = 1

      ab(kl+ku+1-1+5,3*m-5) = h/2
      ab(kl+ku+1-1+4,3*m-4) = -1
      ab(kl+ku+1-1+2,3*m-2) = h/2
      ab(kl+ku+1-1+1,3*m-1) = 1

      ab(kl+ku+1,3*m-0) = 1

c Assign values to b for this iteration
c Gotta fix this for the new system

      b(1,1) = 0
      b(2,1) = 0
      b(3,1) = 0

c Fix this
      do 40 i = 2, m-1
         b(3*i-2,1) = -u(0,i-2) + 2*u(0,i-1) - u(0,i)
     &                       + (h*u(0,i-1))**2
     &                       + h/2*(u(0,i)-u(0,i-2))*v(0,i-1)
     &                       + (h**2)/2*(w(0,i-1)**2 - wc(0)**2)
         b(3*i-1,1) = -(v(0,i-1) - v(0,i-2))
     &               - h/2*(u(0,i-1)+u(0,i-2))
         b(3*i-0,1) = -w(0,i-2) + 2*w(0,i-1) - w(0,i)
     &               + h/2*v(0,i-1)*(w(0,i) - w(0,i-2))
 40   continue

      b(3*m-2,1) = 0
      b(3*m-1,1) = -v(0,m-1) + v(0,m-2) - h/2*(u(0,m-1)+u(0,m-2))
      b(3*m-0,1) = 0

      return
      end

c End setupZeroOrder

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c setupHighOrder - fill AB and B with the coefficients from the finite
c  difference equations, for the nonlinear equations of order k. These
c  equations should be linear, so this only need to be called once per
c  order
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine setupHighOrder(ab, b, h, kl, ku, u, v, w, wc, order,
     &                          m, k)
      integer kl, ku, order, k, m
      double precision ab(2*kl+ku+1, 3*m), b(3*m, 1), h,
     &     u(0:order, 0:m-1), v(0:order, 0:m-1), w(0:order, 0:m-1),
     &     wc(0:order)
c Local variables
      integer i, j, r
      double precision temp, temp2, fac, a
      parameter(a = 2)

c Zero out AB
      do 10 j = 1, 3*m
         do 20 i = kl+1, 2*kl+ku+1
            ab(i,j) = 0
 20      continue
 10   continue

c Assign values to AB
c Recall that AB(kl+ku+1+i-j,j) = A(i,j)
c  This compresses data by stacking diagonals on top of each other, with some
c  column offset

c Setup the lower boundary of AB
      ab(kl+ku+1,1) = 1
      ab(kl+ku+1,2) = 1
      ab(kl+ku+1,3) = 1

c Setup the middle grid points of AB
      do 30 i = 2, m-1

         ab(kl+ku+1-2+5,3*i-5) = (a**(-k-0.5))
     &                                      + h/2*(a**(-0.5))*v(0,i-1)
         ab(kl+ku+1-2+4,3*i-4) = 0
         ab(kl+ku+1-2+3,3*i-3) = 0
         ab(kl+ku+1-2+2,3*i-2) = -2*(a**(-k-0.5))
     &                                       - 2*(h**2)*(k+1)*u(0,i-1)
         ab(kl+ku+1-2+1,3*i-1) = -(a**(-0.5))*h/2*(u(0,i-0)-u(0,i-2))
         ab(kl+ku+1-2-0,3*i-0) = -(a**(-1))*(h**2)*w(0,i-1)
         ab(kl+ku+1-2-1,3*i+1) = (a**(-k-0.5))-(a**(-0.5))*h/2*v(0,i-1)
         ab(kl+ku+1-2-2,3*i+2) = 0
         ab(kl+ku+1-2-3,3*i+3) = 0

         ab(kl+ku+1-1+5,3*i-5) = (2*k+1)*h/2
         ab(kl+ku+1-1+4,3*i-4) = -1*(a**(-0.5))
         ab(kl+ku+1-1+3,3*i-3) = 0
         ab(kl+ku+1-1+2,3*i-2) = (2*k+1)*h/2
         ab(kl+ku+1-1+1,3*i-1) = 1*(a**(-0.5))
         ab(kl+ku+1-1-0,3*i-0) = 0
         ab(kl+ku+1-1-1,3*i+1) = 0
         ab(kl+ku+1-1-2,3*i+2) = 0
         ab(kl+ku+1-1-3,3*i+3) = 0

         ab(kl+ku+1-0+5,3*i-5) = 0
         ab(kl+ku+1-0+4,3*i-4) = 0
         ab(kl+ku+1-0+3,3*i-3) = (a**(-k)) + h/2*v(0,i-1)
         ab(kl+ku+1-0+2,3*i-2) = 0
         ab(kl+ku+1-0+1,3*i-1) = -h/2*(w(0,i-0)-w(0,i-2))
         ab(kl+ku+1-0-0,3*i-0) = -2*(a**(-k))
     &                               - (a**(0.5))*(h**2)*2*k*u(0,i-1)
         ab(kl+ku+1-0-1,3*i+1) = 0
         ab(kl+ku+1-0-2,3*i+2) = 0
         ab(kl+ku+1-0-3,3*i+3) = a**(-k) - h/2*v(0,i-1)

 30   continue

c Setup the upper boundary of AB
      ab(kl+ku+1,3*m-2) = 1

      ab(kl+ku+1-1+5,3*m-5) = (2*k+1)*h/2
      ab(kl+ku+1-1+4,3*m-4) = -1*(a**(-0.5))
      ab(kl+ku+1-1+2,3*m-2) = (2*k+1)*h/2
      ab(kl+ku+1-1+1,3*m-1) = 1*(a**(-0.5))

      ab(kl+ku+1,3*m-0) = 1

c Assign values to B

      b(1,1) = 0
      b(2,1) = 0
      b(3,1) = 0

c The 'u' part
      do 40 i = 2, m-1
         temp = 0
         if(k .GT. 1) then
            do 50 j = 1, (k-1)
               temp = temp + h/2*v(j,i-1)*(u(k-j,i-0)-u(k-j,i-2))
     &                                                    *(a**(-0.5))
     &                + (h**2)*(2*(k-j)+1)*u(j,i-1)*u(k-j,i-1)
     &                + (h**2)/2*(w(j,i-1)*w(k-j,i-1) - wc(j)*wc(k-j))
     &                                                    *(a**(-1))

 50         continue
         endif
         temp = temp - (h**2)*wc(0)*wc(k)*(a**(-1))

         do 60 j = 1, k
            temp2 = 0
            do 70 r = 0, (k-j)
               temp2 = temp2 + w(r,i-1)*w(k-j-r,i-1) - wc(r)*wc(k-j-r)
 70         continue
            temp = temp + (h**2)*((-1)**j)/dble(2*fac(2*j+1))*temp2
     &                                                      *(a**(-1))
 60      continue
         b(3*i-2,1) = temp

c The 'v' part
         b(3*i-1,1) = 0

c The 'w' part
         temp = 0
         if(k .GT. 1) then
            do 80 j = 1, (k-1)
               temp = temp + 2*(h**2)*(k-j)*u(j,i-1)*w(k-j,i-1)
     &                                                   *(a**(0.5))
     &                    + h/2*v(j,i-1)*(w(k-j,i-0)-w(k-j,i-2))
 80         continue
         endif
         b(3*i-0,1) = temp

 40   continue

      b(3*m-2,1) = 0
      b(3*m-1,1) = 0
      b(3*m-0,1) = wc(k)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c integrateU - integrates u's kth order term and puts it in q
c  --> m must be odd since I'm using a composite Simpson's rule
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine integrateU(u, q, h, k, order , m)
      integer k, order, m
      double precision u(0:order, 0:m-1), q(0:order), h
c Local variables
      integer i
      double precision xi0, xi1, xi2

      xi0 = u(k, 0) + u(k, m-1)
      xi1 = 0
      xi2 = 0

      do 10 i = 1, (m-3)/2
         xi2 = xi2 + u(k, 2*i)
         xi1 = xi1 + u(k, 2*i-1)
 10   continue
      xi1 = xi1 + u(k, m-2)

      q(k) = h/3*(xi0 + 2*xi2 + 4*xi1)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c getWc - calculates the kth order wc term (core flow) based on k-1 terms
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getWc(q, wc, k, order)
      integer k, order
      double precision q(0:order), wc(0:order)
c Local variables
      integer i
      double precision temp, fac, a
      parameter(a = 2)

      if(k .EQ. 1) then
         wc(k) = -1/(2*q(0))*(a**(-2*k-0.5))
      else
         temp = ((-1)**k)/dble(fac(2*k-1))*(a**(-2*k-0.5))
         do 10 i = 1, k-1
            temp = temp - (2*(k-i)*q(i)*wc(k-i)
     &       - ((-1)**k)/dble(fac(2*i+1)*fac(2*k-2*i-1)))
     &                                              *(a**(-2*k-0.5))
 10      continue
         wc(k) = temp/(2*k*q(0))
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c fac - returns factorial of k
c  Dunno if it handles 0! = 1
c
c Hah... major problems then this factorial gets too big. Need to switch to
c  float
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function fac(k)
      integer k
c Local variables
      integer i

      fac = 1

      do 10 i = 1, k
         fac = fac*i
 10   continue

      return
      end
