      program pipe

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c pipe.f
c Start Date: 31.1.2002
c Author: Marshall Ward
c
c For Part III Maths Essay, "Fast Flow in Curved Pipes: a Conundrum?" set by
c Stephen Cowley, at the University of Cambridge, UK.
c
c This program solves the equations for fluid flow in a curved pipe as the
c Dean number De->infinity. The equations to be solved are for the components
c of the fluid velocity field within the boundary layer, which the pipe is
c assumed to exhibit. The system is assumed to be symmetric for all
c cross-sections.
c
c They are written as functions of x and y. x is the angle along a
c cross-section of the pipe wall, with x=0 on the outer end and x=Pi (or -Pi)
c on the inner end. y scales like the distance from the wall. Identically,
c y = 1 - r/r0, where r is the radius from the center and r0 is some
c appropriate scaling parameter.
c
c The velocities in the boundary layer are written as {u, v, w}. u is the
c velocity along a cross-section of the pipe wall (along x). v is the velocity
c perpendicular to the wall face (along y). w is the flow into the pipe, and
c perpendicular to the cross section.
c
c As De->infinity, the u and v components of the inner core flow become zero,
c but the cross-sectional flow, w, is still nonzero. In this limit, w has no y
c dependence, but retains some x-dependence. This core flow velocity is
c labeled as wc.
c
c u, v, w, and wc are determined by a set of 3 PDEs and boundary conditions at
c y=0 and y->infinity. They are expanded in powers of x: u is in odd powers of
c x while v, and w are in even powers. The "order" to which I expand is
c determined by how many powers I keep. For example, for "zeroth order", I
c keep u up to x^1, with v and w up to x^0. "First order" keeps u up to x^3,
c and v and w up to x^2. After expanding in x, I'm left with three of ODEs
c in y, different for each order. I then solve these numerically, using the
c boundary conditions.
c
c I'll be using a finite difference method to solve the ODEs for each order.
c The resulting matrix should be block tridiagonal. I'll (hopefully!) use a
c method provided in LAPACK to solve this.
c
c There are three differential equations, sufficient for determining u, v, and
c w. wc is determined by an integral, which is obtained by solving the
c equations for inviscid flow in the core, with De->infinity. In this limit,
c there is no y-dependence, so wc's coefficients are all constant in the
c x-expansion.
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
c  inf   - approximate value for upper bound y->infinity (about y=10)
c  h     - step size = (inf - 0)/(m-1)
c
c  doublereal - velocity components
c  --------------------------------
c  u     - tangential parallel component of BL velocity. See intro notes.
c  v     - (negative) radial component of BL velocity. See intro notes.
c  w     - cross-sectional component of BL velocity. See intro notes.
c  wc    - cross-sectional component of core velocity. See intro notes.
c  wc0   - lowest order value of wc in x-expansion
c  q     - integral of u (to find wc0)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer order, m, n, iterate
      parameter (order = 0, m = 501, n = 3*m, iterate = 20)

      integer kl, ku, nrhs, ldab, ipiv(n), ldb, info
      parameter (kl = 5, ku = 5, nrhs = 1, ldab = 2*kl+ku+1, ldb = n)

      double precision inf, h, wc0
      parameter (inf = 1.4D1, wc0 = 1.763D0, h = inf/dble(m-1))

      double precision u(0:order, 0:m-1), v(0:order, 0:m-1),
     & w(0:order, 0:m-1), wc(0:order), q(0:order)
      double precision ab(ldab, n), b(ldb, nrhs)

c Local variables
      integer i, j

c Initialize parameters
      call initVelocity(u, v, w, wc, wc0, order, m, h)

c      do 50 i = 0, (m-1)
c         write(*,*) u(0,i)
c 50   continue

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
         call setupHighOrder(ab, b, h, kl, ku, u, v, w, wc, order, m, i)
         call dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)

          do 40 j = 0, m-1
            u(i,j) = b(3*j+1,1)
            v(i,j) = b(3*j+2,1)
            w(i,j) = b(3*j+3,1)
 40      continue

 30   continue

c Output results
c      do 60 i = 0, 10
c         write(*,*) u(0,i*(m-1)/10), v(0,i*(m-1)/10), w(0,i*(m-1)/10)
c 60   continue

      do 60 i = 0, (m-1)
          write(*,*) u(0,i)
 60   continue

      stop
      end

c End pipe

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c initVelocity - dumps zeros in u and v, and ramps w from 0 to wc0 as best
c  guess at zeroth order. all higher orders have zero b.c.'s.
c
c  These could be generalized, but fortrans a fuck
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine initVelocity(u, v, w, wc, wc0, k, m, h)
      integer m
      double precision u(0:k, 0:m-1), v(0:k, 0:m-1), w(0:k, 0:m-1),
     & wc(0:k), wc0, h
c Local variables
      integer i, j, r
      double precision node(0:10), d(11), e(10), f(11,1),
     & a(0:9), b(0:9), c(0:9)
      integer n, nrhs, ldb, info
      parameter(n = 11, nrhs = 1, ldb = 11)

c Fix boundary conditions and make initial guesses

c r is the number of points in each subinterval, x=j to x=j+1
      r = 1/h

      node(0)  = 0
      node(1)  = 0.629348
      node(2)  = 0.448817
      node(3)  = 0.196591
      node(4)  = 6.844524e-2
      node(5)  = 2.108094e-2
      node(6)  = 6.029014e-3
      node(7)  = 1.628765e-3
      node(8)  = 4.081907e-4
      node(9)  = 8.280066e-5
      node(10) = 0

c Setup matrix for spline coefficients
      d(1) = 2
      e(1) = 1
      f(1,1) = 6*((node(1)-node(0)) - 1.327501e+00)
      do 10 i = 2, 10
         d(i) = 4
         e(i) = 1
         f(i,1) = 6*(node(i-2)-2*node(i-1)+node(i))
 10   continue
      d(11) = 2
      f(11,1) = 6*(-4.060965e-05 - (node(10)-node(9))

c Solve for curvatures
      call dptsv(n, nrhs, d, e, f, ldb, info)

c Calculate coefficients
      do 20 i = 0, 9
         a(i) = (f(i+2,1)-f(i+1,1))/dble(6)
         b(i) = f(i+1)/dble(2)
         c(i) = (node(i+1)-node(i)) - (2*f(i+1,1)+f(i+2,1))/dble(6)
 20   continue

c Determine values for u
      do 60 i = 0, 9
         do 60 j = 0, r
            u(0,i*r+j) = node(i) + j*h*(node(i+1)-node(i))
 70      continue
 60   continue

      do 80 i = 10*r, m-1
         u(0,i) = 0
 80   continue

      do 90 i = 0, m-1
         v(0,i) = i*(-1.48318744D0)/dble(m-1)
         w(0,i) = i*(1.763D0)/dble(m-1)
 90   continue

      wc(0) = 1.763D0

      return
      end

c End initVelocity

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c setupZeroOrder - fill AB and B with the coefficients from the finite
c  difference equations. This should be called at the beginning of every
c  iteration.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine setupZeroOrder(ab, b, h, kl, ku, u, v, w, wc, k, m)
      integer kl, ku, k, m
      double precision ab(2*kl+ku+1, 3*m), b(3*m, 1), h,
     & u(0:k, 0:m-1), v(0:k, 0:m-1), w(0:k, 0:m-1), wc(0:k)
c Local variables
      integer i, j

c Some unassigned elements need to be zero. This step could use a lot of
c  improvement, but it's probably not worth the effort
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
         b(3*i-1,1) = -v(0,i-1) + v(0,i-2)
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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c setupHighOrder - fill AB and B with the coefficients from the finite
c  difference equations, for the nonlinear equations of order k. These
c  equations should be linear, so this only need to be called once per order
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine setupHighOrder(ab, b, h, kl, ku, u, v, w, wc, order,
     &                          m, k)
      integer kl, ku, order, k, m
      double precision ab(2*kl+ku+1, 3*m), b(3*m, 1), h,
     &     u(0:order, 0:m-1), v(0:order, 0:m-1), w(0:order, 0:m-1),
     &     wc(0:order)
c Local variables
      integer i, j, r, fac
      double precision temp, temp2

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

         ab(kl+ku+1-2+5,3*i-5) = 1 + h/2*v(0,i-1)
         ab(kl+ku+1-2+4,3*i-4) = 0
         ab(kl+ku+1-2+3,3*i-3) = 0
         ab(kl+ku+1-2+2,3*i-2) = -2 - 2*(h**2)*(k+1)*u(0,i-1)
         ab(kl+ku+1-2+1,3*i-1) = -h/2*(u(0,i-0)-u(0,i-2))
         ab(kl+ku+1-2-0,3*i-0) = -(h**2)*w(0,i-1)
         ab(kl+ku+1-2-1,3*i+1) = 1 - h/2*v(0,i-1)
         ab(kl+ku+1-2-2,3*i+2) = 0
         ab(kl+ku+1-2-3,3*i+3) = 0

         ab(kl+ku+1-1+5,3*i-5) = (2*k+1)*h/2
         ab(kl+ku+1-1+4,3*i-4) = -1
         ab(kl+ku+1-1+3,3*i-3) = 0
         ab(kl+ku+1-1+2,3*i-2) = (2*k+1)*h/2
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
         ab(kl+ku+1-0-0,3*i-0) = -2 - (h**2)*2*k*u(0,i-1)
         ab(kl+ku+1-0-1,3*i+1) = 0
         ab(kl+ku+1-0-2,3*i+2) = 0
         ab(kl+ku+1-0-3,3*i+3) = 1 - h/2*v(0,i-1)

 30   continue

c Setup the upper boundary of AB
c Might need to use a backward difference on v...
      ab(kl+ku+1,3*m-2) = 1

      ab(kl+ku+1-1+5,3*m-5) = (2*k+1)*h/2
      ab(kl+ku+1-1+4,3*m-4) = -1
      ab(kl+ku+1-1+2,3*m-2) = (2*k+1)*h/2
      ab(kl+ku+1-1+1,3*m-1) = 1

      ab(kl+ku+1,3*m-0) = 1

c Assign values to B

      b(1,1) = 0
      b(2,1) = 0
      b(3,1) = 0

c The 'u' part
      do 40 i = 2, m-1
         temp = 0
         if(k .GT. 1) then
            do 50 j = 1, k-1
               temp = temp + h/2*v(j,i-1)*(u(k-j,i-0)-u(k-j,i-2))
     &                + (h**2)*(2*(k-j)+1)*u(j,i-1)*u(k-j,i-1)
     &                + (h**2)/2*(w(j,i-1)*w(k-j,i-1) - wc(j)*wc(k-j))
 50         continue
         endif
         temp = temp - (h**2)*wc(0)*wc(k)

         do 60 j = 1, k
            temp2 = 0
            do 70 r = 0, k-j
               temp2 = temp2 + w(r,i-1)*w(k-j-r,i-1) - wc(r)*wc(k-j-r)
 70         continue
            temp = temp + (h**2)*((-1)**j)/(2*fac(2*j+1))*temp2
            b(3*i-2,1) = temp
 60      continue

c The 'v' part
         b(3*i-1,1) = 0

c The 'w' part
         temp = 0
         if(k .GT. 1) then
            do 80 j = 1, k-1
               temp = temp + 2*(h**2)*(k-j)*u(j,i-1)*w(k-j,i-1)
     &                  + h/2*v(j,i-1)*(w(k-j,i-0)-w(k-j,i-2))
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
      double precision q(0:k), wc(0:k)
c Local variables
      integer i, fac
      double precision temp

      if(k .EQ. 1) then
         wc(k) = -1/(2*q(0))
      else
         temp = ((-1)**k)/dble(fac(2*k-1))
         do 10 i = 1, k-1
            temp = temp - (2*(k-i)*q(i)*wc(k-i)
     &                    - ((-1)**k)/dble(fac(2*i+1)*fac(2*k-2*i-1)))
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer function fac(k)
      integer k
c Local variables
      integer i

      fac = 1

      do 10 i = 1, k
         fac = fac*i
 10   continue

      return
      end
