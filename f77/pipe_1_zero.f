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
c  order - Maximum order of calculation. See introductory notes (line 33?)
c  m     - number of points on the grid (including endpoints)
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
c  inf   - approximate value for y->infinity (about y=10)
c  h     - step size = (inf - 0)/(m-1)
c
c  doublereal - velocity components
c  --------------------------------
c  u     - tangential parallel component of BL velocity. See intro notes.
c  v     - (negative) radial component of BL velocity. See intro notes.
c  w     - cross-sectional component of BL velocity. See intro notes.
c  wc    - cross-sectional component of core velocity. See intro notes.
c  wc0   - lowest order value of wc in x-expansion
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer order, m, n, iterate
      parameter (order = 0, m = 20, n = 3*m, iterate = 500)

      integer kl, ku, nrhs, ldab, ipiv(n), ldb, info
      parameter (kl = 5, ku = 5, nrhs = 1, ldab = 2*kl+ku+1, ldb = n)

      double precision inf, h, wc0
      parameter (inf = 1.0D1, wc0 = 1.763D0, h = inf/dble(m-1))

      double precision u(0:order, 0:m-1), v(0:order, 0:m-1),
     & w(0:order, 0:m-1), wc(0:order)
      double precision ab(ldab, n), b(ldb, nrhs)

c Local variables
      integer i, j

c Initialize parameters
      call initVelocity(u, v, w, wc, wc0, order, m)

c Solve lowest order - nonlinear problem, iterative
      do 10 j = 1, iterate
         call setupZeroOrder(ab, b, h, kl, ku, u, v, w, wc, k, m)

         call dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
         do 20 i = 1, m
           u(0,i-1) = u(0,i-1) + b(3*i-2,1)
           v(0,i-1) = v(0,i-1) + b(3*i-1,1)
           w(0,i-1) = w(0,i-1) + b(3*i-0,1)
   20    continue
   10 continue

c for i = 1..order
c  solve ith order problem - linear, fast!

c Output results
      do 30 i = 0, m-1
         write(*,*) u(0,i)
   30 continue

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

      subroutine initVelocity(u, v, w, wc, wc0, k, m)
      integer order, m
      double precision u(0:k, 0:m-1), v(0:k, 0:m-1), w(0:k, 0:m-1),
     & wc(0:k), wc0
c Local variables
      integer i, j

c Fix boundary conditions and make initial guesses

      u(0,0) =  0.0000
      u(0,1) =  0.4906
      u(0,2) =  0.6316
      u(0,3) =  0.5637
      u(0,4) =  0.4182
      u(0,5) =  0.2761
      u(0,6) =  0.1683
      u(0,7) =  0.0971
      u(0,8) =  0.0538
      u(0,9) =  0.0290
      u(0,10) = 0.0153
      u(0,11) = 0.0079
      u(0,12) = 0.0040
      u(0,13) = 0.0020
      u(0,14) = 0.0010
      u(0,15) = 0.0005
      u(0,16) = 0.0002
      u(0,17) = 0.0001
      u(0,18) = 0.0000
      u(0,19) = 0.0000

      v(0,0) =   0.0000
      v(0,1) =  -0.1469
      v(0,2) =  -0.4546
      v(0,3) =  -0.7751
      v(0,4) =  -1.0347
      v(0,5) =  -1.2163
      v(0,6) =  -1.3315
      v(0,7) =  -1.3999
      v(0,8) =  -1.4386
      v(0,9) =  -1.4598
      v(0,10) = -1.4711
      v(0,11) = -1.4770
      v(0,12) = -1.4800
      v(0,13) = -1.4815
      v(0,14) = -1.4823
      v(0,15) = -1.4827
      v(0,16) = -1.4828
      v(0,17) = -1.4829
      v(0,18) = -1.4829
      v(0,19) = -1.4829

      w(0,0) =  0.0000
      w(0,1) =  0.4748
      w(0,2) =  0.9117
      w(0,3) =  1.2566
      w(0,4) =  1.4879
      w(0,5) =  1.6231
      w(0,6) =  1.6950
      w(0,7) =  1.7308
      w(0,8) =  1.7480
      w(0,9) =  1.7561
      w(0,10) = 1.7598
      w(0,11) = 1.7615
      w(0,12) = 1.7623
      w(0,13) = 1.7627
      w(0,14) = 1.7629
      w(0,15) = 1.7629
      w(0,16) = 1.7630
      w(0,17) = 1.7630
      w(0,18) = 1.7630
      w(0,19) = 1.7630

      wc(0) = 1.7630

c Fill higher order entries with zeros (i.e. "guess" zero)
      do 20 i = 1, k
         do 30 j = 0, m-1
            u(i,j) = 0
            v(i,j) = 0
            w(i,j) = 0
   30    continue
         wc(i) = 0
   20 continue

      return
      end

c End initVelocity

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c setupZeroOrder - fill AB and B with the coefficients from the finite
c  difference equations. This should be called at the beginning of every
c  iteration.
c
c Right now this is only for the iterative algorithm on the lowest order
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
   20    continue
   10 continue

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

   30 continue

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
   40 continue

      b(3*m-2,1) = 0
      b(3*m-1,1) = -v(0,m-1) + v(0,m-2) - h/2*(u(0,m-1)+u(0,m-2))
      b(3*m-0,1) = 0

      return
      end

c End constructSystem