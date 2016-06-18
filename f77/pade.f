      program pade

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c pade - reformulate the Taylor coefficients of f to Padé coefficients
c
c variables:
c  a : the function's Taylor coefficients
c  q : the lower polynomial
c  p : the upper polynomial
c  nn: the order of a
c  m : the order of q
c  n : the order of p
c  (nn = m+n)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer nn, m, n input, outp, outq
      parameter (nn = 100, m = 50, n = 50, input = 97, outp = 98,
     &           outq = 99)

      double precision a(0:nn), q(0:m), p(0:n), b(nn+1,nn+1)

      integer i, j, k
      double precision temp, xm

c *** Load Taylor Coefficients ***

      open(input, FILE='W_taylor', STATUS='OLD')

      do 400 i = 0, nn
         read(input,500) a(i)
 400  continue
 500  format((D24.4))

c *** Calculuate Padé coefficients ***

c Set initial values
      q(0) = 1
      p(0) = a(0)

c Construct matrix
      do 20 i = 1, nn

         do 30 j = 1, i-1
            if(j .LE. n) then
               b(i,j) = 0
            endif
 30      continue

         if(i .LE. n) then
            b(i,i) = 1
         endif

         do 40 j = i+1, nn
            b(i,j) = 0
 40      continue

         do 50 j = 1, i
            if(j .LE. m) then
               b(i,n+j) = -a(i-j)
            endif
 50      continue

         do 60 j = n+i+1, nn
            b(i,j) = 0
 60      continue

         b(i,nn+1) = a(i)

 20   continue

c Solve linear system
      do 70 i = n+1, nn-1

         k = i
         do 80 j = i, nn
            if ( abs(b(j,i)) .GT. abs(b(k,i)) ) then
               k = j
            endif
 80      continue

         if(b(k,i) .EQ. 0) then
            write(*,*) 'Singular!!!'
            stop
         endif

         if(k .NE. i) then
            do 90 j = i, nn+1
               temp = b(i,j)
               b(i,j) = b(k,j)
               b(k,j) = temp
 90         continue
         endif

         do 100 j = i+1, nn
            xm = b(j,i) / b(i,i)

            do 110 k = i+1, nn+1
               b(j,k) = b(j,k) - xm*b(i,k)
 110        continue

            b(i,j) = 0

 100     continue

 70   continue

      if(b(nn,nn) .EQ. 0) then
         write(*,*) 'Singular!!!'
         stop
      endif

c Calculate q's and p's
      if(m .GT. 0) then
         q(m) = b(nn,nn+1) / b(nn,nn)
      endif

      do 120 i = nn-1, n+1, -1
         temp = 0
         do 130 j = i+1, nn
            temp = temp + b(i,j)*q(j-n)
 130     continue
         q(i-n) = (b(i,nn+1) - temp) / b(i,i)
 120  continue

      do 140 i = n, 1, -1
         temp = 0
         do 150 j = n+1, nn
            temp = temp + b(i,j)*q(j-n)
 150     continue
         p(i) = b(i,nn+1) - temp
 140  continue

c *** Output Pad'e coefficients ***
c ARGH do these as files, I don't have time to cut and paste

      open(outp, FILE='p_out', STATUS='NEW')

      do 170 i = 0, n
         write(outp,*) p(i)
 170  continue

      close(outp)

      open(outq, FILE='q_out', STATUS='NEW')

       do 180 i = 0, m
         write(outq,*) q(i)
 180  continue

      close(outq)


      stop
      end