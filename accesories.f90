!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

        SUBROUTINE WALLCLOCK_LIMIT (total_min)

        character*8  :: str1
        character*1  :: str2
        character*1  :: str3
        character*1  :: str4
        character*1  :: str5
        integer      :: h

        str1 = 'a'

        open (14, file='RUN', status='unknown')

        do while (str1.ne.'#BSUB -W')  

         read (14,901) str1, str2, str3, str4, str5
901      format(A, TR1, A, A, TR1, A, A)

         h2 = ichar(str2)-48
         h1 = ichar(str3)-48
         m2 = ichar(str4)-48
         m1 = ichar(str5)-48

         total_min = real(h2*600+h1*60+m2*10+m1)

        end do

        close (14)

        RETURN

        END SUBROUTINE WALLCLOCK_LIMIT

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

        Function ranff(idum,iy,iv)
!
!       Minimal random number generator of Park and Miller with
!       Bays-Durham shuffle and added safegards, see Numerical Recipes
!
        Integer idum, IA,IM,IQ,IR,NTAB,NDIV
        Real ranff,AM,EPS,RNMX
        Parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,   &
             NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        INTEGER j,k,iv(NTAB),iy

        if (idum.le.0 .or. iy.eq.0) then
                idum=max(-idum,1)
                do j=NTAB+8,1,-1
                   k=idum/IQ
                   idum=IA*(idum-k*IQ)-IR*k
                   if (idum.lt.0) idum=idum+IM
                   if (j.le.NTAB) iv(j)=idum
                enddo
                iy=iv(1)
        endif
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if(idum.lt.0) idum=idum+IM
        j=1+iy/NDIV
        iy=iv(j)
        iv(j)=idum
        ranff=min(AM*iy,RNMX)

        return
        END

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

!RANDOM
        Function ranpp(idum,iy1,iv1)
!
!       Minimal random number generator of Park and Miller with
!       Bays-Durham shuffle and added safegards, see Numerical Recipes
!
        Integer idum, IA,IM,IQ,IR,NTAB,NDIV
        Real ranpp,AM,EPS,RNMX
        Parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,   &
             NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
        INTEGER j,k,iv1(NTAB),iy1

        if (idum.le.0 .or. iy1.eq.0) then
                idum=max(-idum,1)
                do j=NTAB+8,1,-1
                   k=idum/IQ
                   idum=IA*(idum-k*IQ)-IR*k
                   if (idum.lt.0) idum=idum+IM
                   if (j.le.NTAB) iv1(j)=idum
                enddo
                iy1=iv1(1)
        endif
        k=idum/IQ
        idum=IA*(idum-k*IQ)-IR*k
        if(idum.lt.0) idum=idum+IM
        j=1+iy1/NDIV
        iy1=iv1(j)
        iv1(j)=idum
        ranpp=min(AM*iy1,RNMX)

        return
        END

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

        Function ranff2(idum,iv,iy,idum2)
        INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
        REAL ranff2,AM,EPS,RNMX
        PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1)
        PARAMETER (IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211)
        PARAMETER (IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7)
        PARAMETER (RNMX=1.-EPS)
                                                                                                                                                                      
!Long period (> 2*10^18) random number generator of L'Ecuyer with Bays-Durham
!shuffle and added safeguards. Returns a uniform random deviate between 0.0 and 1.0
!(exclusive of the endpoint values). Call with idum a negative integer to initialize;
!thereafter, do not alter idum between successive deviates in a sequence.
!RNMX should approximate the largest floating value that is less than 1.
                                                                                                                                                                      
         INTEGER idum2,j,k,iv(NTAB),iy
                                                                                                                                                                      
         if (idum.le.0) then
            idum=max(-idum,1)
            idum2=idum
            do j=NTAB+8,1,-1
               k=idum/IQ1
               idum=IA1*(idum-k*IQ1)-k*IR1
               if (idum.lt.0) idum=idum+IM1
               if (j.le.NTAB) iv(j)=idum
            enddo
            iy=iv(1)
         endif
                                                                                                                                                                      
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        k=idum2/IQ2
        idum2=IA2*(idum2-k*IQ2)-k*IR2
        if (idum2.lt.0) idum2=idum2+IM2
        j=1+iy/NDIV
        iy=iv(j)-idum2
        iv(j)=idum
        if(iy.lt.1)iy=iy+IMM1
        ranff2=min(AM*iy,RNMX)
        return
        END
         
!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

        SUBROUTINE hpsort (arr, n)

        INTEGER :: n,M,NSTACK
        INTEGER :: arr(n)
        PARAMETER (M=7,NSTACK=50)
        INTEGER :: i,ir,j,jstack,k,l,istack(NSTACK)
        REAL a,temp

        jstack=0
        l=1
        ir=n

1       if(ir-l.lt.M) then

          do j=l+1,ir
            a=arr(j)
            do i=j-1,l,-1
              if(arr(i).le.a) goto 2
              arr(i+1)=arr(i)
            enddo
            i=l-1
2           arr(i+1)=a
          enddo
          if(jstack.eq.0) return
          ir=istack(jstack)
          l=istack(jstack-1)
          jstack=jstack-2

        else

          k=(l+ir)/2
          temp=arr(k)
          arr(k)=arr(l+1)
          arr(l+1)=temp

          if(arr(l).gt.arr(ir))then
            temp=arr(l)
            arr(l)=arr(ir)
            arr(ir)=temp
          endif

          if(arr(l+1).gt.arr(ir))then
             temp=arr(l+1)
             arr(l+1)=arr(ir)
             arr(ir)=temp
          endif

          if(arr(l).gt.arr(l+1))then
             temp=arr(l)
             arr(l)=arr(l+1)
             arr(l+1)=temp
          endif
          i=l+1
          j=ir
          a=arr(l+1)

3       continue

        i=i+1

        if(arr(i).lt.a)goto 3
4       continue
        j=j-1

        if(arr(j).gt.a)goto 4

        if(j.lt.i)goto 5

        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        jstack=jstack+2

        if(jstack.gt.NSTACK)pause
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif

        endif
        goto 1

        END

!HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
