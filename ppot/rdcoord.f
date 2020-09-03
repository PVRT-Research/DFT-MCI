CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c read coordinates

      subroutine rd(fname,n,xyz,iat)
      implicit real*8 (a-h,o-z)
      dimension xyz(3,n),iat(n),xx(10)
      character*128 line
      character*(*) fname

      ich=142
      open(unit=ich,file=fname)
      read(ich,'(a)')line
      do i=1,n
         read(ich,*) xyz(1:3,i)
      enddo

      rewind ich
      read(ich,'(a)')line
      do i=1,n
         read(ich,'(a)')line
         call elem(line,j)
         iat(i)=j
      enddo

      close(ich)
      end

      subroutine rd0(fname,n)
      implicit real*8 (a-h,o-z)
      dimension xx(10)
      character*128 line
      character*(*) fname

      ich=142
      open(unit=ich,file=fname)
      n=0
 100  read(ich,'(a)',end=200)line
         if(index(line,'$redu').ne.0)goto 200
         if(index(line,'$user').ne.0)goto 200
         call readl(line,xx,nn)
         if(nn.ne.3) goto 100
         n=n+1
      goto 100
 200  continue

      close(ich)
      return
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE ELEM(KEY1, NAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*(*) KEY1
      CHARACTER*2 ELEMNT(107),E

      DATA ELEMNT/'h ','he',
     1 'li','be','b ','c ','n ','o ','f ','ne',
     2 'na','mg','al','si','p ','s ','cl','ar',
     3 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',
     4 'zn','ga','ge','as','se','br','kr',
     5 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',
     6 'cd','in','sn','sb','te','i ','xe',
     7 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     8 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',
     9 'au','hg','tl','pb','bi','po','at','rn',
     1 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx',
     2 'fm','md','cb','xx','xx','xx','xx','xx'/

      nat=0
      e='  '
      k=1
      DO J=1,len(key1)
         if (k.gt.2)exit
         N=ICHAR(key1(J:J))
         if(n.ge.ichar('A') .and. n.le.ichar('Z') )then
            e(k:k)=char(n+ICHAR('a')-ICHAR('A'))
            k=k+1
         endif
         if(n.ge.ichar('a') .and. n.le.ichar('z') )then
            e(k:k)=key1(j:j)
            k=k+1
         endif
      enddo

      DO I=1,107
         if(e.eq.elemnt(i))then
            NAT=I
            RETURN
         ENDIF
      ENDDO

      end

C     *****************************************************************         

      FUNCTION ASYM(I)
      CHARACTER*2 ASYM
      CHARACTER*2 ELEMNT(107), AS
      DATA ELEMNT/'h ','he',
     1 'li','be','b ','c ','n ','o ','f ','ne',
     2 'na','mg','al','si','p ','s ','cl','ar',
     3 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',
     4 'zn','ga','ge','as','se','br','kr',
     5 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',
     6 'cd','in','sn','sb','te','i ','xe',
     7 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',
     8 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',
     9 'au','hg','tl','pb','bi','po','at','rn',
     1 'fr','ra','ac','th','pa','u ','np','pu','am','cm','bk','cf','xx',
     2 'fm','md','cb','xx','xx','xx','xx','xx'/
      AS=ELEMNT(I)
      CALL UPPER(AS)
      ASYM=AS
      if(i.eq.103) asym='XX'
      RETURN
      END

      SUBROUTINE UPPER(AS)
      CHARACTER*2 AS
      NSP=ICHAR(' ')
      ND=ICHAR('A')-ICHAR('a')
      DO 10 I=1,2
         J=ICHAR(AS(I:I))
         IF(J.NE.NSP)J=J+ND
         AS(I:I)=CHAR(J)
  10  CONTINUE
      RETURN
      END

