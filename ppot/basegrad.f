      subroutine basegrad(n,echo,max_elem,iz,coord,rscal,qscal,e,g)  
      implicit none
      integer n,max_elem,maxat,iz(*)
      real*8 coord(3,*),e,g(3,*)
      real*8 autokcal,rscal,qscal     
      real*8 fi,fj,ff,rf,r
c cut-off radii for all element pairs
      real*8 r0ab(max_elem,max_elem),autoang
      parameter (autoang =0.52917726d0)
      parameter (autokcal=627.509541d0)
      logical echo
      real*8 r0
      integer i,j

      g(1:3,1:n)=0
      e=0

      if (echo) then
      write(*,*)'---------------------------------------------------'        
      write(*,*)'                   SRB correction'        
      write(*,*)'---------------------------------------------------'        
      end if
     
      call setr0ab(max_elem,autoang,r0ab)

c paramter of the method are rscal and qscal
c     rscal=0.7
c     qscal=0.03

      do i=1,n-1
         do j=i+1,n
            if(iz(i).lt.1.or.iz(i).gt.18) cycle
            if(iz(j).lt.1.or.iz(j).gt.18) cycle
            r=sqrt( (coord(1,i)-coord(1,j))**2
     .             +(coord(2,i)-coord(2,j))**2
     .             +(coord(3,i)-coord(3,j))**2)
            r0=rscal*r0ab(iz(i),iz(j))**0.75d0   
            fi=float(iz(i))
            fj=float(iz(j))
            ff=-(fi*fj)**1.5d0
c           write(*,*) i,j,-qscal*0.1*fj*fj,r0
c           if(iz(i).gt.10)fi=float(iz(i))-10.
c           if(iz(j).gt.10)fj=float(iz(j))-10.
c           e=e-(fi*fj)/(r0+r**3)
            e=e+ff*exp(-r0*r)
            rf=qscal/r
            g(1,i)=g(1,i)-ff*r0*(coord(1,i)-coord(1,j))*exp(-r0*r)*rf
            g(1,j)=g(1,j)-ff*r0*(coord(1,j)-coord(1,i))*exp(-r0*r)*rf
            g(2,i)=g(2,i)-ff*r0*(coord(2,i)-coord(2,j))*exp(-r0*r)*rf
            g(2,j)=g(2,j)-ff*r0*(coord(2,j)-coord(2,i))*exp(-r0*r)*rf
            g(3,i)=g(3,i)-ff*r0*(coord(3,i)-coord(3,j))*exp(-r0*r)*rf
            g(3,j)=g(3,j)-ff*r0*(coord(3,j)-coord(3,i))*exp(-r0*r)*rf
         enddo
      enddo   

      e=e*qscal

      if (echo) then

      write(*,*)'parameters:'
      write(*,'(''  s     :'',f10.3)') qscal
      write(*,'(''  gamma :'',f10.3)') rscal
      write(*,*)''
      write(*,'('' E_SRB /kcal,au:'',f11.4,f12.8)') e*autokcal,e
      write(*,*)'---------------------------------------------------'
      end if

      end
