c 6hands UHF --------------------------------------------------
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      parameter(nx=50,ny=50,t1=1.0d0,t2=1.0d0)
      dimension eng(nx*ny)
      
      pai=dacos(-1.0d0)
      istep=0
      do ky=-ny/2+1,ny/2
         do kx=-nx/2+1,nx/2
         istep=istep+1
         akx=pai*dble(kx)*2.0d0/dble(nx)
         aky=pai*dble(ky)*2.0d0/dble(ny)
c      write(6,10) dble(kx)*2.0d0/dble(nx),dble(ky)*2.0d0/dble(ny),
c     &-2.0d0*t1*(dcos(akx)+dcos(aky))
c     &-2.0d0*t2*(dcos(akx+aky))
      eng(istep)=-2.0d0*t1*(dcos(akx)+dcos(aky))
     &-2.0d0*t2*(dcos(akx+aky))
         end do
      end do
      
      do made=nx*ny-1,1,-1    
         do i=1,made
            if(eng(i).gt.eng(i+1)) then
               w=eng(i)
               eng(i)=eng(i+1)
               eng(i+1)=w
            end if
         end do
      end do
      
      do i=1,nx*ny
         write(6,*) dble(i)/dble(nx*ny),eng(i)
      end do
      
10    format(F15.7,F15.7,F15.7)      
      stop
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc













      
      
      
      
      





