      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 JOBZ  
      character*1 UPLO
      parameter(nn=30,m2=nn,t1=1.0d0,t2=1.0d0,u=13.0d0,dmagint=0.5d0)
      dimension h2(2,2)
      dimension hk1(nn,nn,2),dm(2*nn),dm2(2*nn),deigen(2*nn)
      dimension hk(2*nn,2*nn)
      dimension hhf(2*nn,2*nn)
      dimension f(2*nn,m2),f0(2*nn,m2)
      dimension g(2*nn,2*nn)
      dimension sqreal(nn)
      dimension sign(nn),WORK2(3*2*nn-1)
c shift site numbers cccccccccccccccccccccccccccc
      dimension ihop1(nn,4),ihop2(nn,2)
      dimension lsitex(nn),lsitey(nn)
      dimension lsitenum(10,10)
      dimension lsub(nn)
      dimension ip(10),im(10)
      dimension ldx(nn),ldy(nn)
ccccccccccccccccccccccccccccccccccccccccccccccccc
      data ITYPE/1/,JOBZ/'V'/,UPLO/'U'/
      data ixran/934757/
      
      open(unit=11,file="120_HFA_t2_1_u13_m.txt"
     &,status='unknown')
      open(unit=66,file="ENG_120_HFA_t2_1_u13_m.txt"
     &,status='unknown')
      open(unit=77,file="DBLE_120_HFA_t2_1_u13_m.txt"
     &,status='unknown')
c      open(unit=16,file="Real_120_Mag_t2_1_u13_m.txt"
c     &,status='unknown')

      energy2=100.0d0
      pai=dacos(-1.0d0)

      nn2=nn*2
      
      
      call shift(nn,ihop1,ihop2,lsitex,lsitey,lsub)

      
c Matrix of kinetic energy ccccccccccccccccccccccccccccccc
      hk=0.0d0
      do i=1,nn
         do j=1,2
      hk(i,ihop1(i,j))=hk(i,ihop1(i,j))
     &-t1*dcos(2.0d0*pai/3.0d0)
      hk(i,ihop1(i,j)+nn)=hk(i,ihop1(i,j)+nn)
     &+t1*dsin(2.0d0*pai/3.0d0)
      hk(i+nn,ihop1(i,j))=hk(i+nn,ihop1(i,j))
     &-t1*dsin(2.0d0*pai/3.0d0)
      hk(i+nn,ihop1(i,j)+nn)=hk(i+nn,ihop1(i,j)+nn)
     &-t1*dcos(2.0d0*pai/3.0d0)
         end do
         
         do j=3,4
      hk(i,ihop1(i,j))=hk(i,ihop1(i,j))
     &-t1*dcos(-2.0d0*pai/3.0d0)
      hk(i,ihop1(i,j)+nn)=hk(i,ihop1(i,j)+nn)
     &+t1*dsin(-2.0d0*pai/3.0d0)
      hk(i+nn,ihop1(i,j))=hk(i+nn,ihop1(i,j))
     &-t1*dsin(-2.0d0*pai/3.0d0)
      hk(i+nn,ihop1(i,j)+nn)=hk(i+nn,ihop1(i,j)+nn)
     &-t1*dcos(-2.0d0*pai/3.0d0)
         end do


         j=1
      hk(i,ihop2(i,j))=hk(i,ihop2(i,j))
     &-t2*dcos(4.0d0*pai/3.0d0)
      hk(i,ihop2(i,j)+nn)=hk(i,ihop2(i,j)+nn)
     &+t2*dsin(4.0d0*pai/3.0d0)
      hk(i+nn,ihop2(i,j))=hk(i+nn,ihop2(i,j))
     &-t2*dsin(4.0d0*pai/3.0d0)
      hk(i+nn,ihop2(i,j)+nn)=hk(i+nn,ihop2(i,j)+nn)
     &-t2*dcos(4.0d0*pai/3.0d0)
         
         j=2
      hk(i,ihop2(i,j))=hk(i,ihop2(i,j))
     &-t2*dcos(-4.0d0*pai/3.0d0)
      hk(i,ihop2(i,j)+nn)=hk(i,ihop2(i,j)+nn)
     &+t2*dsin(-4.0d0*pai/3.0d0)
      hk(i+nn,ihop2(i,j))=hk(i+nn,ihop2(i,j))
     &-t2*dsin(-4.0d0*pai/3.0d0)
      hk(i+nn,ihop2(i,j)+nn)=hk(i+nn,ihop2(i,j)+nn)
     &-t2*dcos(-4.0d0*pai/3.0d0)

      end do
cccccccccccccccccccccccccccccccccccccccccccccccccccc


       
c make 120 degree field cccccccccccccccccccccccccccccccccccccccccccc
      do i=1,nn
         n=lsitex(i)+lsitey(i)
         if(mod(n,2).eq.0) then
            sign(i)=1.0d0
         else
            sign(i)=1.0d0
         end if
      end do
      
      do i=1,nn
         dmag=dmagint*sign(i)
         dm(i)=(1.0d0+dmag)/2.0d0
         dm(i+nn)=(1.0d0-dmag)/2.0d0
      end do
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      write(6,*) 
      write(6,*) 'U/t=',u
      write(6,*) 't2/t1=',t2/t1


      istep2=0
c self consistent cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      do 113
      istep2=istep2+1
c for up spin state ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      hhf=hk
      do 8 i=1,nn
         hhf(i,i)=hhf(i,i)+u*dm(i+nn)
         hhf(i+nn,i+nn)=hhf(i+nn,i+nn)+u*dm(i)
8     continue

      LWORK=3*nn2-1
      call DSYEV(JOBZ,UPLO,nn2,hhf,nn2,deigen,WORK2,LWORK,INFO)

      do inum=1,m2
         do i=1,nn2
            f(i,inum)=hhf(i,inum)
         end do
      end do
      

c make green fnc.ccccccccccccccccccccccccccccccccccccccccccccccccc
c make g matrix
      do j=1,nn2
         do i=1,nn2
         g(i,j)=0.0d0
            do k=1,m2
            g(i,j)=g(i,j)+f(i,k)*f(j,k)
            end do
         end do
      end do
c end make gmatrix
c end make green fnc.ccccccccccccccccccccccccccccccccccccccccccccc      

      do i=1,nn2
         dm(i)=g(i,i)
      end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc














      icount=0
      do 151 i=1,nn2
         if(abs(dm(i)-dm2(i)).le.1.0d-10) icount=icount+1
151   continue
      write(6,*) istep2,dm(1)
c      write(6,*) nn2-icount

      if(icount.eq.nn2) then
         energy=0.0d0
         do i=1,m2
            energy=energy+deigen(i)
         end do
         
         double=0.0d0
         do i=1,nn
            double=double+g(i,i)*g(i+nn,i+nn)
c     &                   -g(i+nn,i)*g(i,i+nn)
            energy=energy-u*g(i,i)*g(i+nn,i+nn)
c            write(16,*) lsitex(i),lsitey(i),g(i,i)-g(i+nn,i+nn)
         end do
         double=double/dble(nn)
            if(energy-energy2.le.1.0d-10)then
               write(66,*) energy
               write(77,*) double
               write(6,*) 'eng=',energy,double,u
               energy2=energy
                  write(11,*) 1
                  f0=f
                  f=0.0d0
                  do j=1,m2
                     do i=1,nn
                        f(i,j)=f(i,j)
     &                 +dcos((lsub(i)-1)*2.0d0*pai/3.0d0)*f0(i,j)
                        f(i,j)=f(i,j)
     &                 -dsin((lsub(i)-1)*2.0d0*pai/3.0d0)*f0(i+nn,j)
                        f(i+nn,j)=f(i+nn,j)
     &                 +dsin((lsub(i)-1)*2.0d0*pai/3.0d0)*f0(i,j)
                        f(i+nn,j)=f(i+nn,j)
     &                 +dcos((lsub(i)-1)*2.0d0*pai/3.0d0)*f0(i+nn,j)
c                        write(6,*) 'sublattice',i,lsub(i)-1
                     end do
                  end do
                  do inum=1,m2
                     do i=1,nn2
                        write(11,*) f(i,inum)
                     end do
                  end do
                  
c make green fnc.ccccccccccccccccccccccccccccccccccccccccccccccccc
c make g matrix
      do j=1,nn2
      do i=1,nn2
         g(i,j)=0.0d0
         do k=1,m2
         g(i,j)=g(i,j)+f(i,k)*f(j,k)
         end do
      end do
      end do
c end make gmatrix
c end make green fnc.ccccccccccccccccccccccccccccccccccccccccccccc      
               do i=1,nn
c            write(16,*) lsitex(i),lsitey(i),g(i,i)-g(i+nn,i+nn)
               end do
               goto 90
           else
               goto 90
           end if
      end if


         dm2=dm


113   continue
c end self consistent cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c90    continue
      
      

90    stop
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc










