c 6hands initial WF --------------------------------------------------
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 JOBZ  
      character*1 UPLO
      parameter(nn=32,m=nn/2,t1=1.0d0,t2=1.0d0,u=6.0d0)
      dimension h2(2,2)
      dimension hk1(nn,nn,2),dm(2*nn),dm2(2*nn)
      dimension hk(nn,nn),deigenup(nn),deigendown(nn)
      dimension hhf(nn,nn)
      dimension fup(nn,m),fdown(nn,m),IPIV(m)
      dimension gup(nn,nn),gdown(nn,nn),WORK(m)
      dimension sqreal(nn)
      dimension sign(nn),WORK2(3*nn-1)
c shift site numbers cccccccccccccccccccccccccccc
      dimension ihop1(nn,4),ihop2(nn,2)
      dimension lsitex(nn),lsitey(nn)
      dimension lsub(nn)
ccccccccccccccccccccccccccccccccccccccccccccccccc
      data ITYPE/1/,JOBZ/'V'/,UPLO/'U'/
      data ixran/499031049/,limconv/100000/,limstep/50/
      data dlimdm/1.0d-5/
      
      open(unit=11,file="UHFA_t2_1_u6_m.txt"
     &,status='unknown')
      pai=dacos(-1.0d0)

      m2=m*2
      nn2=nn*2
      LWORK=3*nn-1
      dm2=0.0d0
      
      call shift(nn,ihop1,ihop2,lsitex,lsitey,lsub)
      
      write(6,*) 'U/t=',u
      write(6,*) 't2/t1=',t2/t1


c Matrix of kinetic energy ccccccccccccccccccccccccccccccc
      hk=0.0d0
      do i=1,nn
         do j=1,2
            hk(i,ihop1(i,j))=hk(i,ihop1(i,j))-t1
         end do
         
         do j=3,4
            hk(i,ihop1(i,j))=hk(i,ihop1(i,j))-t1
         end do


         j=1
            hk(i,ihop2(i,j))=hk(i,ihop2(i,j))-t2
         
         j=2
            hk(i,ihop2(i,j))=hk(i,ihop2(i,j))-t2

      end do
cccccccccccccccccccccccccccccccccccccccccccccccccccc


       
      istep=0
c make randum field cccccccccccccccccccccccccccccccccccccccccccc
      ixran2=ixran

c reset run2
      r=ran2(-ixran)
c end reset run2

      do i=1,nn
         sign(i)=-1.0d0
      end do
      
      do
      isite=int(dble(nn)*ran2(ixran))+1
      sign(isite)=1.0d0
         icount=0
         do j=1,nn
            if(sign(j).gt.0) icount=icount+1
            if(icount.eq.nn/2) goto 100
         end do
      end do
      
100   do i=1,nn
         dmag=ran2(ixran)*sign(i)
         dm(i)=(1.0d0+dmag)/2.0d0
         dm(i+nn)=(1.0d0-dmag)/2.0d0
c      write(6,*) sign(i)
      end do
c      stop
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




c self consistent cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      istep0=0
      do 113
      istep0=istep0+1
c for up spin state ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      hhf=hk
      do i=1,nn
         hhf(i,i)=hhf(i,i)+u*dm(i+nn)
      end do

      call DSYEV(JOBZ,UPLO,nn,hhf,nn,deigenup,WORK2,LWORK,INFO)

      do inum=1,m
         do i=1,nn
            fup(i,inum)=hhf(i,inum)
         end do
      end do
      
c make green fnc.ccccccccccccccccccccccccccccccccccccccccccccccccc
c make g matrix
      do j=1,nn
      do i=1,nn
         gup(i,j)=0.0d0
         do k=1,m
         gup(i,j)=gup(i,j)+fup(i,k)*fup(j,k)
         end do
      end do
      end do
c end make gmatrix
c end make green fnc.ccccccccccccccccccccccccccccccccccccccccccccc      

      do i=1,nn
         dm(i)=gup(i,i)
      end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




c for down spin state ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      hhf=hk
      do i=1,nn
         hhf(i,i)=hhf(i,i)+u*dm(i)
      end do

      call DSYEV(JOBZ,UPLO,nn,hhf,nn,deigendown,WORK2,LWORK,INFO)

      do inum=1,m
         do i=1,nn
            fdown(i,inum)=hhf(i,inum)
         end do
      end do
      
c make green fnc.ccccccccccccccccccccccccccccccccccccccccccccccccc
c make g matrix
      do j=1,nn
      do i=1,nn
         gdown(i,j)=0.0d0
         do k=1,m
         gdown(i,j)=gdown(i,j)+fdown(i,k)*fdown(j,k)
         end do
      end do
      end do
c end make gmatrix
c end make green fnc.ccccccccccccccccccccccccccccccccccccccccccccc      

      do i=1,nn
         dm(i+nn)=gdown(i,i)
      end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc










      icount=0
      do i=1,nn2
         if(abs(dm(i)-dm2(i)).le.dlimdm) icount=icount+1
      end do

      if(icount.eq.nn2) then
         energy=0.0d0
         do i=1,m
            energy=energy+deigenup(i)+deigendown(i)
         end do
         
         double=0.0d0
         do i=1,nn
            double=double+gup(i,i)*gdown(i,i)
            energy=energy-u*gup(i,i)*gdown(i,i)
         end do
         double=double/dble(nn)
               istep=istep+1
               write(6,*) istep,energy,double,ixran2
               energy2=energy
                  write(11,*) 1,ixran2
                  do inum=1,m
                     do i=1,nn
                        write(11,*) fup(i,inum)
                     end do
                  end do
                  do inum=1,m
                     do i=1,nn
                        write(11,*) fdown(i,inum)
                     end do
                  end do
                  if(istep.eq.limstep) stop
                  goto 90
      end if


         dm2=dm

      if(istep0.eq.limconv) then
         write(6,*) 'do not converge'
         stop
      end if
      
113   continue
c end self consistent cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      

90    stop
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc













      
      
      
      
      





