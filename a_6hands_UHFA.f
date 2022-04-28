c 6hands UHF --------------------------------------------------
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 JOBZ  
      character*1 UPLO
      parameter(nn=24,m=nn/2,t1=1.0d0,t2=1.0d0,u=2.0d0
     &,limW2=(3*nn-1)*100)
      dimension dm(2*nn),dm2(2*nn)
      dimension hk(nn,nn),deigenup(nn),deigendown(nn)
      dimension hhf(nn,nn)
      dimension fup(nn,m),fdown(nn,m)
      dimension gup(nn,nn),gdown(nn,nn)
      dimension sign(nn),WORK2(limW2)
c shift site numbers cccccccccccccccccccccccccccc
      dimension ihop1(nn,4),ihop2(nn,2)
      dimension lsitex(nn),lsitey(nn)
      dimension lsub(nn)
ccccccccccccccccccccccccccccccccccccccccccccccccc
      data JOBZ/'V'/,UPLO/'U'/,iWF/0/
      
      if(iWF.eq.1) then
      open(unit=11,file="UHFA_t2_1_u2_m.txt"
     &,status='unknown')
      end if
      open(unit=66,file="ENG_UHFA_t2_1_u2_m.txt"
     &,status='unknown')
      open(unit=77,file="DBLE_UHFA_t2_1_u2_m.txt"
     &,status='unknown')

      call shift(nn,ihop1,ihop2,lsitex,lsitey,lsub)
      

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
      
      hhf=hk
      LWORK2=3*nn-1
      call DSYEV(JOBZ,UPLO,nn,hhf,nn,deigenup,WORK2,LWORK2,INFO)
      LWORK2=int(WORK2(1))
      if(LWORK2.ge.limW2) then
         write(6,*) 'lack of WORK space',LWORK2,limW2
         stop
      end if
      write(6,*) 'WORK space',LWORK2,limW2
      
      call main(nn,m,t1,t2,u,dm,dm2,hk,deigenup,deigendown
     &,hhf,fup,fdown,gup,gdown,sign,WORK2,LWORK2,ihop1,ihop2
     &,lsitex,lsitey,lsub,iWF)
      
      stop
      end 
      
      
      
      
      
      
      subroutine main(nn,m,t1,t2,u,dm,dm2,hk,deigenup,deigendown
     &,hhf,fup,fdown,gup,gdown,sign,WORK2,LWORK2,ihop1,ihop2
     &,lsitex,lsitey,lsub,iWF)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 JOBZ  
      character*1 UPLO
      dimension dm(2*nn),dm2(2*nn)
      dimension hk(nn,nn),deigenup(nn),deigendown(nn)
      dimension hhf(nn,nn)
      dimension fup(nn,m),fdown(nn,m)
      dimension gup(nn,nn),gdown(nn,nn)
      dimension sign(nn),WORK2(LWORK2)
c shift site numbers cccccccccccccccccccccccccccc
      dimension ihop1(nn,4),ihop2(nn,2)
      dimension lsitex(nn),lsitey(nn)
      dimension lsub(nn)
ccccccccccccccccccccccccccccccccccccccccccccccccc
      data JOBZ/'V'/,UPLO/'U'/
      data ixran/562356/,limconv/100000/,limstep/50/
      data dlimdm/1.0d-5/,dlimeng/1.0d-7/
      energy2=100.0d0
      
      WORK2=0.0d0
      energy2=100.0d0

      m2=m*2
      nn2=nn*2
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
      do 90
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

      call DSYEV(JOBZ,UPLO,nn,hhf,nn,deigenup,WORK2,LWORK2,INFO)

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

      call DSYEV(JOBZ,UPLO,nn,hhf,nn,deigendown,WORK2,LWORK2,INFO)

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
            if(energy-energy2.le.dlimeng)then
               istep=istep+1
               write(66,*) istep,energy,ixran2
               write(77,*) istep,double,ixran2
               write(6,*) istep,energy,double,ixran2
               energy2=energy
      if(iWF.eq.1) then
                  write(11,*) 'step',istep,ixran2
                  write(11,*) 1
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
      end if
                  if(istep.eq.limstep) stop
               goto 90
           else
               goto 90
           end if
      end if


         dm2=dm

      if(istep0.eq.limconv) then
         write(6,*) 'do not converge'
         stop
      end if
      
113   continue
c end self consistent cccccccccccccccccccccccccccccccccccccccccccccccccccccc
90    continue
      
      
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc













      
      
      
      
      





