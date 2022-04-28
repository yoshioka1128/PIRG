c 6hands U0 energy (shift check) ----------------------------------------------
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 JOBZ  
      character*1 UPLO
      parameter(nn=16,m=nn/2,t1=1.0d0,t2=1.0d0)
      dimension hk(nn,nn),deigen(nn)
      dimension WORK2(3*nn-1)
c shift site numbers cccccccccccccccccccccccccccc
      dimension ihop1(nn,4),ihop2(nn,2)
      dimension lsitex(nn),lsitey(nn)
      dimension lsub(nn)
ccccccccccccccccccccccccccccccccccccccccccccccccc
      data ITYPE/1/,JOBZ/'V'/,UPLO/'U'/
      data limconv/100000/,limstep/50/
      data dlimdm/1.0d-5/
      
      open(unit=11,file="UHFA_R_t2_1_u8_m.txt"
     &,status='unknown')
      pai=dacos(-1.0d0)

      m2=m*2
      nn2=nn*2
      LWORK=3*nn-1
      dm2=0.0d0
      
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
      LWORK=3*nn-1
      call DSYEV(JOBZ,UPLO,nn,hk,nn,deigen,WORK2,LWORK,INFO)
      eng=0.0d0
      do i=1,nn/2
         eng=eng+deigen(i)
      end do
      
      write(6,*) 'U=0 energy',eng*2.0d0
       
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc













      
      
      
      
      





