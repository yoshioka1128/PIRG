c 6hands read WF and calculate the Physical quantities -----------------
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 JOBZ  
      character*1 UPLO
      
      parameter(nn=32,m=nn/2,u=6.0d0,t1=1.0d0,t2=1.0d0
     &,iqpon=1,ispin=0,Lmax=500)

      dimension fai(nn,2*m,Lmax),faiengv(2*nn,2*m,Lmax)
c schmidt
      dimension aeigen(m),beigen(m),cc(m),fai2(nn,m*2)      
      dimension a(m,m),b(m,m),WORK2(3*Lmax-1)
c size depend !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
c shift site numbers cccccccccccccccccccccccccccc
      dimension ihop1(nn,4),ihop2(nn,2)
      dimension lsitex(nn),lsitey(nn)
      dimension lsub(nn)
c end size depend !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      pai=dacos(-1.0d0)

      ihop1=0.0d0
      ihop2=0.0d0
      lsitex=0
      lsitey=0
      
      call shift(nn,ihop1,ihop2,lsitex,lsitey,lsub)
      
       read(5,*) iLdim
      write(6,*) 'U/t=',u
      write(6,*) 't2/t1=',t2/t1
      write(6,*) 'dim. =',iLdim
c read H.F.A.result 
      do l=1,iLdim
         do inum=1,m
            do i=1,nn
               read(5,*) fai(i,inum,l)
            end do
         end do
      end do
      do l=1,iLdim
         do inum=1,m
            do i=1,nn
               read(5,*) fai(i,inum+m,l)
            end do
         end do
      end do
      
      if(iLdim.ne.1) then
      do l=1,iLdim
      call schmidt(nn,m,iLdim,fai,fai2,l,a,b,aeigen,beigen,cc)
      end do
      end if
      
      
      do l=1,iLdim
      do j=1,m
      do i=1,nn
         faiengv(i,j,l)=fai(i,j,l)
         faiengv(i+nn,j,l)=0.0d0
         faiengv(i,j+m,l)=0.0d0
         faiengv(i+nn,j+m,l)=fai(i,j+m,l)
      end do
      end do
      end do
      
      LWORK2=3*Lmax-1
      call engv(u,t1,t2,iLdim,ispin,iqpon
     &,faiengv,ihop1,ihop2,lsitex,lsitey,lsub,WORK2,LWORK2)
      
      stop
      end
c end PIRG checker board
c----------------------------------------------------------------
      
      
      
      
      
      
      
      






















   
      
