c 6hands PIRG fix L change U/t process ------------------------------------
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 JOBZ  
      character*1 UPLO
      
      parameter(nn=32,m=nn/2,u=8.0d0,t1=1.0d0,t2=1.0d0
     &,Lmax=10,limW=m*100,limW2=(Lmax*3-1)*100
c set up -----------------------------------------------------
     &,iRGtype=0,ichmod=1,ratelim1=0.1d0,ratelim2=0.05d0
     &,delu=0.5d0/u,limit=100000
c engv -----------------------------------------------------
     &,iengv=1,iqpon=1,ispin=0,iWF=0)

c RG type iRGtype=0: independent RG
c RG type iRGtype=1: depend RG

c change RG mode ichmod=0: no chanege
c change RG mode ichmod=1: change

c iengv=0: no engv process
c iengv=1: engv process


      dimension gup(nn,nn,Lmax),gdown(nn,nn,Lmax)
      dimension gup2(nn,nn,Lmax,2),gdown2(nn,nn,Lmax,2)
      dimension ginup(Lmax),gindown(Lmax)
      dimension ginup2(Lmax,2),gindown2(Lmax,2)
      
      dimension h(Lmax,Lmax),hin(Lmax,Lmax)
      dimension hold(Lmax,Lmax),hinold(Lmax,Lmax)
      dimension hsub(Lmax,2),hinsub(Lmax,2)
      
      dimension fai(nn,2*m,Lmax),faiengv(2*nn,2*m,Lmax)
      dimension faiold(nn,2*m)
      
      dimension WORK2(limW2),hgeed(Lmax,Lmax),hingeed(Lmax,Lmax)
      dimension deigen(Lmax),IWORK(Lmax)
      
      dimension lurg(Lmax,Lmax),iss(2)
      dimension gsubu(nn,nn),gsubd(nn,nn)
      dimension bufu(m,nn),bufd(m,nn)
      dimension a(m,m),b(m,m),IPIV(m),WORK(limW)
      dimension alpha1(-1:1)
      dimension pirgmp(2),isubL(Lmax)
      dimension isite(nn),idim0(Lmax)
c schmidt
      dimension stfai(nn,2*m,Lmax),stfaiold(nn,2*m)
      dimension aeigen(m),beigen(m),cc(m),fai2(nn,m*2)
c kinetic2      
      dimension hk(nn,nn),hku(nn,nn)
      dimension ee(nn),hnqei(nn),hkL0(nn,nn)
      dimension huup(nn),hudown(nn)
      dimension huup2(nn,nn),hudown2(nn,nn)
      
      
      
c size depend !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
c shift site numbers cccccccccccccccccccccccccccc
      dimension ihop1(nn,4),ihop2(nn,2)
      dimension lsitex(nn),lsitey(nn)
      dimension lsub(nn)
c end size depend !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      






      data ixran/9862/
      data ITYPE/1/,UPLO/'U'/
      
      ixran=ixran*int(u)
      r=ran2(-ixran)
      
      open(unit=80,file="eng_list_u6_m.txt",status='unknown')
      
      
      illup=0
      illdown=0
      
      
      pai=dacos(-1.0d0)
      engold=100.0d0
      
      huup=0.0d0
      hudown=0.0d0
      huup2=0.0d0
      hudown2=0.0d0
      alpha1=0.0d0
      ginup2=0.0d0
      gindown2=0.0d0
      gup2=0.0d0
      gdown2=0.0d0
      
      fai2=0.0d0
      fai=0.0d0

      h=0.0d0
      hin=0.0d0
      hsub=0.0d0
      hinsub=0.0d0
      


      ihop1=0.0d0
      ihop2=0.0d0
      lsitex=0
      lsitey=0
      
      call shift(nn,ihop1,ihop2,lsitex,lsitey,lsub)
      
      iss(1)=1
      iss(2)=-1

c make Urg dim. ---------------------------------------
      do lr=1,Lmax
         idif=0
         do ll=1,Lmax
            if(ll.eq.lr) then
               idif=1
            end if
               lurg(ll,lr)=ll+idif
         end do
      end do
c end make Urg dim. ----------------------------


      read(5,*) iLdim
      write(6,*) 'start_dim',iLdim
      if(iLdim.ne.Lmax) then
         write(6,*) 'dimension conflict'
         stop
      end if
c read H.F.A.result 
      do l=1,iLdim
         do inum=1,m
            do i=1,nn
               read(5,*) stfai(i,inum,l)
            end do
         end do
      end do
      do l=1,iLdim
         do inum=1,m
            do i=1,nn
               read(5,*) stfai(i,inum+m,l)
            end do
         end do
      end do
c end read First W.F.

      fai=stfai

      do l=1,Lmax
      call schmidt(nn,m,Lmax,fai,fai2,l,a,b,aeigen,beigen,cc)
      end do
      
      LWORK=m
      do lr=1,Lmax
         do ll=1,lr
         call mkgreen(ll,lr,nn,m,Lmax,IPIV,WORK,LWORK
     &  ,bufu,bufd,a,b,fai,gup,gdown,ginup,gindown)
         call mkH(ll,lr,nn,u,t1,t2,Lmax
     &  ,ihop1,ihop2,h,hin,gup,gdown,ginup,gindown)
         end do      
      end do
      
      iLdim=Lmax
      LWORK2=3*iLdim-1
      call geed0(Lmax,iLdim,h,hin
     &,hgeed,hingeed,deigen,eigen,WORK2,LWORK2,'N',label)
      pirgeng=eigen
c end Expand dimension -----------------------------------------------

      write(6,*) 'START eng',pirgeng,Lmax
      
      if(ichmod.eq.0) then
      icoarse=0
      ratelim=ratelim2
      else
      icoarse=1
      ratelim=ratelim1
      end if
      
      write(80,*)
      write(80,*) 0,pirgeng,iLdim,iLdim

      if(WORK(1).ge.limW) then
         write(6,*) 'lack of WORK space',WORK(1),limW
         stop
      end if
      write(6,*) 'WORK space',WORK(1),limW
      LWORK=int(WORK(1))
      if(WORK2(1).ge.limW2) then
         write(6,*) 'lack of WORK2 space',WORK2(1),limW2
         stop
      end if
      write(6,*) 'WORK2 space',WORK2(1),limW2
      LWORK2=int(WORK2(1))
c -----------------------------------------------------------------


























c--------------------------------------------------------------------      
c--------------------------------------------------------------------      
c PIRG step --------------------------------------------------------- 
c MAIN      
      delta=delu
      call sht(u,delta,alpha1)
      call kinetic(nn,m,t1,t2,delta/2.0d0,hkL0,hku
     &            ,ihop1,ihop2,hnqei,ee)
     
      write(6,*) 
      write(6,*) 'U/t=',u
      write(6,*) 't2/t1=',t2/t1
      write(6,*) 'dim. =',iLdim
      write(6,*) 'delta=',delta
      
      if(iRGtype.eq.0) then
         write(6,*) 'independent RG'
      else
         write(6,*) 'depend RG'
      end if
      
      
      
      istep=0
c################################################################
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c################################################################
      do 36 
         istep=istep+1





c Global update in PIRG step ----------------------------------------------  
      if(icoarse.eq.1) then
c select L order ---------------------------------------------    
      idim0=0
         i=0
         do
         i=i+1
            ixran2=int(ran2(ixran)*dble(iLdim))+1
            icount=0
               do j=1,i
                  if(idim0(j).eq.ixran2) icount=icount+1
               end do
            if(icount.eq.0) then
               idim0(i)=ixran2
            else
               i=i-1
            end if
            if(i.eq.iLdim) goto 171
         end do
c end select L order ---------------------------------------------    
      
171   ipstep=0

      do 25 lite2=1,iLdim
      
      lite=idim0(lite2)
      
c select site order ---------------------------------------------    
      isite=0
      i=0
      do
      i=i+1
         ixran2=int(ran2(ixran)*dble(nn))+1
         icount=0
            do j=1,i
               if(isite(j).eq.ixran2) icount=icount+1
            end do
         if(icount.eq.0) then
            isite(i)=ixran2
         else
            i=i-1
         end if
      if(i.eq.nn) goto 780
      end do
c end select site order ------------------------------------------    

780   engold=pirgeng
      do j=1,2*m
         do i=1,nn
            faiold(i,j)=fai(i,j,lite)
         end do
      end do
      do j=1,2*m
         do i=1,nn
            stfaiold(i,j)=stfai(i,j,lite)
         end do
      end do
      do ll=1,lite
         hold(ll,lite)=h(ll,lite)
         hinold(ll,lite)=hin(ll,lite)
      end do
      do ll=lite,iLdim
         hold(lite,ll)=h(lite,ll)
         hinold(lite,ll)=hin(lite,ll)
      end do

c first Hk/2 term-----------------------------------------------      

c preserve updated data ---------------------------------------
      do j=1,m
         do i=1,nn
            fai(i,j,lite)=0.0d0
            fai(i,j+m,lite)=0.0d0
               do k=1,nn
                  fai(i,j,lite)=fai(i,j,lite)
     &           +hkL0(i,k)*stfai(k,j,lite)
                  fai(i,j+m,lite)=fai(i,j+m,lite)
     &           +hkL0(i,k)*stfai(k,j+m,lite)
               end do
         end do
      end do


      call schmidt(nn,m,Lmax,fai,fai2,lite,a,b,aeigen,beigen,cc)
      do ll=1,iLdim
      call mkgreen(ll,lite,nn,m,Lmax,IPIV,WORK,LWORK
     &            ,bufu,bufd,a,b,fai,gup,gdown,ginup,gindown)
      call mkH(ll,lite,nn,u,t1,t2,Lmax
     &        ,ihop1,ihop2,h,hin,gup,gdown,ginup,gindown)
      end do      
      do ll=lite,iLdim
         h(lite,ll)=h(ll,lite)
         hin(lite,ll)=hin(ll,lite)
      end do

c end preserve updated data ----------------------------------

      call geed0(Lmax,iLdim,h,hin
     &,hgeed,hingeed,deigen,eigen,WORK2,LWORK2,'N',label)
      pirgeng=eigen
c end first Hk/2 term-----------------------------------------------

c second Hu*Hk/2 term -----------------------------------------------
      do  ijmk=1,nn
      ijm=isite(ijmk)
      call Urg(nn,m,u,t1,t2,Lmax,iLdim,lite,ijm,alpha1,h,hin
     &,ihop1,ihop2,hsub,hinsub,hgeed,hingeed,WORK2,LWORK2,iss
     &,gup,gdown,ginup,gindown,gsubu,gsubd,pirgeng,lurg
     &,gup2,gdown2,ginup2,gindown2,pirgmp,huup,hudown,deigen
     &,ipstepu,icoarse)
      end do
         
      do i=1,nn
         do j=1,nn
            huup2(i,j)=0.0d0
            hudown2(i,j)=0.0d0
            do k=1,nn
               huup2(i,j)=huup2(i,j)
     &        +hkL0(i,k)*huup(k)*hkL0(k,j)
               hudown2(i,j)=hudown2(i,j)
     &        +hkL0(i,k)*hudown(k)*hkL0(k,j)
            end do
         end do
      end do
c preserve updated data ---------------------------------------
      do j=1,m
         do i=1,nn
            fai(i,j,lite)=0.0d0
            fai(i,j+m,lite)=0.0d0
               do k=1,nn
                  fai(i,j,lite)=fai(i,j,lite)
     &           +huup2(i,k)*stfai(k,j,lite)
                  fai(i,j+m,lite)=fai(i,j+m,lite)
     &           +hudown2(i,k)*stfai(k,j+m,lite)
               end do
         end do
      end do

      do j=1,m*2
         do i=1,nn
            stfai(i,j,lite)=fai(i,j,lite)
         end do
      end do

      call schmidt(nn,m,Lmax,fai,fai2,lite,a,b,aeigen,beigen,cc)
      do ll=1,iLdim
      call mkgreen(ll,lite,nn,m,Lmax,IPIV,WORK,LWORK
     &            ,bufu,bufd,a,b,fai,gup,gdown,ginup,gindown)
      call mkH(ll,lite,nn,u,t1,t2,Lmax
     &        ,ihop1,ihop2,h,hin,gup,gdown,ginup,gindown)
      end do      
      do ll=lite,iLdim
         h(lite,ll)=h(ll,lite)
         hin(lite,ll)=hin(ll,lite)
      end do
c end preserve updated data ----------------------------------
      call geed0(Lmax,iLdim,h,hin
     &,hgeed,hingeed,deigen,eigen,WORK2,LWORK2,'N',label)
      pirgeng=eigen

      if((engold-pirgeng).lt.0) then
      pirgeng=engold
      do j=1,m*2
         do i=1,nn
         fai(i,j,lite)=faiold(i,j)
         stfai(i,j,lite)=stfaiold(i,j)
         end do
      end do   
      do ll=1,lite
         h(ll,lite)=hold(ll,lite)
         hin(ll,lite)=hinold(ll,lite)
      end do
      do ll=lite,iLdim
         h(lite,ll)=hold(lite,ll)
         hin(lite,ll)=hinold(lite,ll)
      end do
      else
         ipstep=ipstep+1
      end if
c end second Hu*Hk/2 term -----------------------------------------------

cc-------------------------------------------------------------------------  
      
25    continue      
c end Global update in PIRG step ---------------------------------------------  
















c Local update in PIRG step -------------------------------------------------  
      else if(icoarse.eq.0) then


c --------------------------------------------------------------------
c first Hk term -------------------------------------------------
c select L order ---------------------------------------------    
      idim0=0
         i=0
         do
         i=i+1
            ixran2=int(ran2(ixran)*dble(iLdim))+1
            icount=0
               do j=1,i
                  if(idim0(j).eq.ixran2) icount=icount+1
               end do
            if(icount.eq.0) then
               idim0(i)=ixran2
            else
               i=i-1
            end if
            if(i.eq.iLdim) goto 199
         end do
c end select L order ---------------------------------------------    

199   ipstepk=0
      do 303 lite2=1,iLdim
      
      lite=idim0(lite2)
      
      engold=pirgeng
      do j=1,2*m
         do i=1,nn
            faiold(i,j)=fai(i,j,lite)
         end do
      end do
      do j=1,2*m
         do i=1,nn
            stfaiold(i,j)=stfai(i,j,lite)
         end do
      end do
      do ll=1,lite
         hold(ll,lite)=h(ll,lite)
         hinold(ll,lite)=hin(ll,lite)
      end do
      do ll=lite,iLdim
         hold(lite,ll)=h(lite,ll)
         hinold(lite,ll)=hin(lite,ll)
      end do
      
c preserve updated data ---------------------------------------
      do j=1,m
         do i=1,nn
            fai(i,j,lite)=0.0d0
            fai(i,j+m,lite)=0.0d0
               do k=1,nn
                  fai(i,j,lite)=fai(i,j,lite)
     &           +hkL0(i,k)*stfai(k,j,lite)
                  fai(i,j+m,lite)=fai(i,j+m,lite)
     &           +hkL0(i,k)*stfai(k,j+m,lite)
               end do
         end do
      end do

      do j=1,m*2
         do i=1,nn
            stfai(i,j,lite)=fai(i,j,lite)
         end do
      end do


      call schmidt(nn,m,Lmax,fai,fai2,lite,a,b,aeigen,beigen,cc)
      do ll=1,iLdim
      call mkgreen(ll,lite,nn,m,Lmax,IPIV,WORK,LWORK
     &            ,bufu,bufd,a,b,fai,gup,gdown,ginup,gindown)
      call mkH(ll,lite,nn,u,t1,t2,Lmax
     &        ,ihop1,ihop2,h,hin,gup,gdown,ginup,gindown)
      end do      
      do ll=lite,iLdim
         h(lite,ll)=h(ll,lite)
         hin(lite,ll)=hin(ll,lite)
      end do
c end preserve updated data ----------------------------------

      call geed0(Lmax,iLdim,h,hin
     &,hgeed,hingeed,deigen,eigen,WORK2,LWORK2,'N',label)
      pirgeng=eigen

      if((engold-pirgeng).lt.0) then
      pirgeng=engold
      do j=1,m*2
         do i=1,nn
         fai(i,j,lite)=faiold(i,j)
         end do
      end do   
      do j=1,m*2
         do i=1,nn
         stfai(i,j,lite)=stfaiold(i,j)
         end do
      end do   
      do ll=1,lite
         h(ll,lite)=hold(ll,lite)
         hin(ll,lite)=hinold(ll,lite)
      end do
      do ll=lite,iLdim
         h(lite,ll)=hold(lite,ll)
         hin(lite,ll)=hinold(lite,ll)
      end do
      else
         ipstepk=ipstepk+1
      end if

cc-------------------------------------------------------------------------  
      
303   continue
c end first Hk term -------------------------------------------      








c --------------------------------------------------------------------
c coulomb energy RG -------------------------------------------------
c select L order ---------------------------------------------    
      idim0=0
         i=0
         do
         i=i+1
            ixran2=int(ran2(ixran)*dble(iLdim))+1
            icount=0
               do j=1,i
                  if(idim0(j).eq.ixran2) icount=icount+1
               end do
            if(icount.eq.0) then
               idim0(i)=ixran2
            else
               i=i-1
            end if
            if(i.eq.iLdim) goto 17
         end do
c end select L order ---------------------------------------------    
      
17    ipstepu=0
      do 24 lite2=1,iLdim
      
      lite=idim0(lite2)
      
c select site order ---------------------------------------------    
      isite=0
      i=0
      do
      i=i+1
         ixran2=int(ran2(ixran)*dble(nn))+1
         icount=0
            do j=1,i
               if(isite(j).eq.ixran2) icount=icount+1
            end do
         if(icount.eq.0) then
            isite(i)=ixran2
         else
            i=i-1
         end if
      if(i.eq.nn) goto 815
      end do
c end select site order ------------------------------------------    

815   do ll=1,iLdim
      call mkgreen(ll,lite,nn,m,Lmax,IPIV,WORK,LWORK
     &            ,bufu,bufd,a,b,fai,gup,gdown,ginup,gindown)
      end do
      
      engold=pirgeng
      
c Hu term -----------------------------------------------
      do  ijmk=1,nn
      ijm=isite(ijmk)
      call Urg(nn,m,u,t1,t2,Lmax,iLdim,lite,ijm,alpha1,h,hin
     &,ihop1,ihop2,hsub,hinsub,hgeed,hingeed,WORK2,LWORK2,iss
     &,gup,gdown,ginup,gindown,gsubu,gsubd,pirgeng,lurg
     &,gup2,gdown2,ginup2,gindown2,pirgmp,huup,hudown,deigen
     &,ipstepu,icoarse)
      end do
         
c      huup2=0.0d0
c      hudown2=0.0d0
c      do i=1,nn
c         huup2(i,i)=huup(i)
c         hudown2(i,i)=hudown(i)
c      end do
c preserve updated data ---------------------------------------
      do j=1,m
         do i=1,nn
            fai(i,j,lite)=0.0d0
            fai(i,j+m,lite)=0.0d0
                  fai(i,j,lite)=fai(i,j,lite)
     &           +huup(i)*stfai(i,j,lite)
                  fai(i,j+m,lite)=fai(i,j+m,lite)
     &           +hudown(i)*stfai(i,j+m,lite)
         end do
      end do

      do j=1,m*2
         do i=1,nn
            stfai(i,j,lite)=fai(i,j,lite)
         end do
      end do

      call schmidt(nn,m,Lmax,fai,fai2,lite,a,b,aeigen,beigen,cc)
      do ll=1,iLdim
      call mkgreen(ll,lite,nn,m,Lmax,IPIV,WORK,LWORK
     &            ,bufu,bufd,a,b,fai,gup,gdown,ginup,gindown)
      call mkH(ll,lite,nn,u,t1,t2,Lmax
     &        ,ihop1,ihop2,h,hin,gup,gdown,ginup,gindown)
      end do      
      do ll=lite,iLdim
         h(lite,ll)=h(ll,lite)
         hin(lite,ll)=hin(ll,lite)
      end do
c end preserve updated data ----------------------------------
      call geed0(Lmax,iLdim,h,hin
     &,hgeed,hingeed,deigen,eigen,WORK2,LWORK2,'N',label)
      pirgeng=eigen
      
cc-------------------------------------------------------------------------  
      
24    continue
c end coulomb energy RG -------------------------------------------------
c------------------------------------------------------------------------











c --------------------------------------------------------------------
c second Hk term -------------------------------------------------
c select L order ---------------------------------------------    
      idim0=0
         i=0
         do
         i=i+1
            ixran2=int(ran2(ixran)*dble(iLdim))+1
            icount=0
               do j=1,i
                  if(idim0(j).eq.ixran2) icount=icount+1
               end do
            if(icount.eq.0) then
               idim0(i)=ixran2
            else
               i=i-1
            end if
            if(i.eq.iLdim) goto 19
         end do
c end select L order ---------------------------------------------    

19    ipstepk=0
      do 30 lite2=1,iLdim
      
      lite=idim0(lite2)
      
      engold=pirgeng
      do j=1,2*m
         do i=1,nn
            faiold(i,j)=fai(i,j,lite)
         end do
      end do
      do j=1,2*m
         do i=1,nn
            stfaiold(i,j)=stfai(i,j,lite)
         end do
      end do
      do ll=1,lite
         hold(ll,lite)=h(ll,lite)
         hinold(ll,lite)=hin(ll,lite)
      end do
      do ll=lite,iLdim
         hold(lite,ll)=h(lite,ll)
         hinold(lite,ll)=hin(lite,ll)
      end do
      
c preserve updated data ---------------------------------------
      do j=1,m
         do i=1,nn
            fai(i,j,lite)=0.0d0
            fai(i,j+m,lite)=0.0d0
               do k=1,nn
                  fai(i,j,lite)=fai(i,j,lite)
     &           +hkL0(i,k)*stfai(k,j,lite)
                  fai(i,j+m,lite)=fai(i,j+m,lite)
     &           +hkL0(i,k)*stfai(k,j+m,lite)
               end do
         end do
      end do

      do j=1,m*2
         do i=1,nn
            stfai(i,j,lite)=fai(i,j,lite)
         end do
      end do

      call schmidt(nn,m,Lmax,fai,fai2,lite,a,b,aeigen,beigen,cc)
      do ll=1,iLdim
      call mkgreen(ll,lite,nn,m,Lmax,IPIV,WORK,LWORK
     &            ,bufu,bufd,a,b,fai,gup,gdown,ginup,gindown)
      call mkH(ll,lite,nn,u,t1,t2,Lmax
     &        ,ihop1,ihop2,h,hin,gup,gdown,ginup,gindown)
      end do      
      do ll=lite,iLdim
         h(lite,ll)=h(ll,lite)
         hin(lite,ll)=hin(ll,lite)
      end do
c end preserve updated data ----------------------------------

      call geed0(Lmax,iLdim,h,hin
     &,hgeed,hingeed,deigen,eigen,WORK2,LWORK2,'N',label)
      pirgeng=eigen

      if((engold-pirgeng).lt.0) then
      pirgeng=engold
      do j=1,m*2
         do i=1,nn
            fai(i,j,lite)=faiold(i,j)
         end do
      end do   
      do j=1,m*2
         do i=1,nn
            stfai(i,j,lite)=stfaiold(i,j)
         end do
      end do   
      do ll=1,lite
         h(ll,lite)=hold(ll,lite)
         hin(ll,lite)=hinold(ll,lite)
      end do
      do ll=lite,iLdim
         h(lite,ll)=hold(lite,ll)
         hin(lite,ll)=hinold(lite,ll)
      end do
      else
         ipstepk=ipstepk+1
      end if
      
      
c end third Hk term -------------------------------------------      
30    continue
c Local update in PIRG step ------------------------------------------  
      end if







      
      
c      ipstep=int((dble(ipstepk)+dble(ipstepu)/dble(nn))/2.0d0)
      if(icoarse.eq.0) then
      ipstep=int(dble(ipstepu)/dble(nn))
      end if
      
      write(6,*) istep,pirgeng,ipstep,iLdim
      write(80,*) istep,pirgeng,ipstep,iLdim
      
      


      if(ipstep.le.int(dble(iLdim)*ratelim)
     &.or.istep.eq.limit) then
      
      if(icoarse.eq.1.and. ichmod.eq.1) then
      write(6,*) 'change RG mode'
      write(80,*) 'change RG mode'
      icoarse=0
      ratelim=ratelim2
      goto 36
c####################################################################
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c####################################################################
      else
      write(6,*) 'converged'
      end if
      
      if(iWF.eq.1) then
      open(unit=86,file="WF_u6.txt",status='unknown')
      write(86,*) iLdim
      do l=1,iLdim
         do inum=1,m
            do i=1,nn
               write(86,*) stfai(i,inum,l)
            end do
         end do
      end do
      do l=1,iLdim
         do inum=1,m
            do i=1,nn
               write(86,*) stfai(i,inum+m,l)
            end do
         end do
      end do
      end if
      
      if(iengv.eq.1) then
      hold=h
      hinold=hin
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
      
      call engv(u,t1,t2,iLdim,ispin,iqpon
     &,faiengv,ihop1,ihop2,lsitex,lsitey,lsub,WORK2,LWORK2)
      h=hold
      hin=hinold
      end if

      
      
      
      write(6,*) 'stop dim.=',iLdim
      write(6,*) 

      goto 600
      end if
      

36    end do      
      
600   stop
      end
c end PIRG checker board
c----------------------------------------------------------------
      
      
      
      
      
      
      
      






















   
      
