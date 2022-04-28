c PIRG checker board --------------------------------------------------
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 JOBZ  
      character*1 UPLO
      
      parameter(nn=36,m=nn/2,u=6.0d0,t1=1.0d0,t2=1.0d0
     &,istdim=1,Lmax=500,iWORKlim=100000
c set up -----------------------------------------------------
     &,iRGtype=0,ierror=0
     &,ichmod=1,ratelim1=0.1d0,ratelim2=0.05d0
     &,delexp=5.0d-4/u,delu=0.5d0/u,limit=100000
c engv -----------------------------------------------------
     &,lqp=12,iengv=1,iqpon=1,ispin=0,iWF=0
c for additional alignment
     &,nx=4,ny=4)

c RG type iRGtype=0: independent RG
c RG type iRGtype=1: depend RG

c change RG mode ichmod=0: no chanege
c change RG mode ichmod=1: change

c error check ierror=0: error check and pass
c error check ierror=1: error check and stop

c iengv=0: no engv process
c iengv=1: engv process


      dimension gup(nn,nn,Lmax),gdown(nn,nn,Lmax)
      dimension gup2(nn,nn,Lmax,2),gdown2(nn,nn,Lmax,2)
      dimension ginup(Lmax),gindown(Lmax)
      dimension ginup2(Lmax,2),gindown2(Lmax,2)
      
      dimension h(Lmax,Lmax),hin(Lmax,Lmax)
      dimension hold(Lmax,Lmax),hinold(Lmax,Lmax)
      dimension hold2(Lmax,Lmax),hinold2(Lmax,Lmax)
      dimension hsub(Lmax,2),hinsub(Lmax,2)
      
      dimension fai0(nn,2*m)
      dimension fai(nn,2*m,Lmax),fai2(nn,2*m,Lmax)
      dimension fupdown(nn,2),faiold(nn,2*m),faiold2(nn,2*m,Lmax)
      
      dimension stexpu(nn,nn,Lmax),stexpd(nn,nn,Lmax)
      dimension oldexpu(nn,nn),oldexpd(nn,nn)
      dimension sup(nn,nn),sdown(nn,nn)
      
      dimension WORK2(iWORKlim),hgeed(Lmax,Lmax),hingeed(Lmax,Lmax)
      dimension deigen(Lmax),EERRBD(Lmax),ZERRBD(Lmax),RCONDZ(Lmax)
      dimension IWORK(Lmax)
      
      dimension cron(2*nn,2*nn),lurg(Lmax,Lmax),iss(2)
      dimension gsubu(nn,nn),gsubd(nn,nn)
      dimension a(m,m),b(m,m),IPIV(m),WORK(m)
      dimension alpha1(-1:1),alpha0(-1:1)
      dimension Wup(Lmax),pirgmp(2),isubL(Lmax)
      dimension isite(nn),idim0(Lmax),iRGdim(Lmax)
      dimension aold(m,m),bold(m,m)
      
c kinetic2      
      dimension hk(nn,nn),hk1(nn,nn)
     &,hk2(nn),ee(nn),hnqei(nn)
     &,hkL0(nn,nn)
      dimension huup(nn),hudown(nn)
      dimension huup2(nn,nn),hudown2(nn,nn)
      
      
      
c size depend !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
c shift site numbers cccccccccccccccccccccccccccc
      dimension ihop1(nn,4),ihop2(nn,2)
      dimension lsitex(nn),lsitey(nn)
      dimension lsitenum(10,10)
      dimension lsub(nn)
      dimension ip(10),im(10)
      dimension ldx(nn),ldy(nn)
      dimension lkx(22),lky(22),lknum(-2:3,0:3)
c end size depend !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      


c engv dimension !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension gqp(2*nn,2*nn),ginqp(Lmax)
      dimension gdelqp(2*nn,2*nn)
      dimension agqp(2*nn,2*nn),agdelqp(2*nn,2*nn),aqp(2*m,2*m)
      dimension aqpold(2*m,2*m)
      
      dimension faiqpL(2*nn,2*m,Lmax),faiqpR(2*nn,2*m,lqp)
      dimension IPIVqp(2*m),WORKqp(2*m)
      dimension hqp(Lmax,Lmax),hinqp(Lmax,Lmax)
      dimension wfin(Lmax,Lmax),delcro(2*nn,2*nn)
      
      dimension sq(nn,nn),cdens(nn,nn)
      
      dimension hspin(2*nn,2*nn,lqp),w(lqp),x(lqp),cof(lqp)
      
      
      dimension subh2(Lmax,Lmax)
      dimension subkin1(Lmax,Lmax,2),subkin2(Lmax,Lmax,2)
      dimension subkin3(Lmax,Lmax,2),subhdble(Lmax,Lmax,2)

      
      
      dimension subssnn1(Lmax,Lmax,2),subssnn2(Lmax,Lmax,2)
      dimension subssnn3(Lmax,Lmax,2)
      

c additional alignment ---------------------------------------
      dimension hnij(2*nn,2*nn)

      dimension subnk(22,Lmax,Lmax)
      dimension subsq(22,Lmax,Lmax)
      dimension subnq(22,Lmax,Lmax)
      
      
      dimension hnk(22),hsq(22),hnq(22)
c end engv dimension !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      data ixran/9862/
      data ITYPE/1/,UPLO/'U'/
      
      ixran=ixran*int(u)
      
      open(unit=80,file="eng_list_u6_m.txt",status='unknown')
      iRGdim=0
            
c      do i2=1,10
c      iRGdim(i2)=50*i2
c      end do
       iRGdim(1)=100
     

c      iRGdim=0
c      iRGdim(1)=10

      
      illup=0
      illdown=0
      
      ifresh=0
      label=0
      label2=0
      
      pai=dacos(-1.0d0)
      pieng2=100.0d0
      engold=100.0d0
      pieng4=100.0d0
      
      stockE=0.0d0
      huup=0.0d0
      hudown=0.0d0
      huup2=0.0d0
      hudown2=0.0d0
      alpha1=0.0d0
      alpha0=0.0d0
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
      lsite=0
      lsitex=0
      lsitey=0
      
      call shift(nn,ihop1,ihop2,lsitex,lsitey,lsitenum
     &,lsub,ip,im,ldx,ldy)
      
c set wave number ----------------------------------------------------
      lkx=0
      lky=0
      lknum=0
      do i=0,3,1
         lknum(i,0)=i+1
      end do
      do i=-2,3,1
         lknum(i,1)=i+7
      end do
      do i=-2,3,1
         lknum(i,2)=i+13
      end do
      do i=-2,3,1
         lknum(i,3)=i+19
      end do
      
      do iy=0,3,1
         do ix=-2,3,1
         if(lknum(ix,iy).ne.0) then
            lkx(lknum(ix,iy))=ix
            lky(lknum(ix,iy))=iy
         end if
         end do
      end do
      
c end set wave number ---------------------------------------------------


c make cronecker delta ---------------------------------------
         cron=0.0d0
         delcro=0.0d0
      do j=1,nn
      do i=1,nn
         
      if(i.eq.j) then
         cron(i,j)=1.0d0
         cron(i+nn,j+nn)=1.0d0
      end if
      
      do k=1,4
      if(i.eq.ihop1(j,k)) then
         delcro(i,j)=delcro(i,j)-t1
      end if
      end do
      
      do k=1,2
      if(i.eq.ihop2(j,k)) then
         delcro(i,j)=delcro(i,j)-t2
      end if
      end do
      
      do k=1,4
      if(i.eq.ihop1(j,k)) then
         delcro(i+nn,j+nn)=delcro(i+nn,j+nn)-t1
      end if
      end do
      
      do k=1,2
      if(i.eq.ihop2(j,k)) then
         delcro(i+nn,j+nn)=delcro(i+nn,j+nn)-t2
      end if
      end do
      
      
      end do
      end do
c end make cronecker delta -----------------------------------
      
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


       do 
       read(5,*) iLdim
       if(iLdim.eq.istdim) then
          goto 8989
       end if
       end do
8989   write(6,*) 'start_dim',iLdim
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
      
      if(iLdim.eq.1) then
      l=1
         do inum=1,m*2
            do i=1,nn
               fai0(i,inum)=fai(i,inum,l)
            end do
         end do
      end if
c end read First W.F.












c Expand dimension ---------------------------------------------------
      do 600 iRGstep=1,Lmax
      if(iRGdim(iRGstep).eq.0) goto 600
      
361   if(iRGtype.eq.1.and.iRGstep.ge.2) then
         iexpdim=iRGdim(iRGstep-1)+1
      else
         iexpdim=1
      end if

      delta=delexp
      call sht(u,delta,alpha1)
      call kinetic(nn,m,t1,t2,delta/2.0d0,hkL0,hk1,hk2
     &,ihop1,ihop2,hnqei,ee)
      
      write(6,*) 'LOCAL delta=',delta
      
      
      
      do 360 iLdim=iexpdim,iRGdim(iRGstep)
c select SH field ------------------------------------------------
      do ijm=1,nn
      if(ran2(ixran).lt.1.0d0/2.0d0) then
         is=1
      else
         is=-1
      end if
      huup(ijm)=dexp(alpha1(is))
      hudown(ijm)=dexp(alpha1(-is))
      end do
      
      do j=1,nn
      do i=1,nn
         stexpu(i,j,iLdim)=0.0d0
         stexpd(i,j,iLdim)=0.0d0
         do k=1,nn
            stexpu(i,j,iLdim)=stexpu(i,j,iLdim)
     &     +hkL0(i,k)*huup(k)*hkL0(k,j)
            stexpd(i,j,iLdim)=stexpd(i,j,iLdim)
     &     +hkL0(i,k)*hudown(k)*hkL0(k,j)
      end do
      end do
      end do
c end select SH field ------------------------------------------------

      
      do inum=1,m
         do i=1,nn
            fai(i,inum,iLdim)=0.0d0
            fai(i,inum+m,iLdim)=0.0d0
               do j=1,nn
                  fai(i,inum,iLdim)=fai(i,inum,iLdim)
     &           +stexpu(i,j,iLdim)*fai0(j,inum)
                  fai(i,inum+m,iLdim)=fai(i,inum+m,iLdim)
     &           +stexpd(i,j,iLdim)*fai0(j,inum+m)
               end do
         end do
      end do
      
      call schmidt(nn,m,Lmax,fai,fai2,fupdown,iLdim)
      do ll=1,iLdim
         call mkgreen(ll,iLdim,nn,m,Lmax,IPIV,WORK
     &  ,a,b,fai,gup,gdown,ginup,gindown)
         call mkH(ll,iLdim,nn,u,t1,t2,Lmax
     &  ,ihop1,ihop2,h,hin,gup,gdown,ginup,gindown)
      end do      
      
      LWORK=3*iLdim-1
      call geed0(Lmax,iLdim,label,h,hin
     &,hgeed,hingeed,deigen,Wup,eigen,WORK2,LWORK)
      pirgeng=eigen
      
      if(iRGtype.eq.1) then
         write(6,*) 'Local',iLdim,pirgeng
      end if
      
      if(label.eq.1) then
         label=0
         goto 361
      end if
      
360   end do      
c end Expand dimension -----------------------------------------------

      write(6,*) 'START eng',pirgeng,iRGdim(iRGstep)
      
      
      if(ichmod.eq.0) then
      icoarse=0
      ratelim=ratelim2
      else
      icoarse=1
      ratelim=ratelim1
      end if
      
      iLdim=iRGdim(iRGstep)
      write(80,*)
      write(80,*) 0,pirgeng,iLdim,iLdim

      write(6,*) 'WORK space',WORK2(1)
      LWORK=int(WORK2(1))
      if(LWORK.gt.iWORKlim) then
         write(6,*) 'over iWORKlim',iWORKlim,LWORK
         stop
      end if
c -----------------------------------------------------------------


























c--------------------------------------------------------------------      
c--------------------------------------------------------------------      
c PIRG step --------------------------------------------------------- 
c MAIN      
      delta=delu
      call sht(u,delta,alpha1)
      call kinetic(nn,m,t1,t2,delta/2.0d0,hkL0,hk1,hk2
     &            ,ihop1,ihop2,hnqei,ee)
     
      write(6,*) 
      write(6,*) 'U/t=',u
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


c coarse grained RG -------------------------------------------------  
      if(icoarse.eq.1) then
c select L order ---------------------------------------------    
      idim0=0
         i=0
         do
         i=i+1
            ixran2=int(ran2(ixran)*dble(iLdim))+1
            icount=0
               do j=2,i
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
      do j=1,nn
         do i=1,nn
            oldexpu(i,j)=stexpu(i,j,lite)
            oldexpd(i,j)=stexpd(i,j,lite)
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
      huup2=hkL0
      hudown2=hkL0

c preserve updated data ---------------------------------------
      call mkfai(nn,m,Lmax,lite,fai,fai0,stexpu,stexpd
     &          ,oldexpu,oldexpd,huup2,hudown2)
      call schmidt(nn,m,Lmax,fai,fai2,fupdown,lite)
      do ll=1,iLdim
      call mkgreen(ll,lite,nn,m,Lmax,IPIV,WORK
     &            ,a,b,fai,gup,gdown,ginup,gindown)
      call mkH(ll,lite,nn,u,t1,t2,Lmax
     &        ,ihop1,ihop2,h,hin,gup,gdown,ginup,gindown)
      end do
      do ll=lite,iLdim
         h(lite,ll)=h(ll,lite)
         hin(lite,ll)=hin(ll,lite)
      end do

c      write(*,*) "B"
c end preserve updated data ----------------------------------

      call geed0(Lmax,iLdim,label,h,hin
     &,hgeed,hingeed,deigen,Wup,eigen,WORK2,LWORK)
      pirgeng=eigen
c end first Hk/2 term-----------------------------------------------


c      write(*,*) "C"

c second Hu*Hk/2 term -----------------------------------------------
      do  ijmk=1,nn
      ijm=isite(ijmk)
      call Urg(nn,m,u,t1,t2,Lmax,iLdim,lite,ijm,alpha1,h,hin
     &,ihop1,ihop2,hsub,hinsub,hgeed,hingeed,Wup,WORK2,LWORK,iss
     &,gup,gdown,ginup,gindown,gsubu,gsubd,pirgeng,label,lurg
     &,gup2,gdown2,ginup2,gindown2,pirgmp,huup,hudown,deigen
     &,ipstepu,icoarse)
      end do
         
c      write(*,*) "CC"
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
      call mkfai(nn,m,Lmax,lite,fai,fai0,stexpu,stexpd
     &          ,oldexpu,oldexpd,huup2,hudown2)
      call schmidt(nn,m,Lmax,fai,fai2,fupdown,lite)
      do ll=1,iLdim
      call mkgreen(ll,lite,nn,m,Lmax,IPIV,WORK
     &            ,a,b,fai,gup,gdown,ginup,gindown)
      call mkH(ll,lite,nn,u,t1,t2,Lmax
     &        ,ihop1,ihop2,h,hin,gup,gdown,ginup,gindown)
      end do      
      do ll=lite,iLdim
         h(lite,ll)=h(ll,lite)
         hin(lite,ll)=hin(ll,lite)
      end do
c end preserve updated data ----------------------------------
      call geed0(Lmax,iLdim,label,h,hin
     &,hgeed,hingeed,deigen,Wup,eigen,WORK2,LWORK)
      pirgeng=eigen

      if((engold-pirgeng).lt.0) then
      pirgeng=engold
      do j=1,m*2
         do i=1,nn
         fai(i,j,lite)=faiold(i,j)
         end do
      end do   
      do j=1,nn
         do i=1,nn
         stexpu(i,j,lite)=oldexpu(i,j)
         stexpd(i,j,lite)=oldexpd(i,j)
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

cc error check -------------------------------------------------------
      if(label.eq.1) then
      label=0

      write(6,*) 'ERROR',lite

         write(6,*) 'back'
         if(ierror.eq.1) stop
         
         goto 361
      
      end if
cc end error check --------------------------------------------------------
cc-------------------------------------------------------------------------  
      
25    continue      
c end coarse grained RG -------------------------------------------------  

















      
      
      
      
      
      
      
c minute RG -------------------------------------------------  
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
               do j=2,i
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
      
c second Hk term ---------------------------------------------------      
      engold=pirgeng
      do j=1,2*m
         do i=1,nn
            faiold(i,j)=fai(i,j,lite)
         end do
      end do
      do j=1,nn
         do i=1,nn
            oldexpu(i,j)=stexpu(i,j,lite)
            oldexpd(i,j)=stexpd(i,j,lite)
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
      
      huup2=hkL0
      hudown2=hkL0
c preserve updated data ---------------------------------------
      call mkfai(nn,m,Lmax,lite,fai,fai0,stexpu,stexpd
     &          ,oldexpu,oldexpd,huup2,hudown2)
      call schmidt(nn,m,Lmax,fai,fai2,fupdown,lite)
      do ll=1,iLdim
      call mkgreen(ll,lite,nn,m,Lmax,IPIV,WORK
     &            ,a,b,fai,gup,gdown,ginup,gindown)
      call mkH(ll,lite,nn,u,t1,t2,Lmax
     &        ,ihop1,ihop2,h,hin,gup,gdown,ginup,gindown)
      end do      
      do ll=lite,iLdim
         h(lite,ll)=h(ll,lite)
         hin(lite,ll)=hin(ll,lite)
      end do
c end preserve updated data ----------------------------------

      call geed0(Lmax,iLdim,label,h,hin
     &,hgeed,hingeed,deigen,Wup,eigen,WORK2,LWORK)
      pirgeng=eigen

      if((engold-pirgeng).lt.0) then
      pirgeng=engold
      do j=1,m*2
         do i=1,nn
         fai(i,j,lite)=faiold(i,j)
         end do
      end do   
      do j=1,nn
         do i=1,nn
         stexpu(i,j,lite)=oldexpu(i,j)
         stexpd(i,j,lite)=oldexpd(i,j)
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

cc error check -------------------------------------------------------
      if(label.eq.1) then
      label=0

      write(6,*) '1st Kinetic term ERROR',lite

         write(6,*) 'back'
         if(ierror.eq.1) stop
         
         goto 361
      
      end if
cc end error check --------------------------------------------------------
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
               do j=2,i
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
      call mkgreen(ll,lite,nn,m,Lmax,IPIV,WORK
     &            ,a,b,fai,gup,gdown,ginup,gindown)
      end do
      
      engold=pirgeng
      do j=1,nn
         do i=1,nn
            oldexpu(i,j)=stexpu(i,j,lite)
            oldexpd(i,j)=stexpd(i,j,lite)
         end do
      end do
      
c first Hu term -----------------------------------------------
      do  ijmk=1,nn
      ijm=isite(ijmk)
      call Urg(nn,m,u,t1,t2,Lmax,iLdim,lite,ijm,alpha1,h,hin
     &,ihop1,ihop2,hsub,hinsub,hgeed,hingeed,Wup,WORK2,LWORK,iss
     &,gup,gdown,ginup,gindown,gsubu,gsubd,pirgeng,label,lurg
     &,gup2,gdown2,ginup2,gindown2,pirgmp,huup,hudown,deigen
     &,ipstepu,icoarse)
      end do
         
      huup2=0.0d0
      hudown2=0.0d0
      do i=1,nn
         huup2(i,i)=huup(i)
         hudown2(i,i)=hudown(i)
      end do
c preserve updated data ---------------------------------------
      call mkfai(nn,m,Lmax,lite,fai,fai0,stexpu,stexpd
     &          ,oldexpu,oldexpd,huup2,hudown2)
      call schmidt(nn,m,Lmax,fai,fai2,fupdown,lite)
      do ll=1,iLdim
      call mkgreen(ll,lite,nn,m,Lmax,IPIV,WORK
     &            ,a,b,fai,gup,gdown,ginup,gindown)
      call mkH(ll,lite,nn,u,t1,t2,Lmax
     &        ,ihop1,ihop2,h,hin,gup,gdown,ginup,gindown)
      end do      
      do ll=lite,iLdim
         h(lite,ll)=h(ll,lite)
         hin(lite,ll)=hin(ll,lite)
      end do
c end preserve updated data ----------------------------------
      call geed0(Lmax,iLdim,label,h,hin
     &,hgeed,hingeed,deigen,Wup,eigen,WORK2,LWORK)
      pirgeng=eigen
      
cc error check -------------------------------------------------------
      if(label.eq.1) then
      label=0

      write(6,*) 'Coulomb term ERROR',lite

         write(6,*) 'back'
         if(ierror.eq.1) stop
         
         goto 361
      
      end if
cc end error check --------------------------------------------------------
cc-------------------------------------------------------------------------  
      
24    continue
c end coulomb energy RG -------------------------------------------------
c------------------------------------------------------------------------








c --------------------------------------------------------------------
c third Hk term -------------------------------------------------
c select L order ---------------------------------------------    
      idim0=0
         i=0
         do
         i=i+1
            ixran2=int(ran2(ixran)*dble(iLdim))+1
            icount=0
               do j=2,i
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
      
c second Hk term ---------------------------------------------------      
      engold=pirgeng
      do j=1,2*m
         do i=1,nn
            faiold(i,j)=fai(i,j,lite)
         end do
      end do
      do j=1,nn
         do i=1,nn
            oldexpu(i,j)=stexpu(i,j,lite)
            oldexpd(i,j)=stexpd(i,j,lite)
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
      
      huup2=hkL0
      hudown2=hkL0
c preserve updated data ---------------------------------------
      call mkfai(nn,m,Lmax,lite,fai,fai0,stexpu,stexpd
     &          ,oldexpu,oldexpd,huup2,hudown2)
      call schmidt(nn,m,Lmax,fai,fai2,fupdown,lite)
      do ll=1,iLdim
      call mkgreen(ll,lite,nn,m,Lmax,IPIV,WORK
     &            ,a,b,fai,gup,gdown,ginup,gindown)
      call mkH(ll,lite,nn,u,t1,t2,Lmax
     &        ,ihop1,ihop2,h,hin,gup,gdown,ginup,gindown)
      end do      
      do ll=lite,iLdim
         h(lite,ll)=h(ll,lite)
         hin(lite,ll)=hin(ll,lite)
      end do
c end preserve updated data ----------------------------------

      call geed0(Lmax,iLdim,label,h,hin
     &,hgeed,hingeed,deigen,Wup,eigen,WORK2,LWORK)
      pirgeng=eigen

      if((engold-pirgeng).lt.0) then
      pirgeng=engold
      do j=1,m*2
         do i=1,nn
         fai(i,j,lite)=faiold(i,j)
         end do
      end do   
      do j=1,nn
         do i=1,nn
         stexpu(i,j,lite)=oldexpu(i,j)
         stexpd(i,j,lite)=oldexpd(i,j)
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
      
cc error check -------------------------------------------------------
      if(label.eq.1) then
      label=0

      write(6,*) '3rd Kinetic term ERROR',lite

         write(6,*) 'back'
         if(ierror.eq.1) stop
         
         goto 361
      
      end if
cc end error check --------------------------------------------------------
cc-------------------------------------------------------------------------  
      
c end third Hk term -------------------------------------------      
30    continue
c end minute RG -------------------------------------------------  
      end if







      
      
c      ipstep=int((dble(ipstepk)+dble(ipstepu)/dble(nn))/2.0d0)
      if(icoarse.eq.0) then
      ipstep=int(dble(ipstepu)/dble(nn))
      end if
      
      
      write(6,*) istep,pirgeng,ipstep,iLdim
      write(80,*) istep,pirgeng,ipstep,iLdim

cccccccc KOGA cccccccccccc  
c      if(istep==2) stop
cccccccccccccccccccccccccc
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
      
      if(iLdim.eq.Lmax.and.iWF.eq.1) then
      open(unit=86,file="WF_u6.txt",status='unknown')
      write(86,*) iLdim
      do l=1,iLdim
         do inum=1,m
            do i=1,nn
               write(86,*) fai(i,inum,l)
            end do
         end do
      end do
      do l=1,iLdim
         do inum=1,m
            do i=1,nn
               write(86,*) fai(i,inum+m,l)
            end do
         end do
      end do
      end if
      
      if(iengv.eq.1) then
      hold=h
      hinold=hin
      call engv(nx,ny,nn,m,u,t1,t2,Lmax,iLdim,lqp,ispin,iqpon
     &,fai,deigen,WORK2,hgeed,hingeed,h,hin,cron      
     &,Wup,isubL,ihop1,ihop2,lsitex,lsitey,lsite,lsub,ip,im
     &,ldx,ldy,gqp,ginqp
     &,gdelqp,agqp,agdelqp,aqp,faiqpL,faiqpR,IPIVqp,WORKqp,hqp,hinqp
     &,wfin,delcro,sq,cdens,hspin
     &,w,x,cof,subh2,subkin1,subkin2,subkin3,subhdble
     &,subssnn1,subssnn2,subssnn3
     &,hnij,subnk,subsq,subnq,hnk,hsq,hnq
     &,lkx,lky,lknum)
      h=hold
      hin=hinold
      end if
      
      
      
      
      write(6,*) 'stop dim.=',iLdim
      write(6,*) 
      write(6,*)
      
      goto 600
      end if
      

   
36    end do      
600   end do      
      
      stop
      end
c end PIRG checker board
c----------------------------------------------------------------
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mkfai(nn,m,Lmax,lite,fai,fai0,stexpu,stexpd
     &                ,oldexpu,oldexpd,huup2,hudown2)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      
      dimension fai(nn,2*m,Lmax),fai0(nn,2*m)
      dimension stexpu(nn,nn,Lmax),stexpd(nn,nn,Lmax)
      dimension oldexpu(nn,nn),oldexpd(nn,nn)
      dimension huup2(nn,nn),hudown2(nn,nn)

      do j=1,nn
         do i=1,nn
            stexpu(i,j,lite)=0.0d0
            stexpd(i,j,lite)=0.0d0
            do k=1,nn
               stexpu(i,j,lite)=stexpu(i,j,lite)
     &        +huup2(i,k)*oldexpu(k,j)
               stexpd(i,j,lite)=stexpd(i,j,lite)
     &        +hudown2(i,k)*oldexpd(k,j)
            end do
         end do
      end do
         
      do inum=1,m
         do i=1,nn
            fai(i,inum,lite)=0.0d0
            fai(i,inum+m,lite)=0.0d0
               do j=1,nn
                  fai(i,inum,lite)=fai(i,inum,lite)
     &           +stexpu(i,j,lite)*fai0(j,inum)
                  fai(i,inum+m,lite)=fai(i,inum+m,lite)
     &           +stexpd(i,j,lite)*fai0(j,inum+m)
               end do
         end do
      end do
      

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine echeck(Lmax,N,label,h,hin
     &,A,B,W,EERRBD,ZERRBD,RCONDZ,eigen,WORK,LWORK,IWORK)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      
      dimension h(Lmax,Lmax),hin(Lmax,Lmax)
      dimension A(N,N),B(N,N),EERRBD(N),ZERRBD(N),RCONDZ(N)
      dimension W(N),WORK(LWORK),IWORK(N)
      character*1 UPLO
      data UPLO/'U'/
      
      LDA=N
      LDB=N
      do lr=1,N
      do ll=1,lr
         A(ll,lr)=h(ll,lr)
         B(ll,lr)=hin(ll,lr)
      end do
      end do
      
      
c LAPACK BOOK p.85      
      EPSMCH=DLAMCH('E')
      ITYPE=1
      ANORM=DLANSY('1',UPLO,N,A,LDA,WORK,LWORK,INFO)
      BNORM=DLANSY('1',UPLO,N,B,LDB,WORK,LWORK,INFO)
      CALL DSYGV(ITYPE,'V',UPLO,N,A,LDA,B,LDB,W,WORK,LWORK,INFO)
      IF(INFO.GT.0.AND.INFO.LE.N) THEN
         PRINT *, 'DSYGV did not converge'
      ELSE IF(INFO.GT.N) THEN
         PRINT *, 'B not positive definite'
      ELSE IF(N.GT.0) THEN
         CALL DTRCON('1',UPLO,'N',N,B,LDB,RCONDB,WORK,IWORK,INFO)
         RCONDB=MAX(RCONDB,EPSMCH)
c         CALL DDGSNA('Eigenvectors',N,N,W,RCONDX,INFO)
         DO 10 I=1,N
            EERRBD(I)=(EPSMCH/RCONDB**2)*(ANORM/BNORM+ABS(W(I)))
c            ZERRBD(I)=(EPSMCH/RCONDB**3)*((ANORM/BNORM)/RCONDZ(I)
c     &                +(ABS(W(I))/RCONDZ(I))*RCONDB)
10       CONTINUE
      END IF
      
      
      
      write(6,*) 'error check',W(1),EERRBD(1)
      
      
      
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      








ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c kinetic energy matrix 2
      subroutine kinetic(nn,m,t1,t2,delta,hk,hk1,hk2
     &,ihop1,ihop2,hnqei,ee)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      
      dimension hk(nn,nn)
      dimension hk1(nn,nn)
      dimension hk2(nn)
      dimension hnqei(nn),ee(nn)
      
      dimension ihop1(nn,4),ihop2(nn,2)
      
      m2=m*2


c Matrix of kinetic energy
         hk=0.0d0
         hk1=0.0d0
         hk2=0.0d0

c Matrix of kinetic energy ccccccccccccccccccccccccccccccc
      do 7 i=1,nn
      do k=1,4
      hk(i,ihop1(i,k))=hk(i,ihop1(i,k))+t1*delta
      end do
      do k=1,2
      hk(i,ihop2(i,k))=hk(i,ihop2(i,k))+t2*delta
      end do
7     continue 
cccccccccccccccccccccccccccccccccccccccccccccccccccc


      ndim=nn
      call tred2(hk,ndim,ndim,hnqei,ee)
      call tqli(hnqei,ee,ndim,ndim,hk)
      

      do i=1,nn
         hk2(i)=dexp(hnqei(i))
         do j=1,nn
            hk1(j,i)=hk(i,j)
         end do
      end do

      hk=0.0d0
      do l=1,nn
      do i=1,nn
      do j=1,nn
         hk(i,l)=hk(i,l)+hk1(j,i)*hk2(j)*hk1(j,l)
      end do
      end do
      end do
      

c      write(6,*) 'reset dim.itehk=',ite,delta

      return
      end
c end kinetc energy matrix
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
















ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine sht(u,delta,alpha1)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      dimension alpha1(-1:1)
      
      b=dsqrt(dtanh(delta*u/4.0d0))
      a=datanh(b)
c      write(6,*) a
      
      alpha1(-1)=2.0d0*a*dble(-1)-delta*u/2.0d0
      alpha1(1)=2.0d0*a*dble(1)-delta*u/2.0d0
      alpha1(0)=0.0d0
c      write(6,*) alpha1(-1),alpha1(1)
c      stop
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function datanh(x)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      datanh=dlog(dsqrt((1.0d0+x)/(1.0d0-x)))
      return
      end 
cccccccccccccccccccccccccccccccccccccccccccccccccc











ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine schmidt(nn,m,Lmax,fai,fai2,fupdown,lite)

      implicit real*8(a-h,o-z)
      implicit integer*4(i-n) 

      dimension fai(nn,2*m,Lmax),fai2(nn,2*m,Lmax)
      dimension fupdown(nn,2)
      
      
c reset first vector ccccccccccccccccc
      j=1 
         resetu=0.0d0
         resetd=0.0d0
         do i=1,nn
            resetu=resetu+fai(i,j,lite)**2
            resetd=resetd+fai(i,j+m,lite)**2
         end do
         resetu = 1.0d0/dsqrt(resetu)
         resetd = 1.0d0/dsqrt(resetd)
         do i=1,nn
            fai(i,j,lite)=fai(i,j,lite)*resetu
            fai(i,j+m,lite)=fai(i,j+m,lite)*resetd
         end do
ccccccccccccccccccccccccccccccccccccccccc
      fai2=fai

      do 1 j=2,m
      
      fupdown=0.0d0

      do 3 k=1,j-1
      
      dinneru=0.0d0
      dinnerd=0.0d0
         do i2=1,nn
            dinneru=dinneru+fai2(i2,j,lite)*fai(i2,k,lite)
            dinnerd=dinnerd+fai2(i2,j+m,lite)*fai(i2,k+m,lite)
         end do
         
         do i3=1,nn
            fupdown(i3,1)=fupdown(i3,1)-dinneru*fai(i3,k,lite)
            fupdown(i3,2)=fupdown(i3,2)-dinnerd*fai(i3,k+m,lite)
         end do
         
3     continue      
      
         do i2=1,nn
            fai(i2,j,lite)=fai2(i2,j,lite)+fupdown(i2,1)
            fai(i2,j+m,lite)=fai2(i2,j+m,lite)+fupdown(i2,2)
         end do
         
c reset no.j2 vector cccccccccccccccccc     
         resetu=0.0d0
         resetd=0.0d0
         do i=1,nn
            resetu=resetu+fai(i,j,lite)**2
            resetd=resetd+fai(i,j+m,lite)**2
         end do

         resetu = 1.0d0/dsqrt(resetu)
         resetd = 1.0d0/dsqrt(resetd)

         do i=1,nn
            fai(i,j,lite)=fai(i,j,lite)*resetu
            fai(i,j+m,lite)=fai(i,j+m,lite)*resetd
         end do
ccccccccccccccccccccccccccccc      
1     continue




      return
      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc











ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mkgreen(ll,lr,nn,m,Lmax,IPIV,WORK
     &,a,b,fai,gup,gdown,ginup,gindown)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n) 
      character*1 UPLO 
      
      dimension a(m,m),b(m,m)
      dimension fai(nn,2*m,Lmax)
      dimension gup(nn,nn,Lmax),gdown(nn,nn,Lmax)
      dimension ginup(Lmax),gindown(Lmax)
      
      dimension IPIV(m),WORK(m)

      dimension bufu(m,nn),bufd(m,nn)

      do jnum=1,m
      do inum=1,m
         a(inum,jnum)=0.0d0
         b(inum,jnum)=0.0d0
         do i=1,nn
            a(inum,jnum)=a(inum,jnum)
     &     +fai(i,inum,ll)*fai(i,jnum,lr)
            b(inum,jnum)=b(inum,jnum)
     &     +fai(i,inum+m,ll)*fai(i,jnum+m,lr)
         end do
      end do
      end do
c end make a matrix
      call DGETRF(m,m,a,m,IPIV,INFO)
      if(INFO.ne.0) then
         write(6,*) 'LU_failure  INFO=',INFO
         stop
      end if
      ginup(ll)=1.0d0
      do j=1,m
         ginup(ll)=ginup(ll)*PIV(IPIV(j),j)*a(j,j)
      end do
      call DGETRI(m,a,m,IPIV,WORK,m,INFO)
      if(INFO.ne.0) then
         write(6,*) 'inverse_failure up spin INFO=',INFO
         stop
      end if
      

      call DGETRF(m,m,b,m,IPIV,INFO)
      if(INFO.ne.0) then
         write(6,*) 'LU_failure  INFO=',INFO
         stop
      end if
      gindown(ll)=1.0d0
      do j=1,m
          gindown(ll)=gindown(ll)*PIV(IPIV(j),j)*b(j,j)
      end do
      call DGETRI(m,b,m,IPIV,WORK,m,INFO)
      if(INFO.ne.0) then
         write(6,*) 'inverse_failure down spin INFO=',INFO
         stop
      end if


c make g matrix

      do j=1,nn
      do i=1,m
         bufu(i,j)=0.0d0
         bufd(i,j)=0.0d0
         do l=1,m
            bufu(i,j)=bufu(i,j)+a(i,l)*fai(j,l,ll)
            bufd(i,j)=bufd(i,j)+b(i,l)*fai(j,l+m,ll)
         end do
      end do
      end do

      do j=1,nn
      do i=1,nn
         gup(i,j,ll)=0.0d0
         gdown(i,j,ll)=0.0d0
         do k=1,m
            gup(i,j,ll)=gup(i,j,ll)+fai(i,k,lr)*bufu(k,j)
            gdown(i,j,ll)=gdown(i,j,ll)+fai(i,k+m,lr)*bufd(k,j)
         end do
      end do
      end do
      return
      end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc









ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mkH(ll,lr,nn,u,t1,t2,Lmax
     &,ihop1,ihop2,h,hin,gup,gdown,ginup,gindown)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n) 
      character*1 UPLO 
      
      dimension h(Lmax,Lmax),hin(Lmax,Lmax)
      dimension gup(nn,nn,Lmax),gdown(nn,nn,Lmax)
      dimension ginup(Lmax),gindown(Lmax)
      
      dimension ihop1(nn,4),ihop2(nn,2)

      hin(ll,lr)=gindown(ll)*ginup(ll)
      h(ll,lr)=0.0d0
      do i=1,nn
      
      do k=1,4
      h(ll,lr)=h(ll,lr)
     &-t1*(gup(i,ihop1(i,k),ll)+gdown(i,ihop1(i,k),ll))
      end do
      do k=1,2
      h(ll,lr)=h(ll,lr)
     &-t2*(gup(i,ihop2(i,k),ll)+gdown(i,ihop2(i,k),ll))
      end do
      h(ll,lr)=h(ll,lr)+u*gup(i,i,ll)*gdown(i,i,ll)
      
      end do
      h(ll,lr)=h(ll,lr)*hin(ll,lr)




      return
      end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc















ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine geed0(Lmax,iLdim,label,h,hin
     &,hgeed,hingeed,deigen,Wup,eigen,WORK2,LWORK)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 JOBZ,JOBZ2
      character*1 UPLO
      dimension h(Lmax,Lmax),hin(Lmax,Lmax)
      dimension hgeed(iLdim,iLdim),hingeed(iLdim,iLdim)
      dimension deigen(iLdim),Wup(iLdim),WORK2(LWORK)
      
      data ITYPE/1/,JOBZ/'N'/,JOBZ2/'N'/,UPLO/'U'/
      
      
      
      
      
c ED ccccccccccccccccccccccccccccccccccccccccccc      
      do lr=1,iLdim
      do ll=1,lr
         hgeed(ll,lr)=h(ll,lr)
         hingeed(ll,lr)=hin(ll,lr)
      end do
      end do
      
      call DSYGV(ITYPE,JOBZ,UPLO,iLdim,hgeed,iLdim,hingeed,iLdim
     &,deigen,WORK2,LWORK,INFO)
      if(INFO.gt.0.and.INFO.le.iLdim) then
         write(6,*) 'DSYGV did not converge'
         label=1
         deigen(1)=100.0d0
      else if(INFO.gt.N) then
         write(6,*) 'B not positive definite'
         label=1
         deigen(1)=100.0d0
      end if
      eigen=deigen(1)
cccccccccccccccccccccccccccccccccccccccccccccccc      


      return
2000  format(a,6f15.5)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc










ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine geedvarg2(nn,m,Lmax,iLdim,label,isubL
     &,h,hin,hgeed,hingeed,deigen,Wup,eigen
     &,WORK2,LWORK)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 JOBZ 
      character*1 JOBZ2 
      character*1 UPLO
c end QP      
      dimension isubL(iLdim)
      dimension h(Lmax,Lmax),hin(Lmax,Lmax)
      dimension hgeed(iLdim,iLdim),hingeed(iLdim,iLdim)
      dimension deigen(iLdim),Wup(iLdim),WORK2(LWORK)
      



      
      data ITYPE/1/,JOBZ/'V'/,JOBZ2/'N'/,UPLO/'U'/
      
      
      do i=1,iLdim
         isubL(i)=i
      end do 
      
      
      
c ED ccccccccccccccccccccccccccccccccccccccccccc      
      do lr=1,iLdim
      do ll=1,lr
         hgeed(ll,lr)=h(ll,lr)
         hingeed(ll,lr)=hin(ll,lr)
      end do
      end do
      
      call DSYGV(ITYPE,JOBZ,UPLO,iLdim,hgeed,iLdim,hingeed,iLdim
     &,deigen,WORK2,LWORK,INFO)
      if(INFO.ne.0.and.iLdim.ne.2) then
         write(6,*) iLdim,'GEED_failure st2 INFO=',INFO
         label=1
         
      deigen(1)=100.0d0
         
         
      end if
      eigen=deigen(1)
cccccccccccccccccccccccccccccccccccccccccccccccc      

c Weight sort ccccccccccccccccccccccccccccccccccccccccccccccccccc    
      do made=iLdim-1,1,-1
         do i=1,made
            if(abs(hgeed(i,1)).lt.abs(hgeed(i+1,1))) then
               
               W2=hgeed(i,1)
               iii=isubL(i)
                              
               hgeed(i,1)=hgeed(i+1,1)
               isubL(i)=isubL(i+1)
               
               hgeed(i+1,1)=W2
               isubL(i+1)=iii
               
            end if
         end do
      end do
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      
      return
2000  format(a,6f15.5)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc














ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine geedv(Lmax,iLdim,label,h,hin
     &,hgeed,hingeed,deigen,Wup,eigen,WORK2,LWORK)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 JOBZ 
      character*1 JOBZ2 
      character*1 UPLO
      dimension h(Lmax,Lmax),hin(Lmax,Lmax)
      dimension hgeed(iLdim,iLdim),hingeed(iLdim,iLdim)
      dimension deigen(iLdim),Wup(iLdim),WORK2(LWORK)
      
      data ITYPE/1/,JOBZ/'V'/,JOBZ2/'N'/,UPLO/'U'/
      
      
      
      
      
c ED ccccccccccccccccccccccccccccccccccccccccccc      
      do lr=1,iLdim
      do ll=1,lr
         hgeed(ll,lr)=h(ll,lr)
         hingeed(ll,lr)=hin(ll,lr)
      end do
      end do
      
      call DSYGV(ITYPE,JOBZ,UPLO,iLdim,hgeed,iLdim,hingeed,iLdim
     &,deigen,WORK2,LWORK,INFO)
      if(INFO.ne.0.and.iLdim.ne.2) then
         write(6,*) iLdim,'GEED_failure st2 INFO=',INFO
         label=1
         
c      call DSYEV(JOBZ2,UPLO,iLdim,hingeed,iLdim,Wup,WORK2,LWORK,INFO)
c      write(6,*) 
c      write(6,*) 'p define check list'
c      do i=1,iLdim
c         write(6,*) i,Wup(i)
c      end do
c         stop
      deigen(1)=100.0d0
         
         
      end if
      eigen=deigen(1)
cccccccccccccccccccccccccccccccccccccccccccccccc      





      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

































ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine Urg(nn,m,u,t1,t2,Lmax,iLdim,lite,ijm,alpha1,h,hin
     &,ihop1,ihop2,hsub,hinsub,hgeed,hingeed,Wup,WORK2,LWORK,iss
     &,gup,gdown,ginup,gindown,gsubu,gsubd,pirgeng,label,lurg
     &,gup2,gdown2,ginup2,gindown2,pirgmp,huup,hudown,deigen
     &,ipstepu,icoarse)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 UPLO
      
      dimension alpha1(-1:1)
      dimension hsub(iLdim,2),hinsub(iLdim,2)
      dimension h(Lmax,Lmax),hin(Lmax,Lmax)
      dimension hgeed(iLdim,iLdim),hingeed(iLdim,iLdim)
      dimension Wup(iLdim),WORK2(3*iLdim-1),deigen(iLdim)
      dimension gup(nn,nn,Lmax),gdown(nn,nn,Lmax)
      dimension ginup(Lmax),gindown(Lmax)
      dimension gsubu(nn,nn),gsubd(nn,nn)
      dimension gup2(nn,nn,Lmax,2),gdown2(nn,nn,Lmax,2)
      dimension ginup2(Lmax,2),gindown2(Lmax,2)
      dimension pirgmp(2),huup(nn),hudown(nn)
      
      dimension ihop1(nn,4),ihop2(nn,2),lurg(Lmax,Lmax),iss(2)
      do 1140 is=1,2
      is2=iss(is)
      
cccccccccccccccccccccccccccccccc
      ex1=dexp(alpha1(is2))-1.0d0
      ex2=dexp(alpha1(-is2))-1.0d0
cccccccccccccccccccccccccccccccc
!$omp parallel do private(ll,co1,co2)
      do 400 llsub=1,iLdim-1
      ll=lurg(llsub,lite)
      
      ginup2(ll,is)=(1.0d0
     &+ex1*gup(ijm,ijm,ll))*ginup(ll)
      gindown2(ll,is)=(1.0d0
     &+ex2*gdown(ijm,ijm,ll))*gindown(ll)

cccccccccccccccccccccccccccccccc
      co1=ex1/(1.0d0+gup(ijm,ijm,ll)*ex1)
      co2=ex2/(1.0d0+gdown(ijm,ijm,ll)*ex2)
cccccccccccccccccccccccccccccccc
      do j=1,nn
         do i=1,nn
      gup2(i,j,ll,is)=gup(i,j,ll)
     &-gup(i,ijm,ll)*gup(ijm,j,ll)*co1

      gdown2(i,j,ll,is)=gdown(i,j,ll)
     &-gdown(i,ijm,ll)*gdown(ijm,j,ll)*co2
         end do
      end do
      do j=1,nn
         gup2(ijm,j,ll,is)=gup2(ijm,j,ll,is)*(ex1+1.0d0)
         gdown2(ijm,j,ll,is)=gdown2(ijm,j,ll,is)*(ex2+1.0d0)
      end do
400   continue


      ll=lite
      ginsubu=(1.0d0
     &+ex1*gup(ijm,ijm,ll))*ginup(ll)
      ginsubd=(1.0d0
     &+ex2*gdown(ijm,ijm,ll))*gindown(ll)

cccccccccccccccccccccccccccccccc
      co1=ex1/(1.0d0+gup(ijm,ijm,ll)*ex1)
      co2=ex2/(1.0d0+gdown(ijm,ijm,ll)*ex2)
cccccccccccccccccccccccccccccccc
      do j=1,nn
         do i=1,nn
      gsubu(i,j)=gup(i,j,ll)
     &-gup(i,ijm,ll)*gup(ijm,j,ll)*co1

      gsubd(i,j)=gdown(i,j,ll)
     &-gdown(i,ijm,ll)*gdown(ijm,j,ll)*co2
         end do
      end do
      do j=1,nn
         gsubu(ijm,j)=gsubu(ijm,j)*(ex1+1.0d0)
         gsubd(ijm,j)=gsubd(ijm,j)*(ex2+1.0d0)
      end do
      ginup2(ll,is)=(1.0d0
     &+ex1*gsubu(ijm,ijm))*ginsubu
      gindown2(ll,is)=(1.0d0
     &+ex2*gsubd(ijm,ijm))*ginsubd


cccccccccccccccccccccccccccccccc
      co1=ex1/(1.0d0+gsubu(ijm,ijm)*ex1)
      co2=ex2/(1.0d0+gsubd(ijm,ijm)*ex2)
cccccccccccccccccccccccccccccccc
      do j=1,nn
         do i=1,nn
      gup2(i,j,ll,is)=gsubu(i,j)
     &-gsubu(i,ijm)*gsubu(ijm,j)*co1
      gdown2(i,j,ll,is)=gsubd(i,j)
     &-gsubd(i,ijm)*gsubd(ijm,j)*co2
         end do
      end do
      do i=1,nn  
         gup2(i,ijm,ll,is)=gup2(i,ijm,ll,is)*(ex1+1.0d0)
         gdown2(i,ijm,ll,is)=gdown2(i,ijm,ll,is)*(ex2+1.0d0)
      end do

c make Hamiltonial --------------------------------------------------
!$omp parallel do private(i1,i2) 
      do ll=1,iLdim
      hinsub(ll,is)=gindown2(ll,is)*ginup2(ll,is)
      hsub(ll,is)=0.0d0
      
      do i=1,nn
      
      do k=1,4
      hsub(ll,is)=hsub(ll,is)
     &-t1*(gup2(i,ihop1(i,k),ll,is)+gdown2(i,ihop1(i,k),ll,is))
      end do
      do k=1,2
      hsub(ll,is)=hsub(ll,is)
     &-t2*(gup2(i,ihop2(i,k),ll,is)+gdown2(i,ihop2(i,k),ll,is))
      end do
      hsub(ll,is)=hsub(ll,is)
     &+u*gup2(i,i,ll,is)*gdown2(i,i,ll,is)
      
      end do
      
      hsub(ll,is)=hsub(ll,is)*hinsub(ll,is)
c end make Hamiltonial ------------------------------------------------


        i1=ll
        i2=lite
        if(ll>lite)then
        i1=lite
        i2=ll
        end if
         h(i1,i2)=hsub(ll,is)
         hin(i1,i2)=hinsub(ll,is)
      
        end do
      call geed0(Lmax,iLdim,label,h,hin
     &,hgeed,hingeed,deigen,Wup,eigen,WORK2,LWORK)
      pirgmp(is)=eigen
1140  continue

c no change
      if(pirgmp(1).ge.pirgeng.and.pirgmp(2).ge.pirgeng
     &.and.icoarse.eq.0) then
      huup(ijm)=1.0d0
      hudown(ijm)=1.0d0
      else
      
      
      ipstepu=ipstepu+1
c minus SH field stable
      if((pirgmp(2).lt.pirgmp(1))) then
      huup(ijm)=dexp(alpha1(-1))
      hudown(ijm)=dexp(alpha1(1))
      
!$omp parallel do private(i1,i2)   
      do l=1,iLdim
        i1=l
        i2=lite
        if(l>lite) then
        i1=lite
        i2=l
        end if
        h(i1,i2)=hsub(l,2)
        hin(i1,i2)=hinsub(l,2)

         ginup(l)=ginup2(l,2)
         gindown(l)=gindown2(l,2)
         do j=1,nn
            do i=1,nn
               gup(i,j,l)=gup2(i,j,l,2)
               gdown(i,j,l)=gdown2(i,j,l,2)
            end do
         end do
      end do
      pirgeng=pirgmp(2)

c plus SH field stable
      else
      huup(ijm)=dexp(alpha1(1))
      hudown(ijm)=dexp(alpha1(-1))

!$omp parallel do private(i1,i2)
      do l=1,iLdim
        i1=l
        i2=lite
        if(l>lite) then
        i1=lite
        i2=l
        end if
        h(i1,i2)=hsub(l,1)
        hin(i1,i2)=hinsub(l,1)
         ginup(l)=ginup2(l,1)
         gindown(l)=gindown2(l,1)
         do j=1,nn
            do i=1,nn
               gup(i,j,l)=gup2(i,j,l,1)
               gdown(i,j,l)=gdown2(i,j,l,1)
            end do
         end do
      end do
      pirgeng=pirgmp(1)
      end if
      
      end if
      
20    return
      end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
















      
      
      
      
      
      SUBROUTINE choldc(a,n,np,p)
      integer*4 n,np
      REAL*8 a(np,np),p(n)
      integer*4 i,j,k
      REAL*8 sum
      do 13 i=1,n
        do 12 j=i,n
          sum=a(i,j)
          do 11 k=i-1,1,-1
            sum=sum-a(i,k)*a(j,k)
11        continue
          if(i.eq.j)then
c          write(6,*) 'sum',sum
            if(sum.le.0.d0)pause 'choldc failed'
            p(i)=dsqrt(sum)
          else
            a(j,i)=sum/p(i)
          endif
          if(i.eq.j) then
c          write(6,*) p(i)
          else
c          write(6,*) j,i,a(j,i)
          end if
12      continue
13    continue
      return
      END
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      SUBROUTINE tred2(a,n,np,d,e)
      integer*4 n,np
      REAL*8 a(np,np),d(np),e(np)
      integer*4 i,j,k,l
      REAL*8 f,g,h,hh,scale
      do 18 i=n,2,-1
        l=i-1
        h=0.d0
        scale=0.d0
        if(l.gt.1)then
          do 11 k=1,l
            scale=scale+abs(a(i,k))
11        continue
          if(scale.eq.0.d0)then
            e(i)=a(i,l)
          else
            do 12 k=1,l
              a(i,k)=a(i,k)/scale
              h=h+a(i,k)**2
12          continue
            f=a(i,l)
            g=-sign(dsqrt(h),f)
            e(i)=scale*g
            h=h-f*g
            a(i,l)=f-g
            f=0.d0
            do 15 j=1,l
C     Omit following line if finding only eigenvalues
              a(j,i)=a(i,j)/h
              g=0.d0
              do 13 k=1,j
                g=g+a(j,k)*a(i,k)
13            continue
              do 14 k=j+1,l
                g=g+a(k,j)*a(i,k)
14            continue
              e(j)=g/h
              f=f+e(j)*a(i,j)
15          continue
            hh=f/(h+h)
            do 17 j=1,l
              f=a(i,j)
              g=e(j)-hh*f
              e(j)=g
              do 16 k=1,j
                a(j,k)=a(j,k)-f*e(k)-g*a(i,k)
16            continue
17          continue
          endif
        else
          e(i)=a(i,l)
        endif
        d(i)=h
18    continue
C     Omit following line if finding only eigenvalues.
      d(1)=0.d0
      e(1)=0.d0
      do 24 i=1,n
C     Delete lines from here ...
        l=i-1
        if(d(i).ne.0.d0)then
          do 22 j=1,l
            g=0.d0
            do 19 k=1,l
              g=g+a(i,k)*a(k,j)
19          continue
            do 21 k=1,l
              a(k,j)=a(k,j)-g*a(k,i)
21          continue
22        continue
        endif
C     ... to here when finding only eigenvalues.
        d(i)=a(i,i)
C     Also delete lines from here ...
        a(i,i)=1.d0
        do 23 j=1,l
          a(i,j)=0.d0
          a(j,i)=0.d0
23      continue
C     ... to here when finding only eigenvalues.
24    continue
      return
      END
      
      
      
      
      
      
      SUBROUTINE tqli(d,e,n,np,z)
      integer*4 n,np
      REAL*8 d(np),e(np),z(np,np)
CU    USES pythag
      integer*4 i,iter,k,l,m
      REAL*8 b,c,dd,f,g,p,r,s,pythag
      do 11 i=2,n
        e(i-1)=e(i)
11    continue
      e(n)=0.d0
      do 15 l=1,n
        iter=0
1       do 12 m=l,n-1
          dd=abs(d(m))+abs(d(m+1))
          if (abs(e(m))+dd.eq.dd) goto 2
12      continue
        m=n
2       if(m.ne.l)then
          if(iter.eq.30)pause 'too many iterations in tqli'
          iter=iter+1
          g=(d(l+1)-d(l))/(2.d0*e(l))
          r=pythag(g,1.d0)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.d0
          c=1.d0
          p=0.d0
          do 14 i=m-1,l,-1
            f=s*e(i)
            b=c*e(i)
            r=pythag(f,g)
            e(i+1)=r
            if(r.eq.0.d0)then
              d(i+1)=d(i+1)-p
              e(m)=0.d0
              goto 1
            endif
            s=f/r
            c=g/r
            g=d(i+1)-p
            r=(d(i)-g)*s+2.d0*c*b
            p=s*r
            d(i+1)=g+p
            g=c*r-b
C     Omit lines from here ...
            do 13 k=1,n
              f=z(k,i+1)
              z(k,i+1)=s*z(k,i)+c*f
              z(k,i)=c*z(k,i)-s*f
13          continue
C     ... to here when finding only eigenvalues.
14        continue
          d(l)=d(l)-p
          e(l)=g
          e(m)=0.d0
          goto 1
        endif
15    continue
      return
      END


      FUNCTION PIV(i,j)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
         if(i.ne.j) then
            PIV=-1.0d0
         else
            PIV=1.0d0
         end if
      return
      END


      FUNCTION pythag(a,b)
      REAL*8 a,b,pythag
      REAL*8 absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*dsqrt(1.d0+(absb/absa)**2)
      else
        if(absb.eq.0.d0)then
          pythag=0.d0
        else
          pythag=absb*dsqrt(1.d0+(absa/absb)**2)
        endif
      endif
      return
      END
      
      
      
      
      
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
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
      ran2=min(AM*iy,RNMX)
      return
      END
      
      
      
      
      
      
      
      

      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

















ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dimer(nn,ll,lr,i,j,lqp,iqp
     &,gqp,agqp,ldx,ldy,awick)

      implicit real*8(a-h,o-z)
      implicit integer*4(i-n) 

      dimension gqp(nn*2,nn*2),agqp(nn*2,nn*2)
      dimension ldx(nn),ldy(nn)

      awick=0.0d0
      if(ll.eq.1.and.lr.eq.1.and.i.eq.1.and.j.eq.1) then
       end if
c SiuSidSjuSjd ------------------------------------------------
      i1=i
      i2=i+nn
      j1=ldx(i)+nn
      j2=ldx(i)
      k1=j
      k2=j+nn
      l1=ldy(j)+nn
      l2=ldy(j)
      awick=awick
     &+0.25d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)
     
c SidSiuSjdSju
      i1=i+nn
      i2=i
      j1=ldx(i)
      j2=ldx(i)+nn
      k1=j+nn
      k2=j
      l1=ldy(j)
      l2=ldy(j)+nn
      awick=awick
     &+0.25d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

c SiuSidSjdSju
      i1=i
      i2=i+nn
      j1=ldx(i)+nn
      j2=ldx(i)
      k1=j+nn
      k2=j
      l1=ldy(j)
      l2=ldy(j)+nn
      awick=awick
     &+0.25d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

c SidSiuSjuSjd
      i1=i+nn
      i2=i
      j1=ldx(i)
      j2=ldx(i)+nn
      k1=j
      k2=j+nn
      l1=ldy(j)+nn
      l2=ldy(j)
      awick=awick
     &+0.25d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)
c --------------------------------------------------------------     
      
      
      
c SuSdSzSz ------------------------------------------------------      
      i1=i
      i2=i+nn
      j1=ldx(i)+nn
      j2=ldx(i)
      k1=j
      k2=j
      l1=ldy(j)
      l2=ldy(j)
      awick=awick
     &+0.125d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)
      
      i1=i
      i2=i+nn
      j1=ldx(i)+nn
      j2=ldx(i)
      k1=j+nn
      k2=j+nn
      l1=ldy(j)
      l2=ldy(j)
      awick=awick
     &-0.125d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)
      
      i1=i
      i2=i+nn
      j1=ldx(i)+nn
      j2=ldx(i)
      k1=j
      k2=j
      l1=ldy(j)+nn
      l2=ldy(j)+nn
      awick=awick
     &-0.125d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)
      
      i1=i
      i2=i+nn
      j1=ldx(i)+nn
      j2=ldx(i)
      k1=j+nn
      k2=j+nn
      l1=ldy(j)+nn
      l2=ldy(j)+nn
      awick=awick
     &+0.125d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)
c -----------------------------------------------------------



c SdSuSzSz ----------------------------------------------------      
      i1=i+nn
      i2=i
      j1=ldx(i)
      j2=ldx(i)+nn
      k1=j
      k2=j
      l1=ldy(j)
      l2=ldy(j)
      awick=awick
     &+0.125d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)
      
      i1=i+nn
      i2=i
      j1=ldx(i)
      j2=ldx(i)+nn
      k1=j+nn
      k2=j+nn
      l1=ldy(j)
      l2=ldy(j)
      awick=awick
     &-0.125d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i+nn
      i2=i
      j1=ldx(i)
      j2=ldx(i)+nn
      k1=j
      k2=j
      l1=ldy(j)+nn
      l2=ldy(j)+nn
      awick=awick
     &-0.125d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i+nn
      i2=i
      j1=ldx(i)
      j2=ldx(i)+nn
      k1=j+nn
      k2=j+nn
      l1=ldy(j)+nn
      l2=ldy(j)+nn
      awick=awick
     &+0.125d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)
c --------------------------------------------------------------------

c SzSzSuSd ----------------------------------------------------------
      i1=i
      i2=i
      j1=ldx(i)
      j2=ldx(i)
      k1=j
      k2=j+nn
      l1=ldy(j)+nn
      l2=ldy(j)
      awick=awick
     &+0.125d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)
      
      i1=i+nn
      i2=i+nn
      j1=ldx(i)
      j2=ldx(i)
      k1=j
      k2=j+nn
      l1=ldy(j)+nn
      l2=ldy(j)
      awick=awick
     &-0.125d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i
      i2=i
      j1=ldx(i)+nn
      j2=ldx(i)+nn
      k1=j
      k2=j+nn
      l1=ldy(j)+nn
      l2=ldy(j)
      awick=awick
     &-0.125d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i+nn
      i2=i+nn
      j1=ldx(i)+nn
      j2=ldx(i)+nn
      k1=j
      k2=j+nn
      l1=ldy(j)+nn
      l2=ldy(j)
      awick=awick
     &+0.125d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)
c ----------------------------------------------------------------



c SzSzSdSu
      i1=i
      i2=i
      j1=ldx(i)
      j2=ldx(i)
      k1=j+nn
      k2=j
      l1=ldy(j)
      l2=ldy(j)+nn
      awick=awick
     &+0.125d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i+nn
      i2=i+nn
      j1=ldx(i)
      j2=ldx(i)
      k1=j+nn
      k2=j
      l1=ldy(j)
      l2=ldy(j)+nn
      awick=awick
     &-0.125d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i
      i2=i
      j1=ldx(i)+nn
      j2=ldx(i)+nn
      k1=j+nn
      k2=j
      l1=ldy(j)
      l2=ldy(j)+nn
      awick=awick
     &-0.125d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i+nn
      i2=i+nn
      j1=ldx(i)+nn
      j2=ldx(i)+nn
      k1=j+nn
      k2=j
      l1=ldy(j)
      l2=ldy(j)+nn
      awick=awick
     &+0.125d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)
c ------------------------------------------------------------



c SzSzSzSz ---------------------------------------------------
      i1=i
      i2=i
      j1=ldx(i)
      j2=ldx(i)
      k1=j
      k2=j
      l1=ldy(j)
      l2=ldy(j)
      awick=awick
     &+0.0625d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i+nn
      i2=i+nn
      j1=ldx(i)
      j2=ldx(i)
      k1=j
      k2=j
      l1=ldy(j)
      l2=ldy(j)
      awick=awick
     &-0.0625d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i
      i2=i
      j1=ldx(i)+nn
      j2=ldx(i)+nn
      k1=j
      k2=j
      l1=ldy(j)
      l2=ldy(j)
      awick=awick
     &-0.0625d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i
      i2=i
      j1=ldx(i)
      j2=ldx(i)
      k1=j+nn
      k2=j+nn
      l1=ldy(j)
      l2=ldy(j)
      awick=awick
     &-0.0625d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i
      i2=i
      j1=ldx(i)
      j2=ldx(i)
      k1=j
      k2=j
      l1=ldy(j)+nn
      l2=ldy(j)+nn
      awick=awick
     &-0.0625d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i+nn
      i2=i+nn
      j1=ldx(i)+nn
      j2=ldx(i)+nn
      k1=j
      k2=j
      l1=ldy(j)
      l2=ldy(j)
      awick=awick
     &+0.0625d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i+nn
      i2=i+nn
      j1=ldx(i)
      j2=ldx(i)
      k1=j+nn
      k2=j+nn
      l1=ldy(j)
      l2=ldy(j)
      awick=awick
     &+0.0625d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i+nn
      i2=i+nn
      j1=ldx(i)
      j2=ldx(i)
      k1=j
      k2=j
      l1=ldy(j)+nn
      l2=ldy(j)+nn
      awick=awick
     &+0.0625d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i
      i2=i
      j1=ldx(i)+nn
      j2=ldx(i)+nn
      k1=j+nn
      k2=j+nn
      l1=ldy(j)
      l2=ldy(j)
      awick=awick
     &+0.0625d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i
      i2=i
      j1=ldx(i)+nn
      j2=ldx(i)+nn
      k1=j
      k2=j
      l1=ldy(j)+nn
      l2=ldy(j)+nn
      awick=awick
     &+0.0625d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i
      i2=i
      j1=ldx(i)
      j2=ldx(i)
      k1=j+nn
      k2=j+nn
      l1=ldy(j)+nn
      l2=ldy(j)+nn
      awick=awick
     &+0.0625d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i+nn
      i2=i+nn
      j1=ldx(i)+nn
      j2=ldx(i)+nn
      k1=j+nn
      k2=j+nn
      l1=ldy(j)
      l2=ldy(j)
      awick=awick
     &-0.0625d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i+nn
      i2=i+nn
      j1=ldx(i)+nn
      j2=ldx(i)+nn
      k1=j
      k2=j
      l1=ldy(j)+nn
      l2=ldy(j)+nn
      awick=awick
     &-0.0625d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i+nn
      i2=i+nn
      j1=ldx(i)
      j2=ldx(i)
      k1=j+nn
      k2=j+nn
      l1=ldy(j)+nn
      l2=ldy(j)+nn
      awick=awick
     &-0.0625d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i
      i2=i
      j1=ldx(i)+nn
      j2=ldx(i)+nn
      k1=j+nn
      k2=j+nn
      l1=ldy(j)+nn
      l2=ldy(j)+nn
      awick=awick
     &-0.0625d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)

      i1=i+nn
      i2=i+nn
      j1=ldx(i)+nn
      j2=ldx(i)+nn
      k1=j+nn
      k2=j+nn
      l1=ldy(j)+nn
      l2=ldy(j)+nn
      awick=awick
     &+0.0625d0*(
c a term
     &+(gqp(i1,i2)*gqp(j1,j2)
     &+gqp(i1,j2)*agqp(i2,j1))
     &*(gqp(k1,k2)*gqp(l1,l2)
     &+gqp(k1,l2)*agqp(k2,l1))
c b term
     &+(gqp(i1,i2)*gqp(j1,k2)
     &+gqp(i1,k2)*agqp(i2,j1))
     &*(agqp(j2,k1)*gqp(l1,l2)
     &-agqp(j2,l1)*gqp(k1,l2))
c c term
     &+(gqp(i1,i2)*gqp(j1,l2)
     &+gqp(i1,l2)*agqp(i2,j1))
     &*(agqp(j2,k1)*agqp(k2,l1)
     &+agqp(j2,l1)*gqp(k1,k2))
c d term
     &+(-gqp(i1,j2)*gqp(j1,k2)
     &+gqp(i1,k2)*gqp(j1,j2))
     &*(agqp(i2,k1)*gqp(l1,l2)
     &-agqp(i2,l1)*gqp(k1,l2))
c e term
     &+(-gqp(i1,j2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,j2))
     &*(agqp(i2,k1)*agqp(k2,l1)
     &+agqp(i2,l1)*gqp(k1,k2))
c f term
     &+(-gqp(i1,k2)*gqp(j1,l2)
     &+gqp(i1,l2)*gqp(j1,k2))
     &*(-agqp(i2,k1)*agqp(j2,l1)
     &+agqp(i2,l1)*agqp(j2,k1))
     &)
     






      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc











ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mkHqp(ll,lr,nn,u,t1,t2,Lmax
     &,ihop1,ihop2,h,hin,gqp,ginqp)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n) 
      character*1 UPLO 
      
c QP    
      dimension gqp(2*nn,2*nn),ginqp(Lmax)
c end QP      
      dimension ihop1(nn,4),ihop2(nn,2)
      dimension h(Lmax,Lmax),hin(Lmax,Lmax)
      

      hin(ll,lr)=ginqp(ll)
      h(ll,lr)=0.0d0
      do i=1,nn
      
      do k=1,4
      h(ll,lr)=h(ll,lr)
     &-t1*(gqp(i,ihop1(i,k))+gqp(i+nn,ihop1(i,k)+nn))
      end do
      do k=1,2
      h(ll,lr)=h(ll,lr)
     &-t2*(gqp(i,ihop2(i,k))+gqp(i+nn,ihop2(i,k)+nn))
      end do
      h(ll,lr)=h(ll,lr)
     &+u*gqp(i,i)*gqp(i+nn,i+nn)
     &-u*gqp(i+nn,i)*gqp(i,i+nn)
      end do
      h(ll,lr)=h(ll,lr)*hin(ll,lr)


      return
      end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mkgreenqp(ll,nn,m,Lmax,lqp,iqp,IPIVqp,WORKqp
     &,aqp,faiqpL,faiqpR,gqp,ginqp)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n) 
      character*1 UPLO 
      
c QP    
      dimension gqp(2*nn,2*nn),ginqp(Lmax)
      dimension aqp(2*m,2*m)
      dimension faiqpL(2*nn,2*m,Lmax)
      dimension faiqpR(2*nn,2*m,lqp)
      dimension IPIVqp(2*m),WORKqp(2*m)

        dimension buf(2*m,2*nn)
c end QP      
      
      
      m2=2*m
!$omp parallel do       
      do inum=1,m2
      do jnum=1,m2
      aqp(inum,jnum)=0.0d0
         do i=1,2*nn
         aqp(inum,jnum)=aqp(inum,jnum)
     &+faiqpL(i,inum,ll)*faiqpR(i,jnum,iqp)
         end do
      end do
      end do

c end make a matrix
      call DGETRF(m2,m2,aqp,m2,IPIVqp,INFO)
      if(INFO.ne.0) then
         write(6,*) 'LU_failure  INFO=',INFO
         stop
      end if
      ginqp(ll)=1.0d0
      do j=1,m2
         ginqp(ll)=ginqp(ll)*PIV(IPIVqp(j),j)*aqp(j,j)
      end do
      
c make g matrix
      call DGETRI(m2,aqp,m2,IPIVqp,WORKqp,m2,INFO)
      if(INFO.ne.0) then
         write(6,*) 'inverse_failure  INFO=',INFO
         stop
      end if
   
c      do j=1,nn*2
c      do i=1,nn*2
c         gqp(i,j)=0.0d0
c         do k=1,m2
c         do l=1,m2
c            gqp(i,j)=gqp(i,j)
c     &     +faiqpR(i,k,iqp)*aqp(k,l)*faiqpL(j,l,ll)
c         end do
c         end do
c      end do
c      end do
!$omp parallel do private(s)
        do j=1,nn*2
        do i=1,m2
        s=0.0d0
        do l=1,m2
        s = s + aqp(i,l)*faiqpL(j,l,ll)
        end do
        buf(i,j)=s
        end do
        end do  

!$omp parallel do private(s)
        do j=1,nn*2
        do i=1,nn*2
        s=0.0d0
        do k=1,m2
        s =s + faiqpR(i,k,iqp)*buf(k,j)
        end do
        gqp(i,j)=s
        end do  
        end do

    


      return
      end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mkgreenqp0(ll,nn,m,Lmax,lqp,iqp,IPIVqp,WORKqp
     &,aqp,aqpold,faiqpL,faiqpR,gqp,ginqp)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n) 
      character*1 UPLO 
      
c QP    
      dimension gqp(2*nn,2*nn),ginqp(Lmax)
      dimension aqp(2*m,2*m),aqpold(2*m,2*m)
      dimension faiqpL(2*nn,2*m,Lmax)
      dimension faiqpR(2*nn,2*m,lqp)
      dimension IPIVqp(2*m),WORKqp(2*m)
c end QP      
      
      
      m2=2*m
      
      do inum=1,m2
      do jnum=1,m2
      aqp(inum,jnum)=0.0d0
         do i=1,2*nn
         aqp(inum,jnum)=aqp(inum,jnum)
     &+faiqpL(i,inum,ll)*faiqpR(i,jnum,iqp)
         end do
      end do
      end do
      aqpold=aqp
c end make a matrix
      call DGETRF(m2,m2,aqp,m2,IPIVqp,INFO)
      if(INFO.ne.0) then
         write(6,*) 'LU_failure  INFO=',INFO
         stop
      end if
      ginqp(ll)=1.0d0
      do j=1,m2
         ginqp(ll)=ginqp(ll)*PIV(IPIVqp(j),j)*aqp(j,j)
      end do
      
c make g matrix
      call DGETRI(m2,aqp,m2,IPIVqp,WORKqp,m2,INFO)
      if(INFO.ne.0) then
         write(6,*) 'inverse_failure  INFO=',INFO
         stop
      end if
   

      do j=1,nn*2
      do i=1,nn*2
         gqp(i,j)=0.0d0
         do k=1,m2
         do l=1,m2
            gqp(i,j)=gqp(i,j)
     &     +faiqpR(i,k,iqp)*aqp(k,l)*faiqpL(j,l,ll)
         end do
         end do
      end do
      end do

      return
      end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc







ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c spin rotation matrix
      subroutine spinrt(nn,m,hspin,lqp,cof,w,x,ispin)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      
      dimension hspin(nn*2,nn*2,lqp)
      
      dimension x(lqp),w(lqp),cof(lqp)
      dimension rambda(nn*2),spiny(nn*2,nn*2)
      dimension spiny2(nn*2,nn*2)
      m2=m*2
      
      
      x1=-1.0d0
      x2=1.0d0
      call gauleg(x1,x2,x,w,lqp)
c      do i=1,lqp
c      write(6,*) x(i),w(i),i
c      end do 
c      stop



      do 6 l=1,lqp
      do 6 nr=1,nn*2
      do 6 nl=1,nn*2
         spiny(nl,nr)=0.0d0
         hspin(nl,nr,l)=0.0d0
6     continue

c make Sy matrix cccccccccccccccccccccccccccccccccccccccccccc      
      do 7 i=1,nn
      spiny(i,i+nn)=spiny(i,i+nn)+1.0d0
      spiny(i+nn,i)=spiny(i+nn,i)-1.0d0
      
      spiny2(i,i+nn)=spiny2(i,i+nn)+1.0d0
      spiny2(i+nn,i)=spiny2(i+nn,i)-1.0d0
7     continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      



      do 10 ii=1,lqp
      beta=dacos(x(ii))
ccccccccccccccccccccccccc
c      beta=dacos(1.0d0)
c      w(ii)=1.0d0
cccccccccccccccccccccccccc
      do 8 nr=1,nn*2
      do 8 nl=1,nn*2
      hspin(nl,nr,ii)=hspin(nl,nr,ii)
     &+spiny(nl,nr)*dsin(beta/2.0d0)
      if(nl.eq.nr) then
      hspin(nl,nr,ii)=hspin(nl,nr,ii)
     &+dcos(beta/2.0d0)
      end if
8     continue
      
      
c      do 9 nr=1,nn*2
c      do 9 nl=1,nn*2
c         write(6,*) nl,nr,hspin(nl,nr,1),spiny(nl,nr)
c9     continue
c         write(6,*) hspin
10    continue       
      
      if(ispin.eq.0) then
         do i=1,lqp
            cof(i)=w(i)/2.0d0
         end do
      else if(ispin.eq.1) then
         do i=1,lqp
            cof(i)=x(i)*w(i)*1.5d0
         end do
      end if
      
      return
      end
c end spin energy matrix
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc






cccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER*4 n
      REAL*8 x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER*4 i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

































ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Exact grand state energy 4by4 u_0 -32 (for one spin)
c PIRG checker board
      subroutine engv(nx,ny,nn,m,u,t1,t2,Lmax,iLdim,lqp,ispin,iqpon
     &,fai,deigen,WORK2,hgeed,hingeed,h,hin,cron      
     &,Wup,isubL,ihop1,ihop2,lsitex,lsitey,lsite,lsub,ip,im
     &,ldx,ldy,gqp,ginqp
     &,gdelqp,agqp,agdelqp,aqp,faiqpL,faiqpR,IPIVqp,WORKqp,hqp,hinqp
     &,wfin,delcro,sq,cdens,hspin
     &,w,x,cof,subh2,subkin1,subkin2,subkin3,subhdble
     &,subssnn1,subssnn2,subssnn3
     &,hnij,subnk,subsq,subnq,hnk,hsq,hnq
     &,lkx,lky,lknum)



      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 UPLO,JOBZ
      
      
     
     
     
      dimension fai(nn,2*m,Lmax)
      dimension deigen(Lmax)
      dimension WORK2(3*Lmax-1),hgeed(Lmax,Lmax),hingeed(Lmax,Lmax)
      
      dimension h(Lmax,Lmax),hin(Lmax,Lmax)
      dimension cron(2*nn,2*nn)
      
      dimension Wup(Lmax),isubL(Lmax)




c size depend !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
c shift site numbers cccccccccccccccccccccccccccc
      dimension ihop1(nn,4),ihop2(nn,2)
      dimension lsitex(nn),lsitey(nn)
      dimension lsitenum(10,10)
      dimension lsub(nn)
      dimension ip(10),im(10)
      dimension ldx(nn),ldy(nn)
      dimension lkx(22),lky(22),lknum(-2:3,0:3)
ccccccccccccccccccccccccccccccccccccccccccccccccc
c end size depend !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      





c engv dimension !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension gqp(2*nn,2*nn),ginqp(Lmax)
      dimension gdelqp(2*nn,2*nn)
      dimension agqp(2*nn,2*nn),agdelqp(2*nn,2*nn),aqp(2*m,2*m)
      
      dimension faiqpL(2*nn,2*m,Lmax),faiqpR(2*nn,2*m,lqp)
      dimension IPIVqp(2*m),WORKqp(2*m)
      dimension hqp(Lmax,Lmax),hinqp(Lmax,Lmax)
      dimension wfin(Lmax,Lmax),delcro(2*nn,2*nn)
      
      dimension sq(nn,nn),cdens(nn,nn)
      dimension hspin(2*nn,2*nn,lqp),w(lqp),x(lqp),cof(lqp)
      
      
      dimension subh2(Lmax,Lmax)
      dimension subkin1(Lmax,Lmax,2),subkin2(Lmax,Lmax,2)
      dimension subkin3(Lmax,Lmax,2),subhdble(Lmax,Lmax,2)

      
      dimension subssnn1(Lmax,Lmax,2),subssnn2(Lmax,Lmax,2)
      dimension subssnn3(Lmax,Lmax,2)
      

c additional alignment ---------------------------------------
      dimension hnij(2*nn,2*nn)

      dimension subnk(22,Lmax,Lmax)
      dimension subsq(22,Lmax,Lmax)
      dimension subnq(22,Lmax,Lmax)
      
      
      dimension hnk(22),hsq(22),hnq(2)
c end engv dimension !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




      
      data ITYPE/1/,JOBZ/'N'/,UPLO/'U'/
      if(iqpon.eq.0) then
      
      
      open(unit=100,file="t2_1_Eng_vs_Engv2_u6_m.txt"
     &,status='unknown')
      open(unit=101,file="t2_1_Eng_u6_m.txt"
     &,status='unknown')
      open(unit=102,file="t2_1_Dble_u6_m.txt"
     &,status='unknown')
c nearest neighbor spin ccccccccccccccccccccccccccccc
      open(unit=103,file="t2_1_nnspin_sq_u6_m.txt"
     &,status='unknown')
      open(unit=104,file="t2_1_nnspin_diag_u6_m.txt"
     &,status='unknown')
      open(unit=105,file="t2_1_nnspin_tot_u6_m.txt"
     &,status='unknown')
c kenetic energy cccccccccccccccccccccccccccccccccccccccccc
      open(unit=106,file="t2_1_nnhop_sq_u6_m.txt"
     &,status='unknown')
      open(unit=107,file="t2_1_nnhop_diag_u6_m.txt"
     &,status='unknown')
      open(unit=108,file="t2_1_nnhop_tot_u6_m.txt"
     &,status='unknown')
     
     
c sq cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(unit=201,file="t2_1_Sq_1_u6_m.txt"
     &,status='unknown')
      open(unit=202,file="t2_1_Sq_2_u6_m.txt"
     &,status='unknown')
      open(unit=203,file="t2_1_Sq_3_u6_m.txt"
     &,status='unknown')
      open(unit=204,file="t2_1_Sq_4_u6_m.txt"
     &,status='unknown')
      open(unit=205,file="t2_1_Sq_5_u6_m.txt"
     &,status='unknown')
      open(unit=206,file="t2_1_Sq_6_u6_m.txt"
     &,status='unknown')
      open(unit=207,file="t2_1_Sq_7_u6_m.txt"
     &,status='unknown')
      open(unit=208,file="t2_1_Sq_8_u6_m.txt"
     &,status='unknown')
      open(unit=209,file="t2_1_Sq_9_u6_m.txt"
     &,status='unknown')
      open(unit=210,file="t2_1_Sq_10_u6_m.txt"
     &,status='unknown')
      open(unit=211,file="t2_1_Sq_11_u6_m.txt"
     &,status='unknown')
      open(unit=212,file="t2_1_Sq_12_u6_m.txt"
     &,status='unknown')
      open(unit=213,file="t2_1_Sq_13_u6_m.txt"
     &,status='unknown')
      open(unit=214,file="t2_1_Sq_14_u6_m.txt"
     &,status='unknown')
      open(unit=215,file="t2_1_Sq_15_u6_m.txt"
     &,status='unknown')
      open(unit=216,file="t2_1_Sq_16_u6_m.txt"
     &,status='unknown')
      open(unit=217,file="t2_1_Sq_17_u6_m.txt"
     &,status='unknown')
      open(unit=218,file="t2_1_Sq_18_u6_m.txt"
     &,status='unknown')
      open(unit=219,file="t2_1_Sq_19_u6_m.txt"
     &,status='unknown')
      open(unit=220,file="t2_1_Sq_20_u6_m.txt"
     &,status='unknown')
      open(unit=221,file="t2_1_Sq_21_u6_m.txt"
     &,status='unknown')
      open(unit=222,file="t2_1_Sq_22_u6_m.txt"
     &,status='unknown')

c nq cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(unit=301,file="t2_1_Nq_1_u6_m.txt"
     &,status='unknown')
      open(unit=302,file="t2_1_Nq_2_u6_m.txt"
     &,status='unknown')
      open(unit=303,file="t2_1_Nq_3_u6_m.txt"
     &,status='unknown')
      open(unit=304,file="t2_1_Nq_4_u6_m.txt"
     &,status='unknown')
      open(unit=305,file="t2_1_Nq_5_u6_m.txt"
     &,status='unknown')
      open(unit=306,file="t2_1_Nq_6_u6_m.txt"
     &,status='unknown')
      open(unit=307,file="t2_1_Nq_7_u6_m.txt"
     &,status='unknown')
      open(unit=308,file="t2_1_Nq_8_u6_m.txt"
     &,status='unknown')
      open(unit=309,file="t2_1_Nq_9_u6_m.txt"
     &,status='unknown')
      open(unit=310,file="t2_1_Nq_10_u6_m.txt"
     &,status='unknown')
      open(unit=311,file="t2_1_Nq_11_u6_m.txt"
     &,status='unknown')
      open(unit=312,file="t2_1_Nq_12_u6_m.txt"
     &,status='unknown')
      open(unit=313,file="t2_1_Nq_13_u6_m.txt"
     &,status='unknown')
      open(unit=314,file="t2_1_Nq_14_u6_m.txt"
     &,status='unknown')
      open(unit=315,file="t2_1_Nq_15_u6_m.txt"
     &,status='unknown')
      open(unit=316,file="t2_1_Nq_16_u6_m.txt"
     &,status='unknown')
      open(unit=317,file="t2_1_Nq_17_u6_m.txt"
     &,status='unknown')
      open(unit=318,file="t2_1_Nq_18_u6_m.txt"
     &,status='unknown')
      open(unit=319,file="t2_1_Nq_19_u6_m.txt"
     &,status='unknown')
      open(unit=320,file="t2_1_Nq_20_u6_m.txt"
     &,status='unknown')
      open(unit=321,file="t2_1_Nq_21_u6_m.txt"
     &,status='unknown')
      open(unit=322,file="t2_1_Nq_22_u6_m.txt"
     &,status='unknown')

c nk cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(unit=401,file="t2_1_nk_1_u6_m.txt"
     &,status='unknown')
      open(unit=402,file="t2_1_nk_2_u6_m.txt"
     &,status='unknown')
      open(unit=403,file="t2_1_nk_3_u6_m.txt"
     &,status='unknown')
      open(unit=404,file="t2_1_nk_4_u6_m.txt"
     &,status='unknown')
      open(unit=405,file="t2_1_nk_5_u6_m.txt"
     &,status='unknown')
      open(unit=406,file="t2_1_nk_6_u6_m.txt"
     &,status='unknown')
      open(unit=407,file="t2_1_nk_7_u6_m.txt"
     &,status='unknown')
      open(unit=408,file="t2_1_nk_8_u6_m.txt"
     &,status='unknown')
      open(unit=409,file="t2_1_nk_9_u6_m.txt"
     &,status='unknown')
      open(unit=410,file="t2_1_nk_10_u6_m.txt"
     &,status='unknown')
      open(unit=411,file="t2_1_nk_11_u6_m.txt"
     &,status='unknown')
      open(unit=412,file="t2_1_nk_12_u6_m.txt"
     &,status='unknown')
      open(unit=413,file="t2_1_nk_13_u6_m.txt"
     &,status='unknown')
      open(unit=414,file="t2_1_nk_14_u6_m.txt"
     &,status='unknown')
      open(unit=415,file="t2_1_nk_15_u6_m.txt"
     &,status='unknown')
      open(unit=416,file="t2_1_nk_16_u6_m.txt"
     &,status='unknown')
      open(unit=417,file="t2_1_nk_17_u6_m.txt"
     &,status='unknown')
      open(unit=418,file="t2_1_nk_18_u6_m.txt"
     &,status='unknown')
      open(unit=419,file="t2_1_nk_19_u6_m.txt"
     &,status='unknown')
      open(unit=420,file="t2_1_nk_20_u6_m.txt"
     &,status='unknown')
      open(unit=421,file="t2_1_nk_21_u6_m.txt"
     &,status='unknown')
      open(unit=422,file="t2_1_nk_22_u6_m.txt"
     &,status='unknown')
      
      
      
      
      else if(ispin.eq.0) then
      open(unit=100,file="QP_S0_t2_1_Eng_vs_Engv2_u6_m.txt"
     &,status='unknown')
      open(unit=101,file="QP_S0_t2_1_Eng_u6_m.txt"
     &,status='unknown')
      open(unit=102,file="QP_S0_t2_1_Dble_u6_m.txt"
     &,status='unknown')
c nearest neighbor spin ccccccccccccccccccccccccccccc
      open(unit=103,file="QP_S0_t2_1_nnspin_sq_u6_m.txt"
     &,status='unknown')
      open(unit=104,file="QP_S0_t2_1_nnspin_diag_u6_m.txt"
     &,status='unknown')
      open(unit=105,file="QP_S0_t2_1_nnspin_tot_u6_m.txt"
     &,status='unknown')
c kenetic energy cccccccccccccccccccccccccccccccccccccccccc
      open(unit=106,file="QP_S0_t2_1_nnhop_sq_u6_m.txt"
     &,status='unknown')
      open(unit=107,file="QP_S0_t2_1_nnhop_diag_u6_m.txt"
     &,status='unknown')
      open(unit=108,file="QP_S0_t2_1_nnhop_tot_u6_m.txt"
     &,status='unknown')
     

c sq cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(unit=201,file="QP_S0_t2_1_Sq_1_u6_m.txt"
     &,status='unknown')
      open(unit=202,file="QP_S0_t2_1_Sq_2_u6_m.txt"
     &,status='unknown')
      open(unit=203,file="QP_S0_t2_1_Sq_3_u6_m.txt"
     &,status='unknown')
      open(unit=204,file="QP_S0_t2_1_Sq_4_u6_m.txt"
     &,status='unknown')
      open(unit=205,file="QP_S0_t2_1_Sq_5_u6_m.txt"
     &,status='unknown')
      open(unit=206,file="QP_S0_t2_1_Sq_6_u6_m.txt"
     &,status='unknown')
      open(unit=207,file="QP_S0_t2_1_Sq_7_u6_m.txt"
     &,status='unknown')
      open(unit=208,file="QP_S0_t2_1_Sq_8_u6_m.txt"
     &,status='unknown')
      open(unit=209,file="QP_S0_t2_1_Sq_9_u6_m.txt"
     &,status='unknown')
      open(unit=210,file="QP_S0_t2_1_Sq_10_u6_m.txt"
     &,status='unknown')
      open(unit=211,file="QP_S0_t2_1_Sq_11_u6_m.txt"
     &,status='unknown')
      open(unit=212,file="QP_S0_t2_1_Sq_12_u6_m.txt"
     &,status='unknown')
      open(unit=213,file="QP_S0_t2_1_Sq_13_u6_m.txt"
     &,status='unknown')
      open(unit=214,file="QP_S0_t2_1_Sq_14_u6_m.txt"
     &,status='unknown')
      open(unit=215,file="QP_S0_t2_1_Sq_15_u6_m.txt"
     &,status='unknown')
      open(unit=216,file="QP_S0_t2_1_Sq_16_u6_m.txt"
     &,status='unknown')
      open(unit=217,file="QP_S0_t2_1_Sq_17_u6_m.txt"
     &,status='unknown')
      open(unit=218,file="QP_S0_t2_1_Sq_18_u6_m.txt"
     &,status='unknown')
      open(unit=219,file="QP_S0_t2_1_Sq_19_u6_m.txt"
     &,status='unknown')
      open(unit=220,file="QP_S0_t2_1_Sq_20_u6_m.txt"
     &,status='unknown')
      open(unit=221,file="QP_S0_t2_1_Sq_21_u6_m.txt"
     &,status='unknown')
      open(unit=222,file="QP_S0_t2_1_Sq_22_u6_m.txt"
     &,status='unknown')

c nq cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(unit=301,file="QP_S0_t2_1_Nq_1_u6_m.txt"
     &,status='unknown')
      open(unit=302,file="QP_S0_t2_1_Nq_2_u6_m.txt"
     &,status='unknown')
      open(unit=303,file="QP_S0_t2_1_Nq_3_u6_m.txt"
     &,status='unknown')
      open(unit=304,file="QP_S0_t2_1_Nq_4_u6_m.txt"
     &,status='unknown')
      open(unit=305,file="QP_S0_t2_1_Nq_5_u6_m.txt"
     &,status='unknown')
      open(unit=306,file="QP_S0_t2_1_Nq_6_u6_m.txt"
     &,status='unknown')
      open(unit=307,file="QP_S0_t2_1_Nq_7_u6_m.txt"
     &,status='unknown')
      open(unit=308,file="QP_S0_t2_1_Nq_8_u6_m.txt"
     &,status='unknown')
      open(unit=309,file="QP_S0_t2_1_Nq_9_u6_m.txt"
     &,status='unknown')
      open(unit=310,file="QP_S0_t2_1_Nq_10_u6_m.txt"
     &,status='unknown')
      open(unit=311,file="QP_S0_t2_1_Nq_11_u6_m.txt"
     &,status='unknown')
      open(unit=312,file="QP_S0_t2_1_Nq_12_u6_m.txt"
     &,status='unknown')
      open(unit=313,file="QP_S0_t2_1_Nq_13_u6_m.txt"
     &,status='unknown')
      open(unit=314,file="QP_S0_t2_1_Nq_14_u6_m.txt"
     &,status='unknown')
      open(unit=315,file="QP_S0_t2_1_Nq_15_u6_m.txt"
     &,status='unknown')
      open(unit=316,file="QP_S0_t2_1_Nq_16_u6_m.txt"
     &,status='unknown')
      open(unit=317,file="QP_S0_t2_1_Nq_17_u6_m.txt"
     &,status='unknown')
      open(unit=318,file="QP_S0_t2_1_Nq_18_u6_m.txt"
     &,status='unknown')
      open(unit=319,file="QP_S0_t2_1_Nq_19_u6_m.txt"
     &,status='unknown')
      open(unit=320,file="QP_S0_t2_1_Nq_20_u6_m.txt"
     &,status='unknown')
      open(unit=321,file="QP_S0_t2_1_Nq_21_u6_m.txt"
     &,status='unknown')
      open(unit=322,file="QP_S0_t2_1_Nq_22_u6_m.txt"
     &,status='unknown')

c nk cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(unit=401,file="QP_S0_t2_1_nk_1_u6_m.txt"
     &,status='unknown')
      open(unit=402,file="QP_S0_t2_1_nk_2_u6_m.txt"
     &,status='unknown')
      open(unit=403,file="QP_S0_t2_1_nk_3_u6_m.txt"
     &,status='unknown')
      open(unit=404,file="QP_S0_t2_1_nk_4_u6_m.txt"
     &,status='unknown')
      open(unit=405,file="QP_S0_t2_1_nk_5_u6_m.txt"
     &,status='unknown')
      open(unit=406,file="QP_S0_t2_1_nk_6_u6_m.txt"
     &,status='unknown')
      open(unit=407,file="QP_S0_t2_1_nk_7_u6_m.txt"
     &,status='unknown')
      open(unit=408,file="QP_S0_t2_1_nk_8_u6_m.txt"
     &,status='unknown')
      open(unit=409,file="QP_S0_t2_1_nk_9_u6_m.txt"
     &,status='unknown')
      open(unit=410,file="QP_S0_t2_1_nk_10_u6_m.txt"
     &,status='unknown')
      open(unit=411,file="QP_S0_t2_1_nk_11_u6_m.txt"
     &,status='unknown')
      open(unit=412,file="QP_S0_t2_1_nk_12_u6_m.txt"
     &,status='unknown')
      open(unit=413,file="QP_S0_t2_1_nk_13_u6_m.txt"
     &,status='unknown')
      open(unit=414,file="QP_S0_t2_1_nk_14_u6_m.txt"
     &,status='unknown')
      open(unit=415,file="QP_S0_t2_1_nk_15_u6_m.txt"
     &,status='unknown')
      open(unit=416,file="QP_S0_t2_1_nk_16_u6_m.txt"
     &,status='unknown')
      open(unit=417,file="QP_S0_t2_1_nk_17_u6_m.txt"
     &,status='unknown')
      open(unit=418,file="QP_S0_t2_1_nk_18_u6_m.txt"
     &,status='unknown')
      open(unit=419,file="QP_S0_t2_1_nk_19_u6_m.txt"
     &,status='unknown')
      open(unit=420,file="QP_S0_t2_1_nk_20_u6_m.txt"
     &,status='unknown')
      open(unit=421,file="QP_S0_t2_1_nk_21_u6_m.txt"
     &,status='unknown')
      open(unit=422,file="QP_S0_t2_1_nk_22_u6_m.txt"
     &,status='unknown')
      


     
     
     
      else if(ispin.eq.1) then
      open(unit=100,file="QP_S1_t2_1_Eng_vs_Engv2_u6_m.txt"
     &,status='unknown')
      open(unit=101,file="QP_S1_t2_1_Eng_u6_m.txt"
     &,status='unknown')
      open(unit=102,file="QP_S1_t2_1_Dble_u6_m.txt"
     &,status='unknown')
c nearest neighbor spin ccccccccccccccccccccccccccccc
      open(unit=103,file="QP_S1_t2_1_nnspin_sq_u6_m.txt"
     &,status='unknown')
      open(unit=104,file="QP_S1_t2_1_nnspin_diag_u6_m.txt"
     &,status='unknown')
      open(unit=105,file="QP_S1_t2_1_nnspin_tot_u6_m.txt"
     &,status='unknown')
c kenetic energy cccccccccccccccccccccccccccccccccccccccccc
      open(unit=106,file="QP_S1_t2_1_nnhop_sq_u6_m.txt"
     &,status='unknown')
      open(unit=107,file="QP_S1_t2_1_nnhop_diag_u6_m.txt"
     &,status='unknown')
      open(unit=108,file="QP_S1_t2_1_nnhop_tot_u6_m.txt"
     &,status='unknown')
     
     
c sq cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(unit=201,file="QP_S1_t2_1_Sq_1_u6_m.txt"
     &,status='unknown')
      open(unit=202,file="QP_S1_t2_1_Sq_2_u6_m.txt"
     &,status='unknown')
      open(unit=203,file="QP_S1_t2_1_Sq_3_u6_m.txt"
     &,status='unknown')
      open(unit=204,file="QP_S1_t2_1_Sq_4_u6_m.txt"
     &,status='unknown')
      open(unit=205,file="QP_S1_t2_1_Sq_5_u6_m.txt"
     &,status='unknown')
      open(unit=206,file="QP_S1_t2_1_Sq_6_u6_m.txt"
     &,status='unknown')
      open(unit=207,file="QP_S1_t2_1_Sq_7_u6_m.txt"
     &,status='unknown')
      open(unit=208,file="QP_S1_t2_1_Sq_8_u6_m.txt"
     &,status='unknown')
      open(unit=209,file="QP_S1_t2_1_Sq_9_u6_m.txt"
     &,status='unknown')
      open(unit=210,file="QP_S1_t2_1_Sq_10_u6_m.txt"
     &,status='unknown')
      open(unit=211,file="QP_S1_t2_1_Sq_11_u6_m.txt"
     &,status='unknown')
      open(unit=212,file="QP_S1_t2_1_Sq_12_u6_m.txt"
     &,status='unknown')
      open(unit=213,file="QP_S1_t2_1_Sq_13_u6_m.txt"
     &,status='unknown')
      open(unit=214,file="QP_S1_t2_1_Sq_14_u6_m.txt"
     &,status='unknown')
      open(unit=215,file="QP_S1_t2_1_Sq_15_u6_m.txt"
     &,status='unknown')
      open(unit=216,file="QP_S1_t2_1_Sq_16_u6_m.txt"
     &,status='unknown')
      open(unit=217,file="QP_S1_t2_1_Sq_17_u6_m.txt"
     &,status='unknown')
      open(unit=218,file="QP_S1_t2_1_Sq_18_u6_m.txt"
     &,status='unknown')
      open(unit=219,file="QP_S1_t2_1_Sq_19_u6_m.txt"
     &,status='unknown')
      open(unit=220,file="QP_S1_t2_1_Sq_20_u6_m.txt"
     &,status='unknown')
      open(unit=221,file="QP_S1_t2_1_Sq_21_u6_m.txt"
     &,status='unknown')
      open(unit=222,file="QP_S1_t2_1_Sq_22_u6_m.txt"
     &,status='unknown')

c nq cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(unit=301,file="QP_S1_t2_1_Nq_1_u6_m.txt"
     &,status='unknown')
      open(unit=302,file="QP_S1_t2_1_Nq_2_u6_m.txt"
     &,status='unknown')
      open(unit=303,file="QP_S1_t2_1_Nq_3_u6_m.txt"
     &,status='unknown')
      open(unit=304,file="QP_S1_t2_1_Nq_4_u6_m.txt"
     &,status='unknown')
      open(unit=305,file="QP_S1_t2_1_Nq_5_u6_m.txt"
     &,status='unknown')
      open(unit=306,file="QP_S1_t2_1_Nq_6_u6_m.txt"
     &,status='unknown')
      open(unit=307,file="QP_S1_t2_1_Nq_7_u6_m.txt"
     &,status='unknown')
      open(unit=308,file="QP_S1_t2_1_Nq_8_u6_m.txt"
     &,status='unknown')
      open(unit=309,file="QP_S1_t2_1_Nq_9_u6_m.txt"
     &,status='unknown')
      open(unit=310,file="QP_S1_t2_1_Nq_10_u6_m.txt"
     &,status='unknown')
      open(unit=311,file="QP_S1_t2_1_Nq_11_u6_m.txt"
     &,status='unknown')
      open(unit=312,file="QP_S1_t2_1_Nq_12_u6_m.txt"
     &,status='unknown')
      open(unit=313,file="QP_S1_t2_1_Nq_13_u6_m.txt"
     &,status='unknown')
      open(unit=314,file="QP_S1_t2_1_Nq_14_u6_m.txt"
     &,status='unknown')
      open(unit=315,file="QP_S1_t2_1_Nq_15_u6_m.txt"
     &,status='unknown')
      open(unit=316,file="QP_S1_t2_1_Nq_16_u6_m.txt"
     &,status='unknown')
      open(unit=317,file="QP_S1_t2_1_Nq_17_u6_m.txt"
     &,status='unknown')
      open(unit=318,file="QP_S1_t2_1_Nq_18_u6_m.txt"
     &,status='unknown')
      open(unit=319,file="QP_S1_t2_1_Nq_19_u6_m.txt"
     &,status='unknown')
      open(unit=320,file="QP_S1_t2_1_Nq_20_u6_m.txt"
     &,status='unknown')
      open(unit=321,file="QP_S1_t2_1_Nq_21_u6_m.txt"
     &,status='unknown')
      open(unit=322,file="QP_S1_t2_1_Nq_22_u6_m.txt"
     &,status='unknown')

c nk cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(unit=401,file="QP_S1_t2_1_nk_1_u6_m.txt"
     &,status='unknown')
      open(unit=402,file="QP_S1_t2_1_nk_2_u6_m.txt"
     &,status='unknown')
      open(unit=403,file="QP_S1_t2_1_nk_3_u6_m.txt"
     &,status='unknown')
      open(unit=404,file="QP_S1_t2_1_nk_4_u6_m.txt"
     &,status='unknown')
      open(unit=405,file="QP_S1_t2_1_nk_5_u6_m.txt"
     &,status='unknown')
      open(unit=406,file="QP_S1_t2_1_nk_6_u6_m.txt"
     &,status='unknown')
      open(unit=407,file="QP_S1_t2_1_nk_7_u6_m.txt"
     &,status='unknown')
      open(unit=408,file="QP_S1_t2_1_nk_8_u6_m.txt"
     &,status='unknown')
      open(unit=409,file="QP_S1_t2_1_nk_9_u6_m.txt"
     &,status='unknown')
      open(unit=410,file="QP_S1_t2_1_nk_10_u6_m.txt"
     &,status='unknown')
      open(unit=411,file="QP_S1_t2_1_nk_11_u6_m.txt"
     &,status='unknown')
      open(unit=412,file="QP_S1_t2_1_nk_12_u6_m.txt"
     &,status='unknown')
      open(unit=413,file="QP_S1_t2_1_nk_13_u6_m.txt"
     &,status='unknown')
      open(unit=414,file="QP_S1_t2_1_nk_14_u6_m.txt"
     &,status='unknown')
      open(unit=415,file="QP_S1_t2_1_nk_15_u6_m.txt"
     &,status='unknown')
      open(unit=416,file="QP_S1_t2_1_nk_16_u6_m.txt"
     &,status='unknown')
      open(unit=417,file="QP_S1_t2_1_nk_17_u6_m.txt"
     &,status='unknown')
      open(unit=418,file="QP_S1_t2_1_nk_18_u6_m.txt"
     &,status='unknown')
      open(unit=419,file="QP_S1_t2_1_nk_19_u6_m.txt"
     &,status='unknown')
      open(unit=420,file="QP_S1_t2_1_nk_20_u6_m.txt"
     &,status='unknown')
      open(unit=421,file="QP_S1_t2_1_nk_21_u6_m.txt"
     &,status='unknown')
      open(unit=422,file="QP_S1_t2_1_nk_22_u6_m.txt"
     &,status='unknown')
      
      
      
      
      end if




      m2=2*m

      pai=dacos(-1.0d0)
      


      if(iqpon.eq.1) then
         write(6,*) '+QP ON','    TOTAL S=',ispin
      else
         write(6,*) '+QP OFF'
      end if


!$omp parallel do
      do ll=1,Lmax
      do jm=1,m
      do i=1,nn
         faiqpL(i,jm,ll)=fai(i,jm,ll)
         faiqpL(i+nn,jm,ll)=0.0d0
         faiqpL(i,jm+m,ll)=0.0d0
         faiqpL(i+nn,jm+m,ll)=fai(i,jm+m,ll)
      end do
      end do
      end do


      
      call spinrt(nn,m,hspin,lqp,cof,w,x,ispin)
      if(iqpon.eq.0) cof(1)=1.0d0
      hqp=0.0d0
      hinqp=0.0d0
      
      subhdble=0.0d0
      subh2=0.0d0
      
      subnk=0.0d0
      subsq=0.0d0
      subnq=0.0d0
      subkin1=0.0d0
      subkin2=0.0d0
      subkin3=0.0d0
      subssnn1=0.0d0
      subssnn2=0.0d0
      subssnn3=0.0d0


c initial eng ---------------------------------------------------
      do 1811 lr=1,iLdim
      
      faiqpR=0.0d0
      if(iqpon.eq.1) then
c QP Rotation ----------------------------------------------------
      lqpend=lqp
!$omp parallel do
      do iqp=1,lqp
      do j=1,m*2
         do i=1,nn*2
            do l=1,nn*2
      faiqpR(i,j,iqp)=faiqpR(i,j,iqp)+hspin(i,l,iqp)*faiqpL(l,j,lr)  
            end do  
         end do
      end do
      end do
c -------------------------------------------------------------
      else
c Normal ----------------------------------------------------
      lqpend=1
      iqp=1
      do j=1,m*2
         do i=1,nn*2
      faiqpR(i,j,iqp)=faiqpL(i,j,lr)  
         end do
      end do
c -------------------------------------------------------------
      end if
         do 1911 ll=1,lr
            do 1711 iqp=1,lqpend
      call mkgreenqp(ll,nn,m,Lmax,lqp,iqp,IPIVqp,WORKqp
     &,aqp,faiqpL,faiqpR,gqp,ginqp)
      
      do i=1,nn*2
         do j=1,nn*2
      agqp(i,j)=cron(j,i)-gqp(j,i)
         end do
      end do

!$omp parallel do
      do j=1,nn
         do i=1,nn
      gdelqp(i,j)=0.0d0
      gdelqp(i+nn,j+nn)=0.0d0
      gdelqp(i+nn,j)=0.0d0
      gdelqp(i,j+nn)=0.0d0
      do k=1,4
      gdelqp(i,j)=gdelqp(i,j)
     &-t1*gqp(i,ihop1(j,k))
      end do
      do k=1,2
      gdelqp(i,j)=gdelqp(i,j)
     &-t2*gqp(i,ihop2(j,k))
      end do
      
      do k=1,4
      gdelqp(i+nn,j+nn)=gdelqp(i+nn,j+nn)
     &-t1*gqp(i+nn,ihop1(j,k)+nn)
      end do
       do k=1,2
      gdelqp(i+nn,j+nn)=gdelqp(i+nn,j+nn)
     &-t2*gqp(i+nn,ihop2(j,k)+nn)
      end do
      
      do k=1,4
      gdelqp(i+nn,j)=gdelqp(i+nn,j)
     &-t1*gqp(i+nn,ihop1(j,k))
      end do
       do k=1,2
      gdelqp(i+nn,j)=gdelqp(i+nn,j)
     &-t2*gqp(i+nn,ihop2(j,k))
      end do
      
      do k=1,4
      gdelqp(i,j+nn)=gdelqp(i,j+nn)
     &-t1*gqp(i,ihop1(j,k)+nn)
      end do
       do k=1,2
      gdelqp(i,j+nn)=gdelqp(i,j+nn)
     &-t2*gqp(i,ihop2(j,k)+nn)
      end do
      
         end do
      end do
      
      do i=1,nn*2
         do j=1,nn*2
      agdelqp(i,j)=delcro(j,i)-gdelqp(j,i)
         end do
      end do







      call mkHqp(ll,lr,nn,u,t1,t2,Lmax
     &,ihop1,ihop2,h,hin,gqp,ginqp)
      hqp(ll,lr)=hqp(ll,lr)+h(ll,lr)*cof(iqp)
      hinqp(ll,lr)=hinqp(ll,lr)+hin(ll,lr)*cof(iqp)
      
c      write(6,*) hqp(ll,lr)
      
      do j=1,nn
c double occupancy 
      subhdble(ll,lr,1)=subhdble(ll,lr,1)
     &+(gqp(j,j)*gqp(j+nn,j+nn)
     &+gqp(j,j+nn)*agqp(j,j+nn))*hin(ll,lr)*cof(iqp)
      end do
      
      
            
      
      h2=0.0d0
      hnij=0.0d0
      cdens=0.0d0
      sq=0.0d0

!$omp parallel do reduction(+:h2) private(i1,i2,j1,j2,awick)     
      do 1951 j=1,nn
      do 1952 i=1,nn
c H*H term -----------------------------------------------------
c make hk^2 term ccccccccccccccccccccccccccccccccccccccccc      
      h2=h2+(
     &+(gdelqp(i,i)+gdelqp(i+nn,i+nn))
     &*(gdelqp(j,j)+gdelqp(j+nn,j+nn))
cc up spin ccccccccc
     &+gdelqp(i,j)*agdelqp(i,j)
     &+gdelqp(i,j+nn)*agdelqp(i,j+nn)
cc down spin cccccccc
     &+gdelqp(i+nn,j)*agdelqp(i+nn,j)
     &+gdelqp(i+nn,j+nn)*agdelqp(i+nn,j+nn)
     &)
     
c make Hu^2 term ccccccccccccccccccccccccccccccccccccc
      h2=h2
c a
     &+(u**2)*((gqp(i,i)*gqp(i+nn,i+nn)
     &+gqp(i,i+nn)*agqp(i,i+nn))
     &*(gqp(j,j)*gqp(j+nn,j+nn)
     &+gqp(j,j+nn)*agqp(j,j+nn))
c b
     &+(gqp(i,i)*gqp(i+nn,j)
     &+gqp(i,j)*agqp(i,i+nn))
     &*(agqp(i+nn,j)*gqp(j+nn,j+nn)
     &-agqp(i+nn,j+nn)*gqp(j,j+nn))
c c
     &+(gqp(i,i)*gqp(i+nn,j+nn)
     &+gqp(i,j+nn)*agqp(i,i+nn))
     &*(agqp(i+nn,j)*agqp(j,j+nn)
     &+agqp(i+nn,j+nn)*gqp(j,j))
c d    
     &+(-gqp(i,i+nn)*gqp(i+nn,j)
     &+gqp(i,j)*gqp(i+nn,i+nn))
     &*(agqp(i,j)*gqp(j+nn,j+nn)
     &-agqp(i,j+nn)*gqp(j,j+nn))
c e
     &+(-gqp(i,i+nn)*gqp(i+nn,j+nn)
     &+gqp(i,j+nn)*gqp(i+nn,i+nn))
     &*(agqp(i,j)*agqp(j,j+nn)
     &+agqp(i,j+nn)*gqp(j,j))
c f
     &+(-gqp(i,j)*gqp(i+nn,j+nn)
     &+gqp(i,j+nn)*gqp(i+nn,j))
     &*(-agqp(i,j)*agqp(i+nn,j+nn)
     &+agqp(i,j+nn)*agqp(i+nn,j))
     &)
c make Hk*Hu term ccccccccccccccccccccccccccccccccccccccccccccc
c 1
      h2=h2
     &+u*(
     &(gdelqp(i,i)+gdelqp(i+nn,i+nn))
     &*(gqp(j,j)*gqp(j+nn,j+nn)
     &+gqp(j,j+nn)*agqp(j,j+nn))
c 2
     &+gqp(i,j)*(
     &agdelqp(i,j)*gqp(j+nn,j+nn)
     &-agdelqp(i,j+nn)*gqp(j,j+nn))
c 3
     &+gqp(i+nn,j)*(
     &+agdelqp(i+nn,j)*gqp(j+nn,j+nn)
     &-agdelqp(i+nn,j+nn)*gqp(j,j+nn))
c 4
     &+gqp(i,j+nn)*(
     &agdelqp(i,j)*agqp(j,j+nn)
     &+agdelqp(i,j+nn)*gqp(j,j))
c 5
     &+gqp(i+nn,j+nn)*(
     &agdelqp(i+nn,j)*agqp(j,j+nn)
     &+agdelqp(i+nn,j+nn)*gqp(j,j))
     &)
c make Hu*Hk term ccccccccccccccccccccccccccccccccccccccccccccc
c 1
      h2=h2
     &+u*(
     &(gqp(i,i)*gqp(i+nn,i+nn)
     &+gqp(i,i+nn)*agqp(i,i+nn))
     &*(gdelqp(j,j)+gdelqp(j+nn,j+nn))
c 2
     &+agqp(i+nn,j)*(
     &gqp(i,i)*gdelqp(i+nn,j)
     &+gdelqp(i,j)*agqp(i,i+nn))
c 3
     &+agqp(i+nn,j+nn)*(
     &gqp(i,i)*gdelqp(i+nn,j+nn)
     &+gdelqp(i,j+nn)*agqp(i,i+nn))
c 4
     &+agqp(i,j)*(
     &-gqp(i,i+nn)*gdelqp(i+nn,j)
     &+gdelqp(i,j)*gqp(i+nn,i+nn))
c 5
     &+agqp(i,j+nn)*(
     &-gqp(i,i+nn)*gdelqp(i+nn,j+nn)
     &+gdelqp(i,j+nn)*gqp(i+nn,i+nn))
     &)
c end H*H term -----------------------------------------------------




c nk ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      hnij(i,j)=hnij(i,j)
     &+gqp(i,j)+gqp(i+nn,j+nn)
     &+gqp(j,i)+gqp(j+nn,i+nn)

c cdens (nq) ccccccccccccccccccccccccccccccccccccccccccccc
      cdens(i,j)=cdens(i,j)
     &+gqp(i,j)*agqp(i,j)
     &+gqp(i,j+nn)*agqp(i,j+nn)
     &+gqp(i+nn,j)*agqp(i+nn,j)
     &+gqp(i+nn,j+nn)*agqp(i+nn,j+nn)

c Sq cccccccccccccccccccccccccccccccccccccccccccccccccccc     
      awick=0.0d0
      i1=i
      i2=i+nn
      j1=j+nn
      j2=j
      awick=awick
     &+0.5d0
     &*(gqp(i1,i2)*gqp(j1,j2)+gqp(i1,j2)*agqp(i2,j1))
      i1=i+nn
      i2=i
      j1=j
      j2=j+nn
      awick=awick
     &+0.5d0
     &*(gqp(i1,i2)*gqp(j1,j2)+gqp(i1,j2)*agqp(i2,j1))
      i1=i
      i2=i
      j1=j
      j2=j
      awick=awick
     &+0.25d0
     &*(gqp(i1,i2)*gqp(j1,j2)+gqp(i1,j2)*agqp(i2,j1))
      i1=i+nn
      i2=i+nn
      j1=j
      j2=j
      awick=awick
     &-0.25d0
     &*(gqp(i1,i2)*gqp(j1,j2)+gqp(i1,j2)*agqp(i2,j1))
      i1=i
      i2=i
      j1=j+nn
      j2=j+nn
      awick=awick
     &-0.25d0
     &*(gqp(i1,i2)*gqp(j1,j2)+gqp(i1,j2)*agqp(i2,j1))
      i1=i+nn
      i2=i+nn
      j1=j+nn
      j2=j+nn
      awick=awick
     &+0.25d0
     &*(gqp(i1,i2)*gqp(j1,j2)+gqp(i1,j2)*agqp(i2,j1))
      sq(i,j)=sq(i,j)+awick


1952  continue
1951  continue
c end H*H term cccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subh2(ll,lr)=subh2(ll,lr)+h2*hin(ll,lr)*cof(iqp)

      do i=1,nn
      
      do k=1,4
      subssnn1(ll,lr,1)=subssnn1(ll,lr,1)
     &+sq(i,ihop1(i,k))*hin(ll,lr)*cof(iqp)
      subkin1(ll,lr,1)=subkin1(ll,lr,1)
     &-(gqp(i,ihop1(i,k))+gqp(i+nn,ihop1(i,k)+nn))
     &*hin(ll,lr)*cof(iqp)
      end do
      do k=1,2
      subssnn1(ll,lr,1)=subssnn1(ll,lr,1)
     &+sq(i,ihop2(i,k))*hin(ll,lr)*cof(iqp)
      subkin1(ll,lr,1)=subkin1(ll,lr,1)
     &-t2*(gqp(i,ihop2(i,k))+gqp(i+nn,ihop2(i,k)+nn))
     &*hin(ll,lr)*cof(iqp)
      end do
      
      do k=1,4
      subssnn2(ll,lr,1)=subssnn2(ll,lr,1)
     &+sq(i,ihop1(i,k))*hin(ll,lr)*cof(iqp)
      subkin2(ll,lr,1)=subkin2(ll,lr,1)
     &-(gqp(i,ihop1(i,k))+gqp(i+nn,ihop1(i,k)+nn))
     &*hin(ll,lr)*cof(iqp)
      end do
      
      do k=1,2
      subssnn3(ll,lr,1)=subssnn3(ll,lr,1)
     &+sq(i,ihop2(i,k))*hin(ll,lr)*cof(iqp)
      subkin3(ll,lr,1)=subkin3(ll,lr,1)
     &-(gqp(i,ihop2(i,k))+gqp(i+nn,ihop2(i,k)+nn))
     &*hin(ll,lr)*cof(iqp)
      end do      
      end do
      
      
      
      


c nk,Nq,Sq -------------------------------------------------
!$omp parallel do private(akx,aky,jx,jy,ix,iy)
      do ik=1,22
      akx=2.0d0*pai*dble(lkx(ik))/dble(6.0d0)
      aky=2.0d0*pai*dble(lky(ik))/dble(6.0d0)
         do j=1,nn
         jx=lsitex(j)
         jy=lsitey(j)

         do i=1,nn
         ix=lsitex(i)
         iy=lsitey(i)
         
      subnk(ik,ll,lr)=subnk(ik,ll,lr)
     &+hnij(i,j)*dcos(akx*dble(ix-jx)+aky*dble(iy-jy))
     &*hin(ll,lr)*cof(iqp)
     
      subnq(ik,ll,lr)=subnq(ik,ll,lr)
     &+cdens(i,j)*dcos(akx*dble(ix-jx)+aky*dble(iy-jy))
     &*hin(ll,lr)*cof(iqp)
     
      subsq(ik,ll,lr)=subsq(ik,ll,lr)
     &+sq(i,j)*dcos(akx*dble(ix-jx)+aky*dble(iy-jy))
     &*hin(ll,lr)*cof(iqp)
     
         end do
         end do
      end do
c end nk,Nq,Sq --------------------------------------------------------      

1711        continue
1911     continue
1811  continue
c end initial eng -------------------------------------------------
      
      
      
      
      




c subspace EXT --------------------------------------------------
      
      LWORK=3*iLdim-1
      call geedv(Lmax,iLdim,label,hqp,hinqp
     &,hgeed,hingeed,deigen,Wup,eigen,WORK2,LWORK)
      pirgeng=eigen
      hinre=0.0d0
      do lr=1,iLdim
         do ll=1,lr
            if(ll.ne.lr) then
         hinre=hinre+2.0d0*hgeed(ll,1)*hgeed(lr,1)*hinqp(ll,lr)
            else
         hinre=hinre+hgeed(ll,1)*hgeed(lr,1)*hinqp(ll,lr)
            end if
         end do
      end do
      
      write(6,*) 'dimension=',iLdim,'eng=',pirgeng
      
      hh2=0.0d0
      hdble=0.0d0
      hkin1=0.0d0
      hkin2=0.0d0
      hkin3=0.0d0
      
      hssnn1=0.0d0
      hssnn2=0.0d0
      hssnn3=0.0d0
      
      
      
      
c additional alignment -------------------------------------------
      hnk=0.0d0
      hsq=0.0d0
      hnq=0.0d0
c -------------------------------------------------------------------
!$omp parallel do private(cweight)
!reduction(+:hh2,hdble,hkin1,hkin2,hkin3,hssnn1,hssnn2,hssnn3)      
      do lr=1,iLdim
         do ll=1,lr
      if(ll.ne.lr) then
         cweight=2.0d0*hgeed(ll,1)*hgeed(lr,1)
      else
         cweight=hgeed(ll,1)*hgeed(lr,1)
      end if
      hh2=hh2+subh2(ll,lr)*cweight
      hdble=hdble+subhdble(ll,lr,1)*cweight
      hkin1=hkin1+subkin1(ll,lr,1)*cweight
      hkin2=hkin2+subkin2(ll,lr,1)*cweight
      hkin3=hkin3+subkin3(ll,lr,1)*cweight
      
      
      hssnn1=hssnn1+subssnn1(ll,lr,1)*cweight
      hssnn2=hssnn2+subssnn2(ll,lr,1)*cweight
      hssnn3=hssnn3+subssnn3(ll,lr,1)*cweight
      
      
      
c additional term ------------------------------------------------      
      do k=1,22   
         hnk(k)=hnk(k)+subnk(k,ll,lr)*cweight
         hsq(k)=hsq(k)+subsq(k,ll,lr)*cweight
         hnq(k)=hnq(k)+subnq(k,ll,lr)*cweight
      end do
c --------------------------------------------------------------------
         end do
      end do


      
      hh2=hh2/hinre
      
      cmod=-u*(dble(m)-dble(nn)/4.0d0)
      engva=(hh2+2.0d0*cmod*pirgeng+cmod**2)
     &/(pirgeng**2+2.0d0*cmod*pirgeng+cmod**2)-1.0d0
      engva2=hh2/pirgeng/pirgeng-1.0d0
        
        
      hdble=hdble/hinre/dble(nn)
      
      hkin1=hkin1/hinre/dble(nn)/6.0d0
      hkin2=hkin2/hinre/dble(nn)/4.0d0
      hkin3=hkin3/hinre/dble(nn)/2.0d0


      hssnn1=hssnn1/hinre/dble(nn)/6.0d0
      hssnn2=hssnn2/hinre/dble(nn)/4.0d0
      hssnn3=hssnn3/hinre/dble(nn)/2.0d0
      
      
      
c additional term --------------------------------------------------      
      do k=1,22
      hnk(k)=hnk(k)*0.25d0/dble(nn)/hinre
      hsq(k)=hsq(k)/dble(nn)/3.0d0/hinre
      hnq(k)=hnq(k)/dble(nn)/hinre
      end do
c -----------------------------------------------------------------
      

      write(100,*) engva2,pirgeng,iLdim
      write(101,*) engva,pirgeng,iLdim
      write(102,*) engva,hdble,iLdim
      write(103,*) engva,hssnn2,iLdim
      write(104,*) engva,hssnn3,iLdim
      write(105,*) engva,hssnn1,iLdim
      write(106,*) engva,hkin2,iLdim
      write(107,*) engva,hkin3,iLdim
      write(108,*) engva,hkin1,iLdim
      
      write(6,*)
      write(6,*) 'total spin'
      write(6,*) hsq(1)*3.0d0*dble(nn)
     
      do k=1,22
         write(400+k,*) engva,hnk(k),iLdim
         write(300+k,*) engva,hnq(k),iLdim
         write(200+k,*) engva,hsq(k),iLdim
      end do
      
2222  format(a,F25.18,F25.18,F25.18,F25.18)
      return
      end
c end PIRG checker board
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
      
      
      
      












cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine shift(nn,ihop1,ihop2,lsitex,lsitey,lsitenum
     &,lsub,ip,im,ldx,ldy)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
c shift site numbers cccccccccccccccccccccccccccc
      dimension ihop1(nn,4),ihop2(nn,2)
      dimension lsitex(nn),lsitey(nn)
      dimension lsitenum(10,10)
      dimension lsub(nn)
      dimension ip(10),im(10)
      dimension ldx(nn),ldy(nn)
ccccccccccccccccccccccccccccccccccccccccccccccccc
      
      
      
      
      lsitex=0.0d0
      lsitey=0.0d0
      lsitenum=0.0d0
c size depend !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
c shift site nubers ccccccccccccccccccccccccccccccc
      if(nn.eq.36) then
      do i=1,10
         ip(i)=i+1
         im(i)=i-1
      end do
c set site num.
      inum=0
      do iy=1,6
         do ix=1,6
            inum=inum+1
            lsitenum(ix+1,iy+1)=inum
         end do
      end do
      
      do ix=1,10
      do iy=1,10
         if(lsitenum(ix,iy).ne.0) then
            lsitex(lsitenum(ix,iy))=ix
            lsitey(lsitenum(ix,iy))=iy
         end if
      end do
      end do
      
      do ix=1,6
         lsitenum(ix+1,1)=lsitenum(ix+1,7)
         lsitenum(ix+1,8)=lsitenum(ix+1,2)
      end do
      do iy=1,6
         lsitenum(1,iy+1)=lsitenum(7,iy+1)
         lsitenum(8,iy+1)=lsitenum(2,iy+1)
      end do
      
      lsitenum(1,1)=36
      lsitenum(8,1)=31
      lsitenum(1,8)=6
      lsitenum(8,8)=1
      
      
c      do j=1,10
c      do i=1,10
c      write(6,*) i,j,lsitenum(i,j)
c      end do
c      end do
c      stop
c      do i=1,nn
c         write(6,*) i,lsitex(i),lsitey(i)
c      end do
c      stop

      
      do i=1,nn
      
      ix=lsitex(i)
      iy=lsitey(i)
         ihop1(i,1)=lsitenum(ip(ix),iy)
         ihop1(i,2)=lsitenum(ix,ip(iy))
         ihop1(i,3)=lsitenum(im(ix),iy)
         ihop1(i,4)=lsitenum(ix,im(iy))
         
         ihop2(i,1)=lsitenum(ip(ix),ip(iy))
         ihop2(i,2)=lsitenum(im(ix),im(iy))
         
      end do


c      do i=1,nn
c         do k=1,4
c         write(6,*) i,ihop1(i,k)
c         end do
c         write(6,*)
c         do k=1,2
c         write(6,*) i,ihop2(i,k)
c         end do
c         write(6,*)
c      end do
c      stop


c set sublattice num.
       do isub=1,6
          i=isub
          do j=1,6
             i=ihop1(ihop1(i,3),2)
             lsub(i)=isub
          end do
       end do

c      do i=1,nn
c      write(6,*) i,lsub(i)
c      end do
c      stop



      end if      

      return
      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc











      
      
