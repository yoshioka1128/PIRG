c 6hands PIRG --------------------------------------------------
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 JOBZ  
      character*1 UPLO
      
      parameter(nn=36,m2=nn,u=9.0d0,t1=1.0d0,t2=1.0d0
     &,istdim=1,Lmax=500,limW=m2*100,limW2=(3*Lmax-1)*100
c set up -----------------------------------------------------
     &,iRGtype=0,ichmod=1,ratelim1=0.1d0,ratelim2=0.05d0
     &,delexp=5.0d-4/u,delu=0.5d0/u,limit=100000
c engv -----------------------------------------------------
     &,iengv=1,iqpon=0,ispin=0,iWF=0)

c RG type iRGtype=0: independent RG
c RG type iRGtype=1: depend RG

c change RG mode ichmod=0: no chanege
c change RG mode ichmod=1: change

c iengv=0: no engv process
c iengv=1: engv process


      dimension g(2*nn,2*nn,Lmax),gup2(2*nn,2*nn,Lmax,2)
      dimension gin(Lmax),g2(2*nn,2*nn,Lmax,2),gin2(Lmax,2)
      
      dimension h(Lmax,Lmax),hin(Lmax,Lmax)
      dimension hold(Lmax,Lmax),hinold(Lmax,Lmax)
      dimension hsub(Lmax,2),hinsub(Lmax,2)
      
      dimension fai(2*nn,m2,Lmax),faiold(2*nn,m2),fai0(2*nn,m2)
      
      dimension WORK2(limW2),hgeed(Lmax,Lmax),hingeed(Lmax,Lmax)
      dimension deigen(Lmax)
      
      dimension lurg(Lmax,Lmax),iss(2)
      dimension gsub(2*nn,2*nn)
      dimension buf(m2,2*nn)
      dimension a(m2,m2),IPIV(m2),WORK(limW)
      dimension alpha1(-1:1)
      dimension pirgmp(2),isubL(Lmax)
      dimension isite(nn),idim0(Lmax),iRGdim(Lmax)
c schmidt
      dimension stfai(2*nn,m2,Lmax),stfaiold(2*nn,m2)
      dimension aeigen(m2),cc(m2),fai2(2*nn,m2)
c kinetic2      
      dimension hk(2*nn,2*nn),hku(2*nn,2*nn)
      dimension ee(2*nn),hnqei(2*nn),hkL0(2*nn,2*nn)
      dimension hu(2*nn),hu2(2*nn,2*nn)
      
      
      
c size depend !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
c shift site numbers cccccccccccccccccccccccccccc
      dimension ihop1(nn,4),ihop2(nn,2)
      dimension lsitex(nn),lsitey(nn)
      dimension lsub(nn)
c end size depend !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      






      data ixran/9862/
      data ITYPE/1/,UPLO/'U'/
      
      ixran=ixran*int(u)
      
      open(unit=80,file="eng_list_u9_m.txt",status='unknown')
      
      nn2=nn*2
      
      iRGdim=0
      do i2=1,10
      iRGdim(i2)=50*i2
      end do
c      iRGdim=0
c      iRGdim(1)=10
      
      
      illup=0
      illdown=0
      
      
      pai=dacos(-1.0d0)
      engold=100.0d0
      
      hu=0.0d0
      hu2=0.0d0
      alpha1=0.0d0
      gin2=0.0d0
      g2=0.0d0
      
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


       do 
       read(5,*) iLdim
       if(iLdim.eq.istdim) then
          goto 8989
       end if
       end do
8989   write(6,*) 'start_dim',iLdim
c read H.F.A.result 
      do l=1,iLdim
         do j=1,m2
            do i=1,nn2
               read(5,*) fai(i,j,l)
            end do
         end do
      end do
      

c trun local quantization axis
      do l=1,iLdim
         do j=1,m2
            do i=1,nn2
               fai2(i,j)=fai(i,j,l)
               fai(i,j,l)=0.0d0
            end do
         end do

         do j=1,m2
            do i=1,nn
               fai(i,j,l)=fai(i,j,l)
     &        +dcos(-(lsub(i)-1)*2.0d0*pai/3.0d0)*fai2(i,j)
               fai(i,j,l)=fai(i,j,l)
     &        -dsin(-(lsub(i)-1)*2.0d0*pai/3.0d0)*fai2(i+nn,j)
               fai(i+nn,j,l)=fai(i+nn,j,l)
     &        +dsin(-(lsub(i)-1)*2.0d0*pai/3.0d0)*fai2(i,j)
               fai(i+nn,j,l)=fai(i+nn,j,l)
     &        +dcos(-(lsub(i)-1)*2.0d0*pai/3.0d0)*fai2(i+nn,j)
            end do
         end do
      end do
c end trun local quantization axis


      if(iLdim.eq.1) then
      l=1
         do inum=1,m2
            do i=1,nn2
               fai0(i,inum)=fai(i,inum,l)
            end do
         end do
      end if
c end read First W.F.








c Expand dimension ---------------------------------------------------
      do 600 iRGstep=1,Lmax
      if(iRGdim(iRGstep).eq.0) goto 600
      
      ierror=0
361   if(iRGtype.eq.1.and.iRGstep.ge.2) then
         iexpdim=iRGdim(iRGstep-1)+1
      else
         iexpdim=1
      end if

      delta=delexp
      call sht(u,delta,alpha1)
      call kinetictri(pai,nn,t1,t2,delta/2.0d0,hkL0,hku
     &,ihop1,ihop2,hnqei,ee)
      
      r=ran2(-ixran)
      write(6,*) 'Expansion delta=',delta,ixran
      write(80,*) 'Expansion delta=',delta,ixran

      
      do 360 iLdim=iexpdim,iRGdim(iRGstep)
c select SH field ------------------------------------------------
      do ijm=1,nn
      if(ran2(ixran).lt.0.5d0) then
         is=1
      else
         is=-1
      end if
      hu(ijm)=dexp(alpha1(is))
      hu(ijm+nn)=dexp(alpha1(-is))
      end do
      
      do j=1,nn2
      do i=1,nn2
         hu2(i,j)=0.0d0
         do k=1,nn2
         hu2(i,j)=hu2(i,j)+hkL0(i,k)*hu(k)*hkL0(k,j)
      end do
      end do
      end do
c end select SH field ------------------------------------------------

      
      do inum=1,m2
         do i=1,nn2
            fai(i,inum,iLdim)=0.0d0
               do j=1,nn2
                  fai(i,inum,iLdim)=fai(i,inum,iLdim)
     &           +hu2(i,j)*fai0(j,inum)
               end do
         end do
      end do
      
      do j=1,m2
         do i=1,nn2
            stfai(i,j,iLdim)=fai(i,j,iLdim)
         end do
      end do

      call schmidttri(nn2,m2,Lmax,fai,fai2,iLdim,a,aeigen,cc)
      
      LWORK=m2
      do ll=1,iLdim
         call mkgreentri(ll,iLdim,nn2,m2,Lmax,IPIV,WORK,LWORK
     &  ,buf,a,fai,g,gin)
         call mkHtri(pai,ll,iLdim,nn,u,t1,t2,Lmax
     &  ,ihop1,ihop2,h,hin,g,gin)
      end do      
      
      label=0
      LWORK2=3*iLdim-1
      call geed0(Lmax,iLdim,h,hin
     &,hgeed,hingeed,deigen,eigen,WORK2,LWORK2,'N',label)
      pirgeng=eigen
      if(label.eq.1) then
         goto 361
      end if
      
      if(iRGtype.eq.1) then
         write(6,*) 'Local',iLdim,pirgeng
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
      call kinetictri(pai,nn,t1,t2,delta/2.0d0,hkL0,hku
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
      do j=1,m2
         do i=1,nn2
            faiold(i,j)=fai(i,j,lite)
         end do
      end do
      do j=1,m2
         do i=1,nn2
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
      do j=1,m2
         do i=1,nn2
            fai(i,j,lite)=0.0d0
               do k=1,nn2
                  fai(i,j,lite)=fai(i,j,lite)
     &           +hkL0(i,k)*stfai(k,j,lite)
               end do
         end do
      end do

      call schmidttri(nn2,m2,Lmax,fai,fai2,lite,a,aeigen,cc)
      do ll=1,iLdim
      call mkgreentri(ll,lite,nn2,m2,Lmax,IPIV,WORK,LWORK
     &            ,buf,a,fai,g,gin)
      call mkHtri(pai,ll,lite,nn,u,t1,t2,Lmax
     &        ,ihop1,ihop2,h,hin,g,gin)
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
      call Urgtri(pai,nn,m2,u,t1,t2,Lmax,iLdim,lite,ijm,alpha1
     &,h,hin,ihop1,ihop2,hsub,hinsub,hgeed,hingeed,WORK2,LWORK2
     &,iss,g,gin,gsub,pirgeng,lurg,g2,gin2,pirgmp,hu,deigen
     &,ipstepu,icoarse)
      end do

      do i=1,nn2
         do j=1,nn2
            hu2(i,j)=0.0d0
            do k=1,nn2
               hu2(i,j)=hu2(i,j)
     &        +hkL0(i,k)*hu(k)*hkL0(k,j)
            end do
         end do
      end do
c preserve updated data ---------------------------------------
      do j=1,m2
         do i=1,nn2
            fai(i,j,lite)=0.0d0
               do k=1,nn2
                  fai(i,j,lite)=fai(i,j,lite)
     &           +hu2(i,k)*stfai(k,j,lite)
               end do
         end do
      end do

      do j=1,m2
         do i=1,nn2
            stfai(i,j,lite)=fai(i,j,lite)
         end do
      end do

      call schmidttri(nn2,m2,Lmax,fai,fai2,lite,a,aeigen,cc)
      do ll=1,iLdim
      call mkgreentri(ll,lite,nn2,m2,Lmax,IPIV,WORK,LWORK
     &            ,buf,a,fai,g,gin)
      call mkHtri(pai,ll,lite,nn,u,t1,t2,Lmax
     &        ,ihop1,ihop2,h,hin,g,gin)
      end do      
      do ll=lite,iLdim
         h(lite,ll)=h(ll,lite)
         hin(lite,ll)=hin(ll,lite)
      end do
c end preserve updated data ----------------------------------
      label=0
      call geed0(Lmax,iLdim,h,hin
     &,hgeed,hingeed,deigen,eigen,WORK2,LWORK2,'N',label)
      if(label.eq.1) then
         ierror=ierror+1
         write(6,*) 'back',ierror
      isign=0
      do lr=1,iLdim
         do ll=1,lr
            if(hin(ll,lr).le.0.0d0) then
               isign=isign+1
            end if
         end do
      end do
      tsign=dble(iLdim*(iLdim-1))/2.0d0
      write(6,*) 'minus sign',dble(isign)/tsign
      write(80,*) 'minus sign',dble(isign)/tsign
         if(ierror.eq.10) stop
         goto 361
      end if
      
      pirgeng=eigen

      if((engold-pirgeng).lt.0) then
      pirgeng=engold
      do j=1,m2
         do i=1,nn2
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

      label=0
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
      do j=1,m2
         do i=1,nn2
            faiold(i,j)=fai(i,j,lite)
         end do
      end do
      do j=1,m2
         do i=1,nn2
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
      do j=1,m2
         do i=1,nn2
            fai(i,j,lite)=0.0d0
               do k=1,nn2
                  fai(i,j,lite)=fai(i,j,lite)
     &           +hkL0(i,k)*stfai(k,j,lite)
               end do
         end do
      end do

      do j=1,m2
         do i=1,nn2
            stfai(i,j,lite)=fai(i,j,lite)
         end do
      end do


      call schmidttri(nn2,m2,Lmax,fai,fai2,lite,a,aeigen,cc)
      do ll=1,iLdim
      call mkgreentri(ll,lite,nn2,m2,Lmax,IPIV,WORK,LWORK
     &            ,buf,a,fai,g,gin)
      call mkHtri(pai,ll,lite,nn,u,t1,t2,Lmax
     &        ,ihop1,ihop2,h,hin,g,gin)
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
      do j=1,m2
         do i=1,nn2
         fai(i,j,lite)=faiold(i,j)
         end do
      end do   
      do j=1,m2
         do i=1,nn2
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
      call mkgreentri(ll,lite,nn2,m2,Lmax,IPIV,WORK,LWORK
     &            ,buf,a,fai,g,gin)
      end do
      
      engold=pirgeng
      
c Hu term -----------------------------------------------
      do  ijmk=1,nn
      ijm=isite(ijmk)
      call Urgtri(pai,nn,m2,u,t1,t2,Lmax,iLdim,lite,ijm,alpha1
     &,h,hin,ihop1,ihop2,hsub,hinsub,hgeed,hingeed,WORK2,LWORK2
     &,iss,g,gin,gsub,pirgeng,lurg,g2,gin2,pirgmp,hu,deigen
     &,ipstepu,icoarse)
      end do
         
c preserve updated data ---------------------------------------
      do j=1,m2
         do i=1,nn2
            fai(i,j,lite)=0.0d0
                  fai(i,j,lite)=fai(i,j,lite)
     &           +hu(i)*stfai(i,j,lite)
         end do
      end do

      do j=1,m2
         do i=1,nn2
            stfai(i,j,lite)=fai(i,j,lite)
         end do
      end do

      call schmidttri(nn2,m2,Lmax,fai,fai2,lite,a,aeigen,cc)
      do ll=1,iLdim
      call mkgreentri(ll,lite,nn2,m2,Lmax,IPIV,WORK,LWORK
     &            ,buf,a,fai,g,gin)
      call mkHtri(pai,ll,lite,nn,u,t1,t2,Lmax
     &        ,ihop1,ihop2,h,hin,g,gin)
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
      do j=1,m2
         do i=1,nn2
            faiold(i,j)=fai(i,j,lite)
         end do
      end do
      do j=1,m2
         do i=1,nn2
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
      do j=1,m2
         do i=1,nn2
            fai(i,j,lite)=0.0d0
               do k=1,nn2
                  fai(i,j,lite)=fai(i,j,lite)
     &           +hkL0(i,k)*stfai(k,j,lite)
               end do
         end do
      end do

      do j=1,m2
         do i=1,nn2
            stfai(i,j,lite)=fai(i,j,lite)
         end do
      end do

      call schmidttri(nn2,m2,Lmax,fai,fai2,lite,a,aeigen,cc)
      do ll=1,iLdim
      call mkgreentri(ll,lite,nn2,m2,Lmax,IPIV,WORK,LWORK
     &            ,buf,a,fai,g,gin)
      call mkHtri(pai,ll,lite,nn,u,t1,t2,Lmax
     &        ,ihop1,ihop2,h,hin,g,gin)
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
      do j=1,m2
         do i=1,nn2
            fai(i,j,lite)=faiold(i,j)
         end do
      end do   
      do j=1,m2
         do i=1,nn2
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
      if(label.eq.1) then
         ierror=ierror+1
         write(6,*) 'back',ierror
      isign=0
      do lr=1,iLdim
         do ll=1,lr
            if(hin(ll,lr).le.0.0d0) then
               isign=isign+1
            end if
         end do
      end do
      tsign=dble(iLdim*(iLdim-1))/2.0d0
      write(6,*) 'minus sign',dble(isign)/tsign
      write(80,*) 'minus sign',dble(isign)/tsign
         if(ierror.eq.10) stop
         goto 361
      end if
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
      isign=0
      do lr=1,iLdim
         do ll=1,lr
            if(hin(ll,lr).le.0.0d0) then
               isign=isign+1
            end if
         end do
      end do
      tsign=dble(iLdim*(iLdim-1))/2.0d0
      write(6,*) 'minus sign',dble(isign)/tsign
      write(80,*) 'minus sign',dble(isign)/tsign
      write(80,*)
c      write(6,*) 
c      do lr=1,iLdim
c         do ll=1,lr
c            write(6,*) hin(ll,lr)
c         end do
c      end do
      end if
      

c retrun quantization axis      
      do l=1,iLdim
         do j=1,m2
            do i=1,nn2
            fai2(i,j)=fai(i,j,l)
            fai(i,j,l)=0.0d0
            end do
         end do
         
         do j=1,m2
            do i=1,nn
               fai(i,j,l)=fai(i,j,l)
     &        +dcos((lsub(i)-1)*2.0d0*pai/3.0d0)*fai2(i,j)
               fai(i,j,l)=fai(i,j,l)
     &        -dsin((lsub(i)-1)*2.0d0*pai/3.0d0)*fai2(i+nn,j)
               fai(i+nn,j,l)=fai(i+nn,j,l)
     &        +dsin((lsub(i)-1)*2.0d0*pai/3.0d0)*fai2(i,j)
               fai(i+nn,j,l)=fai(i+nn,j,l)
     &        +dcos((lsub(i)-1)*2.0d0*pai/3.0d0)*fai2(i+nn,j)
c               write(6,*) 'sublattice',i,lsub(i)-1
            end do
         end do
      end do
c end retrun quantization axis      



      if(iWF.eq.1) then
      open(unit=86,file="WF_u9_m.txt",status='unknown')
      write(86,*) iLdim
      do l=1,iLdim
         do inum=1,m2
            do i=1,nn2
               write(86,*) stfai(i,inum,l)
            end do
         end do
      end do
      end if
      
      if(iengv.eq.1) then
      hold=h
      hinold=hin
      call engv(u,t1,t2,iLdim,ispin,iqpon
     &,fai,ihop1,ihop2,lsitex,lsitey,lsub,WORK2,LWORK2)
      h=hold
      hin=hinold
      end if

      
      
      
      write(6,*) 'stop dim.=',iLdim
      write(6,*) 

      goto 600
      end if
      

36    end do      
600   end do      
      
      stop
      end
c end PIRG checker board
c----------------------------------------------------------------
      
      
      
      
      
      
      
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine Urgtri(pai,nn,m2,u,t1,t2,Lmax,iLdim,lite,ijm,alpha1
     &,h,hin,ihop1,ihop2,hsub,hinsub,hgeed,hingeed,WORK2,LWORK2
     &,iss,g,gin,gsub,pirgeng,lurg,g2,gin2,pirgmp,hu,deigen
     &,ipstepu,icoarse)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 UPLO
      
      dimension alpha1(-1:1)
      dimension hsub(Lmax,2),hinsub(Lmax,2)
      dimension h(Lmax,Lmax),hin(Lmax,Lmax)
      dimension hgeed(Lmax,Lmax),hingeed(Lmax,Lmax)
      dimension WORK2(LWORK2),deigen(Lmax)
      dimension g(2*nn,2*nn,Lmax),gin(Lmax)
      dimension gsub(2*nn,2*nn)
      dimension g2(2*nn,2*nn,Lmax,2),gin2(Lmax,2)
      dimension pirgmp(2),hu(2*nn)
      
      dimension ihop1(nn,4),ihop2(nn,2)
      dimension lurg(Lmax,Lmax),iss(2)
      
      nn2=2*nn
      
      do 1140 is=1,2
      is2=iss(is)
      do 400 llsub=1,iLdim-1
      ll=lurg(llsub,lite)
c up spin
      ginsub=(1.0d0
     &+(dexp(alpha1(is2))-1.0d0)*g(ijm,ijm,ll))*gin(ll)
      do j=1,nn2
         do i=1,nn2
      gsub(i,j)=g(i,j,ll)
     &-g(i,ijm,ll)*(dexp(alpha1(is2))-1.0d0)*g(ijm,j,ll)
     &/(1.0d0+g(ijm,ijm,ll)*(dexp(alpha1(is2))-1.0d0))
         end do
      end do
      do j=1,nn2
      gsub(ijm,j)=gsub(ijm,j)*dexp(alpha1(is2))
      end do
c down spin
      gin2(ll,is)=(1.0d0
     &+(dexp(alpha1(-is2))-1.0d0)*gsub(ijm+nn,ijm+nn))*ginsub
      do j=1,nn2
         do i=1,nn2
      g2(i,j,ll,is)=gsub(i,j)
     &-gsub(i,ijm+nn)*(dexp(alpha1(-is2))-1.0d0)*gsub(ijm+nn,j)
     &/(1.0d0+gsub(ijm+nn,ijm+nn)*(dexp(alpha1(-is2))-1.0d0))
         end do
      end do
      do j=1,nn2
      g2(ijm+nn,j,ll,is)=g2(ijm+nn,j,ll,is)*dexp(alpha1(-is2))
      end do
400   continue

      ll=lite
c first
c up spin
      ginsub=(1.0d0
     &+(dexp(alpha1(is2))-1.0d0)*g(ijm,ijm,ll))*gin(ll)
      do j=1,nn2
         do i=1,nn2
      gsub(i,j)=g(i,j,ll)
     &-g(i,ijm,ll)*(dexp(alpha1(is2))-1.0d0)*g(ijm,j,ll)
     &/(1.0d0+g(ijm,ijm,ll)*(dexp(alpha1(is2))-1.0d0))
         end do
      end do
      do j=1,nn2
      gsub(ijm,j)=gsub(ijm,j)*dexp(alpha1(is2))
      end do
c down spin
      gin2(ll,is)=(1.0d0
     &+(dexp(alpha1(-is2))-1.0d0)*gsub(ijm+nn,ijm+nn))*ginsub
      do j=1,nn2
         do i=1,nn2
      g2(i,j,ll,is)=gsub(i,j)
     &-gsub(i,ijm+nn)*(dexp(alpha1(-is2))-1.0d0)*gsub(ijm+nn,j)
     &/(1.0d0+gsub(ijm+nn,ijm+nn)*(dexp(alpha1(-is2))-1.0d0))
         end do
      end do
      do j=1,nn2
      g2(ijm+nn,j,ll,is)=g2(ijm+nn,j,ll,is)*dexp(alpha1(-is2))
      end do

c second
c up spin
      ginsub=(1.0d0
     &+(dexp(alpha1(is2))-1.0d0)*g2(ijm,ijm,ll,is))*gin2(ll,is)
      do j=1,nn2
         do i=1,nn2
      gsub(i,j)=g2(i,j,ll,is)
     &-g2(i,ijm,ll,is)*(dexp(alpha1(is2))-1.0d0)*g2(ijm,j,ll,is)
     &/(1.0d0+g2(ijm,ijm,ll,is)*(dexp(alpha1(is2))-1.0d0))
         end do
      end do
      do i=1,nn2
      gsub(i,ijm)=gsub(i,ijm)*dexp(alpha1(is2))
      end do
c down spin
      gin2(ll,is)=(1.0d0
     &+(dexp(alpha1(-is2))-1.0d0)*gsub(ijm+nn,ijm+nn))*ginsub
      do j=1,nn2
         do i=1,nn2
      g2(i,j,ll,is)=gsub(i,j)
     &-gsub(i,ijm+nn)*(dexp(alpha1(-is2))-1.0d0)*gsub(ijm+nn,j)
     &/(1.0d0+gsub(ijm+nn,ijm+nn)*(dexp(alpha1(-is2))-1.0d0))
         end do
      end do
      do i=1,nn2
      g2(i,ijm+nn,ll,is)=g2(i,ijm+nn,ll,is)*dexp(alpha1(-is2))
      end do


c make Hamiltonial --------------------------------------------------
      do ll=1,iLdim
      hinsub(ll,is)=gin2(ll,is)
      hsub(ll,is)=0.0d0
      
      do i=1,nn
      
         do j=1,2
      hsub(ll,is)=hsub(ll,is)-t1*(
     &+dcos(2.0d0*pai/3.0d0)*g2(i,ihop1(i,j),ll,is)
     &-dsin(2.0d0*pai/3.0d0)*g2(i,ihop1(i,j)+nn,ll,is)
     &+dsin(2.0d0*pai/3.0d0)*g2(i+nn,ihop1(i,j),ll,is)
     &+dcos(2.0d0*pai/3.0d0)*g2(i+nn,ihop1(i,j)+nn,ll,is))
         end do
         
         do j=3,4
      hsub(ll,is)=hsub(ll,is)-t1*(
     &+dcos(-2.0d0*pai/3.0d0)*g2(i,ihop1(i,j),ll,is)
     &-dsin(-2.0d0*pai/3.0d0)*g2(i,ihop1(i,j)+nn,ll,is)
     &+dsin(-2.0d0*pai/3.0d0)*g2(i+nn,ihop1(i,j),ll,is)
     &+dcos(-2.0d0*pai/3.0d0)*g2(i+nn,ihop1(i,j)+nn,ll,is))
         end do

         j=1
      hsub(ll,is)=hsub(ll,is)-t2*(
     &+dcos(4.0d0*pai/3.0d0)*g2(i,ihop2(i,j),ll,is)
     &-dsin(4.0d0*pai/3.0d0)*g2(i,ihop2(i,j)+nn,ll,is)
     &+dsin(4.0d0*pai/3.0d0)*g2(i+nn,ihop2(i,j),ll,is)
     &+dcos(4.0d0*pai/3.0d0)*g2(i+nn,ihop2(i,j)+nn,ll,is))
         
         j=2
      hsub(ll,is)=hsub(ll,is)-t2*(
     &+dcos(-4.0d0*pai/3.0d0)*g2(i,ihop2(i,j),ll,is)
     &-dsin(-4.0d0*pai/3.0d0)*g2(i,ihop2(i,j)+nn,ll,is)
     &+dsin(-4.0d0*pai/3.0d0)*g2(i+nn,ihop2(i,j),ll,is)
     &+dcos(-4.0d0*pai/3.0d0)*g2(i+nn,ihop2(i,j)+nn,ll,is))


      hsub(ll,is)=hsub(ll,is)
     &+u*g2(i,i,ll,is)*g2(i+nn,i+nn,ll,is)
     &-u*g2(i+nn,i,ll,is)*g2(i,i+nn,ll,is)
      
      end do
      
      hsub(ll,is)=hsub(ll,is)*hinsub(ll,is)
      end do
c end make Hamiltonial ------------------------------------------------



      do ll=1,lite
         h(ll,lite)=hsub(ll,is)
         hin(ll,lite)=hinsub(ll,is)
      end do
      
      do ll=lite,iLdim
         h(lite,ll)=hsub(ll,is)
         hin(lite,ll)=hinsub(ll,is)
      end do
      call geed0(Lmax,iLdim,h,hin
     &,hgeed,hingeed,deigen,eigen,WORK2,LWORK2,'N',label)
     
      pirgmp(is)=eigen
1140  continue
      
c no change
      if(pirgmp(1).ge.pirgeng.and.pirgmp(2).ge.pirgeng
     &.and.icoarse.eq.0) then
      hu(ijm)=1.0d0
      hu(ijm+nn)=1.0d0
      else
      
      
      ipstepu=ipstepu+1
c minus SH field stable
      if((pirgmp(2).lt.pirgmp(1))) then
      hu(ijm)=dexp(alpha1(-1))
      hu(ijm+nn)=dexp(alpha1(1))
      
      do ll=1,lite
         h(ll,lite)=hsub(ll,2)
         hin(ll,lite)=hinsub(ll,2)
      end do
      do ll=lite,iLdim
         h(lite,ll)=hsub(ll,2)
         hin(lite,ll)=hinsub(ll,2)
      end do
      
      do l=1,iLdim
         gin(l)=gin2(l,2)
         do j=1,nn2
            do i=1,nn2
               g(i,j,l)=g2(i,j,l,2)
            end do
         end do
      end do
      pirgeng=pirgmp(2)

c plus SH field stable
      else
      hu(ijm)=dexp(alpha1(1))
      hu(ijm+nn)=dexp(alpha1(-1))

      do ll=1,lite
         h(ll,lite)=hsub(ll,1)
         hin(ll,lite)=hinsub(ll,1)
      end do
      do ll=lite,iLdim
         h(lite,ll)=hsub(ll,1)
         hin(lite,ll)=hinsub(ll,1)
      end do

      do l=1,iLdim
         gin(l)=gin2(l,1)
         do j=1,nn2
            do i=1,nn2
               g(i,j,l)=g2(i,j,l,1)
            end do
         end do
      end do
      pirgeng=pirgmp(1)
      end if
      
      end if
      
      return
      end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc








ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mkgreentri(ll,lr,nn2,m2,Lmax,IPIV,WORK,LWORK
     &,buf,a,fai,g,gin)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n) 
      character*1 UPLO 
      
      dimension g(nn2,nn2,Lmax),gin(Lmax)
      dimension a(m2,m2),buf(m2,nn2)
      dimension fai(nn2,m2,Lmax)
      dimension IPIV(m2),WORK(LWORK)
      
      do inum=1,m2
      do jnum=1,m2
      a(inum,jnum)=0.0d0
         do i=1,nn2
            a(inum,jnum)=a(inum,jnum)
     &     +fai(i,inum,ll)*fai(i,jnum,lr)
         end do
      end do
      end do
c end make a matrix
      call DGETRF(m2,m2,a,m2,IPIV,INFO)
      if(INFO.ne.0) then
         write(6,*) 'LU_failure  INFO=',INFO
         stop
      end if
      gin(ll)=1.0d0
      do j=1,m2
         gin(ll)=gin(ll)*PIV(IPIV(j),j)*a(j,j)
      end do
      
      call DGETRI(m2,a,m2,IPIV,WORK,LWORK,INFO)
      if(INFO.ne.0) then
         write(6,*) 'inverse_failure up spin INFO=',INFO
         stop
      end if

c make g matrix
c fast code by Koga san
      do j=1,nn2
      do i=1,m2
         buf(i,j)=0.0d0
         do l=1,m2
            buf(i,j)=buf(i,j)+a(i,l)*fai(j,l,ll)
         end do
      end do
      end do
   
      do j=1,nn2
      do i=1,nn2
         g(i,j,ll)=0.0d0
         do k=1,m2
            g(i,j,ll)=g(i,j,ll)+fai(i,k,lr)*buf(k,j)
         end do
      end do
      end do

c late code
c      do j=1,nn2
c      do i=1,nn2
c         g(i,j,ll)=0.0d0
c         do k=1,m2
c         do l=1,m2
c            g(i,j,ll)=g(i,j,ll)
c     &     +fai(i,k,lr)*a(k,l)*fai(j,l,ll)
c         end do
c         end do
c      end do
c      end do


      return
      end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc








ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine schmidttri(nn2,m2,Lmax,fai,fai2,lite,a
     &,aeigen,cc)

      implicit real*8(a-h,o-z)
      implicit integer*4(i-n) 

      dimension fai(nn2,m2,Lmax),fai2(nn2,m2)
      dimension a(m2,m2),aeigen(m2),cc(m2)

c inner product      
      do jnum=1,m2
      do inum=1,m2
         a(inum,jnum)=0.0d0
         do i=1,nn2
            a(inum,jnum)=a(inum,jnum)
     &     +fai(i,inum,lite)*fai(i,jnum,lite)
         end do
      end do
      end do
      

      call tred2(a,m2,m2,aeigen,cc)
      call tqli(aeigen,cc,m2,m2,a)

      
      do j=1,m2
         do i=1,m2
            a(i,j)=a(i,j)/dsqrt(aeigen(j))
         end do
      end do
      
       
      do j=1,m2
         do i=1,nn2
            fai2(i,j)=fai(i,j,lite)
         end do
      end do


      
      do j=1,m2
         do i=1,nn2
            fai(i,j,lite)=0.0d0
            do k=1,m2
               fai(i,j,lite)=fai(i,j,lite)
     &                       +fai2(i,k)*a(k,j)
            end do
         end do
      end do



      return
      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc










ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c kinetic energy matrix
      subroutine kinetictri(pai,nn,t1,t2,delta,hk,hku
     &,ihop1,ihop2,hnqei,ee)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      
      dimension hk(nn*2,nn*2)
      dimension hku(nn*2,nn*2)
      dimension hnqei(nn*2),ee(nn*2)
      
      dimension ihop1(nn,4),ihop2(nn,2)
      
      nn2=nn*2

c Matrix of kinetic energy
         hk=0.0d0
         hku=0.0d0

c Matrix of kinetic energy ccccccccccccccccccccccccccccccc
c      do 7 i=1,nn
c      do k=1,4
c      hk(i,ihop1(i,k))=
c     &hk(i,ihop1(i,k))+t1*delta
c      end do
c      do k=1,2
c      hk(i,ihop2(i,k))=
c     &hk(i,ihop2(i,k))+t2*delta
c      end do
c
c      do k=1,4
c      hk(i+nn,ihop1(i,k)+nn)=
c     &hk(i+nn,ihop1(i,k)+nn)+t1*delta
c      end do
c      do k=1,2
c      hk(i+nn,ihop2(i,k)+nn)=
c     &hk(i+nn,ihop2(i,k)+nn)+t2*delta
c      end do
c7     continue 
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c Matrix of kinetic energy ccccccccccccccccccccccccccccccc
      do i=1,nn
         do j=1,2
      hk(i,ihop1(i,j))=hk(i,ihop1(i,j))
     &+t1*dcos(2.0d0*pai/3.0d0)*delta
      hk(i,ihop1(i,j)+nn)=hk(i,ihop1(i,j)+nn)
     &-t1*dsin(2.0d0*pai/3.0d0)*delta
      hk(i+nn,ihop1(i,j))=hk(i+nn,ihop1(i,j))
     &+t1*dsin(2.0d0*pai/3.0d0)*delta
      hk(i+nn,ihop1(i,j)+nn)=hk(i+nn,ihop1(i,j)+nn)
     &+t1*dcos(2.0d0*pai/3.0d0)*delta
         end do
         
         do j=3,4
      hk(i,ihop1(i,j))=hk(i,ihop1(i,j))
     &+t1*dcos(-2.0d0*pai/3.0d0)*delta
      hk(i,ihop1(i,j)+nn)=hk(i,ihop1(i,j)+nn)
     &-t1*dsin(-2.0d0*pai/3.0d0)*delta
      hk(i+nn,ihop1(i,j))=hk(i+nn,ihop1(i,j))
     &+t1*dsin(-2.0d0*pai/3.0d0)*delta
      hk(i+nn,ihop1(i,j)+nn)=hk(i+nn,ihop1(i,j)+nn)
     &+t1*dcos(-2.0d0*pai/3.0d0)*delta
         end do


         j=1
      hk(i,ihop2(i,j))=hk(i,ihop2(i,j))
     &+t2*dcos(4.0d0*pai/3.0d0)*delta
      hk(i,ihop2(i,j)+nn)=hk(i,ihop2(i,j)+nn)
     &-t2*dsin(4.0d0*pai/3.0d0)*delta
      hk(i+nn,ihop2(i,j))=hk(i+nn,ihop2(i,j))
     &+t2*dsin(4.0d0*pai/3.0d0)*delta
      hk(i+nn,ihop2(i,j)+nn)=hk(i+nn,ihop2(i,j)+nn)
     &+t2*dcos(4.0d0*pai/3.0d0)*delta
         
         j=2
      hk(i,ihop2(i,j))=hk(i,ihop2(i,j))
     &+t2*dcos(-4.0d0*pai/3.0d0)*delta
      hk(i,ihop2(i,j)+nn)=hk(i,ihop2(i,j)+nn)
     &-t2*dsin(-4.0d0*pai/3.0d0)*delta
      hk(i+nn,ihop2(i,j))=hk(i+nn,ihop2(i,j))
     &+t2*dsin(-4.0d0*pai/3.0d0)*delta
      hk(i+nn,ihop2(i,j)+nn)=hk(i+nn,ihop2(i,j)+nn)
     &+t2*dcos(-4.0d0*pai/3.0d0)*delta

      end do
ccccccccccccccccccccccccccccccccccccccccccccccccccc


      call tred2(hk,nn2,nn2,hnqei,ee)
      call tqli(hnqei,ee,nn2,nn2,hk)
      

      do i=1,nn2
         do j=1,nn2
            hku(j,i)=hk(i,j)
         end do
      end do

      hk=0.0d0
      do l=1,nn2
      do i=1,nn2
      do j=1,nn2
         hk(i,l)=hk(i,l)+hku(j,i)*dexp(hnqei(j))*hku(j,l)
      end do
      end do
      end do
      

c      write(6,*) 'reset dim.itehk=',ite,delta

      return
      end
c end kinetc energy matrix
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mkHtri(pai,ll,lr,nn,u,t1,t2,Lmax
     &,ihop1,ihop2,h,hin,g,gin)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n) 
      character*1 UPLO 
      
      dimension h(Lmax,Lmax),hin(Lmax,Lmax)
      dimension g(2*nn,2*nn,Lmax),gin(Lmax)
      
      dimension ihop1(nn,4),ihop2(nn,2)
      nn2=2*nn
      

      
      hin(ll,lr)=gin(ll)
      h(ll,lr)=0.0d0
      do i=1,nn
      
         do j=1,2
      h(ll,lr)=h(ll,lr)-t1*(
     &+dcos(2.0d0*pai/3.0d0)*g(i,ihop1(i,j),ll)
     &-dsin(2.0d0*pai/3.0d0)*g(i,ihop1(i,j)+nn,ll)
     &+dsin(2.0d0*pai/3.0d0)*g(i+nn,ihop1(i,j),ll)
     &+dcos(2.0d0*pai/3.0d0)*g(i+nn,ihop1(i,j)+nn,ll))
         end do
         
         do j=3,4
      h(ll,lr)=h(ll,lr)-t1*(
     &+dcos(-2.0d0*pai/3.0d0)*g(i,ihop1(i,j),ll)
     &-dsin(-2.0d0*pai/3.0d0)*g(i,ihop1(i,j)+nn,ll)
     &+dsin(-2.0d0*pai/3.0d0)*g(i+nn,ihop1(i,j),ll)
     &+dcos(-2.0d0*pai/3.0d0)*g(i+nn,ihop1(i,j)+nn,ll))
         end do

         j=1
      h(ll,lr)=h(ll,lr)-t2*(
     &+dcos(4.0d0*pai/3.0d0)*g(i,ihop2(i,j),ll)
     &-dsin(4.0d0*pai/3.0d0)*g(i,ihop2(i,j)+nn,ll)
     &+dsin(4.0d0*pai/3.0d0)*g(i+nn,ihop2(i,j),ll)
     &+dcos(4.0d0*pai/3.0d0)*g(i+nn,ihop2(i,j)+nn,ll))
         
         j=2
      h(ll,lr)=h(ll,lr)-t2*(
     &+dcos(-4.0d0*pai/3.0d0)*g(i,ihop2(i,j),ll)
     &-dsin(-4.0d0*pai/3.0d0)*g(i,ihop2(i,j)+nn,ll)
     &+dsin(-4.0d0*pai/3.0d0)*g(i+nn,ihop2(i,j),ll)
     &+dcos(-4.0d0*pai/3.0d0)*g(i+nn,ihop2(i,j)+nn,ll))


      h(ll,lr)=h(ll,lr)
     &+u*g(i,i,ll)*g(i+nn,i+nn,ll)
     &-u*g(i+nn,i,ll)*g(i,i+nn,ll)
      end do
      h(ll,lr)=h(ll,lr)*hin(ll,lr)

      return
      end 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   
      
