* ----------------
* site    states
* 4        4900
* 5       63504
* 6      853776
* 7    11778624  57 (a) and (b)
* 8   165636900  65 (a) and (b)
* 9  2363904400
* ----------------
      implicit none
      integer num,MAX_ND
      parameter(num=8,MAX_ND=165636900)
      
      integer n_up,n_dw
      integer nd,need(2,165636900),lv(2,0:4**8)

      integer h_pt(-1:65,165636900)
      double precision U

      double precision q1(165636900),q2(165636900),p(165636900)
      double precision w1(165636900),w2(165636900),w3(165636900)
      double precision a(500),b(500),ee,s,ret

      integer i,j

      n_up=num
      n_dw=num

    
      call ndnd(num,n_up,n_dw,nd,need,lv)      
      call init_H(num,nd,need,lv,h_pt)

c      open(555,status='unknown',file='double_u6_m.txt')
      open(556,status='unknown',file='energy_u6_m.txt')

c      do U=2.0d0,2.1d0,0.2d0
      U=6.0d0

         j=0
         do i=1,nd
            if(j<h_pt(0,i)) j=h_pt(0,i)
         end do
c         write(6,*) j
         if(j>65) then
            write(6,*) j," > 65"
            stop
         end if
 
         s=0.0d0
         do i=1,nd
            q2(i)=1.0d0/dble(i)
            s=s+q2(i)*q2(i)
         end do
         
         s=1.0d0/sqrt(s)
         do i=1,nd
            q2(i)=s*q2(i)
         end do
         
         
         call lanczos1(nd,h_pt,q1,q2,p,a,b,ee,i,U)
         write(556,*) U,ee/dble(num*2) 
         
c         call cg(nd,h_pt,U,q1,q2,p,w1,w2,w3,ee)

c         call dd(num,nd,need,q1,ret)
c         write(555,*) U,ret

c      end do
      stop
      end

      subroutine dd(num,nd,need,vec,s)
      implicit none
      integer num,nd,need(2,nd),ddd(0:3)
      double precision vec(nd),s

      integer i,j,k

      data ddd/0,0,0,1/

      s=0.0d0
      do i=1,nd
         k=0
         do j=0,num-1
            k=k+ddd(mod(need(1,i)/4**j,4))
            k=k+ddd(mod(need(2,i)/4**j,4))
         end do
         s=s+vec(i)*vec(i)*dble(k)
      end do

      s=s/dble(2*num)
      return
      end

      subroutine cg(nd,h_pt,U,vec,b,bb,r,p,ap,val)
      implicit none ! q1 : ground state
      integer nd,h_pt(-1:65,165636900)
      double precision b(165636900),bb(165636900),vec(165636900)
      double precision r(165636900),p(165636900),ap(165636900),val,s
      double precision alpha,beta,U
      integer i,j,jj,iter

      do i=1,nd
         r(i)=1.0d0/dble(i)
         p(i)=1.0d0/dble(i)
      end do

      jj=0
      iter=0

 1    call small_H(nd,h_pt,p,ap,U)

      do j=1,nd
         ap(j)=ap(j)-val*p(j)
      end do

      alpha=0.0d0
      s=0.0d0

      do j=1,nd
         alpha=alpha+p(j)*r(j)
         s=s+p(j)*ap(j)
      end do

      alpha=alpha/s

      beta=0.0d0
      do j=1,nd
         vec(j)=vec(j)+alpha*p(j)
         r(j)=r(j)-alpha*ap(j)
         beta=beta-r(j)*ap(j)
      end do
      beta=beta/s

      do j=1,nd
         p(j)=r(j)+beta*p(j)
      end do
      jj=jj+1
      iter=iter+1
      if(jj.le.9) goto 1

      s=0.0d0
      do i=1,nd
         s=s+vec(i)*vec(i)
      end do

      s=1.0d0/sqrt(s)
      do i=1,nd
         b(i)=s*vec(i)
      end do

      call small_H(nd,h_pt,b,bb,U)
      s=0.0d0


      do i=1,nd
         s=s+b(i)*bb(i)
      end do

      if(iter.gt.1000) then
         write(6,*) 'Error : CG method', abs(s-val)

         do i=1,nd
            vec(i)=b(i)
         end do
         return
c         stop
      end if
      if(abs(s-val).gt.1d-5) then
         jj=0
         goto 1
      end if
      write(6,*) iter,s-val
      if(abs(s-val).gt.1d-12) goto 1

c      write(6,*) 'success cg : ',iter,'iteration'
      do i=1,nd
         vec(i)=b(i)
      end do

      return
      end

      subroutine lanczos1(nd,h_pt,q1,q2,p,a,b,ee,iter,U)
      implicit none ! q2 : initial normalized vector
      integer nd,h_pt(-1:65,165636900)
      double precision ee
      double precision q1(165636900),q2(165636900),p(165636900)

      integer i,j,iter,count
      double precision s,ss,U
      double precision a(500),b(500),x(500),y(500)

      call small_H(nd,h_pt,q2,p,U)

      count=0
      ss=0.0d0

      i=0
 10   i=i+1
      count=count+1

      s=0.0d0
      do j=1,nd
         q1(j)=q2(j)
         s=s+q1(j)*p(j)
      end do
      a(i)=s

      s=0.0d0
      do j=1,nd
         q2(j)=p(j)-a(i)*q1(j)
         s=s+q2(j)*q2(j)
      end do

      b(i)=sqrt(s)

      do j=1,nd
         q2(j)=q2(j)/b(i)
      end do

      call small_H(nd,h_pt,q2,p,U)

      do j=1,nd
         p(j)=p(j)-b(i)*q1(j)
      end do
      if(i.ge.450 .and.mod(count,10).eq.0) then
c      if(mod(count,10).eq.0) then
c      if(mod(count,1).eq.0) then
         do j=1,i
            x(j)=a(j)
            y(j+1)=b(j)
         end do
         call tqli(x,y,i,i)
         call sort(i,x)

         write(6,*) i,x(1)
         if(abs(x(1)-ss).gt.1d-10) then
            ss=x(1)
            count=0
            goto 10
         end if
      else
         goto 10
      end if

      iter=i
      ee=x(1)

      return
      end


      subroutine small_H(nd,h_pt,vec1,vec2,U)
      implicit none
      integer i,j,k,nd
      integer h_pt(-1:65,165636900)
      double precision s,vec1(165636900),vec2(165636900),U

!$omp parallel do private(i,j,s)
      do i=1,nd
         s=U*h_pt(-1,i)*vec1(i)
         do j=2,h_pt(0,i)
            k=h_pt(j,i)
            if(k>0) then
               s=s+vec1(k)
            else
               s=s-vec1(-k)
            end if
c            s=s+h_val(j,i)*vec1(h_pt(j,i))
         end do
         vec2(i)=s
      end do

      return
      end

      subroutine get(num,m,n,need,ret)
      implicit none
      integer need(2,165636900),m,ret,num,n

      if(m<=num) then
         ret=need(1,n)/4**(num-m)
      else
         ret=need(2,n)/4**(2*num-m)
      end if
      ret=mod(ret,4)
      return 
      end

c            m1=2*num-datint(1,i)
c            m2=2*num-datint(2,i)

c            k1=mod(need(j)/4**m1,4)



      subroutine init_H(num,nd,need,lv,h_pt)
      implicit none
      integer num,nd,need(2,165636900),lv(2,0:4**8)
      integer h_pt(-1:65,165636900)

      integer mini(0:15,2),buf
      integer i,j,k,l,n,ja,jb,k1,k2,m1,m2,st,ed,i1,i2,infu,infd
      integer up(0:3),down(0:3)

      integer datint(2,48)
      integer n_int
      n_int = 48
      
      data datint/
     &1,2, 2,3, 3,4, 4,1,
     &5,6, 6,7, 7,8, 8,5,
     &9,10, 10,11, 11,12, 12,9,
     &13,14, 14,15, 15,16, 16,13,
     
     &1,5, 5,9, 9,13, 13,1,
     &2,6, 6,10, 10,14, 14,2,
     &3,7, 7,11, 11,15, 15,3,
     &4,8, 8,12, 12,16, 16,4,
     
     &1,6, 6,11, 11,16, 16,1,
     &2,7, 7,12, 12,13, 13,2,
     &3,8, 8,9, 9,14, 14,3,
     &4,5, 5,10, 10,15, 15,4/


      data up/0,1,0,1/
      data down/0,0,1,1/
      data mini/
     &0,4,0,6,1,0,3,0, 0,12,0,14,9,0,11,0,
     &0,0,8,9,0,0,12,13, 2,3,0,0,6,7,0,0/

!$omp parallel do private(i,s,j,k)
      do i=1,nd
         h_pt(-1,i)=0
         h_pt(0,i)=1
         h_pt(1,i)=i
c         s=0.0d0
c         do j=1,num*2
c            k=mod(need(i)/4**(2*num-j),4)
c            if( k.eq.3 ) then
c               s=s+U
c               h_pt(-1,i)=h_pt(-1,i)+1
c            end if
c         end do

         do j=1,num
            k=mod(need(1,i)/4**(num-j),4)
            if(k.eq.3) h_pt(-1,i)=h_pt(-1,i)+1
            k=mod(need(2,i)/4**(num-j),4)
            if(k.eq.3) h_pt(-1,i)=h_pt(-1,i)+1
         end do

      end do

!$omp parallel do private(i,j,k,l,m,n,k1,k2,ja,jb,info)


      do j=1,nd

         do i=1,n_int

            m1=datint(1,i)
            m2=datint(2,i)

c            call get(num,m1,j,need,k1)
c            call get(num,m2,j,need,k2)

      if(m1<=num) then
         k1=need(1,j)/4**(num-m1)
      else
         k1=need(2,j)/4**(2*num-m1)
      end if
      k1=mod(k1,4)

      if(m2<=num) then
         k2=need(1,j)/4**(num-m2)
      else
         k2=need(2,j)/4**(2*num-m2)
      end if
      k2=mod(k2,4)


            n=k1*4+k2

            st = datint(1,i)
            ed = datint(2,i)
            if(datint(2,i)<datint(1,i)) then
               st = datint(2,i)
               ed = datint(1,i)
            endif

c            info=0
            infu=0
            infd=0
            do l=st+1,ed-1
c               call get(num,l,j,need,buf)

               if(l<=num) then
                  buf=need(1,j)/4**(num-l)
               else
                  buf=need(2,j)/4**(2*num-l)
               end if
               buf=mod(buf,4)

               infu=infu+up(buf)
               infd=infd+down(buf)
            end do
* --- up spin --- *
            k=mini(n,1)
            if(k.ne.0) then

               i1=need(1,j)
               i2=need(2,j)

               if(m1<=num) then
                  i1=i1-k1*4**(num-m1)+(k/4)*4**(num-m1)
               else
                  i2=i2-k1*4**(2*num-m1)+(k/4)*4**(2*num-m1)
               end if

               if(m2<=num) then
                  i1=i1-k2*4**(num-m2)+mod(k,4)*4**(num-m2)
               else
                  i2=i2-k2*4**(2*num-m2)+mod(k,4)*4**(2*num-m2)
               end if

               ja=lv(1,i1)
               jb=lv(2,i2)
c               m=need(j)-k1*4**m1-k2*4**m2
c     &+(k/4)*4**m1+mod(k,4)*4**m2
c               ja=lv(1,m/4**num)
c               jb=lv(2,mod(m,4**num))

c               info=0
c               do l=st+1,ed-1
c                  info=info+up(mod(need(j)/4**(2*num-l),4))
c               end do

               infu=infu+1
               h_pt(0,j)=h_pt(0,j)+1
               h_pt(h_pt(0,j),j)=(ja+jb)*(-1)**infu
c               h_val(h_pt(0,j),j)=-1.0d0*(-1)**info
            end if

* --- down spin --- *
            k=mini(n,2)
            if(k.ne.0) then

               i1=need(1,j)
               i2=need(2,j)

               if(m1<=num) then
                  i1=i1-k1*4**(num-m1)+(k/4)*4**(num-m1)
               else
                  i2=i2-k1*4**(2*num-m1)+(k/4)*4**(2*num-m1)
               end if

               if(m2<=num) then
                  i1=i1-k2*4**(num-m2)+mod(k,4)*4**(num-m2)
               else
                  i2=i2-k2*4**(2*num-m2)+mod(k,4)*4**(2*num-m2)
               end if

               ja=lv(1,i1)
               jb=lv(2,i2)

c               m=need(j)-k1*4**m1-k2*4**m2
c     &+(k/4)*4**m1+mod(k,4)*4**m2
c               ja=lv(1,m/4**num)
c               jb=lv(2,mod(m,4**num))
               
c               infd=0
c               do l=st+1,ed-1
c                  info=info+down(mod(need(j)/4**(2*num-l),4))
c               end do

               infd =infd+1
               h_pt(0,j)=h_pt(0,j)+1
               h_pt(h_pt(0,j),j)=(ja+jb)*(-1)**infd
c               h_val(h_pt(0,j),j)=-1.0d0*(-1)**info
            end if
         end do
      end do

      return
      end

      subroutine ndnd(num,n_up,n_dw,nd,need,lv)
      implicit none
      integer i,j,ja,jb
      integer num,n_up,n_dw
      integer up(0:3),down(0:3)
      integer m,n,nd,need(2,165636900),lv(2,0:4**8)

*    0 empty 1 up 2 down 3 double

      data up/0,1,0,1/
      data down/0,0,1,1/

      nd=0

      do i=0,4**num-1
         m=0
         n=0
         do j=0,num-1
            m=m+up(mod(i/4**j,4))
            n=n+down(mod(i/4**j,4))
         end do
         lv(1,i)=m
         lv(2,i)=n
      end do
      do i=0,4**num-1
         do j=0,4**num-1

            if(lv(1,i)+lv(1,j).eq.n_up
     &           .and. lv(2,i)+lv(2,j).eq.n_dw) then
               nd=nd+1
               need(1,nd)=i
               need(2,nd)=j
c               need(nd)=i*4**num+j
            end if
         end do
      end do


      m=-1
      ja=0
      jb=0
      do j=1,nd
c         i=need(j)
c         if(i/4**num.eq.m) then
         if(need(1,j).eq.m) then
            ja=ja+1
         else
            jb=ja+jb
            ja=1
c            m=i/4**num
            m=need(1,j)
         end if
c         lv(1,i/4**num)=jb
c         lv(2,mod(i,4**num))=ja
         lv(1,need(1,j))=jb
         lv(2,need(2,j))=ja
      end do

* --- check --- *

      do j=1,nd
c         i=need(j)
c         if(j.ne.lv(1,i/4**num)+
c     &        lv(2,mod(i,4**num))) write(6,*) 'error'

         if(j.ne.lv(1,need(1,j))+lv(2,need(2,j))) write(6,*) 'error'
      end do

      return
      end

      
      SUBROUTINE sort(n,arr)
      INTEGER n,M,NSTACK
      DOUBLE PRECISION arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
      DOUBLE PRECISION a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,1,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=0
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(l))then
          temp=arr(l+1)
          arr(l+1)=arr(l)
          arr(l)=temp
        endif
        i=l+1
        j=ir
        a=arr(l)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END

      SUBROUTINE tqli(d,e,n,np)
      INTEGER n,np
      DOUBLE PRECISION  d(np),e(np)
CU    USES pythag
      INTEGER i,iter,l,m
      DOUBLE PRECISION  b,c,dd,f,g,p,r,s,pythag
      do 11 i=2,n
        e(i-1)=e(i)
11    continue
      e(n)=0.0d0
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
          g=(d(l+1)-d(l))/(2.0*e(l))
          r=pythag(g,1.0d0)
          g=d(m)-d(l)+e(l)/(g+sign(r,g))
          s=1.0d0
          c=1.0d0
          p=0.0d0
          do 14 i=m-1,l,-1
            f=s*e(i)
            b=c*e(i)
            r=pythag(f,g)
            e(i+1)=r
            if(r.eq.0.0d0)then
              d(i+1)=d(i+1)-p
              e(m)=0.0d0
              goto 1
            endif
            s=f/r
            c=g/r
            g=d(i+1)-p
            r=(d(i)-g)*s+2.0d0*c*b
            p=s*r
            d(i+1)=g+p
            g=c*r-b
14        continue

          d(l)=d(l)-p
          e(l)=g
          e(m)=0.0d0
          goto 1
        endif
15    continue
      return
      END

      FUNCTION pythag(a,b)
      DOUBLE PRECISION  a,b,pythag
      DOUBLE PRECISION  absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(1.0d0+(absb/absa)**2)
      else
        if(absb.eq.0.0d0)then
          pythag=0.0d0
        else
          pythag=absb*sqrt(1.0d0+(absa/absb)**2)
        endif
      endif
      return
      END
      
