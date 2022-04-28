ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c PIRG Checkerboard Lattice 24site
      subroutine engv(u,t1,t2,iLdim,ispin,iqpon
     &,faiqpL,ihop1,ihop2,lsitex,lsitey,lsub,WORK2,LWORK2)



      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 UPLO,JOBZ
      parameter(nn=24,m=nn/2,Lmax=500,lqp=12)
      
      
      
      dimension deigen(Lmax)
      dimension WORK2(LWORK2),hgeed(Lmax,Lmax),hingeed(Lmax,Lmax)
      
      dimension h(Lmax,Lmax),hin(Lmax,Lmax)
      dimension cron(2*nn,2*nn)
      




c size depend !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
c shift site numbers cccccccccccccccccccccccccccc
      dimension ihop1(nn,4),ihop2(nn,2)
      dimension lsitex(nn),lsitey(nn)
      dimension lsub(nn)
c Omit lines from here ... 
cc wave number ----------------------------------------
c      dimension lkx(nx*ny/2+nx/2+1)
c      dimension lky(nx*ny/2+nx/2+1)
c      dimension lknum(-nx/2+1:nx/2,0:ny/2)
c ... to here when finding only real space physical quantities.

c engv dimension !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dimension gqp(2*nn,2*nn),ginqp(Lmax)
      dimension gdelqp(2*nn,2*nn)
      dimension agqp(2*nn,2*nn),agdelqp(2*nn,2*nn),aqp(2*m,2*m)

      dimension bufqp(m*2,nn*2)
      
      dimension faiqpL(2*nn,2*m,Lmax),faiqpR(2*nn,2*m,lqp)
      dimension IPIVqp(2*m),WORKqp(2*m)
      dimension hqp(Lmax,Lmax),hinqp(Lmax,Lmax)
      dimension wfin(Lmax,Lmax),delcro(2*nn,2*nn)
      
      dimension sq(nn,nn)
      dimension hspin(2*nn,2*nn,lqp),w(lqp),x(lqp),cof(lqp)
      
      
      dimension subh2(Lmax,Lmax)
      dimension subkin1(Lmax,Lmax),subkin2(Lmax,Lmax)
      dimension subkin3(Lmax,Lmax),subhdble(Lmax,Lmax)

      
      dimension subssnn1(Lmax,Lmax),subssnn2(Lmax,Lmax)
      dimension subssnn3(Lmax,Lmax)
      
      dimension subtsq(Lmax,Lmax)

c Omit lines from here ... 
c      dimension uk(nx*ny/2+nx/2+1,2,2)
c      dimension hnij(nn,nn),cdens(nn,nn)
c
c      dimension subnku(nx*ny/2+nx/2+1,Lmax,Lmax)
c      dimension subnkl(nx*ny/2+nx/2+1,Lmax,Lmax)
c      
c      dimension subsq11(nx*ny/2+nx/2+1,Lmax,Lmax)
c      dimension subsq22(nx*ny/2+nx/2+1,Lmax,Lmax)
c      dimension subsqre(nx*ny/2+nx/2+1,Lmax,Lmax)
c      dimension subsqim(nx*ny/2+nx/2+1,Lmax,Lmax)
c
c      dimension subnq11(nx*ny/2+nx/2+1,Lmax,Lmax)
c      dimension subnq22(nx*ny/2+nx/2+1,Lmax,Lmax)
c      dimension subnqre(nx*ny/2+nx/2+1,Lmax,Lmax)
c      dimension subnqim(nx*ny/2+nx/2+1,Lmax,Lmax)
c      
c      
c      dimension hnku(nx*ny/2+nx/2+1),hnkl(nx*ny/2+nx/2+1)
c      dimension hsq11(nx*ny/2+nx/2+1),hsq22(nx*ny/2+nx/2+1)
c      dimension hsqre(nx*ny/2+nx/2+1),hsqim(nx*ny/2+nx/2+1)
c      dimension hsq12(nx*ny/2+nx/2+1)
c      dimension hnq11(nx*ny/2+nx/2+1),hnq22(nx*ny/2+nx/2+1)
c      dimension hnqre(nx*ny/2+nx/2+1),hnqim(nx*ny/2+nx/2+1)
c      dimension hnq12(nx*ny/2+nx/2+1)
c ... to here when finding only real space physical quantities.

c plaquette dimer ----------------------------------
c shift sitenumbers plaquette dimer-------------------
      dimension ldx(nn),ldy(nn)
      dimension ldpx(nn),ldpy(nn)
      dimension ldmx(nn),ldmy(nn)
      dimension psign(0:10)
      dimension ipshift(nn,10,0:10)
      dimension ipnum(0:10)
c ----------------------------------------------------
      dimension pla0(0:10),pla1(0:10)
      dimension hpvbc1(0:10),diag1(0:10),hdiag1(0:10)
      dimension dimxx(nn,nn),dimyy(nn,nn)
      dimension dimxy(nn,nn)
      dimension dim2xy(nn,nn)
      dimension subdimxx1(Lmax,Lmax),subdimxx2(Lmax,Lmax)
      dimension subdimyy1(Lmax,Lmax),subdimyy2(Lmax,Lmax)
      dimension subdimxy1(Lmax,Lmax)
      
      dimension subdim2xy1(Lmax,Lmax)
      
      dimension subpvbc1(Lmax,Lmax,0:10)
      dimension subdiag1(Lmax,Lmax,0:10)
c ------------------------------------------------------




      
      data ITYPE/1/,JOBZ/'N'/,UPLO/'U'/
      
      if(iqpon.eq.0) then
      
      
      open(unit=100,file="t2_1_Eng_vs_Engv2_u6_N24_m.txt"
     &,status='unknown')
      open(unit=101,file="t2_1_Eng_u6_N24_m.txt"
     &,status='unknown')
      open(unit=102,file="t2_1_Dble_u6_N24_m.txt"
     &,status='unknown')
c nearest neighbor spin ccccccccccccccccccccccccccccc
      open(unit=103,file="t2_1_nnspin_sq_u6_N24_m.txt"
     &,status='unknown')
      open(unit=104,file="t2_1_nnspin_diag_u6_N24_m.txt"
     &,status='unknown')
      open(unit=105,file="t2_1_nnspin_tot_u6_N24_m.txt"
     &,status='unknown')
c kenetic energy cccccccccccccccccccccccccccccccccccccccccc
      open(unit=106,file="t2_1_nnhop_sq_u6_N24_m.txt"
     &,status='unknown')
      open(unit=107,file="t2_1_nnhop_diag_u6_N24_m.txt"
     &,status='unknown')
      open(unit=108,file="t2_1_nnhop_tot_u6_N24_m.txt"
     &,status='unknown')
      
      open(unit=109,file="t2_1_Dxx_1_u6_N24_m.txt"
     &,status='unknown')
      open(unit=110,file="t2_1_Dxx_2_u6_N24_m.txt"
     &,status='unknown')
      open(unit=111,file="t2_1_Dyy_1_u6_N24_m.txt"
     &,status='unknown')
      open(unit=112,file="t2_1_Dyy_2_u6_N24_m.txt"
     &,status='unknown')
      open(unit=113,file="t2_1_Dxy_1_u6_N24_m.txt"
     &,status='unknown')
     
      open(unit=114,file="t2_1_C-VBC_1_u6_N24_m.txt"
     &,status='unknown')
      
      open(unit=115,file="t2_1_PLA0-tot_u6_N24_m.txt"
     &,status='unknown')
      open(unit=116,file="t2_1_PLA0-0_u6_N24_m.txt"
     &,status='unknown')
      open(unit=117,file="t2_1_PLA0-1_u6_N24_m.txt"
     &,status='unknown')
      open(unit=118,file="t2_1_PLA0-2_u6_N24_m.txt"
     &,status='unknown')
      open(unit=119,file="t2_1_PLA0-3_u6_N24_m.txt"
     &,status='unknown')
      open(unit=120,file="t2_1_PLA0-4_u6_N24_m.txt"
     &,status='unknown')
      open(unit=121,file="t2_1_PLA0-5_u6_N24_m.txt"
     &,status='unknown')
     
      open(unit=122,file="t2_1_DIAG0-tot_u6_N24_m.txt"
     &,status='unknown')
      open(unit=123,file="t2_1_DIAG0-0_u6_N24_m.txt"
     &,status='unknown')
      open(unit=124,file="t2_1_DIAG0-1_u6_N24_m.txt"
     &,status='unknown')
      open(unit=125,file="t2_1_DIAG0-2_u6_N24_m.txt"
     &,status='unknown')
      open(unit=126,file="t2_1_DIAG0-3_u6_N24_m.txt"
     &,status='unknown')
      open(unit=127,file="t2_1_DIAG0-4_u6_N24_m.txt"
     &,status='unknown')
      open(unit=128,file="t2_1_DIAG0-5_u6_N24_m.txt"
     &,status='unknown')

c Omit lines from here ... 
cc sq cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      open(unit=201,file="t2_1_Sq_max_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=202,file="t2_1_Sq_max_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=203,file="t2_1_Sq_max_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=204,file="t2_1_Sq_max_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=205,file="t2_1_Sq_max_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=206,file="t2_1_Sq_max_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=207,file="t2_1_Sq_max_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=208,file="t2_1_Sq_max_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=209,file="t2_1_Sq_max_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=210,file="t2_1_Sq_max_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=211,file="t2_1_Sq_max_11_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=212,file="t2_1_Sq_min_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=213,file="t2_1_Sq_min_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=214,file="t2_1_Sq_min_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=215,file="t2_1_Sq_min_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=216,file="t2_1_Sq_min_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=217,file="t2_1_Sq_min_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=218,file="t2_1_Sq_min_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=219,file="t2_1_Sq_min_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=220,file="t2_1_Sq_min_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=221,file="t2_1_Sq_min_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=222,file="t2_1_Sq_min_11_u6_N24_m.txt"
c     &,status='unknown')
c
cc nq cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      open(unit=301,file="t2_1_Nq_max_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=302,file="t2_1_Nq_max_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=303,file="t2_1_Nq_max_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=304,file="t2_1_Nq_max_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=305,file="t2_1_Nq_max_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=306,file="t2_1_Nq_max_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=307,file="t2_1_Nq_max_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=308,file="t2_1_Nq_max_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=309,file="t2_1_Nq_max_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=310,file="t2_1_Nq_max_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=311,file="t2_1_Nq_max_11_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=312,file="t2_1_Nq_min_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=313,file="t2_1_Nq_min_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=314,file="t2_1_Nq_min_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=315,file="t2_1_Nq_min_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=316,file="t2_1_Nq_min_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=317,file="t2_1_Nq_min_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=318,file="t2_1_Nq_min_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=319,file="t2_1_Nq_min_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=320,file="t2_1_Nq_min_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=321,file="t2_1_Nq_min_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=322,file="t2_1_Nq_min_11_u6_N24_m.txt"
c     &,status='unknown')
c
cc nk cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      open(unit=401,file="t2_1_nk_max_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=402,file="t2_1_nk_max_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=403,file="t2_1_nk_max_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=404,file="t2_1_nk_max_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=405,file="t2_1_nk_max_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=406,file="t2_1_nk_max_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=407,file="t2_1_nk_max_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=408,file="t2_1_nk_max_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=409,file="t2_1_nk_max_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=410,file="t2_1_nk_max_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=411,file="t2_1_nk_max_11_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=412,file="t2_1_nk_min_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=413,file="t2_1_nk_min_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=414,file="t2_1_nk_min_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=415,file="t2_1_nk_min_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=416,file="t2_1_nk_min_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=417,file="t2_1_nk_min_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=418,file="t2_1_nk_min_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=419,file="t2_1_nk_min_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=420,file="t2_1_nk_min_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=421,file="t2_1_nk_min_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=422,file="t2_1_nk_min_11_u6_N24_m.txt"
c     &,status='unknown')
c ... to here when finding only real space physical quantities.
      
      
      else if(ispin.eq.0) then
      open(unit=100,file="QP_S0_t2_1_Eng_vs_Engv2_u6_N24_m.txt"
     &,status='unknown')
      open(unit=101,file="QP_S0_t2_1_Eng_u6_N24_m.txt"
     &,status='unknown')
      open(unit=102,file="QP_S0_t2_1_Dble_u6_N24_m.txt"
     &,status='unknown')
c nearest neighbor spin ccccccccccccccccccccccccccccc
      open(unit=103,file="QP_S0_t2_1_nnspin_sq_u6_N24_m.txt"
     &,status='unknown')
      open(unit=104,file="QP_S0_t2_1_nnspin_diag_u6_N24_m.txt"
     &,status='unknown')
      open(unit=105,file="QP_S0_t2_1_nnspin_tot_u6_N24_m.txt"
     &,status='unknown')
c kenetic energy cccccccccccccccccccccccccccccccccccccccccc
      open(unit=106,file="QP_S0_t2_1_nnhop_sq_u6_N24_m.txt"
     &,status='unknown')
      open(unit=107,file="QP_S0_t2_1_nnhop_diag_u6_N24_m.txt"
     &,status='unknown')
      open(unit=108,file="QP_S0_t2_1_nnhop_tot_u6_N24_m.txt"
     &,status='unknown')
     
      open(unit=109,file="QP_S0_t2_1_Dxx_1_u6_N24_m.txt"
     &,status='unknown')
      open(unit=110,file="QP_S0_t2_1_Dxx_2_u6_N24_m.txt"
     &,status='unknown')
      open(unit=111,file="QP_S0_t2_1_Dyy_1_u6_N24_m.txt"
     &,status='unknown')
      open(unit=112,file="QP_S0_t2_1_Dyy_2_u6_N24_m.txt"
     &,status='unknown')
      open(unit=113,file="QP_S0_t2_1_Dxy_1_u6_N24_m.txt"
     &,status='unknown')
     
      open(unit=114,file="QP_S0_t2_1_C-VBC_1_u6_N24_m.txt"
     &,status='unknown')
      
      open(unit=115,file="QP_S0_t2_1_PLA0-tot_u6_N24_m.txt"
     &,status='unknown')
      open(unit=116,file="QP_S0_t2_1_PLA0-0_u6_N24_m.txt"
     &,status='unknown')
      open(unit=117,file="QP_S0_t2_1_PLA0-1_u6_N24_m.txt"
     &,status='unknown')
      open(unit=118,file="QP_S0_t2_1_PLA0-2_u6_N24_m.txt"
     &,status='unknown')
      open(unit=119,file="QP_S0_t2_1_PLA0-3_u6_N24_m.txt"
     &,status='unknown')
      open(unit=120,file="QP_S0_t2_1_PLA0-4_u6_N24_m.txt"
     &,status='unknown')
      open(unit=121,file="QP_S0_t2_1_PLA0-5_u6_N24_m.txt"
     &,status='unknown')
     
      open(unit=122,file="QP_S0_t2_1_DIAG0-tot_u6_N24_m.txt"
     &,status='unknown')
      open(unit=123,file="QP_S0_t2_1_DIAG0-0_u6_N24_m.txt"
     &,status='unknown')
      open(unit=124,file="QP_S0_t2_1_DIAG0-1_u6_N24_m.txt"
     &,status='unknown')
      open(unit=125,file="QP_S0_t2_1_DIAG0-2_u6_N24_m.txt"
     &,status='unknown')
      open(unit=126,file="QP_S0_t2_1_DIAG0-3_u6_N24_m.txt"
     &,status='unknown')
      open(unit=127,file="QP_S0_t2_1_DIAG0-4_u6_N24_m.txt"
     &,status='unknown')
      open(unit=128,file="QP_S0_t2_1_DIAG0-5_u6_N24_m.txt"
     &,status='unknown')

c Omit lines from here ... 
cc sq cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      open(unit=201,file="QP_S0_t2_1_Sq_max_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=202,file="QP_S0_t2_1_Sq_max_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=203,file="QP_S0_t2_1_Sq_max_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=204,file="QP_S0_t2_1_Sq_max_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=205,file="QP_S0_t2_1_Sq_max_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=206,file="QP_S0_t2_1_Sq_max_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=207,file="QP_S0_t2_1_Sq_max_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=208,file="QP_S0_t2_1_Sq_max_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=209,file="QP_S0_t2_1_Sq_max_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=210,file="QP_S0_t2_1_Sq_max_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=211,file="QP_S0_t2_1_Sq_max_11_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=212,file="QP_S0_t2_1_Sq_min_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=213,file="QP_S0_t2_1_Sq_min_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=214,file="QP_S0_t2_1_Sq_min_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=215,file="QP_S0_t2_1_Sq_min_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=216,file="QP_S0_t2_1_Sq_min_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=217,file="QP_S0_t2_1_Sq_min_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=218,file="QP_S0_t2_1_Sq_min_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=219,file="QP_S0_t2_1_Sq_min_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=220,file="QP_S0_t2_1_Sq_min_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=221,file="QP_S0_t2_1_Sq_min_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=222,file="QP_S0_t2_1_Sq_min_11_u6_N24_m.txt"
c     &,status='unknown')
c
cc nq cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      open(unit=301,file="QP_S0_t2_1_Nq_max_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=302,file="QP_S0_t2_1_Nq_max_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=303,file="QP_S0_t2_1_Nq_max_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=304,file="QP_S0_t2_1_Nq_max_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=305,file="QP_S0_t2_1_Nq_max_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=306,file="QP_S0_t2_1_Nq_max_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=307,file="QP_S0_t2_1_Nq_max_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=308,file="QP_S0_t2_1_Nq_max_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=309,file="QP_S0_t2_1_Nq_max_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=310,file="QP_S0_t2_1_Nq_max_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=311,file="QP_S0_t2_1_Nq_max_11_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=312,file="QP_S0_t2_1_Nq_min_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=313,file="QP_S0_t2_1_Nq_min_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=314,file="QP_S0_t2_1_Nq_min_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=315,file="QP_S0_t2_1_Nq_min_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=316,file="QP_S0_t2_1_Nq_min_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=317,file="QP_S0_t2_1_Nq_min_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=318,file="QP_S0_t2_1_Nq_min_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=319,file="QP_S0_t2_1_Nq_min_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=320,file="QP_S0_t2_1_Nq_min_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=321,file="QP_S0_t2_1_Nq_min_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=322,file="QP_S0_t2_1_Nq_min_11_u6_N24_m.txt"
c     &,status='unknown')
c
cc nk cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      open(unit=401,file="QP_S0_t2_1_nk_max_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=402,file="QP_S0_t2_1_nk_max_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=403,file="QP_S0_t2_1_nk_max_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=404,file="QP_S0_t2_1_nk_max_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=405,file="QP_S0_t2_1_nk_max_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=406,file="QP_S0_t2_1_nk_max_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=407,file="QP_S0_t2_1_nk_max_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=408,file="QP_S0_t2_1_nk_max_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=409,file="QP_S0_t2_1_nk_max_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=410,file="QP_S0_t2_1_nk_max_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=411,file="QP_S0_t2_1_nk_max_11_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=412,file="QP_S0_t2_1_nk_min_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=413,file="QP_S0_t2_1_nk_min_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=414,file="QP_S0_t2_1_nk_min_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=415,file="QP_S0_t2_1_nk_min_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=416,file="QP_S0_t2_1_nk_min_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=417,file="QP_S0_t2_1_nk_min_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=418,file="QP_S0_t2_1_nk_min_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=419,file="QP_S0_t2_1_nk_min_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=420,file="QP_S0_t2_1_nk_min_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=421,file="QP_S0_t2_1_nk_min_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=422,file="QP_S0_t2_1_nk_min_11_u6_N24_m.txt"
c     &,status='unknown')
c ... to here when finding only real space physical quantities.
     
     
     
      else if(ispin.eq.1) then
      open(unit=100,file="QP_S1_t2_1_Eng_vs_Engv2_u6_N24_m.txt"
     &,status='unknown')
      open(unit=101,file="QP_S1_t2_1_Eng_u6_N24_m.txt"
     &,status='unknown')
      open(unit=102,file="QP_S1_t2_1_Dble_u6_N24_m.txt"
     &,status='unknown')
c nearest neighbor spin ccccccccccccccccccccccccccccc
      open(unit=103,file="QP_S1_t2_1_nnspin_sq_u6_N24_m.txt"
     &,status='unknown')
      open(unit=104,file="QP_S1_t2_1_nnspin_diag_u6_N24_m.txt"
     &,status='unknown')
      open(unit=105,file="QP_S1_t2_1_nnspin_tot_u6_N24_m.txt"
     &,status='unknown')
c kenetic energy cccccccccccccccccccccccccccccccccccccccccc
      open(unit=106,file="QP_S1_t2_1_nnhop_sq_u6_N24_m.txt"
     &,status='unknown')
      open(unit=107,file="QP_S1_t2_1_nnhop_diag_u6_N24_m.txt"
     &,status='unknown')
      open(unit=108,file="QP_S1_t2_1_nnhop_tot_u6_N24_m.txt"
     &,status='unknown')
     
      open(unit=109,file="QP_S1_t2_1_Dxx_1_u6_N24_m.txt"
     &,status='unknown')
      open(unit=110,file="QP_S1_t2_1_Dxx_2_u6_N24_m.txt"
     &,status='unknown')
      open(unit=111,file="QP_S1_t2_1_Dyy_1_u6_N24_m.txt"
     &,status='unknown')
      open(unit=112,file="QP_S1_t2_1_Dyy_2_u6_N24_m.txt"
     &,status='unknown')
      open(unit=113,file="QP_S1_t2_1_Dxy_1_u6_N24_m.txt"
     &,status='unknown')
     
      open(unit=114,file="QP_S1_t2_1_C-VBC_1_u6_N24_m.txt"
     &,status='unknown')
      
      open(unit=115,file="QP_S1_t2_1_PLA0-tot_u6_N24_m.txt"
     &,status='unknown')
      open(unit=116,file="QP_S1_t2_1_PLA0-0_u6_N24_m.txt"
     &,status='unknown')
      open(unit=117,file="QP_S1_t2_1_PLA0-1_u6_N24_m.txt"
     &,status='unknown')
      open(unit=118,file="QP_S1_t2_1_PLA0-2_u6_N24_m.txt"
     &,status='unknown')
      open(unit=119,file="QP_S1_t2_1_PLA0-3_u6_N24_m.txt"
     &,status='unknown')
      open(unit=120,file="QP_S1_t2_1_PLA0-4_u6_N24_m.txt"
     &,status='unknown')
      open(unit=121,file="QP_S1_t2_1_PLA0-5_u6_N24_m.txt"
     &,status='unknown')
     
      open(unit=122,file="QP_S1_t2_1_DIAG0-tot_u6_N24_m.txt"
     &,status='unknown')
      open(unit=123,file="QP_S1_t2_1_DIAG0-0_u6_N24_m.txt"
     &,status='unknown')
      open(unit=124,file="QP_S1_t2_1_DIAG0-1_u6_N24_m.txt"
     &,status='unknown')
      open(unit=125,file="QP_S1_t2_1_DIAG0-2_u6_N24_m.txt"
     &,status='unknown')
      open(unit=126,file="QP_S1_t2_1_DIAG0-3_u6_N24_m.txt"
     &,status='unknown')
      open(unit=127,file="QP_S1_t2_1_DIAG0-4_u6_N24_m.txt"
     &,status='unknown')
      open(unit=128,file="QP_S1_t2_1_DIAG0-5_u6_N24_m.txt"
     &,status='unknown')
     
c Omit lines from here ... 
cc sq cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      open(unit=201,file="QP_S1_t2_1_Sq_max_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=202,file="QP_S1_t2_1_Sq_max_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=203,file="QP_S1_t2_1_Sq_max_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=204,file="QP_S1_t2_1_Sq_max_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=205,file="QP_S1_t2_1_Sq_max_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=206,file="QP_S1_t2_1_Sq_max_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=207,file="QP_S1_t2_1_Sq_max_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=208,file="QP_S1_t2_1_Sq_max_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=209,file="QP_S1_t2_1_Sq_max_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=210,file="QP_S1_t2_1_Sq_max_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=211,file="QP_S1_t2_1_Sq_max_11_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=212,file="QP_S1_t2_1_Sq_min_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=213,file="QP_S1_t2_1_Sq_min_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=214,file="QP_S1_t2_1_Sq_min_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=215,file="QP_S1_t2_1_Sq_min_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=216,file="QP_S1_t2_1_Sq_min_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=217,file="QP_S1_t2_1_Sq_min_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=218,file="QP_S1_t2_1_Sq_min_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=219,file="QP_S1_t2_1_Sq_min_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=220,file="QP_S1_t2_1_Sq_min_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=221,file="QP_S1_t2_1_Sq_min_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=222,file="QP_S1_t2_1_Sq_min_11_u6_N24_m.txt"
c     &,status='unknown')
c
cc nq cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      open(unit=301,file="QP_S1_t2_1_Nq_max_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=302,file="QP_S1_t2_1_Nq_max_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=303,file="QP_S1_t2_1_Nq_max_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=304,file="QP_S1_t2_1_Nq_max_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=305,file="QP_S1_t2_1_Nq_max_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=306,file="QP_S1_t2_1_Nq_max_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=307,file="QP_S1_t2_1_Nq_max_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=308,file="QP_S1_t2_1_Nq_max_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=309,file="QP_S1_t2_1_Nq_max_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=310,file="QP_S1_t2_1_Nq_max_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=311,file="QP_S1_t2_1_Nq_max_11_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=312,file="QP_S1_t2_1_Nq_min_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=313,file="QP_S1_t2_1_Nq_min_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=314,file="QP_S1_t2_1_Nq_min_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=315,file="QP_S1_t2_1_Nq_min_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=316,file="QP_S1_t2_1_Nq_min_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=317,file="QP_S1_t2_1_Nq_min_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=318,file="QP_S1_t2_1_Nq_min_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=319,file="QP_S1_t2_1_Nq_min_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=320,file="QP_S1_t2_1_Nq_min_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=321,file="QP_S1_t2_1_Nq_min_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=322,file="QP_S1_t2_1_Nq_min_11_u6_N24_m.txt"
c     &,status='unknown')
c
cc nk cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      open(unit=401,file="QP_S1_t2_1_nk_max_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=402,file="QP_S1_t2_1_nk_max_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=403,file="QP_S1_t2_1_nk_max_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=404,file="QP_S1_t2_1_nk_max_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=405,file="QP_S1_t2_1_nk_max_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=406,file="QP_S1_t2_1_nk_max_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=407,file="QP_S1_t2_1_nk_max_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=408,file="QP_S1_t2_1_nk_max_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=409,file="QP_S1_t2_1_nk_max_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=410,file="QP_S1_t2_1_nk_max_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=411,file="QP_S1_t2_1_nk_max_11_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=412,file="QP_S1_t2_1_nk_min_1_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=413,file="QP_S1_t2_1_nk_min_2_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=414,file="QP_S1_t2_1_nk_min_3_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=415,file="QP_S1_t2_1_nk_min_4_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=416,file="QP_S1_t2_1_nk_min_5_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=417,file="QP_S1_t2_1_nk_min_6_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=418,file="QP_S1_t2_1_nk_min_7_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=419,file="QP_S1_t2_1_nk_min_8_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=420,file="QP_S1_t2_1_nk_min_9_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=421,file="QP_S1_t2_1_nk_min_10_u6_N24_m.txt"
c     &,status='unknown')
c      open(unit=422,file="QP_S1_t2_1_nk_min_11_u6_N24_m.txt"
c     &,status='unknown')
c ... to here when finding only real space physical quantities.
      
      
      end if

      call shiftpd(nn,ihop1,ihop2,ldx,ldy
     &         ,ldpx,ldpy,ldmx,ldmy,psign,ipshift,ipnum)

c Omit lines from here ... 
cc set wave number ----------------------------------------------------
c      lkx=0
c      lky=0
c      lknum=0
c      do i=0,nx/2,1
c         lknum(i,0)=i+1
c      end do
c      do j=1,ny/2
c         do i=-nx/2+1,nx/2,1
c            lknum(i,j)=i+1+j*nx
c         end do
c      end do
c      
c      do iy=0,ny/2,1
c         do ix=-nx/2+1,nx/2,1
c         if(lknum(ix,iy).ne.0) then
c            lkx(lknum(ix,iy))=ix
c            lky(lknum(ix,iy))=iy
c         end if
c         end do
c      end do
c      
cc end set wave number ---------------------------------------------------
c ... to here when finding only real space physical quantities.

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
      

      m2=2*m

      pai=dacos(-1.0d0)
      


      if(iqpon.eq.1) then
         write(6,*) '+QP ON','    TOTAL S=',ispin
      else
         write(6,*) '+QP OFF'
      end if

c Omit lines from here ... 
cc make unitary matrix cccccccccccccccccccccccccccccccccccc
c      do k=1,nx*ny/2+nx/2+1
c      kx=lkx(k)
c      ky=lky(k)
c
c      akx=2.0d0*pai*dble(kx)/dble(nx)
c      aky=2.0d0*pai*dble(ky)/dble(ny)
c      
c         h11=-t2*2.0d0*dcos(aky)
c         h22=-t2*2.0d0*dcos(akx)
c         h12=-t1*4.0d0*dcos(akx/2.0d0)*dcos(aky/2.0d0)
c         h21=-t1*4.0d0*dcos(akx/2.0d0)*dcos(aky/2.0d0)
c      
c      d=(h11-h22)**2+4.0d0*h12*h21
c      
c      uk(k,1,1)=(1.0d0-(h11-h22)/dsqrt(d))/2.0d0
c      uk(k,1,2)=(1.0d0+(h11-h22)/dsqrt(d))/2.0d0
cc remove singlarity cccccccccccccccccccccccccccccccccccccc
c      sign=-1.0d0
c      if(uk(k,1,1).lt.1.0d-16) then
c         uk(k,1,1)=0.0d0
c      else if(uk(k,1,1).eq.1.0d0) then
c         uk(k,1,1)=1.0d0
c      else 
c         uk(k,1,1)=dsqrt(uk(k,1,1))  
c      end if   
c      
c      if(uk(k,1,2).lt.1.0d-16) then
c         uk(k,1,2)=0.0d0
c      else if(uk(k,1,2).eq.1.0d0) then
c         uk(k,1,2)=sign*1.0d0
c      else 
c         uk(k,1,2)=sign*dsqrt(uk(k,1,2))
c      end if   
c      
c      uk(k,2,1)=-uk(k,1,2)
c      uk(k,2,2)=uk(k,1,1)
c      
c      if(kx.eq.nx/2.and.ky.eq.ny/2) then
c         uk(k,1,1)=1.0/dsqrt(2.0d0)
c         uk(k,2,1)=1.0/dsqrt(2.0d0)
c         uk(k,1,2)=-1.0/dsqrt(2.0d0)
c         uk(k,2,2)=1.0/dsqrt(2.0d0)
c      end if
c      end do
cc end remove singlarity ccccccccccccccccccccccccccccccccccccccc      
cc end make unitary matrix ccccccccccccccccccccccccccccccccccc
cc ------------------------------------------------------------------
c ... to here when finding only real space physical quantities.



      
      call spinrt(nn,m,hspin,lqp,cof,w,x,ispin)
      if(iqpon.eq.0) cof(1)=1.0d0
      
      hqp=0.0d0
      hinqp=0.0d0
      
      subhdble=0.0d0
      subh2=0.0d0
      

      subnku=0.0d0
      subnkl=0.0d0

      subtsq=0.0d0
       
      subsq11=0.0d0
      subsq22=0.0d0
      subsqre=0.0d0
      subsqim=0.0d0
      
      subnq11=0.0d0
      subnq22=0.0d0
      subnqre=0.0d0
      subnqim=0.0d0


      subkin1=0.0d0
      subkin2=0.0d0
      subkin3=0.0d0
      subssnn1=0.0d0
      subssnn2=0.0d0
      subssnn3=0.0d0
      subdimxx1=0.0d0
      subdimxx2=0.0d0
      subdimyy1=0.0d0
      subdimyy2=0.0d0
      subdimxy1=0.0d0
      
      subdim2xy1=0.0d0
      
      subpvbc1=0.0d0
      subdiag1=0.0d0

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
     &,aqp,faiqpL,faiqpR,gqp,ginqp,bufqp)
      
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
      subhdble(ll,lr)=subhdble(ll,lr)
     &+(gqp(j,j)*gqp(j+nn,j+nn)
     &+gqp(j,j+nn)*agqp(j,j+nn))*hin(ll,lr)*cof(iqp)
      end do
      
      
            
      dimxx=0.0d0
      dimyy=0.0d0
      dimxy=0.0d0
      
      dim2xy=0.0d0
      pla1=0.0d0
      diag1=0.0d0
      
      h2=0.0d0
      hnij=0.0d0
      cdens=0.0d0
      sq=0.0d0
      
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

c dimer ----------------------------------------------------------------
      call dimer(nn,ll,lr,i,j,lqp,iqp
     &,gqp,agqp,ldx,ldy,awick)
      dimxy(i,j)=dimxy(i,j)+awick
      call dimer(nn,ll,lr,j,i,lqp,iqp
     &,gqp,agqp,ldy,ldx,awick)
      dimxy(i,j)=dimxy(i,j)+awick

      call dimer(nn,ll,lr,i,j,lqp,iqp
     &,gqp,agqp,ldx,ldx,awick)
       dimxx(i,j)=dimxx(i,j)+awick

      call dimer(nn,ll,lr,i,j,lqp,iqp
     &,gqp,agqp,ldy,ldy,awick)
       dimyy(i,j)=dimyy(i,j)+awick

c C-VBC pattern -------------------------------------------------      
      if(lsub(i).eq.1) then
      
      if(lsub(j).eq.1) then
      call dimer(nn,ll,lr,i,j,lqp,iqp
     &,gqp,agqp,ldpy,ldpy,awick)
      dim2xy(i,j)=dim2xy(i,j)+awick
      else
      call dimer(nn,ll,lr,i,j,lqp,iqp
     &,gqp,agqp,ldpy,ldpx,awick)
      dim2xy(i,j)=dim2xy(i,j)+awick
      end if
      
      else
      
      if(lsub(j).eq.1) then
      call dimer(nn,ll,lr,i,j,lqp,iqp
     &,gqp,agqp,ldpx,ldpy,awick)
      dim2xy(i,j)=dim2xy(i,j)+awick
      else
      call dimer(nn,ll,lr,i,j,lqp,iqp
     &,gqp,agqp,ldpx,ldpx,awick)
      dim2xy(i,j)=dim2xy(i,j)+awick
      end if
      
      end if
c C-VBC pattern -----------------------------------------------------      
      
      
c Omit lines from here ... 
cc nk ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      hnij(i,j)=hnij(i,j)
c     &+gqp(i,j)+gqp(i+nn,j+nn)
c     &+gqp(j,i)+gqp(j+nn,i+nn)
c
cc cdens (nq) ccccccccccccccccccccccccccccccccccccccccccccc
c      cdens(i,j)=cdens(i,j)
c     &+gqp(i,j)*agqp(i,j)
c     &+gqp(i,j+nn)*agqp(i,j+nn)
c     &+gqp(i+nn,j)*agqp(i+nn,j)
c     &+gqp(i+nn,j+nn)*agqp(i+nn,j+nn)
c ... to here when finding only real space physical quantities.

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



c P-VBC ----------------------------------------------------------
      do 1444 i=1,nn
      if(lsub(i).eq.1) then
         
      pla0=0.0d0
      
      do nnsite=0,10
         do ihand=1,10
            if(ipshift(i,ihand,nnsite).ne.0) then
               j=ipshift(i,ihand,nnsite)
      call cPLA(nn,Lmax,ll,lr,i,j,lqp,iqp,pla0,nnsite
     &,psign,gqp,agqp,ldx,ldy)
            end if
            
         end do
      end do
     
      do l5=0,5
         pla1(l5)=pla1(l5)+pla0(l5)
      end do

c DIAG
      else  

      pla0=0.0d0
      
      do nnsite=0,10
         do ihand=1,10
            if(ipshift(i,ihand,nnsite).ne.0) then
               j=ipshift(i,ihand,nnsite)
      call cPLA(nn,Lmax,ll,lr,i,j,lqp,iqp,pla0,nnsite
     &,psign,gqp,agqp,ldx,ldy)
            end if
            
         end do
      end do
      
      do l5=0,5
         diag1(l5)=diag1(l5)+pla0(l5)
      end do
      end if
1444  continue
c end P-VBC ------------------------------------------------------






      do i=1,nn
      do k=1,4
      subssnn1(ll,lr)=subssnn1(ll,lr)
     &+sq(i,ihop1(i,k))*hin(ll,lr)*cof(iqp)
      subkin1(ll,lr)=subkin1(ll,lr)
     &-(gqp(i,ihop1(i,k))+gqp(i+nn,ihop1(i,k)+nn))
     &*hin(ll,lr)*cof(iqp)
      end do
      do k=1,2
      subssnn1(ll,lr)=subssnn1(ll,lr)
     &+sq(i,ihop2(i,k))*hin(ll,lr)*cof(iqp)
      subkin1(ll,lr)=subkin1(ll,lr)
     &-t2*(gqp(i,ihop2(i,k))+gqp(i+nn,ihop2(i,k)+nn))
     &*hin(ll,lr)*cof(iqp)
      end do
      
      do k=1,4
      subssnn2(ll,lr)=subssnn2(ll,lr)
     &+sq(i,ihop1(i,k))*hin(ll,lr)*cof(iqp)
      subkin2(ll,lr)=subkin2(ll,lr)
     &-(gqp(i,ihop1(i,k))+gqp(i+nn,ihop1(i,k)+nn))
     &*hin(ll,lr)*cof(iqp)
      end do
      
      do k=1,2
      subssnn3(ll,lr)=subssnn3(ll,lr)
     &+sq(i,ihop2(i,k))*hin(ll,lr)*cof(iqp)
      subkin3(ll,lr)=subkin3(ll,lr)
     &-(gqp(i,ihop2(i,k))+gqp(i+nn,ihop2(i,k)+nn))
     &*hin(ll,lr)*cof(iqp)
      end do      
      end do
      
      do j=1,nn
         do i=1,nn
         ix=lsitex(i)
         iy=lsitey(i)
         jx=lsitex(j)
         jy=lsitey(j)
      subdimxx1(ll,lr)=subdimxx1(ll,lr)+(
     &+dimxx(i,j)*(-1.0d0)**(abs(ix-jx))
     &)*hin(ll,lr)*cof(iqp)
      
      subdimxx2(ll,lr)=subdimxx2(ll,lr)+(
     &+dimxx(i,j)*(-1.0d0)**(abs(ix-jx)+abs(iy-jy))
     &)*hin(ll,lr)*cof(iqp)
     
      subdimyy1(ll,lr)=subdimyy1(ll,lr)+(
     &+dimyy(i,j)*(-1.0d0)**(abs(iy-jy))
     &)*hin(ll,lr)*cof(iqp)
     
      subdimyy2(ll,lr)=subdimyy2(ll,lr)+(
     &+dimyy(i,j)*(-1.0d0)**(abs(ix-jx)+abs(iy-jy))
     &)*hin(ll,lr)*cof(iqp)
     
     
      subdimxy1(ll,lr)=subdimxy1(ll,lr)+(
     &-dimxy(i,j)*(-1.0d0)**(abs(ix)+abs(jy))
     &)*hin(ll,lr)*cof(iqp)
     
c C-VBC pattern
      subdim2xy1(ll,lr)=subdim2xy1(ll,lr)+(
     &+dim2xy(i,j)*(-1.0d0)**(abs(iy-jy))
     &)*hin(ll,lr)*cof(iqp)
         end do
      end do

c P-VBC pattern     
      do l5=0,5
      subpvbc1(ll,lr,l5)=subpvbc1(ll,lr,l5)
     &+pla1(l5)*hin(ll,lr)*cof(iqp)
      subdiag1(ll,lr,l5)=subdiag1(ll,lr,l5)
     &+diag1(l5)*hin(ll,lr)*cof(iqp)
      end do


c total S -------------------------------------------------
      do j=1,nn
      do i=1,nn
      subtsq(ll,lr)=subtsq(ll,lr)
     &+sq(i,j)*hin(ll,lr)*cof(iqp)
      end do
      end do
c --------------------------------------------------------------      

c Omit lines from here ... 
cc nk,Nq,Sq -------------------------------------------------
c      do k=1,nx*ny/2+nx/2+1
c      kx=lkx(k)
c      ky=lky(k)
c      akx=2.0d0*pai*dble(kx)/dble(nx)
c      aky=2.0d0*pai*dble(ky)/dble(ny)
c      
c      do j=1,nx*ny*2
c         ixr=j-(j-1)/nx*nx
c         iyr=1+(j-1)/nx
c         if(j.le.nx*ny) then
c            pxr=0.0d0
c            pyr=0.0d0
c            mr=1
c         else if(j.le.nx*ny*2) then
c            iyr=iyr-ny
c            pxr=0.5d0
c            pyr=-0.5d0
c            mr=2
c         end if
c
c      do i=1,nx*ny*2
c         ixl=i-(i-1)/nx*nx
c         iyl=1+(i-1)/nx
c         if(i.le.nx*ny) then
c            pxl=0.0d0
c            pyl=0.0d0
c            ml=1
c         else if(i.le.nx*ny*2) then
c            iyl=iyl-ny
c            pxl=0.5d0
c            pyl=-0.5d0
c            ml=2
c         end if
c
c      subnku(k,ll,lr)=subnku(k,ll,lr)
c     &+hnij(i,j)*uk(k,ml,2)*uk(k,mr,2)*dcos(
c     &akx*(dble(ixl)+pxl-dble(ixr)-pxr)
c     &+aky*(dble(iyl)+pyl-dble(iyr)-pyr))*hin(ll,lr)*cof(iqp)
c     
c      subnkl(k,ll,lr)=subnkl(k,ll,lr)
c     &+hnij(i,j)*uk(k,ml,1)*uk(k,mr,1)*dcos(
c     &akx*(dble(ixl)+pxl-dble(ixr)-pxr)
c     &+aky*(dble(iyl)+pyl-dble(iyr)-pyr))*hin(ll,lr)*cof(iqp)
c     
c      if(i.le.nx*ny.and.j.le.nx*ny) then
c      subnq11(k,ll,lr)=subnq11(k,ll,lr)
c     &+cdens(i,j)*dcos(
c     &akx*(dble(ixl)-dble(ixr))
c     &+aky*(dble(iyl)-dble(iyr)))*hin(ll,lr)*cof(iqp)
c     
c      subsq11(k,ll,lr)=subsq11(k,ll,lr)
c     &+sq(i,j)*dcos(
c     &akx*(dble(ixl)-dble(ixr))
c     &+aky*(dble(iyl)-dble(iyr)))*hin(ll,lr)*cof(iqp)
c     
c      else if(i.ge.nx*ny+1.and.j.ge.nx*ny+1) then
c      subnq22(k,ll,lr)=subnq22(k,ll,lr)
c     &+cdens(i,j)*dcos(
c     &akx*(dble(ixl)-dble(ixr))
c     &+aky*(dble(iyl)-dble(iyr)))*hin(ll,lr)*cof(iqp)
c     
c      subsq22(k,ll,lr)=subsq22(k,ll,lr)
c     &+sq(i,j)*dcos(
c     &akx*(dble(ixl)-dble(ixr))
c     &+aky*(dble(iyl)-dble(iyr)))*hin(ll,lr)*cof(iqp)
c     
c      else if(i.le.nx*ny.and.j.ge.nx*ny+1) then
c      subnqre(k,ll,lr)=subnqre(k,ll,lr)
c     &+cdens(i,j)*dcos(
c     &akx*(dble(ixl)-dble(ixr))
c     &+aky*(dble(iyl)-dble(iyr)))*hin(ll,lr)*cof(iqp)
c      subnqim(k,ll,lr)=subnqim(k,ll,lr)
c     &+cdens(i,j)*dsin(
c     &akx*(dble(ixl)-dble(ixr))
c     &+aky*(dble(iyl)-dble(iyr)))*hin(ll,lr)*cof(iqp)
c     
c      subsqre(k,ll,lr)=subsqre(k,ll,lr)
c     &+sq(i,j)*dcos(
c     &akx*(dble(ixl)-dble(ixr))
c     &+aky*(dble(iyl)-dble(iyr)))*hin(ll,lr)*cof(iqp)
c      subsqim(k,ll,lr)=subsqim(k,ll,lr)
c     &+sq(i,j)*dsin(
c     &akx*(dble(ixl)-dble(ixr))
c     &+aky*(dble(iyl)-dble(iyr)))*hin(ll,lr)*cof(iqp)
c     
c      end if
c      end do
c      end do
c      end do
cc --------------------------------------------------------------      
c ... to here when finding only real space physical quantities.

1711        continue
1911     continue
1811  continue
c end initial eng -------------------------------------------------
      
      
      
      




c subspace EXT --------------------------------------------------
      
      call geed0(Lmax,iLdim,hqp,hinqp
     &,hgeed,hingeed,deigen,eigen,WORK2,LWORK2,'V',label)
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
      
      hdimxx1=0.0d0
      hdimxx2=0.0d0
      hdimyy1=0.0d0
      hdimyy2=0.0d0
      hdimxy1=0.0d0
      
      hdim2xx1=0.0d0
      hdim2xx2=0.0d0
      hdim2yy1=0.0d0
      hdim2yy2=0.0d0
      hdim2xy1=0.0d0

      hpvbc1=0.0d0
      hdiag1=0.0d0

c additional alignment -------------------------------------------
      hnku=0.0d0
      hnkl=0.0d0
      
      tsq=0.0d0
      
      hsq11=0.0d0
      hsq22=0.0d0
      hsq12=0.0d0
      hsqre=0.0d0
      hsqim=0.0d0
      
      hnq11=0.0d0
      hnq22=0.0d0
      hnq12=0.0d0
      hnqre=0.0d0
      hnqim=0.0d0
c --------------------------------------------------------------------
      
      do lr=1,iLdim
         do ll=1,lr
      if(ll.ne.lr) then
         cweight=2.0d0*hgeed(ll,1)*hgeed(lr,1)
      else
         cweight=hgeed(ll,1)*hgeed(lr,1)
      end if
      hh2=hh2+subh2(ll,lr)*cweight
      hdble=hdble+subhdble(ll,lr)*cweight
      hkin1=hkin1+subkin1(ll,lr)*cweight
      hkin2=hkin2+subkin2(ll,lr)*cweight
      hkin3=hkin3+subkin3(ll,lr)*cweight
      
      hssnn1=hssnn1+subssnn1(ll,lr)*cweight
      hssnn2=hssnn2+subssnn2(ll,lr)*cweight
      hssnn3=hssnn3+subssnn3(ll,lr)*cweight
      
      hdimxx1=hdimxx1+subdimxx1(ll,lr)*cweight
      hdimxx2=hdimxx2+subdimxx2(ll,lr)*cweight
      hdimyy1=hdimyy1+subdimyy1(ll,lr)*cweight
      hdimyy2=hdimyy2+subdimyy2(ll,lr)*cweight 
      hdimxy1=hdimxy1+subdimxy1(ll,lr)*cweight
      
      hdim2xy1=hdim2xy1+subdim2xy1(ll,lr)*cweight
      
      do l5=0,5
      hpvbc1(l5)=hpvbc1(l5)+subpvbc1(ll,lr,l5)*cweight
      hdiag1(l5)=hdiag1(l5)+subdiag1(ll,lr,l5)*cweight
      end do
      
      tsq=tsq+subtsq(ll,lr)*cweight


c Omit lines from here ... 
c      do k=1,nx*ny/2+nx/2+1
c      hnku(k)=hnku(k)+subnku(k,ll,lr)*cweight
c      hnkl(k)=hnkl(k)+subnkl(k,ll,lr)*cweight
c      
c      hsq11(k)=hsq11(k)+subsq11(k,ll,lr)*cweight
c      hsq22(k)=hsq22(k)+subsq22(k,ll,lr)*cweight
c      hsqre(k)=hsqre(k)+subsqre(k,ll,lr)*cweight
c      hsqim(k)=hsqim(k)+subsqim(k,ll,lr)*cweight
c      
c      hnq11(k)=hnq11(k)+subnq11(k,ll,lr)*cweight
c      hnq22(k)=hnq22(k)+subnq22(k,ll,lr)*cweight
c      hnqre(k)=hnqre(k)+subnqre(k,ll,lr)*cweight
c      hnqim(k)=hnqim(k)+subnqim(k,ll,lr)*cweight
c      end do
c ... to here when finding only real space physical quantities.
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
      
      hdimxx1=hdimxx1/hinre/dble(nn)/dble(nn)
      hdimxx2=hdimxx2/hinre/dble(nn)/dble(nn)
      hdimyy1=hdimyy1/hinre/dble(nn)/dble(nn)
      hdimyy2=hdimyy2/hinre/dble(nn)/dble(nn)
      hdimxy1=hdimxy1/hinre/dble(nn)/dble(nn)*0.5d0
      
      hdim2xy1=hdim2xy1/hinre/dble(nn)/dble(nn)
      
      pvbc1tot=0.0d0
      diag1tot=0.0d0
      
      do l5=0,5
      pvbc1tot=pvbc1tot+hpvbc1(l5)
      diag1tot=diag1tot+hdiag1(l5)
      end do
      pvbc1tot=pvbc1tot/hinre/dble(nn/2)/dble(nn/2)
      diag1tot=diag1tot/hinre/dble(nn/2)/dble(nn/2)
      
      
      do i=1,5
      if(ipnum(i).ne.0) then
         hpvbc1(i)=hpvbc1(i)/hinre/dble(ipnum(i))
         hdiag1(i)=hdiag1(i)/hinre/dble(ipnum(i))
      end if
      end do
      
      
      do l5=0,5
         hpvbc1(l5)=hpvbc1(l5)/dble(nn/2)
         hdiag1(l5)=hdiag1(l5)/dble(nn/2)
      end do
      
      tsq=tsq/hinre

c Omit lines from here ... 
c      do k=1,nx*ny/2+nx/2+1
c      hnku(k)=hnku(k)*0.25d0/dble(nx*ny)/hinre
c      hnkl(k)=hnkl(k)*0.25d0/dble(nx*ny)/hinre
c      
c      hsq11(k)=hsq11(k)/dble(nx*ny)/3.0d0/hinre
c      hsq22(k)=hsq22(k)/dble(nx*ny)/3.0d0/hinre
c      hsq12(k)
c     &=dsqrt(hsqre(k)**2+hsqim(k)**2)/dble(nx*ny)/3.0d0/hinre
c      
c      hnq11(k)=hnq11(k)/dble(nx*ny)/hinre
c      hnq22(k)=hnq22(k)/dble(nx*ny)/hinre
c      hnq12(k)
c     &=dsqrt(hnqre(k)**2+hnqim(k)**2)/dble(nx*ny)/hinre
c      end do
c ... to here when finding only real space physical quantities.
      

      write(100,*) engva2,pirgeng,iLdim
      write(101,*) engva,pirgeng,iLdim
      write(102,*) engva,hdble,iLdim
      write(103,*) engva,hssnn2,iLdim
      write(104,*) engva,hssnn3,iLdim
      write(105,*) engva,hssnn1,iLdim
      write(106,*) engva,hkin2,iLdim
      write(107,*) engva,hkin3,iLdim
      write(108,*) engva,hkin1,iLdim
      
      write(109,*) engva,hdimxx1,iLdim
      write(110,*) engva,hdimxx2,iLdim
      write(111,*) engva,hdimyy1,iLdim
      write(112,*) engva,hdimyy2,iLdim
      write(113,*) engva,hdimxy1,iLdim
      
      write(114,*) engva,hdim2xy1,iLdim

      write(115,*) engva,pvbc1tot,iLdim
      write(116,*) engva,hpvbc1(0),iLdim
      write(117,*) engva,-hpvbc1(1),iLdim
      write(118,*) engva,hpvbc1(2),iLdim
      write(119,*) engva,hpvbc1(3),iLdim
      write(120,*) engva,-hpvbc1(4),iLdim
      write(121,*) engva,hpvbc1(5),iLdim
      

      write(122,*) engva,diag1tot,iLdim
      write(123,*) engva,hdiag1(0),iLdim
      write(124,*) engva,-hdiag1(1),iLdim
      write(125,*) engva,hdiag1(2),iLdim
      write(126,*) engva,hdiag1(3),iLdim
      write(127,*) engva,-hdiag1(4),iLdim
      write(128,*) engva,hdiag1(5),iLdim


      write(6,*)
      write(6,*) 'total spin',tsq

c Omit lines from here ... 
c      write(6,*)
c      write(6,*) 'total spin2',
c     &(hsq11(1)+hsq22(1)-2.0d0*hsq12(1))*3.0d0*dble(nx*ny)
c      a=0.0d0
c      kdel=nx*ny/2+nx/2+1
c      do k=1,kdel
c      
c      sqmax=(hsq11(k)+hsq22(k)
c     &+dsqrt((hsq11(k)-hsq22(k))**2
c     &+4.0d0*hsq12(k)**2))*0.5d0
c      sqmin=(hsq11(k)+hsq22(k)
c     &-dsqrt((hsq11(k)-hsq22(k))**2
c     &+4.0d0*hsq12(k)**2))*0.5d0
c      
c
c      hnqmax=(hnq11(k)+hnq22(k)
c     &+dsqrt((hnq11(k)-hnq22(k))**2
c     &+4.0d0*hnq12(k)**2))*0.5d0
c      hnqmin=(hnq11(k)+hnq22(k)
c     &-dsqrt((hnq11(k)-hnq22(k))**2
c     &+4.0d0*hnq12(k)**2))*0.5d0
c
c         write(400+k,*) engva,hnku(k),iLdim
c         write(400+kdel+k,*) engva,hnkl(k),iLdim
c         if(k.eq.1.or.lkx(k).eq.nx/2.and.lky(k).eq.0
c     &                          .or.lky(k).eq.ny/2) then
c            a=a+hnku(k)+hnkl(k)
c         else
c            a=a+(hnku(k)+hnkl(k))*2.0d0
c         end if
c         write(300+k,*) engva,hnqmax,iLdim
c         write(300+kdel+k,*) engva,hnqmin,iLdim
c         write(200+k,*) engva,sqmax,iLdim
c         write(200+kdel+k,*) engva,sqmin,iLdim
c      end do
c      write(6,*) 'particle number',a*2.0d0/dble(nn)
c ... to here when finding only real space physical quantities.



      return
      end
c end PIRG checker board
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
      
      
      
      
      
      
      
      
      
      
      
      
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine shift(nn,ihop1,ihop2,lsitex,lsitey,lsub)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
c shift site numbers cccccccccccccccccccccccccccc
      dimension ihop1(nn,4),ihop2(nn,2)
      dimension lsitex(nn),lsitey(nn),lsub(nn)
      dimension lsitenum(10,10)
      dimension ip(10),im(10)
ccccccccccccccccccccccccccccccccccccccccccccccccc
      
      
      
      lsitenum=0.0d0
      do i=1,10
         ip(i)=i+1
         im(i)=i-1
      end do
c set site num.
      lsitenum(2,2)=1
      lsitenum(3,2)=2
      lsitenum(4,2)=3
      lsitenum(5,2)=4
      lsitenum(6,2)=5
      lsitenum(7,2)=6

      lsitenum(2,3)=7
      lsitenum(3,3)=8
      lsitenum(4,3)=9
      lsitenum(5,3)=10
      lsitenum(6,3)=11
      lsitenum(7,3)=12
      
      lsitenum(2,4)=13
      lsitenum(3,4)=14
      lsitenum(4,4)=15
      lsitenum(5,4)=16
      lsitenum(6,4)=17
      lsitenum(7,4)=18
      
      lsitenum(2,5)=19
      lsitenum(3,5)=20
      lsitenum(4,5)=21
      lsitenum(5,5)=22
      lsitenum(6,5)=23
      lsitenum(7,5)=24
      
      do ix=1,10
      do iy=1,10
         if(lsitenum(ix,iy).ne.0) then
            lsitex(lsitenum(ix,iy))=ix
            lsitey(lsitenum(ix,iy))=iy
         end if
      end do
      end do
      
      lsitenum(1,1)=24
      lsitenum(2,1)=19
      lsitenum(3,1)=20
      lsitenum(4,1)=21
      lsitenum(5,1)=22
      lsitenum(6,1)=23
      lsitenum(7,1)=24
      lsitenum(8,1)=19
      
      lsitenum(1,6)=6
      lsitenum(2,6)=1
      lsitenum(3,6)=2
      lsitenum(4,6)=3
      lsitenum(5,6)=4
      lsitenum(6,6)=5
      lsitenum(7,6)=6
      lsitenum(8,6)=1
      
      lsitenum(1,2)=6
      lsitenum(1,3)=12
      lsitenum(1,4)=18
      lsitenum(1,5)=24
      
      lsitenum(8,2)=1
      lsitenum(8,3)=7
      lsitenum(8,4)=13
      lsitenum(8,5)=19
      
      
c      do j=1,10
c      do i=1,10
c      if(lsitenum(i,j).ne.0) then
c      write(6,*) i,j,lsitenum(i,j)
c      end if
c      end do
c      end do
c      do i=1,nn
c         write(6,*) i,lsitex(i),lsitey(i)
c      end do
c      stop

c set sublattice num.
c !!! note numbering of 1 or 2
      do i=1,nn
         if(mod(lsitex(i)+lsitey(i),2).eq.0) then
            lsub(i)=2
         else
            lsub(i)=1
         end if
      end do
c      do i=1,nn
c      write(6,*) i,lsub(i)
c      end do
c      stop


      
      do i=1,nn
      
      ix=lsitex(i)
      iy=lsitey(i)
         ihop1(i,1)=lsitenum(ip(ix),iy)
         ihop1(i,2)=lsitenum(ix,ip(iy))
         ihop1(i,3)=lsitenum(im(ix),iy)
         ihop1(i,4)=lsitenum(ix,im(iy))
         
      if(lsub(i).eq.2) then
         ihop2(i,1)=lsitenum(ip(ix),ip(iy))
         ihop2(i,2)=lsitenum(im(ix),im(iy))
      else
         ihop2(i,1)=lsitenum(im(ix),ip(iy))
         ihop2(i,2)=lsitenum(ip(ix),im(iy))
      end if
         
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
c shift site nubers 2 ccccccccccccccccccccccccccccccc

      return
      end
      













c --------------------------------------------------------------------------
      subroutine shiftpd(nn,ihop1,ihop2,ldx,ldy
     &         ,ldpx,ldpy,ldmx,ldmy,psign,ipshift,ipnum)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      dimension ihop1(nn,4),ihop2(nn,2)
      
      dimension ldx(nn),ldy(nn)
      dimension ldpx(nn),ldpy(nn)
      dimension ldmx(nn),ldmy(nn)
      dimension psign(0:10)
      dimension ipshift(nn,10,0:10)
      dimension ipnum(0:10)
c set sublattice num.


c shift site nubers 2 ccccccccccccccccccccccccccccccc
      do i=1,nn
         ldx(i)=ihop1(i,1)
         ldy(i)=ihop1(i,2)
      end do
      
      do i=1,nn
         ldpx(i)=ihop1(ihop1(i,1),2)
         ldpy(i)=ihop1(ihop1(i,3),2)
         
         ldmx(i)=ihop1(ihop1(i,3),4)
         ldmy(i)=ihop1(ihop1(i,1),4)
      end do
     
c      do i=1,nn
c         write(6,*) i,ldpx(i)
c         write(6,*) i,ldpy(i)
c         write(6,*) i,ldmx(i)
c         write(6,*) i,ldmy(i)
c         write(6,*)
c      end do
c      stop
c end shift site nubers 2 ccccccccccccccccc

c end shift site nubers ccccccccccccccccccccccccccccc





c shift plaquette numbers cccccccccccccccccccccccccccccccccccccc
c PLA sign ------------------------------------------------------
      ipshift=0
      psign=-1.0d0

      do nnsite=0,5
      if(nnsite.eq.0.or.nnsite.eq.2
     &.or.nnsite.eq.3.or.nnsite.eq.5) then
      psign(nnsite)=1.0d0
      end if
      end do
c ---------------------------------------------------------------
      do i=1,nn
      
         ipshift(i,1,0)=i
         ipnum(0)=1
         
         ipshift(i,1,1)=ldpx(i)
         ipshift(i,2,1)=ldpy(i)
         ipshift(i,3,1)=ldmx(i)
         ipshift(i,4,1)=ldmy(i)
         ipnum(1)=4
         
         ipshift(i,1,2)=ldpy(ldpx(i))
         ipshift(i,2,2)=ldpy(ldmx(i))
         ipshift(i,3,2)=ldmy(ldpx(i))
         ipnum(2)=3
         
         ipshift(i,1,3)=ldpx(ldpx(i))
         ipshift(i,2,3)=ldpy(ldpy(i))
         ipnum(3)=2
         
         ipshift(i,1,4)=ldmy(ldpx(ldpx(i)))
         ipshift(i,2,4)=ldmy(ldmy(ldpx(i)))
         ipnum(4)=2
      end do
      
c      do i=1,nn
c      do nnsite=0,10
c         if(ipnum(nnsite).ne.0) then
c         write(6,*) ipnum(nnsite)
c         end if
c         do ihand=1,10
c            if(ipshift(i,ihand,nnsite).ne.0) then
c               write(6,*) i,ipshift(i,ihand,nnsite)
c            end if
c         end do
c      end do
c      end do
c      stop
c end shift plaquette numbers ccccccccccccccccccccccccccccccccc


      return
      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc





      
      
      
     
     
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




ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cPLA(nn,Lmax,ll,lr,i,j,lqp,iqp,pla0,nnsite
     &,psign,gqp,agqp,ldx,ldy)

      implicit real*8(a-h,o-z)
      implicit integer*4(i-n) 

      dimension gqp(nn*2,nn*2),agqp(nn*2,nn*2)
      dimension ldx(nn),ldy(nn)
      dimension pla0(0:10)
      dimension psign(0:10)
      
      ii=i
      jj=j
      call dimer(nn,ll,lr,ii,jj,lqp,iqp
     &,gqp,agqp,ldx,ldx,awick)
      pla0(nnsite)=pla0(nnsite)+psign(nnsite)*awick
      ii=i
      jj=j
      call dimer(nn,ll,lr,ii,jj,lqp,iqp
     &,gqp,agqp,ldx,ldy,awick)
      pla0(nnsite)=pla0(nnsite)+psign(nnsite)*awick
      ii=i
      jj=ldy(j)
      call dimer(nn,ll,lr,ii,jj,lqp,iqp
     &,gqp,agqp,ldx,ldx,awick)
      pla0(nnsite)=pla0(nnsite)+psign(nnsite)*awick
      ii=i
      jj=ldx(j)
      call dimer(nn,ll,lr,ii,jj,lqp,iqp
     &,gqp,agqp,ldx,ldy,awick)
      pla0(nnsite)=pla0(nnsite)+psign(nnsite)*awick
      
      ii=i
      jj=j
      call dimer(nn,ll,lr,ii,jj,lqp,iqp
     &,gqp,agqp,ldy,ldx,awick)
      pla0(nnsite)=pla0(nnsite)+psign(nnsite)*awick
      ii=i
      jj=j
      call dimer(nn,ll,lr,ii,jj,lqp,iqp
     &,gqp,agqp,ldy,ldy,awick)
      pla0(nnsite)=pla0(nnsite)+psign(nnsite)*awick
      ii=i
      jj=ldy(j)
      call dimer(nn,ll,lr,ii,jj,lqp,iqp
     &,gqp,agqp,ldy,ldx,awick)
      pla0(nnsite)=pla0(nnsite)+psign(nnsite)*awick
      ii=i
      jj=ldx(j)
      call dimer(nn,ll,lr,ii,jj,lqp,iqp
     &,gqp,agqp,ldy,ldy,awick)
      pla0(nnsite)=pla0(nnsite)+psign(nnsite)*awick
      
      ii=ldy(i)
      jj=j
      call dimer(nn,ll,lr,ii,jj,lqp,iqp
     &,gqp,agqp,ldx,ldx,awick)
      pla0(nnsite)=pla0(nnsite)+psign(nnsite)*awick
      ii=ldy(i)
      jj=j
      call dimer(nn,ll,lr,ii,jj,lqp,iqp
     &,gqp,agqp,ldx,ldy,awick)
      pla0(nnsite)=pla0(nnsite)+psign(nnsite)*awick
      ii=ldy(i)
      jj=ldy(j)
      call dimer(nn,ll,lr,ii,jj,lqp,iqp
     &,gqp,agqp,ldx,ldx,awick)
      pla0(nnsite)=pla0(nnsite)+psign(nnsite)*awick
      ii=ldy(i)
      jj=ldx(j)
      call dimer(nn,ll,lr,ii,jj,lqp,iqp
     &,gqp,agqp,ldx,ldy,awick)
      pla0(nnsite)=pla0(nnsite)+psign(nnsite)*awick
      
      ii=ldx(i)
      jj=j
      call dimer(nn,ll,lr,ii,jj,lqp,iqp
     &,gqp,agqp,ldy,ldx,awick)
      pla0(nnsite)=pla0(nnsite)+psign(nnsite)*awick
      ii=ldx(i)
      jj=j
      call dimer(nn,ll,lr,ii,jj,lqp,iqp
     &,gqp,agqp,ldy,ldy,awick)
      pla0(nnsite)=pla0(nnsite)+psign(nnsite)*awick
      ii=ldx(i)
      jj=ldy(j)
      call dimer(nn,ll,lr,ii,jj,lqp,iqp
     &,gqp,agqp,ldy,ldx,awick)
      pla0(nnsite)=pla0(nnsite)+psign(nnsite)*awick
      ii=ldx(i)
      jj=ldx(j)
      call dimer(nn,ll,lr,ii,jj,lqp,iqp
     &,gqp,agqp,ldy,ldy,awick)
      pla0(nnsite)=pla0(nnsite)+psign(nnsite)*awick
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc






ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c kinetic energy matrix
      subroutine kinetic(nn,m,t1,t2,delta,hk,hku
     &,ihop1,ihop2,hnqei,ee)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      
      dimension hk(nn,nn)
      dimension hku(nn,nn)
      dimension hnqei(nn),ee(nn)
      
      dimension ihop1(nn,4),ihop2(nn,2)
      
      m2=m*2


c Matrix of kinetic energy
         hk=0.0d0
         hku=0.0d0

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


      call tred2(hk,nn,nn,hnqei,ee)
      call tqli(hnqei,ee,nn,nn,hk)
      

      do i=1,nn
         do j=1,nn
            hku(j,i)=hk(i,j)
         end do
      end do

      hk=0.0d0
      do l=1,nn
      do i=1,nn
      do j=1,nn
         hk(i,l)=hk(i,l)+hku(j,i)*dexp(hnqei(j))*hku(j,l)
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
      subroutine schmidt(nn,m,Lmax,fai,fai2,lite,a,b
     &,aeigen,beigen,cc)

      implicit real*8(a-h,o-z)
      implicit integer*4(i-n) 

      dimension fai(nn,2*m,Lmax),fai2(nn,2*m)
      dimension a(m,m),b(m,m)
      dimension aeigen(m),beigen(m),cc(m)

c inner product      
!$omp parallel do
      do jnum=1,m
      do inum=1,m
         a(inum,jnum)=0.0d0
         b(inum,jnum)=0.0d0
         do i=1,nn
            a(inum,jnum)=a(inum,jnum)
     &     +fai(i,inum,lite)*fai(i,jnum,lite)
            b(inum,jnum)=b(inum,jnum)
     &     +fai(i,inum+m,lite)*fai(i,jnum+m,lite)
         end do
      end do
      end do
      

      call tred2(a,m,m,aeigen,cc)
      call tqli(aeigen,cc,m,m,a)

      call tred2(b,m,m,beigen,cc)
      call tqli(beigen,cc,m,m,b)
      
      
      do j=1,m
         do i=1,m
            a(i,j)=a(i,j)/dsqrt(aeigen(j))
         end do
      end do
      do j=1,m
         do i=1,m
            b(i,j)=b(i,j)/dsqrt(beigen(j))
         end do
      end do
      
       
      do j=1,m*2
         do i=1,nn
            fai2(i,j)=fai(i,j,lite)
         end do
      end do


      
!$omp parallel do      
      do j=1,m
         do i=1,nn
            fai(i,j,lite)=0.0d0
            do k=1,m
               fai(i,j,lite)=fai(i,j,lite)
     &                       +fai2(i,k)*a(k,j)
               fai(i,j+m,lite)=fai(i,j+m,lite)
     &                       +fai2(i,k+m)*b(k,j)
            end do
         end do
      end do



      return
      end 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc









ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mkgreen(ll,lr,nn,m,Lmax,IPIV,WORK,LWORK
     &,bufu,bufd,a,b,fai,gup,gdown,ginup,gindown)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n) 
      character*1 UPLO 
      
      dimension a(m,m),b(m,m)
      dimension fai(nn,2*m,Lmax)
      dimension gup(nn,nn,Lmax),gdown(nn,nn,Lmax)
      dimension ginup(Lmax),gindown(Lmax)
      
      dimension IPIV(m),WORK(LWORK)

      dimension bufu(m,nn),bufd(m,nn)
      
!$omp parallel do
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
      call DGETRI(m,a,m,IPIV,WORK,LWORK,INFO)
      if(INFO.ne.0) then
         write(6,*) 'inverse_failure INFO=',INFO
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
      call DGETRI(m,b,m,IPIV,WORK,LWORK,INFO)
      if(INFO.ne.0) then
         write(6,*) 'inverse_failure down spin INFO=',INFO
         stop
      end if


c make g matrix
c fast code by Koga san
!$omp parallel do
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

!$omp parallel do
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

c late code
c      do j=1,nn
c      do i=1,nn
c         gup(i,j,ll)=0.0d0
c         gdown(i,j,ll)=0.0d0
c         do k=1,m
c         do l=1,m
c            gup(i,j,ll)=gup(i,j,ll)
c     &     +fai(i,k,lr)*a(k,l)*fai(j,l,ll)
c            gdown(i,j,ll)=gdown(i,j,ll)
c     &     +fai(i,k+m,lr)*b(k,l)*fai(j,l+m,ll)
c         end do
c         end do
c      end do
c      end do

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
      subroutine Urg(nn,m,u,t1,t2,Lmax,iLdim,lite,ijm,alpha1,h,hin
     &,ihop1,ihop2,hsub,hinsub,hgeed,hingeed,WORK2,LWORK2,iss
     &,gup,gdown,ginup,gindown,gsubu,gsubd,pirgeng,lurg
     &,gup2,gdown2,ginup2,gindown2,pirgmp,huup,hudown,deigen
     &,ipstepu,icoarse)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 UPLO
      
      dimension alpha1(-1:1)
      dimension hsub(Lmax,2),hinsub(Lmax,2)
      dimension h(Lmax,Lmax),hin(Lmax,Lmax)
      dimension hgeed(Lmax,Lmax),hingeed(Lmax,Lmax)
      dimension WORK2(LWORK2),deigen(Lmax)
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
        if(ll>lite) then
        i1=lite
        i2=ll
        end if
 
         h(i1,i2)=hsub(ll,is)
         hin(i1,i2)=hinsub(ll,is)
      end do
      
      call geed0(Lmax,iLdim,h,hin
     &,hgeed,hingeed,deigen,eigen,WORK2,LWORK2,'N',label)
     
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
      
      return
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







ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mkgreenqp(ll,nn,m,Lmax,lqp,iqp,IPIVqp,WORKqp
     &,aqp,faiqpL,faiqpR,gqp,ginqp,bufqp)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n) 
      character*1 UPLO 
      
c QP    
      dimension gqp(2*nn,2*nn),ginqp(Lmax)
      dimension aqp(2*m,2*m)
      dimension faiqpL(2*nn,2*m,Lmax)
      dimension faiqpR(2*nn,2*m,lqp)
      dimension IPIVqp(2*m),WORKqp(2*m)
      dimension bufqp(2*m,2*nn)
c end QP      
      
      
      m2=2*m
      nn2=2*nn
      
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
   
c make g matrix
c fast code by Koga san
!$omp parallel do
      do j=1,nn2
      do i=1,m2
         bufqp(i,j)=0.0d0
         do l=1,m2
            bufqp(i,j)=bufqp(i,j)+aqp(i,l)*faiqpL(j,l,ll)
         end do
      end do
      end do

!$omp parallel do
      do j=1,nn2
      do i=1,nn2
         gqp(i,j)=0.0d0
         do k=1,m2
            gqp(i,j)=gqp(i,j)
     &                  +faiqpR(i,k,iqp)*bufqp(k,j)
         end do
      end do
      end do

c late code
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



!$omp parallel do
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



!$omp parallel do private(beta)
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










ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine geed0(Lmax,iLdim,h,hin
     &,hgeed,hingeed,deigen,eigen,WORK2,LWORK2,JOBZ,label)
      implicit real*8(a-h,o-z)
      implicit integer*4(i-n)
      character*1 JOBZ 
      character*1 UPLO
      dimension h(Lmax,Lmax),hin(Lmax,Lmax)
      dimension hgeed(iLdim,iLdim),hingeed(iLdim,iLdim)
      dimension deigen(iLdim),WORK2(LWORK2)
      
      data ITYPE/1/,UPLO/'U'/
      
c ED ccccccccccccccccccccccccccccccccccccccccccc      
      do lr=1,iLdim
      do ll=1,lr
         hgeed(ll,lr)=h(ll,lr)
         hingeed(ll,lr)=hin(ll,lr)
      end do
      end do
      
      call DSYGV(ITYPE,JOBZ,UPLO,iLdim,hgeed,iLdim,hingeed,iLdim
     &,deigen,WORK2,LWORK2,INFO)
      if(INFO.gt.0.and.INFO.le.iLdim) then
         write(6,*) 'DSYGV did not converge'
         label=1
c         stop
      else if(INFO.gt.iLdim) then
         write(6,*) 'B not positive definite'
         label=1
c         stop
      end if
      eigen=deigen(1)
cccccccccccccccccccccccccccccccccccccccccccccccc      





      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc















