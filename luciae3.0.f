C**                                                                    **
C**                            LUCIAE version 3.0                          **
C**
C**                    FINISHED IN  March 1998                                                     **
C**                                                                          **
C**          AUTHORS:   SA BEN-HAO AND TAI AN

C**                                           **
C**                                                                          **
C**                                                                          **
C**    Notice:                                                               **
C**     LUCIAE version 3.0 subroutine packages must 
C**    be compiled together with Fritiof 7.02R      **
C**    ARIADNE 4.02R, JETSET 7.4R, and PYTHIA 5.5.  Be sure you have the      **
C**    proper versions of these LUND programs.                               **
C**                                                                          **





C******************************** FRINGEB  *******************************


C**************************************************************************
c This program, LUCIAE 3.0 contains Firecracker Model and rescattering
C Model which handle collective effects in heavy-ion collisions-----
C multi-gluon emission and rescattering of final particles in nuclear 
Cenvironment.

C******************************** FRINGEB  *******************************

C******************************** FRINGEB  *******************************

      SUBROUTINE FRINGEB

C........................This routine administrates one complete event

      PARAMETER (KSZJ=40000,KSZ1=30,KSZ2=300,KSZ4=150)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FCRSOUT/MCL(10),RCL(20),INEL(400),NST(2,KSZ4),ACL(2,KSZ4)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRCONT2/ICT(10),ICTT(10)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/FRCLUS1/NCL(2,KSZ4),CLM(2,KSZ4),NSL(2,KSZ4)
      COMMON/FRCLUS2/ICL(2,KSZ4),JCL(2,KSZ4),NC(2)
      common/ctllist/nctl,noinel(400),nctl0
      COMMON/TMP/NMEM,MSTU24
      COMMON/LCARSI/IARSIG
      COMMON/ROPE/itime,akapa(5),parj1,parj2,parj3,parj21
      common/sa12/psa(5),ptai(5),clorenp,clorent
      COMMON/LUDAT3/MDCY(500,3),MDME(2000,2),BRAT(2000),KFDP(2000,5)
      SAVE /FRINTN0/,/FRINTN1/,/FRINTN3/,/FRPARA1/,/FRCLUS1/,/LCARSI/
     >    ,/FRCONT2/,/LUJETS/,/ctllist/,/FRCLUS2/,/sa12/,/ROPE/
      SAVE MDCYT1,MDCYT2, MDCYT3,MDCYT4,MDCYT5,MDCYT6,MDCYT7,MDCYT8,
     >    MDCYT9,MDCYT10, MDCYT11,MDCYT12,MDCYT13,MDCYT14,MDCYT15,
     >  MDCYT16,
     >    MDCYT17,MDCYT18, MDCYT19,MDCYT20,MDCYT21,MDCYT22, MDCYT23,
     >    nkapa,avkapa1,avkapa2,avkapa3,avkapa4
        dimension KSP(2,5),PSP(2,5)
C*************************************************************************
c..A rope consists of strings
c..A cluster consists of nucleons that could later on stretch to be  strings
c..or set back to be nucleons if excitation are too small. 
C..2 in array means projectile (1) and target (2)
c..NC(2) number of clusters in projectile and target(MCL(1)MCL(2))
c***********************************************************************
        DO L=1,10
        MCL(L)=0
        ENDDO
        NCTL0=0
        itime=0
        do L=1,5
        akapa(L)=0.
        psa(L)=0
        enddo   
       NFR(1) = NFR(1) + 1
       IF(NFR(1).EQ.1)THEN      
C when KFR(25)=1 PARJ(1),PARJ(2) etc. take tuned values for pp at 200 GeV/c
c Now these values are stored in PARJ1,PARJ2. They will send back to
c PARJ(1),PARJ(2) etc. after this event is finished.
        PARJ21=PARJ(21)
        PARJ1=PARJ(1)
        PARJ2=PARJ(2)
        PARJ3=PARJ(3)    
        avkapa1=0
        avkapa2=0 
        avkapa3=0
        avkapa4=0 
        do L=1,5
        ptai(L)=0
        enddo           
        nkapa=0
        IOP(21)=0  
        DO 55 L=1,400
        noinel(L)=0
55      CONTINUE
        MDCYT1=MDCY(LUCOMP(111),1)
        MDCYT2=MDCY(LUCOMP(310),1) 
        MDCYT3=MDCY(LUCOMP(3122),1)
        MDCYT4=MDCY(LUCOMP(-3122),1)
        MDCYT5=MDCY(LUCOMP(3212),1)
        MDCYT6=MDCY(LUCOMP(-3212),1)
        MDCYT7=MDCY(LUCOMP(3112),1)
        MDCYT8=MDCY(LUCOMP(-3112),1)
        MDCYT9=MDCY(LUCOMP(3222),1)
        MDCYT10=MDCY(LUCOMP(-3222),1)
        MDCYT11=MDCY(LUCOMP(3312),1)
        MDCYT12=MDCY(LUCOMP(-3312),1)
        MDCYT13=MDCY(LUCOMP(3322),1)
        MDCYT14=MDCY(LUCOMP(-3322),1)
        MDCYT15=MDCY(LUCOMP(3334),1)
        MDCYT16=MDCY(LUCOMP(-3334),1)
        MDCYT17=MDCY(LUCOMP(1114),1)
        MDCYT18=MDCY(LUCOMP(2114),1)
        MDCYT19=MDCY(LUCOMP(2214),1)
        MDCYT20=MDCY(LUCOMP(2224),1)
        MDCYT21=MDCY(LUCOMP(213),1)
        MDCYT22=MDCY(LUCOMP(-213),1)
        MDCYT23=MDCY(LUCOMP(113),1)
        ENDIF

        MDCY(LUCOMP(111),1) = 0
        MDCY(LUCOMP(310),1) = 0
        MDCY(LUCOMP(3122),1) = 0
        MDCY(LUCOMP(-3122),1) = 0
        MDCY(LUCOMP(3212),1) = 0
        MDCY(LUCOMP(-3212),1) = 0
        MDCY(LUCOMP(3112),1) = 0
        MDCY(LUCOMP(-3112),1) = 0
        MDCY(LUCOMP(3222),1) = 0
        MDCY(LUCOMP(-3222),1) = 0
        MDCY(LUCOMP(3312),1) = 0
        MDCY(LUCOMP(-3312),1) = 0
        MDCY(LUCOMP(3322),1) = 0
        MDCY(LUCOMP(-3322),1) = 0
        MDCY(LUCOMP(3334),1) = 0
        MDCY(LUCOMP(-3334),1) = 0
c       MDCY(LUCOMP(1114),1)=0
c       MDCY(LUCOMP(2114),1)=0
c       MDCY(LUCOMP(2214),1)=0
c       MDCY(LUCOMP(2224),1)=0
c       MDCY(LUCOMP(213),1)=0
c       MDCY(LUCOMP(-213),1)=0
c       MDCY(LUCOMP(113),1)=0
c do not allow  lambda,lambdaba,K_S0 etc. to decay before rescattering,
c but Delta and rho are decayed because of their short life-times.
2       IOP(16)=0     
        IOP(19)=0
        IARSIG=0
c-All binary collisions will be proceeded with four-momentum stored 
c-  in COMMON/FRINTN1/ and RPS information in COMMON/FRINTN2  
       CALL FRCOLLS

C-calculating incoming energy and charge for checking four-momentum conservation
c-in the end
        IF(KFR(23).EQ.1)CALL FRARCHE(0)

          NMEM=0
          MSTU24=0      
        IF(KFR(15).EQ.0.OR.IOP(2).LT.2)THEN
        N=0
        GOTO 400
        ELSE
        ENDIF
C..when KFR(15)=0 or there is only one collision,then no firecracker effect
        CALL FRCLUEX(0)
C..store information of wounded nucleons.Later on event record is divided
c..into one without RPS(soft part) and another with RPS (hard part)
c..nucleons without hard PRS gluons go through radiation
        CALL FRCLUST 
        CALL FRCRGLU
C-sort out nucleons with hard collision 
        CALL FRSORSH
c..nucleons with hard PRS gluons go through radiation
        IOP(24)=0
400     CALL FRINGEH
        IF(IOP(24).EQ.1)GOTO 2

C..store useful information in IOP()

        IF(IOP(2).GE.2)THEN
        IF(KFR(15).EQ.1)THEN
        MCL(1)=NC(1)
        MCL(2)=NC(2)
        DO 56 L=1,2
        DO 56 L1=1,MCL(L)
        NST(L,L1)=NSL(L,L1)
        ACL(L,L1)=CLM(L,L1)
56      CONTINUE
C..restore original information about wounded nucleons
        CALL FRCLUEX(1)    
        ELSE
        MCL(1)=0
        MCL(2)=0
        DO 58 L=1,2
        NST(L,1)=0
        ACL(L,1)=0.
58      CONTINUE        
        ENDIF
        ELSE
        MCL(1)=1
        MCL(2)=1
        DO 59 L=1,2
        NST(L,1)=1
        ACL(L,1)=1.
59      CONTINUE
        ENDIF
        MCL(3)=IOP(19)  
C....TO ADD ONTO LUJETS THE COLOUR NEUTRAL PARTICLES THAT MAY 
C....HAVE BEEN PRODUCED FROM PARTON-PARTON PROCESSES: 
      CALL FRFILHW
      IF(IARSIG.EQ.1)then
c       WRITE(MSTU(11),*)'Warning:An error might be found in Ariande'
c       CALL LULIST(2)
         GOTO 2
        endif
      IF(N.GE.KSZJ-2) CALL FRMGOUT(0,1,
     > 'LUJETS array size KSZJ must be expanded',float(N),float(KSZJ),
     >  0.,0.,0.)
        DO 670 L=1,N
        IF(K(L,2).EQ.21)P(L,5)=0.
c       WRITE(MSTU(11),*)'Warning:A error due to precision is found, the event
c     & was thrown away'
        IF(ABS(P(L,1)).GT.10000.OR.ABS(P(L,2)).GT.10000.OR.
     &  P(L,4).LT.0.0)GOTO 2
670     CONTINUE        
      IF(KFR(1).EQ.1) THEN
  
        MSTJ21=MSTJ(21)
        MSTJ(21)=0
      CALL  LUEXEC
       if(KFR(25).eq.1)then
        PARJ(21)=PARJ21
        PARJ(1)=PARJ1
        PARJ(2)=PARJ2
        PARJ(3)=PARJ3  
        endif     
      IF(MSTU(24).EQ.4) MSTU24=1
        IF(KFR(19).EQ.0.OR.KFR(21).EQ.1)     CALL LUEDIT(1)
        MSTJ(21)=MSTJ21
      CALL  LUEXEC

        IF(KFR(19).EQ.0.OR.KFR(21).EQ.1)     CALL LUEDIT(1)
        IF(KFR(13).GE.4) THEN 
      CALL FREDITD()
        ENDIF
      ENDIF
C... Regenerate event in case of 'infinite-loop error' in Jetset:
      IF(MSTU24.GT.0)then
        mstu(23)=0
         GOTO 2
        endif
Ccheck particle list before rescattering
        DO 660 L=1,N
        IF(ABS(P(L,1)).GT.10000.OR.ABS(P(L,2)).GT.10000.OR.
     &  P(L,4).LT.0.0)THEN
c       WRITE(MSTU(11),*)'Warning:A error due to precision is found, the event
c     & was thrown away','EVENT=',NFR(1)
        GOTO 2
        ENDIF
660     CONTINUE
c       itime is the number of strings in the current event     
        akapa(1)=akapa(1)/max(1,itime)
        akapa(2)=akapa(2)/max(1,itime)
        akapa(3)=akapa(3)/max(1,itime)
        akapa(4)=akapa(4)/max(1,itime)
        akapa(5)=akapa(5)/max(1,itime)
        IF(itime.ne.0)nkapa=nkapa+1
        avkapa1=avkapa1+akapa(4)
        avkapa2=avkapa2+akapa(2)
        avkapa3=avkapa3+akapa(5)
        avkapa4=avkapa4+akapa(3)
C....Record the number of N-N collisions:
      IOP(2) = IOP(2)-ICT(3)-ICT(8)-ICT(10)
      NFR(3) = NFR(3) + IOP(2)

       NFR(4) = NFR(4) + IOP(13)
       IF(IOP(13).GE.1) NFR(5) = NFR(5)+1
      
C-----RECORDING OF IMPACT PARAMETER AND COUNTING OF SPECTATOR PROTONS  

      IOP(11)=IOP(4)
      IOP(12)=IOP(6)
      DO 200 L=1, 2
      DO 200 J=1,IOP(8+L)
      IF(IDN(L,J).EQ.2212) IOP(10+L)=IOP(10+L)-1
  200 CONTINUE

      IF(IDN(1,1).NE.2212.AND.IDN(1,1).NE.2112) IOP(11)=0 

C.....Add the nuclei spectators onto the event record:
      DO 300 L=1,2
      IF(IOP(3+2*(L-1))-IOP(8+L).GE.1) THEN
      N = N+1
      P(N,1)=PPA(L,1)
      P(N,2)=PPA(L,2)
      P(N,3)=0.5*(PPA(L,4)-PPA(L,3))
      P(N,4)=0.5*(PPA(L,4)+PPA(L,3))
      P(N,5)=PPA(L,5)
      K(N,1) = 1
      K(N,2) = (10000+IOP(10+L))*(-1)**(L-1)
      ENDIF
300   CONTINUE
C....Dump out the event for inspection when error occurs in FRMGOUT:
      IF(IOP(16).GE.1) CALL LULIST(2)

      IF(KFR(14).EQ.1) CALL FRCHKEP(1)
C proceed rescattering  of final particles
        IF(KFR(21).EQ.1.AND.KFR(1).EQ.1)CALL FRRESEXE
        IF(NFR(1).EQ.1)CALL FRVALUE(0)          
        MDCY(LUCOMP(111),1)=MDCYT1
        MDCY(LUCOMP(310),1) =MDCYT2
        MDCY(LUCOMP(3122),1)=MDCYT3
        MDCY(LUCOMP(-3122),1) =MDCYT4
        MDCY(LUCOMP(3212),1)=MDCYT5
        MDCY(LUCOMP(-3212),1) =MDCYT6
        MDCY(LUCOMP(3112),1)=MDCYT7
        MDCY(LUCOMP(-3112),1) =MDCYT8
        MDCY(LUCOMP(3222),1)=MDCYT9
        MDCY(LUCOMP(-3222),1) =MDCYT10
        MDCY(LUCOMP(3312),1) =MDCYT11
        MDCY(LUCOMP(-3312),1) =MDCYT12
        MDCY(LUCOMP(3322),1) =MDCYT13
        MDCY(LUCOMP(-3322),1) =MDCYT14
        MDCY(LUCOMP(3334),1) =MDCYT15
        MDCY(LUCOMP(-3334),1) =MDCYT16
        MDCY(LUCOMP(1114),1)=MDCYT17
        MDCY(LUCOMP(2114),1)=MDCYT18
        MDCY(LUCOMP(2214),1)=MDCYT19
        MDCY(LUCOMP(2224),1)=MDCYT20
        MDCY(LUCOMP(213),1)=MDCYT21
        MDCY(LUCOMP(-213),1)=MDCYT22
        MDCY(LUCOMP(113),1)=MDCYT23
c output some physics quantities which may be interesting
        RCL(1)=avkapa1/max(1,nkapa)
        RCL(2)=avkapa2/max(1,nkapa)
        RCL(3)=avkapa3/max(1,nkapa)
        RCL(4)=avkapa4/max(1,nkapa)             
        MCL(4)=NCTL0
        RCL(5)=ptai(1)
        RCL(6)=ptai(2)
        RCL(7)=ptai(3)
        RCL(8)=ptai(4)
        DO 552 L=1,400
        INEL(L)=noinel(L)
552     CONTINUE
cstore spectators
        LSP=0
        DO 570 L=N-1,N
        IF(ABS(K(L,2)).GE.10000)THEN
        LSP=LSP+1
        DO 550 L1=1,5
        KSP(LSP,L1)=K(L,L1)
        PSP(LSP,L1)=P(L,L1)
550     CONTINUE
        ENDIF
570     CONTINUE
        N=N-LSP
        IF(KFR(1).EQ.1) THEN
        CALL  LUEXEC
        IF(MSTU(24).EQ.4)THEN
        WRITE(MSTU(11),*)'infinite loop occur at event=',NFR(1)
        call lulist(1)
        endif
c       CALL  LUEDIT(1)
        ENDIF
C put spectators back
        LSP0=0
        DO 560 L=N+1,N+LSP
        LSP0=LSP0+1
        DO 560 L1=1,5
        K(L,L1)=KSP(LSP0,L1)
        P(L,L1)=PSP(LSP0,L1)
560     CONTINUE
        N=N+LSP0
c--check four-momentum and charge conservation
        IF(KFR(23).EQ.1)CALL FRARCHE(1)
Clast check
        DO 640 L=1,N
        IF(ABS(P(L,1)).GT.10000.OR.ABS(P(L,2)).GT.10000.OR.
     &  P(L,4).LT.0.0)THEN
c       WRITE(MSTU(11),*)'Warning:A error due to precision is found, the event
c     & was thrown away'
        GOTO 2
        ENDIF
640     CONTINUE
      RETURN
      END 

C******************************** END FRINGEB  ***************************
        subroutine strtension(IP,np)
      PARAMETER (KSZJ=40000,KSZ1=30)
      COMMON/ROPE/itime,akapa(5),parj1,parj2,parj3,parj21
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      SAVE /ROPE/
        VFR26=1.15      
c calculate the effective string tension of each string and then
c calculate the new PARJ(1),PARJ(2) etc. in the current event.
C  VFR26=1.15 the effective string tension calculated by FRITIOF 
c       for the pp at 200 GeV/c 
c when this subroutine is called in JETSET the string has been boosted
c to its own CMS frame. The whole string takes up item n+1 to n+np
c in the particle list.
      toteng=0.0
      toten=0.0
      totglukt=0.0
      pmax=0.
        do i=n+1,n+np
        toten=toten+p(i,4)
        pp=sqrt(p(i,1)**2+p(i,2)**2)
        if(k(i,2).eq.21.and.pp.gt.VFR(25))
     &  toteng=toteng+log(pp/VFR(25))   
        if(k(i,2).eq.21.and.pp.gt.pmax)pmax=pp
        enddo
        if(pmax.gt.VFR(25))totglukt=totglukt+log(pmax/VFR(25))  
        ss=log(toten/VFR(25))+toteng
        VFR24=(1.-(totglukt/ss))**(-VFR(24))    
        akapa(1)=akapa(1)+VFR24
        PARJ(21)=PARJ21*((VFR24/VFR26)**(0.5))
        PARJ(1)=PARJ1**(VFR26/VFR24)
        PARJ(2)=PARJ2**(VFR26/VFR24)
        PARJ(3)=PARJ3**(VFR26/VFR24)
        akapa(2)=akapa(2)+parj(2)
        akapa(3)=akapa(3)+parj(21)
        akapa(4)=akapa(4)+parj(1)
        akapa(5)=akapa(5)+parj(3)
        itime=itime+1
        return
        end
C*********************** SUBROUTINE FRARIAD ***********************

      SUBROUTINE FRARROP

C..Fritiof interface to Ariadne_4.02r.  LUJETS entries from IOP(17) to N
C..are copied to Ariadne event record ARJETX, and after emission is done
C..partons are copied back onto LUJETS.


      PARAMETER (KSZJ=40000,KSZ1=30,KSZ2=300)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      COMMON/ARJETX/NO,KO(1000,5),PO(1000,5),VO(1000,5)

      SAVE /FRINTN0/,/LUJETS/,/ARJETX/


      DO 100 I=IOP(17),N
        
      NO=NO+1

      DO 100 LO=1,5



      KO(NO,LO) = K(I,LO)
100   PO(NO,LO) = P(I,LO)

        RETURN
        END


C*********************** END FRARIAD *****************************
      SUBROUTINE FRCLUST
      PARAMETER (KSZJ=40000,KSZ1=30, KSZ2=300,KSZ4=150,PI=3.1415926)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)
      COMMON/FRCLUS1/NCL(2,KSZ4),CLM(2,KSZ4),NSL(2,KSZ4)
      COMMON/FRCLUS2/ICL(2,KSZ4),JCL(2,KSZ4),NC(2)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRTEMP1/NOICF(2,KSZ2,100)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
        DOUBLE PRECISION PTS(2,KSZ2,5),P1,P2,PH1,PH2,
     &  Y1,Y2,PTST(2,5),PTSY(2,KSZ2,5)
        INTEGER NN(2,KSZ2),NNX(2,KSZ2),WN(2),NOIC(2,KSZ2,100),
     &  NOICT(100),NOICTA(KSZ2,100)
        DIMENSION PTCUT(2)
      SAVE /FRINTN0/,/FRINTN1/,/FRPARA1/,/FRCLUS1/,/FRCLUS2/,/FRTEMP1/
C*******************************************************************************
c..NN(2,KSZ2) number of nucleons in KSZ2_th entry
c..NOIC(2,KSZ2,J) original entry number of j_th nucleon in KSZ2_th entry
C********************************************************************************
c.. calculating average value of Pt and take VFR(17)*sqrt(<Pt>^2) to be size of cluster 
c::in A-A,in p-A,only one cluster is formed.Soft energy is used to pick up cluster
c.. collisions when KFR(18)=1.Otherwise size of cluster = VFR(18)
        DO 31 L=1,2
        DO 31 L2=1,IOP(8+L)     
        NN(L,L2)=1
        NOIC(L,L2,NN(L,L2))=L2  
        DO 31 L3=1,5
        PTS(L,L2,L3)=PPS(L,L2,L3)
31      CONTINUE
        DO 29 L=1,2
        APT=0
        IF(IOP(3).EQ.1)THEN
        PTCUT(L)=1.E+10
        ELSE
        IF(IOP(8+L).GT.0)THEN
        IF(KFR(18).EQ.1)THEN
        DO 30 L2=1,IOP(8+L)
        APT=SQRT(PTS(L,L2,1)**2+PTS(L,L2,2)**2)+APT
30      CONTINUE
        PTCUT(L)=VFR(17)*APT/IOP(8+L)
        ELSE
        PTCUT(L)=VFR(18)
        ENDIF
        ELSE
        PTCUT(L)=0
        ENDIF
        ENDIF
        NC(L)=0
        WN(L)=IOP(8+L)
29      CONTINUE

        DO 33 L=1,2
        IF(IOP(8+L).EQ.1)THEN
        NC(L)=IOP(8+L)
        NCL(L,1)=IOP(8+L)
        NOICF(L,1,1)=IOP(8+L)
        GOTO 33
        ELSEIF(IOP(8+L).EQ.0)THEN
        NC(L)=IOP(8+L)
        NCL(L,1)=IOP(8+L)
C..matrix index >0
        NOICF(L,1,1)=1
        GOTO 33
        ELSE
        ENDIF
c..cluster nucleons in which 'distance' between nucleons is smaller than ptcut
32      PTMIN=1.E+10
        DO 34 L2=1,IOP(8+L)-1
        PT1=DSQRT(PTS(L,L2,2)**2+PTS(L,L2,1)**2)
        PH1=0.5*(PTS(L,L2,4)-PTS(L,L2,3))
        P1=DSQRT(PT1**2+PH1**2)
        DO 35 L3=L2+1,IOP(8+L)
        PT2=DSQRT(PTS(L,L3,2)**2+PTS(L,L3,1)**2)
        PH2=0.5*(PTS(L,L3,4)-PTS(L,L3,3))
        P2=DSQRT(PT2**2+PH2**2)
        Y1=0.5*DLOG(DMAX1(1.D-10,(P1+PH1))/DMAX1(1.D-10,(P1-PH1)))
        Y2=0.5*DLOG(DMAX1(1.D-10,(P2+PH2))/DMAX1(1.D-10,(P2-PH2)))
        FI1=ULANGL(REAL(PTS(L,L2,1)),REAL(PTS(L,L2,2)))
        IF(FI1.LT.0.)FI1=2*PI+FI1
        FI2=ULANGL(REAL(PTS(L,L3,1)),REAL(PTS(L,L3,2)))
        IF(FI2.LT.0.)FI2=2*PI+FI2
        FI12=ABS(FI1-FI2)       
        RD=SQRT(REAL((Y1-Y2))**2+(MIN(FI12,2*PI-FI12))**2)      
        PT=MIN(PT1,PT2)*RD
        IF(PT.LT.PTMIN)THEN
        PTMIN=PT
        KT1=L2
        KT2=L3
        ENDIF   
35      CONTINUE
34      CONTINUE
        IF(PTMIN.GT.PTCUT(L))THEN
        DO 387 JI=1,IOP(8+L)
        NC(L)=NC(L)+1
        IF(NC(L).GT.150) THEN
        WRITE(MSTU(11),*)'number of clusters > 150,stop runing'
        STOP
        ENDIF
        NCL(L,NC(L))=NN(L,JI)
        IF(NCL(L,NC(L)).GT.100) THEN
        WRITE(MSTU(11),*)'number of nucleons in a cluster > 
     &  100,stop runing'
        STOP
        ENDIF
        DO 72 I1=1,NN(L,JI)
72      NOICF(L,NC(L),I1)=NOIC(L,JI,I1)
387     CONTINUE
        GOTO 33
        ELSE
c..update event record for next grouping,here one entry could 
C..be a pseduparticle consisting of many nucleons.The newest 
C..pseduparticle is always put into the first entry
        DO 36 L4=1,4
        PTST(L,L4)=PTS(L,KT1,L4)+PTS(L,KT2,L4)
36      CONTINUE
        I2=1
        NNT=NN(L,KT1)+NN(L,KT2)
        DO 42 I1=1,NN(L,KT1)
42      NOICT(I1)=NOIC(L,KT1,I1)
        DO 43 I1=1,NN(L,KT2)
43      NOICT(I1+NN(L,KT1))=NOIC(L,KT2,I1)
        ENDIF
        DO 37 I1=1,IOP(8+L)
        IF(I1.EQ.KT1.OR.I1.EQ.KT2)GOTO 37
        I2=I2+1
        NNX(L,I2)=NN(L,I1)
        DO 73 I3=1,NN(L,I1)
73      NOICTA(I2,I3)=NOIC(L,I1,I3)
        DO 38 L4=1,4
        PTSY(L,I2,L4)=PTS(L,I1,L4)
38      CONTINUE
37      CONTINUE
        DO 57 L4=1,4
        PTS(L,1,L4)=PTST(L,L4)
57      CONTINUE
        NN(L,1)=NNT
        DO 75 I3=1,NN(L,1)
75      NOIC(L,1,I3)=NOICT(I3)
        IOP(8+L)=IOP(8+L)-1
        DO 60 I1=2,IOP(8+L)
        DO 59 L4=1,4
        PTS(L,I1,L4)=PTSY(L,I1,L4)
59      CONTINUE
        NN(L,I1)=NNX(L,I1)
        DO 77 I3=1,NN(L,I1)
77      NOIC(L,I1,I3)=NOICTA(I1,I3)
60      CONTINUE
        IF(IOP(8+L).GT.1)THEN
        GOTO 32
        ELSEIF(IOP(8+L).EQ.1)THEN
        NC(L)=NC(L)+1
        NCL(L,NC(L))=NN(L,1)
        DO 78 I1=1,NN(L,1)
78      NOICF(L,NC(L),I1)=NOIC(L,1,I1)
        ELSE
        ENDIF   
33      CONTINUE
        IOP(9)=WN(1)
        IOP(10)=WN(2)
c..nucleons are ordered group by group in event record
        CALL FRARORD
        RETURN
        END
C ************************************************************
          SUBROUTINE  FRARORD

      PARAMETER (KSZJ=40000,KSZ1=30, KSZ2=300,KSZ4=150)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)
      COMMON/FRCLUS1/NCL(2,KSZ4),CLM(2,KSZ4),NSL(2,KSZ4)
      COMMON/FRCLUS2/ICL(2,KSZ4),JCL(2,KSZ4),NC(2)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRTEMP1/NOICF(2,KSZ2,100)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      COMMON/FRINTN4/KFEND(2,KSZ2,2)
      SAVE /FRINTN0/,/FRINTN1/,/FRCLUS1/,/FRCLUS2/,/FRINTN4/,/FRINTN3/,
     &   /FRTEMP1/
      DIMENSION JC(2),PTSY(2,KSZ2,5),PTSX(2,KSZ2,5),PTS(2,KSZ2,5)
        DIMENSION IDNT(2,KSZ2),FMNT(2,KSZ2),KFENDT(2,KSZ2,2)
c..purpose: event record is reproduced with nucleons in the same group 
c..being placed together
c..ICL(2,KSZ4)entry number of the first nucleon in KSZ4_th group
c..JCL(2,KSZ4)entry number of the last  nucleon in KSZ4_th group
c..ICL(2,KSZ4) and JCL(2,KSZ4) correspond to PPS(2,KSZ2,5),PPH(2,KSZ2,5) etc
c..which index KSZ2 starts from 1 for both target and projectile
        JC(1)=0
        JC(2)=0
        DO 301 L=1,2
        ICL(L,1)=1
        JCL(L,1)=NCL(L,1)
        IF(NC(L).EQ.0)GOTO 301
        DO 300 I1=1,NC(L)
        IF(I1.GT.1)THEN
        ICL(L,I1)=JCL(L,I1-1)+1
        JCL(L,I1)=ICL(L,I1)+NCL(L,I1)-1
        ENDIF
        DO 300 I2=1,NCL(L,I1)
        JC(L)=JC(L)+1
        IF(NOICF(L,I1,I2).EQ.0)then
        WRITE(MSTU(11),*)'an error may have taken place 
     &  in subroutine FRARORD ,some information is outputed'
        DO 10 J1=1,2
        WRITE(MSTU(11),*)'NC(L)=',NC(J1),IOP(8+J1)
        DO 10 J2=1,NC(J1)
        WRITE(MSTU(11),*)'NCL(J1,J2)=',NCL(J1,J2)
        DO 10 J3=1,NCL(J1,J2)
        WRITE(MSTU(11),*)'NOICF(J1,J2,J3)=',NOICF(J1,J2,J3)
10      CONTINUE
        ENDIF
        DO 310 I3=1,5   
        PTSX(L,JC(L),I3)=PPS(L,NOICF(L,I1,I2),I3)
        PTSY(L,JC(L),I3)=PPH(L,NOICF(L,I1,I2),I3)
        PTS(L,JC(L),I3)=PPSY(L,NOICF(L,I1,I2),I3)
310     CONTINUE

        IDNT(L,JC(L))=IDN(L,NOICF(L,I1,I2))
        FMNT(L,JC(L))=FMN(L,NOICF(L,I1,I2))

        DO 234 J1=1,2
        KFENDT(L,JC(L),J1)=KFEND(L,NOICF(L,I1,I2),J1)
234     CONTINUE
300     CONTINUE
301     CONTINUE
        DO 330 L=1,2
        DO 330 I1=1,IOP(8+L)
        DO 340 I3=1,5   
        PPS(L,I1,I3)=PTSX(L,I1,I3)
        PPH(L,I1,I3)=PTSY(L,I1,I3)
        PPSY(L,I1,I3)=PTS(L,I1,I3)
340     CONTINUE
        IDN(L,I1)=IDNT(L,I1)
        FMN(L,I1)=FMNT(L,I1)
        DO 231 J1=1,2
231     KFEND(L,I1,J1)=KFENDT(L,I1,J1)
330     CONTINUE
        CALL FRARORDT
        RETURN     
        END
C***********************************************************************
C************************************************************************************
        SUBROUTINE FRARCHE(IG)
         PARAMETER (KSZJ=40000,KSZ1=30)
        COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
        COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa12/psa(5),ptai(5),clorenp,clorent
        SAVE /sa12/
        COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
         COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
          COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
        DIMENSION PS(2,6)
        SAVE /LUJETS/,/FRINTN0/,/LUDAT2/
        SAVE PS
C--check four-momentum and charge conservation
        IF(IG.EQ.0)THEN
C-...Sum up total momentum, energy and charge for collision system at 
c    beginning
      DO 50 I=1,2
      DO 50 J=1,6
50    PS(I,J)=0.
      DO 60 I=1,2
      PS(1,1)=PS(1,1)+PLI0(I,1)*IOP(2*I+1)
      PS(1,2)=PS(1,2)+PLI0(I,2)*IOP(2*I+1)
      PS(1,3)=PS(1,3)+0.5*(PLI0(I,4)-PLI0(I,3))*IOP(2*I+1)
      PS(1,4)=PS(1,4)+0.5*(PLI0(I,4)+PLI0(I,3))*IOP(2*I+1)
      PS(1,6)=PS(1,6)+IOP(4+2*(I-1))
60      CONTINUE
        ELSE
        NCO=0
C...Check that momentum, energy and charge were conserved.
      DO 58 I=1,N
      IF(K(I,1).LE.0.OR.K(I,1).GT.10) GOTO 58
        KC=LUCOMP(K(I,2))
        IF(KC.EQ.0)GOTO 57
        KQ=KCHG(KC,2)*ISIGN(1,K(I,2))
C..check if there is coloured particles on decayed strings
        IF(KQ.NE.0)NCO=NCO+1    
57      DO 56 J=1,4
56    PS(2,J)=PS(2,J)+P(I,J)
        IF(ABS(K(I,2)).GE.10000)GOTO 58
      PS(2,6)=PS(2,6)+PLU(I,6)
58    CONTINUE
C..add charges of spectators
c..accuracy parameter used in Jetset for checking four-momentum conservation
c is increased by a factor of 5

        PS(2,6)=PS(2,6)+IOP(11)+IOP(12)+psa(5)
        
      PDEV=(ABS(PS(2,1)-PS(1,1))+ABS(PS(2,2)-PS(1,2))+ABS(PS(2,3)-
     &PS(1,3))+ABS(PS(2,4)-PS(1,4)))/(1.+ABS(PS(2,4))+ABS(PS(1,4)))
      IF(PDEV.GT.5.0*PARU(11))THEN
         IOP(21)=IOP(21)+1
        IF(IOP(21).EQ.MSTU(22)+1)THEN
        WRITE(MSTU(11),650) MSTU(22)
        ELSEIF(IOP(21).LE.MSTU(22))THEN
                CALL LULIST(1)
         WRITE(MSTU(11),*)'EVENT=',NFR(1)
         WRITE(MSTU(11),*)'four-momentum was not conserved'
         WRITE(MSTU(11),*)'initial four-momentum',ps(1,1),ps(1,2),
     &          ps(1,3),ps(1,4)
         WRITE(MSTU(11),*)'final four-momentum',ps(2,1),ps(2,2),
     &          ps(2,3),ps(2,4)
         WRITE(MSTU(11),*)'annihilation four-momentum=',psa(1),
     &  psa(2),psa(3),psa(4),psa(5)

        ELSE
        ENDIF
        ENDIF
        IF(NCO.NE.0.AND.KFR(1).EQ.1)THEN
         WRITE(MSTU(11),*)'EVENT=',NFR(1),'NCP=',NCO
         WRITE(MSTU(11),*)'colour particles found on decayed 
     &        strings error is serious,list last event and stop runing'
        CALL LULIST(1)
        STOP
        ENDIF
       IF(ABS(PS(2,6)-PS(1,6)).GT.0.1)THEN
         IOP(21)=IOP(21)+1
        IF(IOP(21).EQ.MSTU(22)+1)THEN
        WRITE(MSTU(11),650) MSTU(22)
        ELSEIF(IOP(21).LE.MSTU(22))THEN
                CALL LULIST(1)
         WRITE(MSTU(11),*)'EVENT=',NFR(1)
                WRITE(MSTU(11),*)'charge was not conserved'
         WRITE(MSTU(11),*)'initial net charge',PS(1,6)
         WRITE(MSTU(11),*)'final net charge',PS(2,6)
         WRITE(MSTU(11),*)'annihilation charge=',psa(1),
     &  psa(2),psa(3),psa(4),psa(5)

        ELSE
        ENDIF
        ENDIF
650     FORMAT(/,2X,'charge or four-momentum conservation
     &  was violated',I3,'times,no warning any more hereafter',/)
        ENDIF
      RETURN
      END       

C***********************************************************************
      SUBROUTINE FRSORSH
      PARAMETER (KSZ1=30, KSZ2=300)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)
      COMMON/FRINTN2/NHP(2),IHQP(2,KSZ2),KHP(2,KSZ2,100,5),
     >   PHP(2,KSZ2,100,5)   
      COMMON/FRINTN4/KFEND(2,KSZ2,2)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)    
        SAVE /FRINTN1/,/FRINTN2/,/FRINTN3/,/FRINTN4/,/FRINTN0/

        DIMENSION NHN(2)
c-purpose:sorting out nucleons with RPS 

        DO 340 L=1,2
        NHN(L)=0
        DO 330 I1=1,IOP(8+L)
        IF(IHQP(L,I1).GT.0)THEN
        NHN(L)=NHN(L)+1
        DO 310 I3=1,5   
        PPS(L,NHN(L),I3)=PPS(L,I1,I3)
        PPH(L,NHN(L),I3)=PPH(L,I1,I3)
                PPSY(L,NHN(L),I3)=PPSY(L,I1,I3)
310     CONTINUE
        IHQP(L,NHN(L))=IHQP(L,I1)
        IDN(L,NHN(L))=IDN(L,I1)
        FMN(L,NHN(L))=FMN(L,I1)
        DO 235 J1=1,2
235     KFEND(L,NHN(L),J1)=KFEND(L,I1,J1)       
        DO 350 I4=1,IHQP(L,I1)
        DO 350 I3=1,5           
        KHP(L,NHN(L),I4,I3)=KHP(L,I1,I4,I3)
        PHP(L,NHN(L),I4,I3)=PHP(L,I1,I4,I3)
350     CONTINUE
        ENDIF
330     CONTINUE
        IOP(8+L)=NHN(L)
340     CONTINUE
        RETURN
        END
C******************************************************************

      SUBROUTINE FRCOLLS

C........................This routine administrates collission of hadrons

      PARAMETER (KSZ1=30,KSZ2=300)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      COMMON/FRINTN4/KFEND(2,KSZ2,2)
      SAVE /FRINTN0/,/FRINTN3/,/FRINTN4/



C.....Randomly order protons and neutrons:
      CALL FRHILDN

C.....Create the nuclei and calculate number of collisions
      CALL FRANGAN

C.....Fix the end flavors for all the strings
      DO 5 L=1,2 
      DO 5 II=1, IOP(3+2*(L-1))
   5  CALL FRBELEO(KFEND(L,II,1),KFEND(L,II,2),IDN(L,II))

C.....Generate the masses and momenta after the collisions
C.....Generate the masses and momenta after the collisions

 10   CALL FRRINGO
        RETURN
        END
C******************************************************************************
        SUBROUTINE FRARCLU(IOO)

      PARAMETER (KSZJ=40000)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      COMMON/ARJETX/NO,KO(1000,5),PO(1000,5),VO(1000,5)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUJETS/,/ARJETX/



c.purpose interface between Fritiof and Ariande

        CALL AREXEC 
c..put colour singlet particle onto event record
        IOO=0
        DO 509 IO=1,N
                KC=LUCOMP(K(IO,2))
        KQ=KCHG(KC,2)*ISIGN(1,K(IO,2))
        IF(KQ.NE.0)GOTO 509
        IOO=IOO+1
       DO 405 LO=1,5
       KO(NO+IOO,LO) = K(IO,LO)
405    PO(NO+IOO,LO) = P(IO,LO)
509    CONTINUE
        NO=NO+IOO
        N=0
      DO 400 IO=1,NO
      IF(KO(IO,1).GE.11) GOTO 400
      N=N+1
      DO 450 LO=1,5
      K(N,LO) = KO(IO,LO)
450   P(N,LO) = PO(IO,LO)
400   CONTINUE


        RETURN
        END

C************************************************************************

      SUBROUTINE FRCLUEX(I)
      PARAMETER (KSZ2=300,KSZ1=30)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)
      COMMON/FRINTN2/NHP(2),IHQP(2,KSZ2),KHP(2,KSZ2,100,5),
     & PHP(2,KSZ2,100,5)
      COMMON/FRINTN3/IDN(2,KSZ2),FMN(2,KSZ2),NUC(2,3000)
      COMMON/FRINTN3T/IDNT(2,KSZ2),FMNT(2,KSZ2)
      COMMON/FRINTN4/KFEND(2,KSZ2,2)
      COMMON/FRINTN4T/KFENDT(2,KSZ2,2)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRINTN1T/PPST(2,KSZ2,5),PPHT(2,KSZ2,5),PPSYT(2,KSZ2,5)
      COMMON/FRINTN2T/IHQPT(2,KSZ2),KHPT(2,KSZ2,40,5),
     & PHPT(2,KSZ2,40,5)
        DIMENSION IOPT(2)
        SAVE /FRINTN1T/,/FRINTN2T/,/FRINTN3T/,/FRINTN4T/,
     &  /FRINTN1/,/FRINTN3/,/FRINTN2/,/FRINTN4/
        SAVE  IOPT      
C-I=0 store original event record into temperary one
c-I=1 restore original event
        IF(I.EQ.0)THEN
        DO 100 L1=1,2
        IOPT(L1)=IOP(8+L1)      
        DO 100 L2=1,IOP(8+L1)
        IDNT(L1,L2)=IDN(L1,L2)
        FMNT(L1,L2)=FMN(L1,L2)
        IHQPT(L1,L2)=IHQP(L1,L2)
        DO 232 J1=1,2
        KFENDT(L1,L2,J1)=KFEND(L1,L2,J1)
232     CONTINUE
        DO 100 L4=1,IHQP(L1,L2)
        DO 100 L3=1,5
        PPST(L1,L2,L3)=PPS(L1,L2,L3)
        PPHT(L1,L2,L3)=PPH(L1,L2,L3)
        PPSYT(L1,L2,L3)=PPSY(L1,L2,L3)
        KHPT(L1,L2,L4,L3)=KHP(L1,L2,L4,L3)
        PHPT(L1,L2,L4,L3)=PHP(L1,L2,L4,L3)
100     CONTINUE
        ELSE
        DO 200 L1=1,2
        IOP(8+L1)       =IOPT(L1)
        DO 200 L2=1,IOP(8+L1)
        IHQP(L1,L2)=IHQPT(L1,L2)
        IDN(L1,L2)=IDNT(L1,L2)
        FMN(L1,L2)=FMNT(L1,L2)
        DO 234 J1=1,2
        KFEND(L1,L2,J1)=KFENDT(L1,L2,J1)
234     CONTINUE
        DO 200 L4=1,IHQPT(L1,L2)
        DO 200 L3=1,5
        PPS(L1,L2,L3)=PPST(L1,L2,L3)
        PPH(L1,L2,L3)=PPHT(L1,L2,L3)
        PPSY(L1,L2,L3)=PPSYT(L1,L2,L3)
        KHP(L1,L2,L4,L3)=KHPT(L1,L2,L4,L3)
        PHP(L1,L2,L4,L3)=PHPT(L1,L2,L4,L3)
200     CONTINUE
        ENDIF
        RETURN
        END
C*************************************************************************************
        
      SUBROUTINE FRINGEH

C........................This routine administrates radiation of nucleons with RPS
      PARAMETER (KSZJ=40000,KSZ1=30,KSZ2=300,KSZ4=150)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
        COMMON /TMP/ NMEM,MSTU24
       SAVE /FRINTN0/,/FRPARA1/,/LUJETS/

        KFRT=KFR(15)
        KFR(15)=0  
        NSTR = 0
        IOP(15)=0
      DO 390 L=1, 2
      DO 390 J=1,IOP(8+L)
      CALL FRTORST(L,J)
        IF(IOP(24).EQ.1)RETURN
      NSTR = NSTR+1
      NQG = N-NMEM
               IF( (KFR(1).EQ.1.AND.KFR(13).GE.1).AND.
     >             (NSTR.GT.100.OR.NQG.GT.(KSZJ-N)/10) ) THEN
      CALL LUEXEC
      IF(MSTU(24).EQ.4) MSTU24=1
 
          IF(KFR(13).LE.3) THEN
        CALL LUEDIT(KFR(13))
          ELSEIF(KFR(13).GE.4) THEN 
        CALL FREDITD()
          ENDIF
      NMEM=N
      NSTR=0
               ENDIF
390   CONTINUE
        KFR(15)=KFRT    

        RETURN
        END
C****************************************************************************************

        SUBROUTINE FRCRGLU
C........................This routine administrates radiation of firecracker gluons
      PARAMETER (KSZJ=40000,KSZ1=30,KSZ2=300,KSZ4=150,KSZ5=100)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRINTN1/PPS(2,KSZ2,5),PPH(2,KSZ2,5),PPSY(2,KSZ2,5),PPA(2,5)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      COMMON/FRCLUS1/NCL(2,KSZ4),CLM(2,KSZ4),NSL(2,KSZ4)
      COMMON/FRCLUS2/ICL(2,KSZ4),JCL(2,KSZ4),NC(2)
      COMMON/FRINTN2/NHP(2),IHQP(2,KSZ2),KHP(2,KSZ2,100,5),
     >   PHP(2,KSZ2,100,5)
       COMMON/HARD/IHQPS(2,KSZ2)
        COMMON/FRROPE1/NS,SIM2(KSZ5),RIM2
      COMMON/ARROPD/NR,RPT2(100),PTD2(100),IFT(100),PTEN2(2),NPR,NAC
      COMMON /ARDAT1/ PARA(40),MSTA(40)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
        COMMON/ARJETX/NO,KO(1000,5),PO(1000,5),VO(1000,5)
       COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)               
        DIMENSION KT(2000,5),PT(2000,5),VT(2000,5),KCS(300)
        DIMENSION KT0(1000,5),PT0(1000,5),KSS(300),KHH(300)
      COMMON/FCRSOUT/MCL(10),RCL(20),INEL(400),NST(2,KSZ4),ACL(2,KSZ4)
      SAVE /FRINTN0/,/FRINTN1/,/LUJETS/,/FRCLUS1/,/FRROPE1/,/HARD/,
     & /FRCLUS2/,/ARDAT1/,/FRPARA1/,/FRINTN2/,/ARROPD/,/ARJETX/
c..change event as if all evnets were soft 
C************************************************************************************
c..NS number of string in the rope proceeding radiation
c..NCL(2,L)  number of nucleons in the L-th cluster
C..NSL(2,L)  number of strings in the L-th cluster,(=NCL(2,L) if no colour
C  singlet in the cluster
c..CLM(2,L) inavriant mass of the L-th rope
C.if no colour singlet in a cluster,NS = NCL(2,L)
C**************************************************************************************
c..cluster ALL participants into clusters.Now event record is ordered cluster by 
c..cluster.Information of cluster is stored in /FRCLUS2/.ICL(2,KSZ4) is entry number
c.. of the first nucleon of KSZ4-th cluster while JCL(2,KSZ4) is  entry number
c.. of the last nucleon of KSZ4-th cluster,four-momentum of strings is used for
C..grouping

c cluster participant neucleons from projectile and target separately by their
c distance in transverse momentum space,update event record (/FRINTN1/,
c /FRINTN2/,/FRINTN3/,/FRINTN4/)cluster by cluster.
        DO 7 L=1,2
        DO 310 I1=1,IOP(8+L)
        IHQPS(L,I1)=IHQP(L,I1)
        IHQP(L,I1)=0
310     CONTINUE
7       CONTINUE
c::c through FRTORST(LH,J) event record in COMMON/FRINTN1[PPS(),PPH()]....is changed
c..into COMMON/LUJETS [K(),P()...].At the moment,only soft part of string is delt with
c..since we assume only soft excitation energy is used for emitting firecracker gluons

        NB=0
        N0=0

        KFR2=KFR(2)
        MSTA17=MSTA(17)
      DO 150 LH=1, 2
      DO 150 LJ=1,NC(LH)
        NO=0
        N=0
         IOP(15) = 0
        JJ=1
        I0=0
      DO 100 J=ICL(LH,LJ),JCL(LH,LJ)
        NX=N
        KCS(J)=0
        KSS(JJ)=0
        KHH(JJ)=0
        CALL FRTORST(LH,J)
c to judge if a participant nucleon is extended to be a string or 
c remains to be a nucleon
        IF((N-NX).EQ.2.AND.IHQPS(LH,J).EQ.0)THEN
        KSS(JJ)=1
        JJ=JJ+1
        ELSEIF((N-NX).EQ.1)THEN
        KCS(J)=1
        I0=I0+1
        KHH(I0)=0
        IF(IHQPS(LH,J).GT.0)KHH(I0)=1
        ELSE
        JJ=JJ+1
        ENDIF
100   CONTINUE
c event record in  COMMON/LUJETS [K(),P()...has changed into COMMON
c./ARJETX/[KO(),PO()...] for being used in Ariadne.here all event record
c.  in a cluster are sent into Ariadne as a whole
c:: pick firecracker gluons from clusters only and assign them to string 
c::if kinematically possible
c::NCN number of colour singlet in the cluster
c::initialize so that only fire cracker gluon is emited
        KFR(2)=0
        MSTA(17)=0
        NCN=0
c--deal with firecracker gluon radiation from clusters
        CALL FRARCLU(NCN)
c--event record  COMMON/LUJETS/ is obtained after a cluster finishes 
c-radiation and event record is temperarily stored in KT(),PT() and VT()

        IO=0
c::store firecracker gluon into KO() and PO()
        DO 145 L=1,N
        IF(K(L,2).EQ.21)THEN
        IO=IO+1
        DO 146 L1=1,5
        KO(IO,L1)=K(L,L1)
        PO(IO,L1)=P(L,L1)
146     CONTINUE
        ENDIF
145     CONTINUE
        kg=io
        IT=0
c::store colour-singlet nucleon from soft string into KT0 and PT0,entry of a hard 
c::string which is shrunk into a nucleon is still unchanged
        DO 135 IO=1,NCN 
        IF(KHH(NCN-IO+1).EQ.1)GOTO 135
        IT=IT+1
        if(((N0+IT).gt.1000).or.((N-IO+1).gt.20000).or.((N-IO+1).le.0))
     &  write(MSTU(11),*)'segment fault','N0,IT,NCN,N',N0,IT,NCN,N

        DO 1350 L1=1,5  
        KT0(N0+IT,L1)=K(N-IO+1,L1)
        PT0(N0+IT,L1)=P(N-IO+1,L1)
1350    CONTINUE
135     CONTINUE
        N0=N0+IT
        IO=0
        kg1=0
        IOO=1
        DO 140 L=ICL(LH,LJ),JCL(LH,LJ)
        MS=L-ICL(LH,LJ)+1-IO
        IF(KCS(L).EQ.1)THEN
        IO=IO+1
c:::put a firecracker gluon into a hard string being last hard gluon entry,later on the 
c::gluon either is looked upon as a hard gluon if it is hardest of all or a bremsstrahlung
c::gluon
        ELSEIF(IHQPS(LH,L).GT.0.AND.IFT(MS).EQ.1)THEN
        kg1=kg1+1
        IHQPS(LH,L)=IHQPS(LH,L)+1
        PPS(LH,L,1)=PPS(LH,L,1)-PO(IOO,1)
        PPS(LH,L,2)=PPS(LH,L,2)-PO(IOO,2)
        PPS(LH,L,3)=PPS(LH,L,3)-(PO(IOO,4)-PO(IOO,3))
        PPS(LH,L,4)=PPS(LH,L,4)-(PO(IOO,4)+PO(IOO,3))
        PPS(LH,L,5)=0.0
        KHP(LH,L,IHQPS(LH,L),1)=21
        KHP(LH,L,IHQPS(LH,L),2)=21
        KHP(LH,L,IHQPS(LH,L),3)=LH
        KHP(LH,L,IHQPS(LH,L),4)=0
        KHP(LH,L,IHQPS(LH,L),5)=0
        DO 136 L1=1,5   
        PPH(LH,L,L1)=PPSY(LH,L,L1)-PPS(LH,L,L1)
        PHP(LH,L,IHQPS(LH,L),L1)=PO(IOO,L1)
136     CONTINUE
        LCOM=1
        LFIN=IHQPS(LH,L)
c check if the fircracker gluon is the hardest on the string
        PHD0=SQRT(PHP(LH,L,LFIN,1)**2+PHP(LH,L,LFIN,2)**2)
        DO LHD=1,IHQPS(LH,L)-1
        PHD=SQRT(PHP(LH,L,LHD,1)**2+PHP(LH,L,LHD,2)**2)
        IF(PHD0.LT.PHD)LCOM=0
        ENDDO
        IF(LCOM.EQ.1)MCL(5)=MCL(5)+1
C MCL(5) is the number of gluons which have been accepted as hard
c gluons on a string.
        IOO=IOO+1
        ELSEIF(IHQPS(LH,L).EQ.0.AND.IFT(MS).EQ.1)THEN
        IOO=IOO+1
        ELSE
        ENDIF
140     CONTINUE        
        IOO=1
        IO=0
        kg2=0
c::store soft string and their fire-cracker gluon into KO,PO
        DO 155 L=1,N-NCN
        IF(K(L,1).EQ.2.AND.KSS(IOO).EQ.1.AND.K(L,2).NE.21)THEN
        IF(IFT(IOO).EQ.1)THEN
        kg2=kg2+1
        MA=3
        ELSE
        MA=2
        ENDIF
        DO 156 L2=1,MA
        IO=IO+1
        DO 156 L1=1,5
        KO(IO,L1)=K(L2+L-1,L1)
        PO(IO,L1)=P(L2+L-1,L1)
156     CONTINUE
        IOO=IOO+1
        ELSEIF(K(L,1).EQ.2.AND.KSS(IOO).EQ.0.AND.K(L,2).NE.21)THEN
        IOO=IOO+1
        ELSE    
        ENDIF
155     CONTINUE
        DO 130 L=1,IO
        DO 130 L1=1,5
        KT(NB+L,L1)=KO(L,L1)
        PT(NB+L,L1)=PO(L,L1)
        VT(NB+L,L1)=VO(L,L1)                                                                     
130     CONTINUE
        IOP(19)=NAC+IOP(19)     
        NSL(LH,LJ)=NS
        CLM(LH,LJ)=SQRT(RIM2)   
        NB=NB+IO
        if(kg.ne.(kg1+kg2))THEN
c :: kg is nuber of gluons from firecracker.kg1 number of gluons for hard
c::string and kg2 for soft string.NAC=kg
         write(MSTU(11),*)'gluons from firecracker are not properly 
     &   put onto hard and soft strings,some information is outputed'
        write(MSTU(11),*)'kg,kg1,kg2,NAC',kg,kg1,kg2,NAC
        ENDIF
150     CONTINUE
        KFR(2)=KFR2
        MSTA(17)=MSTA17
        DO 190 L=1,2
        DO 190 I1=1,IOP(8+L)
        IHQP(L,I1)=IHQPS(L,I1)
190     CONTINUE
        CALL FRINGES(KT,PT,VT,KT0,PT0,NB,N0)
        RETURN
        END


C****************************************************************************
        SUBROUTINE FRINGES(KT,PT,VT,KT0,PT0,NB,N0)
C........................This routine administrates radiation of nucleons without RPS
C copy strings back to an old entry after all clusters have finished firecracker
c..radiation
      PARAMETER (KSZJ=40000,KSZ1=30)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      COMMON/ARJETX/NO,KO(1000,5),PO(1000,5),VO(1000,5)
       COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)               
        DIMENSION KT(2000,5),PT(2000,5),VT(2000,5)
        DIMENSION KT0(1000,5),PT0(1000,5)
        save /ARJETX/
        KFR15=KFR(15)
        KFR(15)=0
        NO=0
        N=0
        DO 180 L=1,NB
        NO=NO+1
        DO 170 L1=1,5
        KO(NO,L1)=KT(L,L1)
        PO(NO,L1)=PT(L,L1)
        VO(NO,L1)=VT(L,L1)  
170     CONTINUE
        IF(KT(L,1).EQ.1)THEN
        DO I=1,NO
         IF(KO(I,2).EQ.21) PO(I,5)=0.0
        ENDDO
        IF(KFR(2).EQ.1)CALL AREXEC 
      DO 400 IO=1,NO
      IF(KO(IO,1).GE.11) GOTO 400
      N=N+1
      DO 450 LO=1,5
      K(N,LO) = KO(IO,LO)
450   P(N,LO) = PO(IO,LO)
400   CONTINUE
        NO=0
        ENDIF
180     CONTINUE
        KFR(15)=KFR15

c...insert colour-singlet nucleons into event record
      DO 300 IO=1,N0
      DO 300 LO=1,5
      K(N+IO,LO) = KT0(IO,LO)
      P(N+IO,LO) = PT0(IO,LO)
300     CONTINUE
        N=N+N0
        RETURN
        END

C*********************************************************************

C************************************************************************
          SUBROUTINE  FRARORDT

      PARAMETER (KSZJ=40000,KSZ1=30, KSZ2=300,KSZ4=150)

      COMMON/FRCLUS1/NCL(2,KSZ4),CLM(2,KSZ4),NSL(2,KSZ4)
      COMMON/FRCLUS2/ICL(2,KSZ4),JCL(2,KSZ4),NC(2)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/FRTEMP1/NOICF(2,KSZ2,100)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
      COMMON/FRINTN2/NHP(2),IHQP(2,KSZ2),KHP(2,KSZ2,100,5),
     >   PHP(2,KSZ2,100,5)

         DIMENSION KHPT(250,25,5),IHQPT(2,KSZ2),
     &   PHPT(250,25,5),JC(2)
      SAVE /FRINTN0/,/FRCLUS1/,/FRCLUS2/,/FRINTN2/,
     &  /FRTEMP1/  

c..purpose: event record is reproduced with nucleons in the same group 
c..being placed together
c..ICL(2,KSZ4)entry number of the first nucleon in KSZ4_th group
c..JCL(2,KSZ4)entry number of the last  nucleon in KSZ4_th group
        JC(1)=0
        JC(2)=0
        DO 400 L=1,2
        IF(NC(L).EQ.0)GOTO 400
        DO 300 I1=1,NC(L)

        DO 300 I2=1,NCL(L,I1)
        JC(L)=JC(L)+1
        IF(NOICF(L,I1,I2).EQ.0)then
        WRITE(MSTU(11),*)'an error may have taken place 
     &          in subroutine FRARORDT,some information is outputed'
        DO 10 J1=1,2
        WRITE(MSTU(11),*)'NC(L)=',NC(J1),IOP(8+J1)
        DO 10 J2=1,NC(J1)
        WRITE(MSTU(11),*)'NCL(J1,J2)=',NCL(J1,J2)
        DO 10 J3=1,NCL(J1,J2)
        WRITE(MSTU(11),*)'NOICF(J1,J2,J3)=',NOICF(J1,J2,J3)
10      CONTINUE
        ENDIF
        IEND=IHQP(L,NOICF(L,I1,I2))
        DO 100 I3=1,IEND
        DO 150 I4=1,5
        KHPT(JC(L),I3,I4)=KHP(L,NOICF(L,I1,I2),I3,I4)
        PHPT(JC(L),I3,I4)=PHP(L,NOICF(L,I1,I2),I3,I4)
150     CONTINUE
100     CONTINUE
        IHQPT(L,JC(L))=IHQP(L,NOICF(L,I1,I2))

300     CONTINUE
        DO 330 I1=1,IOP(8+L)
        IEND=IHQPT(L,I1)
        DO 250 I3=1,IEND
        DO 200 I4=1,5
        KHP(L,I1,I3,I4)=KHPT(I1,I3,I4)
        PHP(L,I1,I3,I4)=PHPT(I1,I3,I4)
200     CONTINUE
250     CONTINUE
        IHQP(L,I1)=IHQPT(L,I1)

330     CONTINUE
400     CONTINUE
        RETURN     
        END
C***********************************************************************
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      SUBROUTINE FRARSOR
      PARAMETER (KSZ1=30,KSZ5=100)
      IMPLICIT DOUBLE PRECISION (D)

C...starts from biggest invariant masse
C***********************************************************************
C NP(L) is numbers of partons in Lth string.
C SIM2(L) is invariant mass squared of Lth string
C RIM2    is invariant mass squared of colour rope 
C NS is the numbers of strings that will be fragmented
C NPA is the total numbers of partons 
C**********************************************************************   

      COMMON/ARJETX/ N,K(1000,5),P(1000,5),V(1000,5)
      SAVE /ARJETX/
      COMMON/FRROPE1/NS,SIM2(KSZ5),RIM2
        SAVE /FRROPE1/
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      SAVE /LUDAT1/,/LUDAT2/
      DIMENSION DPS(5),ROP(5),NP(100),PS1(5),PS2(5)
         I=0
        NPA=0
        NS=0
        DO 40 L=1,100
        NP(L)=0
40      CONTINUE
        DO 60 L=1,5      
        ROP(L)=0
60      CONTINUE
C...Reset counters.Identify parton system finding starting point of a string
80     KQSUM=0
        NS=NS+1 
        DO 100 L=1,5      
        DPS(L)=0
100     CONTINUE
120     I=I+1
        IF(I.GT.N) THEN
        NS=NS-1
        GOTO 180
        ENDIF
      IF(K(I,1).NE.1.AND.K(I,1).NE.2.AND.K(I,1).NE.41) GOTO 120
      KC=LUCOMP(K(I,2))
      IF(KC.EQ.0) GOTO 120
      KQ=KCHG(KC,2)*ISIGN(1,K(I,2))
      IF(KQ.EQ.0) GOTO 120
C...Starting point has been found .Take copy of partons to be considered. 
C...Check flavour sum.
      NP(NS)=NP(NS)+1
      DO 140 J=1,5
140   IF(J.NE.4) DPS(J)=DPS(J)+P(I,J)
      DPS(4)=DPS(4)+SQRT(DBLE(P(I,1))**2+DBLE(P(I,2))**2+
     &DBLE(P(I,3))**2+DBLE(P(I,5))**2)
      IF(KQ.NE.2) KQSUM=KQSUM+KQ
      IF(K(I,1).EQ.41) THEN
        KQSUM=KQSUM+2*KQ
      ENDIF
C...  loop back going through other partons on a string.
      IF(K(I,1).EQ.2.OR.K(I,1).EQ.41) GOTO 120
        SIM2(NS)=DPS(4)**2-DPS(1)**2-DPS(2)**2-DPS(3)**2
        SIM2(NS)=MAX(0.,SIM2(NS))
C...  Another end of a string has been found.check if the string is colour
C... white.
      IF(KQSUM.NE.0) THEN
        CALL LUERRM(12,'(FRARSOR:) unphysical flavour combination')
        IF(MSTU(21).GE.1) RETURN
      ENDIF
        NPA=NPA+NP(NS)
      DO 160 J=1,5
160   ROP(J)=ROP(J)+DPS(J)
        GOTO 80
180     IF(NS.LE.1)GOTO 430
        DO 310 L1=1,5
        PS1(L1)=0
        PS2(L1)=0
310     CONTINUE
        
C...add a colourrope as first dipole in ARJETX

330     DO 340 L=N,1,-1
        DO 340 L1=1,5
        K(L+2,L1)=K(L,L1)
        P(L+2,L1)=P(L,L1)
        V(L+2,L1)=V(L,L1)    
c..calculate four-momentum of two ends of colourrope
        IF(K(L,1).EQ.2)PS1(L1)=PS1(L1)+P(L,L1)
        IF(K(L,1).EQ.1)PS2(L1)=PS2(L1)+P(L,L1)                                                                 
340     CONTINUE
        DO 360 L=1,5
        P(1,L)=PS1(L)
        P(2,L)=PS2(L)
360     CONTINUE
        N=N+2
430     RIM2=ROP(4)**2-ROP(3)**2-ROP(2)**2-ROP(1)**2
        RIM2=MAX(0.,RIM2)
        RETURN
        END
C...****************************************************************
      SUBROUTINE FRARGLU
      PARAMETER (KSZJ=40000,KSZ1=30, KSZ2=300,KSZ5=100)
      COMMON/FRROPE1/NS,SIM2(KSZ5),RIM2
      COMMON/ARROPD/NR,RPT2(100),PTD2(100),IFT(100),PTEN2(2),NPR,NAC
      COMMON/ARDAT1/ PARA(40),MSTA(40)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/ARJETX/N,K(1000,5),P(1000,5),V(1000,5)
      SAVE /ARDAT1/,/ARROPD/,/ARJETX/

C...Decide phase spaces of all dipoles for checking if a gluon from rope
c.. can be given to a dipole
c**********************************************************************
c PETM2 maximum Pt square in phase space for an extended dipole
c PPTM2 maximum Pt square in phase space for a point-like dipole
C PTEN2() are not used actually
C NAC number of gluons from fire-cracker that are accepted by
c string

      DO 50 I=1,NS
       NER1=MOD(K(2*I-1+2,4),10)
       NER3=MOD(K(2*I+2,4),10)
        S=SIM2(I)
      IF(NER1.GT.0.AND.NER3.EQ.0) PETM2=((0.25*S*(PARA(10+NER1)**
     $                          PARA(10)))**(2.0/(2.0+PARA(10))))
      IF(NER1.EQ.0.AND.NER3.GT.0) PETM2=((0.25*S*(PARA(10+NER3)**
     $                          PARA(10)))**(2.0/(2.0+PARA(10))))
      IF(NER1.GT.0.AND.NER3.GT.0) PETM2=((0.25*S*((PARA(10+NER1)*
     $           PARA(10+NER3))**PARA(10)))**(1.0/(1.0+PARA(10))))
        PPTM2=0.25*SIM2(I)
        PETM2=MIN(PPTM2,PETM2)
        PETM2=MAX(PARA(3)**2,PETM2)
        IF(KFR(17).EQ.1)PETM2=PARA(3)**2
        IF(I.EQ.1)PTEN2(1)=0.25*SIM2(I)
        IF(I.EQ.NS)PTEN2(2)=PETM2
70      IF(RPT2(I).LT.PPTM2.AND.RPT2(I).GT.PETM2)THEN
        NAC=NAC+1
        IFT(I)=1
        PTD2(I)=RPT2(I)
        ELSE
        IFT(I)=0
        PTD2(I)=0.
        ENDIF
50      CONTINUE
        RETURN 
        END


                 SUBROUTINE FRRESEXE
c       Fritiof + rescattering.
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       Rescattering part of this program was written by Wang Zhong-Qi
c        on 12/05/91,improved by Sa Ben-Hao on 02/05/93, 
c        improved further to be able to connect with Fritiof 7.0 by
c        Sa Ben-Hao on 10/03/93,bring to Lund and improved again by 
c        Sa Ben-Hao and Tai An in collaboration with Bo Andersson and 
c        G. Gustafson on 11/22/93.
c       final modification is finished in Apr.,94
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      PARAMETER (KSZJ=40000,KSZ1=30,nsize=100000)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
        DOUBLE PRECISION DB(4),DBB
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      SAVE /FRPARA1/,/FRINTN0/,/LUJETS/
        common/sa10/csnn,cspin,cskn,cspipi,rcsit,ifram
        common/sa12/psa(5),ptai(5),clorenp,clorent
         save /sa12/    
        dimension lc(nsize,5),tc(nsize),tw(nsize)
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C THIS ROUTINE ADMINISTATE THE COMPLETE RESCATTERING
C OF A EVENT. THE RESCATTERING IS DONE IN THE CMS FRAME
C OF TWO NUCLEI.
c       ifram = 0 for fixed target exp., = 1 for collider exp.
c       cspipi (fm^2) is total cross section of pion + pion
c       sig (fm^2) is the cross section of pion + pion to kaon + kaon
c       cspin (fm^2) is total cross section of pion + nucleon interaction
c       cskn (fm^2) is total cross section of kaon + nucleon interaction
c       csnn (fm^2) is total cross section of n + n interaction
c       rcsit : ratio of total cross section to inelastic
c       kfmax : the maximum # of particles with different kf code tagged 
c        in the program.
c       kfaco(i) : KF code of i-th kind of particle among kfmax
c       numb(i) : cumulated # of particles cooresponding to the kfaco(i).
c       last line of i-th kind of particle of flavor kfaco(i) in event record
c       disbe(i,j) : minimum allowable distance between particles
c               kfaco(i) & kfaco(j).
c       c17(i,1-3) : coordinates of the particle i with origin at the center 
c       of the target nucleus at current time of the colliding system
c       tlco(i,1-4) the space-time coordinates of a particle at its own 
c       freeze-out time. They are not used because it is hard to defined
c       freeze-out in the program 
c       tp(i) : time of particle started from collision of two nuclei
c       ishp(i)=1 if i-th particle inside the simulated volume
c              =0 if i-th particle outside the simulated volume
c       tau(i) : formation time of particle i.
c       p(i,1-4) : four momentum of particle i
c       p(i,5) : mass of particle i
c       isinel(i) = 0 without i-th special inelastic process
c                 = 1 with i-th special inelastic process
c       nctl0    number of initial collision pairs
c       noinel(400) frequency of rescattering occur through certain channel
c       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c       0. is the hard distance between two pions
c       0.5 is the hard distance between two nucleons
c       0. is the hard distance between pion and nucleon
c       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        iframt=ifram

        ntest=1
        CALl sysini(ntest)
c       set value of parameters at the begining,ntest=0 rescattering effect
c       is not considered for very peripheral collision,before rescattering
c       Ks0,lamda and lamdaba etc. are not decayed and Ks0,KL0 have been changed
c       into K0,K0ba
        if(ntest.eq.0)return
        ij=1
        i=NFR(1)
c the rescattering is performed in the CMS frame of the colliding system.
        if(KFR(28).eq.0)ifram=1
cboost LUJETS into CMS of collision system. 
        P(N,5)=SQRT(MAX(-P(N,1)**2-P(N,2)**2-P(N,3)**2+P(N,4)**2,0.0))
        P(N-1,5)=SQRT(MAX(-P(N-1,1)**2-P(N-1,2)**2-P(N-1,3)**2
     &  +P(N-1,4)**2,0.0))
        DB(1)=PLI0(1,1)*IOP(3)+PLI0(2,1)*IOP(5)
        DB(2)=PLI0(1,2)*IOP(3)+PLI0(2,2)*IOP(5)
        DB(3)=0.5*(PLI0(1,4)-PLI0(1,3))*IOP(3)+0.5*(PLI0(2,4)-
     &  PLI0(2,3))*IOP(5)
        DB(4)=0.5*(PLI0(1,4)+PLI0(1,3))*IOP(3)+0.5*(PLI0(2,4)+
     &  PLI0(2,3))*IOP(5)
        DB(1)=-DB(1)/DB(4)
        DB(2)=-DB(2)/DB(4)
        DB(3)=-DB(3)/DB(4)
        DBB=1.D0-DB(1)*DB(1)-DB(2)*DB(2)-DB(3)*DB(3)
        cloren=1.D0/DSQRT(MAX(1.D-5,DBB))
        EP=0.5*(PLI0(1,4)+PLI0(1,3))*IOP(3)
        PP=0.5*(PLI0(1,4)-PLI0(1,3))*IOP(3)
        ET=0.5*(PLI0(2,4)+PLI0(2,3))*IOP(5)
        PT=0.5*(PLI0(2,4)-PLI0(2,3))*IOP(5)
        EPCM=cloren*(EP+DB(3)*PP)
        ETCM =cloren*(ET+DB(3)*PT)
        clorenp=EPCM/(IOP(3)*0.938)
        clorent=ETCM/(IOP(5)*0.938)
        if(ifram.eq.0)then
        clorent=1
        vb=(PLI0(1,4)-PLI0(1,3))/(PLI0(1,4)+PLI0(1,3))
        clorenp=1./sqrt(MAX(1.D-5,1.-vb**2))
        endif
        if(KFR(28).eq.0)CALL LUDBRB(1,N,0.0,0.0,DB(1),DB(2),DB(3))

c*********************************************************************
        call even_init(lc,tc,tw,i)
        time=0.
c initialize for proceeding rescattering among final particles and 
c spectator nucleons. The spectator nucleons have given time and 
c coordinates , four-momentum and put into the particle list one by one
        call rescat(i,ij,time)

C start rescattering


c add the psa() as the last item in order to boost 
        DO LL=1,4
        P(N+1,LL)=psa(LL)
        ENDDO
cboost LUJETS back to lab system.
        P(N,5)=SQRT(MAX(-P(N,1)**2-P(N,2)**2-P(N,3)**2+P(N,4)**2,0.0))
        P(N-1,5)=SQRT(MAX(-P(N-1,1)**2-P(N-1,2)**2-P(N-1,3)**2
     &  +P(N-1,4)**2,0.0))
        if(KFR(28).eq.0)CALL LUDBRB(1,N+1,0.0,0.0,-DB(1),-DB(2),-DB(3))
        ifram=iframt
        DO LL=1,4
        psa(LL)=P(N+1,LL)
        ENDDO
c       DO LL=1,4
c       P(N+1,LL)=0.
c       ENDDO
c put the spectator nucleons without experiencing inelastic collision back to
c one entry
        call RESCSET2(time)
c       N=N+1
c       P(N,5)=SQRT(P(N,4)**2-P(N,1)**2-P(N,2)**2-P(N,3)**2)
c       K(N,1)=1
c       K(N,2)=20000
c       K(N,3)=0
c       K(N,4)=0
c       K(N,5)=0
c       N=N-1
        return  
        end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        subroutine sysini(ntest)
c       initiate the rescattering program. give the initial values to
c       the quantities needed in the calculation 
        parameter (KSZ1=30,KSZJ=40000)
        COMMON/LUCIDAT1/KFACOT(50),DISDET(50),ISINELT(400)
        common/sa5/kfmax,kfaco(50),numb(50),disbe(50,50)
        common/count/isinel(400)
        COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/LUCIDAT2/KFMAXT,PARAM(20),WEIGH(400)
        COMMON/INDATA/NAP,NZP,NAT,NZT,BP
        common/sa6/kfmaxi,nwhole,nna
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn
     c  ,rnt,rnp,rao,rou0
        common/sa10/csnn,cspin,cskn,cspipi,rcsit,ifram
        COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
        SAVE /LUCIDAT1/,/sa5/,/LUCIDAT2/,/count/,/INDATA/,/sa6/,
     &        /papr/,/sa10/,/LUJETS/,/FRINTN0/
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c       edipi : interaction distance between two pions.
c       epin : interaction distance between pion and nucleon.
c       ekn : interaction distance between kaon and nucleon.
c       ecsnn : interaction distance between two nucleons.
c       nap (nzp) : # of nucleons (protons) in projectile nucleus.
c       nat (nzt) : # of nucleons (protons) in target nucleus.
c       t0 : average proper formation time.
c       ddt : time accuracy
c       dep : the accuracy in four momentum conservation
c       rou0 : normal nucleon density.
c       rao : enlarged factor for the radius of target nucleus. 
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        nap=IOP(3)
        nat=IOP(5)
        anat=nat
        anap=nap
c       in the program the x-sections are given in a unit of fm*fm.
        csnn=PARAM(1)*0.1
        cspin=PARAM(2)*0.1
        cskn=PARAM(3)*0.1
        cspipi=PARAM(4)*0.1
        rou0=PARAM(11)
c       Considering the nucleus as a sphere with radii rnt for target
c        and rnp for projectile.
        rnt=(3.*anat/(4.*3.1415926*rou0))**(0.33333)
        rnp=(3.*anap/(4.*3.1415926*rou0))**(0.33333)
        bp=aop(2)
c       bp : the impact parameter of current event.
        rnpt=rnp+rnt
C       if collision is too peripheric,go back
        if(bp.ge.0.95*rnpt)then
        if(abs(k(n,2)).ge.10000)
     &  k(n,3)=iop(5)-iop(10)-(abs(k(n,2))-10000)
        if(abs(k(n-1,2)).ge.10000)
     &  k(n-1,3)=iop(3)-iop(9)-(abs(k(n-1,2))-10000)
        ntest=0
C k(n,3) and k(n-1,3) are the neutron numbers of the spectator clusters.
        RETURN  
        endif

c given the collision distance of two colliding particles.
        edipi=sqrt(cspipi/3.1416)
        epin=sqrt(cspin/3.1416)
        ekn=sqrt(cskn/3.1416)
        ecsnn=sqrt(csnn/3.1416)
C       set initial values to some quantities
        nzp=IOP(4)
        nzt=IOP(6)
        sig=PARAM(5)*0.1
        rcsit=PARAM(6)
        t0=PARAM(7)
        dep=PARAM(9)
        ddt=PARAM(8)
        rao=PARAM(10)
        kfmax=KFMAXT
        do i=1,50
        kfaco(i)=KFACOT(i)
        enddo
        do j=1,400
        isinel(j)=ISINELT(j)
        enddo
        do i=1,50
        do j=1,50
        disbe(i,j)=0.
        enddo
        enddo
        do j=1,kfmax
c 1 to 4 are p,n,pbar,nbar. 26-29 are four Delta resonances.
        disbe(1,j)=DISDET(j)
        disbe(2,j)=DISDET(j)
        disbe(3,j)=DISDET(j)
        disbe(4,j)=DISDET(j)
        disbe(26,j)=DISDET(j)
        disbe(27,j)=DISDET(j)
        disbe(28,j)=DISDET(j)
        disbe(29,j)=DISDET(j)
        enddo
400     do i=1,49
        do j=i+1,50
        disbe(j,i)=disbe(i,j)
        enddo
        enddo
        kfmaxi=kfmax
C change Ks,KL from Fritiof into K0,K0ba
        do j=1,n
        kf=k(j,2)
        if(kf.eq.130 .or. kf.eq.310)then
        rrlu=rlu(1)
        k(j,2)=311
        if(rrlu.gt.0.5)k(j,2)=-311
        endif
        enddo

        return
        end
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


        subroutine rescat(jjj,iii,time)
c       This subroutine is used to simulate the rescattering. 
        parameter(KSZJ=40000,KSZ1=30)
        parameter(nsize=100000)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/ctllist/nctl,noinel(400),nctl0
          save /ctllist/
        common/sa5/kfmax,kfaco(50),numb(50),disbe(50,50)
        common/sa6/kfmaxi,nwhole,nna
        common/tai/NT,KT(500,5),PT(500,5),VT(500,5)
        common/sa12/psa(5),ptai(5),clorenp,clorent
        common/sa4/tau(kszj),tlco(kszj,4)
         save /sa12/
        double precision b(3)
        common/sa10/csnn,cspin,cskn,cspipi,rcsit,ifram
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        dimension pi(4),pj(4),pii(4),pjj(4)
        dimension noinelt(400)
        integer winel
        do i=1,400
        noinelt(i)=0
        enddo   
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       lc(i,1) and lc(i,2) are,respectively,the order # in particle  
c       list of colliding particles in i-th collision pair.
c       lc(i,3) and lc(i,4) are the flavor codes of scattered particles
c       in i-th collision.
c       lc(i,5) identify the different inelastic processes.
c       tc(i) is the collision time of i-th colli.
c       tw(i) is the cross section ratio of (i-th inelas.)/tot
c::nctl0    number of initial collision pairs
c ptai(1),ptai(2),ptai(3),ptai(4) is the number of elastic scattering,
c inelastic scattering, inelastic scattering but treated as elastic
cscattering and total rescattering in the agregate.

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C       dump those particles which are not involved in the rescattering
c       process into a temporatory array KT etc.
        NT=N-numb(kfmax)
        DO I=1,NT
        DO J=1,5
        KT(I,J)=K(numb(kfmax)+I,J)
        PT(I,J)=P(numb(kfmax)+I,J)
        VT(I,J)=V(numb(kfmax)+I,J)
        ENDDO
        ENDDO
        n=numb(kfmax)
        time=0.
        call copl(time)
        call ctlcre(lc,tc,tw)
c       create initial colli. list and  fill up lc(i,1-2),tc(i).
        nctl0=nctl
        iii=1
10      if(iii.eq.1)goto 1000
        call copl(time)
1000    call find(icp,tcp,lc,tc,tw)
c       find out the binary colli. with minimum collsion time
        if(icp.eq.0) goto 100
c       icp=0 means the collision list is empty
        time=tcp
        l=lc(icp,1)
        l1=lc(icp,2)
        ilo=0
        pi(4)=p(l,4)
        pj(4)=p(l1,4)
        do i=1,3
        pi(i)=p(l,i)
        pj(i)=p(l1,i)
        b(i)=(pi(i)+pj(i))/(pi(4)+pj(4))
        enddo
        call lorntz(ilo,b,pi,pj)
C       boost to CMS frame
        ss=pi(4)+pj(4)
500     continue
        call his(time,lc,tc,tw,istop)
c       perform classical Newton motion
        if(istop.eq.1)goto 100
c       istop=1 means all particles have get out of considered volume
        ww=rcsit
        mflag=0 
c       assume the cross section ratio of (ela.)/tot =1- rcsit
        rrlu=rlu(1)
        if(rrlu.gt.ww)then
        ptai(1)=ptai(1)+1
        winel=0
        goto 640
        endif
700     winel=1
C       two particles with four-momentum pi and pj in CMS frame and flavor
C       k(l,2),k(l1,2) may go through inelastic reaction
        call coinel(l,l1,ss,b,pi,pj,icp,pii,pjj,lc,tc,tw,winel,ik3,ik4)

C       decided if one channel of inelastic reaction is occured
C       winel=0 means that the inelastic reaction has not been really happened
C       because either CMS energy is too small or none of channels
C       that we are interested in has taken place. If inelastic 
c       collision happens ik3 and ik4 are new flavors of the produced 
c       particles from the inelastic process. pii and pjj are four-momentum 
c       of two produced particles in Lab frame and the two colliding 
c       particles are  still at positions labeled by l and l1 in the 
c       particle list.
        mflag=1 
640     if(winel.ne.0)then
        ptai(2)=ptai(2)+1
c       treat the inelastic collision processe.
        icp5=lc(icp,5)
        noinelt(icp5)=noinelt(icp5)+1

c       if(icp5.eq.397 .or. icp5.eq.398 .or. icp5.eq.399.or.icp5.eq.400)
c     &  then


        call updpli(l,l1,icp,ss,pii,pjj,lc,tc,tw,winel,time,icp5)
c l and l1 are now the postions of the two produced particles in the 
c particle list.
        l=lc(icp,1)
        l1=lc(icp,2)
c       it uses to update particle list after inelastic collision,
c        truncates collision list correspondingly and update collision
c        list for the case of n(nbar) annihilation
c       if(icp5.eq.397 .or. icp5.eq.398 .or. icp5.eq.399.or.icp5.eq.400)
c     & goto 300
c  this part of statements are left here for testing the program
        if(jjj.eq.-391.or.jjj.eq.-1052)then
        write(mstu(11),*)'inelastic',l01,l02,k01,k02,'after',l,l1,
     &  k(l,2),k(l1,2),icp5,'iii=',iii,'n=',n,kfmax,nwhole-n,'jjj=',jjj
        do i=1,kfmax
        write(mstu(11),*)i,numb(i)
        enddo
        call lulist(1)
        endif
        k(l1,3)=0
        k(l,3)=0

        else
        if(mflag.eq.1)ptai(3)=ptai(3)+1
        call coelas(l,l1,ss,pi,pj)

C       calculate four-momentum of two particles after elastic 
c reaction pi and pj in CMS frame

        call updple(l,l1,b,pi,pj,time)



c       update the particle list for elastic scattering,
c       pi and pj have been boosted back to Lab fram .

        endif

400     call updatl(l,l1,time,lc,tc,tw,winel,n)
c       update the collision list update the concerned variables in particle list
c       after both an inelastic and elastic reaction.
        ptai(4)=ptai(4)+1
300     iii=iii+1

        if(iii.gt.100*(nctl0+100))then
        write(mstu(11),*)'infinite loop may have happened in
     &         subroutine rescat'
        stop 'infinite loop occurs'
        endif
        goto 10
100     do i=1,400
        noinel(i)=noinelt(i)+noinel(i)
        enddo

        DO I=1,NT
        DO J=1,5
        K(N+I,J)=KT(I,J)
        P(N+I,J)=PT(I,J)
        V(N+I,J)=VT(I,J)
        ENDDO
        ENDDO
        n=nna
1200    return
        end
C************************************************************************
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


        subroutine even_init(lc,tc,tw,ijk)
c       initiate the event after fritiof
        parameter (mcludi=40000)
      PARAMETER (KSZJ=40000,KSZ1=30)
        parameter(nsize=100000)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      common/sa5/kfmax,kfaco(50),numb(50),disbe(50,50)
        common/sa6/kfmaxi,nwhole,nna
        common/sa10/csnn,cspin,cskn,cspipi,rcsit,ifram
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
      COMMON/INDATA/NAP,NZP,NAT,NZT,BP
        COMMON/LUCIDAT2/KFMAXT,PARAM(20),WEIGH(400)
        SAVE/LUCIDAT2/
        common/wz/c17(mcludi,3),ishp(mcludi),tp(mcludi)
     c   ,coor(3)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn
     c  ,rnt,rnp,rao,rou0
        common/ctllist/nctl,noinel(400),nctl0
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        dimension pp(4),pt(4),ppp(250,5),ppt(250,5)
c       pp : the rest four momentum of projectile nucleus
c       pt : the rest four momentum of target nucleus
c       ppp : four momentum of projectile spectator nucleons
c       ppt : four momentum of target spectator nucleons
        nctl=0
        do i=1,nsize
        do j=1,5
        lc(i,j)=0
        enddo
        tc(i)=0.
        tw(i)=0.
        enddo
        do i=1,KSZJ
        k(i,3)=0
        enddo
        kfmax=kfmaxi
        do i=1,50
        numb(i)=0
        enddo
        nztr=iop(12)
        natr=nat-iop(10)
        nzpr=iop(11)
        napr=nap-iop(9)
        natp=natr+napr
        nztp=nztr+nzpr
c       iop(9) (iop(10)) : wounded (i.e. participant) projectile (target)
c        nucleons.
c       iop(11) (iop(12)) : projectile (target) spectator protons.
c       natr (nztr) : the # of targer spectator nucleons (protons)
c       napr (nzpr) : the # of projectile spectator nucleons (protons)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        do i=1,250
        do j=1,5
        ppp(i,j)=0.
        ppt(i,j)=0.
        enddo
        enddo
        do i=1,4
        pp(i)=0.
        pt(i)=0.
        enddo
c       read the momenta of target & projetile spectator from fritiof 
        if(natr.eq.0 .and. napr.eq.0)goto 2000
        if(natr.ne.0 .and. napr.ne.0)then
        do i=1,4
        pp(i)=p(n-1,i)
        pt(i)=p(n,i)
        enddo
        goto 1200
        endif
        if(natr.ne.0 .and. napr.eq.0)then
        do i=1,4
        pt(i)=p(n,i)
        enddo
        goto 1200
        endif
        if(natr.eq.0 .and. napr.ne.0)then
        do i=1,4
        pp(i)=p(n,i)
        enddo
        endif
c        width of target & projectile spectator momentum 
c        distribution is set to be 40 Mev
1200    if(napr.ne.0)then
        pp3=pp(3)/napr
        temp=PARAM(12)
        pmaxp=3.*temp
        endif
        if(natr.ne.0)then
        pt3=pt(3)/natr
        temtc=PARAM(12)
        temtft=PARAM(12)
        temtf=temtft/(3.*0.93828)       
        pmaxt=3.*temtc
        endif
c       sample the momenta of projectile & target nucleons and keep four 
c        momentum conservation,respectively
        if(napr.eq.0)goto 1100
        call spemd1(napr,ppp,pp3,temp,pmaxp,pp,1,jj)
1100    if(natr.eq.0)goto 2000
        if(ifram.eq.1)call spemd1(natr,ppt,pt3,temtc,pmaxt,pt,2,jj)
        if(ifram.eq.0)call spemd2(natr,ppt,temtf,pt,2,jj)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       filter out those kind of particles wanted to study from fritiof
2000    continue
        iii=0
        jjj=0
        do i=1,kfmax
        kf=kfaco(i)
        do j=iii+1,n
        call ord(jjj,j,kf)
        enddo
        iii=jjj
        numb(i)=jjj
        enddo
        lll=numb(1)+1
        mmm=numb(kfmax)
        m9=numb(9)
        m10=numb(10)
        do i=mmm+1,n
        ishp(i)=0
        tp(i)=0.
        enddo
c       make the order of proton,neutrons,pba,nba,pi+,pi-,pi0,k-,k0-,sigma0,
c       sigma-,sigma+,sigma0ba,sigma-ba,sigma+ba,lamda,lamdaba,k0,k+,
c       cascade-,cascade-ba,cascade0,cascade0ba,omega-,
c       omega+,Delta-,Delta0,Delta+,Delta++,rho+,rho-,rho0. altogether 32 particles.
        do i=n,lll,-1
        ii=i+natp
        ishp(ii)=ishp(i)
        tp(ii)=tp(i)
        do jj=1,5
        k(ii,jj)=k(i,jj)
        p(ii,jj)=p(i,jj)
        v(ii,jj)=v(i,jj)
        enddo
        enddo
        do i=2,kfmax
        numb(i)=numb(i)+natp
        enddo
        mmm=numb(kfmax)

        nwhole=n+natp
c leaving space in event record for puting spector nueleons in 
        ipro=numb(1)
        do i=1,kfmax
        if(i.eq.1)then
        nnn=numb(i)
        goto 200
        endif
        nnn=numb(i)-numb(i-1)
        if(i.eq.2)nnn=nnn-natp
200     call arrse(i,natp,nnn,bp)
c       'arrse'is used to dis. particles with flavor i from fritiof  
c        in the tube dug by projectile in target
        enddo

        if(nat.eq.1 .and. nap.eq.1)goto 600
        if(nztr.gt.natr)nztr=natr
        if(nzpr.gt.napr)nzpr=napr
c       initiate the target spectator protons
        if(nztr.eq.0)goto 900
        call dsp(2,nztr,ipro,ppt)
c       dis. the target spectator protons outside the  
c        tube dug by projectile in target & inside the target.
        call arrnp(ipro,nztr,nat,nap,bp,rnt,rnp,nztr,natr,
     &  nzpr,napr,1)
c       initiate the projectile spectator protons
900     if(nzpr.eq.0)goto 1000
        nmm=ipro+nztr
        call dsp(1,nzpr,nmm,ppp)
c       dis. the projectile spectator protons outside the  
c        tube dug by projectile in target & inside the projectile.
        call arrnp(nmm,nzpr,nat,nap,bp,rnp,rnt,nztr,natr,
     &  nzpr,napr,2)

c       initiate the target spectator neutrons
1000    if(natr-nztr.eq.0)goto 700
        nmm=ipro+nztp
        call dsn(2,natr-nztr,nztr,nmm,ppt)
c       dis. the target spectator neotrons outside the  
c        tube dug by projectile in target & inside the target.
        call arrnp(nmm,natr-nztr,nat,nap,bp,rnt,rnp,nztr,natr,
     &  nzpr,napr,3)
c       initiate the projectile spectator neutrons
700     if(napr-nzpr.eq.0)goto 600
        nmm=ipro+nztp+(natr-nztr)
        call dsn(1,napr-nzpr,nzpr,nmm,ppp)
c       dis. the projectile spectator neotrons outside the  
c        tube dug by projectile in target & inside the projectile.
        call arrnp(nmm,napr-nzpr,nat,nap,bp,rnp,rnt,nztr,natr,
     &  nzpr,napr,4)
600     numb(1)=numb(1)+nztp
c*********************************************************************************
        if(abs(k(nwhole,2)).ge.10000.and.abs(k(nwhole-1,2)).
     &       ge.10000)then
        n=nwhole-2
        elseif(abs(k(nwhole,2)).ge.10000)then
        n=nwhole-1
        else
        endif
        nna=n
c*********************************************************************************

1300    return
        end
C******************************************************************************
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

        subroutine spemd1(np,pp,p3,tem,pmax,ps,ii,jj)
c       sample the momenta of spectator nucleons in moving nucleus and
c        keep four momentum conservation
c       np : the # of spectator nucleons
c       ps : the rest four momentum of nucleus
c       p3 = ps(3)/np
c       pp : four momentum of spectator nucleon
c       tem,pmax : the parameters in two dimension Gaussian distribution
        dimension pp(250,5),ps(4) 
        call tdgaus(tem,pmax,np-1,pp)
        px=0.
        py=0.
        pz=0.
        do i=1,np-1
        pp(i,3)=p3
        px=px+pp(i,1)
        py=py+pp(i,2)
        pz=pz+p3
        pp(i,4)=sqrt(0.880+pp(i,1)**2+pp(i,2)**2+p3**2)
        pp(i,5)=0.938
        enddo
        pp(np,1)=ps(1)-px
        pp(np,2)=ps(2)-py
        pp(np,3)=ps(3)-pz
        pp(np,4)=sqrt(0.880+pp(np,1)**2+pp(np,2)**2+pp(np,3)**2)
        pp(np,5)=0.938
200     format(4(1x,f8.3))
        call conse(np,pp,ps,ii,jj)
        return
        end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


      subroutine tdgaus(v,pmax,np,pp)
c.... 2-d Gaussian distribution with width v, i.e., e^(-p2/v)dp2, 0<p2<pmax
c.... set pmax < 0 if pmax should be infinity.
c.... np : the total # of particles wanted to sample their transverse momenta
        dimension pp(250,5)
        do 30 i=1,np
        p2 = 0
      if(v.le.1.e-8)return 
      if(pmax.lt.1.E-9)return
      if(pmax.lt.0)then
        a = 1.
        goto 10
        endif
        aa=-pmax/v
        if(aa.lt.-70)then
        a=1.
        goto 10
        endif
      a = 1. - exp(aa)
10    p2 = -v*log(max(1.e-20,1. - a*rlu(1)))
      if(p2.LT.0.)goto 10
        ps=sqrt(p2)
        fi=2.*3.1415926*rlu(1)
        pp(i,1)=ps*cos(fi)
        pp(i,2)=ps*sin(fi)
30      continue
      return
      end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        subroutine spemd2(np,pp,tem,ps,ii,jj)
c       sample the momenta of spectator nucleons in rest nucleus and
c        keep four momentum conservation
c       np : the # of spectator nucleons
c       ps : the rest four momentum of nucleus
c       pp : four momentum of spectator nucleon
c       tem : the parameter in Boltzmann distribution
        dimension pp(250,5),ps(4) 
        px=0.
        py=0.
        pz=0.
        do l=1,np-1
100     se1=rlu(1)
        se1=se1*se1
        se2=rlu(1)
        se2l=log(se2)
        udg=-2.718*se2*se2l
        if(se1.gt.udg)goto 100
        xf=-3.*0.938*tem*se2l
        pf=sqrt(xf)
        cta=2.*rlu(1)-1.
        sta=sqrt(1.-cta*cta)
        fi=2.*3.1415926*rlu(1)
        pp(l,1)=pf*sta*cos(fi)
        pp(l,2)=pf*sta*sin(fi)
        pp(l,3)=pf*cta
        pp(l,4)=sqrt(0.880+pp(l,1)**2+pp(l,2)**2+pp(l,3)**2)
        pp(l,5)=0.938
        px=px+pp(l,1)
        py=py+pp(l,2)
        pz=pz+pp(l,3)
        enddo
        pp(np,1)=ps(1)-px
        pp(np,2)=ps(2)-py
        pp(np,3)=ps(3)-pz
        pp(np,4)=sqrt(0.880+pp(np,1)**2+pp(np,2)**2+pp(np,3)**2)
        pp(np,5)=0.938
200     format(4(1x,f8.3))
        call conse(np,pp,ps,ii,jj)
        return
        end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        subroutine conse(np,pp,ps,ii,jj)
c       keep four momentum conservation
c       np : the # of spectator nucleons
c       ps : the rest four momentum of nucleus
c       pp : four momentum of spectator nucleon
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn
     c  ,rnt,rnp,rao,rou0
        COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        dimension pp(250,5),ps(4),ff(250),pxyz(3),arp(3)
        ps4=ps(4)
        do i=1,3
        pxyz(i)=0.
        enddo
        jj=0
100     es=0.
        do i=1,np
        es=es+pp(i,4)
        enddo
        fr=es/ps4
        if(abs(1.-fr) .le. dep)goto 200
        do i=1,np
        ppm=pp(i,4)/0.938
        ppf=ppm/fr
        ff(i)=sqrt(abs(ppf*ppf-1.)/(ppm*ppm-1.))
        do j=1,3
        ppp=ff(i)*pp(i,j)
        pp(i,j)=ppp
        pxyz(j)=pxyz(j)+ppp
        enddo
        enddo
        do i=1,3
        arp(i)=abs(1.-pxyz(i)/ps(i))
        pxyz(i)=pxyz(i)-ps(i)
        enddo
        if(abs(1.-fr).le.dep .and.arp(1).le.dep .and. arp(2).le.dep  
     c   .and. arp(3).le.dep) goto 200
        do i=1,3
        pxyz(i)=pxyz(i)/np
        enddo
        do i=1,np
        do j=1,3
        pp(i,j)=pp(i,j)-pxyz(j)
        enddo
        pp(i,4)=sqrt(0.880+pp(i,1)**2+pp(i,2)**2+pp(i,3)**2)
c       0.880 = 0.938*0.938
        enddo
        jj=jj+1
        if(jj.eq.4000)then
        write(MSTU(11),*)'infinitive loop may occur in subroutine
     &        conse(),which means four-momentum conservation to 
     &         accuracy needed is hard to be achieved,check value 
     &  PARAM(9)'
        return
        endif
        goto 100
200     return
        end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


        subroutine ord(ipi,j,kf)
c       make order for particles with flavor code kf.
c       j : the particle wanted to order
c       ipi : j-th particle should order after ipi
        parameter(mcludi=40000)
        parameter(KSZJ=40000)
        COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/wz/c17(mcludi,3),ishp(mcludi),tp(mcludi)
     c   ,coor(3)
        dimension pp(5)
        ik=k(j,2)
        if(ik.eq.kf)then
        ipi=ipi+1
        kk=k(ipi,2)
        k1=k(ipi,1)
        k3=k(ipi,3)
        do jj=1,5
        pp(jj)=p(ipi,jj)
        enddo
        k(ipi,2)=k(j,2)
        k(ipi,1)=k(j,1)
        k(ipi,3)=k(j,3)
        do jj=1,5
        p(ipi,jj)=p(j,jj)
        enddo
        ishp(ipi)=1
        tp(ipi)=0.
        k(j,2)=kk
        k(j,1)=k1
        k(j,3)=k3
        do jj=1,5
        p(j,jj)=pp(jj)
        enddo
        endif
        return
        end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        subroutine arrse(i,natp,nnn,bp)
c       It is to arrange particles with flavor i from Fritiof in the 
c         tube dug by projectile in target
c       nnn : the # of i-th particles from Fritiof wanted to arrange
        parameter (mcludi=40000,KSZJ=40000,KSZ1=30)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
        common/sa12/psa(5),ptai(5),clorenp,clorent
        common/sa5/kfmax,kfaco(50),numb(50),disbe(50,50)
        common/wz/c17(mcludi,3),ishp(mcludi),tp(mcludi)
     c   ,coor(3)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn
     c  ,rnt,rnp,rao,rou0
        b=bp/rnt
        kkk=0
        rnpt=rnp+rnt
        if(rnpt-bp.le.0.8)kkk=1
c kkk=1 and kk1=1 non-overlapping condition is abolished
        dis=disbe(i,i)
        il=0
        ll=0
        if(i.eq.1)goto 10
        ll=numb(i-1)
        if(i.eq.2)ll=ll+natp
        il=ll
10      do 53 j=1,nnn   
        il=il+1
        iii=0
        kk1=0
54      iii=iii+1
        if(((rnpt-bp).gt.0.8 .and. (rnpt-bp).le.1.6) .and. iii.gt.2000)
     &   kk1=1
        if(iii.eq.10000)then
       write(mstu(11),*)'difficult to arrange produced particles in 
     &  ubroutine arrse(),infinitive loop may occur'
        endif
        x=1.-2.*rlu(1)
        y=1.-2.*rlu(1)
        z=1.-2.*rlu(1)
        rr=sqrt(x*x+y*y+z*z)
        if(rr.gt.1) goto 54
        r1=rnt*sqrt(x*x+(b-y)*(b-y))
        if(r1.gt.rnp)goto 54

        c17(il,1)=x*rnt
        c17(il,2)=y*rnt
        c17(il,3)=z*rnt

        if(il.eq.ll+1) goto 52
        if(kkk.eq.1)goto 52
        if(kk1.eq.1)goto 52
c       non-overlapping demand between current particle and 
c        already distributed same kf particles.
        do j1=ll+1,il-1
        dcc=sqrt((c17(j1,1)-c17(il,1))**2+(c17(j1,2)-c17(il,2))**2
     c  +(c17(j1,3)-c17(il,3))**2)
        if(dcc.lt.dis) goto 54
        enddo
52      continue
        if(i.eq.1)goto 55
        if(kkk.eq.1) goto 55
        if(kk1.eq.1)goto 55
c       non-overlapping demand between current particle and 
c        already distributed different kf particles.
        do 100 jj=1,i-1
        diss=disbe(i,jj)
        if(jj.eq.1)then
        lll=1
        goto 200
        endif
        lll=numb(jj-1)+1
        if(jj.eq.2)lll=lll+natp
200     mmm=numb(jj)
        do 218 i11=lll,mmm
        dcc=sqrt((c17(il,1)-c17(i11,1))**2+(c17(il,2)-c17(i11,2))**2
     c  +(c17(il,3)-c17(i11,3))**2)
        if(dcc.lt.diss) goto 54
218     continue
100     continue
Ctai************************************
55      if(i.eq.1 .or. i.eq.2)then
        kil=k(il,2)
        tau(il)=t0*p(il,4)/pmas(lucomp(kil),1)
        goto 56 
        endif
Ctai************************************
        kil=k(il,2)
        tau(il)=t0*p(il,4)/pmas(lucomp(kil),1)
56      ishp(il)=1
        c17(il,3)=c17(il,3)/clorent
        tp(il)=0.
53      continue
        return
        end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


        subroutine dsp(ijk,nzr,nnu,pp)
c       initiate the spectator protons
c       nzr : the # of protons
c       nnu : the protons are ordered after nnu
        parameter (mcludi=40000)
        parameter(KSZJ=40000)
        common/wz/c17(mcludi,3),ishp(mcludi),tp(mcludi)
     c   ,coor(3)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        dimension pp(250,5)
        do i=1,nzr
        l=i+nnu
        k(l,2)=2212
c give K(L,3)=33 for spectator protons of the projectile and 
c K(L,3)=55 for spectator protons of the target.
        k(l,1)=1
        k(l,3)=33
        if(ijk.eq.2)k(l,3)=55
        p(l,1)=pp(i,1)
        p(l,2)=pp(i,2)
        p(l,3)=pp(i,3)
        p(l,4)=pp(i,4)
        p(l,5)=pp(i,5)
        ishp(l)=1
        tau(l)=0.
        tp(l)=0.
        enddo
        return
        end

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


        subroutine dsn(ijk,nur,nz,nnu,pp)
c       initiate the spectator neutrons
c       nur : the # of neutrons; nz : the # of corresponding protons
c       nnu : the neutrons are ordered after nnu
        parameter (mcludi=40000)
        parameter(KSZJ=40000)
        common/wz/c17(mcludi,3),ishp(mcludi),tp(mcludi)
     c   ,coor(3)
         COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        dimension pp(250,5)
        do i=1,nur
        j=i+nz
        l=i+nnu
        k(l,2)=2112
        k(l,1)=1
cc give K(L,3)=33 for spectator neutrons of the projectile and 
c K(L,3)=55 for spectator neutrons of the target.
        k(l,3)=33
        if(ijk.eq.2)k(l,3)=55
        p(l,1)=pp(j,1)
        p(l,2)=pp(j,2)
        p(l,3)=pp(j,3)
        p(l,4)=pp(j,4)
        p(l,5)=pp(j,5)
        ishp(l)=1
        tau(l)=0.
        tp(l)=0.
        enddo
300     return
        end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


        subroutine arrnp(nnu,naw,nat,nap,bp,rn1,rn2,nztr,natr,nzpr,napr,
     &   jjj)
c       arrange the spectator nucleons outside the tube
c        dug by projectile in target and inside target (projectile)
c       naw : the # of spectator nucleons wanted to arrange
c       nnu : the arranged nucleons should be ordered after nnu
c       for target (projectile) spectator rn1 = rnt (rnp)
c       for target (projectile) spectator rn2 = rnp (rnt)
c       y direction is directed from center of target to center of projectile
c       z is along beam direction
        parameter (mcludi=40000,KSZJ=40000,KSZ1=30)
        common/wz/c17(mcludi,3),ishp(mcludi),tp(mcludi)
     c   ,coor(3)
        COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
        common/sa4/tau(kszj),tlco(kszj,4)
          COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa5/kfmax,kfaco(50),numb(50),disbe(50,50)
        common/sa12/psa(5),ptai(5),clorenp,clorent
        dis=disbe(1,1)
        b=bp/rn1
c       rn1 : the unit of radius in this subroutine
        kkk=0
cKKK=1 non-overlapping condition is abolished
        rnbp=rn1+bp
        ll=nnu
        il=nnu
        do 59 j=1,naw
        il=il+1
        iii=0
54      iii=iii+1
        if(iii.eq.10000)then
        write(mstu(11),*)'difficult to arrange spectator nucleons in 
     &         subroutine arrnp(),infinitive loop may have occured'
        endif
        x=1.-2.*rlu(1)
        y=1.-2.*rlu(1)
        z=1.-2.*rlu(1)
        rr=sqrt(x*x+y*y+z*z)
        if(rr.gt.1) goto 54
        yy=y
        if(jjj.eq.2 .or. jjj.eq.4)yy=-y
        r1=rn1*sqrt(x*x+(b-yy)*(b-yy))
        if(bp.lt.1.e-5 .and. nat.eq.nap)goto 55
        if(rnbp.le.rn2)goto 55
c       it happens when projectile is smaller than target
c       at above case there should not be projectile spectator nucleon indeed
        if(rnbp.gt.rn2 .and.rnbp.le.rn2+0.8 .and. iii.gt.2000)then
        kkk=1
c luciae3.0s should the condition if(r1.lt.rn2)goto 54 be used here ??
        goto 55
        endif
        if(r1.lt.rn2)goto 54
55      c17(il,1)=x*rn1
        c17(il,2)=y*rn1
        c17(il,3)=z*rn1
        if(jjj.eq.2 .or. jjj.eq.4)c17(il,2)=c17(il,2)+b*rn1
        if(il.eq.ll+1) goto 52
        if(kkk.eq.1) goto 52
c       non-overlapping demand between spectator nucleon and already 
c        distributed particles.
        do j1=ll+1,il-1
        dcc=sqrt((c17(j1,1)-c17(il,1))**2+(c17(j1,2)-c17(il,2))**2
     c  +(c17(j1,3)-c17(il,3))**2)
        if(dcc.lt.dis) goto 54
        enddo
52      continue
c       non-overlapping demand between current specator particle and already 
c        distributed ,belonging to same nucleus, spectator particles.
500     if (jjj.eq.1 .or. jjj.eq. 2)goto 53
        if(jjj.eq.3)then
c       the current distributed spectator particles are target neutrons
        if(nztr.eq.0)goto 53 
        mmm=nnu-nzpr
        lll=mmm-nztr+1
        goto 300
        endif
c       the distributed current spectator particles are projectile neutrons
        if(nzpr.eq.0)goto 53 
        mmm=nnu-(natr-nztr)
        lll=mmm-nzpr+1
300     do 400 i11=lll,mmm      
        dcc=sqrt((c17(il,1)-c17(i11,1))**2+(c17(il,2)-c17(i11,2))**2
     c  +(c17(il,3)-c17(i11,3))**2)
        if(dcc.lt.dis) goto 54
400     continue
        
53      if(jjj.eq.1.or.jjj.eq.3)then
        c17(il,3)=c17(il,3)/clorent
        else
        c17(il,3)=c17(il,3)/clorenp
        endif
59      continue
        return
        end
C**************************************************************************
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        subroutine ctlcre(lc,tc,tw)
c       creat the initial collision list  
        parameter (KSZJ=40000,KSZ1=30)
        parameter(nsize=100000)
          COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
        COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa5/kfmax,kfaco(50),numb(50),disbe(50,50)
        common/ctllist/nctl,noinel(400),nctl0
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        time=0.
        nctl=1
        m1=numb(1)
        m2=numb(2)
        m3=numb(3)
        m4=numb(4)
        m7=numb(7)
        m19=numb(19)
        m32=numb(32)
        m25=numb(25)
        m29=numb(29)
        if(KFR(22).eq.1)then
        do 10 l=1,m32
        do 10 l1=1,m32
        if(l.le.m2 .and. l1.le.m2.and.l1.lt.l)goto 10
c no double counting for NN
        if((l1.gt.m2 .and. l1.le.m4).and.l.gt.m25)goto 10
c delta and rho do not scatter with Nbar        
        if(l.le.m4.and. l1.gt.m2)goto 10
c NNbr and NN are needed
c       if((l1.gt.m2 .and. l1.le.m4) .and. l.eq.l1)goto 10
        if((l.gt.m7.and.l1.gt.m25).or.
     &  ((l.gt.m4.and.l.le.m7).and.l1.gt.m29))goto 10
c delta and rho do not collid with the strange particles and rho do not
c collid with pi
        if((k(l1,3).eq.33.or.k(l1,3).eq.55).and.
     &         (k(l,3).eq.33.or.k(l,3).eq.55))goto 10
C       only collisions between produced nucleons and spectator nucleons 
c       are considered

c only pi scatter with detal 
        if(l.gt.m25.and.l1.gt.m25)goto 10
c delta and rho do not scatter with each other

        iflag=0
        call rsfilt(l,l1,iflag,0)
        if(iflag.eq.0)goto 10
        if(nctl.gt.nsize)then
        write(mstu(11),*)'size of array "nsize" needs to be extended
     &         error is serious,stop running'
         stop 
        endif
        tc(nctl)=0.
        call tcolij(l,l1,time,nctl,lc,tc,tw)
        if(tc(nctl).gt.1.0e-7)nctl=nctl+1
10      continue

        else
        do 20 l=1,m25
        do 20 l1=1,m2 
        iflag=0
cfilter pairs that are of interest,this choice is put here for the case that we
c are only interested in certain rescattering processes
        call rsfilt(l,l1,iflag,0)
        if(iflag.eq.0)goto 20
        if(nctl.gt.nsize) stop 100000
        tc(nctl)=0.
        call tcolij(l,l1,time,nctl,lc,tc,tw)
        if(tc(nctl).gt.1.0e-7)nctl=nctl+1
20      continue
        endif
        if(tc(nctl).le.1.0e-7)nctl=nctl-1
        do i=nctl+1,nsize
        do m=1,5
        lc(i,m)=0
        enddo
        tc(i)=0.
        tw(i)=0.
        enddo
        return
        end
C*******************************************************************
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


        subroutine rsfilt(l,l1,iflag,id)
c        subroutine rsfilt plays the role of first range filter 
c         and guarantees that the collision list is composed due to the
c         entrance channels of selected inelastic reactions
c        subroutine intdis plays the role of second range filter
c        both of them filter those collision pairs that are not of our concern
        parameter (KSZJ=40000,KSZ1=30,mcludi=40000)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa5/kfmax,kfaco(50),numb(50),disbe(50,50)
        common/wz/c17(mcludi,3),ishp(mcludi),tp(mcludi)
     c   ,coor(3)
        m2=numb(2)
        m4=numb(4)
        kl=k(l,2)
        kl1=k(l1,2)
        if(KFR(22).eq.1)then
        if(l.eq.l1) goto 10
        if(ishp(l).eq.0.or.ishp(l1).eq.0) goto 10
c       constraints on the direct reactions
        if(kl.eq.211 .and. (kl1.eq.-211 .or. kl1.eq.111 .or.
     &   abs(kl1).eq.2212 .or. abs(kl1).eq.2112 .or. kl1.eq.
     &   3112 .or. kl1.eq.-3122 .or. kl1.eq.-3222 .or. kl1
     &   .eq.-3212 .or. kl1.eq.3212 .or. kl1.eq.3122 .or. kl1.eq.3312
     &   .or. kl1.eq.-3322))goto 11
        if(kl1.eq.211 .and. (kl.eq.-211 .or. kl.eq.111 .or.
     &   abs(kl).eq.2212 .or. abs(kl).eq.2112 .or. kl.eq.
     &   3112 .or. kl.eq.-3122 .or. kl.eq.-3222 .or. kl
     &   .eq.-3212 .or. kl.eq.3212 .or. kl.eq.3122 .or. kl.eq.3312
     &   .or. kl.eq.-3322))goto 11
        if(kl.eq.-211 .and. (kl1.eq.111 .or.
     &   abs(kl1).eq.2212 .or. abs(kl1).eq.2112 .or. kl1.eq.
     &   -3112 .or. kl1.eq.3122 .or. kl1.eq.3222 .or. kl1
     &   .eq.3212 .or. kl1.eq.-3212 .or. kl1.eq.-3122 .or. kl1.eq.
     &   -3312 .or. kl1.eq.3322))goto 11
        if(kl1.eq.-211 .and. (kl.eq.111 .or.
     &   abs(kl).eq.2212 .or. abs(kl).eq.2112 .or. kl.eq.
     &   -3112 .or. kl.eq.3122 .or. kl.eq.3222 .or. kl
     &   .eq.3212 .or. kl.eq.-3212 .or. kl.eq.-3122 .or. kl.eq.
     &   -3312 .or. kl.eq.3322))goto 11
        if(kl.eq.111 .and. (kl1.eq.111 .or. abs(kl1).eq.2212
     &   .or. abs(kl1).eq.2112 .or. abs(kl1).eq.3122 .or.
     &   abs(kl1).eq.3112 .or. abs(kl1).eq.3212 .or. abs(kl1).eq.3222
     &   .or. abs(kl1).eq.3312 .or. abs(kl1).eq.3322))goto 11
        if(kl1.eq.111 .and. (kl.eq.111 .or. abs(kl).eq.2212
     &   .or. abs(kl).eq.2112 .or. abs(kl).eq.3122 .or.
     &   abs(kl).eq.3112 .or. abs(kl).eq.3212 .or. abs(kl).eq.3222
     &   .or. abs(kl).eq.3312 .or. abs(kl).eq.3322))goto 11
        if(kl.eq.321 .and. (kl1.eq.-2212 .or. kl1.eq.-2112 .or. 
     &   kl1.eq.-3122 .or. kl1.eq.-3222 .or. kl1.eq.-3112
     &   .or. kl1.eq.-3212 .or.kl1.eq.-3312 .or.kl1.eq.-3322))goto 11
        if(kl1.eq.321 .and. (kl.eq.-2212 .or. kl.eq.-2112 .or. 
     &   kl.eq.-3122 .or. kl.eq.-3222 .or. kl.eq.-3112
     &   .or. kl.eq.-3212 .or.kl.eq.-3312 .or.kl.eq.-3322))goto 11
        if(kl.eq.-321 .and. (kl1.eq.2212 .or. kl1.eq.2112 .or. 
     &   kl1.eq.3122 .or. kl1.eq.3222 .or. kl1.eq.3112
     &   .or. kl1.eq.3212 .or.kl1.eq.3312 .or. kl1.eq.3322))goto 11
        if(kl1.eq.-321 .and. (kl.eq.2212 .or. kl.eq.2112 .or. 
     &   kl.eq.3122 .or. kl.eq.3222 .or. kl.eq.3112
     &   .or. kl.eq.3212 .or.kl.eq.3312 .or. kl.eq.3322))goto 11
        if(kl.eq.311 .and. (kl1.eq.-2212 .or. kl1.eq.-2112 .or. kl1.eq. 
     &   -3122 .or. kl1.eq.-3212 .or. kl1.eq.-3112 .or.kl1.eq.-3222
     &   .or. kl1.eq.-3312 .or. kl1.eq.-3322))goto 11
        if(kl1.eq.311 .and. (kl.eq.-2212 .or. kl.eq.-2112 .or. kl.eq. 
     &   -3122 .or. kl.eq.-3212 .or. kl.eq.-3112 .or.kl.eq.-3222
     &   .or. kl.eq.-3312 .or. kl.eq.-3322))goto 11
        if(kl.eq.-311 .and. (kl1.eq.2212 .or. kl1.eq.2112 .or. kl1.eq. 
     &   3122 .or. kl1.eq.3212 .or. kl1.eq.3112 .or.kl1.eq.3222
     &   .or.kl1.eq.3312 .or.kl1.eq.3322))goto 11
        if(kl1.eq.-311 .and. (kl.eq.2212 .or. kl.eq.2112 .or. kl.eq. 
     &   3122 .or. kl.eq.3212 .or. kl.eq.3112 .or.kl.eq.3222
     &   .or.kl.eq.3312 .or.kl.eq.3322))goto 11
        if(kl.eq.2112.and.(kl1.eq.2112.or.kl1.eq.2212))goto 11
        if(kl.eq.2212.and.(kl1.eq.2112.or.kl1.eq.2212))goto 11
c       constraints on the annihilation reactions
        if(kl.eq.-3212 .and. (kl1.eq.2212 .or. kl1.eq.2112))goto 11
        if(kl.eq.-3122 .and. (kl1.eq.2212.or. kl1.eq.2112))goto 11
        if(kl1.eq.-3212 .and. (kl.eq.2212 .or. kl.eq.2112))goto 11
        if(kl1.eq.-3122 .and. (kl.eq.2212.or. kl.eq.2112))goto 11
        if(kl.eq.-2212 .and. (kl1.eq.2212 .or. kl1.eq.2112))goto 11
        if(kl1.eq.-2212 .and. (kl.eq.2212 .or. kl.eq.2112))goto 11
        if(kl.eq.-2112 .and. (kl1.eq.2212.or. kl1.eq.2112))goto 11
        if(kl1.eq.-2112 .and. (kl.eq.2212.or. kl.eq.2112))goto 11
c       constraints on the reverse reactions
        if(kl.eq.211 .and. (kl1.eq.3112 .or. abs(kl1).eq.3122 .or. kl1
     &   .eq.-3222 .or. abs(kl1).eq.3212 .or. abs(kl1).eq.3312 .or.
     &   abs(kl1).eq.3322 .or. abs(kl1).eq.3334.or.kl1.eq.1114
     &    .or.kl1.eq.2114.or.kl1.eq.2214))goto 11
        if(kl1.eq.211 .and. (kl.eq.3112 .or. abs(kl).eq.3122 .or. kl
     &   .eq.-3222 .or. abs(kl).eq.3212 .or. abs(kl).eq.3312 .or.
     &   abs(kl).eq.3322 .or. abs(kl).eq.3334.or.kl.eq.1114
     &    .or.kl.eq.2114.or.kl.eq.2214))goto 11
        if(kl.eq.-211 .and. (kl1.eq.3222 .or. abs(kl1).eq.3122 .or. kl1
     &   .eq.-3112 .or. abs(kl1).eq.3212 .or. abs(kl1).eq.3312 .or.
     &   abs(kl1).eq.3322 .or. abs(kl1).eq.3334.or.
     &   kl1.eq.2224.or.kl1.eq.2114.or.kl1.eq.2214))goto 11
        if(kl1.eq.-211 .and. (kl.eq.3222 .or. abs(kl).eq.3122 .or. kl
     &   .eq.-3112 .or. abs(kl).eq.3212 .or. abs(kl).eq.3312 .or.
     &   abs(kl).eq.3322 .or. abs(kl).eq.3334.or.
     &   kl.eq.2224.or.kl.eq.2114.or.kl.eq.2214))goto 11
        if(kl.eq.111 .and. (abs(kl1).eq.3112 .or. abs(kl1).eq.3122 .or. 
     &   abs(kl1).eq.3222 .or. abs(kl1).eq.3212 .or. abs(kl1).eq.3312
     &   .or. abs(kl1).eq.3322 .or. abs(kl1).eq.3334.or.kl1.eq.2224.
     &    or.kl1.eq.2114.or.kl1.eq.2214.or.kl1.eq.1114))goto 11
        if(kl1.eq.111 .and. (abs(kl).eq.3112 .or. abs(kl).eq.3122 .or. 
     &   abs(kl).eq.3222 .or. abs(kl).eq.3212 .or. abs(kl).eq.3312
     &   .or. abs(kl).eq.3322 .or. abs(kl).eq.3334.or.kl.eq.2224.
     &    or.kl.eq.2114.or.kl.eq.2214.or.kl.eq.1114))goto 11
        if(kl.eq.321 .and. (kl1.eq.-321 .or. kl1.eq.-311 .or. kl1.eq.
     &   3222 .or. kl1.eq.3212 .or. kl1.eq.3112 .or. kl1.eq.3122 .or. 
     &   kl1.eq.3312 .or. kl1.eq.3322 .or. kl1.eq.3334))goto 11
        if(kl1.eq.321 .and. (kl.eq.-321 .or. kl.eq.-311 .or. kl.eq.
     &   3222 .or. kl.eq.3212 .or. kl.eq.3112 .or. kl.eq.3122 .or. 
     &   kl.eq.3312 .or. kl.eq.3322 .or. kl.eq.3334))goto 11
        if(kl.eq.-321 .and. (kl1.eq.311 .or. kl1.eq.-3222 .or. kl1.eq.
     &   -3212 .or. kl1.eq.-3112 .or. kl1.eq.-3122 .or. kl1.eq.-3312 
     &  .or. kl1.eq.-3322 .or. kl1.eq.-3334))goto 11
        if(kl1.eq.-321 .and. (kl.eq.311 .or. kl.eq.-3222 .or. kl.eq.
     &   -3212 .or. kl.eq.-3112 .or. kl.eq.-3122 .or. kl.eq.-3312
     &   .or. kl.eq.-3322 .or. kl.eq.-3334))goto 11
        if(kl.eq.311 .and. (kl1.eq.-311 .or. kl1.eq.3222 .or. kl1.eq.
     &   3212 .or. kl1.eq.3112 .or. kl1.eq.3122 .or. 
     &   kl1.eq.3312 .or. kl1.eq.3322 .or. kl1.eq.3334))goto 11
        if(kl1.eq.311 .and. (kl.eq.-311 .or. kl.eq.3222 .or. kl.eq.
     &   3212 .or. kl.eq.3112 .or. kl.eq.3122 .or. 
     &   kl.eq.3312 .or. kl.eq.3322 .or. kl.eq.3334))goto 11
        if(kl.eq.-311 .and. (kl1.eq.-3222 .or. kl1.eq.-3212 .or. kl1
     &   .eq.-3112 .or. kl1.eq.-3122 .or. kl1.eq.-3312 .or.
     &   kl1.eq.-3322 .or. kl1.eq.-3334))goto 11
        if(kl1.eq.-311 .and. (kl.eq.-3222 .or. kl.eq.-3212 .or. kl
     &   .eq.-3112 .or. kl.eq.-3122 .or. kl.eq.-3312 .or.
     &   kl.eq.-3322 .or. kl.eq.-3334))goto 11
        if(kl1.eq.2212.and.(kl.eq.1114.or.kl.eq.2114.or.
     &       kl.eq.2214.or.abs(kl).eq.213.or.kl.eq.113))goto 11
        if(kl1.eq.2112.and.(kl.eq.2224.or.kl.eq.2114.or.
     &       kl.eq.2214.or.abs(kl).eq.213.or.kl.eq.113))goto 11
        if(kl.eq.2212.and.(kl1.eq.1114.or.kl1.eq.2114.or.
     &       kl1.eq.2214.or.abs(kl1).eq.213.or.kl1.eq.113))goto 11
        if(kl.eq.2112.and.(kl1.eq.2224.or.kl1.eq.2114.or.
     &       kl1.eq.2214.or.abs(kl1).eq.213.or.kl1.eq.113))goto 11
        goto 10
11      iflag=1
10      continue
        else
cThis part related to KFR(22) is not used now
        if(l.eq.l1)goto 20
        if (l.gt.m2 .and. l.le.m4)goto 20
        k3l=k(l,3)
        k3l1=k(l1,3)
        if(l.gt.m4 .and. (k3l1.ne.55.or. k3l1.ne. 44))goto 20
        if(l.le.m2 .and. k3l.eq.k3l1)goto 20
        if(l.le.m2 .and. (k3l.eq.55 .and. k3l1.eq.44))goto 20
        if(l.le.m2 .and. (k3l1.eq.55 .and. k3l.eq.44))goto 20
        iflag=1
20      continue
        endif
        return
        end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        subroutine tcolij(l,l1,time,icp,lc,tc,tw)
c       It is used to calculate the collision time & fill
c        up lc(i,1-2),tc(i). 
        parameter (mcludi=40000)
        parameter (KSZJ=40000,KSZ1=30)
        parameter(nsize=100000)
        COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
        common/sa10/csnn,cspin,cskn,cspipi,rcsit,ifram
        COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        common/sa12/psa(5),ptai(5),clorenp,clorent
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn
     c  ,rnt,rnp,rao,rou0
        common/wz/c17(mcludi,3),ishp(mcludi),tp(mcludi)
     c   ,coor(3)
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
        double precision b(3),bta
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        dimension dr(3),db(3),pi(4),pj(4),vi(3),vj(3)
        dimension ri(4),rj(4),rfi(4),rfj(4)
         pel=p(l,4)
        pel1=p(l1,4)
        pi(4)=pel
        pj(4)=pel1
        do i=1,3
        pi(i)=p(l,i)
        pj(i)=p(l1,i)
        b(i)=(pi(i)+pj(i))/(pi(4)+pj(4))
        enddo
        ilo=0
        call lorntz(ilo,b,pi,pj)
c       perform Lorentz transf. to CMS frame for momentum.
        bta=dsqrt(b(1)**2+b(2)**2+b(3)**2)
cif boost is too violent,put particles on mass shell by hand.
        if(bta.gt.0.99999d+0)then
        bmi=pmas(lucomp(k(l,2)),1)
        bmj=pmas(lucomp(k(l1,2)),1)
        pi(4)=sqrt(bmi**2+pi(1)**2+pi(2)**2+pi(3)**2+pi(4)**2)
        pj(4)=sqrt(bmj**2+pj(1)**2+pj(2)**2+pj(3)**2+pj(4)**2)
        endif
        ss=pi(4)+pj(4)
        IF(KFR(26).EQ.0)THEN    
c do not pair into the collision list if the threshold is too small.
        if(((abs(k(l,2)).eq.211.or.k(l,2).eq.111).and.
     &  (abs(k(l1,2)).eq.211.or.k(l1,2).eq.111)).and.ss.le.
     &  2.*pmas(lucomp(321),1))goto 10
        if(((abs(k(l,2)).eq.211.or.k(l,2).eq.111).and.
     &  (k(l1,2).eq.2112.or.k(l1,2).eq.2212)).and.ss.le.
     &  (pmas(lucomp(111),1)+1.232))goto 10
        if(((abs(k(l1,2)).eq.211.or.k(l1,2).eq.111).and.
     &  (k(l,2).eq.2112.or.k(l,2).eq.2212)).and.ss.le.
     &  (pmas(lucomp(111),1)+1.232))goto 10
        ELSE
        ENDIF

        do i=1,3
        ri(i)=c17(l,i)
        rj(i)=c17(l1,i)
        enddo
        ri(4)=time
        rj(4)=time
        call lorntz(ilo,b,ri,rj)
c       perform Lorentz transf. to CMS frame for coordinate.
        rb=0.
        bb=0.
        rr=0.
        kflag=0
        rtai=0.
        do ik=1,3
        vi(ik)=pi(ik)/pi(4)
        vj(ik)=pj(ik)/pj(4)
        enddo
        do i=1,3
        rfi(i)=c17(l,i)+(tau(l)-time)*(p(l,i)/p(l,4))
        rfj(i)=c17(l1,i)+(tau(l1)-time)*(p(l1,i)/p(l1,4))
        enddo
        rfi(4)=tau(l)
        rfj(4)=tau(l1)
        call lorntz(ilo,b,rfi,rfj)
c       gamli=p(l,4)/p(l,5)
c       gamlj=p(l1,4)/p(l1,5)
        ctaui=rfi(4)
        ctauj=rfj(4)
        tcol=ctaui
        if(ctaui.lt.ctauj)tcol=ctauj
        do ik=1,3
        db(ik)=(vi(ik)-vj(ik))*tcol
        dr(ik)=ri(ik)-rj(ik)-(vi(ik)*ri(4)-vj(ik)*rj(4))+db(ik)
        rtai=rtai+dr(ik)*dr(ik) 
        enddo
        dot=0.
        do ik=1,3
        dot=dr(ik)*pi(ik)+dot
        enddo
c       dot=-1
        if(dot.ge.0.)then
        kflag=1
        if(tcol.le.ri(4) )goto 10
        if(tcol.le.rj(4) )goto 10
        else
        rtai=0.
        do ik=1,3
        dr(ik)=ri(ik)-rj(ik)-(vi(ik)*ri(4)-vj(ik)*rj(4))
        db(ik)=vi(ik)-vj(ik)
        rb=rb+dr(ik)*db(ik)
        bb=bb+db(ik)*db(ik)
        rr=rr+dr(ik)*dr(ik)
        enddo
        if(bb .le. 1.e-10)goto 10
        tcol=0.-rb/bb
        if(tcol-ri(4) .le. 0.0)goto 10
        if(tcol-rj(4) .le. 0.0)goto 10
        if(tcol-ctaui .le. 0.)goto 10
        if(tcol-ctauj .le. 0.)goto 10   
        do ik=1,3
        dr(ik)=ri(ik)-rj(ik)-(vi(ik)*ri(4)-vj(ik)*rj(4))+tcol*db(ik)
        rtai=rtai+dr(ik)*dr(ik)
        enddo
c       for collision occurs,time must one step ahead
cwhen collision happens,particles should already be produced    
        sg1=rr+tcol*rb
        endif
        sg=rtai
c       if(sg.lt.0.)then
c       dmin=0.
c       tcol=-rr/rb
c       goto 20
c       endif
        dmin=sqrt(sg)
20      call intdis(l,l1,ss,rsig)
c       'intdis' : to calculate the interaction distance between 
c        particles l & l1.
        if(dmin.gt.rsig)goto 10
c       distance between the two particles should be smaller than rsig
        do ik=1,3
        ri(ik)=ri(ik)+vi(ik)*(tcol-ri(4))
        rj(ik)=rj(ik)+vj(ik)*(tcol-rj(4))
        enddo
c       move along Newton trajectory in CMS
        ri(4)=tcol
        rj(4)=tcol
        ilo=1
        call lorntz(ilo,b,ri,rj)
c       transform back to Lab.
        tcol1=ri(4)
        tcol2=rj(4)
        if(kflag.eq.0)then
        if(tcol1-tau(l).lt.0.) goto 10
        if(tcol2-tau(l1).lt.0.) goto 10
        else
        if(tcol1-tau(l).lt.-1.E-4) goto 10
        if(tcol2-tau(l1).lt.-1.E-4) goto 10
        endif
        if(ri(4).gt.rj(4)) ri(4)=rj(4)
        tcol=ri(4)
        if(tcol.le.time)goto 10
c       collision happens in the further
c       if(ifram.eq.0)coor(3)=coor(3)+rnt
        do i=1,3
        ri(i)=c17(l,i)+p(l,i)*(tcol-time)/pel-coor(i)
        rj(i)=c17(l1,i)+p(l1,i)*(tcol-time)/pel1-coor(i)
        enddo
        rri=sqrt(ri(1)*ri(1)+ri(2)*ri(2)+ri(3)*ri(3))
        rrj=sqrt(rj(1)*rj(1)+rj(2)*rj(2)+rj(3)*rj(3))
c the  rnt in rao*max(rnt,rnp)+rnt is due to the fact that
C we could not know the postion of the mass-center in the future.
        rrr=rao*rnt
c       if(ifram.eq.0)rrr=rao*max(rnt,rnp)+rnt
c       if(abs(k(l,2)).gt.1000.or.abs(k(l1,2)).gt.1000)rrr=1.E+10*rrr
c we do not apply any restriction to the interacting region in this version.
        rrr=1000000.0*rrr
c       if(rri.gt.rrr)goto 10
c       if(rrj.gt.rrr)goto 10
C       particles under considerarion are still within interaction region
C       when the collision happens
        tc(icp)=tcol
        lc(icp,1)=l
        lc(icp,2)=l1
10      return
        end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


        subroutine intdis(l,l1,ss,rsig)
c       to calculate the interaction distance between particles l 
c        and l1.
c       It plays also the role of second range filter
c       Only the collisions of (pion)(pion),(pion)n (and
c        inverse),(pion)y,kn and ky are considered here
c       y : refers to lambda or sigma
        parameter (KSZJ=40000)
                COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn
     c  ,rnt,rnp,rao,rou0
        rsig=0.
        kl=k(l,2)
        kl1=k(l1,2)
        if(abs(kl).eq.2212 .or. abs(kl).eq.2112)idpl=1
        if(abs(kl).eq.211 .or. kl.eq.111)idpl=3
        if(abs(kl).eq.321 .or. abs(kl).eq.311)idpl=4
        if(abs(kl).eq.3212 .or. abs(kl).eq.3112 .or. abs(kl).eq.3222
     c   .or. abs(kl).eq.3122 .or. abs(kl).eq.3312 .or. abs(kl).eq.
     c   3322 .or. abs(kl).eq.3334)idpl=5
        if(abs(kl).eq.213 .or. kl.eq.113)idpl=6
        if(kl.eq.1114 .or. kl.eq.2114.or.kl.eq.2214 .or. kl.eq.2224)
     &       idpl=7

        if(abs(kl1).eq.2212 .or. abs(kl1).eq.2112)idpl1=1
        if(abs(kl1).eq.211 .or. kl1.eq.111)idpl1=3
        if(abs(kl1).eq.321 .or. abs(kl1).eq.311)idpl1=4
        if(abs(kl1).eq.3212 .or. abs(kl1).eq.3112 .or. abs(kl1)
     c  .eq.3222 .or. abs(kl1).eq.3122 .or. abs(kl1).eq.3312
     c  .or. abs(kl1).eq.3322 .or. abs(kl1).eq.3334)idpl1=5
        if(abs(kl1).eq.213 .or. kl1.eq.113)idpl1=6
        if(kl1.eq.1114 .or. kl1.eq.2114.or.kl1.eq.2214 .or. kl1.eq.2224)
     c   idpl1=7

        if(idpl.eq.1 .and. idpl1.eq.1)rsig=ecsnn
        if(idpl.eq.3 .and. idpl1.eq.3)rsig=edipi
        if(idpl.eq.1 .and. idpl1.eq.3)rsig=epin
        if(idpl.eq.3 .and. idpl1.eq.1)rsig=epin
        if(idpl.eq.3 .and. idpl1.eq.5)rsig=epin
        if(idpl.eq.5 .and. idpl1.eq.3)rsig=epin
c       assume the total cross section of (pion)y,((pion)cascade) and  
c        ((pion)omiga) = (pion)n
        if(idpl.eq.4 .and. (idpl1.eq.1 .or. idpl1.eq.5))rsig=ekn
        if((idpl.eq.1 .or. idpl.eq.5) .and. idpl1.eq.4)rsig=ekn
        if(idpl.eq.4 .and. idpl1.eq.4)rsig=edipi

c       assume the total cross section of ky (k cascade) and (k omiga)
c        = kn
        if(idpl.eq.1 .and. idpl1.eq.6)rsig=epin
        if(idpl.eq.1 .and. idpl1.eq.7)rsig=ecsnn
        if(idpl.eq.3 .and. idpl1.eq.7)rsig=epin
        if(idpl1.eq.1 .and. idpl.eq.6)rsig=epin
        if(idpl1.eq.1 .and. idpl.eq.7)rsig=ecsnn
        if(idpl1.eq.3 .and. idpl.eq.7)rsig=epin
        if(idpl.eq.1 .and. idpl1.eq.5)rsig=ecsnn
        if(idpl.eq.5 .and. idpl1.eq.1)rsig=ecsnn
        return
        end


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        subroutine his(t1,lc,tc,tw,istop)
c       classical Newton motion
        parameter (mcludi=40000)
        parameter (KSZJ=40000)
        parameter(nsize=100000)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        common/sa12/psa(5),ptai(5),clorenp,clorent
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn
     c  ,rnt,rnp,rao,rou0
        common/ctllist/nctl,noinel(400),nctl0
        common/wz/c17(mcludi,3),ishp(mcludi),tp(mcludi)
     c   ,coor(3)
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        istop=1
        in=0
        do 200 i=1,n
c       if(t1.le.tau(i))goto 100
C       do move particles which have not produced
        if(ishp(i).eq.1) goto 10
        pp4=p(i,4)
        do j=1,3
        vp=p(i,j)/pp4
        c17(i,j)=c17(i,j)+vp*(t1-tp(i)) 
        enddo
        in=in+1
        goto 100
10      aa=0.
        pp4=p(i,4)
c       due to the fast speed of bayons, we could not use a limited interaction
c       region
        r0=rao*rnt
c we do not use any restriction to the interacting region in this version
        r0=1000000.0*r0
c       if(abs(k(i,2)).gt.1000)r0=1.E+10*r0
        do j=1,3
        vp=p(i,j)/pp4
        c17(i,j)=c17(i,j)+vp*(t1-tp(i))
        aa=aa+(c17(i,j)-coor(j))**2
        enddo
        aa=sqrt(aa)
        tp(i)=t1
c release the restriction to the interacting region
        r0=aa+100.
        if(aa.lt.r0) goto 100
cif  freeze-out occurs deduct the distance between the last collision and now
        ishp(i)=0
        do il=1,nctl
        if(lc(il,1).eq.i.or.lc(il,2).eq.i) tc(il)=0.
        enddo
100     continue
c       tp(i)=t1
200     continue
        if(in.eq.n) return
        istop=0
        return
        end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


        subroutine find(icp,tcp,lc,tc,tw)
        parameter(nsize=100000)
c       find out the binary collision with minimum collision time
        common/ctllist/nctl,noinel(400),nctl0
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        icp=0
        tcp=20000.
        do i=1,nctl
        if(tc(i).le.1.0e-7) goto 241
        if(tcp.lt.tc(i))  goto 241
        icp=i
        tcp=tc(i)
241     continue
        enddo
        return
        end


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine coinel(l,l1,ss,b,pi,pj,icp,pii,pjj,lc,tc,tw,winel,
     &   ik1,ik2)

c       treat the collision process.
c       the inelastic processes considered are:

c        strangeness production reactions:
c        1. pion+ + pion- to k+ + k-
c        2. pion+ + pion- to k0 + k0-
c        3. pion+ + pion0 to k+ + k0-
c        4. pion- + pion0 to k- + k0
c        5. pion0 + pion0 to k+ + k-
c        6. pion0 + pion0 to k0 + k0-

c        7. pion+ + p to k+ + sigma+
c        8. pion+ + n to k+ + sigma0
c        9. pion+ + n to k+ + lambda
c        10. pion+ + n to k0 + sigma+
c        11. pion- + p to k+ + sigma-
c        12. pion- + p to k0 + lambda
c        13. pion- + p to k0 + sigma0
c        14. pion- + n to k0 + sigma-
c        15. pion0 + p to k+ + sigma0
c        16. pion0 + p to k+ + lambda
c        17. pion0 + p to k0 + sigma+
c        18. pion0 + n to k+ + sigma-
c        19. pion0 + n to k0 + lambda
c        20. pion0 + n to k0 + sigma0

c        21. pion+ + pba to k0- + lambdaba
c        22. pion+ + pba to k0- + sigma0ba
c        23. pion+ + pba to k- + sigma-ba
c        24. pion+ + nba to k0- + sigma-ba
c        25. pion- + pba to k- + sigma+ba
c        26. pion- + nba to k- + lambdaba
c        27. pion- + nba to k- + sigma0ba
c        28. pion- + nba to k0- + sigma+ba
c        29. pion0 + pba to k- + lambdaba
c        30. pion0 + pba to k- + sigma0ba
c        31. pion0 + pba to k0- + sigma+ba
c        32. pion0 + nba to k- + sigma-ba
c        33. pion0 + nba to k0- + lambdaba
c        34. pion0 + nba to k0- + sigma0ba

c        35. pion+ + sigma- to k+ + cascade-
c        36. pion- + lambda to k0 + cascade- 
c        37. pion- + sigma+ to k+ + cascade- 
c        38. pion- + sigma0 to k0 + cascade- 
c        39. pion0 + lambda to k+ + cascade- 
c        40. pion0 + sigma- to k0 + cascade- 
c        41. pion0 + sigma0 to k+ + cascade- 
c        42. pion+ + lambdaba to k0- + cascade-ba 
c        43. pion+ + sigma+ba to k- + cascade-ba 
c        44. pion+ + sigma0ba to k0- + cascade-ba 
c        45. pion- + sigma-ba to k- + cascade-ba 
c        46. pion0 + lambdaba to k- + cascade-ba 
c        47. pion0 + sigma-ba to k0- + cascade-ba 
c        48. pion0 + sigma0ba to k- + cascade-ba
 
c        strangeness exchange reactions:
c        49. k- + p to pion0 + lambda
c        50. k- + p to pion0 + sigma0
c        51. k- + p to pion- + sigma+
c        52. k- + p to pion+ + sigma-
c        53. k- + p to k+ + cascade- (strangeness production)
c        54. k- + n to pion- + sigma0
c        55. k- + n to pion- + lambda
c        56. k- + n to pion0 + sigma-
c        57. k- + n to k0 + cascade- (strangeness production)

c        58. k0- + p to pion+ + lambda
c        59. k0- + p to pion+ + sigma0
c        60. k0- + p to pion0 + sigma+
c        61. k0- + n to pion+ + sigma-
c        62. k0- + n to pion- + sigma+
c        63. k0- + n to pion0 + sigma0
c        64. k0- + n to pion0 + lambda
c        65. k0- + n to k+ + cascade- (strangeness production)

c        66. k+ + p- to pion0 + lambda-
c        67. k+ + p- to pion0 + sigma0-
c        68. k+ + p- to pion+ + sigma+ba
c        69. k+ + p- to pion- + sigma-ba
c        70. k+ + p- to k- + cascade-ba (strangeness production)

c        71. k+ + n- to pion+ + sigma0-
c        72. k+ + n- to pion+ + lambda-
c        73. k+ + n- to pion0 + sigma-ba
c        74. k+ + n- to k0- + cascade-ba(strangeness production)

c        75. k0 + p- to pion- + lambda-
c        76. k0 + p- to pion- + sigma0-
c        77. k0 + p- to pion0 + sigma+ba
c        78. k0 + n- to pion- + sigma-ba
c        79. k0 + n- to pion+ + sigma+ba
c        80. k0 + n- to pion0 + sigma0-
c        81. k0 + n- to pion0 + lambda-
c        82. k0 + n- to k- + cascade-ba (strangeness production)

c        83. k- + lambda to pion0 + cascade-
c        84. k- + sigma+ to pion+ + cascade-
c        85. k- + sigma- to pion- + cascade-
c        86. k- + sigma0 to pion0 + cascade-
c        87. k0- + lambda to pion+ + cascade-
c        88. k0- + sigma0 to pion+ + cascade-
c        89. k0- + sigma- to pion0 + cascade-
c        90. k+ + lambda- to pion0 + cascade-ba
c        91. k+ + sigma+ba to pion- + cascade-ba
c        92. k+ + sigma-ba to pion+ + cascade-ba
c        93. k+ + sigma0- to pion0 + cascade-ba
c        94. k0 + lambda- to pion- + cascade-ba
c        95. k0 + sigma0- to pion- + cascade-ba
c        96. k0 + sigma-ba to pion0 + cascade-ba
c        97 pion+ + sigma- to k0 + cascade0
c        98 pion+ + sigma0 to k+ + cascade0
c        99 pion+ + lambda0 to k+ + cascade0
c        100    pion- + sigma+ to k0 + cascade0
c        101    pion0 + sigma+ to k+ + cascade0
c        102    pion0 + sigma0 to k0 + cascade0
c        103    pion0 + lambda to k0 + cascade0
c        104    pion+ + sigma+ba to k0- + cascade0-
c        105    pion- + sigma-ba to k0- + cascade0-
c        106    pion- + sigma0- to k- + cascade0-
c        107    pion- + lambda- to k- + cascade0-
c        108    pion0 + sigma+- to k- + cascade0-
c        109    pion0 + sigma0- to k0- + cascade0-
c        110    pion0 + lambda- to k0- + cascade0-
c        111    k- + sigma+ to pion0 + cascade0
c        112    k- + sigma0 to pion0 + cascade-
c     113  k- + lambda to pion- + cascade0
c     114  k0- + sigma+ to pion+ + cascade0
c     115  k0- + sigma- to pion- + cascade0
c     116  k0- + sigma0 to pion0 + cascade0
c     117  k0- + lambda to pion0 + cascade0
c     118  k+ + sigma+ba to pion0 + cascade0-
c     119  k+ + sigma0- to pion+ + cascade0-
c     120  k+ + lambda- to pion+ + cascade0-
c     121  k+ + cascade-ba to pion+ + omiga-ba
c     122  k0 + sigma-ba to pion+ + cascade0ba
c     123       k0 + sigma0- to pion0 + cascade0-
c     124       k0 + lambda- to pion0 + cascade0ba
c     125       k- + p to k0 + cascade0
c     126       k0- + p to k+ + cascade0
c     127       k0- + n to k0 + cascade0
c     128       k+ + p- to k0- + cascade0ba
c     129       k0 + p- to k- + cascade0-
c     130       k0 + n- to k0- + cascade0-
c     131       pion+ + cascade- to k+ + omiga-
c     132       pion0 + cascade- to k0 + omiga-
c     133       pion- + cascade-ba to k- + omiga-ba
c     134       pion0 + cascade-ba to k0- + omiga-ba
c     135       pion- + cascade0 to k0 + omiga-
c     136       pion0 + cascade0 to k+ + omiga-
c     137       pion+ + cascade0- to k0- + omiga-ba
c     138       pion0 + cascade0- to k- + omiga-ba
c     139     k- + cascade- to pion- + omiga-
c     140     k0- + cascade- to pion0 + omiga-
c     141     k- + cascade0 to pion0 + omiga-
c     142     k0- + cascade0 to pion+ + omiga-
c     143     k+ + cascade-ba to pion+ + omiga-ba
c     144     k0 + cascade-ba to pion0 + omiga-ba
c     145     k+ + cascade0- to pion0 + omiga-ba
c     146     k0 + cascade0- to pion- + omiga-ba
c     147       pion- + p to delta- + pion+
c     148       pion- + p to rho0 + n
c     149       pion- + p to rho- + p
c     150       pion- + p to delta+ + pion-
c     151       pion- + p to delta0 + pion0
c     152       pion- + n to delta- + pion0
c     153       pion- + n to rho- + n
c     154       pion- + n to delta0 + pion-
c     155       pion+ + p to delta++ + pion0
c     156       pion+ + p to delta+ + pion+
c     157       pion+ + p to rho+ + p
c     158       pion+ + n to delta++ + pion-
c     159       pion+ + n to delta0 + pion+
c     160       pion+ + n to delta+ + pion0
c     161       pion+ + n to rho0 + p
c     162       pion+ + n to rho+ + n
c     163       pion0 + p to delta0 + pion+
c     164       pion0 + p to delta++ + pion-
c     165       pion0 + p to rho+ + n
c     166       pion0 + p to rho0 + p
c     167       pion0 + p to delta+ + pion0
c     168       pion0 + n to delta+ + pion-
c     169       pion0 + n to delta- + pion+
c     170       pion0 + n to delta0 + pion0
c     171       pion0 + n to rho0 + n
c     172       pion0 + n to rho- + p
c     173       p + p to delta+ + p
c     174       p + p to delta++ + n
c     175       p + n to delta+ + n
c     176       p + n to delta0 + p
c     177       n + n to delta0 + n
c     178       n + n to delta- + p
c     there are reverse reactions after '201'
c               201  k+ + k- to pion+ + pion-
c               202  k+ + k- to pion0 + pion0 
c               203  k+ + k0- to pion+ + pion0
c               204  k- + k0 to pion- + pion0
c               205  k0 + k0- to pion+ + pion-
c               206  k0 + k0- to pion0 + pion0
c               207       k+ + sigma+ to pion+ + p
c               208       k+ + sigma- to pion- + p
c               209       k+ + sigma- to pion0 + n
c               210       k+ + sigma0 to pion+ + n
c               211       k+ + sigma0 to pion0 + p
c               212       k+ + lambda0 to pion+ + n
c               213       k+ + lambda0 to pion0 + p
c               214       k+ + sigma+ to pion+ + n
c               215       k+ + sigma+ to pion0 + p
c               216       k0 + sigma- to pion- +n 
c               217       k0 + sigma0 to pion- + p
c               218       k0 + sigma0 to pion0 + n
c               219       k0 + lambda0 to pion- + p
c               220       k0 + lambda0 to pion0 + n
c               221       k- + sigma+ba to pion- + pba
c               222       k- + sigma-ba to pion+ + pba
c               223       k- + sigma-ba to pion0 + nba
c               224       k- + sigma0ba to pion- + nba
c               225       k- + sigma0ba to pion0 + pba
c               226       k- + lambda0ba to pion- + nba
c               227       k- + lambda0ba to pion0 + pba
c               228       k0ba + sigma+ba to pion- + nba
c               229       k0ba + sigma+ba to pion0 + pba
c               230       k0- + sigma-ba to pion+ + nba
c       231 k0- + sigma0- to pion+ + pba
c       232 k0- + sigma0- to pion0 + nba
c       233 k0- + lambda0- to pion+ + pba
c       234 k0- + lambda0- to pion0 + nba
c       235 k++ cascade- to pi+ sigma-
c       236 k+ + cascade- to pion- + sigma+
c       237 k+ + cascade- to pion0 + lambda0
c       238 k+ + cascade- to pion0 + sigma0
c       239 k- + cascade-ba to pion+ + sigma+ba
c       240 k- + cascade-ba to pion- + sigma-ba
c       241 k- + cascade-ba to pion0 + lambda0-
c       242 k- + cascade-ba to pion0 + sigma0-
c       243 k0 + cascade- to pion- + lambda0
c       244 k0 + cascade- to pion- + sigma0
c       245 k0 + cascade- to pion0 + sigma-
c       246 k0- + cascade-ba to pion+ + lambda0-
c       247 k0- + cascade-ba to pion+ + sigma0-
c       248 k0- + cascade-ba to pion0 + sigma-ba
c       249 pion+ + sigma- to k- + p
c       250 pion+ + sigma- to k0- + n
c       251 pion+ + sigma0 to k0- + p
c       252 pion+ + lambda0 to k0- + p
c       253 k+ + cascade- to k- + p
c       254 pion- + sigma+ to k- + p
c       255 pion- + sigma+ to k0- + n
c       256 pion- + sigma0 to k- + n
c       257 k0 + cascade- to k- + n
c       258 pion- + lambda0 to k- + n
c       259 pion0 + sigma+ to k0- + p
c       260 pion0 + sigma- to k- + n
c       261 pion0 + sigma0 to k- + p
c       262 pion0 + sigma0 to k0- + n
c       263 pion0 + lambda0 to k- + p
c       264 pion0 + lambda0 to k0- + n
c       265 k+ + cascade- to k0- + n
c       266 pion+ + sigma+ba to k+ + pba
c       267 pion+ + sigma+ba to k0 + nba
c       268 pion+ + sigma0- to k+ + nba
c       269 pion+ + lambda0- to k+ + nba
c       270 k- + cascade-ba to k+ + pba
c       271 pion- + sigma-ba to k+ + pba
c       272 pion- + sigma-ba to k0 + nba
c       273 pion- + sigma0- to k0 + pba
c       274 k0- + cascade-ba to k+ + nba
c       275 pion- + lambda0- to k0 + pba
c       276 pion0 + sigma+- to k0 + pba
c       277 pion0 + sigma-ba to k+ + nba
c       278 pion0 + sigma0- to k+ + pba
c       279 pion0 + sigma0- to k0 + nba
c       280 pion0 + slambda0- to k+ + pba
c       281 pion0 + lambda0- to k0 + nba
c       282 k- + cascade-ba to k0 + nba
c       283 pion+ + cascade- to k- + sigma+
c       284 pion+ + cascade- to k0- + lambda0
c       285 pion+ + cascade- to k0- + sigma0
c       286 pion- + cascade- to k- + sigma-
c       287 pion0 + cascade- to k- + lambda0
c       288 pion0 + cascade- to k- + gigma0
c       289 pion0 + cascade- to k0- + sigma-
c       290 pion+ + cascade-ba to k+ + sigma-ba
c       291 pion- + cascade-ba to k+ + sigma+-
c       292 pion- + cascade-ba to k0 + lambda0-
c       293 pion- + cascade-ba to k0 + sigma0-
c       294 pion0 + cascade-ba to k+ + lambda0-
c       295 pion0 + cascade-ba to k+ + gigma0-
c       296 pion0 + cascade-ba to k0 + sigma-ba
c       297 k0 + cascade0 to pion+ + sigma-
c       298 k+ + cascade0 to pion+ + sigma0
c       299 k+ + cascade0 to pion+ + lambda
c       300 k0 + cascade0 to pion- + sigma+
c       301 k+ + cascade0 to pion0 + sigma+
c       302 k0 + cascade0 to pion0 + sigma0
c       303 k0 + cascade0 to pion0 + lambda
c       304 k-0 + cascade0- to pion+ + sigma-ba
c       305 k0- + cascade0- to pion- + sigma-ba
c       306 k- + cascade0- to pion- + sigma0-
c       307 k- + cascade0- to pion- + lambda-
c       308 k- + cascade0- to pion0 + sigma+ba
c       309 k0- + cascade0- to pion0 + sigma0-
c       310 k0- + cascade0- to pion0 + lambda-
c       311 pion0 + cascade0 to k- + sigma+
c       312 pion- + cascade0 to k- + sigma0
c       313 pion- + cascade0 to k- + lambda
c       314 pion+ + cascade0 to k0- + sigma+
c       315 pion- + cascade0 to k0- + sigma-
c       316 pion0 + cascade0 to k0- + sigma0
c       317 pion0 + cascade0 to k0- + lambda
c       318 pion0 + cascade0- to k+ + sigma+ba
c       319 pion+ + cascade0- to k+ + sigma0-
c       320 pion+ + cascade0- to k+ + lambda-
c       321 pion- + cascade0- to k0 + sigma+ba
c       322 pion+ + cascade0- to k0 + sigma-ba
c       323 pion0 + cascade0- to k0 + sigma0-
c       324 pion0 + cascade0- to k0 + lambda-
c       325 k0 + cascade0 to k- + p 
c       326 k+ + cascade0 to k0- + p
c       327 k0 + cascade0 to k0- + n
c       328 k0- + cascade0- to k+ + p- 
c       329 k- + cascade0- to k0 + p-
c       330 k0- + cascade0- to k0 + n-
c       331 k+ + omiga- to pion+ + cascade-
c       332 k0 + omiga- to pion0 + cascade-
c       333 k- + omiga-ba to pion- + cascade-ba
c       334 k0- + omiga-ba to pion0 + cascade-ba
c       335 k0 + omiga- to pion- + cascade0
c       336 k+ + omiga- to pion0 + cascade0
c       337 k0- + omiga-ba to pion+ + cascade0-
c       338 k- + omiga-ba to pion0 + cascade0-
c       339 pion- + omiga- to k- + cascade-
c       340 pion0 + omiga- to k0- + cascade-
c       341 pion0 + omiga- to k- + cascade0
c       342 pion+ + omiga- to k0- + cascade0
c       343 pion+ + omiga-ba to k+ + cascade-ba
c       344 pion0 + omiga-ba to k0 + cascade-ba
c       345 pion0 + omiga-ba to k+ + cascade0-
c       346 pion- + omiga-ba to k0 + cascade0-
c       347 pion+ + delta- to pion- + p
c       348 pion+ + delta- to pion0 + n
c       349     pion+ + delta0 to pion+ + n
c       350     pion+ + delta0 to pion0 + p
c       351     pion+ + delta+ to pion+ + p
c       352     pion0 + delta++ to pion+ p
c       353     pion0 + delta+ to pion0 + p
c       354     pion0 + delta+ to pion+ + n
c       355     pion0 + delta0 to pion0 + n
c       356     pion0 + delta0 to pion- + p
c       357     pion0 + delta- to pion- + n
c       358     pion- + delta++ to pion0 + p
c       359     pion- + delta++ to pion+ + n
c       360     pion- + delta+ to pion- + p
c       361     pion- + delta+ to pion0 + n
c       362     pion- + delta0 to pion- + n
c       363     rho0 + n to pion- + p
c       364     rho0 + n to pion0 + n
c       365     rho- + n to pion- + n
c       366     rho+ + n to pion0 + p
c       367     rho+ + n to pion+ + n
c       368     rho0 + p to pion0 + p
c       369     rho0 + p to pion+ + n
c       370     rho- + p to pion0 + n
c       371     rho- + p to pion- + p
c       372     rho+ + p to pion+ + p
c       373     delta++ + n to p + p
c       374     delta+ + n  to p + n
c       375     delta+ + p  to p + p
c       376     delta0 + p  to p + n
c       377     delta0 + n  to n + n
c       378     delta- + p  to n + n
c       393     lambda- + p to K*+ + omiga
c       394     lambda- + n to K*0 + omiga
c       395     sigma0- + p to K*+ + omiga
c       396     sigma0- + n to K*0 + omiga
c       397     p- + p to rho0 + omiga
c       398     p- + n to rho- + omiga
c       399     n- + p to rho+ + omiga
c       400     n- + n to rho0 + omiga
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       increase the reactions involving cascade0 and omiga;27/04/95
      PARAMETER (KSZJ=40000,KSZ1=30)
        parameter(nsize=100000)
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
        COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/sa10/csnn,cspin,cskn,cspipi,rcsit,ifram
c       common/sa11/an(20,5,20),bn(20),san(20,5,20),sbn(20)
        double precision b(3)
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        dimension pi(4),pj(4),pii(4),pjj(4)
        integer winel
        jjj=0
        ww=tw(icp)
        kl=k(l,2)
        kl1=k(l1,2)
        am01=pmas(lucomp(kl),1)
        am02=pmas(lucomp(kl1),1)
c the following statements guarantees the line number of pions and kions 
c correspond to the first position of subroutine prod()
c       tw : the cross section ratio of (special inela.)/tot
        if(abs(kl).eq.211 .or. abs(kl).eq.111 .or. abs(kl).eq.321 
     &   .or. abs(kl).eq.311.or. abs(kl).eq.213.or. abs(kl).eq.113)then
        idpl=1
        else
        idpl=3
        endif
c put the meson of a colliding pair in the first position in prod(). But if
Cthe collision occurs between two baryons (or two mesons) the position is uncertain.
        if(idpl.eq.1)call prod(l,l1,kl,kl1,ss,icp,lc,tc,tw)
        if(idpl.eq.3)call prod(l1,l,kl1,kl,ss,icp,lc,tc,tw)     
        ik1=lc(icp,3)
        ik2=lc(icp,4)
        w1=tw(icp)/rcsit
c       1/rcsit : the cross section ratio of tot/(inela.)
        ww=1.
        if(rlu(1).gt.w1)then
        winel=0
c       treated as elastic then
        return
        endif
        icp5=lc(icp,5)
c       if(icp5.eq.397 .or. icp5.eq.398 .or. icp5.eq.399.or.icp5.
c     & eq.400)goto 400
100     fi1=atan2(pi(2),pi(1))
        cta1=atan2(sqrt(pi(1)**2+pi(2)**2),pi(3))
        cfi1=cos(fi1)
        sfi1=sin(fi1)
        ccta1=cos(cta1)
        scta1=sin(cta1)
        am1=pmas(lucomp(ik1),1)
        am2=pmas(lucomp(ik2),1)
        pp=(ss*ss-(am1+am2)**2)*(ss*ss-(am1-am2)**2)/(4.*ss*ss)
        if(pp.lt.0.)pp=1.e-10
        pp=sqrt(pp)
        pii(4)=(ss*ss+am1**2-am2**2)/(2.*ss)
c       energy of one particle (between two) after scattering
        
        fis=2.*3.1415926*rlu(1)
        cfis=cos(fis)
        sfis=sin(fis)
        if(KFR(27).eq.1)then
        call coinelas(am01,am1,ss,pi,pp,cctas)
        else
        cctas=2*rlu(1)-1.
        endif
        sctas=sqrt(1.-cctas*cctas)
111     continue
        pii(1)=cfi1*(ccta1*sctas*cfis+scta1*cctas)-sfi1*sctas*sfis
        pii(2)=sfi1*(ccta1*sctas*cfis+scta1*cctas)+cfi1*sctas*sfis
        pii(3)=ccta1*cctas-scta1*sctas*cfis
        pii(1)=pp*pii(1)
        pii(2)=pp*pii(2)
        pii(3)=pp*pii(3)
        do i=1,3 
        pjj(i)=0.-pii(i)
        enddo
        pjj(4)=ss-pii(4)
        if(pii(4).lt.0. .or. pjj(4).lt.0.)then
        write(mstu(11),*)'error may happen here in subroutine
     &         coinel(),energy is negative','chann=',icp5
        write(mstu(11),*)'pi,pj,pii,pjj=',pi(4),pj(4),pii(4),pjj(4)   
        endif   
        ilo=1
        call lorntz(ilo,b,pii,pjj)
        if(pii(4).lt.0. .or. pjj(4).lt.0.)then
        write(mstu(11),*)'error may happen here in subroutine
     &         coinel(),energy is negative','chann=',icp5
        endif   
c       if(jjj.eq.1)goto 300
400     return
        end


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine coinelas(am01,am1,eij,pi,pp,cctas)
c       perform inelastic scattering
        dimension pi(4)
        double precision abt,abu
        COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
c       d=3.65*(eij-am1-am2)
c       if(d.lt.1.e-10)return
c       pt=0.2
c       a=min(10.3,1./(1.12*pt)/(1.12*pt))
c       d6=d**6
c       b=d6*a/(1.+d6)
        b=10.3
c       if(b.lt.1.e-20)then
c       b=1.e-20
c       endif
        pm2=pi(1)**2+pi(2)**2+pi(3)**2
        pm=sqrt(pm2)
        em=sqrt(pm*pm+am01*am01)
        em1=sqrt(pp*pp+am1*am1)
        tmin=am01**2+am1**2-2*(em*em1+pm*pp)
        tmax=am01**2+am1**2-2*(em*em1-pm*pp)
        if(tmin.gt.tmax)write(mstu(11),*)tmin,tmax
        cc=rlu(1)

        abt=dexp(dmax1(-7.0D2,dble(b*tmin)))
        abu=dexp(dmax1(-7.0D2,dble(b*tmax)))

        tt1=dlog(cc*abu+(1.-cc)*abt)

        tt=tt1/b
        cctas=(0.5*(tt-am01**2-am1**2)+em*em1)/(pm*pp)

        if(abs(cctas).gt.1)cctas=sign(1.,cctas)
        return
        end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


        subroutine prod(l,l1,kl,kl1,ss,icp,lc,tc,tw)
        parameter(nsize=100000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       l must be larger than numb(4) (pion,k-,etc.)  
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/LUCIDAT2/KFMAXT,PARAM(20),WEIGH(400)
        common/count/isinel(400)
        COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn
     c  ,rnt,rnp,rao,rou0       
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        common/sa10/csnn,cspin,cskn,cspipi,rcsit,ifram
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        integer fact1,fact2
        para13=PARAM(13)*0.1
        PARAM14=PARAM(14)*0.1
c       p- p annihilation cross section is assumed to be equal to 0.8*(
c        total inelastic cross section)
        ioo=0
        ilo=1
        ilo1=1
        ilo2=1
        ilo3=1
        ilo4=1
        ilo5=1
        ilo6=1
        ilo7=1
        ilo8=1
        ilo9=1

        nchargei=plu(l,6)+plu(l1,6)
c p-p ------------------------>
        if(kl.eq. 2212.and. kl1.eq.2212)then    
        call ppdelta(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)   
        if(ioo.eq.0)goto 13
        goto 10
cp+n ----------------------->
        elseif((kl.eq. 2212.and. kl1.eq.2112).or.
     &          (kl.eq. 2112.and. kl1.eq.2212))then    
        call pndelta(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
        if(ioo.eq.0)goto 13
        goto 10
c n+n ------->
        elseif(kl.eq. 2112.and. kl1.eq.2112)then  
        call nndelta(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
        if(ioo.eq.0)goto 13
        goto 10
        endif

c       pion+ + pion- 
        if((kl.eq.211 .and. kl1.eq.-211)
     c   .or.(kl.eq.-211 .and. kl1.eq.211))then
        if(isinel(1).eq.0)then
        fact1=0
        goto 101
        endif
        fact1=1
        ik3=-321
        ik4=321
        ic3=1
        call spipi(ik3,ik4,ss,ilo1)
        if(ilo1.eq.0)fact1=0.
101     if(isinel(2).eq.0)then
        fact2=0
        goto 102
        endif
        fact2=1
        ik5=-311
        ik6=311
        ic5=2
        call spipi(ik5,ik6,ss,ilo2)
        if(ilo2.eq.0)fact2=0.
102     fact=fact1+fact2
        if(fact1.eq.0. .and. fact2.eq.0.)goto 13
c       if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=ic3
        if(rlu(1).gt.fact1/fact)then
        lc(icp,3)=ik5
        lc(icp,4)=ik6
        lc(icp,5)=ic5
        endif
        tw(icp)=fact*sig/cspipi
        goto 10
        endif

c       pion+ + pion0
        if((kl.eq.211 .and. kl1.eq.111)
     c   .or.(kl.eq.111 .and. kl1.eq.211))then
        if(isinel(3).eq.0)goto 13
        ik3=-311
        ik4=321
        call spipi(ik3,ik4,ss,ilo1)
        if(ilo1.eq.0)goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=3
        tw(icp)=sig/cspipi
        goto 10
        endif
        
c       pion- + pion0
        if((kl.eq.-211 .and. kl1.eq.111)
     c   .or.(kl.eq.111 .and. kl1.eq.-211))then
        if(isinel(4).eq.0)goto 13
        ik3=311
        ik4=-321
        call spipi(ik3,ik4,ss,ilo1)
        if(ilo1.eq.0)goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=4
        tw(icp)=sig/cspipi
        goto 10
        endif

c       pion0 + pion0 
        if(kl.eq.111 .and. kl1.eq.111)then
        if(isinel(5).eq.0)then
        fact1=0.
        goto 103
        endif
        fact1=1.
        ik3=-321
        ik4=321
        ic3=5
        call spipi(ik3,ik4,ss,ilo1)
        if(ilo1.eq.0)fact1=0.
103     if(isinel(6).eq.0)then
        fact2=0.
        goto 104
        endif
        fact2=1.
        ik5=-311
        ik6=311
        ic5=6
        call spipi(ik5,ik6,ss,ilo2)
        if(ilo2.eq.0)fact2=0.
104     fact=fact1+fact2
        if(fact1.eq.0. .and. fact2.eq.0.)goto 13
c       if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=ic3
        if(rlu(1).gt.fact1/fact)then
        lc(icp,3)=ik5
        lc(icp,4)=ik6
        lc(icp,5)=ic5
        endif
        tw(icp)=fact*sig/cspipi
        goto 10
        endif

c       pion+ + p 

500     if(kl.eq.211 .and. kl1.eq.2212)then
        call pip2(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
        if(ioo.eq.0)goto 13
        goto 10
        endif

c       pion+ + n 
        if(kl.eq.211 .and. kl1.eq.2112)then
        call pin1(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
        if(ioo.eq.0)goto 13
        goto 10
        endif


c       pion- + p
        if(kl.eq.-211 .and. kl1.eq.2212)then
        call pip1(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
        if(ioo.eq.0)goto 13
        goto 10
        endif

c       pion- + n 
        if(kl.eq.-211 .and. kl1.eq.2112)then
        call pin3(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
        if(ioo.eq.0)goto 13
        goto 10
        endif

c       pion0 + p
        if(kl.eq.111 .and. kl1.eq.2212)then  
        call pip3(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
        if(ioo.eq.0)goto 13
        goto 10
        endif

c       pion0 + n 
        if(kl.eq.111 .and. kl1.eq.2112)then
        call pin2(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
        if(ioo.eq.0)goto 13
        goto 10
        endif

c       pion+ + pba 
        if(kl.eq.211 .and. kl1.eq.-2212)then
        if(isinel(21).eq.0)then
        si1=0.
        goto 319
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)
c       cross section of pion+ + pba to k0- + lambdaba
c       cross section of pion+ + pba to k0- + lambdaba is assumed to be 
c        equal to pion- + p to k0 + lambda = s1724
319     if(isinel(22).eq.0)then
        si2=0.
        goto 320
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)
c       cross section of pion+ + pba to k0- + sigma0-
320     if(isinel(23).eq.0)then
        si3=0.
        goto 321
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3112),1)
        si3=s1724(ss,ilo3,0,the)
c       cross section of pion+ + pba to k- + sigma-ba
321     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6)goto 13
        si12=si1+si2
        sit=si12+si3
        s1=si1/sit
        s2=si12/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-311
        ik2=-3122
        ic=21
c       pion+ + pba to k0- + lambda-
        goto 322
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-311
        ik2=-3212
        ic=22
c       pion+ + pba to k0- + sigma0-
        goto 322
        endif
        ik1=-321
        ik2=-3112
        ic=23
c       pion+ + pba to k- + sigma-ba
322     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif

c       pion+ + nba to k0- + sigma-ba
        if(kl.eq.211 .and. kl1.eq.-2112)then
        if(isinel(24).eq.0)goto 13
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3112),1)
        tw(icp)=s1724(ss,ilo,0,the)/cspin/10.
c       cross section of pion+ + nba to k0- + sigma-ba is assumed to be 
c        equal to pion- + n to k + y (isotropic averaged)
        if(ilo.eq.0)goto 13
        lc(icp,3)=-311
        lc(icp,4)=-3112
        lc(icp,5)=24
        goto 10
        endif

c       pion- + pba to k- + sigma+ba
        if(kl.eq.-211 .and. kl1.eq.-2212)then
        if(isinel(25).eq.0)goto 13
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3222),1)
        tw(icp)=s1724(ss,ilo,0,the)/cspin/10.
c       cross section of pion- + pba to k- + sigma+ba is assumed to be 
c        equal to pion+ + p to k + y (isotropic averaged)
        if(ilo.eq.0)goto 13
        lc(icp,3)=-321
        lc(icp,4)=-3222
        lc(icp,5)=25
        goto 10
        endif

c       pion- + nba
        if(kl.eq.-211 .and. kl1.eq.-2112)then
        if(isinel(26).eq.0)then
        si1=0.
        goto 323
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)
c       cross section of pion- + nba to k- + y- is assumed to be 
c        equal to pion+ + n to k + y
323     if(isinel(27).eq.0)then
        si2=0.
        goto 324
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)
324     if(isinel(28).eq.0)then
        si3=0.
        goto 325
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3222),1)
        si3=s1724(ss,ilo3,0,the)
325     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6)goto 13
        si12=si1+si2
        sit=si12+si3
        s1=si1/sit
        s2=si12/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-321
        ik2=-3122
        ic=26
c       pion- + nba to k- + lambdaba
        goto 326
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-321
        ik2=-3212
        ic=27
c       pion- + nba to k- + sigma0ba
        goto 326
        endif
        ik1=-311
        ik2=-3222
        ic=28
c       pion-+ nba to k0- + sigma+ba
326     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif

c       pion0 + pba
        if(kl.eq.111 .and. kl1.eq.-2212)then
        if(isinel(29).eq.0)then
        si1=0.
        goto 327
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)
327     if(isinel(30).eq.0)then
        si2=0.
        goto 328
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)
328     if(isinel(31).eq.0)then
        si3=0.
        goto 329
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3222),1)
        si3=s1724(ss,ilo3,0,the)
329     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6)goto 13
        si12=si1+si2
        sit=si12+si3
        s1=si1/sit
        s2=si12/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-321
        ik2=-3122
        ic=29
c       pion0 + pba to k- + lambda-
        goto 330
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-321
        ik2=-3212
        ic=30
c       pion0 + pba to k- + sigma0-
        goto 330
        endif
        ik1=-311
        ik2=-3222
        ic=31
c       pion0 + pba to k0- + sigma+ba
330     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif

c       pion0 + nba
        if(kl.eq.111 .and. kl1.eq.-2112)then
        if(isinel(32).eq.0)then
        si1=0.
        goto 331
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3112),1)
        si1=s1724(ss,ilo1,0,the)
331     if(isinel(33).eq.0)then
        si2=0.
        goto 332
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3122),1)
        si2=s1724(ss,ilo2,0,the)
332     if(isinel(34).eq.0)then
        si3=0.
        goto 333
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3212),1)
        si3=s1724(ss,ilo3,0,the)
333     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6)goto 13
        si12=si1+si2
        sit=si12+si3
        s1=si1/sit
        s2=si12/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-321
        ik2=-3112
        ic=32
c       pion0 + nba to k- + sigma-ba
        goto 334
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-311
        ik2=-3122
        ic=33
c       pion0 + nba to k0- + lambdaba
        goto 334
        endif
        ik1=-311
        ik2=-3212
        ic=34
c       pion0 + nba to k0- + sigma0ba
334     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif

c       pion+ + sigma- 
        if(kl.eq.211 .and. kl1.eq.3112)then
        if(isinel(35).eq.0)then
        si1=0.
        goto 701
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
c       cross section of pion + y to kaon + cascade is assumed to be 
c        equal to pion + n to kaon + y,but take the different of
c        threshold energy into account

701     if(isinel(97).eq.0)then
        si2=0.
        goto 660
        endif
        the=pmas(lucomp(311),1)+pmas(lucomp(3322),1)
        si2=s1724(ss,ilo2,0,the)

660     if(isinel(249).eq.0)then
        si3=0.
        goto 661
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(2212),1)
        ik1=-321
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10*s1724(ss,ilo3,0,the)*fac
661     if(isinel(250).eq.0)then
        si4=0.
        goto 702
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(2112),1)
        ik1=-311
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10*s1724(ss,ilo4,0,the)*fac
702     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &  si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=rlu(1)     
         if(rlus.le.s1)then

        ik1=321
        ik2=3312
        ic=35
c       pion+ + sigma- to k+ + cascade-
        goto 703
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=311
        ik2=3322
        ic=97
c       pion+ + sigma- to k0 + cascade0
        goto 703
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=-321
        ik2=2212
        ic=249
c       pion+ + sigma- to k- + p
        goto 703
        endif
        ik1=-311
        ik2=2112
        ic=250
c       pion+ + sigma- to k0- + n
703     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion- + sigma-ba
        if(kl.eq.-211 .and. kl1.eq.-3112)then
        if(isinel(45).eq.0)then
        si1=0.
        goto 7011
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
c       cross section of pion + y to kaon + cascade is assumed to be 
c        equal to pion + n to kaon + y,but take the different of
c        threshold energy into account

7011    if(isinel(105).eq.0)then
        si2=0.
        goto 6601
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)

6601    if(isinel(271).eq.0)then
        si3=0.
        goto 6611
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(-2212),1)
        ik1=321
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10*s1724(ss,ilo3,0,the)*fac
6611    if(isinel(272).eq.0)then
        si4=0.
        goto 7021
        endif
        the=pmas(lucomp(311),1)+pmas(lucomp(-2112),1)
        ik1=311
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10*s1724(ss,ilo4,0,the)*fac
7021    if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.
     &  and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.
     &  and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=rlu(1)     
         if(rlus.le.s1)then

        ik1=-321
        ik2=-3312
        ic=45
c       pion- + sigma-ba to k- + cascade-ba
        goto 7031
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-311
        ik2=-3322
        ic=105
c       pion- + sigma-ba to k0- + cascade0-
        goto 7031
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=321
        ik2=-2212
        ic=271
c       pion- + sigma-ba to k+ + pba
        goto 7031
        endif
        ik1=311
        ik2=-2112
        ic=272
c       pion- + sigma-ba to k0 + nba
7031    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif





c       pion- + sigma+ 
        if(kl.eq.-211 .and. kl1.eq.3222)then
        if(isinel(37).eq.0)then
        si1=0.
        goto 704
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
 
704     if(isinel(100).eq.0)then
        si2=0.
        goto 663
        endif
        the=pmas(lucomp(311),1)+pmas(lucomp(3322),1)
        si2=s1724(ss,ilo2,0,the)

663        if(isinel(254).eq.0)then
        si3=0.
        goto 664
        endif
        ik1=-321
        ik2=2212
        the=pmas(lucomp(-321),1)+pmas(lucomp(2212),1)
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10*s1724(ss,ilo3,0,the)*fac
664    if(isinel(255).eq.0)then
        si4=0.
        goto 705
        endif
        ik1=-311
        ik2=2112
        the=pmas(lucomp(-311),1)+pmas(lucomp(2112),1)
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10*s1724(ss,ilo4,0,the)*fac
705     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.
     &  and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=321
        ik2=3312
        ic=37
c       pion- + sigma+ to k+ + cascade-
        goto 706
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=311
        ik2=3322
        ic=100
c       pion- + sigma+ to k0 + cascade0
        goto 706
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=-321
        ik2=2212
        ic=254
c       pion- + sigma+ to k- + p
        goto 706
        endif
        ik1=-311
        ik2=2112
        ic=255
c       pion- + sigma+ to k0- + n
706     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

cc      pion+ + sigma+bar
        if(kl.eq.211 .and. kl1.eq.-3222)then
        if(isinel(43).eq.0)then
        si1=0.
        goto 7041
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
 
7041     if(isinel(104).eq.0)then
        si2=0.
        goto 6631
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)

6631        if(isinel(266).eq.0)then
        si3=0.
        goto 6641
        endif
        ik1=321
        ik2=-2212
        the=pmas(lucomp(321),1)+pmas(lucomp(-2212),1)
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10*s1724(ss,ilo3,0,the)*fac
6641    if(isinel(267).eq.0)then
        si4=0.
        goto 7051
        endif
        ik1=311
        ik2=-2112
        the=pmas(lucomp(311),1)+pmas(lucomp(-2112),1)
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10*s1724(ss,ilo4,0,the)*fac
7051     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.
     &       and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-321
        ik2=-3312
        ic=43
c       pion+ + sigma+bar to k- + cascade-bar
        goto 7061
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-311
        ik2=-3322
        ic=104
c       pion+ + sigma+bar to k0- + cascade0bar
        goto 7061
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=321
        ik2=-2212
        ic=266
c       pion+ + sigma+bar to k+ + pbar
        goto 7061
        endif
        ik1=311
        ik2=-2112
        ic=267
c       pion+ + sigma+bar to k0 + nbar
7061     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif
c       pion- + sigma0 
        if(kl.eq.-211 .and. kl1.eq.3212)then
        if(isinel(38).eq.0)then
        si1=0.
        goto 1681
        endif
        the=pmas(lucomp(311),1)+pmas(lucomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
1681     if(isinel(256).eq.0)then
        si2=0.
        goto 1682
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(2112),1)
       ik3=-321
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10.*fac
1682     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=311
        ik2=3312
        ic=38
c   c   pion- + sigma0 to k0 + cascade-
        goto 683
        endif
        ik1=-321
        ik2=2112
        ic=256
c    c       pion- + sigma0 to k- + n
683     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif


c       pion+ + sigma0 bar
        if(kl.eq.211 .and. kl1.eq.-3212)then
        if(isinel(44).eq.0)then
        si1=0.
        goto 1781
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
1781     if(isinel(268).eq.0)then
        si2=0.
        goto 1782
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(-2112),1)
       ik3=321
        ik4=-2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10.*fac
1782     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-311
        ik2=-3312
        ic=44
c   c   pion+ + sigma0bar to k0- + cascade-bar
        goto 6831
        endif
        ik1=321
        ik2=-2112
        ic=268
c    c       pion+ + sigma0bar to k+ + nbar
6831     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif


c       pion0 + sigma0 
        if(kl.eq.111 .and. kl1.eq.3212)then
        if(isinel(41).eq.0)then
        si1=0.
        goto 710
        endif   
        the=pmas(lucomp(321),1)+pmas(lucomp(3312),1)
        si1=s1724(ss,ilo1,0,the)

710     if(isinel(102).eq.0)then
        si2=0.
        goto 666
        endif
        the=pmas(lucomp(311),1)+pmas(lucomp(3322),1)
        si2=s1724(ss,ilo2,0,the)

666        if(isinel(261).eq.0)then
        si3=0.
        goto 667
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(2212),1)
        ik1=-321
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10*s1724(ss,ilo3,0,the)*fac
667         if(isinel(262).eq.0)then
        si4=0.
        goto 711
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(2112),1)
        ik1=-311
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10*s1724(ss,ilo4,0,the)*fac
711   if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.
     &  and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=321
        ik2=3312
        ic=41
c       pion0 + sigma0 to k+ + cascade-
        goto 712
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=311
        ik2=3322
        ic=102
c       pion0 + sigma0 to k0 + cascade0
        goto 712
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=-321
        ik2=2212
        ic=261
c      pion0 + sigma0 to k- + p
        goto 712
        endif
       ik1=-311
        ik2=2112
        ic=262
c       pion0 + sigma0 to k0- + n
712     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif
c       pion0 + sigma0bar
        if(kl.eq.111 .and. kl1.eq.-3212)then
        if(isinel(48).eq.0)then
        si1=0.
        goto 1910
        endif   
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)

1910     if(isinel(109).eq.0)then
        si2=0.
        goto 1666
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)

1666        if(isinel(278).eq.0)then
        si3=0.
        goto 1667
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(-2212),1)
        ik1=321
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si3=10*s1724(ss,ilo3,0,the)*fac
1667        if(isinel(279).eq.0)then
        si4=0.
        goto 5711
        endif
        the=pmas(lucomp(311),1)+pmas(lucomp(-2112),1)
        ik1=311
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si4=10*s1724(ss,ilo4,0,the)*fac
5711    if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.
     &  and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-321
        ik2=-3312
        ic=48
c       pion0 + sigma0bar to k- + cascade-bar
        goto 7125
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-311
        ik2=-3322
        ic=109
c       pion0 + sigma0bar to k0- + cascade0bar
        goto 7125
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=321
        ik2=-2212
        ic=278
c      pion0 + sigma0bar to k+ + pbar
        goto 7125
        endif
       ik1=311
        ik2=-2112
        ic=279
c       pion0 + sigma0bar to k0 + nbar
7125     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif







c       pion+ + sigma0 
        if(kl.eq.211 .and. kl1.eq.3212)then
        if(isinel(98).eq.0)then
        si1=0.
        goto 2681
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(3322),1)
        si1=s1724(ss,ilo1,0,the)
2681     if(isinel(251).eq.0)then
        si2=0.
        goto 2682
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(2212),1)
       ik3=-311
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10.*fac
2682     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=321
        ik2=3322
        ic=98
c       pion+ + sigma0 to k+ + cacade0
        goto 1683
        endif
        ik1=-311
        ik2=2212
        ic=251
c       pion+ + sigma0 to k0- + p
1683     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif
c       pion- + sigma0-
        if(kl.eq.-211 .and. kl1.eq.-3212)then
        if(isinel(106).eq.0)then
        si1=0.
        goto 2781
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3322),1)
        si1=s1724(ss,ilo1,0,the)
2781     if(isinel(273).eq.0)then
        si2=0.
        goto 2782
        endif
        the=pmas(lucomp(311),1)+pmas(lucomp(-2212),1)
       ik3=311
        ik4=-2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10.*fac
2782     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-321
        ik2=-3322
        ic=106
c       pion- + sigma0- to k- + cacade0-
        goto 1783
        endif
        ik1=311
        ik2=-2212
        ic=273
c       pion- + sigma0- to k0 + pba
1783     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif

c       pion0 + sigma+
        if(kl.eq.111 .and. kl1.eq.3222)then
        if(isinel(101).eq.0)then
        si1=0.
        goto 3681
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(3322),1)
        si1=s1724(ss,ilo1,0,the)
3681     if(isinel(259).eq.0)then
        si2=0.
        goto 3682
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(2212),1)
       ik3=-311
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10.*fac
3682     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=321
        ik2=3322
        ic=101
c       pion0 + sigma+ to k+ + cascade0
        goto 7683
        endif
        ik1=-311
        ik2=2212
        ic=259
c       pion0 + sigma+ to k0- + p
7683     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif
c       pion0 + sigma+bar
        if(kl.eq.111 .and. kl1.eq.-3222)then
        if(isinel(108).eq.0)then
        si1=0.
        goto 3781
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3322),1)
        si1=s1724(ss,ilo1,0,the)
3781     if(isinel(276).eq.0)then
        si2=0.
        goto 3782
        endif
        the=pmas(lucomp(311),1)+pmas(lucomp(-2212),1)
       ik3=311
        ik4=-2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10.*fac
3782     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-321
        ik2=-3322
        ic=108
c       pion0 + sigma+bar to k- + cascade0bar
        goto 6783
        endif
        ik1=311
        ik2=-2212
        ic=276
c       pion0 + sigma+bar to k0 + pbar
6783     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif

c       pion0 + sigma-
        if(kl.eq.111 .and. kl1.eq.3112)then
        if(isinel(40).eq.0)then
        si1=0.
        goto 4681
        endif
        the=pmas(lucomp(311),1)+pmas(lucomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
4681     if(isinel(260).eq.0)then
        si2=0.
        goto 4682
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(2112),1)
       ik3=-321
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10.*fac
4682     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=311
        ik2=3312
        ic=40
c       pion0 + sigma- to k0+cascade-
        goto 2683
        endif
        ik1=-321
        ik2=2112
        ic=260
c     pion0 + sigma- to k- + n
2683     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif
c       pion0 + sigma-bar
        if(kl.eq.111 .and. kl1.eq.-3112)then
        if(isinel(47).eq.0)then
        si1=0.
        goto 7681
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
7681     if(isinel(277).eq.0)then
        si2=0.
        goto 7682
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(-2112),1)
       ik3=321
        ik4=-2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &  1.,1.,0.,0.5,0.5,0.5,0.,0.5,1.)
        si2=s1724(ss,ilo2,0,the)*10.*fac
7682     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-311
        ik2=-3312
        ic=47
c       pion0 + sigma-bar to k0- +cascade-bar
        goto 7688
        endif
        ik1=321
        ik2=-2112
        ic=277
c     pion0 + sigma-bar to k+ + nbar
7688     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cspin/10.
        goto 10
        endif





 
c       k- + p
        if(kl.eq.-321 .and. kl1.eq.2212)then
        if(isinel(49).eq.0)then
        si1=0.
        goto 407
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(3122),1)
        si1=s1724(ss,ilo1,0,the)
c       cross section of k + n to pion + y is assumed to be equal to 
c        ten times of pion + n to k + y
407     if(isinel(50).eq.0)then
        si2=0.
        goto 408
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(3212),1)
        si2=s1724(ss,ilo2,0,the)
408     if(isinel(51).eq.0)then
        si3=0.
        goto 411
        endif
        the=pmas(lucomp(-211),1)+pmas(lucomp(3222),1)
        si3=s1724(ss,ilo3,0,the)
411     if(isinel(52).eq.0)then
        si4=0.
        goto 431
        endif
        the=pmas(lucomp(211),1)+pmas(lucomp(3112),1)
        si4=s1724(ss,ilo4,0,the)
431     if(isinel(53).eq.0)then
        si5=0.
        goto 412
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(3312),1)
        si5=s1724(ss,ilo5,0,the)/10.
c       cross section of kaon + n to kaon + cascade is assumed to be 
c        equal to pion + n to kaon + y,but take the different of
c        threshold energy into account

412     if(isinel(125).eq.0)then
        si6=0.
        goto 725
        endif
        the=pmas(lucomp(311),1)+pmas(lucomp(3322),1)
          si6=s1724(ss,ilo6,0,the)/10.

725     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0 
     c    .and. ilo5.eq.0 .and. ilo6.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6 .and. si5.eq.1.e-6 .and. si6.eq.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=3122
        ic=49
c       k- + p to pion0 + lambda
        goto 414
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=111
        ik2=3212
        ic=50
c       k- + p to pion0 + sigma0
        goto 414
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then 
        ik1=-211
        ik2=3222
        ic=51
c       k- + p to pion- + sigma+
        goto 414
        endif
        if(rlu1.gt.s3 .and. rlu1.le.s4)then 
        ik1=211
        ik2=3112
        ic=52
c       k- + p to pion+ + sigma-
        goto 414
        endif
        if(rlu1.gt.s4 .and. rlu1.le.s5)then
        ik1=321
        ik2=3312
        ic=53
c       k- + p to k+ + cascade-
        goto 414
        endif
        ik1=311
        ik2=3322
        ic=125
c       k- + p to k0 + cascade0
414     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn
        goto 10
        endif

c       k- + n
        if(kl.eq.-321 .and. kl1.eq.2112)then
        if(isinel(54).eq.0)then
        si1=0.
        goto 415
        endif
        the=pmas(lucomp(-211),1)+pmas(lucomp(3212),1)
        si1=s1724(ss,ilo1,0,the)
415     if(isinel(55).eq.0)then
        si2=0.
        goto 416
        endif
        the=pmas(lucomp(-211),1)+pmas(lucomp(3122),1)
        si2=s1724(ss,ilo2,0,the)
416     if(isinel(56).eq.0)then
        si3=0.
        goto 417
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(3112),1)
        si3=s1724(ss,ilo3,0,the)
417     if(isinel(57).eq.0)then
        si4=0.
        goto 432
        endif
        the=pmas(lucomp(311),1)+pmas(lucomp(3312),1)
        si4=s1724(ss,ilo4,0,the)/10.

432     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     $   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=3212
        ic=54
c       k- + n to pion- + sigma0
        goto 418
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-211
        ik2=3122
        ic=55
c       k- + n to pion- + lambda
        goto 418
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=111
        ik2=3112
        ic=56
c       k- + n to pion0 + sigma-
        goto 418
        endif
        ik1=311
        ik2=3312
        ic=57
c       k- + n to k0 + cascade-
418     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn
        goto 10
        endif

c       k0- + p
        if(kl.eq.-311 .and. kl1.eq.2212)then
        if(isinel(58).eq.0)then
        si1=0.
        goto 419
        endif
        the=pmas(lucomp(211),1)+pmas(lucomp(3122),1)
        si1=s1724(ss,ilo1,0,the)
419     if(isinel(59).eq.0)then
        si2=0.
        goto 420
        endif
        the=pmas(lucomp(211),1)+pmas(lucomp(3212),1)
        si2=s1724(ss,ilo2,0,the)
420     if(isinel(60).eq.0)then
        si3=0.
        goto 421
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(3222),1)
        si3=s1724(ss,ilo3,0,the)
421     if(isinel(126).eq.0)then
        si4=0.
        goto 726
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(3322),1)
        si4=s1724(ss,ilo4,0,the)/10.

726     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.si4.lt.
     &  1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        sit=si13+si4
        s1=si1/sit
        s2=si12/sit
        s3=si13/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=3122
        ic=58
c       k0- + p to pion+ + lambda
        goto 422
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=211
        ik2=3212
        ic=59
c       k0- + p to pion+ + sigma0
        goto 422
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=111
        ik2=3222
        ic=60
c       k0- + p to pion0 + sigma+
        goto 422
        endif
        ik1=321
        ik2=3322
        ic=126
c       k0- + p to k+ + cascade0
422     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0- + n 
        if(kl.eq.-311 .and. kl1.eq.2112)then
        if(isinel(61).eq.0)then
        si1=0.
        goto 423
        endif
        the=pmas(lucomp(211),1)+pmas(lucomp(3112),1)
        si1=s1724(ss,ilo1,0,the)
423     if(isinel(62).eq.0)then
        si2=0.
        goto 424
        endif
        the=pmas(lucomp(-211),1)+pmas(lucomp(3222),1)
        si2=s1724(ss,ilo2,0,the)
424     if(isinel(63).eq.0)then
        si3=0.
        goto 425
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(3212),1)
        si3=s1724(ss,ilo3,0,the)
425     if(isinel(64).eq.0)then
        si4=0.
        goto 426
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(3122),1)
        si4=s1724(ss,ilo4,0,the)
426     if(isinel(65).eq.0)then
        si5=0.
        goto 433
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(3312),1)
        si5=s1724(ss,ilo5,0,the)/10.
433     if(isinel(127).eq.0)then
        si6=0.
        goto 729
        endif
        the=pmas(lucomp(311),1)+pmas(lucomp(3322),1)
        si6=s1724(ss,ilo6,0,the)/10.
729     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0 
     c   .and. ilo5.eq.0 .and. ilo6.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6.and.si5.lt.1.e-6.and.si6.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=211
        ik2=3112
        ic=61
c       k0- + n to pion+ + sigma-
        goto 427
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-211
        ik2=3222
        ic=62
c       k0- + n to pion- + sigma+
        goto 427
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then 
        ik1=111
        ik2=3212
        ic=63
c       k0- + n to pion0 + sigma0
        goto 427
        endif
        if(rlu1.gt.s3 .and. rlu1.le.s4)then 
        ik1=111
        ik2=3122
        ic=64
c       k0- + n to pion0 + lambda
        goto 427
        endif
        if(rlu1.gt.s4 .and. rlu1.le.s5)then
        ik1=321
        ik2=3312
        ic=65
c       k0- + n to k+ + cascade-
        goto 427
        endif
        ik1=311
        ik2=3322
        ic=127
c       k0- + n to k0 + cascade0
427     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn
        goto 10
        endif

c       k+ + p-
        if(kl.eq.321 .and. kl1.eq.-2212)then
        if(isinel(66).eq.0)then
        si1=0.
        goto 510
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)
c       cross section of k + n- to pion + y- is assumed to be equal to 
c        ten times of pion + n to k + y
510     if(isinel(67).eq.0)then
        si2=0.
        goto 511
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)
511     if(isinel(68).eq.0)then
        si3=0.
        goto 512
        endif
        the=pmas(lucomp(211),1)+pmas(lucomp(-3222),1)
        si3=s1724(ss,ilo3,0,the)
512     if(isinel(69).eq.0)then
        si4=0.
        goto 513
        endif
        the=pmas(lucomp(-211),1)+pmas(lucomp(-3112),1)
        si4=s1724(ss,ilo4,0,the)
513     if(isinel(70).eq.0)then
        si5=0.
        goto 514
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3312),1)
        si5=s1724(ss,ilo5,0,the)/10.

514     if(isinel(128).eq.0)then
        si6=0.
        goto 730
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3322),1)

        si6=s1724(ss,ilo6,0,the)/10.

730     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0 
     c    .and. ilo5.eq.0 .and. ilo6.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6 .and. si5.eq.1.e-6.and.si6.eq.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=-3122
        ic=66
c       k+ + p- to pion0 + lambda-
        goto 515
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=111
        ik2=-3212
        ic=67
c       k+ + p- to pion0 + sigma0-
        goto 515
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then 
        ik1=211
        ik2=-3222
        ic=68
c       k+ + p- to pion+ + sigma+ba
        goto 515
        endif
        if(rlu1.gt.s3 .and. rlu1.le.s4)then 
        ik1=-211
        ik2=-3112
        ic=69
c       k+ + p- to pion- + sigma-ba
        goto 515
        endif
        if(rlu1.gt.s4 .and. rlu1.le.s5)then
        ik1=-321
        ik2=-3312
        ic=70
c       k+ + p- to k- + cascade-ba
        goto 515
        endif
        ik1=-311
        ik2=-3322
        ic=128
c       k+ + p- to k0- + cascade0ba
515     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn
        goto 10
        endif

c       k+ + n-
        if(kl.eq.321 .and. kl1.eq.-2112)then
        if(isinel(71).eq.0)then
        si1=0.
        goto 516
        endif
        the=pmas(lucomp(211),1)+pmas(lucomp(-3212),1)
        si1=s1724(ss,ilo1,0,the)
516     if(isinel(72).eq.0)then
        si2=0.
        goto 517
        endif
        the=pmas(lucomp(211),1)+pmas(lucomp(-3122),1)
        si2=s1724(ss,ilo2,0,the)
517     if(isinel(73).eq.0)then
        si3=0.
        goto 518
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(-3112),1)
        si3=s1724(ss,ilo3,0,the)
518     if(isinel(74).eq.0)then
        si4=0.
        goto 519
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3312),1)
        si4=s1724(ss,ilo4,0,the)/10.

519     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     $   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=-3212
        ic=71
c       k+ + n- to pion+ + sigma0-
        goto 520
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=211
        ik2=-3122
        ic=72
c       k+ + n- to pion+ + lambda-
        goto 520
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=111
        ik2=-3112
        ic=73
c       k+ + n- to pion0 + sigma-ba
        goto 520
        endif
        ik1=-311
        ik2=-3312
        ic=74
c       k+ + n- to k0- + cascade-ba
520     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn
        goto 10
        endif

c       k0 + p-
        if(kl.eq.311 .and. kl1.eq.-2212)then
        if(isinel(75).eq.0)then
        si1=0.
        goto 521
        endif
        the=pmas(lucomp(-211),1)+pmas(lucomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)
521     if(isinel(76).eq.0)then
        si2=0.
        goto 522
        endif
        the=pmas(lucomp(-211),1)+pmas(lucomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)
522     if(isinel(77).eq.0)then
        si3=0.
        goto 523
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(-3222),1)
        si3=s1724(ss,ilo3,0,the)
523     if(isinel(129).eq.0)then
        si4=0.
        goto 731
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3322),1)
        si4=s1724(ss,ilo4,0,the)/10.

731     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &   si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        sit=si13+si4
        s1=si1/sit
        s2=si12/sit
        s3=si13/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=-3122
        ic=75
c       k0 + p- to pion- + lambda-
        goto 524
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-211
        ik2=-3212
        ic=76
c       k0 + p- to pion- + sigma0-
        goto 524
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=111
        ik2=-3222
        ic=77
c       k0 + p- to pion0 + sigma+ba
        goto 524
        endif
        ik1=-321
        ik2=-3322
        ic=129
c       k0 + p- to k- + cascade0-
524     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0 + n-
        if(kl.eq.311 .and. kl1.eq.-2112)then
        if(isinel(78).eq.0)then
        si1=0.
        goto 525
        endif
        the=pmas(lucomp(-211),1)+pmas(lucomp(-3112),1)
        si1=s1724(ss,ilo1,0,the)
525     if(isinel(79).eq.0)then
        si2=0.
        goto 526
        endif
        the=pmas(lucomp(211),1)+pmas(lucomp(-3222),1)
        si2=s1724(ss,ilo2,0,the)
526     if(isinel(80).eq.0)then
        si3=0.
        goto 527
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(-3212),1)
        si3=s1724(ss,ilo3,0,the)
527     if(isinel(81).eq.0)then
        si4=0.
        goto 528
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(-3122),1)
        si4=s1724(ss,ilo4,0,the)
528     if(isinel(82).eq.0)then
        si5=0.
        goto 529
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3312),1)
        si5=s1724(ss,ilo5,0,the)/10.

        if(isinel(130).eq.0)then
        si6=0.
        goto 529
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3322),1)
        si6=s1724(ss,ilo6,0,the)/10.
   
529     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0 
     c   .and. ilo5.eq.0.and. ilo6.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6.and.si5.lt.1.e-6.and.si6.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=-211
        ik2=-3112
        ic=78
c       k0 + n- to pion- + sigma-ba
        goto 530
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=211
        ik2=-3222
        ic=79
c       k0 + n- to pion+ + sigma+ba
        goto 530
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then 
        ik1=111
        ik2=-3212
        ic=80
c       k0 + n- to pion0 + sigma0-
        goto 530
        endif
        if(rlu1.gt.s3 .and. rlu1.le.s4)then 
        ik1=111
        ik2=-3122
        ic=81
c       k0 + n- to pion0 + lambda-
        goto 530
        endif
        if(rlu1.gt.s4 .and. rlu1.le.s5)then
        ik1=-321
        ik2=-3312
        ic=82
c       k0 + n- to k- + cascade-ba
        goto 530
        endif
        ik1=-311
        ik2=-3322
        ic=130
c       k0 + n- to k0- + cascade0-
530     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn
        goto 10
        endif

c       k- + lambda 
        if(kl.eq.-321 .and. kl1.eq.3122)then
        if(isinel(83).eq.0)then
        si1=0.
        goto 732
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
c       cross section of k + y to pion + cascade is assumed to be equal to 
c        k + n to pion + y
732     if(isinel(113).eq.0)then
        si2=0.
        goto 733
        endif
        the=pmas(lucomp(-211),1)+pmas(lucomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
733     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=111
        ik2=3312
        ic=83
c       k- + lambda to pion0 + cascade-
        goto 734
        endif
        ik1=-211
        ik2=3322
        ic=113
c       k- + lambda to pion- + cascade0
734     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k- + sigma+
        if(kl.eq.-321 .and. kl1.eq.3222)then
        if(isinel(84).eq.0)then
        si1=0.
        goto 735
        endif
        the=pmas(lucomp(211),1)+pmas(lucomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
735     if(isinel(111).eq.0)then
        si2=0.
        goto 736
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
736     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=3312
        ic=84
c       k- + sigma+ to pion+ + cascade-
        goto 737
        endif
        ik1=111
        ik2=3322
        ic=111
c       k- + sigma+ to pion0 + cascade0
737     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k- + sigma- to pion- + cascade-
        if(kl.eq.-321 .and. kl1.eq.3112)then
        if(isinel(85).eq.0)goto 13
        the=pmas(lucomp(-211),1)+pmas(lucomp(3312),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=3312
        lc(icp,5)=85
        goto 10
        endif

c       k- + sigma0 
        if(kl.eq.-321 .and. kl1.eq.3212)then
        if(isinel(86).eq.0)then
        si1=0.
        goto 738
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
738     if(isinel(112).eq.0)then
        si2=0.
        goto 739
        endif
        the=pmas(lucomp(-211),1)+pmas(lucomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
739     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=111
        ik2=3312
        ic=86
c       k- + sigma0 to pion0 + cascade-
        goto 740
        endif
        ik1=-211
        ik2=3322
        ic=112
740     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0- + lambda
        if(kl.eq.-311 .and. kl1.eq.3122)then
        if(isinel(87).eq.0)then
        si1=0.
        goto 741
        endif
        the=pmas(lucomp(211),1)+pmas(lucomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
741     if(isinel(117).eq.0)then
        si2=0.
        goto 742
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
742     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=3312
        ic=87
c       k0- + lambda to pion+ + cascade-
        goto 743
        endif
        ik1=111
        ik2=3322
        ic=117
c       k0- + lambda to pion0 + cascade0
743     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0- + sigma0
        if(kl.eq.-311 .and. kl1.eq.3212)then
        if(isinel(88).eq.0)then
        si1=0.
        goto 744
        endif
        the=pmas(lucomp(211),1)+pmas(lucomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
744     if(isinel(116).eq.0)then
        si2=0.
        goto 745
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
745     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=3312
        ic=88
c       k0- + sigma0 to pion+ + cascade-
        goto 746
        endif
        ik1=111
        ik2=3322
        ic=116
c       k0- + sigma0 to pion0 + cascade0
746     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0- + sigma- 
        if(kl.eq.-311 .and. kl1.eq.3112)then
        if(isinel(89).eq.0)then
        si1=0.
        goto 747
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(3312),1)
        si1=s1724(ss,ilo1,0,the)
747     if(isinel(115).eq.0)then
        si2=0.
        goto 748
        endif
        the=pmas(lucomp(-211),1)+pmas(lucomp(3322),1)
        si2=s1724(ss,ilo2,0,the)
748     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=111
        ik2=3312
        ic=89
c       k0- + sigma- to pion0 + cascade-
        goto 749
        endif
        ik1=-211
        ik2=3322
        ic=115
c       k0- + sigma- to pion- + cascade0
749     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k+ + lambda- 
        if(kl.eq.321 .and. kl1.eq.-3122)then
        if(isinel(90).eq.0)then
        si1=0.
        goto 750
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
750     if(isinel(120).eq.0)then
        si2=0.
        goto 751
        endif
        the=pmas(lucomp(211),1)+pmas(lucomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
751     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=111
        ik2=-3312
        ic=90
c       k+ + lambda- to pion0 + cascade-ba
        goto 752
        endif
        ik1=211
        ik2=-3322
        ic=120
c       k+ + lambda- to pion+ + cascade0-
752     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k+ + sigma+ba 
        if(kl.eq.321 .and. kl1.eq.-3222)then
        if(isinel(91).eq.0)then
        si1=0.
        goto 753
        endif
        the=pmas(lucomp(-211),1)+pmas(lucomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
753     if(isinel(118).eq.0)then
        si2=0.
        goto 754
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
754     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=-3312
        ic=91
c       k+ + sigma+ba to pion- + cascade-ba
        goto 755
        endif
        ik1=111
        ik2=-3322
        ic=118
c       k+ + sigma+ba to pion0 + cascade0-
755     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k+ + sigma-ba to pion+ + cascade-ba
        if(kl.eq.321 .and. kl1.eq.-3112)then
        if(isinel(92).eq.0)goto 13
        the=pmas(lucomp(211),1)+pmas(lucomp(-3312),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=-3312
        lc(icp,5)=92
        goto 10
        endif

c       k+ + sigma0ba
        if(kl.eq.321 .and. kl1.eq.-3212)then
        if(isinel(93).eq.0)then
        si1=0.
        goto 756
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
756     if(isinel(119).eq.0)then
        si2=0.
        goto 757
        endif
        the=pmas(lucomp(211),1)+pmas(lucomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
757     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=111
        ik2=-3312
        ic=93
c       k+ + sigma0- to pion0 + cascade-ba
        goto 758
        endif
        ik1=211
        ik2=-3322
        ic=119
c       k+ + sigma0- to pion+ + cascade0-
758     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0+ + lambda-
        if(kl.eq.311 .and. kl1.eq.-3122)then
        if(isinel(94).eq.0)then
        si1=0.
        goto 759
        endif
        the=pmas(lucomp(-211),1)+pmas(lucomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
759     if(isinel(124).eq.0)then
        si2=0.
        goto 760
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
760     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=-3312
        ic=94
c       k0 + lambda- to pion- + cascade-ba
        goto 761
        endif
        ik1=111
        ik2=-3322
        ic=124
c       k0 + lambda- to pion0 + cascade0ba

761     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0 + sigma0-
        if(kl.eq.311 .and. kl1.eq.-3212)then
        if(isinel(95).eq.0)then
        si1=0.
        goto 762
        endif
        the=pmas(lucomp(-211),1)+pmas(lucomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
762     if(isinel(123).eq.0)then
        si2=0.
        goto 763
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
763     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=-3312
        ic=95
c       k0 + sigma0- to pion- + cascade-ba
        goto 764
        endif
        ik1=111
        ik2=-3322
        ic=123
c       k0 + sigma0- to pion0 + cascade0-
764     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0 + sigma-ba
        if(kl.eq.311 .and. kl1.eq.-3112)then
        if(isinel(96).eq.0)then
        si1=0.
        goto 765
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)
765     if(isinel(122).eq.0)then
        si2=0.
        goto 766
        endif
        the=pmas(lucomp(211),1)+pmas(lucomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)
766     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=111
        ik2=-3312
        ic=96
c       k0 + sigma-ba to pion0 + cascade-ba
        goto 767
        endif
        ik1=211
        ik2=-3322
        ic=122
c       k0 + sigma-ba to pion+ + cascade0ba
767     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn
        goto 10
        endif

c       k0- + sigma+ to pion+ + cascade0
        if(kl.eq.-311 .and. kl1.eq.3222)then
        if(isinel(114).eq.0)goto 13
        the=pmas(lucomp(211),1)+pmas(lucomp(3322),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=3322
        lc(icp,5)=114
        goto 10
        endif

c       k0 + sigma+bar to pion- + cascade0ba
        if(kl.eq.311 .and. kl1.eq.-3222)then
        if(isinel(121).eq.0)goto 13
        the=pmas(lucomp(-211),1)+pmas(lucomp(-3322),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=-3322
        lc(icp,5)=121
        goto 10
        endif

c       k+ + cascade-ba to pion+ + omiga-ba
        if(kl.eq.321 .and. kl1.eq.-3312)then
        if(isinel(143).eq.0)goto 13
        the=pmas(lucomp(211),1)+pmas(lucomp(-3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=-3334
        lc(icp,5)=143
        goto 10
        endif

c       k+ + cascade0- to pion0 + omiga-ba
        if(kl.eq.321 .and. kl1.eq.-3322)then
        if(isinel(145).eq.0)goto 13
        the=pmas(lucomp(111),1)+pmas(lucomp(-3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=111
        lc(icp,4)=-3334
        lc(icp,5)=145
        goto 10
        endif

c       k- + cascade- to pion- + omiga-
        if(kl.eq.-321 .and. kl1.eq.3312)then
        if(isinel(139).eq.0)goto 13
        the=pmas(lucomp(-211),1)+pmas(lucomp(3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=3334
        lc(icp,5)=139
        goto 10
        endif
       
c       k- + cascade0 to pion0 + omiga-
        if(kl.eq.-321 .and. kl1.eq.3322)then
        if(isinel(141).eq.0)goto 13
        the=pmas(lucomp(111),1)+pmas(lucomp(3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=111
        lc(icp,4)=3334
        lc(icp,5)=141
        goto 10
        endif

c       k0 + cascade-ba to pion0 + omiga-ba
        if(kl.eq.311 .and. kl1.eq.-3312)then
        if(isinel(144).eq.0)goto 13
        the=pmas(lucomp(111),1)+pmas(lucomp(-3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=111
        lc(icp,4)=-3334
        lc(icp,5)=144
        goto 10
        endif

c       k0 + cascade0- to pion- + omiga-ba
        if(kl.eq.311 .and. kl1.eq.-3322)then
        if(isinel(146).eq.0)goto 13
        the=pmas(lucomp(-211),1)+pmas(lucomp(-3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=-3334
        lc(icp,5)=146
        goto 10
        endif

c       k0- + cascade- to pion0 + omiga-
        if(kl.eq.-311 .and. kl1.eq.3312)then
        if(isinel(140).eq.0)goto 13
        the=pmas(lucomp(111),1)+pmas(lucomp(3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=111
        lc(icp,4)=3334
        lc(icp,5)=140
        goto 10
        endif

c       k0- + cascade0 to pion+ + omiga-
        if(kl.eq.-311 .and. kl1.eq.3322)then
        if(isinel(142).eq.0)goto 13
        the=pmas(lucomp(211),1)+pmas(lucomp(3334),1)
        tw(icp)=s1724(ss,ilo,0,the)/cskn
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=3334
        lc(icp,5)=142
        goto 10
        endif

c       follows are for reverse reactions
c       k+ + k-
        if((kl.eq.321 .and. kl1.eq.-321)
     c   .or.(kl.eq.-321 .and. kl1.eq.321))then
        if(isinel(201).eq.0)then
        fact1=0.
        goto 601
        endif
        ik3=211
        ik4=-211
        ic3=201
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  0.5,0.5,0.,0.,1.,1.,0.,0.,1.)
c k+ + k- to pi+ + pi-
        fact1=fac
        if(ilo1.eq.0)fact1=0.
601     if(isinel(202).eq.0)then
        fact2=0.
        goto 602
        endif
        ik5=111
        ik6=111
        ic5=202
        call srev(kl,kl1,ik5,ik6,ss,ilo2,fac,
     &  0.5,0.5,0.,0.,1.,1.,0.,0.,0.5)
c k+ + k- to pi0 + pi0
        fact2=fac
        if(ilo2.eq.0)fact2=0.
602     fact=fact1+fact2
        if(fact1.eq.0. .and. fact2.eq.0.)goto 13
c       if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=ic3
        if(rlu(1).gt.fact1/fact)then
        lc(icp,3)=ik5
        lc(icp,4)=ik6
        lc(icp,5)=ic5
        endif
        tw(icp)=fact*sig/cspipi
        goto 10
        endif

c       k+ + k0-
        if((kl.eq.321 .and. kl1.eq.-311)
     c   .or.(kl.eq.-311 .and. kl1.eq.321))then
        if(isinel(203).eq.0)goto 13
        ik3=211
        ik4=111
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  0.5,0.5,0.,0.,1.,1.,0.,0.,1.)
c k+ + k0- to pi+ + pi0
        if(ilo1.eq.0)goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=203
        tw(icp)=fac*sig/cspipi
        goto 10
        endif 

c       k- + k0
        if((kl.eq.-321 .and. kl1.eq.311)
     c   .or.(kl.eq.311 .and. kl1.eq.-321))then
        if(isinel(204).eq.0)goto 13
        ik3=-211
        ik4=111
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  0.5,0.5,0.,0.,1.,1.,0.,0.,1.)
c k- + k0 to pi- + pi0
        if(ilo1.eq.0)goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=204
        tw(icp)=fac*sig/cspipi
        goto 10
        endif 

c       k0 + k0-
        if((kl.eq.311 .and. kl1.eq.-311)
     c   .or.(kl.eq.-311 .and. kl1.eq.311))then
        if(isinel(205).eq.0)then
        fact1=0.
        goto 603
        endif
        ik3=211
        ik4=-211
        ic3=205
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  0.5,0.5,0.,0.,1.,1.,0.,0.,1.)
c k0 + k0- to pi- + pi+
        fact1=fac
        if(ilo1.eq.0)fact1=0.
603     if(isinel(206).eq.0)then
        fact2=0.
        goto 604
        endif
        ik5=111
        ik6=111
        ic5=206
        call srev(kl,kl1,ik5,ik6,ss,ilo2,fac,
     &  0.5,0.5,0.,0.,1.,1.,0.,0.,0.5)
c k0 + k0- to pi0 + pi0
        fact2=fac
        if(ilo2.eq.0)fact2=0.
604     fact=fact1+fact2
        if(fact1.eq.0. .and. fact2.eq.0.)goto 13
c       if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        lc(icp,3)=ik3
        lc(icp,4)=ik4
        lc(icp,5)=ic3
        if(rlu(1).gt.fact1/fact)then
        lc(icp,3)=ik5
        lc(icp,4)=ik6
        lc(icp,5)=ic5
        endif
        tw(icp)=fact*sig/cspipi
        goto 10
        endif
  

c       k+ + sigma+ to pion+ + p
        if(kl.eq.321 .and. kl1.eq.3222)then
c       3222 is the flavor code of sigma+
        if(isinel(207).eq.0)goto 13
        the=pmas(lucomp(211),1)+pmas(lucomp(2212),1)
        ww=s1713(ss,ilo,0,the)/cskn/10.
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=2212
        lc(icp,5)=207
        ik3=211
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif


c       k+ + sigma-
        if(kl.eq.321 .and. kl1.eq.3112)then
        if(isinel(208).eq.0)then
        si1=0.
        goto 606
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(-211),1)+pmas(lucomp(2212),1)
        si1=s0715(ss,ilo1,0,the)*fac
606     if(isinel(209).eq.0)then
        si2=0.
        goto 607
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(2112),1)
        si2=s2325(ss,ilo2,0,the)*fac
607     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=2212
        ic=208
c       k+ + sigma- to pion- + p
        goto 609
        endif
        ik1=111
        ik2=2112
        ic=209
c       k+ + sigma- to pion0 + n
609     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

   
   
c       k+ + sigma0
        if(kl.eq.321 .and. kl1.eq.3212)then
        if(isinel(210).eq.0)then
        si1=0.
        goto 610
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(211),1)+pmas(lucomp(2112),1)
        si1=s1724(ss,ilo1,0,the)*fac
610     if(isinel(211).eq.0)then
        si2=0.
        goto 611
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(2212),1)
        si2=s2314(ss,ilo2,0,the)*fac
611     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=2112
        ic=210
c       k+ + sigma0 to pion+ + n
        goto 612
        endif
        ik1=111
        ik2=2212
        ic=211
c       k+ + sigma0 to pion0 + p
612     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif


c       k+ + lambda0
        if(kl.eq.321 .and. kl1.eq.3122)then
        if(isinel(212).eq.0)then
        si1=0.
        goto 613
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(211),1)+pmas(lucomp(2112),1)
        si1=s1727(ss,ilo1,0,the)*fac
613     if(isinel(213).eq.0)then
        si2=0.
        goto 614
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(2212),1)
        si2=s2317(ss,ilo2,0,the)*fac
614     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=2112
        ic=212
c       k+ + lambda0 to pion+ + n
        goto 615
        endif
        ik1=111
        ik2=2212
        ic=213
c       k+ + lambda0 to pion0 + p
615     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif


c       k0 + sigma+
        if(kl.eq.311 .and. kl1.eq.3222)then
        if(isinel(214).eq.0)then
        si1=0.
        goto 616
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(211),1)+pmas(lucomp(2112),1)
        si1=s1724(ss,ilo1,0,the)*fac
616     if(isinel(215).eq.0)then
        si2=0.
        goto 617
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(2212),1)
        si2=s1724(ss,ilo2,0,the)*fac
617     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=2112
        ic=214
c       k+ + sigma+ to pion+ + n
        goto 618
        endif
        ik1=111
        ik2=2212
        ic=215
c       k+ + sigma+ to pion0 + p
618     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k0 + sigma- to pion- +n 
        if(kl.eq.311 .and. kl1.eq.3112)then
c       3112 is the flavor code of sigma-
        if(isinel(216).eq.0)goto 13
        the=pmas(lucomp(-211),1)+pmas(lucomp(2112),1)
        ww=s1724(ss,ilo,0,the)/cskn/10.
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=2112
        lc(icp,5)=216
        ik3=-211
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       k0 + sigma0
        if(kl.eq.311 .and. kl1.eq.3212)then
        if(isinel(217).eq.0)then
        si1=0.
        goto 619
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(-211),1)+pmas(lucomp(2212),1)
        si1=s07123(ss,ilo1,0,the)*fac
619     if(isinel(218).eq.0)then
        si2=0.
        goto 620
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(2112),1)
        si2=s1724(ss,ilo2,0,the)*fac
620     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=2212
        ic=217
c       k0 + sigma0 to pion- + p
        goto 621
        endif
        ik1=111
        ik2=2112
        ic=218
c       k0 + sigma0 to pion0 + n
621     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k0 + lambda0
        if(kl.eq.311 .and. kl1.eq.3122)then
        if(isinel(219).eq.0)then
        si1=0.
        goto 622
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(-211),1)+pmas(lucomp(2212),1)
        si1=s07122(ss,ilo1,0,the)*fac
622     if(isinel(220).eq.0)then
        si2=0.
        goto 623
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(2112),1)
        si2=s1724(ss,ilo2,0,the)*fac
623     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=2212
        ic=219
c       k0 + lambda0 to pion- + p
        goto 624
        endif
        ik1=111
        ik2=2112
        ic=220
c       k0 + lambda0 to pion0 + n
624     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k- + sigma+ba to pion- + pba
        if(kl.eq.-321 .and. kl1.eq.-3222)then
        if(isinel(221).eq.0)goto 13
        the=pmas(lucomp(-211),1)+pmas(lucomp(-2212),1)
        ww=s1724(ss,ilo,0,the)/cskn/10.
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=-2212
        lc(icp,5)=221
        ik3=-211
        ik4=-2212
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       k- + sigma-ba
        if(kl.eq.-321 .and. kl1.eq.-3112)then
        if(isinel(222).eq.0)then
        si1=0.
        goto 625
        endif
        ik1=211
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(211),1)+pmas(lucomp(-2212),1)
        si1=s1724(ss,ilo1,0,the)*fac
625     if(isinel(223).eq.0)then
        si2=0.
        goto 626
        endif
        ik1=111
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(-2112),1)
        si2=s1724(ss,ilo2,0,the)*fac
626     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=-2212
        ic=222
c       k- + sigma-ba to pion+ + pba
        goto 627
        endif
        ik1=111
        ik2=-2112
        ic=223
c       k- + sigma-ba to pion0 + nba
627     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif


c       k- + sigma0ba
        if(kl.eq.-321 .and. kl1.eq.-3212)then
        if(isinel(224).eq.0)then
        si1=0.
        goto 628
        endif
        ik1=-211
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(-211),1)+pmas(lucomp(-2112),1)
        si1=s1724(ss,ilo1,0,the)*fac
628     if(isinel(225).eq.0)then
        si2=0.
        goto 629
        endif
        ik1=111
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(-2212),1)
        si2=s1724(ss,ilo2,0,the)*fac
629     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=-2112
        ic=224
c       k- + sigma0ba to pion- + nba
        goto 630
        endif
        ik1=111
        ik2=-2212
        ic=225
c       k- + sigma0ba to pion0 + pba
630     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k- + lambda0-
        if(kl.eq.-321 .and. kl1.eq.-3122)then
        if(isinel(226).eq.0)then
        si1=0.
        goto 631
        endif
        ik1=-211
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(-211),1)+pmas(lucomp(-2112),1)
        si1=s1724(ss,ilo1,0,the)*fac
631     if(isinel(227).eq.0)then
        si2=0.
        goto 632
        endif
        ik1=111
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(-2212),1)
        si2=s1724(ss,ilo2,0,the)*fac
632     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=-2112
        ic=226
c       k- + lambda0ba to pion- + nba
        goto 633
        endif
        ik1=111
        ik2=-2212
        ic=227
c       k- + lambda0ba to pion0 + pba
633     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k0- + sigma+ba
        if(kl.eq.-311 .and. kl1.eq.-3222)then
        if(isinel(228).eq.0)then
        si1=0.
        goto 634
        endif
        ik1=-211
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(-211),1)+pmas(lucomp(-2112),1)
        si1=s1724(ss,ilo1,0,the)*fac
634     if(isinel(229).eq.0)then
        si2=0.
        goto 635
        endif
        ik1=111
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(-2212),1)
        si2=s1724(ss,ilo2,0,the)*fac
635     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=-2112
        ic=228
c       k0ba + sigma+ba to pion- + nba
        goto 636
        endif
        ik1=111
        ik2=-2212
        ic=229
c       k0ba + sigma+ba to pion0 + pba
636     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k0- + sigma-ba to pion+ + nba
        if(kl.eq.-311 .and. kl1.eq.-3112)then
        if(isinel(230).eq.0)goto 13
        the=pmas(lucomp(211),1)+pmas(lucomp(-2112),1)
        ww=s1724(ss,ilo,0,the)/cskn/10.
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=-2112
        lc(icp,5)=230
        ik3=211
        ik4=-2112
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       k0- + sigma0-
        if(kl.eq.-311 .and. kl1.eq.-3212)then
        if(isinel(231).eq.0)then
        si1=0.
        goto 637
        endif
        ik1=211
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(211),1)+pmas(lucomp(-2212),1)
        si1=s1724(ss,ilo1,0,the)*fac
637     if(isinel(232).eq.0)then
        si2=0.
        goto 638
        endif
        ik1=111
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,1.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(-2112),1)
        si2=s1724(ss,ilo2,0,the)*fac
638     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=-2212
        ic=231
c       k0- + sigma0- to pion+ + pba
        goto 639
        endif
        ik1=111
        ik2=-2112
        ic=232
c       k0- + sigma0- to pion0 + nba
639     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k0- + lambda0-
        if(kl.eq.-311 .and. kl1.eq.-3122)then
        if(isinel(233).eq.0)then
        si1=0.
        goto 640
        endif
        ik1=211
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(211),1)+pmas(lucomp(-2212),1)
        si1=s1724(ss,ilo1,0,the)*fac
640     if(isinel(234).eq.0)then
        si2=0.
        goto 641
        endif
        ik1=111
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,0.,0.,0.5,1.,0.5,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(-2112),1)
        si2=s1724(ss,ilo2,0,the)*fac
641     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=-2212
        ic=233
c       k0- + lambda0- to pion+ + pba
        goto 642
        endif
        ik1=111
        ik2=-2112
        ic=234
c       k0- + lambda0- to pion0 + nba
642     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k+ + cascade-
        if(kl.eq.321 .and. kl1.eq.3312)then
        if(isinel(235).eq.0)then
        si1=0.
        goto 699
        endif
        ik1=211
        ik2=3112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(lucomp(211),1)+pmas(lucomp(3112),1)
        si1=s1724(ss,ilo1,0,the)*fac
699     if(isinel(236).eq.0)then
        si2=0.
        goto 643
        endif
        ik1=-211
        ik2=3222
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(lucomp(-211),1)+pmas(lucomp(3222),1)
        si2=s1724(ss,ilo2,0,the)*fac
643     if(isinel(237).eq.0)then
        si3=0.
        goto 644
        endif
        ik1=111
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(3122),1)
        si3=s1724(ss,ilo3,0,the)*fac
644     if(isinel(238).eq.0)then
        si4=0.
        goto 645
        endif
        ik1=111
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(3212),1)
        si4=s1724(ss,ilo4,0,the)*fac
645     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0
     c   )goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=211
        ik2=3112
        ic=235
c       k+ + cascade- to pion+ + sigma-
        goto 646
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-211
        ik2=3222
        ic=236
c       k+ + cascade- to pion- + sigma+
        goto 646
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=111
        ik2=3122
        ic=237
c       k+ + cascade- to pion0 + lambda0
        goto 646
        endif
        ik1=111
        ik2=3212
        ic=238
c       k+ + cascade- to pion0 + sigma0
646     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn/10.
        goto 10
        endif

c       k- + cascade-ba
        if(kl.eq.-321 .and. kl1.eq.-3312)then
        if(isinel(239).eq.0)then
        si1=0.
        goto 647
        endif
        ik1=211
        ik2=-3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(lucomp(211),1)+pmas(lucomp(-3222),1)
        si1=s1724(ss,ilo1,0,the)*fac
647     if(isinel(240).eq.0)then
        si2=0.
        goto 648
        endif
        ik1=-211
        ik2=-3112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(lucomp(-211),1)+pmas(lucomp(-3112),1)
        si2=s1724(ss,ilo2,0,the)*fac
648     if(isinel(241).eq.0)then
        si3=0.
        goto 649
        endif
        ik1=111
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(-3122),1)
        si3=s1724(ss,ilo3,0,the)*fac
649     if(isinel(242).eq.0)then
        si4=0.
        goto 650
        endif
        ik1=111
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(-3212),1)
        si4=s1724(ss,ilo4,0,the)*fac
650     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0
     c   )goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=211
        ik2=-3222
        ic=239
c       k- + cascade-ba to pion+ + sigma+ba
        goto 651
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-211
        ik2=-3112
        ic=240
c       k- + cascade-ba to pion- + sigma-ba
        goto 651
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=111
        ik2=-3122
        ic=241
c       k- + cascade-ba to pion0 + lambda0-
        goto 651
        endif
        ik1=111
        ik2=-3212
        ic=242
c       k- + cascade-ba to pion0 + sigma0-
651     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn/10.
        goto 10
        endif

c       k0 + cascade-
        if(kl.eq.311 .and. kl1.eq.3312)then
        if(isinel(243).eq.0)then
        si1=0.
        goto 652
        endif
        ik1=-211
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(lucomp(-211),1)+pmas(lucomp(3122),1)
        si1=s1724(ss,ilo1,0,the)*fac
652     if(isinel(244).eq.0)then
        si2=0.
        goto 653
        endif
        ik1=-211
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(lucomp(-211),1)+pmas(lucomp(3212),1)
        si2=s1724(ss,ilo2,0,the)*fac
653     if(isinel(245).eq.0)then
        si3=0.
        goto 654
        endif
        ik1=111
        ik2=3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(3112),1)
        si3=s1724(ss,ilo3,0,the)*fac
654     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        s1=si1/si13
        s2=si12/si13
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=-211
        ik2=3122
        ic=243
c       k0 + cascade- to pion- + lambda0
        goto 655
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-211
        ik2=3212
        ic=244
c       k0 + cascade- to pion- + sigma0
        goto 655
        endif
        ik1=111
        ik2=3112
        ic=245
c       k0 + cascade- to pion0 + sigma-
655     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si13/cskn/10.
        goto 10
        endif

c       k0- + cascade-ba
        if(kl.eq.-311 .and. kl1.eq.-3312)then
        if(isinel(246).eq.0)then
        si1=0.
        goto 656
        endif
        ik1=211
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(lucomp(211),1)+pmas(lucomp(-3122),1)
        si1=s1724(ss,ilo1,0,the)*fac
656     if(isinel(247).eq.0)then
        si2=0.
        goto 657
        endif
        ik1=211
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(lucomp(211),1)+pmas(lucomp(-3212),1)
        si2=s1724(ss,ilo2,0,the)*fac
657     if(isinel(248).eq.0)then
        si3=0.
        goto 658
        endif
        ik1=111
        ik2=-3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(-3112),1)
        si3=s1724(ss,ilo3,0,the)*fac
658     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        s1=si1/si13
        s2=si12/si13
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=211
        ik2=-3122
        ic=246
c       k0- + cascade-ba to pion+ + lambda0-
        goto 659
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=211
        ik2=-3212
        ic=247
c       k0- + cascade-ba to pion+ + sigma0-
        goto 659
        endif
        ik1=111
        ik2=-3112
        ic=248
c       k0- + cascade-ba to pion0 + sigma-ba
659     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si13/cskn/10.
        goto 10
        endif



c       pion+ + lambda0 
        if(kl.eq.211 .and. kl1.eq.3122)then
        if(isinel(99).eq.0)then
        sigma1=0.
        goto 2522
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(3322),1)
        sigma1=s1724(ss,ilo1,0,the)
c   c   pion+ + lambda0 to k+ + cascade
2522            if(isinel(252).eq.0)then
        sigma2=0.
        goto 2523
        endif        
        ik3=-311
             ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &  1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(-311),1)+pmas(lucomp(2212),1)
        sigma2=fac*s1724(ss,ilo2,0,the)*10.

c       cross section of     pion+ + lambda0 to k0- + p
2523    if(ilo1.eq.0.and.ilo2.eq.0) goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6)  goto 13
        ik1=321
        ik2=3322
        ic=99
c       pion+ + lambda0 to k+ + cascade0
        sigm12=sigma1+sigma2
        if(rlu(1).gt.sigma1/sigm12)then
        ik1=-311
        ik2=2212
        ic=252
c       cross section of     pion+ + lambda0 to k0- + p
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/cspin/10.
        goto 10
        endif 

c       pion- + lambda0-
        if(kl.eq.-211 .and. kl1.eq.-3122)then
        if(isinel(107).eq.0)then
        sigma1=0.
        goto 3522
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3322),1)
        sigma1=s1724(ss,ilo1,0,the)
c       pion- + lambda- to k- + cascade0-
3522            if(isinel(275).eq.0)then
        sigma2=0.
        goto 3523
        endif        
        ik3=311
             ik4=-2212
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &  1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(311),1)+pmas(lucomp(-2212),1)
        sigma2=fac*s1724(ss,ilo2,0,the)*10.
c       pion- + lambda0- to k0 + pba

3523    if(ilo1.eq.0.and.ilo2.eq.0) goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6)  goto 13
        ik1=-321
        ik2=-3322
        ic=107
c       pion- + lambda- to k- + cascade0-
        sigm12=sigma1+sigma2
        if(rlu(1).gt.sigma1/sigm12)then
        ik1=311
        ik2=-2212
        ic=275
c       pion- + lambda0- to k0 + pba
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/cspin/10.
        goto 10
        endif 







c       k+ + cascade- to k- + p
        if(kl.eq.321 .and. kl1.eq.3312)then
        if(isinel(253).eq.0)goto 13
        the=pmas(lucomp(-321),1)+pmas(lucomp(2212),1)
        ww=s1724(ss,ilo,0,the)/cskn/10.
        if(ilo.eq.0) goto 13
        lc(icp,3)=-321
        lc(icp,4)=2212
        lc(icp,5)=253
        ik3=-321
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif 



c       k0 + cascade- to k- + n
        if(kl.eq.311 .and. kl1.eq.3312)then
        if(isinel(257).eq.0)goto 13
        the=pmas(lucomp(-321),1)+pmas(lucomp(2112),1)
        ww=s1724(ss,ilo,0,the)/cskn/10.
        if(ilo.eq.0) goto 13
        lc(icp,3)=-321
        lc(icp,4)=2112
        lc(icp,5)=257
        ik3=-321
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif 

c       pion- + lambda0 
        if(kl.eq.-211 .and. kl1.eq.3122)then
        if(isinel(36).eq.0)then
        sigma1=0.
        goto 1522
        endif
        the=pmas(lucomp(311),1)+pmas(lucomp(3312),1)
        sigma1=s1724(ss,ilo1,0,the)
c       pion- + lambda to k0 + cascade-
1522            if(isinel(258).eq.0)then
        sigma2=0.
        goto 1523
        endif        
        ik3=-321
             ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &  1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(-321),1)+pmas(lucomp(2112),1)
        sigma2=fac*s1724(ss,ilo2,0,the)*10.

c       cross section of   pion- + lambda0 to k- + n
1523    if(ilo1.eq.0.and.ilo2.eq.0) goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6)  goto 13
        ik1=311
        ik2=3312
        ic=36
c       pion- + lambda to k0 + cascade-
        sigm12=sigma1+sigma2
        if(rlu(1).gt.sigma1/sigm12)then
        ik1=-321
        ik2=2112
        ic=258
c       cross section of pion- + lambda0 to k- + n
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/cspin/10.
        goto 10
        endif 
c       pion+ + lambda-
        if(kl.eq.211 .and. kl1.eq.-3122)then
        if(isinel(42).eq.0)then
        sigma1=0.
        goto 6522
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3312),1)
        sigma1=s1724(ss,ilo1,0,the)
c       pion+ + lambdaba to k0- + cascade-ba
6522            if(isinel(269).eq.0)then
        sigma2=0.
        goto 6523
        endif        
        ik3=321
             ik4=-2112
        call srev(kl,kl1,ik3,ik4,ss,ilo2,fac,
     &  1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(321),1)+pmas(lucomp(-2112),1)
        sigma2=fac*s1724(ss,ilo2,0,the)*10.

c         pion+ + lambda0- to k+ + nba
6523    if(ilo1.eq.0.and.ilo2.eq.0) goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6)  goto 13
        ik1=-311
        ik2=-3312
        ic=42
cpion+ + lambdaba to k0- + cascade-ba
        sigm12=sigma1+sigma2
        if(rlu(1).gt.sigma1/sigm12)then
        ik1=321
        ik2=-2112
        ic=269
c  pion+ + lambda0- to k+ + nba
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/cspin/10.
        goto 10
        endif 


c       pion0 + lambda0
        if(kl.eq.111 .and. kl1.eq.3122)then
        if(isinel(263).eq.0)then
        si1=0.
        goto 669
        endif
        ik1=-321
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(-321),1)+pmas(lucomp(2212),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
669     if(isinel(264).eq.0)then
        si2=0.
        goto 670
        endif
        ik1=-311
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(-311),1)+pmas(lucomp(2112),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
670      if(isinel(39).eq.0)then
        si3=0.
        goto 1670
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(3312),1)
        si3=s1724(ss,ilo3,0,the)
1670      if(isinel(103).eq.0)then
        si4=0.
        goto 1671
        endif
        the=pmas(lucomp(311),1)+pmas(lucomp(3322),1)
        si4=s1724(ss,ilo4,0,the)
1671    if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.
     &  and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
       si14=si13+si4
        s1=si1/si14
       s2=si12/si14
       s3=si13/si14
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-321
        ik2=2212
        ic=263
c       pion0 + lambda0 to k- + p
        goto 671
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-311
        ik2=2112
        ic=264
c       pion0 + lambda0 to k0- + n
        goto 671
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
      ik1=321
        ik2=3312
        ic=39
c       pion0 + lambda0 to  k+ + cascade-
        goto 671
        endif
      ik1=311
        ik2=3322
        ic=103
c       pion0 + lambda to k0 + cascade0
671     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif
c       pion0 + lambda0-
        if(kl.eq.111 .and. kl1.eq.-3122)then
        if(isinel(280).eq.0)then
        si1=0.
        goto 2669
        endif
        ik1=321
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(321),1)+pmas(lucomp(-2212),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
2669     if(isinel(281).eq.0)then
        si2=0.
        goto 2677
        endif
        ik1=311
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.,0.,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(311),1)+pmas(lucomp(-2112),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
2677      if(isinel(46).eq.0)then
        si3=0.
        goto 2670
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3312),1)
        si3=s1724(ss,ilo3,0,the)
c       pion0 + lambda- to k- + cascade-ba
2670      if(isinel(110).eq.0)then
        si4=0.
        goto 2671
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3322),1)
        si4=s1724(ss,ilo4,0,the)
2671    if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0.and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.
     &  and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
       si14=si13+si4
        s1=si1/si14
       s2=si12/si14
       s3=si13/si14
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=321
        ik2=-2212
        ic=280
c       pion0 + lambda0- to k+ + pbar
        goto 6711
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=311
        ik2=-2112
        ic=281
c       pion0 + lambda0- to k0 + nba
        goto 6711
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
      ik1=-321
        ik2=-3312
        ic=46
c       pion0 + lambda- to k- + cascade-ba
        goto 6711
        endif
      ik1=-311
        ik2=-3322
        ic=110
c          pion0 + lambda- to k0- + cascade0-
6711     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif
c       k+ + cascade- to k0- + n
        if(kl.eq.321 .and. kl1.eq.3312)then
        if(isinel(265).eq.0)goto 13
        the=pmas(lucomp(-311),1)+pmas(lucomp(2112),1)
        ww=s1724(ss,ilo,0,the)/cskn/10.
        if(ilo.eq.0) goto 13
        lc(icp,3)=-311
        lc(icp,4)=2112
        lc(icp,5)=265
        ik3=-311
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif 





c       k- + cascade-ba to k+ + pba
        if(kl.eq.-321 .and. kl1.eq.-3312)then
        if(isinel(270).eq.0)goto 13
        the=pmas(lucomp(321),1)+pmas(lucomp(-2212),1)
        ww=s1724(ss,ilo,0,the)/cskn/10.
        if(ilo.eq.0) goto 13
        lc(icp,3)=321
        lc(icp,4)=-2212
        lc(icp,5)=270
        ik3=321
        ik4=-2212
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif 


c       k0- + cascade-ba to k+ + nba
        if(kl.eq.-311 .and. kl1.eq.-3312)then
        if(isinel(274).eq.0)goto 13
        the=pmas(lucomp(321),1)+pmas(lucomp(-2112),1)
        ww=s1724(ss,ilo,0,the)/cskn/10.
        if(ilo.eq.0) goto 13
        lc(icp,3)=321
        lc(icp,4)=-2112
        lc(icp,5)=274
        ik3=321
        ik4=-2112
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif 

c       k- + cascade-ba to k0 + nba
        if(kl.eq.-321 .and. kl1.eq.-3312)then
        if(isinel(282).eq.0)goto 13
        the=pmas(lucomp(311),1)+pmas(lucomp(-2112),1)
        ww=s1724(ss,ilo,0,the)/cskn/10.
        if(ilo.eq.0) goto 13
        lc(icp,3)=311
        lc(icp,4)=-2112
        lc(icp,5)=282
        ik3=311
        ik4=-2112
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif 

c       pion+ + cascade-
        if(kl.eq.211 .and. kl1.eq.3312)then
        if(isinel(283).eq.0)then
        si1=0.
        goto 684
        endif
        ik1=-321
        ik2=3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(lucomp(-321),1)+pmas(lucomp(3222),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
684     if(isinel(284).eq.0)then
        si2=0.
        goto 685
        endif
        ik1=-311
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.,0.5,0.,0.5,0.5,0.,0.,0.5,1.)
        the=pmas(lucomp(-311),1)+pmas(lucomp(3122),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
685     if(isinel(285).eq.0)then
        si3=0.
        goto 686
        endif
        ik1=-311
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(lucomp(-311),1)+pmas(lucomp(3212),1)
        si3=10*s1724(ss,ilo3,0,the)*fac
686     if(isinel(131).eq.0)then
        si4=0.
        goto 1686
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(3334),1)
        si4=s1724(ss,ilo4,0,the)

1686    if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0.
     &  and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &   si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=-321
        ik2=3222
        ic=283
c       pion+ + cascade- to k- + sigma+
        goto 687
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-311
        ik2=3122
        ic=284
c       pion+ + cascade- to k0- + lambda0
        goto 687
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=-311
        ik2=3212
        ic=285
c       pion+ + cascade- to k0- + sigma0
        goto 687
        endif
       ik1=321
        ik2=3334
        ic=131
c       pion+ + cascade- to k+ + omiga-
687     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion- + cascade- to k- + sigma-
        if(kl.eq.-211 .and. kl1.eq.3312)then
        if(isinel(286).eq.0)goto 13
        the=pmas(lucomp(-321),1)+pmas(lucomp(3112),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo.eq.0) goto 13
        lc(icp,3)=-321
        lc(icp,4)=3112
        lc(icp,5)=286
        ik3=-321
        ik4=3112
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif 

c       pion0 + cascade-
        if(kl.eq.111 .and. kl1.eq.3312)then
        if(isinel(287).eq.0)then
        si1=0.
        goto 688
        endif
        ik1=-321
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.,0.5,0.,0.5,0.5,0.,0.,0.5,1.)
        the=pmas(lucomp(-321),1)+pmas(lucomp(3122),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
688     if(isinel(288).eq.0)then
        si2=0.
        goto 689
        endif
        ik1=-321
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(lucomp(-321),1)+pmas(lucomp(3212),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
689     if(isinel(289).eq.0)then
        si3=0.
        goto 690
        endif
        ik1=-311
        ik2=3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(lucomp(-311),1)+pmas(lucomp(3112),1)
        si3=10*s1724(ss,ilo3,0,the)*fac
690     if(isinel(132).eq.0)then
        si4=0.
        goto 1690
        endif
        the=pmas(lucomp(311),1)+pmas(lucomp(3334),1)            
        si4=s1724(ss,ilo4,0,the)
   
c       pion0 + cascade- to k0 + omiga-
1690    if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0.and.
     &  ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &  si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=-321
        ik2=3122
        ic=287
c       pion0 + cascade- to k- + lambda0
        goto 691
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-321
        ik2=3212
        ic=288
c       pion0 + cascade- to k- + sigma0
        goto 691
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=-311
        ik2=3112
        ic=289
c       pion0 + cascade- to k0- + sigma-
        goto 691
        endif
        ik1=311
        ik2=3334
        ic=132
c       pion0 + cascade- to k0 + omiga-
691     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion+ + cascade-ba to k+ + sigma-ba
        if(kl.eq.211 .and. kl1.eq.-3312)then
        if(isinel(290).eq.0)goto 13
        the=pmas(lucomp(321),1)+pmas(lucomp(-3112),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo.eq.0) goto 13
        lc(icp,3)=321
        lc(icp,4)=-3112
        lc(icp,5)=290
        ik3=321
        ik4=-3112
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)

        tw(icp)=ww*fac
        goto 10
        endif 

c       pion- + cascade-ba
        if(kl.eq.-211 .and. kl1.eq.-3312)then
        if(isinel(291).eq.0)then
        si1=0.
        goto 692
        endif
        ik1=321
        ik2=-3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(lucomp(321),1)+pmas(lucomp(-3222),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
692     if(isinel(292).eq.0)then
        si2=0.
        goto 693
        endif
        ik1=311
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.,0.5,0.,0.5,0.5,0.,0.,0.5,1.)
        the=pmas(lucomp(311),1)+pmas(lucomp(-3122),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
693     if(isinel(293).eq.0)then
        si3=0.
        goto 694
        endif
        ik1=311
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(lucomp(311),1)+pmas(lucomp(-3212),1)
        si3=10*s1724(ss,ilo3,0,the)*fac
694     if(isinel(133).eq.0)then
        si4=0.
        goto 1694
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3334),1)          
        si4=s1724(ss,ilo4,0,the)
1694    if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0.and.
     &  ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &  si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=321
        ik2=-3222
        ic=291
c       pion- + cascade-ba to k+ + sigma+-
        goto 695
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=311
        ik2=-3122
        ic=292
c       pion- + cascade-ba to k0 + lambda0-
        goto 695
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=311
        ik2=-3212
        ic=293
c       pion- + cascade-ba to k0 + sigma0-
        goto 695
        endif
        ik1=-321
        ik2=-3334
        ic=133
c       pion- + cascade-ba to k- + omiga-ba
695     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion0 + cascade-ba
        if(kl.eq.111 .and. kl1.eq.-3312)then
        if(isinel(294).eq.0)then
        si1=0.
        goto 696
        endif
        ik1=321
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.,0.5,0.,0.5,0.5,0.,0.,0.5,1.)
        the=pmas(lucomp(321),1)+pmas(lucomp(-3122),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
696     if(isinel(295).eq.0)then
        si2=0.
        goto 697
        endif
        ik1=321
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(lucomp(321),1)+pmas(lucomp(-3212),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
697     if(isinel(296).eq.0)then
        si3=0.
        goto 698
        endif
        ik1=311
        ik2=-3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  1.,0.5,0.,0.5,0.5,1.,0.,0.5,1.)
        the=pmas(lucomp(311),1)+pmas(lucomp(-3112),1)
        si3=10*s1724(ss,ilo3,0,the)*fac
698     if(isinel(134).eq.0)then
        si4=0.
        goto 1698
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3334),1)          
        si4=s1724(ss,ilo4,0,the)
1698     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0.and.
     &  ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &  si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=321
        ik2=-3122
        ic=294
c       pion0 + cascade-ba to k+ + lambda0-
        goto 700
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=321
        ik2=-3212
        ic=295
c       pion0 + cascade-ba to k+ + sigma0-
        goto 700
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=311
        ik2=-3112
        ic=296
c       pion0 + cascade-ba to k0 + sigma-ba
        goto 700
        endif
       ik1=-311
        ik2=-3334
        ic=134
c       pion0 + cascade-ba to k0- + omiga-ba
700     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       k+ + cascade0 
        if(kl.eq.321 .and. kl1.eq.3322)then
        if(isinel(298).eq.0)then
        si1=0.
        goto 768
        endif
        ik1=211
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(lucomp(211),1)+pmas(lucomp(3212),1)
        si1=s1724(ss,ilo1,0,the)*fac
768     if(isinel(299).eq.0)then
        si2=0.
        goto 769
        endif
        ik1=211
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(lucomp(211),1)+pmas(lucomp(3122),1)
        si2=s1724(ss,ilo2,0,the)*fac
769     if(isinel(301).eq.0)then
        si3=0.
        goto 770
        endif
        ik1=111
        ik2=3222
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(3222),1)
        si3=s1724(ss,ilo3,0,the)*fac
770     if(isinel(326).eq.0)then
        si4=0.
        goto 771
        endif
        ik1=-311
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &  0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(-311),1)+pmas(lucomp(2212),1)
        si4=s1724(ss,ilo4,0,the)*fac
771     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0
     c   )goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=211
        ik2=3212
        ic=298
c       k+ + cascade0 to pion+ + sigma0
        goto 772
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=211
        ik2=3122
        ic=299
c       k+ + cascade0 to pion+ + lambda
        goto 772
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=111
        ik2=3222
        ic=301
c       k+ + cascade0 to pion0 + sigma+
        goto 772
        endif
        ik1=-311
        ik2=2212
        ic=326
c       k+ + cascade0 to k0- + p
772     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn/10.
        goto 10
        endif

c       k- + cascad0-
        if(kl.eq.-321 .and. kl1.eq.-3322)then
        if(isinel(306).eq.0)then
        si1=0.
        goto 773
        endif
        ik1=-211
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(lucomp(-211),1)+pmas(lucomp(-3212),1)
        si1=s1724(ss,ilo1,0,the)*fac
773     if(isinel(307).eq.0)then
        si2=0.
        goto 774
        endif
        ik1=-211
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(lucomp(-211),1)+pmas(lucomp(3122),1)
        si2=s1724(ss,ilo2,0,the)*fac
774     if(isinel(308).eq.0)then
        si3=0.
        goto 775
        endif
        ik1=111
        ik2=-3222
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(-3222),1)
        si3=s1724(ss,ilo3,0,the)*fac
775     if(isinel(329).eq.0)then
        si4=0.
        goto 776
        endif
        ik1=311
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &  0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(311),1)+pmas(lucomp(-2212),1)
        si4=s1724(ss,ilo4,0,the)*fac
776     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0
     c   )goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=-211
        ik2=-3212
        ic=306
c       k- + cascade0- to pion- + sigma0-
        goto 777
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-211
        ik2=-3122
        ic=307
c       k- + cascade0- to pion- + lambda-
        goto 777
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=111
        ik2=-3222
        ic=308
c       k- + cascade0- to pion0 + sigma+ba
        goto 777
        endif
        ik1=311
        ik2=-2212
        ic=329
c       k- + cascade0- to k0 + p-
777     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cskn/10.
        goto 10
        endif

c       k0 + cascad0
        if(kl.eq.311 .and. kl1.eq.3322)then
        if(isinel(297).eq.0)then
        si1=0.
        goto 778
        endif
        ik1=211
        ik2=3112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(lucomp(211),1)+pmas(lucomp(3112),1)
        si1=s1724(ss,ilo1,0,the)*fac
778     if(isinel(300).eq.0)then
        si2=0.
        goto 779
        endif
        ik1=-211
        ik2=3222
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(lucomp(-211),1)+pmas(lucomp(3222),1)
        si2=s1724(ss,ilo2,0,the)*fac
779     if(isinel(302).eq.0)then
        si3=0.
        goto 780
        endif
        ik1=111
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  0.5,0.5,0.,0.5,1.,1.,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(3212),1)
        si3=s1724(ss,ilo3,0,the)*fac
780     if(isinel(303).eq.0)then
        si4=0.
        goto 781
        endif
        ik1=111
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &  0.5,0.5,0.,0.5,1.,0.,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(3122),1)
        si4=s1724(ss,ilo4,0,the)*fac
781     if(isinel(325).eq.0)then
        si5=0.
        goto 782
        endif
        ik1=-321
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo5,fac,
     &  0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(-321),1)+pmas(lucomp(2212),1)
        si5=s1724(ss,ilo5,0,the)*fac
782     if(isinel(327).eq.0)then
        si6=0.
        goto 783
        endif
        ik1=-311
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo6,fac,
     &  0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(-311),1)+pmas(lucomp(2112),1)
        si6=s1724(ss,ilo6,0,the)*fac
783     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0
     c    .and. ilo5.eq.0.and. ilo6.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6.and.si5.lt.1.e-6.and.si6.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=211
        ik2=3112
        ic=297
c       k0 + cascade0 to pion+ + sigma-
        goto 784
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-211
        ik2=3222
        ic=300
c       k0 + cascade0 to pion- + sigma+
        goto 784
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=111
        ik2=3212
        ic=302
c       k0 + cascade0 to pion0 + sigma0
        goto 784
        endif
        if(rlu1.gt.s3 .and. rlu1.le.s4)then
        ik1=111
        ik2=3122
        ic=303
c       k0 + cascade0 to pion0 + lambda
        goto 784
        endif
        if(rlu1.gt.s4 .and. rlu1.le.s5)then
        ik1=-321
        ik2=2212
        ic=325
c       k0 + cascade0 to k- + p 
        goto 784
        endif
        ik1=-311
        ik2=2112
        ic=327
c       k0 + cascade0 to k0- + n
784     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn/10.
        goto 10
        endif

c       k0- + cascad0-
        if(kl.eq.-311 .and. kl1.eq.-3322)then
        if(isinel(304).eq.0)then
        si1=0.
        goto 785
        endif
        ik1=211
        ik2=-3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,  
     &  0.5,0.5,0.,0.5,1.0,1.0,0.,0.5,1.)
        the=pmas(lucomp(211),1)+pmas(lucomp(-3222),1)
        si1=s1724(ss,ilo1,0,the)*fac
785     if(isinel(305).eq.0)then
        si2=0.
        goto 786
        endif
        ik1=-211
        ik2=-3112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,0.5,0.,0.5,1.0,1.0,0.,0.5,1.)
        the=pmas(lucomp(-211),1)+pmas(lucomp(-3112),1)
        si2=s1724(ss,ilo2,0,the)*fac
786     if(isinel(309).eq.0)then
        si3=0.
        goto 787
        endif
        ik1=111
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  0.5,0.5,0.,0.5,1.0,1.0,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(-3212),1)
        si3=s1724(ss,ilo3,0,the)*fac
787     if(isinel(310).eq.0)then
        si4=0.
        goto 788
        endif
        ik1=111
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo4,fac,
     &  0.5,0.5,0.,0.5,1.0,0.,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(-3122),1)
        si4=s1724(ss,ilo4,0,the)*fac
788     if(isinel(328).eq.0)then
        si5=0.
        goto 789
        endif
        ik1=321
        ik2=-2212
        call srev(kl,kl1,ik1,ik2,ss,ilo5,fac,
     &  0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(321),1)+pmas(lucomp(-2212),1)
        si5=s1724(ss,ilo5,0,the)*fac
789     if(isinel(330).eq.0)then
        si6=0.
        goto 790
        endif
        ik1=311
        ik2=-2112
        call srev(kl,kl1,ik1,ik2,ss,ilo6,fac,
     &  0.5,0.5,0.,0.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(311),1)+pmas(lucomp(-2112),1)
        si6=s1724(ss,ilo6,0,the)*fac
790     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0 .and. ilo4.eq.0
     c    .and. ilo5.eq.0.and. ilo6.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6
     c   .and.si4.lt.1.e-6.and.si5.lt.1.e-6.and.si6.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        si15=si14+si5
        si16=si15+si6
        s1=si1/si16
        s2=si12/si16
        s3=si13/si16
        s4=si14/si16
        s5=si15/si16
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=211
        ik2=-3222
        ic=304
c       k-0 + cascade0- to pion+ + sigma-ba
        goto 791
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-211
        ik2=-3112
        ic=305
c       k0- + cascade0- to pion- + sigma-ba
        goto 791
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=111
        ik2=-3212
        ic=309
c       k0- + cascade0- to pion0 + sigma0-
        goto 791
        endif
        if(rlu1.gt.s3 .and. rlu1.le.s4)then
        ik1=111
        ik2=-3122
        ic=310
c       k0- + cascade0- to pion0 + lambda-
        goto 791
        endif
        if(rlu1.gt.s4 .and. rlu1.le.s5)then
        ik1=321
        ik2=-2212
        ic=328
c       k0- + cascade0- to k+ + p- 
        goto 791
        endif
        ik1=311
        ik2=-2112
        ic=330
c       k0- + cascade0- to k0 + n-
791     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si16/cskn/10.
        goto 10
        endif

c       k+ + omiga- 
        if(kl.eq.321 .and. kl1.eq.3334)then
        if(isinel(331).eq.0)then
        si1=0.
        goto 792
        endif
        ik1=211
        ik2=3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(lucomp(211),1)+pmas(lucomp(3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
792     if(isinel(336).eq.0)then
        si2=0.
        goto 793
        endif
        ik1=111
        ik2=3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
793     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=211
        ik2=3312
        ic=331
c       k+ + omiga- to pion+ + cascade-
        goto 794
        endif
        ik1=111
        ik2=3322
        ic=336
c       k+ + omiga- to pion0 + cascade0
794     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k- + omiga-ba
        if(kl.eq.-321 .and. kl1.eq.-3334)then
        if(isinel(333).eq.0)then
        si1=0.
        goto 795
        endif
        ik1=-211
        ik2=-3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(lucomp(-211),1)+pmas(lucomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
795     if(isinel(338).eq.0)then
        si2=0.
        goto 796
        endif
        ik1=111
        ik2=-3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
796     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=-211
        ik2=-3312
        ic=333
c       k- + omiga-ba to pion- + cascade-ba
        goto 797
        endif
        ik1=111
        ik2=-3322
        ic=338
c       k- + omiga-ba to pion0 + cascade0-
797     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k0 + omiga-
        if(kl.eq.311 .and. kl1.eq.3334)then
        if(isinel(332).eq.0)then
        si1=0.
        goto 798
        endif
        ik1=111
        ik2=3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
798     if(isinel(335).eq.0)then
        si2=0.
        goto 799
        endif
        ik1=-211
        ik2=3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(lucomp(-211),1)+pmas(lucomp(3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
799     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=111
        ik2=3312
        ic=332
c       k0 + omiga- to pion0 + cascade-
        goto 800
        endif
        ik1=-211
        ik2=3322
        ic=335
c       k0 + omiga- to pion- + cascade0
800     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       k0- + omiga-ba
        if(kl.eq.-311 .and. kl1.eq.-3334)then
        if(isinel(334).eq.0)then
        si1=0.
        goto 801
        endif
        ik1=111
        ik2=-3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(lucomp(111),1)+pmas(lucomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
801     if(isinel(337).eq.0)then
        si2=0.
        goto 802
        endif
        ik1=211
        ik2=-3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  0.5,0.,0.,1.5,1.0,0.5,0.,0.5,1.)
        the=pmas(lucomp(211),1)+pmas(lucomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
802     if(ilo1.eq.0.and.ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        sit=si1+si2
        s1=si1/sit
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=111
        ik2=-3312
        ic=334
c       k0- + omiga-ba to pion0 + cascade-ba
        goto 803
        endif
        ik1=211
        ik2=-3322
        ic=337
c       k0- + omiga-ba to pion+ + cascade0-
803     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sit/cskn/10.
        goto 10
        endif

c       pion+ + cascade0 to k0- + sigma+
        if(kl.eq.211 .and. kl1.eq.3322)then
        if(isinel(314).eq.0)goto 13
        the=pmas(lucomp(-311),1)+pmas(lucomp(3222),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo.eq.0) goto 13
        lc(icp,3)=-311
        lc(icp,4)=3222
        lc(icp,5)=314
        ik3=-311
        ik4=3222
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       pion+ + cascade0- 
        if(kl.eq.211 .and. kl1.eq.-3322)then
        if(isinel(319).eq.0)then
        si1=0.
        goto 804
        endif
        ik1=321
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(lucomp(321),1)+pmas(lucomp(-3212),1)
        si1=10.*s1724(ss,ilo1,0,the)*fac
804     if(isinel(320).eq.0)then
        si2=0.
        goto 805
        endif
        ik1=321
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.0,0.5,0.,0.5,0.5,0.0,0.,0.5,1.)
        the=pmas(lucomp(321),1)+pmas(lucomp(-3122),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
805     if(isinel(322).eq.0)then
        si3=0.
        goto 806
        endif
        ik1=311
        ik2=-3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(lucomp(311),1)+pmas(lucomp(-3112),1)
        si3=10*s1724(ss,ilo3,0,the)*fac
806    if(isinel(137).eq.0)then
        si4=0.
        goto 1806
        endif
        the=pmas(lucomp(-311),1)+pmas(lucomp(-3334),1)          
        si4=s1724(ss,ilo4,0,the)
1806    if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0.and.
     &   ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &  si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=321
        ik2=-3212
        ic=319
c       pion+ + cascade0- to k+ + sigma0-
        goto 807
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=321
        ik2=-3122
        ic=320
c       pion+ + cascade0- to k+ + lambda-
        goto 807
        endif
       if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=311
        ik2=-3112
        ic=322
c       pion+ + cascade0- to k0 + sigma-ba
        goto 807
        endif
        ik1=-311
        ik2=-3334
        ic=137
cpion+ + cascade0- to k0- + omiga-ba
807     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion- + cascade0 
        if(kl.eq.-211 .and. kl1.eq.3322)then
        if(isinel(312).eq.0)then
        si1=0.
        goto 808
        endif
        ik1=-321
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(lucomp(-321),1)+pmas(lucomp(3212),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
808     if(isinel(313).eq.0)then
        si2=0.
        goto 809
        endif
        ik1=-321
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.0,0.5,0.,0.5,0.5,0.0,0.,0.5,1.)
        the=pmas(lucomp(-321),1)+pmas(lucomp(3122),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
809     if(isinel(315).eq.0)then
        si3=0.
        goto 810
        endif
        ik1=-311
        ik2=3112
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(lucomp(-311),1)+pmas(lucomp(3112),1)
        si3=10*s1724(ss,ilo3,0,the)*fac
810     if(isinel(135).eq.0)then
        si4=0.
        goto 1810
        endif
        the=pmas(lucomp(311),1)+pmas(lucomp(3334),1)    
        si4=s1724(ss,ilo4,0,the)
                
1810    if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0.and.
     &  ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &  si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=-321
        ik2=3212
        ic=312
c       pion- + cascade0 to k- + sigma0
        goto 811
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-321
        ik2=3122
        ic=313
c       pion- + cascade0 to k- + lambda
        goto 811
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=-311
        ik2=3112
        ic=315
c       pion- + cascade0 to k0- + sigma-
        goto 811
        endif
        ik1=311
        ik2=3334
        ic=135
c       pion- + cascade0 to k0 + omiga-
811     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion- + cascade0- to k0 + sigma+ba
        if(kl.eq.-211 .and. kl1.eq.-3322)then
        if(isinel(321).eq.0)goto 13
        the=pmas(lucomp(311),1)+pmas(lucomp(-3222),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo.eq.0) goto 13
        lc(icp,3)=311
        lc(icp,4)=-3222
        lc(icp,5)=321
        ik3=311
        ik4=-3222
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       pion0 + cascade0 
        if(kl.eq.111 .and. kl1.eq.3322)then
        if(isinel(311).eq.0)then
        si1=0.
        goto 812
        endif
        ik1=-321
        ik2=3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(lucomp(-321),1)+pmas(lucomp(3222),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
812     if(isinel(316).eq.0)then
        si2=0.
        goto 813
        endif
        ik1=-311
        ik2=3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(lucomp(-311),1)+pmas(lucomp(3212),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
813     if(isinel(317).eq.0)then
        si3=0.
        goto 814
        endif
        ik1=-311
        ik2=3122
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  1.0,0.5,0.,0.5,0.5,0.0,0.,0.5,1.)
        the=pmas(lucomp(-311),1)+pmas(lucomp(3122),1)
        si3=10*s1724(ss,ilo3,0,the)*fac
814     if(isinel(136).eq.0)then
        si4=0.
        goto 1814
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(3334),1)            
        si4=s1724(ss,ilo4,0,the)
        
1814      if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0.
     &  and.ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &  si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=-321
        ik2=3222
        ic=311
c       pion0 + cascade0 to k- + sigma+
        goto 815
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-311
        ik2=3212
        ic=316
c       pion0 + cascade0 to k0- + sigma0
        goto 815
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=-311
        ik2=3122
        ic=317
c       pion0 + cascade0 to k0- + lambda
        goto 815
        endif
        ik1=321
        ik2=3334
        ic=136
c       pion0 + cascade0 to k+ + omiga-
815     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion0 + cascade0- 
        if(kl.eq.111 .and. kl1.eq.-3322)then
        if(isinel(318).eq.0)then
        si1=0.
        goto 816
        endif
        ik1=321
        ik2=-3222
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(lucomp(321),1)+pmas(lucomp(-3222),1)
        si1=10*s1724(ss,ilo1,0,the)*fac
816     if(isinel(323).eq.0)then
        si2=0.
        goto 817
        endif
        ik1=311
        ik2=-3212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.0,0.5,0.,0.5,0.5,1.0,0.,0.5,1.)
        the=pmas(lucomp(311),1)+pmas(lucomp(-3212),1)
        si2=10*s1724(ss,ilo2,0,the)*fac
817     if(isinel(324).eq.0)then
        si3=0.
        goto 818
        endif
        ik1=311
        ik2=-3122
        call srev(kl,kl1,ik1,ik2,ss,ilo3,fac,
     &  1.0,0.5,0.,0.5,0.5,0.0,0.,0.5,1.)
        the=pmas(lucomp(311),1)+pmas(lucomp(-3122),1)
        si3=10*s1724(ss,ilo3,0,the)*fac
818     if(isinel(138).eq.0)then
        si4=0.
        goto 1818
        endif
        the=pmas(lucomp(-321),1)+pmas(lucomp(-3334),1)          
        si4=s1724(ss,ilo4,0,the)
        
1818    if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0.and.
     &  ilo4.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6.and.si3.lt.1.e-6.and.
     &  si4.lt.1.e-6)goto 13
        si12=si1+si2
        si13=si12+si3
        si14=si13+si4
        s1=si1/si14
        s2=si12/si14
        s3=si13/si14
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=321
        ik2=-3222
        ic=318
c       pion0 + cascade0- to k+ + sigma+ba
        goto 819
        endif
        if(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=311
        ik2=-3212
        ic=323
c       pion0 + cascade0- to k0 + sigma0-
        goto 819
        endif
        if(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=311
        ik2=-3122
        ic=324
c       pion0 + cascade0- to k0 + lambda-
        goto 819
        endif
        ik1=-321
        ik2=-3334
        ic=138
c       pion0 + cascade0- to k- + omiga-ba
819     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si14/cspin/10.
        goto 10
        endif

c       pion+ + omiga- to k0- + cascade0
        if(kl.eq.211 .and. kl1.eq.3334)then
        if(isinel(342).eq.0)goto 13
        the=pmas(lucomp(-311),1)+pmas(lucomp(3322),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo.eq.0) goto 13
        lc(icp,3)=-311
        lc(icp,4)=3322
        lc(icp,5)=342
        ik3=-311
        ik4=3322
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       pion+ + omiga-ba to k+ + cascade-ba
        if(kl.eq.211 .and. kl1.eq.-3334)then
        if(isinel(343).eq.0)goto 13
        the=pmas(lucomp(321),1)+pmas(lucomp(-3312),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo.eq.0) goto 13
        lc(icp,3)=321
        lc(icp,4)=-3312
        lc(icp,5)=343
        ik3=321
        ik4=-3312
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       pion- + omiga- to k- + cascade-
        if(kl.eq.-211 .and. kl1.eq.3334)then
        if(isinel(339).eq.0)goto 13
        the=pmas(lucomp(-321),1)+pmas(lucomp(3312),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo.eq.0) goto 13
        lc(icp,3)=-321
        lc(icp,4)=3312
        lc(icp,5)=339
        ik3=-321
        ik4=3312
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       pion- + omiga-ba to k0 + cascade0-
        if(kl.eq.-211 .and. kl1.eq.-3334)then
        if(isinel(346).eq.0)goto 13
        the=pmas(lucomp(311),1)+pmas(lucomp(-3322),1)
        ww=s1724(ss,ilo,0,the)/cspin
        if(ilo.eq.0) goto 13
        lc(icp,3)=311
        lc(icp,4)=-3322
        lc(icp,5)=346
        ik3=311
        ik4=-3322
        call srev(kl,kl1,ik3,ik4,ss,ilo1,fac,
     &  1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        tw(icp)=ww*fac
        goto 10
        endif

c       pion0 + omiga- 
        if(kl.eq.111 .and. kl1.eq.3334)then
        if(isinel(340).eq.0)then
        si1=0.
        goto 820
        endif
        ik1=-311
        ik2=3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(-311),1)+pmas(lucomp(3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
820     if(isinel(341).eq.0)then
        si2=0.
        goto 821
        endif
        ik1=-321
        ik2=3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(-321),1)+pmas(lucomp(3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
821     if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=-311
        ik2=3312
        ic=340
c       pion0 + omiga- to k0- + cascade-
        goto 822
        endif
        ik1=-321
        ik2=3322
        ic=341
c       pion0 + omiga- to k- + cascade0
822     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si13/cspin
        goto 10
        endif

c       pion0 + omiga-ba 
        if(kl.eq.111 .and. kl1.eq.-3334)then
        if(isinel(344).eq.0)then
        si1=0.
        goto 823
        endif
        ik1=311
        ik2=-3312
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(311),1)+pmas(lucomp(-3312),1)
        si1=s1724(ss,ilo1,0,the)*fac
823     if(isinel(345).eq.0)then
        si2=0.
        goto 824
        endif
        ik1=321
        ik2=-3322
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.0,0.0,0.,1.5,0.5,0.5,0.,0.5,1.)
        the=pmas(lucomp(321),1)+pmas(lucomp(-3322),1)
        si2=s1724(ss,ilo2,0,the)*fac
824     if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=311
        ik2=-3312
        ic=344
c       pion0 + omiga-ba to k0 + cascade-ba
        goto 825
        endif
        ik1=321
        ik2=-3322
        ic=345
c       pion0 + omiga-ba to k+ + cascade0-
825     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin
        goto 10
        endif
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       pion+ + delta+ to pi+ + p
        if((kl.eq.211 .and. kl1.eq.2214).or.
     &  (kl1.eq.211 .and. kl.eq.2214))then
        if(isinel(351).eq.0)goto 13
        ik3=211
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &  1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=2212
        lc(icp,5)=351
        ww=sdelta(ss,ilo1,1,0.)/cspin/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       pion0 + delta++ to pi+ + p
        if((kl.eq.111 .and. kl1.eq.2224).or.
     &  (kl1.eq.111 .and. kl.eq.2224))then
        if(isinel(352).eq.0)goto 13
        ik3=211
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &  1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=2212
        lc(icp,5)=352
        ww=sdelta(ss,ilo1,1,0.)/cspin/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       pion+ + delta- 
        if((kl.eq.211 .and. kl1.eq.1114).or.
     &  (kl1.eq.211 .and. kl.eq.1114))then
        if(isinel(347).eq.0)then
        si1=0.
        goto 1000
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        
        if(ilo1.eq.0)then
        si1=0.
        else
        si1=sdelta(ss,ilo1,1,0.)*WEIGH(19)*fac
        endif
c       cross section of pion+ + delta- to pi- + p

1000    if(isinel(348).eq.0)then
        si2=0.
        goto 1002
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=sdelta(ss,ilo2,1,0.)*WEIGH(19)*fac
        endif
c       cross section of pion+ + delta- to pi0 + n

1002    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=-211
        ik2=2212
        ic=347
c       cross section of pion+ + delta- to pi- + p
        goto 1004
        endif
        ik1=111
        ik2=2112
        ic=348
c       cross section of pion+ + delta- to pi0 + n
1004    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       pion+ + delta0
        if((kl.eq.211 .and. kl1.eq.2114).or.
     &  (kl1.eq.211 .and. kl.eq.2114))then
        if(isinel(349).eq.0)then
        si1=0.
        goto 1006
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        
        if(ilo1.eq.0)then
        si1=0.
        else
        si1=sdelta(ss,ilo1,1,0.)*WEIGH(19)*fac
        endif
c       cross section of pion+ + delta0 to pi+ + n

1006    if(isinel(350).eq.0)then
        si2=0.
        goto 1008
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=sdelta(ss,ilo2,1,0.)*WEIGH(19)*fac
        endif
c       cross section of pion+ + delta0 to pi0 + p

1008    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=211
        ik2=2112
        ic=349
c       cross section of pion+ + delta0 to pi+ + n
        goto 1010
        endif
        ik1=111
        ik2=2212
        ic=350
c       cross section of pion+ + delta0 to pi0 + p
1010    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif   

c       pion0 + delta+
        if((kl.eq.111 .and. kl1.eq.2214).or.
     &  (kl1.eq.111 .and. kl.eq.2214))then
        if(isinel(353).eq.0)then
        si1=0.
        goto 1012
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        
        if(ilo1.eq.0)then
        si1=0.
        else
        si1=sdelta(ss,ilo1,1,0.)*WEIGH(19)*fac
        endif
c       cross section of pion0 + delta+ to pi0 + p

1012    if(isinel(354).eq.0)then
        si2=0.
        goto 1014
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=sdelta(ss,ilo2,1,0.)*WEIGH(19)*fac
        endif
c       cross section of pion0 + delta+ to pi+ + n

1014    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=2212
        ic=353
c       cross section of pion0 + delta+ to pi0 + p
        goto 1016
        endif
        ik1=211
        ik2=2112
        ic=354
c       cross section of pion0 + delta+ to pi+ + n
1016    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif   
c       pion0 + delta0
        if((kl.eq.111 .and. kl1.eq.2114).or.
     &  (kl1.eq.111 .and. kl.eq.2114))then
        if(isinel(355).eq.0)then
        si1=0.
        goto 1018
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        
        if(ilo1.eq.0)then
        si1=0.
        else
        si1=sdelta(ss,ilo1,1,0.)*WEIGH(19)*fac
        endif
c       cross section of pion0 + delta0 to pi0 + n

1018    if(isinel(356).eq.0)then
        si2=0.
        goto 1020
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=sdelta(ss,ilo2,1,0.)*WEIGH(19)*fac
        endif
c       cross section of pion0 + delta0 to pi- + p

1020    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=2112
        ic=355
c       cross section of pion0 + delta0 to pi0 + n
        goto 1022
        endif
        ik1=-211
        ik2=2212
        ic=356
c       cross section of pion0 + delta0 to pi- + p
1022    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif   
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       pion0 + delta- to pi- + n
        if((kl.eq.111 .and. kl1.eq.1114).or.
     &  (kl1.eq.111 .and. kl.eq.1114))then
        if(isinel(357).eq.0)goto 13
        ik3=-211
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &  1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=2112
        lc(icp,5)=357
        ww=sdelta(ss,ilo1,1,0.)/cspin/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       pion- + delta++
        if((kl.eq.-211 .and. kl1.eq.2224).or.
     &  (kl1.eq.-211 .and. kl.eq.2224))then
        if(isinel(358).eq.0)then
        si1=0.
        goto 1024
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        
        if(ilo1.eq.0)then
        si1=0.
        else
        si1=sdelta(ss,ilo1,1,0.)*WEIGH(19)*fac
        endif
c       cross section of pion- + delta++ to pi0 + p

1024    if(isinel(359).eq.0)then
        si2=0.
        goto 1026
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=sdelta(ss,ilo2,1,0.)*WEIGH(19)*fac
        endif
c       cross section of pion- + delta++ to pi+ + n

1026    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=2212
        ic=358
c       cross section of pion- + delta++ to pi0 + p
        goto 1028
        endif
        ik1=211
        ik2=2112
        ic=359
c       cross section of pion- + delta++ to pi+ + n
1028    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif   
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       pion- + delta+
        if((kl.eq.-211 .and. kl1.eq.2214).or.
     &  (kl1.eq.-211 .and. kl.eq.2214))then
        if(isinel(360).eq.0)then
        si1=0.
        goto 1030
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        
        if(ilo1.eq.0)then
        si1=0.
        else
        si1=sdelta(ss,ilo1,1,0.)*WEIGH(19)*fac
        endif
c       cross section of pion- + delta+ to pi- + p

1030    if(isinel(361).eq.0)then
        si2=0.
        goto 1032
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=sdelta(ss,ilo2,1,0.)*WEIGH(19)*fac
        endif
c       cross section of pion- + delta+ to pi0 + n

1032    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=-211
        ik2=2212
        ic=360
c       cross section of pion- + delta+ to pi- + p
        goto 1034
        endif
        ik1=111
        ik2=2112
        ic=361
c       cross section of pion- + delta+ to pi0 + n
1034    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif   
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       pion- + delta0 to pi- + n
        if((kl.eq.-211 .and. kl1.eq.2114).or.
     &  (kl1.eq.-211 .and. kl.eq.2114))then
        if(isinel(362).eq.0)goto 13
        ik3=-211
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &  1.0,1.5,0.,1.5,1.0,0.5,0.,0.5,1.)
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=2112
        lc(icp,5)=362
        ww=sdelta(ss,ilo1,1,0.)/cspin/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       prho0 + n
        if((kl.eq.113 .and. kl1.eq.2112).or.
     &  (kl1.eq.113 .and. kl.eq.2112))then
        if(isinel(363).eq.0)then
        si1=0.
        goto 1036
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        
        if(ilo1.eq.0)then
        si1=0.
        else
        si1=srho(ss,ilo1,1,0.)*WEIGH(19)*fac
        endif
c       cross section of rho0 + n to pi- + p

1036    if(isinel(364).eq.0)then
        si2=0.
        goto 1038
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=srho(ss,ilo2,1,0.)*WEIGH(19)*fac
        endif
c       cross section of rho0 + n to pi0 + n

1038    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=-211
        ik2=2212
        ic=363
c       cross section of rho0 + n to pi- + p
        goto 1040
        endif
        ik1=111
        ik2=2112
        ic=364
c       cross section of rho0 + n to pi0 + n
1040    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif   
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       prho+ + n
        if((kl.eq.213 .and. kl1.eq.2112).or.
     &  (kl1.eq.213 .and. kl.eq.2112))then
        if(isinel(366).eq.0)then
        si1=0.
        goto 1042
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        
        if(ilo1.eq.0)then
        si1=0.
        else
        si1=srho(ss,ilo1,1,0.)*WEIGH(19)*fac
        endif
c       cross section of rho+ + n to pi0 + p

1042    if(isinel(367).eq.0)then
        si2=0.
        goto 1044
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=srho(ss,ilo2,1,0.)*WEIGH(19)*fac
        endif
c       cross section of rho+ + n to pi+ + n

1044    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=2212
        ic=366
c       cross section of rho+ + n to pi0 + p
        goto 1046
        endif
        ik1=211
        ik2=2112
        ic=367
c       cross section of rho+ + n to pi+ + n
1046    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif   
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       rho0 + p
        if((kl.eq.113 .and. kl1.eq.2212).or.
     &  (kl1.eq.113 .and. kl.eq.2212))then
        if(isinel(368).eq.0)then
        si1=0.
        goto 1048
        endif
        ik1=111
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        
        if(ilo1.eq.0)then
        si1=0.
        else
        si1=srho(ss,ilo1,1,0.)*WEIGH(19)*fac
        endif
c       cross section of rho0 + p to pi0 + p

1048    if(isinel(369).eq.0)then
        si2=0.
        goto 1050
        endif
        ik1=211
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=srho(ss,ilo2,1,0.)*WEIGH(19)*fac
        endif
c       cross section of rho0 + p to pi+ + n

1050    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=2212
        ic=368
c       cross section of rho0 + p to pi0 + p
        goto 1052
        endif
        ik1=211
        ik2=2112
        ic=369
c       cross section of rho0 + p to pi+ + n
1052    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif   
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       rho- + p
        if((kl.eq.-213 .and. kl1.eq.2212).or.
     &  (kl1.eq.-213 .and. kl.eq.2212))then
        if(isinel(370).eq.0)then
        si1=0.
        goto 1054
        endif
        ik1=111
        ik2=2112
        call srev(kl,kl1,ik1,ik2,ss,ilo1,fac,
     &  1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        
        if(ilo1.eq.0)then
        si1=0.
        else
        si1=srho(ss,ilo1,1,0.)*WEIGH(19)*fac
        endif
c       cross section of rho- + p to pi0 + n

1054    if(isinel(371).eq.0)then
        si2=0.
        goto 1056
        endif
        ik1=-211
        ik2=2212
        call srev(kl,kl1,ik1,ik2,ss,ilo2,fac,
     &  1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo2.eq.0)then
        si2=0.
        else
c since threshold energy is check in srev() we set the=0. here in
c sdelta()
        si2=srho(ss,ilo2,1,0.)*WEIGH(19)*fac
        endif
c       cross section of rho- + p to pi- + p

1056    if(ilo1.eq.0 .and. ilo2.eq.0)goto 13
        if(si1.lt.1.e-6.and.si2.lt.1.e-6)goto 13
        si12=si1+si2
        s1=si1/si12
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=2112
        ic=370
c       cross section of rho- + p to pi0 + n
        goto 1058
        endif
        ik1=-211
        ik2=2212
        ic=371
c       cross section of rho- + p to pi- + p
1058    lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=si12/cspin/10.
        goto 10
        endif   
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       rho+ + p to pi+ + p
        if((kl.eq.213 .and. kl1.eq.2212).or.
     &  (kl1.eq.213 .and. kl.eq.2212))then
        if(isinel(372).eq.0)goto 13
        ik3=211
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &  1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo.eq.0) goto 13
        lc(icp,3)=211
        lc(icp,4)=2212
        lc(icp,5)=372
        ww=srho(ss,ilo1,1,0.)/cspin/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       rho- + n to pi- + n
        if((kl.eq.-213 .and. kl1.eq.2112).or.
     &  (kl1.eq.-213 .and. kl.eq.2112))then
        if(isinel(365).eq.0)goto 13
        ik3=-211
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &  1.0,0.5,1.0,0.5,1.0,0.5,0.,0.5,1.)
        if(ilo.eq.0) goto 13
        lc(icp,3)=-211
        lc(icp,4)=2112
        lc(icp,5)=365
        ww=srho(ss,ilo1,1,0.)/cspin/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       delta++ + n to p + p
        if((kl.eq.2224 .and. kl1.eq.2112).or.
     &  (kl1.eq.2224 .and. kl.eq.2112))then
        if(isinel(373).eq.0)goto 13
        ik3=2212
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &  1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,0.5)
        if(ilo.eq.0) goto 13
        lc(icp,3)=2212
        lc(icp,4)=2212
        lc(icp,5)=373
        ww=snn(ss,ilo1,1,0.)/csnn/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       delta+ + n to p + n
        if((kl.eq.2214 .and. kl1.eq.2112).or.
     &  (kl1.eq.2214 .and. kl.eq.2112))then
        if(isinel(374).eq.0)goto 13
        ik3=2212
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &  1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,1.0)
        if(ilo.eq.0) goto 13
        lc(icp,3)=2212
        lc(icp,4)=2112
        lc(icp,5)=374
        ww=snn(ss,ilo1,1,0.)/csnn/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       delta+ + p to p + p
        if((kl.eq.2214 .and. kl1.eq.2212).or.
     &  (kl1.eq.2214 .and. kl.eq.2212))then
        if(isinel(375).eq.0)goto 13
        ik3=2212
        ik4=2212
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &  1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,0.5)
        if(ilo.eq.0) goto 13
        lc(icp,3)=2212
        lc(icp,4)=2212
        lc(icp,5)=375
        ww=snn(ss,ilo1,1,0.)/csnn/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       delta0 + p to p + n
        if((kl.eq.2114 .and. kl1.eq.2212).or.
     &  (kl1.eq.2114 .and. kl.eq.2212))then
        if(isinel(376).eq.0)goto 13
        ik3=2212
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &  1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,1.0)
        if(ilo.eq.0) goto 13
        lc(icp,3)=2212
        lc(icp,4)=2112
        lc(icp,5)=376
        ww=snn(ss,ilo1,1,0.)/csnn/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       delta0 + n to n + n
        if((kl.eq.2114 .and. kl1.eq.2112).or.
     &  (kl1.eq.2114 .and. kl.eq.2112))then
        if(isinel(377).eq.0)goto 13
        ik3=2112
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &  1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,0.5)
        if(ilo.eq.0) goto 13
        lc(icp,3)=2112
        lc(icp,4)=2112
        lc(icp,5)=377
        ww=snn(ss,ilo1,1,0.)/csnn/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c       delta- + p to n + n
        if((kl.eq.1114 .and. kl1.eq.2212).or.
     &  (kl1.eq.1114 .and. kl.eq.2212))then
        if(isinel(378).eq.0)goto 13
        ik3=2112
        ik4=2112
        call srev(kl,kl1,ik3,ik4,ss,ilo,fac,
     &  1.5,0.5,1.5,0.5,0.5,0.5,0.5,0.5,0.5)
        if(ilo.eq.0) goto 13
        lc(icp,3)=2112
        lc(icp,4)=2112
        lc(icp,5)=378
        ww=snn(ss,ilo1,1,0.)/csnn/10.
        tw(icp)=ww*fac
        goto 10
        endif
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


c       lambda- p annihilation
c       lambdabar + p -----> K*+ + omiga
        if((kl.eq.-3122.and.kl1.eq.2212)
     c   .or.(kl1.eq.-3122.and.kl.eq.2212))then
        if(isinel(393).eq.0)goto 13
        lc(icp,3)=323
        lc(icp,4)=223
        lc(icp,5)=393
        tw(icp)=PARAM14/csnn
        goto 10
        endif

c       lambda- n annihilation
c       lambdabar + n -----> K*0 + omiga
        if((kl.eq.-3122.and.kl1.eq.2112)
     c   .or.(kl1.eq.-3122.and.kl.eq.2112))then
        if(isinel(394).eq.0)goto 13
        tw(icp)=PARAM14/csnn
c       if(ilo.eq.0) goto 13
        lc(icp,3)=313
        lc(icp,4)=223
        lc(icp,5)=394
        goto 10
        endif

c       sigma0- p annihilation
c        sigma0bar + p -----> K*+ + omiga
        if((kl.eq.-3212.and.kl1.eq.2212)
     c   .or.(kl1.eq.-3212.and.kl.eq.2212))then
        if(isinel(395).eq.0)goto 13
        tw(icp)=PARAM14/csnn
c       if(ilo.eq.0) goto 13
        lc(icp,3)=323
        lc(icp,4)=223
        lc(icp,5)=395
        goto 10
        endif

c       sigma0- n annihilation
c       sigma0bar + n -----> K*0 + omiga
        if((kl.eq.-3212.and.kl1.eq.2112)
     c   .or.(kl1.eq.-3212.and.kl.eq.2112))then
        if(isinel(396).eq.0)goto 13
        tw(icp)=PARAM14/csnn
c       if(ilo.eq.0) goto 13
        lc(icp,3)=313
        lc(icp,4)=223
        lc(icp,5)=396
        goto 10
        endif

c       p- p annihilation
c       pbar p annihilation -----> rho0+omiga
        if((kl.eq.-2212.and.kl1.eq.2212)
     c   .or.(kl1.eq.-2212.and.kl.eq.2212))then
        if(isinel(397).eq.0)goto 13
        tw(icp)=para13/csnn
c       if(ilo.eq.0) goto 13
        lc(icp,3)=113
        lc(icp,4)=223
        lc(icp,5)=397
        goto 10
        endif

c       p- n annihilation
c       pbar n annihilation -----> rho- + omiga
        if((kl.eq.-2212.and.kl1.eq.2112)
     c   .or.(kl1.eq.-2212.and.kl.eq.2112))then
        if(isinel(398).eq.0)goto 13
        tw(icp)=para13/csnn
c       if(ilo.eq.0) goto 13
        lc(icp,3)=-213
        lc(icp,4)=223
        lc(icp,5)=398
        goto 10
        endif

c       n- p annihilation
c       nbar p annihilation -----> rho+ + omiga
        if((kl.eq.-2112.and.kl1.eq.2212)
     c   .or.(kl1.eq.-2112.and.kl.eq.2212))then
        if(isinel(399).eq.0)goto 13
        tw(icp)=para13/csnn
c       if(ilo.eq.0) goto 13
        lc(icp,3)=213
        lc(icp,4)=223
        lc(icp,5)=399
        goto 10
        endif

c       n- n annihilation
        if((kl.eq.-2112.and.kl1.eq.2112)
     c   .or.(kl1.eq.-2112.and.kl.eq.2112))then
c       nbar n annihilation -----> rho0 + omiga
        if(isinel(400).eq.0)goto 13
        tw(icp)=para13/csnn
c       if(ilo.eq.0) goto 13
        lc(icp,3)=113
        lc(icp,4)=223
        lc(icp,5)=400
        goto 10
        endif

        write(mstu(11),*)'nothing is found,kl,kl1=',kl,kl1
13      do m=3,5
        lc(icp,m)=0
        enddo   
        tw(icp)=0.
        goto 999
c10     if(lc(icp,5).ne.397.and.lc(icp,5).ne.398.and.lc(icp,5).ne.399
c     & .and.lc(icp,5).ne.400)then
10      nchargef=LUCHGE(lc(icp,3))+LUCHGE(lc(icp,4))
        if(nchargei.ne.nchargef/3)write(mstu(11),*)'initial,final,
     &  lc(icp,5)=',nchargei,nchargef/3,lc(icp,5)
999     return
        end     

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine spipi(lc3,lc4,ss,ilo)
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
        am1=pmas(lucomp(lc3),1)
        am2=pmas(lucomp(lc4),1)
        ilo=1
        if(ss.lt.(am1+am2))ilo=0
        return
        end



C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine srev(kl,kl1,lc3,lc4,ss,ilo,fac,xii1,xii2,xsi1,xsi2,
     &  xif1,xif2,xsf1,xsf2,pauli)
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/LUCIDAT2/KFMAXT,PARAM(20),WEIGH(400)       
        am1=pmas(lucomp(kl),1)
        am2=pmas(lucomp(kl1),1)
        am3=pmas(lucomp(lc3),1)
        am4=pmas(lucomp(lc4),1)
        am12=am1*am1
        am22=am2*am2
        am32=am3*am3
        am42=am4*am4
        ss2=ss*ss
        ss4=ss2*ss2
        ilo=1
        
        if(ss.lt.(am3+am4))then
        ilo=0
        fac=0.
        else
        pfp=ss4-2.*ss2*(am32+am42)+(am32-am42)*(am32-am42)
        if(pfp.lt.0.)pfp=1.e-10
        pfp=sqrt(pfp)/2./ss
        pip=ss4-2.*ss2*(am12+am22)+(am12-am22)*(am12-am22)
        if(pip.lt.0.)pip=1.e-10
        pip=sqrt(pip)/2./ss
c       phase=pauli*(2*xif1+1)*(2*xif2+1)*(2*xsf1+1)*(2*xsf2+1)/
c     & ((2*xii1+1)*(2*xii2+1)*(2*xsi1+1)*(2*xsi2+1))
        phase=pauli*(2*xsf1+1)*(2*xsf2+1)/((2*xsi1+1)*(2*xsi2+1))
        fac=phase*pfp/pip
        endif
        return
        end





C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine updpli(l,l1,icp,ss,pii,pjj,lc,tc,tw,winel,time,iia)
c       it uses to update particle list after inelastic collision &
c        truncates collision list correspondingly.
        parameter (mcludi=40000)
        parameter(KSZJ=40000)
        parameter(nsize=100000)
        common/wz/c17(mcludi,3),ishp(mcludi),tp(mcludi)
     c   ,coor(3)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/ctllist/nctl,noinel(400),nctl0
        common/sa4/tau(kszj),tlco(kszj,4)
        common/sa5/kfmax,kfaco(50),numb(50),disbe(50,50)
        common/sa6/kfmaxi,nwhole,nna
        common/sa12/psa(5),ptai(5),clorenp,clorent
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn
     c  ,rnt,rnp,rao,rou0
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        dimension pii(4),pjj(4),pp(4)
        integer winel
        kf1=k(l,2)
        kf2=k(l1,2)
        n1=nwhole-n
c       if(iia.eq.397 .or. iia.eq.398 .or. iia.eq.399.or.iia.eq.400)goto 
c     c  1300   
1400    ik1=lc(icp,3)
        ik2=lc(icp,4)
c       put the scattered (produced) particles into particle list
c        (i.e. update particle list in inelastic scattering case)
c        & turncate collision list correspondingly.
        ll=l
        ll1=l1
        kf=ik1
        do i=1,4
        pp(i)=pii(i)
        enddo
        do 500 i=1,2
        nnn=0
        do 600 j=1,kfmax
        if(kf.ne.kfaco(j))goto 900
        jj=numb(j)+1
c       update the particle list.

1000    do m=n+n1,jj,-1
        mm=m+1
        k(mm,2)=k(m,2)
        k(mm,1)=1
        k(mm,3)=k(m,3)
        do m1=1,3
        p(mm,m1)=p(m,m1)
        c17(mm,m1)=c17(m,m1)
        enddo
        do m2=4,5
        p(mm,m2)=p(m,m2)
        enddo
        ishp(mm)=ishp(m)
        tp(mm)=tp(m)
        tau(mm)=tau(m)
        do itlc=1,4
        tlco(mm,itlc)=tlco(m,itlc)
        enddo
        enddo
        if(ll.ge.jj)ll=ll+1
        if(ll1.ge.jj)ll1=ll1+1
c       update the values of lc(m,1-2) with m.eq.ll & .gt.jj.
        do m=1,nctl
        lc1=lc(m,1)
        if(lc1.ge.jj)lc(m,1)=lc1+1
        lc2=lc(m,2)
        if(lc2.ge.jj)lc(m,2)=lc2+1
        enddo
        do m=1,nctl
        lc1=lc(m,1)
        lc2=lc(m,2)
        if(lc1.eq.ll)lc(m,1)=jj
        if(lc2.eq.ll)lc(m,2)=jj
        enddo
c       give proper values to particle jj.
        k(jj,2)=kf
        k(jj,1)=1
        k(jj,3)=0
        tlco(jj,4)=tlco(ll,4)
        do m=1,3
        p(jj,m)=pp(m)
        c17(jj,m)=c17(ll,m)
        tlco(jj,m)=tlco(ll,m)
        enddo
        p(jj,4)=pp(4)
        p(jj,5)=pmas(lucomp(kf),1)
        ishp(jj)=ishp(ll)
        tp(jj)=tp(ll)
Cwe  give a zero formation time for particles produced from
c rescatttering 
c       tau(jj)=t0*p(jj,4)/p(jj,5)
        tau(jj)=time
        lc(icp,i)=jj
        if(nnn.eq.1)goto 300
        do m=j,kfmax
        numb(m)=numb(m)+1
        enddo
300     goto 200
900     if(j.lt.kfmax)goto 600
c       scattered particle with new flavor.
        nnn=1
        jj=numb(kfmax)+1
        kfmax=kfmax+1
        numb(kfmax)=jj
        kfaco(kfmax)=kf
        goto 1000
600     continue
200     continue
        n=n+1
        if(i.eq.2)goto 500
        ll2=ll
        ll=ll1  
        ll1=ll2
        kf=ik2
        do j=1,4
        pp(j)=pjj(j)
        enddo
500     continue
        l=ll1
        l1=ll
1100    continue
c       take out the scattering particles from particle list (i.e.
c        update particle list in inelastic scattering case) &
c        truncate the collision list correspondingly.
1200    kf=kf1
        ll=l
        do 700 i=1,2
        if(ll.eq.n)then
        kfd=numb(kfmax)-numb(kfmax-1)
c       following statements are added by sa on 07/March/97
        if(kfd.eq.0)then   
        numbm=numb(kfmax)
        do i1=1,kfmax
        if(numb(i1).ne.numbm)goto 3000
        i2=i1   
        goto 3001
3000    enddo
3001    do i1=i2,kfmax
        numb(i1)=numb(i1)-1
        enddo
        goto 400
        endif

        if(kfd.eq.1)then
        kfmax=kfmax-1
        goto 400
        endif
        numb(kfmax)=numb(kfmax)-1
400     goto 100
        endif
        do j=ll+1,n+n1
        jj=j-1
        k(jj,2)=k(j,2)
        k(jj,1)=1
        k(jj,3)=k(j,3)
        tlco(jj,4)=tlco(J,4)
        do m=1,3
        p(jj,m)=p(j,m)
        c17(jj,m)=c17(j,m)
        tlco(jj,m)=tlco(J,m)
        enddo
        do m=4,5
        p(jj,m)=p(j,m)
        enddo
        ishp(jj)=ishp(j)
        tp(jj)=tp(j)
        tau(jj)=tau(j)
        enddo
        do m=1,nctl
        lc1=lc(m,1)
        lc2=lc(m,2)
        if(lc1.gt.ll)lc(m,1)=lc1-1
        if(lc2.gt.ll)lc(m,2)=lc2-1
        enddo
        do 800 j=1,kfmax
        if(kf.ne.kfaco(j))goto 800
        do m=j,kfmax
        numb(m)=numb(m)-1
        enddo
        goto 100
800     continue
100     continue
        n=n-1
        if(l1.gt.ll)l1=l1-1
        if(i.eq.2)goto 700
        ll=l1
        kf=kf2
700     continue
c       if(iia.eq.397 .or. iia.eq.398 .or. iia.eq.399.or.iia.eq.400)goto 
c     c  1500
1600    return
        end


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&




        subroutine coelas(ic,jc,eij,pi,pj)
c       perform elastic scattering
        parameter (mcludi=40000)
        parameter (KSZJ=40000)
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
        COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
        dimension pi(4),pj(4)
        double precision abt
        iic=k(ic,2)
        jjc=k(jc,2)
        d=3.65*(eij-pmas(lucomp(iic),1)-pmas(lucomp(jjc),1))
        if(d.lt.1.e-10)return
        pt=0.2
        a=min(10.3,1./(1.12*pt)/(1.12*pt))
        d6=d**6
        b=d6*a/(1.+d6)
        if(b.lt.1.e-20)then
        b=1.e-20
        endif
        pm2=pi(1)**2+pi(2)**2+pi(3)**2
        pm=sqrt(pm2)
        t0=-4.*pm2
        if(abs(t0).lt.1.e-20)then
        cctas=1.
        goto 100
        endif
        cc=rlu(1)
        if(abs(b*t0).lt.0.0001)then
        abt=1.
c       elseif(b*t0.lt.-50.)then
c       abt=0.
        else
        abt=dexp(dmax1(-7.0D2,dble(b*t0)))
        endif
        tt1=dlog(cc+(1.-cc)*abt)
        if(abs(tt1).lt.1.e-30 .and. b.le.1.e-20)then
        cctas=1.
        goto 100
        endif
        tt=tt1/b
        if(abs(tt).lt.1.e-20)then
        cctas=1.
        goto 100
        endif
        cctas=1.-tt*2./t0
        if(abs(cctas).gt.1.)then
        cctas=sign(1.,cctas)
        endif
100     continue
        sctas=sqrt(1.-cctas**2)
        fis=2.*3.1416*rlu(1)
        cfis=cos(fis)
        sfis=sin(fis)
        call rotate(cctas,sctas,cfis,sfis,pm,pi,pj)
        if(pi(4).lt.0..or.pj(4).lt.0.)then
        write(mstu(11),*)'An error might occur in subroutine
     &  coelas()','pi,pj=',pi,pj
        endif
        return
        end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&



        subroutine rotate(cctas,sctas,cfis,sfis,pp3,pi,pj)
c       perform the rotate for elastic scattering
        dimension pi(4),pj(4)
        fi1=atan2(pi(2),pi(1))
        cta1=atan2(sqrt(pi(1)**2+pi(2)**2),pi(3))
        cfi1=cos(fi1)
        sfi1=sin(fi1)
        ccta1=cos(cta1)
        scta1=sin(cta1)
        pi(1)=cfi1*(ccta1*sctas*cfis+scta1*cctas)-sfi1*sctas*sfis
        pi(2)=sfi1*(ccta1*sctas*cfis+scta1*cctas)+cfi1*sctas*sfis
        pi(3)=ccta1*cctas-scta1*sctas*cfis
        pi(1)=pp3*pi(1)
        pi(2)=pp3*pi(2)
        pi(3)=pp3*pi(3)
        do i=1,3
        pj(i)=0.-pi(i)
        enddo
        return
        end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


        subroutine lorntz(ilo,b,pi,pj)
c       It uses to perform Lorentz (or inverse Lorentz) transformation
c       implicit real*8 (a-h,o-z)
        double precision b(3),bb,bbb,gam,ga,pib,pjb
        dimension pi(4),pj(4)
        bb=b(1)*b(1)+b(2)*b(2)+b(3)*b(3)
        bbb=1.-bb
        if(bbb.le.1.d-10)bbb=1.d-10
        gam=1./dsqrt(bbb)
        ga=gam*gam/(gam+1.)
        if(ilo.eq.1) goto 100
c       Lorentz transformation
        pib=pi(1)*b(1)+pi(2)*b(2)+pi(3)*b(3)
        pjb=pj(1)*b(1)+pj(2)*b(2)+pj(3)*b(3)
        do i=1,3
        pi(i)=pi(i)+b(i)*(ga*pib-gam*pi(4))
        pj(i)=pj(i)+b(i)*(ga*pjb-gam*pj(4))
        enddo
        pi(4)=gam*(pi(4)-pib)
        pj(4)=gam*(pj(4)-pjb)
        return
100     continue
c       inverse Lorentz transformation
        pib=pi(1)*b(1)+pi(2)*b(2)+pi(3)*b(3)
        pjb=pj(1)*b(1)+pj(2)*b(2)+pj(3)*b(3)
        do i=1,3
        pi(i)=pi(i)+b(i)*(ga*pib+gam*pi(4))
        pj(i)=pj(i)+b(i)*(ga*pjb+gam*pj(4))
        enddo
        pi(4)=gam*(pi(4)+pib)
        pj(4)=gam*(pj(4)+pjb)
        return
        end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        subroutine updple(ic,jc,b,pi,pj,time)
c       update the particle list for elastic scattering only.
        parameter (mcludi=40000)
        parameter (KSZJ=40000)
        double precision b(3)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
      COMMON/LUDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200) 
        common/sa4/tau(kszj),tlco(kszj,4)
        common/wz/c17(mcludi,3),ishp(mcludi),tp(mcludi)
     c   ,coor(3)
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn
     c  ,rnt,rnp,rao,rou0
        dimension pi(4),pj(4)
        save /LUJETS/
c       k(ic,3)=0
c       k(jc,3)=0
        ilo=1
c       ilo=1 for inverse Lorentz transformation
        call lorntz(ilo,b,pi,pj)
        if(pi(4).lt.0..or.pj(4).lt.0.)then
        write(mstu(11),*)'An error might occur in subroutine
     &  updple()','pi,pj=',pi,pj
        endif
        p(ic,4)=pi(4)
        p(jc,4)=pj(4)
        do i=1,3
        p(ic,i)=pi(i)
        p(jc,i)=pj(i)
        enddo   
        tau(ic)=time
        tau(jc)=time
c       tau(ic)=time+0.5*t0*p(ic,4)/p(ic,5)
c       tau(jc)=time+0.5*t0*p(jc,4)/p(jc,5)
c the spectators which have experienced elastic reactions are given
c k(ic,4)=333 for projectile spectators and k(ic,4)=555 for target 
c       spectators.
        if((k(ic,2).eq.2112.or.k(ic,2).eq.2212).and.(k(ic,3).eq.33))
     &  k(ic,4)=333
        if((k(ic,2).eq.2112.or.k(ic,2).eq.2212).and.(k(ic,3).eq.55))
     &  k(ic,4)=555
        if((k(jc,2).eq.2112.or.k(jc,2).eq.2212).and.(k(jc,3).eq.33))
     &  k(jc,4)=333
        if((k(jc,2).eq.2112.or.k(jc,2).eq.2212).and.(k(jc,3).eq.55))
     &  k(jc,4)=555
        return
        end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


        subroutine updatl(ic,jc,time,lc,tc,tw,winel,n)
c       update the collision time list after either an elastic or an 
c       inelastic interaction
        parameter (KSZ1=30)
        parameter(nsize=100000)
      COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
        common/sa5/kfmax,kfaco(50),numb(50),disbe(50,50)
        common/ctllist/nctl,noinel(400),nctl0
        common/papr/t0,sig,dep,ddt,edipi,epin,ecsnn,ekn
     c  ,rnt,rnp,rao,rou0
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        integer winel

        m32=numb(kfmax)
c       loop over old colliding pairs
        j=0     
        do i=1,nctl
        i1=lc(i,1)
        j1=lc(i,2)
c       ia=(i1-ic)*(j1-jc)*(i1-jc)*(j1-ic)
c       if(ia.eq.0) goto 400
        if(i1.eq.ic .or. i1.eq.jc)goto 400
        if(j1.eq.ic .or. j1.eq.jc)goto 400
        if((tc(i)-time).le.ddt) goto 400
c       through away the pairs which have tc<= time
        j=j+1
        tc(j)=tc(i)
        tw(j)=tw(i)
        do m=1,5
        lc(j,m)=lc(i,m)
        enddo
400     continue
        enddo
        do i=j+1,nctl+1
        tc(i)=0.0
        tw(i)=0.0
        do m=1,5
        lc(i,m)=0
        enddo
        enddo
c       loop over particle list
        nctl=j+1
        if(kfr(22).eq.1)then
        j1=ic
        do ik=1,2
        if(j1.gt.m32)goto 300
        do i=1,m32
        if(nctl.gt.nsize)stop 100000
        i1=i
        iflag=0
        if((i.eq.jc.and.j1.eq.ic).or.(i.eq.ic.and.j1.eq.jc))
     &  goto 100
c do not pair these two particles which was just yielded into a
c new collision pair.
        call rsfilt(j1,i1,iflag,1)
        if(iflag.eq.0)goto 100
        tc(nctl)=0.0
        call tcolij(i1,j1,time,nctl,lc,tc,tw)
        if(tc(nctl).gt.1.0e-7) nctl=nctl+1
100     continue                
        enddo
300     if(ik.eq.2)goto 500
        j1=jc
500     enddo

        else

        j1=ic
        do ik=1,2
        do i=1,m32
        i1=i
        iflag=0
cfilter pairs that are of interest
        call rsfilt(j1,i1,iflag,1)
        if(iflag.eq.0)goto 200
        if(nctl.gt.nsize) stop 100000
        tc(nctl)=0.0
        call tcolij(i1,j1,time,nctl,lc,tc,tw)
        if(tc(nctl).gt.1.e-7) nctl=nctl+1
200     continue                
        enddo
        if(ik.eq.2)goto 600
        j1=jc
600     enddo
        endif
        if(tc(nctl).le.1.e-7) nctl=nctl-1
        do i=nctl+1,nsize
        do m=1,5
        lc(i,m)=0
        enddo
        tc(i)=0.
        tw(i)=0.
        enddo
        return
        end


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        function s0715(ss,ilo,i,the)
        parameter (KSZ1=30)      
        COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
        SAVE /FRPARA1/
c       pion- + p to k+ + sigma-,channel 4
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        if(ss.ge.1.9) goto 10
        IF(KFR(20).EQ.0)THEN
        si=0.25*(1.-0.75*(ss-1.691))
        ELSE
        si=VFR(20)      
        ENDIF
        goto 100
10      IF(KFR(20).EQ.0)THEN
        si=309.1*exp(max(-40.,-3.77*ss))
        ELSE
        si=VFR(20)
        ENDIF   
100     continue
        s0715=si
        return
        end



        function s07122(ss,ilo,i,the)
        parameter (KSZ1=30)      
        COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
        SAVE /FRPARA1/
c       pion- + p to k0 + lambda,channel 5
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        if(ss.ge.1.684) goto 10
        IF(KFR(20).EQ.0)THEN
        si=0.9/0.091*(ss-1.613)
        ELSE
        si=VFR(20)
        ENDIF
        goto 100
10      if(ss.ge.2.1) goto 20
        IF(KFR(20).EQ.0)THEN
        si=436.3*exp(max(-40.,-4.154*ss))
        ELSE
        si=VFR(20)
        ENDIF
        goto 100
20      IF(KFR(20).EQ.0)THEN
        si=0.314*exp(max(-40.,-0.301*ss))
        ELSE
        si=VFR(20)
        ENDIF
100     s07122=si
        return
        end



        function s07123(ss,ilo,i,the)
        parameter (KSZ1=30)      
        COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
        SAVE /FRPARA1/
c       pion- + p to k0 + sigma0,channel 6
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        if(ss.ge.1.722) goto 10
        IF(KFR(20).EQ.0)THEN
        si=10.6*(ss-1.689)
        ELSE
        si=VFR(20)
        ENDIF
        goto 100
10      if(ss.ge.3.) goto 20
        IF(KFR(20).EQ.0)THEN
        si=13.7*exp(max(-40.,-1.92*ss))
        ELSE
        si=VFR(20)
        ENDIF
        goto 100
20      IF(KFR(20).EQ.0)THEN
        si=0.188*exp(max(-40.,-0.611*ss))
        ELSE
        si=VFR(20)
        ENDIF
100     s07123=si
        return
        end


        
        
        function s1724(ss,ilo,i,the)
        parameter (KSZ1=30)      
        COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
        SAVE /FRPARA1/
c       pion+ + n to k+ + sigma0,channel 8
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        IF(KFR(20).EQ.0)THEN
        si=0.25*(s0715(ss,ilo,ii,the)+s07123(ss,ilo,ii,the)+
     c  s1713(ss,ilo,ii,the))
        ELSE
        si=VFR(20)
        ENDIF
100     s1724=si
        return
        end


        function s1727(ss,ilo,i,the)
        parameter (KSZ1=30)      
        COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
        SAVE /FRPARA1/
c       pion+ + n to k+ + lambda,channel 9
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        if(ss.ge.1.684) goto 10
        IF(KFR(20).EQ.0)THEN
        si=0.9*(ss-1.613)/0.091
        ELSE
        si=VFR(20)
        ENDIF
        goto 100
10      if(ss.ge.2.1) goto 20
        IF(KFR(20).EQ.0)THEN
        si=436.3*exp(max(-40.,-4.154*ss))
        ELSE
        si=VFR(20)
        ENDIF
        goto 100
20      IF(KFR(20).EQ.0)THEN
        si=0.314*exp(max(-40.,-0.301*ss))
        ELSE
        si=VFR(20)
        ENDIF
C*SA
c100    s1727=si*0.25
100     s1727=si
        return
        end


        function s1713(ss,ilo,i,the)
        parameter (KSZ1=30)      
        COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
        SAVE /FRPARA1/
c       pion+ + p to k+ + sigma+,channel 7
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        if(ss.ge.1.934) goto 10
        IF(KFR(20).EQ.0)THEN
        si=0.7*(ss-1.683)/0.218
        ELSE
        si=VFR(20)
        ENDIF
        goto 100
10      if(ss.ge.3.) goto 20
        IF(KFR(20).EQ.0)THEN
        si=60.26*exp(max(-40.,-2.31*ss))
        ELSE
        si=VFR(20)
        ENDIF
        goto 100
20      IF(KFR(20).EQ.0)THEN
        si=0.36*exp(max(-40.,-0.605*ss))
        ELSE
        si=VFR(20)
        ENDIF
100     s1713=si
        return
        end


        function s2325(ss,ilo,i,the)
c       pion0 + n to k+ + sigma-,channel 12
        ii=i
        s2325=s1724(ss,ilo,ii,the)
        return
        end


        function s2314(ss,ilo,i,the)
c       pion0 + p to k+ + sigma0,channel 10
        ii=i
        s2314=s1724(ss,ilo,ii,the)
        return
        end


        function s2317(ss,ilo,i,the)
c       pion0 + p to k+ + lambda,channel 11
        ii=i
        s2317=s1727(ss,ilo,ii,the)
        return
        end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine ppdelta(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'proc' to deal with pp -> ...
        parameter(nsize=100000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       l must be pion or k- & l1 might be pion,nucleon or antinucleon
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/count/isinel(400)
        COMMON/LUCIDAT2/KFMAXT,PARAM(20),WEIGH(400)
        SAVE /LUCIDAT2/
        common/sa10/csnn,cspin,cskn,cspipi,rcsit,ifram
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        if(isinel(174).eq.0)then
        sigma1=0.
        goto 522
        endif   
        the=pmas(lucomp(2224),1)+pmas(lucomp(2112),1)
        sigma1=snn(ss,ilo1,0,the)*WEIGH(46)
c       cross section of p + p to delta++  +  n

522     if(isinel(173).eq.0)then
        sigma2=0.
        goto 523
        endif
        the=pmas(lucomp(2214),1)+pmas(lucomp(2212),1)
        sigma2=snn(ss,ilo2,0,the)*WEIGH(45)
c       cross section of p+p to delta+ +p
523     if(ilo1.eq.0.and.ilo2.eq.0) goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6)  goto 13
        ik1=2224
        ik2=2112
        ic=174
c       p+p to delta++ +  n
        sigm12=sigma1+sigma2
        if(rlu(1).gt.sigma1/sigm12)then
        ik1=2214
        ik2=2212
        ic=173
c       cross section of p+p to delta+ +p
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/csnn/10.
        ioo=1
13      return
        end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine pndelta(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'proc' to deal with pn -> ...
        parameter(nsize=100000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       l must be pion or k- & l1 might be pion,nucleon or antinucleon
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/count/isinel(400)
        COMMON/LUCIDAT2/KFMAXT,PARAM(20),WEIGH(400)
        SAVE /LUCIDAT2/
        common/sa10/csnn,cspin,cskn,cspipi,rcsit,ifram
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        if(isinel(175).eq.0)then
        sigma1=0.
        goto 524
        endif   
        the=pmas(lucomp(2214),1)+pmas(lucomp(2112),1)
        sigma1=snn(ss,ilo1,0,the)*WEIGH(47)
c       cross section of p + n to delta+ +  n

524     if(isinel(176).eq.0)then
        sigma2=0.
        goto 525
        endif
        the=pmas(lucomp(2114),1)+pmas(lucomp(2212),1)
        sigma2=snn(ss,ilo2,0,the)*WEIGH(48)
c       cross section of p+n to delta0 +p
525     if(ilo1.eq.0.and.ilo2.eq.0) goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6)goto 13 
        ik1=2214
        ik2=2112
        ic=175
c       p+n to delta+ +  n
        sigm12=sigma1+sigma2
        if(rlu(1).gt.sigma1/sigm12)then
        ik1=2114
        ik2=2212
        ic=176
c       cross section of p+n to delta0 +p
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/csnn/10.
        ioo=1
13      return
        end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine nndelta(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'proc' to deal with nn-> ...
        parameter(nsize=100000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       l must be pion or k- & l1 might be pion,nucleon or antinucleon
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/count/isinel(400)
        COMMON/LUCIDAT2/KFMAXT,PARAM(20),WEIGH(400)
        SAVE /LUCIDAT2/
        common/sa10/csnn,cspin,cskn,cspipi,rcsit,ifram
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        if(isinel(177).eq.0)then
        sigma1=0.
        goto 527
        endif   
        the=pmas(lucomp(2114),1)+pmas(lucomp(2112),1)
        sigma1=snn(ss,ilo1,0,the)*WEIGH(49)
c       cross section of n + n to delta0 +  n

527     if(isinel(178).eq.0)then
        sigma2=0.
        goto 526

        endif
        the=pmas(lucomp(1114),1)+pmas(lucomp(2212),1)
        sigma2=snn(ss,ilo2,0,the)*WEIGH(50)
c       cross section of n+n to delta- +p
526     if(ilo1.eq.0.and.ilo2.eq.0) goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6)goto 13 
        ik1=2114
        ik2=2112
        ic=177
c       n+n to delta0 +  n
        sigm12=sigma1+sigma2
        if(rlu(1).gt.sigma1/sigm12)then
        ik1=1114
        ik2=2212
        ic=178
c       cross section of n+n to delta- +p
        endif
        lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm12/csnn/10.
        ioo=1
13      return
        end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine pip1(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'proc' to deal with pion- + p -> ...
        parameter(nsize=100000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       l must be pion or k- & l1 might be pion,nucleon or antinucleon
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
              COMMON/LUCIDAT2/KFMAXT,PARAM(20),WEIGH(400)
        common/count/isinel(400)
        common/sa10/csnn,cspin,cskn,cspipi,rcsit,ifram
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        if(isinel(11).eq.0)then
        sigma1=0.
        goto 401
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(3112),1)
        sigma1=s0715(ss,ilo1,0,the)*WEIGH(15)
c       cross section of pion- + p to k+ + sigma-
c       cross section is here in the unit of mb
401     if(isinel(12).eq.0)then
        sigma2=0.
        goto 402
        endif
        the=pmas(lucomp(311),1)+pmas(lucomp(3122),1)
        sigma2=s07122(ss,ilo2,0,the)*WEIGH(7)
c       cross section of pion- + p to k0 + lambda
402     if(isinel(13).eq.0)then
        sigma3=0.
        goto 403
        endif   
        the=pmas(lucomp(311),1)+pmas(lucomp(3212),1)

        sigma3=s07123(ss,ilo3,0,the)*WEIGH(16)
c       cross section of pion- + p to k0 + sigma0
403     if(isinel(147).eq.0)then
        sigma4=0.
        goto 203
        endif
c the--threshold energy of a reaction
        the=pmas(lucomp(1114),1)+pmas(lucomp(211),1)
        sigma4=sdelta(ss,ilo4,0,the)*WEIGH(57)
c       cross section of pion- + p to delta- + pi+
203     if(isinel(148).eq.0)then
        sigma5=0.
        goto 204
        endif
        the=pmas(lucomp(113),1)+pmas(lucomp(2112),1)
        sigma5=srho(ss,ilo5,0,the)*WEIGH(63)
c       cross section of pion- + p to rho0 + n
204     if(isinel(149).eq.0)then
        sigma6=0.
        goto 205
        endif
        the=pmas(lucomp(-213),1)+pmas(lucomp(2212),1)
        sigma6=srho(ss,ilo6,0,the)*WEIGH(64)
c       cross section of pion- + p to rho- + p
205     if(isinel(150).eq.0)then
        sigma7=0.
        goto 310
        endif
        the=pmas(lucomp(-211),1)+pmas(lucomp(2214),1)
        sigma7=sdelta(ss,ilo7,0,the)*WEIGH(58)
c       cross section of pion- + p to delta+ + pion-
310     if(isinel(151).eq.0)then
        sigma8=0.
        goto 805
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(2114),1)
        sigma8=sdelta(ss,ilo8,0,the)*WEIGH(58)
c       cross section of pion- + p to delta0 + pion0

805     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0
     &  .and.ilo4.eq.0 .and. ilo5.eq.0 .and. ilo6.eq.0.
     &  and.ilo7.eq.0.and.ilo8.eq.0)goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6 .and. sigma3.lt.1.e-6
     &         .and.sigma4.lt.1.e-6 .and. sigma5.lt.1.e-6 .and. 
     &          sigma6.lt.1.e-6.and. sigma7.lt.1.e-6.
     &          and.sigma8.lt.1.e-6)goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        sigma14=sigma13+sigma4
        sigma15=sigma14+sigma5
        sigma16=sigma15+sigma6
        sigma17=sigma16+sigma7
        sigma18=sigma17+sigma8
        s1=sigma1/sigma18
        s2=sigma12/sigma18
        s3=sigma13/sigma18
        s4=sigma14/sigma18
        s5=sigma15/sigma18
        s6=sigma16/sigma18
        s7=sigma17/sigma18
        rlus=rlu(1)
        if(rlus.le.s1)then
        ik1=321
        ik2=3112
c       3112 is the flavor code of sigma-
        ic=11
c       pion- + p to k+ + sigma-
        goto 416
        endif
        if(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=311
        ik2=3122
c       3122 is the Kf code of lambda
        ic=12
c       pion- + p to k0 + lambda
        goto 416
        endif
        if(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=311
        ik2=3212
c       3212 is the KF code of sigma0
        ic=13
c       pion- + p to k0 + sigma0
        goto 416
        endif
        if(rlus.gt.s3 .and. rlus.le.s4)then
        ik1=211
        ik2=1114
c       1114 is the Kf code of delta-
        ic=147
c       pion- + p to pi+ + delta-
        goto 416
        endif
        if(rlus.gt.s4 .and. rlus.le.s5)then
        ik1=113
        ik2=2112
c       113 is the Kf code of rho0
        ic=148
c       pion- + p to rho0 + n
        goto 416
        endif
        if(rlus.gt.s5 .and. rlus.le.s6)then
        ik1=-213
        ik2=2212
c       -213 is the Kf code of rho-
        ic=149
c       pion- + p to rho- + p
        goto 416
        endif
        if(rlus.gt.s6 .and. rlus.le.s7)then
        ik1=-211
        ik2=2214
c       2214 is the Kf code of delta+
        ic=150
c       pion- + p to delta++pion-
        goto 416
        endif
        ik1=111
        ik2=2114
c       2114 is the Kf code of delta0
        ic=151
c       pion- + p to delta0+pion0
        goto 416
416     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma18/cspin/10.
        ioo=1
13      return
        end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine pin3(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'proc' to deat with pion- + n -> ...
        parameter(nsize=100000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       l must be pion or k- & l1 might be pion,nucleon or antinucleon
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/count/isinel(400)
        COMMON/LUCIDAT2/KFMAXT,PARAM(20),WEIGH(400)
        SAVE /LUCIDAT2/
        common/sa10/csnn,cspin,cskn,cspipi,rcsit,ifram
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        if(isinel(152).eq.0)then
        sigma1=0.
        goto 222
        endif   
        the=pmas(lucomp(111),1)+pmas(lucomp(1114),1)
        sigma1=sdelta(ss,ilo1,0,the)*WEIGH(61)
c       cross section of pion- + n to delta-  +  pi0

222     if(isinel(153).eq.0)then
        sigma2=0.
        goto 223
        endif
        the=pmas(lucomp(-213),1)+pmas(lucomp(2112),1)
        sigma2=srho(ss,ilo2,0,the)*WEIGH(65)
c       cross section of pion- + n to rho- + n
223     if(isinel(154).eq.0)then
        sigma3=0.
        goto 228
        endif
        the=pmas(lucomp(-211),1)+pmas(lucomp(2114),1)
        sigma3=sdelta(ss,ilo3,0,the)*WEIGH(62)
c       cross section of pion- + n to delta0 + pion-
228     if(isinel(14).eq.0)then
        sigma4=0.
        goto 818
        endif
        ilo4=1
        the=pmas(lucomp(311),1)+pmas(lucomp(3112),1)            

        sigma4=s1724(ss,ilo4,0,the)*WEIGH(24)
        
c       cross section of pion- + n to k0 + sigma-
818     if(ilo1.eq.0.and.ilo2.eq.0.and.ilo3.eq.0
     &  .and.ilo4.eq.0) goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6.and. 
     &  sigma3.lt.1.e-6.and. sigma4.lt.1.e-6)goto 13 
        sigm12=sigma1+sigma2
        sigm13=sigm12+sigma3
        sigm14=sigm13+sigma4
        s1=sigma1/sigm14
        s2=sigm12/sigm14
        s3=sigm13/sigm14
        rlu1=rlu(1)
        if(rlu1.le.s1)then
        ik1=111
        ik2=1114
        ic=152
c       pion- + n to delta-  +  pi0
        elseif(rlu1.gt.s1 .and. rlu1.le.s2)then
        ik1=-213
        ik2=2112
        ic=153
c       cross section of pion- + n to rho- + n
        elseif(rlu1.gt.s2 .and. rlu1.le.s3)then
        ik1=-211
        ik2=2114
c       cross section of pion- + n to delta0 + pion-
        ic=154
        else
        ik1=311
        ik2=3112
        ic=14
c       cross section of pion- + n to k0 + sigma-
        endif
810     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigm14/cspin/10.
        ioo=1
13      return
        end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine pin2(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'proc' to deal with pion0 + n to ...
        parameter(nsize=100000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       l must be pion or k- & l1 might be pion,nucleon or antinucleon
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/count/isinel(400)
        COMMON/LUCIDAT2/KFMAXT,PARAM(20),WEIGH(400)
        SAVE /LUCIDAT2/
        common/sa10/csnn,cspin,cskn,cspipi,rcsit,ifram
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        if(isinel(172).eq.0)then
        sigma1=0.
        goto 216
        endif
        ilo1=1
        the=pmas(lucomp(-213),1)+pmas(lucomp(2212),1)           
        sigma1=srho(ss,ilo1,0,the)*WEIGH(72)
c       cross section of pion0 + n to rho- + p
216     if(isinel(168).eq.0)then
        sigma2=0.
        goto 217
        endif
c the--threshold energy of a reaction
        the=pmas(lucomp(-211),1)+pmas(lucomp(2214),1)
        sigma2=sdelta(ss,ilo2,0,the)*WEIGH(59)
c       cross section of pion0 + n to delta+  +  pi-
217     if(isinel(171).eq.0)then
        sigma3=0.
        goto 218
        endif
        the=pmas(lucomp(113),1)+pmas(lucomp(2112),1)
        sigma3=srho(ss,ilo3,0,the)*WEIGH(70)
c       cross section of pion0 + n to rho0 + n
218     if(isinel(169).eq.0)then
        sigma4=0.
        goto 807
        endif
        the=pmas(lucomp(211),1)+pmas(lucomp(1114),1)
        sigma4=sdelta(ss,ilo4,0,the)*WEIGH(60)
c       cross section of pion0 + n to delta- + pion+
807     if(isinel(19).eq.0)then
        sigma5=0.
        goto 811
        endif
        ilo5=1
        the=pmas(lucomp(311),1)+pmas(lucomp(3122),1)

        sigma5=s1724(ss,ilo5,0,the)*WEIGH(10)

c       cross section of pion0 + n to k0 + lambda is assumed
c to be  cross section of pion0 + p to k+ + lambda
811     if(isinel(18).eq.0)then
        sigma6=0.
        goto 816
        endif
        ilo6=1
        the=pmas(lucomp(321),1)+pmas(lucomp(3112),1)            

        sigma6=s2325(ss,ilo6,0,the)*WEIGH(22)
        
c       cross section of pion0 + n to k+ + sigma-
816     if(isinel(20).eq.0)then
        sigma7=0.
        goto 817
        endif
        ilo7=1
        the=pmas(lucomp(311),1)+pmas(lucomp(3212),1)            

        sigma7=s1724(ss,ilo7,0,the)*WEIGH(23)

c       cross section of pion0 + n to k0 + sigma0

817     if(isinel(170).eq.0)then
        sigma8=0.
        goto 827
        endif
        the=pmas(lucomp(111),1)+pmas(lucomp(2114),1)
        sigma8=sdelta(ss,ilo8,0,the)*WEIGH(60)
c       cross section of pion0 + n to delta0 + pion0
827     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0.and. 
     &  ilo4.eq.0.and. ilo5.eq.0
     &  .and. ilo6.eq.0 .and. ilo7.eq.0 .and. ilo8.eq.0)goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6 .and. sigma3.lt.1.e-6
     &          .and. sigma4.lt.1.e-6.and. sigma5.lt.1.e-6
     &  .and. sigma6.lt.1.e-6.and. sigma7.lt.1.e-6.
     &  and. sigma8.lt.1.e-6)goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        sigma14=sigma13+sigma4
        sigma15=sigma14+sigma5
        sigma16=sigma15+sigma6
        sigma17=sigma16+sigma7
        sigma18=sigma17+sigma8
        s1=sigma1/sigma18
        s2=sigma12/sigma18
        s3=sigma13/sigma18
        s4=sigma14/sigma18
        s5=sigma15/sigma18
        s6=sigma16/sigma18
        s7=sigma17/sigma18
        rlus=rlu(1)
        if(rlus.le.s1)then      
        ik1=-213
        ik2=2212
c       -213 is the flavor code of rho-
        ic=172
c       pion0 + n to p + rho-
        elseif(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=-211
        ik2=2214
c       2214 is the flavor code of delta+
        ic=168
c       pion0 + n to delta+  +  pi-
        elseif(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=113
        ik2=2112
c       113 is the Kf code of rho
        ic=171
c       pion0 + n to rho  +  n
        elseif(rlus.gt.s3 .and. rlus.le.s4)then
        ik1=211
        ik2=1114
c       1114 is the Kf code of delta-
        ic=169
        elseif(rlus.gt.s4 .and. rlus.le.s5)then
        ik1=311
        ik2=3122
c       3122 is the Kf code of lambda
c       pion0 + n to k0 + lambda
        ic=19
        elseif(rlus.gt.s5 .and. rlus.le.s6)then
        ik1=321
        ik2=3112
c       3112 is the Kf code of sigma-
        ic=18
        elseif(rlus.gt.s6 .and. rlus.le.s7)then
c       pion0 + n to k+ + sigma-
        ik1=311
        ik2=3212
c       3212 is the Kf code of sigma0
        ic=20
c       pion0 + n to k0 + sigma0
        else
        ik1=111
        ik2=2114
c       2114 is the Kf code of delta0
        ic=170
c       pion0 + n to pi0 + delta0
        endif
219     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma18/cspin/10.
        ioo=1
13      return
        end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine pin1(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'proc' to deal with pion+ + n to ...
        parameter(nsize=100000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       l must be pion or k- & l1 might be pion,nucleon or antinucleon
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/count/isinel(400)
        COMMON/LUCIDAT2/KFMAXT,PARAM(20),WEIGH(400)
        SAVE /LUCIDAT2/
        common/sa10/csnn,cspin,cskn,cspipi,rcsit,ifram
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        if(isinel(8).eq.0)then
        sigma1=0.
        goto 306
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(3212),1)            
        sigma1=s1724(ss,ilo1,0,the)*WEIGH(18)
c       cross section of pion+ + n to k+ + sigma0
306     if(isinel(9).eq.0)then
        sigma2=0.
        goto 207
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(3122),1)    
        sigma2=s1727(ss,ilo2,0,the)*WEIGH(8)
c       cross section of pion+ + n to k+ + lambda
207     if(isinel(158).eq.0)then
        sigma3=0.
        goto 208
        endif
c the--threshold energy of a reaction
        the=pmas(lucomp(-211),1)+pmas(lucomp(2224),1)
        sigma3=sdelta(ss,ilo3,0,the)*WEIGH(55)
c       cross section of pion+ + n to delta++  +  pi-
208     if(isinel(161).eq.0)then
        sigma4=0.
        goto 209
        endif
        the=pmas(lucomp(113),1)+pmas(lucomp(2212),1)
        sigma4=srho(ss,ilo4,0,the)*WEIGH(66)
c       cross section of pion+ + n to rho0 + p
209     if(isinel(162).eq.0)then
        sigma5=0.
        goto 210
        endif
        the=pmas(lucomp(213),1)+pmas(lucomp(2112),1)
        sigma5=srho(ss,ilo5,0,the)*WEIGH(67)
c       cross section of pion+ + n to rho+ + n
210     if(isinel(159).eq.0)then
        sigma6=0.
        goto 806
        endif
        the=pmas(lucomp(211),1)+pmas(lucomp(2114),1)
        sigma6=sdelta(ss,ilo6,0,the)*WEIGH(56)
c       cross section of pion+ + n to delta0 + pion+
806     if(isinel(10).eq.0)then
        sigma7=0.
        goto 814
        endif
        ilo7=1
        the=pmas(lucomp(311),1)+pmas(lucomp(3222),1)            

        sigma7=s1724(ss,ilo7,0,the)*WEIGH(19)
        
c       cross section of pion+ + n to k0 + sigma+

814     if(isinel(160).eq.0)then
        sigma8=0.
        goto 818
        endif
        ilo8=1
        the=pmas(lucomp(111),1)+pmas(lucomp(2214),1)            

        sigma8=sdelta(ss,ilo8,0,the)*WEIGH(19)

c       cross section of pion+ + n to pi0 + delta+

818     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0
     &  .and.ilo4.eq.0 .and. ilo5.eq.0.and.ilo6.eq.0
     &  .and.ilo7.eq.0.and.ilo8.eq.0)goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6 .and. sigma3.lt.1.e-6
     &         .and.sigma4.lt.1.e-6 .and. sigma5.lt.1.e-6.and. 
     &          sigma6.lt.1.e-6
     &  .and.sigma7.lt.1.e-6.and.sigma8.lt.1.e-6)goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        sigma14=sigma13+sigma4
        sigma15=sigma14+sigma5
        sigma16=sigma15+sigma6
        sigma17=sigma16+sigma7
        sigma18=sigma17+sigma8
        s1=sigma1/sigma18
        s2=sigma12/sigma18
        s3=sigma13/sigma18
        s4=sigma14/sigma18
        s5=sigma15/sigma18
        s6=sigma16/sigma18
        s7=sigma17/sigma18
        rlus=rlu(1)
        if(rlus.le.s1)then      
        ik1=321
        ik2=3212
c       3212 is the flavor code of sigma0
        ic=8
c       pion+ + n to k+ + sigma0
        elseif(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=321
        ik2=3122
c       3122 is the flavor code of lambda
        ic=9
c       pion+ + n to k+ + lambda
        elseif(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=-211
        ik2=2224
c       2224 is the Kf code of delta++
        ic=158
c       pion+ + n to pi- + delta++
        elseif(rlus.gt.s3 .and. rlus.le.s4)then
        ik1=113
        ik2=2212
c       113 is the Kf code of rho0
        ic=161
c       pion+ + n to rho0 + p
        elseif(rlus.gt.s4 .and. rlus.le.s5)then
        ik1=213
        ik2=2112
c       213 is the Kf code of rho+
        ic=162
c       pion+ + n to rho+ + n
        elseif(rlus.gt.s5 .and. rlus.le.s6)then
        ik1=211
        ik2=2114
c       2114 is the Kf code of delta0
        ic=159
        elseif(rlus.gt.s6 .and. rlus.le.s7)then
        ik1=311
        ik2=3222
c       3222 is the Kf code of sigma+
        ic=10
c       pion+ + n to k0 + sigma+
        else
        ik1=111
        ik2=2214
c       2214 is the Kf code of delta+
        ic=160
c       pion+ + n to pi0 + delta+
        endif
211     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma18/cspin/10.
        ioo=1
13      return
        end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine pip2(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'proc' to deal with pion+ + p to ...
        parameter(nsize=100000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       l must be pion or k- & l1 might be pion,nucleon or antinucleon
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/count/isinel(400)
        COMMON/LUCIDAT2/KFMAXT,PARAM(20),WEIGH(400)
        SAVE /LUCIDAT2/
        common/sa10/csnn,cspin,cskn,cspipi,rcsit,ifram
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        if(isinel(7).eq.0)then
        sigma1=0.
        goto 212
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(3222),1)
        sigma1=s1713(ss,ilo1,0,the)*WEIGH(17)
c       cross section of pion+ + p to k+ + sigma+
212     if(isinel(155).eq.0)then
        sigma2=0.
        goto 213
        endif
c the--threshold energy of a reaction
        the=pmas(lucomp(111),1)+pmas(lucomp(2224),1)
        sigma2=sdelta(ss,ilo2,0,the)*WEIGH(51)
c       cross section of pion+ + p to delta++  +  pi0
213     if(isinel(157).eq.0)then
        sigma3=0.
        goto 214
        endif
        the=pmas(lucomp(213),1)+pmas(lucomp(2212),1)
        sigma3=srho(ss,ilo3,0,the)*WEIGH(68)
c       cross section of pion+ + p to rho+ + p
214     if(isinel(156).eq.0)then
        sigma4=0.
        goto 804
        endif
        the=pmas(lucomp(211),1)+pmas(lucomp(2214),1)
        sigma4=sdelta(ss,ilo4,0,the)*WEIGH(52)
c       cross section of pion+ + p to delta+ + pion+
804     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0.
     &  and. ilo4.eq.0)goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6 .and. sigma3.lt.1.e-6
     &        .and. sigma4.lt.1.e-6 )goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        sigma14=sigma13+sigma4
        s1=sigma1/sigma14
        s2=sigma12/sigma14
        s3=sigma13/sigma14
        rlus=rlu(1)
        if(rlus.le.s1)then      
        ik1=321
        ik2=3222
c       3222 is the flavor code of sigma+
        ic=7
c       pion+ + p to k+ + sigma+
        elseif(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=111
        ik2=2224
c       2224 is the flavor code of delta++
        ic=155
c       pion+ + p to delta++  +  pi0
        elseif(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=213
        ik2=2212
c       213 is the Kf code of rho+
        ic=157
c       pion+ + p to rho+ + p
        else
        ik1=211
        ik2=2214
c       2214 is the Kf code of delta+
        ic=156
c       pion+ + p to delta+ + pion+
        endif
215     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma14/cspin/10.
        ioo=1
13      return
        end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        function sdelta(ss,ilo,i,the)
        parameter (KSZ1=30)      
        COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
        SAVE /FRPARA1/
c       pion- + p to pi0 + delta0
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        if(ss.ge.1.6941) goto 10
        IF(KFR(20).EQ.0)THEN
        si=-61.127+42.9365*ss
        ELSE
        si=VFR(21)
        ENDIF
        goto 100
10      IF(KFR(20).EQ.0)THEN
        sit=-0.0186959*ss**3+0.310359*ss**2-0.755106*ss+0.565481
        si=1.0/sit
        ELSE
        si=VFR(21)
        ENDIF
100     sdelta=si
        return
        end

        function srho(ss,ilo,i,the)
        parameter (KSZ1=30)      
        COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
        SAVE /FRPARA1/
c       pion- + p to pi0 + delta0
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        if(ss.ge.1.8837) goto 10
        IF(KFR(20).EQ.0)THEN
        si=-23.3607+13.9936*ss
        ELSE
        si=VFR(22)
        ENDIF
        goto 100
10      IF(KFR(20).EQ.0)THEN
        sit=0.331583*ss**3-1.86123*ss**2+3.81364*ss-2.50068
        si=1.0/sit
        ELSE
        si=VFR(22)
        ENDIF
100     srho=si
        return
        end

        function snn(ss,ilo,i,the)
        parameter (KSZ1=30)      
        COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
        SAVE /FRPARA1/
c       n + n to n + delta
c       parameterized x-section is not correct for n-n
        ilo=0
        si=0.
        ii=i
        if(ii.eq.1)goto 30
        if(ss.le.the) goto 100
30      ilo=1
        IF(KFR(20).EQ.0)THEN
        si=20*(ss-2.015)**2/(0.015+(ss-2.015)**2)
        ELSE
        si=VFR(23)
        ENDIF
100     snn=si
        return
        end


C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        subroutine pip3(l,l1,kl,kl1,ss,icp,lc,tc,tw,ioo)
c       a part of 'proc' to deal with pion0 + p to ...
        parameter(nsize=100000)
c       calculate particle production weight and fill up lc(i,3-5),tw(i).
c       l must be pion or k- & l1 might be pion,nucleon or antinucleon
c       tw : the ratio of cross section of (special inela.)/tot
        COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/count/isinel(400)
        COMMON/LUCIDAT2/KFMAXT,PARAM(20),WEIGH(400)
        SAVE /LUCIDAT2/
        common/sa10/csnn,cspin,cskn,cspipi,rcsit,ifram
        dimension lc(nsize,5),tc(nsize),tw(nsize)
        common/iloval/ilo1,ilo2,ilo3,ilo4,ilo5,ilo6,ilo7,ilo8,ilo9
        if(isinel(15).eq.0)then
        sigma1=0.
        goto 308
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(3212),1)
        sigma1=s2314(ss,ilo1,0,the)*WEIGH(20)
c       cross section of pion0 + p to k+ + sigma0
308     if(isinel(16).eq.0)then
        sigma2=0.
        goto 239
        endif
        the=pmas(lucomp(321),1)+pmas(lucomp(3122),1)
        sigma2=s2317(ss,ilo2,0,the)*WEIGH(9)
c       cross section of pion0 + p to k+ + lambda
239     if(isinel(163).eq.0)then
        sigma3=0.
        goto 220
        endif
c the--threshold energy of a reaction
        the=pmas(lucomp(211),1)+pmas(lucomp(2114),1)
        sigma3=sdelta(ss,ilo3,0,the)*WEIGH(53)
c       cross section of pion0 + p to delta0  +  pi+
220     if(isinel(165).eq.0)then
        sigma4=0.
        goto 240
        endif
        the=pmas(lucomp(213),1)+pmas(lucomp(2112),1)
        sigma4=srho(ss,ilo4,0,the)*WEIGH(69)
c       cross section of pion0 + p to rho+ + n
240     if(isinel(164).eq.0)then
        sigma5=0.
        goto 808
        endif
        the=pmas(lucomp(-211),1)+pmas(lucomp(2224),1)
        sigma5=sdelta(ss,ilo5,0,the)*WEIGH(54)
c       cross section of pion0 + p to delta++ + pion-
808     if(isinel(17).eq.0)then
        sigma6=0.
        goto 815
        endif
        ilo6=1
        the=pmas(lucomp(311),1)+pmas(lucomp(3222),1)            

        sigma6=s1724(ss,ilo6,0,the)*WEIGH(21)
        
c       cross section of pion0 + P to k0 + sigma+

815     if(isinel(166).eq.0)then
        sigma7=0.
        goto 835
        endif
        ilo7=1
        the=pmas(lucomp(113),1)+pmas(lucomp(2212),1)            
        sigma7=srho(ss,ilo7,0,the)*WEIGH(71)
c       cross section of pion0 + P to rho0 + p
835     if(isinel(167).eq.0)then
        sigma8=0.
        goto 845
        endif
        ilo8=1
        the=pmas(lucomp(111),1)+pmas(lucomp(2214),1)            
        sigma8=sdelta(ss,ilo8,0,the)*WEIGH(54)
c       cross section of pion0 + P to pi0 + delta+
845     if(ilo1.eq.0 .and. ilo2.eq.0 .and. ilo3.eq.0
     &  .and.ilo4.eq.0.and.ilo5.eq.0.and.ilo6.eq.0
     &  .and.ilo7.eq.0.and.ilo8.eq.0)goto 13
        if(sigma1.lt.1.e-6 .and. sigma2.lt.1.e-6 .and. sigma3.lt.1.e-6
     &         .and.sigma4.lt.1.e-6.and.sigma5.lt.1.e-6.and.sigma6.lt.1.e-6
     &  .and.sigma7.lt.1.e-6.and.sigma8.lt.1.e-6 )goto 13
        sigma12=sigma1+sigma2
        sigma13=sigma12+sigma3
        sigma14=sigma13+sigma4
        sigma15=sigma14+sigma5
        sigma16=sigma15+sigma6
        sigma17=sigma16+sigma7
        sigma18=sigma17+sigma8
        s1=sigma1/sigma18
        s2=sigma12/sigma18
        s3=sigma13/sigma18
        s4=sigma14/sigma18
        s5=sigma15/sigma18
        s6=sigma16/sigma18
        s7=sigma17/sigma18
        rlus=rlu(1)
        if(rlus.le.s1)then      
        ik1=321
        ik2=3212
        ic=15
c       3212 is the flavor code of sigma0
c       pion0 + p to k+ + sigma0
        elseif(rlus.gt.s1 .and. rlus.le.s2)then
        ik1=321
        ik2=3122
        ic=16
c       3122 is the flavor code of lambda
c       pion0 + p to k+ + lambda
        elseif(rlus.gt.s2 .and. rlus.le.s3)then
        ik1=211
        ik2=2114
c       2214 is the Kf code of delta+
        ic=163
c       pion0 + p to pi+ + delta0
        elseif(rlus.gt.s3 .and. rlus.le.s4)then
        ik1=213
        ik2=2112
c       213 is the Kf code of rho+
        ic=165
c       pion0 + p to rho+ + n
        elseif(rlus.gt.s4 .and. rlus.le.s5)then
        ik1=-211
        ik2=2224
c       2224 is the Kf code of delta++
        ic=164
        elseif(rlus.gt.s5 .and. rlus.le.s6)then
        ik1=311
        ik2=3222
c       3222 is the Kf code of sigma+
c       pion0 + p to k0 + sigma+
        ic=17
        elseif(rlus.gt.s6 .and. rlus.le.s7)then
        ik1=113
        ik2=2212
c       113 is the Kf code of rho0
        ic=166
c       pion0 + p to rho0 + p
        else
        ik1=111
        ik2=2214
c       2214 is the Kf code of delta+
c       pion0 + p to pi0 + delta+
        ic=167
        endif
221     lc(icp,3)=ik1
        lc(icp,4)=ik2
        lc(icp,5)=ic
        tw(icp)=sigma18/cspin/10.
        ioo=1
13      return
        end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

        SUBROUTINE RESCSET2(time)
c       put unaffected target (projectile)spectator nucleons back to a single entry
        PARAMETER (KSZJ=40000,KSZ1=30,mcludi=40000)
      COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
        common/sa4/tau(kszj),tlco(kszj,4)
        COMMON/FRPARA1/KFR(KSZ1),VFR(KSZ1)
      COMMON/FRINTN0/PLI0(2,4),AOP(KSZ1),IOP(KSZ1),NFR(KSZ1)
        common/wz/c17(mcludi,3),ishp(mcludi),tp(mcludi)
     c   ,coor(3)
      SAVE /FRPARA1/,/FRINTN0/,/LUDAT1/,/LUJETS/,/WZ/
        real ep(2,5)

        do 44 j=1,2
        do 44 jf=1,5
        ep(j,jf)=0.
44      continue
        jf=0
        jr1=0
        jr2=0
        jrn1=0
        jrn2=0
        if(kfr(24).eq.0)then
        k55=55
        k33=33
        else
        k55=-55
        k33=-33
        endif
        do j=1,n
        if(((k(j,3).eq.k55.or.k(j,3).eq.k33)
     &  .and.(k(j,2).eq.2112.or.k(j,2).eq.2212).and.k(j,4).ne.555.
     &  and.k(j,4).ne.333)
     &   .or.((k(j,3).eq.k55.or.k(j,3).eq.k33)
     &  .and.(k(j,2).eq.2112.or.k(j,2).eq.2212).and.(k(j,4).eq.555.
     &    or.k(j,4).eq.333)
     &  .and.((p(j,4)-0.938).lt.0.026.or.
     &  (p(j,4)-0.938).gt.0.375)))then
cif a particle is spectator nucleon which experiences neither inelastic 
c collisions nor within the grey particle region after elastic reactions 
C it will be put into a cluster if KFR(24)=0.
c(k(j,3)=55 target spectator,(k(j,3)=33 projectile spectator

        if(k(j,2).eq.2212.and.k(j,3).eq.33)jr1=jr1+1
        if(k(j,2).eq.2212.and.k(j,3).eq.55)jr2=jr2+1
        if(k(j,2).eq.2112.and.k(j,3).eq.33)jrn1=jrn1+1
        if(k(j,2).eq.2112.and.k(j,3).eq.55)jrn2=jrn2+1
        do 43 jy=1,5
        if(k(j,3).eq.33)ep(1,jy)=ep(1,jy)+p(j,jy)
        if(k(j,3).eq.55)ep(2,jy)=ep(2,jy)+p(j,jy)
43      continue
        else
C-change K0,K0ba to KL and Ks 
        kf=k(j,2)
        if(kf.eq.311 .or. kf.eq.-311)then
        rrlu=rlu(1)
        k(j,2)=130
        if(rrlu.gt.0.5)k(j,2)=310
        endif   
        jf=jf+1
        k(jf,1)=k(j,1)
        k(jf,2)=k(j,2)
        if(k(j,3).eq.33.or.k(j,3).eq.55)then
        k(jf,3)=k(j,3)
        k(jf,4)=k(j,4)
        k(jf,5)=k(j,5)  
        else
        k(jf,3)=0
        k(jf,4)=0
        k(jf,5)=0
        endif
        tp(jf)=tp(j)
        v(jf,4)=time
        v(jf,5)=0.
        c17(jf,1)=c17(j,1)
        c17(jf,2)=c17(j,2)
        c17(jf,3)=c17(j,3)
        v(jf,1)=c17(jf,1)
        v(jf,2)=c17(jf,2)
        v(jf,3)=c17(jf,3)
        do 40 jp=1,5
        p(jf,jp)=p(j,jp)
40      continue
        endif
        enddo
        if(ep(1,4).gt.1.e-4)then
        do 406 jp=1,5
        p(jf+1,jp)=ep(1,jp)
406     continue
        k(jf+1,1)=1
        k(jf+1,2)=10000+jr1
        k(jf+1,3)=jrn1
        k(jf+1,4)=0
        k(jf+1,5)=0
        jf=jf+1
        endif
        if(ep(2,4).gt.1.e-4)then
        do 407 jp=1,5
        p(jf+1,jp)=ep(2,jp)
407     continue
        k(jf+1,1)=1
        k(jf+1,2)=-10000-jr2
        k(jf+1,3)=jrn2
        k(jf+1,4)=0
        k(jf+1,5)=0
        jf=jf+1
        endif
        n=jf
        IOP(11)=jr1 
        IOP(12)=jr2
        return
        end
C*****************************************************************

c******************************************************
        subroutine copl(tt)
c       calculate the coordinate of the center of mass of
c the non-freeze-out system. The distance of a particle when checked
c if it freezes out or not is given with respect to this origin.
        parameter(mcludi=40000)
        parameter(KSZJ=40000)
                  COMMON/LUJETS/N,K(KSZJ,5),P(KSZJ,5),V(KSZJ,5)
         COMMON/LUDAT2/KCHG(500,3),PMAS(500,4),PARF(2000),VCKM(4,4)
        common/sa4/tau(kszj),tlco(kszj,4)
        common/wz/c17(mcludi,3),ishp(mcludi),tp(mcludi)
     c   ,coor(3)
        do ii=1,3
        coor(ii)=0.
        enddo
        samass=0.
        do 110 ii=1,n
c       if(tau(ii).gt.tt)goto 110
        if(ishp(ii).eq.0)goto 110
        kf=k(ii,2)
        amass=pmas(lucomp(kf),1)
        samass=samass+amass
        do 100 jj=1,3
        coor(jj)=coor(jj)+amass*c17(ii,jj)
100     continue
110     continue
        do ii=1,3
        coor(ii)=coor(ii)/max(0.14,samass)
        enddo
        return
        end
c********************************************************
C***********************DATA LUCIDATA*****************
        BLOCK DATA LUCIDATA
        COMMON/LUCIDAT1/KFACOT(50),DISDET(50),ISINELT(400)
        COMMON/LUCIDAT2/KFMAXT,PARAM(20),WEIGH(400)
        SAVE /LUCIDAT1/,/LUCIDAT2/
        DATA KFACOT/2212,2112,-2212,-2112,211,-211,111,-321,-311,
     &        3212,3112,3222,-3212,-3112,-3222,3122,-3122,311,
     &        321,3312,-3312,3322,-3322,3334,-3334,1114,2114,2214,
     &  2224,213,-213,113,18*0/
        DATA DISDET/0.5,0.5,0.5,0.5,21*0.,0.5,0.5,0.5,0.5,21*0./
        DATA ISINELT/178*1,22*0,178*1,14*0,8*1/
        DATA KFMAXT/32/
        DATA PARAM/40.,25.,35.,10.,3.0,0.85,1.0,0.0,0.1,2.0,
     &  0.16,0.01,28.0,5.6,0,0,0,0,0,0/         
                  DATA WEIGH/400*1.0/

        END
C******************************************************************
C...........Main switches and parameters...........................
C\item[KFACOT] flavor order of considered particles
C  \item[DISDET]  minimum allowable distance between two
C  particles,=0.5 if between two necleons,=0,otherwise
C  \item[ISINELT] switch of inelastic channel
C =0, i-th channel among interesting reactions is conceled,=1,opened
C \item[KFMAXT](D=12) KFMAXT kinds of particles are involved in rescattering
C PARAM(1)(D=40.0mb) totle cross-section of nucleon-nucleon 
C PARAM(2)(D=25.0mb)  totle cross-section of pi-nucleon 
C PARAM(3)(D=35.0mb) totle cross-section of K-nucleon 
C PARAM(4)(D=10.0mb)  totle cross-section of pi-pi
C PARAM(5)(D=3.mb)  cross-section of pi+pi -->K K 
C PARAM(6)(D=0.85) ratio of inelastic cross-section
C                              to  totle cross-section
C PARAM(7)(D=1.0fm)formation time at rest-frame of particles
C PARAM(8)(D=0.0fm)step of time elapse. D=0.0fm means no restriction to the time
c difference between two consecutive collisions.
C PARAM(9)(D=0.1) accuracy of four-momentum conservation
C PARAM(10)(D=2.0) size of effective rescattering region is PARAM(10)*radius of target.
c The condition is not used in this version, which means that rescattering stops when
c the collision pair list is empty.
C origin is center of target nucleus
C PARAM(11)(D=0.16fm^-3) nucleon density of nucleus
C PARAM(12)(D=0.01 GeV^2/c^2) The <Pt^2> for the Gaussian distribution of 
C spectator nucleons 
C PARAM(13)(D=28.0 mb) N Nbar annihilation cross section
C PARAM(14)(D=5.6 mb) N Lambdabar annihilation cross section
C nucleons

C@@@@@@@@@@@@@@@@@@@@@  END  @@@@@@@@@@@@@@@@@@@@@@@@

              


        
