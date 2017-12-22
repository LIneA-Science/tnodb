	implicit real*8 (a-h,o-z)
	CHARACTER*200 fichier_in, fichier_out
	parameter (nmax=100000)
	dimension xx(nmax), yy(nmax)
	dimension  x(nmax),  y(nmax), a(5), yout(nmax)
	integer*2 deg

 1000	format(a)
	zero= 0.d0
  	write(*,*) 'nome do arquivo de entrada: '
	read 1000, fichier_in

	k= 0
	xminmin=  1.d30
	xmaxmax= -1.d30
 1	write(*,*) 'Valores min e max de x: '
	read(*,*)  xmin, xmax
c
	open(unit=22,file=fichier_in,status='old',form='formatted')
	do i= 1, nmax
		read(22,*,err=99,end=99) xxx, yyy
		if (xxx.ge.xmin.and.xxx.le.xmax) then
			k= k+1
			xx(k)= xxx
			yy(k)= yyy
			if (xx(k).lt.xminmin) xminmin = xx(k)
			if (xx(k).gt.xmaxmax) xmaxmax = xx(k)
		endif
	enddo
C
C
99	continue
	if (k.eq.nmax) write(*,*)  'Nem todos os pontos lidos! Verifique o erro!'
	write(*,*)  'fichier ', fichier_in , ' lu ', k , ' points'
	close(22)

	write(*,*)  '1= outro intervalo'
	read(*,*)  iop
	if (iop.eq.1) go to 1

	write(*,*)  'Alisamento em ? pontos: '
	read(*,*)  ipas
	ipas2= ipas/2	

	npt= 0
	do i= 1+ipas2, k-ipas2, ipas
		som= 0.d0
		do j= i-ipas2, i+ipas2
			som= som + yy(j)
		enddo
		npt= npt+1
		x(npt)= xx(i)
		y(npt)= (som/dfloat(2*ipas2+1))
	enddo
	write(*,*)  npt, ' points medios'


	write(*,*)  'grau do polinomio a ajustar:  (<=4)'
	read(*,*)  deg

	call polyfit (x,y,npt,deg, a,yout)
	som= 0.d0
	do i= 1, npt
		som= som + (y(i) - yout(i))**2
	enddo
	som= dsqrt(som/dfloat(npt))
	write(*,*)  'rms= ', som

  	write(*,*) 'nome do arquivo de saida ? (Traco do polinomio)'
	read 1000, fichier_out
	open(unit=23,file=fichier_out,status='unknown',form='formatted')

  	write(*,*) 'nome do arquivo de saida (dados-polinomio)?'
	read 1000, fichier_out
	open(unit=24,file=fichier_out,status='unknown',form='formatted')

  	write(*,*) 'nome do arquivo de saida (dados/polinomio)?'
	read 1000, fichier_out
	open(unit=25,file=fichier_out,status='unknown',form='formatted')

	open(unit=22,file=fichier_in,status='old',form='formatted')
	npt= 0
	do i= 1, nmax
		read(22,*,err=98,end=98) xxx, yyy
		npt= npt + 1
		x(i)= xxx
		y(i)= yyy
	enddo
 98	continue
	close(22)
	nint= 1
	do i= 1, npt, nint
		yout(i)= 0.d0
		do j= 0, deg
			if (j.eq.0.and.x(i).eq.(0.d0)) then
				yout(i)= yout(i) + a(j+1)
			else
				yout(i)= yout(i) + a(j+1)*(x(i)**j)
			endif
 		enddo
		fac= (x(i)-xminmin)*(x(i)-xmaxmax) 
		if (fac.le.zero) then
		 write (23,*) x(i),  sngl(yout(i))
		 write (24,*) x(i),  sngl(y(i)-yout(i))
		 write (25,*) x(i),  sngl(y(i)/yout(i))
		endif
	enddo
	write(*,*) 'Um ponto cada ', nint, ' no arquivo de saida'

	stop
	end

***************************** fin de main **********************************

	subroutine polyfit (x,y,npt,deg, a,yout)
	implicit real*8 (a-h,o-z)
c
	CHARACTER*100 fichier
	parameter (nmax=100000)
	dimension    x(nmax),y(nmax),yout(nmax)
	integer*2 deg
	dimension xx(nmax),yy(nmax)
	dimension sig(nmax)
	dimension u(nmax,5),v(nmax,5),w(5),a(5),p(5)
c
	do i= 1, npt
		xx(i)= x(i)
		yy(i)= y(i)
	enddo
c
	ma= deg+1
	ndata= npt
	mp= ndata
	np= ma
c
	do i= 1,npt
		sig(i)= 1.d0
	enddo
c
      if (ma.eq.2) then
	    call droite(XX,YY,NDATA,A)
      else
      	call SVDFIT(XX,YY,SIG,NDATA,A,MA,U,V,W,MP,NP,CHISQ)
      endif
c
	write(*,*)  'Coeficientes do polinomio (do grau 0 ao maximo)'
c
	do i= 1,ma
	 	write(*,*)  a(i)
		write (10,*) a(i)
	enddo
c
	write(*,*)  ' '
c
	do i= 1, npt
	        call funcs(dble(x(i)),p,ma)
		yout(i)= 0.d0
		do j= 1, ma
			yout(i)= yout(i) + a(j)*p(j)
 		enddo
	enddo
c
c
	return
	end
c ***************************************************************************
c
	subroutine funcs(x,p,np)
c
	implicit real*8 (a-h,o-z)
	dimension p(np)
c
	p(1)= 1.d0
c
	do 11 j= 2, np
		p(j)= p(j-1)*x
 11	continue
c
	return
	end
c
c **************************************************************************
c
      SUBROUTINE SVDFIT(X,Y,SIG,NDATA,A,MA,U,V,W,MP,NP,CHISQ)
      implicit real*8 (a-h,o-z)
      PARAMETER(NMAX=100000,MMAX=5,TOL=1.D-16)
c
      DIMENSION X(NDATA),Y(NDATA),SIG(NDATA),A(MA),V(NP,NP),
     *    U(MP,NP),W(NP),B(NMAX),AFUNC(MMAX)
      DO 12 I=1,NDATA
        CALL funcs(X(I),AFUNC,MA)
        TMP=1.d0/SIG(I)
        DO 11 J=1,MA
          U(I,J)=AFUNC(J)*TMP
11      CONTINUE
        B(I)=Y(I)*TMP
12    CONTINUE
      CALL SVDCMP(U,NDATA,MA,MP,NP,W,V)
      WMAX=0.d0
      DO 13 J=1,MA
        IF(W(J).GT.WMAX)WMAX=W(J)
13    CONTINUE
      THRESH=TOL*WMAX
      DO 14 J=1,MA
        IF(W(J).LT.THRESH)W(J)=0.d0
14    CONTINUE
      CALL SVBKSB(U,W,V,NDATA,MA,MP,NP,B,A)
      CHISQ=0.d0
      DO 16 I=1,NDATA
        CALL funcs(X(I),AFUNC,MA)
        SUM=0.d0
        DO 15 J=1,MA
          SUM=SUM+A(J)*AFUNC(J)
15      CONTINUE
        CHISQ=CHISQ+((Y(I)-SUM)/SIG(I))**2
16    CONTINUE
      RETURN
      END
c
c **************************************************************************
c
      SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
      implicit real*8 (a-h,o-z)
      PARAMETER (NMAX=100)
c
      DIMENSION A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)
      G=0.0d0
      SCALE=0.0d0
      ANORM=0.0d0
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=0.0d0
        S=0.0d0
        SCALE=0.0d0
        IF (I.LE.M) THEN
          DO 11 K=I,M
            SCALE=SCALE+DABS(A(K,I))
11        CONTINUE
          IF (SCALE.NE.0.0d0) THEN
            DO 12 K=I,M
              A(K,I)=A(K,I)/SCALE
              S=S+A(K,I)*A(K,I)
12          CONTINUE
            F=A(I,I)
            G=-SIGN(DSQRT(S),F)
            H=F*G-S
            A(I,I)=F-G
            IF (I.NE.N) THEN
              DO 15 J=L,N
                S=0.0d0
                DO 13 K=I,M
                  S=S+A(K,I)*A(K,J)
13              CONTINUE
                F=S/H
                DO 14 K=I,M
                  A(K,J)=A(K,J)+F*A(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K= I,M
              A(K,I)=SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=0.0d0
        S=0.0d0
        SCALE=0.0d0
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K=L,N
            SCALE=SCALE+DABS(A(I,K))
17        CONTINUE
          IF (SCALE.NE.0.0d0) THEN
            DO 18 K=L,N
              A(I,K)=A(I,K)/SCALE
              S=S+A(I,K)*A(I,K)
18          CONTINUE
            F=A(I,L)
            G=-SIGN(DSQRT(S),F)
            H=F*G-S
            A(I,L)=F-G
            DO 19 K=L,N
              RV1(K)=A(I,K)/H
19          CONTINUE
            IF (I.NE.M) THEN
              DO 23 J=L,M
                S=0.0d0
                DO 21 K=L,N
                  S=S+A(J,K)*A(I,K)
21              CONTINUE
                DO 22 K=L,N
                  A(J,K)=A(J,K)+S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K=L,N
              A(I,K)=SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM=MAX(ANORM,(DABS(W(I))+DABS(RV1(I))))
25    CONTINUE
      DO 32 I=N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.0.0d0) THEN
            DO 26 J=L,N
              V(J,I)=(A(I,J)/A(I,L))/G
26          CONTINUE
            DO 29 J=L,N
              S=0.0d0
              DO 27 K=L,N
                S=S+A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=0.0d0
            V(J,I)=0.0d0
31        CONTINUE
        ENDIF
        V(I,I)=1.0d0
        G=RV1(I)
        L=I
32    CONTINUE
      DO 39 I=N,1,-1
        L=I+1
        G=W(I)
        IF (I.LT.N) THEN
          DO 33 J=L,N
            A(I,J)=0.0d0
33        CONTINUE
        ENDIF
        IF (G.NE.0.0d0) THEN
          G=1.0d0/G
          IF (I.NE.N) THEN
            DO 36 J=L,N
              S=0.0d0
              DO 34 K=L,M
                S=S+A(K,I)*A(K,J)
34            CONTINUE
              F=(S/A(I,I))*G
              DO 35 K=I,M
                A(K,J)=A(K,J)+F*A(K,I)
35            CONTINUE
36          CONTINUE
          ENDIF
          DO 37 J=I,M
            A(J,I)=A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J= I,M
            A(J,I)=0.0d0
38        CONTINUE
        ENDIF
        A(I,I)=A(I,I)+1.0d0
39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,30
          DO 41 L=K,1,-1
            NM=L-1
            IF ((DABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((DABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C=0.0d0
          S=1.0d0
          DO 43 I=L,K
            F=S*RV1(I)
            IF ((DABS(F)+ANORM).NE.ANORM) THEN
              G=W(I)
              H=DSQRT(F*F+G*G)
              W(I)=H
              H=1.0d0/H
              C= (G*H)
              S=-(F*H)
              DO 42 J=1,M
                Y=A(J,NM)
                Z=A(J,I)
                A(J,NM)=(Y*C)+(Z*S)
                A(J,I)=-(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
2         Z=W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.0.0d0) THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
          IF (ITS.EQ.30) STOP 'Sem convergencia em 30 iteracoes'
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0d0*H*Y)
          G=DSQRT(F*F+1.0d0)
          F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
          C=1.0d0
          S=1.0d0
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=DSQRT(F*F+H*H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 NM=1,N
              X=V(NM,J)
              Z=V(NM,I)
              V(NM,J)= (X*C)+(Z*S)
              V(NM,I)=-(X*S)+(Z*C)
45          CONTINUE
            Z=DSQRT(F*F+H*H)
            W(J)=Z
            IF (Z.NE.0.0d0) THEN
              Z=1.0d0/Z
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 NM=1,M
              Y=A(NM,J)
              Z=A(NM,I)
              A(NM,J)= (Y*C)+(Z*S)
              A(NM,I)=-(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L)=0.0d0
          RV1(K)=F
          W(K)=X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
      END
c
c **************************************************************************
c
      SUBROUTINE SVBKSB(U,W,V,M,N,MP,NP,B,X)
      implicit real*8 (a-h,o-z)
      PARAMETER (NMAX=100)
c
      DIMENSION U(MP,NP),W(NP),V(NP,NP),B(MP),X(NP),TMP(NMAX)
      DO 12 J=1,N
        S=0.d0
        IF(W(J).NE.0.d0)THEN
          DO 11 I=1,M
            S=S+U(I,J)*B(I)
11        CONTINUE
          S=S/W(J)
        ENDIF
        TMP(J)=S
12    CONTINUE
      DO 14 J=1,N
        S=0.d0
        DO 13 JJ=1,N
          S=S+V(J,JJ)*TMP(JJ)
13      CONTINUE
        X(J)=S
14    CONTINUE
      RETURN
      END
c
c ************************************************************************
c
	subroutine droite (x,y,ndata,a)
	implicit real*8 (a-h,o-z)
	parameter (nmax=100000)
	dimension a(5), x(nmax), y(nmax)
c
	sx= 0.d0
	sy= 0.d0
	sxy= 0.d0
	sx2= 0.d0
c
	do i=1, ndata
		sx= sx + x(i)
		sy= sy + y(i)
		sxy= sxy + x(i)*y(i)
		sx2= sx2 + x(i)*x(i)
	enddo
c
	denom= dfloat(ndata)*sx2 - sx*sx
	a(1)=  ( sx2*sy - sx*sxy)/denom					! coeff. de degre 0 (cste)
	a(2)=  ( dfloat(ndata)*sxy - sx*sy )/denom		! coeff. de degre 1 (pente)
c
c
	return
	end
c
c **********************************************************************
