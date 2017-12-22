c
c                          Programme ellipse_fit.for
c                          -------------------------
c
c fit des points du limbe per une ellipse. Les parametres libres sont 
c a, ksic, etac, e, P. Les cordes ne sont pas deplacees longitudinalement, ce
c qui semble plus stable que le programme ellipsefit.for. 
c
c Minimisation d'une fonction par methode du simplex. Le sous-programme
c qui deplace le simplex est tire de "Numerical Recipes", voir explications
c detaillees dans le livre, p. 289.
c
c La fonction a minimiser est FUNK (param), ou param est un tableau de NDIM
c valeurs correspondant aux differents parametres a ajuster. 
c
c Les parametres sont:
c
c	parafit(1)= demi grand-axe de l'ellipse 
c	parafit(2)= coordonnee ksi du centre de l'ellipse dans le plan du ciel 
c	parafit(3)= coordonnee eta du centre de l'ellipse dans le plan du ciel 
c   parafit(4)= aplatissement apparent de l'ellipse
c   parafit(5)= angle de position du demi petit-axe
c
c Sorties: 
c fort.51: meilleur limbe
c fort.52: centre du meilleur limbe
c fort.53: latitude, rayon
c fort.54: latitude, contributions au chi2: [(O-C)**2/sigma**2] individuels
c fort.55: i, ecarts radiaux (pour histogramme entre autres)
c fort.56: i, ecarts radiaux, min et max
c fort.57: P, epsilon (incremental)
c fort.58: P, R_eq (incremental)
c fort.59: P, sigma (incremental)
c fort.60: sigma, R_equiv (incremental)
c fort.61: R_equiv, epsilon (incremental)
c fort.62: i, rayon (pour histogramme entre autres)
c fort.63: i, rayons min et max
c fort.64: i, meilleur demi-grand-axe a (et non pas rayon equivalent) pour plot_i_r.i
c fort.65: latitude, rayons min et max

	parameter (MP=20, NP=20, NDIM= 5, NDIM1= 6)
	parameter (NCORDES=1000)
	dimension P(MP,NP)
	dimension parafit(NDIM), dparam(NDIM), Y(NDIM1)
	dimension dparam_err(NDIM)
c
	dimension eta(NCORDES), sigma(NCORDES) 
	dimension phi(2500), dr(2500)	! a utiliser dans correc
	real*4 ksi(NCORDES), ksic
	real*8 aa, lat, rho		! pour etre en accord avec limbe
	character*100 nom
	common/coor/npt,ksi,eta,sigma,B_pla,ifit,iplanete
	common/corr/phi,dr
	pi= 4.e+00*atan(1.e+0)
	rad= pi/180.d0
 100	format(a)
 
	open(unit=57,file='fort.57',status='unknown',form='formatted',position='append')
	open(unit=58,file='fort.58',status='unknown',form='formatted',position='append')
	open(unit=59,file='fort.59',status='unknown',form='formatted',position='append')
	open(unit=60,file='fort.60',status='unknown',form='formatted',position='append')
	open(unit=61,file='fort.61',status='unknown',form='formatted',position='append')

c	iplanete= 1        ! parametres gravitationnels de Jupiter
c	B_pla= 3.33677*rad ! elevation planetocentrique le 10/10/99

c	iplanete= 5        ! parametres gravitationnels de Pluton
c	B_pla= -28.46*rad  ! elevation planetocentrique le 20/07/02

c	iplanete= 6        ! parametres gravitationnels de Titan
c	B_pla= -23.53*rad  ! elevation planetocentrique le 14/11/03

c	iplanete= 7        ! parametres gravitationnels de Charon
c	B_pla= -34.19*rad  ! elevation planetocentrique le 11/07/05

c	iplanete= 8        ! parametres gravitationnels de Titania
c	B_pla= -24.17*rad  ! elevation planetocentrique le 08/09/01
c!	B_pla= +62.00*rad  ! elevation planetocentrique en prenant pour "pole" le point sub_Uranus 

c	iplanete= 9        ! parametres gravitationnels de Ceres
c	B_pla= +5.9381*rad ! elevation planetocentrique le 17/08/10

c	iplanete= 10        ! parametres gravitationnels d'Eris
c	B_pla= +0.0*rad		! elevation planetocentrique inconnue

c	iplanete= 11        ! parametres gravitationnels de Charon
c	B_pla= -46.5759*rad	! elevation planetocentrique le 04/06/11

	iplanete= 12        ! parametres gravitationnels de Chariklo
	B_pla=  00.d0*rad	! d'apres fit des anneaux arccsin(1-0.44414940) 
!	B_pla=  33.77d0*rad	! d'apres fit des anneaux arccsin(1-0.44414940) 

c Dans le cas de Jupiter, on lit les corrections a apporter au limbe, dues
c aux vents zonaux. Il faut alors activer le sous-programme "correc"
c
c	open(unit=9,file='lat_dr.dat',status='old',form='formatted')
c	npt= 0
c	do i= 1, 2500
c		read(9,*,err=97,end=97) phi(i), dr(i)
c		npt= npt + 1
c	enddo
c 97	close(9)
c	write (*,*) npt, ' points lus dans lat_dr.dat'
c
c ..........................................................................
c
c lecture des coordonnees des points d'immersion et d'emersion
c
	write (*,*) 'nom du fichier contenant les coordonnees des points'
	write (*,*) 'd''immersion-emersion ?'
	read (*,*) nom
	open(unit=10,file=nom,status='old',form='formatted')
	npt= 0
	do i= 1, NCORDES
		read(10,*,err=98,end=99) ksi(i), eta(i)
		npt= npt + 1
 98		continue
	enddo
 99	close(10)
	write (*,*) npt, 'points pris en compte'
c
	write (*,*) 'Type de l''ajustement? (1= ellipse, 2= isopotentielle)'
	read (*,*) ifit
c
	write (*,*)'sigma associe a chaque point? (en km, ds la direction RADIALE)'
	do i= 1, npt
		read (*,*) sigma(i)
	enddo
c
c fin de la lecture des coordonnees des points d'immersion et d'emersion
c
c ..........................................................................
c
c initialisation du simplex
c
c
	write (*,*) ' '
	write (*,*) 'Nombre de parametres consideres: ', NDIM
	write (*,*) ' '
c
 	write (*,*) 'Tolerance ? (variation relative minimum de la fonction)'
     	write (*,*) 'au sommet du simplex'
	read (*,*) FTOL
c
c Choix des parametres initiaux et de leur variation.
c Les points P(j,i) forment un matrice de NDIM+1 lignes et de NDIM colonnes.
c Chaque LIGNE est une vecteur de NDIM coordonnees dans l'espace des parametres
c a ajuster. Ce vecteur definit l'un des sommets du simplex, il y a donc NDIM+1
c sommets, chacun indexe par j. On definit d'abord l'un des sommets du simplex, 
c par exemple P(1,i), ou i explore l'ensemble des parametres. On definit 
c ensuite les autres sommets du simplex P(j,i) en prenant P(j,i)= P(1,i), SAUF
c si i= j-1, auquel cas P(j,i)= P(1,i) + dparam(i), i.e. on fait varier la
c (j-1) eme coordonnee de P(j,i) d'une quantite dparam(i), valeur typique de la
c variation du i eme parametre au cours de l'exploration d'ajustement.
c
	write (*,*) 'rayon EQUATORIAL et variation ? (km)'
	read (*,*) P(1,1), dparam(1)
	write (*,*) 'ksi du centre du cercle et variation ? (km)'
	read (*,*) P(1,2), dparam(2)
	write (*,*) 'eta du centre du cercle et variation ? (km)'
	read (*,*) P(1,3), dparam(3)
	write (*,*) 'aplatissement et variation ?' 
	read (*,*) P(1,4), dparam(4)
	write (*,*) 'angle de position et variation ? (deg)' 
	read (*,*) P(1,5), dparam(5)

	write (*,*) 'ATTENTION: dans le cas d''un fit avec isopotentielle,'
	write (*,*) 'il faut affecter une barre d''erreur NULLE a l''aplatissement'
	write (*,*) 'car ce dernier est DEDUIT du rayon equatorial'
C
C ON CALCULE LE NOMBRE DE PARAMETRE LIBRES
C
	NPARAM= 0
	DO I= 1, NDIM
		IF (DPARAM(I).NE.0.D0) NPARAM= NPARAM+1
	ENDDO

c---------------------------------------------------------------
c On explore differentes conditions initiales
c
c	write (*,*) 'Conditions initiales: a_min, a_max et pas?'
c	read (*,*) a_min, a_max, a_pas
c	a_ini= a_min
c	do while(a_ini.le.a_max)
c	P(1,1)= a_ini
c---------------------------------------------------------------

c
c Initialisation des sommets du simplex
c
	do j= 1, NDIM+1
		do i= 1, NDIM
			if(i.ne.(j-1)) P(j,i)= P(1,i)
			if(i.eq.(j-1)) P(j,i)= P(1,i) + dparam(i)
		enddo
	enddo
c
c Calcul de FUNK aux sommets du simplex 
c
	do j= 1, NDIM+1 
		do i= 1, NDIM
			parafit(i)= P(j,i)
		enddo
		y(j)= FUNK (parafit)
	enddo
c
c fin de l'initialisation du simplex
c
c ..........................................................................
c
c Appel du sous-programme de minimisation et resultats finaux
c	
	call amoeba (P,Y,MP,NP,NDIM,FTOL,FUNK,ITER)
	write (*,*) ' '	
	write (*,*) 'Rayon EQUATORIAL (km):',  P(1,1)
	write (*,*) 'Aplatissement:', P(1,4), '(apparent si fit elliptique,'
	write (*,*) 'reel si fit isopotentielle)'
	write (*,*) 'Rayon polaire (km):',  P(1,1)*(1. - P(1,4)), 'Idem'
	write (*,*) 'Centre de l''ellipse (km):', P(1,2), P(1,3)
	write (*,*) 'Angle de position du pole (deg):', P(1,5)
	write (*,*) 'Residu minimum:', y(1)
	write (*,*) '[ sqrt( Som (O-C)**2/(npt*sigma**2) ]'
c
c Calcul des barres d'erreur
c
	write (*,*) 'Nombre de points pris en compte:', npt
	write (*,*) 'Nombre de parametres libres:', nparam 
	nlib= npt-nparam
	write (*,*) 'Nombre de degres de liberte:', nlib
	if (nlib.eq.0) then
	 write (*,*) 'Nbre de points, chi2', npt, float(npt)*(y(1)**2),
     *	'[ Som (O-C)**2/sigma**2 ]'
	 write (*,*) 'chi2 reduit infini'
	else
	 write (*,*) 'Nbre de points, chi2', npt, float(npt)*(y(1)**2),
     *	'[ Som (O-C)**2/sigma**2 ]'
	write (*,*) 'chi2 par degre de liberte:', (float(npt)/float(nlib))*(y(1)**2), '[ Som (O-C)**2/(nlib*sigma**2) ]'
	endif
	
	write (*,*) 'Nombre d''iterations:', ITER
	
c---------------------------------------------------------------
c
c Fin de la boucle d'exploration des conditions initiales
c 
c	write (40,*) P(1,1), y(1)
c	write (41,*) P(1,4), y(1)
c	write (42,*) P(1,4), P(1,1)
c	a_ini= a_ini + a_pas
c	enddo
c---------------------------------------------------------------

	do i= 1, ndim
		parafit(i)= P(1,i)
	enddo
c
c on definit les pas pour explorer les parametres autour du meilleur fit
c
c NB avec un fit isopotentiel, il faut mettre un pas NUL pour l'aplatissement
c car il est IMPOSE par la valeur de a 
c

	write (*,*) 'entrer les', ndim, ' pas pour explorer les parametres '
	do i= 1, ndim
		read (*,*) dparam_err(i)
	enddo
	write (*,*) ' '
c
c	dchi2= 3.53	! definit a combien de sigma on veut les parametres
c	write (*,*) ' '
c	write (*,*) 'Barres d''erreur donnees a dchi2:', dchi2,
c     *	'(68% avec 3 param libres)'
c
c	write (*,*) 'Niveau de confiance de 95.4% avec 3 param libres'
c	dchi2= 8.02	! definit a combien de sigma on veut les parametres
	write (*,*) 'Niveau de confiance de 99.7% avec 3 param libres: dchi2= 14.2 ?'
c	dchi2= 14.2d0
	dchi2= 1.d0
c
c Chi2 est ici somme ( (O-C)/(sigma**2) )**2
c
	chi2= ( (funk(parafit))**2 )*float(npt)
c	chi2max= chi2 + dchi2			! en supposant que le chi2 est representatif
	chi2max= float(npt) + dchi2		! en obligeant le chi2 minimum a etre egal a npt
	write (*,*) 'Barres d''erreur donnees a Dchi2:', dchi2,
     *	'(chi2max=', chi2max, ')'

	do i= 1, ndim					! boucle sur les ndim parametres
		if (dparam_err(i).eq.0.d0) go to 1
		do while (chi2.le.chi2max)
			parafit(i)=   parafit(i) + dparam_err(i)
			chi2= ( (funk(parafit))**2 )*float(npt)
		enddo
 1	    continue
		erreur= parafit(i)-p(1,i)
		write (*,*)'Parametre:',i,'=', p(1,i),' +/-', erreur

		do k= 1, ndim		 ! on re-initialise les param. au minimum
			parafit(k)= p(1,k)
		enddo
		chi2= ( (funk(parafit))**2 )*float(npt)
	enddo
 

c
c On trace le dernier fit, i.e. la derniere ellipse, rho en fonction de theta. 
c theta est l'angle de position du point courant (par rapport a la direction N). Donc, 
c par ex. theta = posi indique la region polaire nord (demi petit-axe apparent).
c NB la definition de theta est differente de celle utilisee dans FUNK, ou theta=
c angle a partir de l'axe des ksi (direction E)
c
	a=       P(1,1)
	ksic=    P(1,2)
	etac=    P(1,3)
	epsilon= P(1,4)
	posi=    P(1,5)*rad
	cosP=    cos(posi)
	sinP=    sin(posi)
	cosB=    cos(B_pla)
	sinB=    sin(B_pla)
	b=       a*(1.d0 - epsilon)

	write (*,*) ' '
	write (*,'(A,f8.4,1x,f8.4)') 'Elevation B et P pole (deg): ', B_pla/rad, P(1,5)
	write (*,*) 'Rayon equivalent apparent:', sqrt(a*b)
	write (*,*) 'Aplatissement et rayon polaire reels'
	write (*,*) '(cas du fit elliptique seulement)'
	epsilon_reel= 1. - (sqrt( (1-epsilon)**2 - sinB**2 )/cosB )
	b_reel= a*(1-epsilon_reel)
	write (*,*) epsilon_reel, b_reel

	theta=		-180.
	pas=		0.25
	thetamax=	180. 
	do while(theta.le.thetamax) ! theta est l'ANGLE DE POSITION du point courant
		thetar= theta*rad

		if (ifit.eq.1) then
		 rho= a*b/sqrt((a*cos(thetar-posi))**2+(b*sin(thetar-posi))**2)
		 xx= rho*sin(thetar) + ksic
		 yy= rho*cos(thetar) + etac
		endif

		if (ifit.eq.2) then
		 write (*,*) ' '
		 write (*,*) 'ATTENTION: Pas de formule encore disponible'
		 write (*,*) 'pour tracer le limbe!'
		 go to 2
		endif

		write (51,*) xx, yy
		theta= theta + pas
	enddo
 2	continue

c centre de l'ombre (deux fois pour pouvoir plotter comme corde)
	write (52,*) sngl(ksic), sngl(etac) 
	write (52,*) sngl(ksic), sngl(etac) 
c
c On trace les points observes (extremites de chaque corde)
c
	dispersion= 0. 
	npt_selec=  0 ; seuil= 1000.d0		! nombre de points selectionnes, ie dont le sigma est < seuil
	do i= 1, npt
		xx= ksi(i)-ksic
		yy= eta(i)-etac
		rr= sqrt(xx**2 + yy**2)
		sin_lat= ( (xx*sinP + yy*cosP)*cosB )/rr
		lat= asin(sin_lat)

		if (ifit.eq.1) then
		 theta= atan(xx/yy)
		 if (yy.le.0.) theta= theta + pi
		 rho= a*b/sqrt((a*cos(theta-posi))**2+(b*sin(theta-posi))**2)
		endif

		if (ifit.eq.2) then
		 aa= a							! aa en double precision
c		 call correc (0., d_r)
c		 aa= aa - dble(d_r)				! aa sans vents zonaux
		 call limbe (aa, lat, rho)
c		 call correc (sngl(lat),d_r)	! effets vents zonaux
c		 rho= rho + d_r
		endif

c ecarts radiaux en km:
		write (53,*) sngl(lat)/rad,  rr				, i
		write (53,*) sngl(lat)/rad,  rr				, i			! on repete pour pouvoir utiliser plot_pla ! 
		write (55,*) i,             (rr-sngl(rho))				! pour faire histogramme de Delta r
		write (55,*) i,             (rr-sngl(rho))				! on repete pour pouvoir utiliser plot_pla ! 
		write (56,*) i,				(rr-sngl(rho)) - sigma(i)	! pour tracer les barres d'erreur
		write (56,*) i,				(rr-sngl(rho)) + sigma(i)	!

		write (62,*) i,              rr							! pour tracer les rayons
		write (62,*) i,              rr							! on repete pour pouvoir utiliser plot_pla ! 
		write (63,*) i,				 rr            - sigma(i)	! pour tracer les barres d'erreur sur (i,Dr)
		write (63,*) i,				 rr            + sigma(i)	!
		write (64,*) i,				 P(1,1)

		write (65,*) sngl(lat)/rad,	 rr - sigma(i)				! pour tracer les barres d'erreur sur (lat,Dr)
		write (65,*) sngl(lat)/rad,	 rr + sigma(i)				!
		
c contributions au chi2:		
		write (54,*) sngl(lat)/rad, ((rr-sngl(rho))/sigma(i))**2, i
c		write (54,*) sngl(lat)/rad, rr-sngl(rho)+sigma(i)
c		write (54,*) sngl(lat)/rad, rr-sngl(rho)-sigma(i)
		if (sigma(i).lt.seuil) then
		 npt_selec= npt_selec + 1
		 dispersion= dispersion + (rr-sngl(rho))**2
		endif
	enddo
	dispersion= sqrt(dispersion/float(npt_selec))
	write (*,'(A,I2,A,f7.3)') 'dispersion sur les ', npt_selec, ' points selectionnes (km)=', dispersion
	R_equiv=  sqrt(a*b)											! rayon equivalent APPARENT
	write (57,*) P(1,5), epsilon	
	write (58,*) P(1,5), R_equiv
	write (59,*) P(1,5), dispersion
	write (60,*) dispersion, R_equiv
	write (61,*) R_equiv, epsilon								! abs(epsilon) si on peut avoir epsilon > 0 quand a < b 
																! (que l'on retrouvre avec a > b et en faisant P ---> P + 90)
        open(81,file='fort.81',status='unknown',action='write',form='formatted',position='append')
	write(81,*) P(1,1), P(1,4), P(1,2), P(1,3), P(1,5), P(1,1)*(1-P(1,4)), (float(npt)/float(nlib))*(y(1)**2), dispersion
        close(81)
	
	stop
	end
c
c				FIN DE MAIN
c
c **********************************************************************
C
C
C
      SUBROUTINE AMOEBA(P,Y,MP,NP,NDIM,FTOL,FUNK,ITER)
      PARAMETER (NMAX=20,ALPHA=1.0,BETA=0.5,GAMMA=2.0,ITMAX=500)
      DIMENSION P(MP,NP),Y(MP),PR(NMAX),PRR(NMAX),PBAR(NMAX)
      MPTS=NDIM+1
      ITER=0
 1    continue
C      write (*,*) 'ITERATION NO=', ITER
C      write (*,*) 'VALEURS DES PARAMETRES'
      DO 10 I= 1,NDIM
C	write (*,*) 'PARAMETRE', I, ' =', P(1,I)
10    CONTINUE 
      ILO=1
      IF(Y(1).GT.Y(2))THEN
        IHI=1
        INHI=2
      ELSE
        IHI=2
        INHI=1
      ENDIF
      DO 11 I=1,MPTS
        IF(Y(I).LT.Y(ILO)) ILO=I
        IF(Y(I).GT.Y(IHI))THEN
          INHI=IHI
          IHI=I
        ELSE IF(Y(I).GT.Y(INHI))THEN
          IF(I.NE.IHI) INHI=I
        ENDIF
11    CONTINUE
      RTOL=2.*ABS(Y(IHI)-Y(ILO))/(ABS(Y(IHI))+ABS(Y(ILO)))
C      write (*,*) 'VALEUR MINIMALE DE LA FONCTION', Y(ILO)
      IF(RTOL.LT.FTOL)RETURN
      IF(ITER.EQ.ITMAX) THEN
         WRITE(*,*) 'Amoeba exceeding maximum iterations.'
         RETURN
      ENDIF
      ITER=ITER+1
      DO 12 J=1,NDIM
        PBAR(J)=0.
12    CONTINUE
      DO 14 I=1,MPTS
        IF(I.NE.IHI)THEN
          DO 13 J=1,NDIM
            PBAR(J)=PBAR(J)+P(I,J)
13        CONTINUE
        ENDIF
14    CONTINUE
      DO 15 J=1,NDIM
        PBAR(J)=PBAR(J)/NDIM
        PR(J)=(1.+ALPHA)*PBAR(J)-ALPHA*P(IHI,J)
15    CONTINUE
      YPR=FUNK(PR)
      IF(YPR.LE.Y(ILO))THEN
        DO 16 J=1,NDIM
          PRR(J)=GAMMA*PR(J)+(1.-GAMMA)*PBAR(J)
16      CONTINUE
        YPRR=FUNK(PRR)
        IF(YPRR.LT.Y(ILO))THEN
          DO 17 J=1,NDIM
            P(IHI,J)=PRR(J)
17        CONTINUE
          Y(IHI)=YPRR
        ELSE
          DO 18 J=1,NDIM
            P(IHI,J)=PR(J)
18        CONTINUE
          Y(IHI)=YPR
        ENDIF
      ELSE IF(YPR.GE.Y(INHI))THEN
        IF(YPR.LT.Y(IHI))THEN
          DO 19 J=1,NDIM
            P(IHI,J)=PR(J)
19        CONTINUE
          Y(IHI)=YPR
        ENDIF
        DO 21 J=1,NDIM
          PRR(J)=BETA*P(IHI,J)+(1.-BETA)*PBAR(J)
21      CONTINUE
        YPRR=FUNK(PRR)
        IF(YPRR.LT.Y(IHI))THEN
          DO 22 J=1,NDIM
            P(IHI,J)=PRR(J)
22        CONTINUE
          Y(IHI)=YPRR
        ELSE
          DO 24 I=1,MPTS
            IF(I.NE.ILO)THEN
              DO 23 J=1,NDIM
                PR(J)=0.5*(P(I,J)+P(ILO,J))
                P(I,J)=PR(J)
23            CONTINUE
              Y(I)=FUNK(PR)
            ENDIF
24        CONTINUE
        ENDIF
      ELSE
        DO 25 J=1,NDIM
          P(IHI,J)=PR(J)
25      CONTINUE
        Y(IHI)=YPR
      ENDIF
      GO TO 1
      END
C
C				FIN DE AMOEBA
C
C **************************************************************************

	function funk (param)
c
c Fonction a minimiser. Cette fonction depend des parametres param(i), ainsi
c que de variables eventuelles qui sont passees dans la fonction via des ordres 
c COMMON.
c
	parameter (nmax=5)									! ATTENTION: nmax doit etre egal a NDIM !
	parameter (NCORDES=1000)
	common/coor/npt,ksi,eta,sigma,B_pla,ifit,iplanete
	dimension param(nmax)
	dimension eta(NCORDES), sigma(NCORDES) 
	real*4 ksi(NCORDES), ksic 
	real*8 aa, lat, rayon ! pour etre compatible avec limbe...

	pi= acos(-1.0)
	rad= pi/180.
	cosB=  cos(B_pla)	
	call parametres(iplanete)

	a=     param(1)
	ksic=  param(2)
	etac=  param(3)
	eps=   param(4)
	posi=  param(5)*rad
	b=     a*(1.-eps)

	cosP=  cos(posi)
	sinP=  sin(posi)

	funk= 0.
c
c Calcul de funk avec un limbe elliptique:
c
	if (ifit.eq.1) then 
	 do i= 1, npt
		x= ksi(i)-ksic
		y= eta(i)-etac
		theta= atan(y/x)
		if (x.le.0.) theta= theta + pi
		rho= a*b/sqrt((a*sin(theta+posi))**2 + (b*cos(theta+posi))**2)
		xe= ksic + rho*cos(theta)
		ye= etac + rho*sin(theta)
		funk= funk + ((ksi(i)-xe)**2 + (eta(i)-ye)**2)/sigma(i)**2
	 enddo

	 funk= sqrt(funk/float(npt))
	 return
	endif
c
c Calcul avec un limbe defini par J2, J4, J6 et q de la planete
c
c Calcul a priori de l'aplatissement:
c
	if (ifit.eq.2) then 
	 aa=  dble(a) 				! aa en double precision
c	 call correc (0., d_r)
c	 aa= aa - dble(d_r)			! aa s'il n'y avait pas de vents

	 lat= dble(90.*rad)	
	 call limbe (aa, lat, rayon)
c	 call correc (sngl(lat),d_r)		! effets vents zonaux
c	 rayon= rayon + d_r
	 param(4)= (a-sngl(rayon))/a

	 do i= 1, npt
		x=     ksi(i) - ksic
		y=     eta(i) - etac
		r_obs= sqrt( x**2 + y**2 )
		sin_lat= ((x*sinP + y*cosP)*cosB)/r_obs
		lat=  dasin(dble(sin_lat))
		call limbe (aa, lat, rayon) 
		r_theo= rayon				! simple precision 
c		call correc (sngl(lat),d_r)		! effet vents zonaux
c		r_theo= r_theo + d_r
		funk= funk + ((r_obs - r_theo)**2)/sigma(i)**2
	 enddo

	 funk= sqrt(funk/float(npt))
	 return
	endif

	end
c
c				FIN DE FUNK
c
c **************************************************************************

	subroutine parametres(iplanete)
        implicit real*8 (a-h,o-z)
        real*8 MT,J(6)
        common /grav/GM,a,J,q
	G = 6.67259d-11 ! Site JPL/Horizon
	MT= 5.9746d24	! masse de la Terre en kg (precision?)
c 
c attention, les unites sont les km (a) et km^3 s^-2 (GM) ci-dessous!
c
c Jupiter (Cf. Table de Bill), valeurs a 1 bar 
c
	if (iplanete.eq.1) then
		GM= (G*MT*317.735d0)/1.d9	! en km^3 s^-2 !
		a = 71492.d0
		q = 0.08918d0
		J(2) = 14697.d-6
		J(4) =  -584.d-6 
		J(6) =    31.d-6
	endif
c
c Saturne (Cf. Table de Bill), valeurs a 1 bar
c
	if (iplanete.eq.2) then
		GM= (G*MT*95.147d0)/1.d9	! en km^3 s^-2 !
		a = 60268.d0
		q = 0.154766d0
		J(2) = 16331.d-6
		J(4) =  -914.d-6 
		J(6) =   108.d-6
	endif
c
c Uranus (Cf. Table de Bill)
c
	if (iplanete.eq.3) then
		GM= (G*MT*14.53d0)/1.d9		! en km^3 s^-2 !
		a = 25559.d0
		q = 0.02951d0
		J(2) = 3516.d-6
		J(4) = -31.9d-6 
		J(6) =    0.d-6		! ???
	endif
c
c Neptune (Cf. Table de Bill)
c
	if (iplanete.eq.4) then
		GM= (G*MT*17.14d0)/1.d9		! en km^3 s^-2 !
		a = 24764.d0
		q = 0.0261d0
		J(2) = 3538.d-6
		J(4) =  -38.d-6 
		J(6) =    0.d-6		! ???
	endif
c
c Pluton (Tholen et al AJ 2008)
c
	if (iplanete.eq.5) then
		GM= 870.3d0		! en km^3 s^-2 !
		a = 1195.d0
		q =    0.d0
		J(2) = 0.d0
		J(4) = 0.d0
		J(6) = 0.d0
	endif
c
c Titan
c
	if (iplanete.eq.6) then
		GM= 8979.6d0		! en km^3 s^-2 !
		a = 2575.d0
		q =    0.d0
		J(2) = 0.d0
		J(4) = 0.d0
		J(6) = 0.d0
	endif
c
c Charon (Bob Jacobson, email mai 2005, et Sicardy Nature 2005)
c
	if (iplanete.eq.7) then
		GM= (G*1.58d21)/1.d9		! en km^3 s^-2 !
		a = 603.d0
		q =    0.d0
		J(2) = 0.d0
		J(4) = 0.d0
		J(6) = 0.d0
	endif
c
c Titania  (de http://ssd.jpl.nasa.gov/?sat_phys_par)
c
	if (iplanete.eq.8) then
		GM= 235.3d0		! en km^3 s^-2 !
		a = 789.d0
		q =    0.d0
		J(2) = 0.d0
		J(4) = 0.d0
		J(6) = 0.d0
	endif

c
c Ceres  (de Carry et al. 2008)
c
	if (iplanete.eq.9) then
		GM= 62.9d0		! en km^3 s^-2 !
		a = 490.d0
		q =    0.d0
		J(2) = 0.d0
		J(4) = 0.d0
		J(6) = 0.d0
	endif

c
c Eris ???
c
	if (iplanete.eq.10) then
		GM= 1114.1d0		! en km^3 s^-2 !
		a = 1111.d0
		q =    0.d0
		J(2) = 0.d0
		J(4) = 0.d0
		J(6) = 0.d0
	endif
c
c Charon (Bob Jacobson, email mai 2005, et Sicardy Nature 2005)
c
	if (iplanete.eq.11) then
		GM= (G*1.58d21)/1.d9		! en km^3 s^-2 !
		a = 603.d0
		q =    0.d0
		J(2) = 0.d0
		J(4) = 0.d0
		J(6) = 0.d0
	endif
c
c Chariklo ???
c
	if (iplanete.eq.12) then
		GM=   1.2d0		! en km^3 s^-2 !
		a =  110.d0
		q =    0.d0
		J(2) = 0.d0
		J(4) = 0.d0
		J(6) = 0.d0
	endif

	return
	end
c---------------------------- Fin de parametres -----------------------------

c Sous-programme qui calcule par iteration le rayon de l'isopotentielle
c de rayon equatorial req, pour la latitude lat, et avec les parametres
c planetaires definis dans le common
c
	subroutine limbe (req, lat,r)
	implicit real*8 (a-h,o-z)
	real*8 lat, J(6)
	common /grav/GM,a,J,q
	pi= dacos(-1.d0)
	rad= pi/180.d0
	itermax= 20
	iter=     0

	theta= 90.d0*rad
	U= -(GM/req)*F(req,theta)	! valeur du potentiel a l'equateur

c calcul de r pour la colatitude theta 

	theta=  pi/2.d0 - lat
	r_0=   -(GM/U)			! estimation initiale du rayon 

	tol= 0.001d0	
 99	r= -(GM*F(r_0,theta)/U) 	! correction
	if (dabs(r-r_0).ge.tol) then
		r_0= r
		iter= iter + 1
		if (iter.ge.itermax) then
		   write (*,*) 'Nombre d''iterations maximales depasse!', iter
		   return
		endif
		go to 99
	endif

	return
	end

c---------------------------- Fin de limbe -----------------------------

c----------------------------------------------------------------------------
	function F(r,theta)
        implicit real*8 (a-h,o-z)
	real*8 mu, J(6)
	common /grav/GM,a,J,q

	mu=   dcos(theta)	
	sin2= 1.d0 - mu**2
	som= 0.d0

	do i= 2, 6, 2
		som = som + J(i)*( (a/r)**i )*Pn(i,mu)
	enddo

	F = 1.d0 - som + (q/2.d0)*( (r/a)**3 )*sin2

	return
	end
c-------------------------- fin de F ----------------------------------------

c----------------------------------------------------------------------------
	function Pn(n,mu)
        implicit real*8 (a-h,o-z)
        real*8 mu
c
c Polynomes de Legendre d'ordres n= 2, 4, 6
c
c Ici, on calcule explicitement les polynomes. On pourrait aussi utiliser la
c relation de recurrence:
c
c (n+1)*P_(n+1)[mu] = (2n+1)*mu*P_n[mu] - nP_(n-1)[mu]
c
c Cf. "Analyse Numerique" , Francis Scheid, p. 134
c
	if (n.eq.2) then
		Pn= ( 3.d0*(mu**2) - 1.d0 )/2.d0
	endif

	if (n.eq.4) then
		Pn= ( 35.d0*(mu**4) - 30.d0*(mu**2) + 3.d0 )/8.d0
	endif

	if (n.eq.6) then
		Pn= ( 231.d0*(mu**6) - 315.d0*(mu**4) + 
     *  	      105.d0*(mu**2) - 5.d0 )/16.d0
	endif

	return
	end
c-------------------------- fin de Pn ---------------------------------------

	subroutine correc (lat, d_r)

	dimension phi(2500), dr(2500)
	real*4 lat
	common/corr/phi,dr
	rad= acos(-1.)/180.
	latd= lat/rad		! car Phi est en degres
c
c Corrections a apporter au limbe dues aux vents zonaux
c ND. dr en km !!!!
c
	do i= 1, 2000
	fac= (latd-phi(i))*(latd-phi(i+1))
	 if (fac.le.0.) then
		coeff= (dr(i+1)-dr(i))/(phi(i+1)-phi(i))
		d_r = dr(i) + coeff*(latd-phi(i))
		return
	 endif
	enddo

	return
	end
c ------------------------- fin de correc -----------------------------------
