c
c Calcule le decalage en temps a appliquer a chaque corde afin de l'aligner
c sur un cercle de centre (xc, yc) donne. NB. le rayon du cercle n'a pas besoin
c d'etre precis. Il est utlise uniquement en sortie pour calculer l'ecart de
c chaque point au cercle de rayon "Rayon" et de centre (xc, yc)
c
c On ressort en ecrivant la nouvelle corde decalee.
c
c NB. fort.99: utile pour ellipse_fit.i (erreur_radiale, station)
c

	real*4 long, mhx, mhy
c	integer*1 iobsa, iobsb ?????????????????????
	character*80 fichier
	character*50 nom(1000)
 1000	format (A)

	call liste(nom)
c
c (Titania: Rayon= 787.7)
c
	write(*,*) 'Rayon de reference?'
	read(*,*) Rayon 

1	write (*,*) 'Fichier a lire ?'
	read 1000, fichier
	open (unit=20,file=fichier,status='old',form='formatted',err=1)

2	write (*,*) 'Fichier de sortie (texte)?'
	read 1000, fichier
	open (unit=21,file=fichier,status='unknown',form='formatted',err=2)
	write(21,*) ' '
	write(21,'(A,1x,f8.3)') ' Rayon de reference:', Rayon
	write(21,*) ' '
	write(21,'(A)') 'Code   Dt     r-r_ref   vrA     vrB      C/A      site'
	write(21,*) ' '

3	write (*,*) 'Fichier de sortie (nouvelle cordes)?'
	read 1000, fichier
	open (unit=22,file=fichier,status='unknown',form='formatted',err=3)

	write (*,*)  'Position supposee du centre du cercle?'
	read (*,*)  xc, yc

	nmax= 1000
	k= 0
	seuil= 1.d-3
	
	do i= 1, nmax
c
c A: 1er point de la corde, B: 2eme point de la corde
c
!		read (20,*,end=99,err=99) iobsa, heurea, xa, ya 
!		read (20,*,end=99,err=99) iobsb, heureb, xb, yb 
		read (20,*,end=99,err=99) xa, ya, heurea, dumb, iobsa
		read (20,*,end=99,err=99) xb, yb, heureb, dumb, iobsb
c
c verifie que les deux points appartiennent bien a la meme corde
c
		if (iobsa.ne.iobsb) write (*,*)  'Erreur, pas le meme iobs'
		long= sqrt( (xa-xb)**2 + (ya-yb)**2 )
		Delta_t= heureb - heurea	! heures en secondes
		vx= (xb-xa)/Delta_t
		vy= (yb-ya)/Delta_t
		V= long/Delta_t				! km/sec
		ux= (xb-xa)/long			! vec(AB)/AB
		uy= (yb-ya)/long
		acu= (xc-xa)*ux + (yc-ya)*uy	! vec(AC).vec(u)
		xh= xa + acu*ux					! H: projection de C sur AB
		yh= ya + acu*uy
		CH= sqrt( (xh-xc)**2 + (yh-yc)**2 )	! plus courte distance du centre a la corde AB
		if ((yh-yc).lt.(0.0)) CH= -CH			! si la corde passe en-dessous de C
		xm= (xa+xb)/2.				! M: milieu de AB
		ym= (ya+yb)/2.
		mhx= xh-xm					! MH
		mhy= yh-ym
		Dl = mhx*ux + mhy*uy		! decalage en distance 
		dt = Dl/V					! decalage en temps
		xa_decal= xa + Dl*ux		! point A decale
		ya_decal= ya + Dl*uy
		ra= sqrt((xa_decal - xc)**2 + (ya_decal - yc)**2)
		xb_decal= xb + Dl*ux		! point B decale
		yb_decal= yb + Dl*uy
		rb= sqrt((xb_decal - xc)**2 + (yb_decal - yc)**2)
		if (abs(rb-ra).gt.seuil) then
		 write (*,*)  'Les deux rayons sont differents!'
		endif
		vra= ((xa_decal - xc)*vx + (ya_decal - yc)*vy)/ra			! vitesse radiale en A
		vrb= ((xb_decal - xc)*vx + (yb_decal - yc)*vy)/rb			! vitesse radiale en B
		write (21,'(I3,1x,f8.3,2x,f8.3,2x,f7.3,1x,f7.3,2x,f7.1,4x,A)') iobsa, dt, ra - Rayon, vra, vrb, CH, nom(iobsa)
		write (22,*) xa_decal, ya_decal, nom(iobsa)
		write (22,*) xb_decal, yb_decal
		write (99,'(A,A)') '1.3     ', nom(iobsa)
		write (99,'(A,A)') '1.3     ', nom(iobsb)
		k= k+1
	enddo
 99	continue

	write (*,*)  k, ' cordes examinees'
	write (21,*) ' '
	write (21,'(I2,1x,A)') k, 'cordes examinees'

	stop
	end
********************************** fin de main ***************************

	subroutine liste(nom)
	character*50 nom(1000)

	nom(51)= 'Salinas'
	nom(52)= 'Cuenca'
	nom(53)= 'Maracaibo'
	nom(54)= 'Bobares '
	nom(55)= 'Caracas'
	nom(56)= 'Arikok, Aruba'
	nom(57)= 'Fort de France'
	nom(58)= 'Tobago'
	nom(59)= 'Trinidad'
	nom(60)= 'Barbados'
	nom(61)= 'Acores'
	nom(62)= 'Oeiras'
	nom(63)= 'Portimao'
	nom(64)= 'COAA Algarve (Portimao)'
	nom(65)= 'Linhaceira'
	nom(66)= 'Alvito'
	nom(67)= 'Granada'
	nom(68)= 'Alcublas'
	nom(69)= 'Bordeaux'
	nom(70)= 'Pic du Midi'
	nom(71)= 'St Maurice (Poitou)'
	nom(72)= 'Pezenas'
	nom(73)= 'Mauguio'
	nom(74)= 'Nimes'
	nom(75)= 'Orfeuilles'
	nom(76)= 'OHP'
	nom(77)= 'Binfield'
	nom(78)= 'Worth Hill, Bournemouth'
	nom(79)= 'Teversham'
	nom(80)= 'Plateau d''Albion '
	nom(81)= 'Bd H. Fabre, Marseille'
	nom(82)= 'Guitalens'
	nom(83)= 'Sabadell (Ardanuy)'
	nom(84)= 'Salon'
	nom(85)= 'Wela, Aruba'
	nom(86)= 'St Maurice Cazevieille'
	nom(87)= 'Marinha Grande'
	nom(88)= 'Almeirim'
	nom(89)= 'Carcavelos'
	nom(90)= 'Setubal'
	nom(91)= 'Alcacer do Sal'
	nom(92)= 'Barcelona 1'
	nom(93)= 'St Savinien'
	nom(94)= 'St Martin de Crau'
	nom(95)= 'Barcelona 2'
	nom(96)= 'St Esteve'
	nom(97)= 'Alella'
	nom(98)= 'Castellon'
	nom(99)= 'Hortoneda'
	nom(100)= 'Barcelona 3'
	nom(101)= 'Sabadell (Casas)'
	nom(102)= 'Esplugues de Llobregat'
	nom(103)= 'Zaragoza'
	nom(104)= 'Calern'
	nom(105)= 'Dax'
	nom(106)= 'Puimichel'
	nom(107)= 'Chatellerault'

	nom(110)= 'George Observatory'
	nom(111)= 'Monterrey'
	
	nom(160)= 'Tomar'
	nom(161)= 'Sabadell (Casas)'
	nom(162)= 'Pic du Midi'
	nom(163)= 'Czarna Bialostocka'
	nom(164)= 'Vitebsk'
	nom(165)= 'Meudon'
	nom(166)= 'Kooriyama'
	nom(167)= 'Ooe'
	nom(168)= 'Kashiwa'
	nom(169)= 'Hitachi'
	nom(170)= 'Chichibu'
	nom(171)= 'Musashino'
	nom(172)= 'Mitaka'

	nom(220)= 'Windhoek'
	nom(221)= 'Tivoli'
	nom(222)= 'HESS'
	nom(223)= 'Hakos'
	nom(224)= 'Kleinbegin'
	nom(225)= 'Sandfontein'
	nom(226)= 'Springbok'
	nom(227)= 'Nuwerus'
	nom(228)= 'Gifberg'
	nom(229)= 'Cederberg'
	nom(230)= 'SAAO Sutherland'
	nom(231)= 'SAAO Cape Town'
	nom(232)= 'Boyden'
	nom(240)= 'Maido'
	nom(241)= 'Les Makes'
	nom(242)= 'Fournaise'

	nom(630)= 'Charon23jun11 Kauai Kekaha'
	nom(631)= 'Charon23jun11 Maui Lahaina'
	nom(632)= 'Charon23jun11 CFHT'
	nom(633)= 'Charon23jun11 Haleakala FTN'
	nom(634)= 'Charon23jun11 Hale-a'
	nom(635)= 'Charon23jun11 Majuro'
	nom(636)= 'Charon23jun11 S Pedro Martir 84'
	nom(637)= 'Charon23jun11 S Pedro Martir 210'
	nom(638)= 'Charon23jun11 Kauai KEASA'
	nom(639)= 'Charon23jun11 xxx'

	nom(700)= '2003VS2 07nov2014 - BA'
	nom(701)= '2003VS2 07nov2014 - BU'
	nom(702)= '2003VS2 07nov2014 - LM'
	nom(703)= '2003VS2 07nov2014 - SR'





	
	nom(900)= '2007Uk126 14nov2014 - Carson (B)'
	nom(901)= '2007Uk126 14nov2014 - Carson (S)'
	nom(902)= '2007Uk126 14nov2014 - Gardnerville'
	nom(903)= '2007Uk126 14nov2014 - Tonopah'
	nom(904)= '2007Uk126 14nov2014 - Yerington'
	nom(905)= '2007Uk126 14nov2014 - IOTA'
	nom(908)= '2007Uk126 14nov2014 - Reno'
	
	
	return
	end
*********************************** fin de liste **************************
