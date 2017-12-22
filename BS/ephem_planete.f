	
c Ephemeride de Pluton pour le 01 juillet 2002 (coordonnees ICRS???, JPL Horizon)
c Option: ephemeride de Charon 01 juillet 2002
c                              20 juillet 2002
c                       Pluton 20 juillet 2002
c                       Pluton 21 aout    2002
c etc...
c 
c fort.20: t, ksi, Delta_km
c fort.21: t, eta
c fort.22: ksi, eta, t, sqrt( ksi**2 + eta**2)
c fort.23: bande planete
c fort.24: bande atmosphere
c fort.25: Draa*cos(delta), Ddec (arcsec), t, de la planete /t a l'etoile
c fort.26: ksi, -eta, t pour etre en accord avec positionv
c fort.27: hhmmsec, Draa*cos(delta), Ddec (arcsec) pour faire une "tracking table" NACO
c (NB offset **geocentrique** et non pas topocentrique, mais difference negligeable en general,
c sinon, voir positionv
c fort.28: coordonnees de l'etoile
c
	implicit real*8 (a-h,o-z)
	dimension t(1000),  at(1000),  dt(1000),  delta(1000)
	dimension xi(1000), eta(1000)
	pi=   dacos(-1.d+0)
	radd= pi/180.d+00
	radh= pi/12.d+00
	arcsec= pi/(180.d0*3600.d0)
	ua= 1.4959787061d+08
	R_Plu=      1223.0d0		! rayon de l'ombre a miocc d'apres 2002 mais distance 18 mars 07 (cahier 4 avril 07)
	R_Cha=       604.0d0
	R_Hyd=        50.0d0
	R_Titania=   788.4d0
	R_atmo=     1557.0d0		! rayon de l'ombre a 1% de chute d'apres 2002 mais distance 18 mars 07 (cahier 4 avril 07)
	R_Tri=      1353.0d0		! rayon Triton JPL
	R_atmo_Tri= 1500.0d0		! cf. Olkin et al. annees 96 ou 97
c
c NB. les distances ci-dessous sont changeables plus bas pour 2007 
c
c	write(*,*)  '1=  Pluton, 01 juillet 2002'
c	write(*,*)  '2=  Charon, 01 juillet 2002'
c	write(*,*)  '3=  Pluton, 20 juillet 2002'
c	write(*,*)  '4=  Charon, 20 juillet 2002'
c	write(*,*)  '5=  Pluton, 09 aout    2002'
c	write(*,*)  '6=  Charon, 09 aout    2002'
c	write(*,*)  '7=  Pluton, 21 aout    2002'
c 	write(*,*)  '8=  Pluton, 07 novembr 2002'

c 	write(*,*)  '9=  Titan,  14 novembr 2003, TYC 1343-1615, JPL'
c 	write(*,*)  '10= Titan,  14 novembr 2003, TYC 1343-1615, BdL'
c 	write(*,*)  '11= Titan,  14 novembr 2003, TYC 1343-1615, BS Tass7'
c 	write(*,*)  '12= Titan,  14 novembr 2003, TYC 1343-1615, BS Tass8'
c 	write(*,*)  '13= Titan,  14 novembr 2003, TYC 1343-1615, JAC'

c 	write(*,*)  '14= Titan,  14 novembr 2003, TYC 1343-1865, JPL'
c 	write(*,*)  '15= Titan,  14 novembr 2003, TYC 1343-1865, BdL'
c 	write(*,*)  '16= Titan,  14 novembr 2003, TYC 1343-1865, BS Tass7'
c 	write(*,*)  '17= Titan,  14 novembr 2003, TYC 1343-1865, BS Tass8'
c 	write(*,*)  '18= Titan,  14 novembr 2003, TYC 1343-1865, JAC'

c 	write(*,*)  '20= Triton, 29 novembr 2003' 
c 	write(*,*)  '21= Tethys, 15 decembr 2002' 
c 	write(*,*)  '22= Jupite, 01 avril   2003' 
c 	write(*,*)  '23= Jupite, 10 Octobre 1999' 

c 	write(*,*)  '24= Titan,  24 aout    2004, HIP 37084, BdL, DE405+TASS'
c 	write(*,*)  '25= Titan,  24 aout    2004, HIP 37084, BdL, VSOP87+TASS'
c 	write(*,*)  '26= Titan,  24 aout    2004, BdL, DE405+TASS, 5sec'

c	write(*,*)  '27= Pluton, 22 mai 2005, JPL DE-0406LE-0406, VLT'
c	write(*,*)  '28= Pluton, 22 mai 2005, JPL DE-0406LE-0406, GEO'
c	write(*,*)  '29= Pluton, 22 mai 2005, JPL DE-0406LE-0406, CFH'

c	write(*,*)  '30= Charon, 11 juillet 2005, JPL DE-0406LE-0406'
c	write(*,*)  '31= Charon, 11 juillet 2005, JPL DE-0413+Stone+JLX'
c	write(*,*)  '32= Charon, 11 juillet 2005, JPL DE-0413+LNA'
	
c	write(*,*)  '42= Charon, 20 septemb 2005, JPL DE-0406LE-0406'
c	write(*,*)  '43= Charon, 26 juillet 2005, JPL DE413'
	
c	write(*,*)  '44= Pluton, 10 avril 2006, JPL DE413 + R. Stone 2004'
c	write(*,*)  '45= Pluton, 10 avril 2006, JPL DE413 + R. Behrend 2006'
c	write(*,*)  '46= Charon, 10 avril 2006, JPL DE413 + R. Behrend 2006'
c	write(*,*)  '47= Pluton, 10 avril 2006, etoile voisine,R. Behrend 2006'

c	write(*,*)  '48= P1,     06 aout  2006, DE413 + R. Behrend 2006'

c	write(*,*)  '49= Pluton, 12 juin  2006, DE413 + R. Stone   2004'
c	write(*,*)  '50= Pluton, 12 juin  2006, DE413 + R. Behrend 2006'
c	write(*,*)  '51= P1,     12 juin  2006, DE413 + R. Behrend 2006'

c	write(*,*)  '60= Pluton, 18 mars  2007, DE413 + R. Behrend ICRF 04 mars 2007'
c	write(*,*)  '61= Pluton, 18 mars  2007, DE413 + M. Assafim ICRF 05 mars 2007'
		
c	write(*,*)  '70= Pluton, 12 mai   2007, DE413 + R. Behrend A ICRF'
c	write(*,*)  '71= Pluton, 12 mai   2007, DE413 + M. Assafim A ICRF'
c	write(*,*)  '72= Pluton, 12 mai   2007, DE413 + R. Behrend B ICRF'
c	write(*,*)  '73= Pluton, 12 mai   2007, DE413 + M. Assafim B ICRF'
		
c	write(*,*)  '80= Pluton, 14 juin  2007, + R. Behrend ICRF'
c	write(*,*)  '81= Pluton, 14 juin  2007, + M. Assafim ICRF'

c	write(*,*)  '82= Pluton, 31 juillet  2007, + R. Behrend ICRF'
c	write(*,*)  '83= Pluton, 31 juillet  2007, + M. Assafim ICRF'

c	write(*,*)  '84= Pluton, 09 juin  2007, + R. Behrend ICRF, predic. 08 sept 07'
c	write(*,*)  '85= Pluton, 14 juin  2007, + R. Behrend ICRF, predic. 16 avr 07'

c	write(*,*)  '86= Titania, 01 aout 2003'

c	write(*,*)  '87= Triton, 21 mai 2008, predic Rio 27 avril 2008'

c	write(*,*)  '88= Pluton, 22 juin 2008, predic Rio 02 juin 2008'
c	write(*,*)  '89= Charon, 22 juin 2008, predic Rio 02 juin 2008'

c	write(*,*)  '90= Pluton, 24 juin 2008, predic Rio 09 juin 2008'

c	write(*,*)  '91= 2002 MS4, 30 juin 2009, predic Rio 25 juin 2009'

c	write(*,*)  '92= Pluton, 25 aout 2008, predic Rio 23 juillet 2008'

c	write(*,*)  '93= Varuna, 19 fevrier 2010, predic Rio 04 fevrier 2010'

c	write(*,*)  '94= Pluton, 14 fevrier 2010, predic Rio 10 nov. 2009 '

c	write(*,*)  '95= Pluton, 19 mai 2010, predic Rio 10 nov. 2009'

c	write(*,*)  '96= Pluton, 04 juin 2010, predic Rio 10 nov. 2009'

c	write(*,*)  '97= Pluton, 04 juillet 2010, predic Rio 02 juillet 2010'

c	write(*,*)  '98= Ceres, 17 aout 2010, predic UCAC2'

c	write(*,*)  '99=  Eris, 06 novembre 2010 (JPL28), star Bill Owen'
c	write(*,*)  '100= Eris, 06 novembre 2010 (JPL29), star Bill Owen'
c	write(*,*)  '101= Eris, 05 novembre 2010 (JPL29 OPD), star Bill Owen'
c	write(*,*)  '102= Eris, 05 novembre 2010 (JPL29 geocentrique), star Bill Owen'

c	write(*,*)  '110= Makemake, 23 avril 2011 (JPL43 geocentrique), star JL Ortiz DDT '
c	write(*,*)  '110= Makemake, 23 avril 2011 (JPL46 geocentrique), star JL Ortiz DDT '

c	write(*,*)  '111= Quaoar, 04 mai 2011 (JPL17 geocentrique), star JL Ortiz 03 May 2011'
	
c	write(*,*)  '112= Pluton, 04 juin 2011, predic Rio'
c	write(*,*)  '113= Charon, 04 juin 2011, predic Rio'
c	write(*,*)  '114= Pluton, 04 juin 2011, predic Rio, ephem MB/DE413'
c	write(*,*)  '115= Charon, 04 juin 2011, predic Rio, ephem MB/DE413'

c	write(*,*)  '120= Pluton, 23 juin 2011, ephem DE413 + offset MB, predic Rio'
c	write(*,*)  '121= Charon, 23 juin 2011, ephem DE413 + offset MB, predic Rio'

c	write(*,*)  '122= Pluton, 27 juin 2011, ephem DE413 + offset MB, predic Rio'
c	write(*,*)  '123= Hydra,  27 juin 2011, ephem DE413 + offset MB, predic Rio'

c	write(*,*)  '124= Pluton, 23 juin 2011, ephem DE413 + PLU021, predic Rio'
c	write(*,*)  '125= Charon, 23 juin 2011, ephem DE413 + PLU021, predic Rio'

c	write(*,*)  '126= Pluton, 14 juin 2012, ephem DE413 + PLU021, predic Rio'

c	write(*,*)  '127= Pluton, 18 juil 2012, ephem DE413 + PLU022, predic Rio'

c	write(*,*)  '128= Varuna, 05 janv 2013, ephem JPL25, etoile Raoul'
c	write(*,*)  '129= Varuna, 06 janv 2013, ephem JPL25, etoile Raoul'
c	write(*,*)  '130= Varuna, 08 janv 2013, ephem JPL25, etoile Raoul'

c	write(*,*)  '140= 2002 KX14, 26 avr 2012, ephem JPL9, etoile papier Alvaro'

c	write(*,*)  '141= Pluton, 04 mai 2013, ephem DE413 + PLU031, etoile IAG Rio (mail FB 06mai13)'

c	write(*,*)  '142= Chariklo, 03 juin 2013, ephem JPL20, etoile OPD Rio (carte FB mai 2013)'

c	write(*,*)  '143= Chariklo, 29 avril 2014, ephem JPL21, etoile OPD Rio (carte FB 17 avril 2014)'

c	write(*,*)  '144= Chariklo, 16 fevrier 2014, ephem JPL21, etoile OPD 08sep13_T60'

c	write(*,*)  '145= Chariklo, 28 juin 2014, ephem JPL22, etoile OPD 30may14_T160'
	
c	write(*,*)  '146= 2007 UK126, 2014 nov 15, ephem JPL7, star OPD'


c	write(*,*)  '147= 2003 AZ84, 2011 jan 08'
c	write(*,*)  '148= 2003 AZ84, 2012 fev 03'
c	write(*,*)  '149= 2003 AZ84, 2014 nov 15'
	write(*,*)  '150= 2003 VS2, 2014 nov 07'

	read(*,*)    iop


	if (iop.eq.1) R_pla= R_Plu
	if (iop.eq.2) then
		R_pla= R_Cha
		R_atmo=R_Cha
	endif
	if (iop.eq.3) R_pla= R_Plu
	if (iop.eq.4) then
		R_pla= R_Cha
		R_atmo=R_Cha
	endif
	if (iop.eq.5) R_pla= R_Plu
	if (iop.eq.6) then
		R_pla= R_Cha
		R_atmo=R_Cha
	endif
        if (iop.eq.7) R_pla= R_Plu
	if (iop.eq.8) R_pla= R_Plu
	if ((iop.ge.9).and.(iop.le.18)) then 
		R_pla=  3100.d0		! en fait, rayon de L'ATMOSPHRE
		R_atmo=  170.d0		! rayon du FLASH CENTRAL
	endif
	if (iop.eq.20) then 
		R_pla=  1352.
		R_atmo= 1552.
	endif
	if (iop.eq.21) then 
		R_pla=  535. 
		R_atmo= 535.
	endif
	if (iop.eq.22) then 
		R_pla=  71700. 
		R_atmo= 71700.
	endif
c
	if (iop.eq.23) then 
		R_pla=  67164.		! polar half-light level
		R_atmo= 67164.
	endif

	if ((iop.ge.24).and.(iop.le.26)) then 
		R_pla=  3100.d0		! en fait, rayon de L'ATMOSPHRE
		R_atmo=  150.d0		! rayon du FLASH CENTRAL
	endif

	if ((iop.ge.27).and.(iop.le.29)) then
			R_pla=  R_Plu 
			R_atmo= R_Plu + 200.
	endif

	if ((iop.ge.30).and.(iop.le.31)) then
		R_pla= R_Cha
		R_atmo=R_Cha + 785.d0	! 785 km= 36 mas= sqrt(20^2 + 30^2), 
					! erreur possible sur ephem.+etoile???
	endif
	
	if (iop.eq.32) then
		R_pla= R_Cha
c		R_atmo=R_Cha + 325.d0	! 325 km= 15 mas= sqrt(10^2 + 11^2), 
					! erreur possible sur ephem.+etoile???
		R_atmo=R_Cha + 0.d0	! post-occultation !
	endif

	if ((iop.ge.42).and.(iop.le.43)) then
		R_pla= R_Cha
		R_atmo=R_Cha 		! a preciser! 
	endif	
	
	if ((iop.ge.44).and.(iop.le.44)) then
			R_pla=  R_Plu 
			R_atmo= R_Plu + 2200.	! erreur typique de 100 mas
	endif

	if ((iop.ge.45).and.(iop.le.47)) then
			R_pla=  R_Plu 
			R_atmo= R_Plu + 220.	! erreur typique de 10 mas
	endif

	if ((iop.ge.48).and.(iop.le.48)) then
			R_pla=   50.		! rayon P1 25-75 km		 
			R_atmo= R_Plu + 500.	! incertitude (?)
	endif

	if ((iop.ge.49).and.(iop.le.50)) then
			R_pla=  R_Plu 
			R_atmo= R_Plu + 250.	! epaisseur typique atmos. detectable (5% chute)
	endif

	if ((iop.ge.51).and.(iop.le.51)) then
			R_pla=  100.		! rayon P1 ?
			R_atmo= 100. + 500.	! incertitude
	endif

	if ((iop.ge.60).and.(iop.le.85)) then
			R_pla=  R_Plu 
			R_atmo= 1557.	!  1% chute si = 2002, mais avec distance 18 mars 07 (cahier 4/4/07) 
c			R_atmo= 1507.	!  2% chute si = 2002, mais avec distance 18 mars 07 (cahier 4/4/07)
c			R_atmo= 1443.	!  5% chute si = 2002, mais avec distance 18 mars 07 (cahier 4/4/07)			
c			R_atmo= 1394.	! 10% chute si = 2002, mais avec distance 18 mars 07 (cahier 4/4/07)
c			R_atmo= 1223.	! 50% chute si = 2002, mais avec distance 18 mars 07 (cahier 4/4/07)
c			R_atmo= 1102.	! 70% chute si = 2002, mais avec distance 18 mars 07 (cahier 4/4/07)
c			R_atmo=  953.	! 90% chute si = 2002, mais avec distance 18 mars 07 (cahier 4/4/07)						
	endif

	if (iop.eq.86) then
			R_pla=  R_Titania
			R_atmo= R_Titania
	endif

	if (iop.eq.87) then
			R_pla=  R_Tri
			R_atmo= R_atmo_Tri
	endif
c
c Pluton et Charon 22 juin 2008:
c
	if (iop.eq.88) then
			R_pla=  R_Plu
			R_atmo= 1557.	!  1% chute si = 2002, mais avec distance 18 mars 07 (cahier 4/4/07) 
	endif
	if (iop.eq.89) then
			R_pla=  R_Cha
			R_atmo= R_Cha 
	endif
c
c Pluton 24 juin 2008:
c 
	if (iop.eq.90) then
			R_pla=  R_Plu
			R_atmo= 1557.	!  1% chute si = 2002, mais avec distance 18 mars 07 (cahier 4/4/07) 
	endif
c
c 2002 MS4, 30 juin 2009:
c 
	if (iop.eq.91) then
			R_pla=  350.
			R_atmo= 350.	!  Stansberry et al, UoA 
	endif
c
c Pluton 25 aout 2008:
c 
	if (iop.eq.92) then
			R_pla=  R_Plu
			R_atmo= 1557.	!  1% chute si = 2002, mais avec distance 18 mars 07 (cahier 4/4/07) 
	endif
c
c Varuna 19 fevrier 2010:
c 
	if (iop.eq.93) then
			R_pla=   500.
			R_atmo=  500.	! approx. 
	endif
c
c Pluto 14 fevrier 2010, 19 mai 2010, 04 juin 2010, 04 juillet 2010:
c 
c 
	if (iop.ge.94.and.iop.le.97) then
			R_pla=   R_Plu
			R_atmo=  1557.	! ~ 1% chute  
	endif

c
c Ceres 17 aout 2010
c 
c 
	if (iop.ge.98.and.iop.le.98) then
			R_pla=    480.
			R_atmo=   480.
	endif
c
c Eris 06 novembre 2010
c 
c 
	if (iop.ge.99.and.iop.le.102) then
			R_pla=   1200.
			R_atmo=  1200.
	endif
c
c Makemake 23 avril 2011
c 
c 
	if (iop.ge.110.and.iop.le.110) then
			R_pla=    750.
			R_atmo=   750.
	endif
c
c Quaoar 04 mai 2011
c 
c 
	if (iop.ge.111.and.iop.le.111) then
			R_pla=    630.
			R_atmo=   630.
	endif
c
c Pluto 04 juin 2011:
c 
c 
	if (iop.ge.112.and.iop.le.112) then
			R_pla=   R_Plu
			R_atmo=  1557.	! ~ 1% chute  
	endif
c
c Charon 04 juin 2011:
c 
c 
	if (iop.ge.113.and.iop.le.113) then
			R_pla=   R_Cha
			R_atmo=  R_Cha
	endif
c
c Pluto 04 juin 2011:
c 
c 
	if (iop.ge.114.and.iop.le.114) then
			R_pla=   R_Plu
			R_atmo=  1557.	! ~ 1% chute  
	endif
c
c Charon 04 juin 2011:
c 
c 
	if (iop.ge.115.and.iop.le.115) then
			R_pla=   R_Cha
			R_atmo=  R_Cha
	endif
c
c Pluto 23 juin 2011 pour DE413 + offset MB:
c 
c 
	if (iop.ge.120.and.iop.le.120) then
			R_pla=   R_Plu
			R_atmo=  1557.	! ~ 1% chute  
	endif
c
c Charon 23 juin 2011 pour DE413 + offset MB:
c 
c 
	if (iop.ge.121.and.iop.le.121) then
			R_pla=   R_Cha
			R_atmo=  R_Cha
	endif
c
c Pluto 27 juin 2011:
c 
c 
	if (iop.ge.122.and.iop.le.122) then
			R_pla=   R_Plu
			R_atmo=  1557.	! ~ 1% chute  
	endif
c
c Hydra 27 juin 2011:
c 
c 
	if (iop.ge.123.and.iop.le.123) then
			R_pla=   R_Hyd
			R_atmo=  R_Hyd
	endif
c
c Pluto 23 juin 2011 pour DE413 + PLU021:
c 
c 
	if (iop.ge.124.and.iop.le.124) then
			R_pla=   R_Plu
			R_atmo=  1557.	! ~ 1% chute  
	endif
c
c Charon 23 juin 2011 pour DE413 + PLU021:
c 
c 
	if (iop.ge.125.and.iop.le.125) then
			R_pla=   R_Cha
			R_atmo=  R_Cha
	endif
c
c Pluto 14 juin 2012 pour DE413 + PLU021 et 18 juillet 2012 DE413 + PLU022
c 
	if (iop.ge.126.and.iop.le.127) then
		R_pla=   R_Plu
		R_atmo=  1557.	! ~ 1% chute  
	endif

c
c Varuna 08 janvier 2013:
c 
	if (iop.ge.128.and.iop.le.130) then
			R_pla=   500.
			R_atmo=  500.	! approx. 
	endif

c
c 2002 KX14 26 avril 2012: rayon equivalent papier Alvaro
c 
	if (iop.ge.140.and.iop.le.140) then
			R_pla=   516.
			R_atmo=  516.	! approx. 
	endif
c
c Pluto 04 mai 2013:
c 
	if (iop.ge.141.and.iop.le.141) then
			R_pla=   R_Plu
			R_atmo=  1557.	! ~ 1% chute  
	endif
c
c Chariklo 03 juin 2013:
c 
	if (iop.ge.142.and.iop.le.142) then
			R_pla=   130.		! approx. !!!
			R_atmo=  130. 
	endif
c
c Chariklo 29 avril 2014:
c 
	if (iop.ge.143.and.iop.le.143) then
			R_pla=   124.		! approx. !!!
			R_atmo=  400.		! rings
	endif
c
c Chariklo 16 fevrier 2014:
c 
	if (iop.ge.144.and.iop.le.144) then
			R_pla=   124.		! approx. !!!
			R_atmo=  400.		! rings
	endif
c
c Chariklo 28 juin 2014:
c 
	if (iop.ge.145.and.iop.le.145) then
			R_pla=   124.		! approx. !!!
			R_atmo=  400.		! rings
	endif


c
c 2007UK126 15 nov 2014:
c 
	if (iop.ge.146.and.iop.le.146) then
			R_pla=   306.		! Value from Michael E. Brown = 612 km (diam.)
			R_atmo=  306.		! Value from P. Santos-Sanz = 599 +- 77 km (diam.)
	endif


c
c 2003AZ84 08 jan 2011:
c 
	if (iop.ge.147.and.iop.le.147) then
			R_pla=   345.		! Value from Michael E. Brown = 612 km (diam.)
			R_atmo=  345.		! Value from P. Santos-Sanz = 599 +- 77 km (diam.)
	endif
c
c 2003AZ84 03 fev 2012:
c 
	if (iop.ge.148.and.iop.le.148) then
			R_pla=   400.		! Value from Michael E. Brown = 612 km (diam.)
			R_atmo=  400.		! Value from P. Santos-Sanz = 599 +- 77 km (diam.)
	endif
c
c 2003AZ84 15 nov 2014:
c 
	if (iop.ge.149.and.iop.le.149) then
			R_pla=   400.		! Value from Michael E. Brown = 612 km (diam.)
			R_atmo=  400.		! Value from P. Santos-Sanz = 599 +- 77 km (diam.)
	endif

cMommert, Michael; Harris, A. W.; Kiss, C.; Pál, A.; Santos-Sanz, P.; Stansberry, J.; Delsanti, A.; Vilenius, E.; Müller, T. G.; Peixinho, N.; Lellouch, E.; Szalai, N.; Henry, F.; Duffard, R.; Fornasier, S.; Hartogh, P.; Mueller, M.; Ortiz, J. L.; Protopapa, S.; Rengel, M.; Thirouin, A. (May 2012). "TNOs are cool: A survey of the trans-Neptunian region—V. Physical characterization of 18 Plutinos using Herschel-PACS observations". Astronomy & Astrophysics. 541: A93. arXiv:1202.3657 ￼. Bibcode:2012A&A...541A..93M. doi:10.1051/0004-6361/201118562.
	if (iop.ge.150.and.iop.le.150) then
			R_pla=   260.		! 
			R_atmo=  260.		! 
	endif




c	
c
c etoile du 01 juillet 2002 (origine ???).
c
c	if ((iop.eq.1).or.(iop.eq.2)) then
c	  ae=    17.d+00 + 01.d+00/60.d+00 + 52.319d+0/3600.d+00
c	  de= -( 12.d+00 + 38.d+00/60.d+00 + 56.50d+00/3600.d+00 )
c	endif
c
c etoile du 01 juillet 2002 (site MIT).
c
	if ((iop.eq.1).or.(iop.eq.2)) then
	  ae=    17.d+00 + 01.d+00/60.d+00 + 52.3141d+0/3600.d+00
	  de= -( 12.d+00 + 38.d+00/60.d+00 + 56.491d+00/3600.d+00 )
	endif

c etoile du 20 juillet 2002 (ICRS ???).
c
c position Papier McDonald Elliot
c
c       if ((iop.eq.3).or.(iop.eq.4)) then
c         ae=    17.d+00 + 00.d+00/60.d+00 + 18.002d+0/3600.d+00
c         de= -( 12.d+00 + 41.d+00/60.d+00 + 41.63d+00/3600.d+00 )
c       endif
c
c position janvier 2002:
c
c	if ((iop.eq.3).or.(iop.eq.4)) then
c	  ae=    17.d+00 + 00.d+00/60.d+00 + 18.006d+0/3600.d+00
c	  de= -( 12.d+00 + 41.d+00/60.d+00 + 41.91d+00/3600.d+00 )
c	endif
c
c position 20 mars 2002:
c
c	if ((iop.eq.3).or.(iop.eq.4)) then
c	  ae=    17.d+00 + 00.d+00/60.d+00 + 18.0235d+0/3600.d+00
c	  de= -( 12.d+00 + 41.d+00/60.d+00 + 41.968d+00/3600.d+00 )
c	endif
c
c position 24 avril 2002:
c
c	if ((iop.eq.3).or.(iop.eq.4)) then
c	  ae=    17.d+00 + 00.d+00/60.d+00 + 18.0235d+0/3600.d+00
c	  de= -( 12.d+00 + 41.d+00/60.d+00 + 41.951d+00/3600.d+00 )
c	endif
c
c position 23 mai 2002: P126a
c
c	if ((iop.eq.3).or.(iop.eq.4)) then
c	  ae=    17.d+00 + 00.d+00/60.d+00 + 18.0199d+0/3600.d+00
c	  de= -( 12.d+00 + 41.d+00/60.d+00 + 41.995d+00/3600.d+00 )
c	endif
c
c position 7 juin 2002: Michel Rapaport
c
c	if ((iop.eq.3).or.(iop.eq.4)) then
c	  ae=    17.d+00 + 00.d+00/60.d+00 + 18.0197d+0/3600.d+00
c	  de= -( 12.d+00 + 41.d+00/60.d+00 + 41.915d+00/3600.d+00 )
c	endif
c position 3 juillet 2002: Michel Rapaport, legeres modifications
c
c	if ((iop.eq.3).or.(iop.eq.4)) then
c	  ae=    17.d+00 + 00.d+00/60.d+00 + 18.0199d+0/3600.d+00
c	  de= -( 12.d+00 + 41.d+00/60.d+00 + 41.918d+00/3600.d+00 )
c	endif
c
c position 25 juin 2002: Agnes Fienga (T55 Pic)
c
c	if ((iop.eq.3).or.(iop.eq.4)) then
c	  ae=    17.d+00 + 00.d+00/60.d+00 + 17.99222d+0/3600.d+00
c	  de= -( 12.d+00 + 41.d+00/60.d+00 + 41.9489d+00/3600.d+00 )
c	endif
c
c position 27 juin 2002: Agnes Fienga (T120 OHP)
c
c	if ((iop.eq.3).or.(iop.eq.4)) then
c	  ae=    17.d+00 + 00.d+00/60.d+00 + 18.0151d+0/3600.d+00
c	  de= -( 12.d+00 + 41.d+00/60.d+00 + 42.034d+00/3600.d+00 )
c	endif
c
c position 5 juillet 2002: Ron Stone (mesures centroide + Da = +0.583"
c                                                         Dd = -0.192"
c	if ((iop.eq.3).or.(iop.eq.4)) then
c	  ae=    17.d+00 + 00.d+00/60.d+00 + 18.0364d+0/3600.d+00
c	  de= -( 12.d+00 + 41.d+00/60.d+00 + 42.081d+00/3600.d+00 )
c	endif
c
c  position 9 juillet 2002: Francois Colas T55 Pic
c
c 	if ((iop.eq.3).or.(iop.eq.4)) then
c 	  ae=    17.d+00 + 00.d+00/60.d+00 + 18.0186d+0/3600.d+00
c 	  de= -( 12.d+00 + 41.d+00/60.d+00 + 42.060d+00/3600.d+00 )
c 	endif
c
c  position 11 juillet 2002: JLX T1M Pic
c
c	if ((iop.eq.3).or.(iop.eq.4)) then
c	  ae=    17.d+00 + 00.d+00/60.d+00 + 18.0102d+0/3600.d+00
c	  de= -( 12.d+00 + 41.d+00/60.d+00 + 41.963d+00/3600.d+00 )
c	endif
c
c  position 25 juillet 2002: site MIT
c
	if ((iop.eq.3).or.(iop.eq.4)) then
	  ae=    17.d+00 + 00.d+00/60.d+00 + 18.03104d+0/3600.d+00
	  de= -( 12.d+00 + 41.d+00/60.d+00 + 42.028d+00/3600.d+00 )
	endif
c
c etoile du 09 aout 2002 (ICRS ???).
c
	if ((iop.eq.5).or.(iop.eq.6)) then
	  ae=    16.d+00 + 59.d+00/60.d+00 + 07.407d+0/3600.d+00
	  de= -( 12.d+00 + 47.d+00/60.d+00 + 20.30d+00/3600.d+00 )
	endif
c
c etoile du 21 aout 2002 (MIT site 16 aout 2002).
c
c	if (iop.eq.7) then
c	  ae=    16.d+00 + 58.d+00/60.d+00 + 49.43132d+0/3600.d+00
c	  de= -( 12.d+00 + 51.d+00/60.d+00 + 31.869d+00/3600.d+00 )
c	endif
c
c etoile du 21 aout 2002 (MIT site 16 aout 2002). Correction le 
c 7 novembre 2002: il y avait une erreur dans la r.a. ci-dessus (?),
c voir cahier du 7 novembre 2002.
c
	if (iop.eq.7) then
	  ae=    16.d+00 + 58.d+00/60.d+00 + 49.431793d+0/3600.d+00
	  de= -( 12.d+00 + 51.d+00/60.d+00 + 31.869d+00/3600.d+00 )
	endif
c
c etoile du 07 novembre 2002 (MIT site)
c
	if (iop.eq.8) then
	  ae=    17.d+00 + 04.d+00/60.d+00 + 01.816d+0/3600.d+00
	  de= -( 13.d+00 + 28.d+00/60.d+00 + 00.13d+00/3600.d+00 )
	endif
c
c etoile du 14 novembre 2003: catalogue Vizier, ICRS 2000 position en 2000?
c (29 octobre 2002)
c
c	if (iop.eq.9) then
c	  ae=    06.d+00 + 55.d+00/60.d+00 + 20.9639d+0/3600.d+00
c	  de= +( 22.d+00 + 06.d+00/60.d+00 + 07.901d+00/3600.d+00 )
c	endif
c
c etoile du 14 novembre 2003 (Titan): TYC 1343-1615-1 vers 00h11
c catalogue Vizier, ICRS 2000 + mvt propre propage du 01 janvier 2000 au
c 14 novembre 2003. Erreurs typiques en alpha et delta: 8 et 12 mas
c (25 novembre 2002):
c
	if ((iop.ge.9).and.(iop.le.13)) then
	  ae=    06.d+00 + 55.d+00/60.d+00 + 20.9672d+0/3600.d+00
	  de= +( 22.d+00 + 06.d+00/60.d+00 + 07.812d+00/3600.d+00 )
	endif
c
c etoile du 14 novembre 2003 (Titan): TYC 1343-1865-1 vers 06h44
c catalogue Vizier, ICRS 2000 + mvt propre propage du 01 janvier 2000 au
c 14 novembre 2003. Erreurs typiques en alpha et delta: 41 et 56 mas
c (25 novembre 2002):
c
c	if ((iop.ge.14).and.(iop.le.18)) then
c	  ae=    06.d+00 + 55.d+00/60.d+00 + 17.7733d+0/3600.d+00
c	  de= +( 22.d+00 + 06.d+00/60.d+00 + 01.168d+00/3600.d+00 )
c	endif
c 
c  Idem, mais avec position Manek (email Rubicon 01/07/03), erreurs 18 mas
c
c	
	if ((iop.ge.14).and.(iop.le.18)) then
	  ae=    06.d+00 + 55.d+00/60.d+00 + 17.7690d+0/3600.d+00
	  de= +( 22.d+00 + 06.d+00/60.d+00 + 01.226d+00/3600.d+00 )
	endif
c
c Etoile du 29 novembre 2003 (Triton), position site MIT Ron Stone FASTT
c (25 novembre 2002):
c
	if (iop.eq.20) then
	  ae=    20.d+00 + 52.d+00/60.d+00 + 45.908d+0/3600.d+00
	  de= -( 17.d+00 + 33.d+00/60.d+00 + 34.27d+00/3600.d+00 )
	endif
c
c Etoile du 15 decembre 2002 (Tethys), TYC 1310-24351-1 Vizier avec 
c mouvement propre pris en compte (12 decembre 2002):
c
	if (iop.eq.21) then
	  ae=    05.d+00 + 41.d+00/60.d+00 + 33.8871d+0/3600.d+00
	  de= +( 22.d+00 + 03.d+00/60.d+00 + 31.214d+00/3600.d+00 )
	endif
c
c Etoile du debut avril 2003 (Jupiter), TYC 1396-00214, voir email
c de JLX du 6 fevrier 2003 (parallaxe annuelle prise en compte)
c
	if (iop.eq.22) then
	  ae=    08.d+00 + 42.d+00/60.d+00 + 42.4731d+0/3600.d+00
	  de= +( 19.d+00 + 05.d+00/60.d+00 + 58.880d+00/3600.d+00 )
	endif
c
c HIP 9369, 10 octobre 1999
c voir mon rapport (parallaxe annuelle prise en compte)
c
	if (iop.eq.23) then
	  ae=    02.d+00 + 00.d+00/60.d+00 + 22.7321d+0/3600.d+00
	  de= +( 10.d+00 + 37.d+00/60.d+00 + 47.706d00/3600.d+00 )
	endif
c
c HIP 37084, 24 aout 2004, catalogue Hipparcos
c calcul 29 juin 2004 (parallaxe annuelle **non** prise en compte)
c
	if (iop.ge.24.and.iop.le.26) then
	  ae=    07.d+00 + 37.d+00/60.d+00 + 13.6462d+0/3600.d+00
	  de= +( 21.d+00 + 22.d+00/60.d+00 + 22.5486d00/3600.d+00 )
	endif
c
c P292, 22 mai 2005
c calcul mars 2004 (cf. classeur) 
c
	if ((iop.ge.27).and.(iop.le.29)) then
	  ae=    17.d+00 + 34.d+00/60.d+00 + 04.8684d+0/3600.d+00
	  de= -( 14.d+00 + 59.d+00/60.d+00 + 08.249d00/3600.d+00 )
	endif
c
c P313.2, Charon 11 juillet 2005
c position 01/09/04, site du MIT
c
	if (iop.eq.30) then
	  ae=    17.d+00 + 28.d+00/60.d+00 + 55.0174d+0/3600.d+00
	  de= -( 15.d+00 + 00.d+00/60.d+00 + 54.750d00/3600.d+00 )
	endif
c
c P313.2, Charon 11 juillet 2005
c position Stone (email Leslie Young 3/9/04) + calculs JLX (email 4/9/04)
c
	if (iop.eq.31) then
	  ae=    17.d+00 + 28.d+00/60.d+00 + 55.0130d+0/3600.d+00
	  de= -( 15.d+00 + 00.d+00/60.d+00 + 54.719d00/3600.d+00 )
	endif
c
c P313.2, Charon 11 juillet 2005
c position LNA (email Assafim, Andrei, da Silva Neto )
c
	if (iop.eq.32) then
	  ae=    17.d+00 + 28.d+00/60.d+00 + 55.01670d+0/3600.d+00
	  de= -( 15.d+00 + 00.d+00/60.d+00 + 54.7260d00/3600.d+00 )
	endif
c       
c P328.1, Charon 20 septembre 2005 
c position 02/09/04, site du MIT.
c Attention!: etoile double !!!
c
c        if (iop.eq.42) then
c          ae=    17.d+00 + 26.d+00/60.d+00 + 25.4567d+0/3600.d+00
c          de= -( 15.d+00 + 21.d+00/60.d+00 + 35.248d00/3600.d+00 )
c        endif
c
c position email L. Young 26/03/05
c
	if (iop.eq.42) then
	  ae=    17.d+00 + 26.d+00/60.d+00 + 25.4610d+0/3600.d+00
	  de= -( 15.d+00 + 21.d+00/60.d+00 + 35.241d00/3600.d+00 )
	endif
c
c C395, Charon 26 juillet 2006
c position D. Herald, email 26 juin 05
c
	if (iop.eq.43) then
	  ae=    17.d+00 + 36.d+00/60.d+00 + 51.725d+0/3600.d+00
	  de= -( 15.d+00 + 47.d+00/60.d+00 + 01.07d00/3600.d+00 )
	endif

c
c P363, Pluton 10 avril 2006
c position R. Stone, cf. document L. Young 23 septembre 2005
c +/- 102 mas en ra
c +/-  75 mas en dec
c
	if (iop.eq.44) then
	  ae=    17.d+00 + 46.d+00/60.d+00 + 06.8197d+0/3600.d+00
	  de= -( 15.d+00 + 46.d+00/60.d+00 + 10.332d00/3600.d+00 )
	endif
	
c
c P363, Pluton (ou Charon) 10 avril 2006
c position R. Behrend, cf. mail 6 mars 2006
c +/- 10 mas en ra
c +/- 10 mas en dec
c
	if (iop.eq.45.or.iop.eq.46) then
	  ae=    17.d+00 + 46.d+00/60.d+00 + 06.880d+0/3600.d+00
	  de= -( 15.d+00 + 46.d+00/60.d+00 + 10.11d00/3600.d+00 )
	endif
c
c VOISINE de P363, Pluton 10 avril 2006
c position R. Behrend, cf. mail 8 mars 2006
c erreur? 
c
	if (iop.eq.47) then
	  ae=    17.d+00 + 46.d+00/60.d+00 + 06.497d+0/3600.d+00
	  de= -( 15.d+00 + 46.d+00/60.d+00 + 11.24d00/3600.d+00 )
	endif
c
c P1 06 aout 2006
c position R. Behrend, cf. email 19 avril 2006
c 
c	if (iop.eq.48) then
c	  ae=    17.d+00 + 36.d+00/60.d+00 + 07.8431d+0/3600.d+00
c	  de= -( 15.d+00 + 49.d+00/60.d+00 + 32.800d0/3600.d+00 )
c	endif
c
c P1 06 aout 2006
c position R. Behrend, cf. email 15 juillet 2006
c 
	if (iop.eq.48) then
	  ae=    17.d+00 + 36.d+00/60.d+00 + 07.8427d+0/3600.d+00
	  de= -( 15.d+00 + 49.d+00/60.d+00 + 32.792d0/3600.d+00 )
	endif
c
c P384.2, Pluton 12 juin 2006
c position R. Stone 2004, cf. document L. Young 
c erreur 19 mas en ra, 24 mas en dec
c
	if (iop.eq.49) then
	  ae=    17.d+00 + 41.d+00/60.d+00 + 12.0940d+0/3600.d+00
	  de= -( 15.d+00 + 41.d+00/60.d+00 + 34.446d0/3600.d+00 )
	endif
c
c P384.2, Pluton 12 juin 2006
c position R. Behrend, cf. email 02 avril 2006
c
	if ((iop.ge.50).and.(iop.le.51)) then
	  ae=    17.d+00 + 41.d+00/60.d+00 + 12.0769d+0/3600.d+00
	  de= -( 15.d+00 + 41.d+00/60.d+00 + 34.485d0/3600.d+00 )
	endif
c
c Pluton 18 mars 2007
c position R. Behrend UCAC2 ---> ICRF et mouvement propre, cahier 21/02/07
c
c	if ((iop.ge.60).and.(iop.le.60)) then
c	  ae=    17.d+00 + 55.d+00/60.d+00 + 05.6969d+0/3600.d+00
c	  de= -( 16.d+00 + 28.d+00/60.d+00 + 34.389d0/3600.d+00 )
c	endif
c
c Pluton 18 mars 2007
c position R. Behrend UCAC2 ---> ICRF, mail (mail 4 mars 2007)
c
	if ((iop.ge.60).and.(iop.le.60)) then
	  ae=    17.d+00 + 55.d+00/60.d+00 + 05.6971d+0/3600.d+00
	  de= -( 16.d+00 + 28.d+00/60.d+00 + 34.360d0/3600.d+00 )
	endif
c
c Pluton 18 mars 2007
c position M. Assafim UCAC2 ---> ICRF et mouvement propre, cahier 21/02/07
c (erreurs 15 et 17 mas en ra et dec)
c
c	if ((iop.ge.61).and.(iop.le.61)) then
c	  ae=    17.d+00 + 55.d+00/60.d+00 + 05.693774d+0/3600.d+00
c	  de= -( 16.d+00 + 28.d+00/60.d+00 + 34.35947d0/3600.d+00 )
c	endif
c
c
c Pluton 18 mars 2007
c position M. Assafim UCAC2 ---> ICRF mail 05 mars 2004
c
	if ((iop.ge.61).and.(iop.le.61)) then
	  ae=    17.d+00 + 55.d+00/60.d+00 + 05.69579d+0/3600.d+00
	  de= -( 16.d+00 + 28.d+00/60.d+00 + 34.363d0/3600.d+00 )
	endif
c
c Pluton 12 mai 2007  A
c position R. Behrend ICRF, cahier 28/01/07
c
	if ((iop.ge.70).and.(iop.le.70)) then
	  ae=    17.d+00 + 53.d+00/60.d+00 + 31.9635d+0/3600.d+00
	  de= -( 16.d+00 + 22.d+00/60.d+00 + 47.017d0/3600.d+00 )
	endif
c
c Pluton 12 mai 2007  A
c position M. Assafim ICRF, cahier 28/01/07
c
	if ((iop.ge.71).and.(iop.le.71)) then
	  ae=    17.d+00 + 53.d+00/60.d+00 + 31.9590d+0/3600.d+00
	  de= -( 16.d+00 + 22.d+00/60.d+00 + 47.031d0/3600.d+00 )
	endif
c
c Pluton 12 mai 2007  B
c position R. Behrend ICRF, cahier 28/01/07
c
	if ((iop.ge.72).and.(iop.le.72)) then
	  ae=    17.d+00 + 53.d+00/60.d+00 + 32.1058d+0/3600.d+00
	  de= -( 16.d+00 + 22.d+00/60.d+00 + 47.365d0/3600.d+00 )
	endif
c
c Pluton 12 mai 2007  B
c position M. Assafim ICRF, cahier 28/01/07
c
	if ((iop.ge.73).and.(iop.le.73)) then
	  ae=    17.d+00 + 53.d+00/60.d+00 + 32.1006d+0/3600.d+00
	  de= -( 16.d+00 + 22.d+00/60.d+00 + 47.367d0/3600.d+00 )
	endif
c
c Pluton 14 juin 2007
c position UCAC2: 17 50 20.744, -16 22 42.18 J2000 + pm pendant 2721 jours, 
c -6.6 +/- 8.0 et 6.4 +/- 8.1 mas en ra et dec: -49.2 mas et +47.7 mas avec 
c incertitude 72 mas depuis 1998.5
c
c position R. Behrend ICRF, cahier 28/01/07:
c
c	if ((iop.ge.80).and.(iop.le.80)) then
c	  ae=    17.d+00 + 50.d+00/60.d+00 + 20.7422d+0/3600.d+00
c	  de= -( 16.d+00 + 22.d+00/60.d+00 + 42.212d0/3600.d+00 )
c	endif
c
c position R. Behrend ICRF, mail 16/04/07 et cahier 28/01/07:
c
c	if ((iop.ge.80).and.(iop.le.80)) then
c	  ae=    17.d+00 + 50.d+00/60.d+00 + 20.7442d+0/3600.d+00
c	  de= -( 16.d+00 + 22.d+00/60.d+00 + 42.214d0/3600.d+00 )
c	endif
c
c position R. Behrend ICRF, mail 02/06/07 et cahier 04/06/07:
c
	if ((iop.ge.80).and.(iop.le.80)) then
	  ae=    17.d+00 + 50.d+00/60.d+00 + 20.7441d+0/3600.d+00
	  de= -( 16.d+00 + 22.d+00/60.d+00 + 42.208d0/3600.d+00 )
	endif
c
c Pluton 14 juin 2007
c position M. Assafim ICRF, cahier 28/01/07:
c
c	if ((iop.ge.81).and.(iop.le.81)) then
c	  ae=    17.d+00 + 50.d+00/60.d+00 + 20.7412d+0/3600.d+00
c	  de= -( 16.d+00 + 22.d+00/60.d+00 + 42.236d0/3600.d+00 )
c	endif
c
c Pluton 14 juin 2007
c position M. Assafim ICRF, mail 22 mai 07:
c
	if ((iop.ge.81).and.(iop.le.81)) then
	  ae=    17.d+00 + 50.d+00/60.d+00 + 20.7392d+0/3600.d+00
	  de= -( 16.d+00 + 22.d+00/60.d+00 + 42.210d0/3600.d+00 )
	endif
c
c Pluton, 31 juillet 07
c position R. Behrend ICRF, mail 02/06/07 et cahier 04/06/07:
c
	if ((iop.ge.82).and.(iop.le.82)) then
	  ae=    17.d+00 + 45.d+00/60.d+00 + 41.9854d+0/3600.d+00
	  de= -( 16.d+00 + 29.d+00/60.d+00 + 31.613d0/3600.d+00 )
	endif
c
c Pluton, 31 juillet 07
c position M. Assafin ICRF, mail 31/05/07 et cahier 01/06/07:
c
c	if ((iop.ge.83).and.(iop.le.83)) then
c	  ae=    17.d+00 + 45.d+00/60.d+00 + 41.98502d+0/3600.d+00
c	  de= -( 16.d+00 + 29.d+00/60.d+00 + 31.647d0/3600.d+00 )
c	endif
c
c Mail Roberto 29 juillet 07 (erratum), avec correction ICRF:
c
	if ((iop.ge.83).and.(iop.le.83)) then
	  ae=    17.d+00 + 45.d+00/60.d+00 + 41.98943d+0/3600.d+00
	  de= -( 16.d+00 + 29.d+00/60.d+00 + 31.639d0/3600.d+00 )
	endif
c
c Pluton, 09 juin 07
c position R. Behrend UCAC2, catalogue 08 septembre 06, mon mail 18/05/07
c
	if ((iop.ge.84).and.(iop.le.84)) then
	  ae=    17.d+00 + 50.d+00/60.d+00 + 50.6517d+0/3600.d+00
	  de= -( 16.d+00 + 22.d+00/60.d+00 + 29.325d0/3600.d+00 )
	endif
c
c Pluton, 09 juin 07
c position R. Behrend UCAC2, catalogue 16 avril 07, mon mail 18/05/07
c
	if ((iop.ge.85).and.(iop.le.85)) then
	  ae=    17.d+00 + 50.d+00/60.d+00 + 50.6535d+0/3600.d+00
	  de= -( 16.d+00 + 22.d+00/60.d+00 + 29.297d0/3600.d+00 )
	endif
c
c Titania, 01 aout 03
c position + pm TYC2: 
c  22:15:54.55762 +/- 30 mas 
c -11:36:56.3395  +/- 30 mas
c
c UCAC2:
c  22:15:54.55433 +/- 20 mas 
c -11:36:56.3584  +/- 13 mas
c
c Nomad:
c  22:15:54.55430 +/- 20 mas 
c -11:36:56.3584  +/- 13 mas
c
c je choisis Nomad:
c
	if ((iop.ge.86).and.(iop.le.86)) then
	  ae=    22.d+00 + 15.d+00/60.d+00 + 54.55430d+0/3600.d+00
	  de= -( 11.d+00 + 36.d+00/60.d+00 + 56.3584d0/3600.d+00 )
	endif
c
c Triton, 21 mai 2008
c position Rio 27 avril 08, cf. mail 27/04/08
c
	if (iop.eq.87) then
	  ae=    21.d+00 + 46.d+00/60.d+00 + 11.049528d+0/3600.d+00
	  de= -( 13.d+00 + 46.d+00/60.d+00 + 45.9484716d0/3600.d+00 )
	endif
c
c Pluton et Charon, 22 juin  2008
c position Rio 02 juin 08 (mail Roberto)
c
	if ((iop.ge.88).and.(iop.le.89)) then
	  ae=    17.d+00 + 58.d+00/60.d+00 + 33.01376909d+00/3600.d+00
	  de= -( 17.d+00 + 02.d+00/60.d+00 + 38.34870554d+00/3600.d+00 )
	endif
c
c Pluton, 24 juin  2008
c position Rio 09 juin 08 (mail Roberto)
c
	if (iop.eq.90) then
	  ae=    17.d+00 + 58.d+00/60.d+00 + 22.392984d+00/3600.d+00
	  de= -( 17.d+00 + 02.d+00/60.d+00 + 49.349328d+00/3600.d+00 )
	endif
c
c 2002 MS4, 30 juin  2009
c position Rio 25 juin 09 (mail Marcelo Assafin)
c
	if (iop.eq.91) then
	  ae=    18.d+00 + 03.d+00/60.d+00 + 37.93133d+00/3600.d+00
	  de= -( 08.d+00 + 28.d+00/60.d+00 + 17.6971d+00/3600.d+00 )
	endif
c
c Pluton, 25 aout 2008
c position Rio 23 juillet 08 (mail Roberto)
c
	if (iop.eq.92) then
	  ae=    17.d+00 + 53.d+00/60.d+00 + 27.103164d+00/3600.d+00
	  de= -( 17.d+00 + 15.d+00/60.d+00 + 27.54630d+00/3600.d+00 )
	endif
c
c Varuna, 19 fevrier 2010
c position Rio 04 fevrier 2010 (mail Roberto)
c
	if (iop.eq.93) then
	  ae=    07.d+00 + 29.d+00/60.d+00 + 22.4714d+00/3600.d+00
	  de= +( 26.d+00 + 07.d+00/60.d+00 + 23.173d+00/3600.d+00 )
	endif
c
c Pluton, 14 fevrier 2010
c position Rio 10 novembre 2009, voir table
c g3_occ_data_pluto_2010_table
c
	if (iop.eq.94) then
	  ae=    18.d+00 + 19.d+00/60.d+00 + 14.3851d+00/3600.d+00
	  de= -( 18.d+00 + 16.d+00/60.d+00 + 42.313d+00/3600.d+00 )
	endif
c
c Pluton, 19 mai 2010
c position Rio 10 novembre 2009, voir table
c g3_occ_data_pluto_2010_table
c
	if (iop.eq.95) then
	  ae=    18.d+00 + 20.d+00/60.d+00 + 16.7941d+00/3600.d+00
	  de= -( 18.d+00 + 11.d+00/60.d+00 + 48.070d+00/3600.d+00 )
	endif
c
c Pluton, 04 juin 2010
c position Rio 10 novembre 2009, voir table
c g3_occ_data_pluto_2010_table
c
	if (iop.eq.96) then
	  ae=    18.d+00  + 18.d+00/60.d+00 + 47.9349d+00/3600.d+00
	  de= -( 18.d+00  + 12.d+00/60.d+00 + 51.794d+00/3600.d+00 )
	endif
c
c Pluton, 04 juillet 2010
c position Rio 02 juillet 2010
c
	if (iop.eq.97) then
	  ae=    18.d+00  + 15.d+00/60.d+00 + 42.1027d+00/3600.d+00
	  de= -( 18.d+00  + 16.d+00/60.d+00 + 41.124d+00/3600.d+00 )
	endif
c
c Ceres,  17 aout 2010
c position UCAC2
c
	if (iop.eq.98) then
	  ae=    17.d+00  + 18.d+00/60.d+00 + 29.008d+00/3600.d+00
	  de= -( 27.d+00  + 26.d+00/60.d+00 + 38.89d+00/3600.d+00 )
	endif

c
c Eris,  06 novembre 2010
c position Bill Owen, mail 4 novembre 2010, rms scatter 11 mas and 9mas in ra and dec, resp.
c
	if (iop.ge.99.and.iop.le.102) then
	  ae=    01.d+00  + 39.d+00/60.d+00 + 09.9392d+00/3600.d+00
	  de= -( 04.d+00  + 21.d+00/60.d+00 + 12.140d+00/3600.d+00 )
	endif
c
c Makemake, 23 avril 2011
c position JL Ortiz, position dans la demande DDT ESO 18 avril 2011
c
	if (iop.ge.110.and.iop.le.110) then
	  ae=    12.d+00  + 36.d+00/60.d+00 + 11.4140d+00/3600.d+00
	  de= +( 28.d+00  + 11.d+00/60.d+00 + 10.606d+00/3600.d+00 )
	endif
c
c Quaoar, 04 mai 2011
c position JL Ortiz, position dans l'update 03 mai 2011
c
	if (iop.ge.111.and.iop.le.111) then
	  ae=    17.d+00  + 28.d+00/60.d+00 + 50.8183d+00/3600.d+00
	  de= -( 15.d+00  + 27.d+00/60.d+00 + 42.657d+00/3600.d+00 )
	endif
c
c Pluton et Charon, 04 juin 2011
c position Rio 10 novembre 2009, voir table
c g3_occ_data_pluto_2011_table
c
	if (iop.ge.112.and.iop.le.115) then
	  ae=    18.d+00  + 27.d+00/60.d+00 + 53.8249d+00/3600.d+00
	  de= -( 18.d+00  + 45.d+00/60.d+00 + 30.725d+00/3600.d+00 )
	endif
c
c Pluton et Charon, 23 juin 2011
c position Rio 10 novembre 2009, voir table
c g3_occ_data_pluto_2011_table
c
	if (iop.ge.120.and.iop.le.121) then
	  ae=    18.d+00  + 25.d+00/60.d+00 + 55.475d+00/3600.d+00
	  de= -( 18.d+00  + 48.d+00/60.d+00 + 07.015d+00/3600.d+00 )
	endif
c
c Pluton et Hydra, 27 juin 2011 pour DE413 + offset MB
c position Rio 10 novembre 2009, voir table
c g3_occ_data_pluto_2011_table
c
	if (iop.ge.122.and.iop.le.123) then
	  ae=    18.d+00  + 25.d+00/60.d+00 + 29.010d+00/3600.d+00
	  de= -( 18.d+00  + 48.d+00/60.d+00 + 47.570d+00/3600.d+00 )
	endif
c
c Pluton et Charon, 23 juin 2011 pour DE413 + PLU021
c position Rio 10 novembre 2009, voir table
c g3_occ_data_pluto_2011_table
c
	if (iop.ge.124.and.iop.le.125) then
	  ae=    18.d+00  + 25.d+00/60.d+00 + 55.475d+00/3600.d+00
	  de= -( 18.d+00  + 48.d+00/60.d+00 + 07.015d+00/3600.d+00 )
	endif
c
c Pluton, 14 juin 2012 pour DE413 + PLU021
c voir table g4_occ_data_pluto_2012_table
c
	if (iop.ge.126.and.iop.le.126) then
		ae=    18.d+00  + 35.d+00/60.d+00 + 48.6931d+00/3600.d+00
		de= -( 19.d+00  + 17.d+00/60.d+00 + 43.617d+00/3600.d+00 )
	endif
c
c Pluton, 18 juillet 2012 pour DE413 + PLU022
c voir meilleure mesure OPD 6 Jul 2012, mail JC 11 jul 12
c
	if (iop.ge.127.and.iop.le.127) then
	 ae=    18.d+00  + 32.d+00/60.d+00 + 14.6720d+00/3600.d+00
	 de= -( 19.d+00  + 24.d+00/60.d+00 + 19.295d+00/3600.d+00 )
	endif
c
c Varuna, 08 janvier 2013
c etoile Raoul 06 et 07 jan 2013, donnees CAHA et Pic 5/6/7 jan 2013
c
	if (iop.ge.128.and.iop.le.130) then
	 ae=    07.d+00  + 49.d+00/60.d+00 + 37.0011d+00/3600.d+00
	 de= +( 26.d+00  + 25.d+00/60.d+00 + 51.979d+00/3600.d+00)
	endif
c
c 2002 KX14, 26 avril 2012
c etoile papier Alvaro, position NOMAD
c
	if (iop.ge.140.and.iop.le.140) then
	 ae=    16.d+00  + 35.d+00/60.d+00 + 04.279d+00/3600.d+00
	 de= -( 22.d+00  + 15.d+00/60.d+00 + 22.75d+00/3600.d+00)
	endif
c
c Pluton, 04 mai 2013 pour DE413 + PLU031
c position IAG mai 2013, voir mail Felipe 6 mai 2013, dossier "hydranixocc"
c
	if (iop.ge.141.and.iop.le.141) then
		ae=    18.d+00  + 47.d+00/60.d+00 + 52.5322d+00/3600.d+00
		de= -( 19.d+00  + 41.d+00/60.d+00 + 24.3738d+00/3600.d+00 )
	endif
c
c Chariklo, 03 juin 2013 pour JPL20
c position OPD mai 2013, voir carte Felipe mai 2013
c
	if (iop.ge.142.and.iop.le.142) then
		ae=    16.d+00  + 56.d+00/60.d+00 + 06.4876d+00/3600.d+00
		de= -( 40.d+00  + 31.d+00/60.d+00 + 30.205d+00/3600.d+00 )
	endif
c
c Chariklo, 29 avril 2014
c position OPD 06 sept. 2013, voir lail/carte Felipe 17 avril 2014
c
	if (iop.ge.143.and.iop.le.143) then
		ae=    17.d+00  + 39.d+00/60.d+00 + 02.1336d+00/3600.d+00
		de= -( 38.d+00  + 52.d+00/60.d+00 + 48.801d+00/3600.d+00 )
	endif
c
c Chariklo, 16 fevrier 2014
c position OPD 08sep13_T60
c
	if (iop.ge.144.and.iop.le.144) then
		ae=    17.d+00  + 35.d+00/60.d+00 + 55.3333d+00/3600.d+00
		de= -( 38.d+00  + 05.d+00/60.d+00 + 17.184d+00/3600.d+00 )
	endif
c
c Chariklo, 28 juin 2014
c position OPD 3nov14
c
	if (iop.ge.145.and.iop.le.145) then
		ae=    17.d+00  + 24.d+00/60.d+00 + 50.3800d+00/3600.d+00
		de= -( 38.d+00  + 41.d+00/60.d+00 + 05.618d+00/3600.d+00 )
	endif
	
c
c 2007 uk 126, 15 nov 2014
c position Ortiz  30may14_T160
c
	if (iop.ge.146.and.iop.le.146) then
		ae=    4.d+00  + 29.d+00/60.d+00 + 30.6265d+00/3600.d+00
		de= -( 0.d+00  + 28.d+00/60.d+00 + 20.742d+00/3600.d+00 )
	endif	
	



c
c 2003AZ84, WFI star position 01/08/2011
c
	if (iop.ge.147.and.iop.le.147) then
	  ae=    07.d+00  + 43.d+00/60.d+00 + 41.8220d+00/3600.d+00
	  de= +( 11.d+00  + 30.d+00/60.d+00 + 23.569d+00/3600.d+00 )
	endif	

c
c 2003 AZ84, Moyenne J.Lecacheux 01/fev/2012 star position, confirmee par 1,6m OPD 17/fev/2012
c
	if (iop.ge.148.and.iop.le.148) then
	  ae=  07.d+00  + 45.d+00/60.d+00 + 54.76965d+00/3600.d+00
	  de=  11.d+00  + 12.d+00/60.d+00 + 43.0933d+00/3600.d+00
	endif


c
c 2003 AZ84
c
	if (iop.ge.149.and.iop.le.149) then
	  ae=  08.d+00  + 03.d+00/60.d+00 + 51.2981d+00/3600.d+00
	  de=  09.d+00  + 57.d+00/60.d+00 + 18.729d+00/3600.d+00
	endif

c
c 2003 VS2
c
	if (iop.ge.150.and.iop.le.150) then
c GAIA
	  ae=  04.d+00  + 48.d+00/60.d+00 + 32.1390d+00/3600.d+00
	  de=  33.d+00  + 58.d+00/60.d+00 + 36.407d+00/3600.d+00
c NOMAD 
c	  ae=  04.d+00  + 48.d+00/60.d+00 + 32.1457d+00/3600.d+00
c	  de=  33.d+00  + 58.d+00/60.d+00 + 36.365d+00/3600.d+00
c UC4 
c	  ae=  04.d+00  + 48.d+00/60.d+00 + 32.134d+00/3600.d+00
c	  de=  33.d+00  + 58.d+00/60.d+00 + 36.57d+00/3600.d+00

	endif




	
	call dms (ae, iahe, iame, asece)	! on sauve la position de l'etoile dans fort.28,
	call dms (de, idhe, idme, dsece)	! pour utiliser dans fit_d2_ksi_eta.f
	write (28,*) iahe, iame, asece
	write (28,*) idhe, idme, dsece

c	write(*,*)  'Coordonnees astrometriques de l''etoile a l''epoque?'
c	write(*,*)  'sec et abs[arcsec] seulement'
c	read(*,*)  asec, dsec
c	ae=    17.d+00 + 01.d+00/60.d+00 + asec/3600.d+00
c	de= -( 12.d+00 + 38.d+00/60.d+00 + dsec/3600.d+00 )
c
c ephemeride de la planete
c
	call ephem (t,at,dt,delta,npt,iop)
c
c
c on utilise les formules:
c 
c xi=  -D*cos(delta)*sin(alpha-alphae)
c eta= -D*{sin(delta-deltae) + 2cos(delta)*sin(deltae)*sin^2[(alpha-alphae)/2]}
c ou: 
c (alpha, delta )= coord. planete
c (alphae,deltae)= coord. etoile
c
	do i= 1, npt
	 cdt= dcos( dt(i)*radd )
	 satmae=  dsin( ( at(i)-ae )*radh )
	 satmae2= (dsin((at(i)-ae)*(radh/2.d+00)))**2
	 sdtmde=  dsin( ( dt(i)-de )*radd )
	 sde= dsin( de*	radd )
	 xi(i)=  -ua*delta(i)*cdt*satmae
	 eta(i)= -ua*delta(i)*( sdtmde + 2.d+00*cdt*sde*satmae2 )
	 eta(i)= -eta(i) ! question d'orientation pour tracer sur la Terre
!	 write(*,*) sngl(t(i)), sngl(xi(i)), sngl(eta(i))
	enddo

	do i= 2, npt-1
		delta_t= ( t(i+1)-t(i-1) )*3600.d0
		vxi = (  xi(i+1) -  xi(i-1) )/delta_t
		veta= ( eta(i+1) - eta(i-1) )/delta_t
		v   = dsqrt( vxi**2 + veta**2 )
c		write(*,*)  sngl(t(i)), sngl(v)
		ux=  vxi/v
		uy= veta/v
		dkm= ua*delta(i)
		write(20,*) t(i),  xi (i), dkm 
		write(21,*) t(i),  eta(i)
		write(22,*) xi(i), eta(i), t(i), dsqrt(xi(i)**2 + eta(i)**2)
		write(23,*) xi(i) - uy*R_pla,  eta(i) + ux*R_pla,  t(i)
		write(24,*) xi(i) - uy*R_atmo, eta(i) + ux*R_atmo, t(i)
c si on veut mvt dans le ciel en arcsec:
		call dms (t(i), ih_tu,im_tu,sec_tu)
		dracosd= xi(i)/(dkm*arcsec)
		ddec   = eta(i)/(dkm*arcsec)
		write(25,*) -sngl(dracosd), sngl(ddec), sngl(t(i))
		write(26,*) xi(i), -eta(i), t(i)					! pour avoir le bon signe de eta
		write(27,'(A,I2,A,I2,A,f5.2,1X,f6.2,1x,f6.2)') 
     *  ':', ih_tu, ':', im_tu, ':', sngl(sec_tu), -sngl(dracosd), sngl(ddec)
	enddo

	do i= npt-1, 2, -1
		delta_t= ( t(i+1)-t(i-1) )*3600.d0
		vxi = (  xi(i+1) -  xi(i-1) )/delta_t
		veta= ( eta(i+1) - eta(i-1) )/delta_t
		v   = dsqrt( vxi**2 + veta**2 )
c		write(*,*)  sngl(t(i)), sngl(v)
		ux=  vxi/v
		uy= veta/v
		write(23,*) xi(i) + uy*R_pla,  eta(i) - ux*R_pla,  t(i)
		write(24,*) xi(i) + uy*R_atmo, eta(i) - ux*R_atmo, t(i)
	enddo
c
c
c
	stop
	end
*********************************** fin de main *****************************

	subroutine ephem (t,at,dt,delta,k,iop)
	implicit real*8 (a-h,o-z)
	dimension t(1000),at(1000),dt(1000),delta(1000)
	character*50 fichier
	character*11 date
	character*1  column

	pi=   dacos(-1.d+0)
	radd= pi/180.d+00
	radh= pi/12.d+00

	if (iop.eq.1)  fichier= 'ephem_pluton.01jul2002'
	if (iop.eq.2)  fichier= 'ephem_charon.01jul2002'
	if (iop.eq.3)  fichier= 'ephem_pluton.20jul2002'
	if (iop.eq.4)  fichier= 'ephem_charon.20jul2002'
	if (iop.eq.5)  fichier= 'ephem_pluton.09aug2002'
	if (iop.eq.6)  fichier= 'ephem_charon.09aug2002' 
	if (iop.eq.7)  fichier= 'ephem_pluton.21aug2002' 
	if (iop.eq.8)  fichier= 'ephem_pluton.07nov2002'

	if  (iop.eq.9) fichier= 'ephem_titan_jpl.14nov2003'
	if (iop.eq.10) fichier= 'ephem_titan_bdl.14nov2003'
	if (iop.eq.11) fichier= 'ephem_titan_bs_tass7.14nov2003'
	if (iop.eq.12) fichier= 'ephem_titan_bs_tass8.14nov2003'
	if (iop.eq.13) fichier= 'ephem_titan_jac.14nov2003'

	if (iop.eq.14) fichier= 'ephem_titan_jpl.14nov2003'
	if (iop.eq.15) fichier= 'ephem_titan_bdl.14nov2003'
	if (iop.eq.16) fichier= 'ephem_titan_bs_tass7.14nov2003'
	if (iop.eq.17) fichier= 'ephem_titan_bs_tass8.14nov2003'
	if (iop.eq.18) fichier= 'ephem_titan_jac.14nov2003'
	
	if (iop.eq.20) fichier= 'ephem_triton.29nov2003'
	if (iop.eq.21) fichier= 'ephem_tethys.15dec2002'
	if (iop.eq.22) fichier= 'ephem_jupiter.01apr2003'
	if (iop.eq.23) fichier= 'ephem_jupiter.10oct1999'
	
	if (iop.eq.24) fichier= 'ephem_titan_bdl.24aug2004'
	if (iop.eq.25) fichier= 'ephem_titan_bdl_vsop87.24aug2004'
	if (iop.eq.26) fichier= 'ephem_titan_bdl_5sec.24aug2004'

	if (iop.eq.27) fichier= 'ephem_pluton_jpl_vlt.22mai2005'
	if (iop.eq.28) fichier= 'ephem_pluton_jpl_geo.22mai2005'
	if (iop.eq.29) fichier= 'ephem_pluton_jpl_cfh.22mai2005'

	if (iop.eq.30) fichier= 'ephem_charon_jpl.11jul2005'
	if (iop.eq.31) fichier= 'ephem_charon_jpl_de413.11jul2005'
	if (iop.eq.32) fichier= 'ephem_charon_jpl_de413.11jul2005'
	
	if (iop.eq.42) fichier= 'ephem_charon_jpl_de413.20sep2005' 

	if (iop.eq.43) fichier= 'ephem_charon_jpl_de413.26jul2006' 
	
	if (iop.eq.44) fichier= 'ephem_pluton_geo.10avril2006'
!	if (iop.eq.45) fichier= 'ephem_pluton.10avril2006'				! ancienne DE413 (2006) 'ephem_pluton_etendue_geo.10avril2006'
	if (iop.eq.45) fichier= 'ephem_pluton_bary.10avril2006'			! si on veut la barycentrique
	if (iop.eq.46) fichier= 'ephem_charon_etendue_geo.10avril2006'
	if (iop.eq.47) fichier= 'ephem_pluton_etendue_geo.10avril2006'
	
	if (iop.eq.48) fichier= 'ephem_pluton.06aout2006'	
	
	if (iop.eq.49) fichier= 'ephem_pluton.12juin2006'
	if (iop.eq.50) fichier= 'ephem_pluton.12juin2006'	
	if (iop.eq.51) fichier= 'ephem_pluton.12juin2006'

	if (iop.eq.60) fichier= 'ephem_pluton.18mars2007'
	if (iop.eq.61) fichier= 'ephem_pluton.18mars2007'
		
	if (iop.eq.70) fichier= 'ephem_pluton.12mai2007'
	if (iop.eq.71) fichier= 'ephem_pluton.12mai2007'
	if (iop.eq.72) fichier= 'ephem_pluton.12mai2007'
	if (iop.eq.73) fichier= 'ephem_pluton.12mai2007'
				
!	if (iop.eq.80) fichier= 'ephem_pluton.14juin2007'
!	if (iop.eq.81) fichier= 'ephem_pluton.14juin2007'
	if (iop.eq.80) fichier= 'ephem_pluton_bary.14juin2007'
	if (iop.eq.81) fichier= 'ephem_pluton_bary.14juin2007'		! si on veut la barycentrique

	if (iop.eq.82) fichier= 'ephem_pluton.31jul2007'
	if (iop.eq.83) fichier= 'ephem_pluton.31jul2007'
			
	if (iop.eq.84) fichier= 'ephem_pluton.09jun2007'
	if (iop.eq.85) fichier= 'ephem_pluton.09jun2007'

	if (iop.eq.86) fichier= 'ephem_titania.01aug2003'

	if (iop.eq.87) fichier= 'ephem_triton.21mai2008'

	if (iop.eq.88) fichier= 'ephem_pluton.22jun2008'
	if (iop.eq.89) fichier= 'ephem_charon.22jun2008'

	if (iop.eq.90) fichier= 'ephem_pluton.24jun2008'

	if (iop.eq.91) fichier= 'ephem_mS4.30jun2009'

	if (iop.eq.92) fichier= 'ephem_pluton.25aug2008'

	if (iop.eq.93) fichier= 'ephem_varuna.19fev2010'

	if (iop.eq.94) fichier= 'ephem_pluton.14fev2010'

	if (iop.eq.95) fichier= 'ephem_pluton.19may2010'

	if (iop.eq.96) fichier= 'ephem_pluton.04jun2010'

	if (iop.eq.97) fichier= 'ephem_pluton.04jul2010'

	if (iop.eq.98) fichier= 'ephem_ceres.17aug2010'

	if (iop.eq.99) fichier=  'ephem_eris_JPL28.06nov2010'
	if (iop.eq.100) fichier= 'ephem_eris_JPL29.06nov2010'
	if (iop.eq.101) fichier= 'ephem_eris_JPL29_OPD.05nov2010'
	if (iop.eq.102) fichier= 'ephem_eris_JPL29_geo.05nov2010'

!	if (iop.eq.110) fichier= 'ephem_makemake_JPL43.23apr11'
	if (iop.eq.110) fichier= 'ephem_makemake_JPL46.23apr11'

	if (iop.eq.111) fichier= 'ephem_quaoar_JPL17.04may11'

	if (iop.eq.112) fichier= 'ephem_pluton.04jun2011'
	if (iop.eq.113) fichier= 'ephem_charon.04jun2011'
	if (iop.eq.114) fichier= 'ephem_pluton_DE413_MB.04jun2011'
	if (iop.eq.115) fichier= 'ephem_charon_DE413_MB.04jun2011'

	if (iop.eq.120) fichier= 'ephem_pluto_DE413_MB.23jun2011'
	if (iop.eq.121) fichier= 'ephem_charon_DE413_MB.23jun2011'

	if (iop.eq.122) fichier= 'ephem_pluto_DE413_MB.27jun2011'
	if (iop.eq.123) fichier= 'ephem_hydra_DE413_MB.27jun2011'

	if (iop.eq.124) fichier= 'ephem_pluto_DE413_PLU021.23jun2011'
	if (iop.eq.125) fichier= 'ephem_charon_DE413_PLU021.23jun2011'

	if (iop.eq.126) fichier= 'ephem_pluto_DE413_PLU021.14jun2012'

	if (iop.eq.127) fichier= 'ephem_pluto_DE413_PLU022.18jul2012'

	if (iop.eq.128) fichier= 'ephem_varuna_jpl25.05jan2013'
	if (iop.eq.129) fichier= 'ephem_varuna_jpl25.06jan2013'
	if (iop.eq.130) fichier= 'ephem_varuna_jpl25.08jan2013'

	if (iop.eq.140) fichier= 'ephem_kx14_jpl9.26avr2012'

	if (iop.eq.141) fichier= 'ephem_pluto_DE413_PLU031.04mai2013'

	if (iop.eq.142) fichier= 'ephem_chariklo_JPL20.03jun2013'

	if (iop.eq.143) fichier= 'ephem_chariklo_JPL21.29apr2014'

	if (iop.eq.144) fichier= 'ephem_chariklo_JPL21.16fev2014'

	if (iop.eq.145) fichier= 'ephem_chariklo_JPL22.28jun2014'
	if (iop.eq.146) fichier= 'ephem_2007UK126_15nov2014'


	if (iop.eq.147) fichier= 'ephem_2003AZ84_08jan2011'
	if (iop.eq.148) fichier= 'ephem_2003AZ84_03fev2012'
	if (iop.eq.149) fichier= 'ephem_2003AZ84_15nov2014'
	if (iop.eq.150) fichier= 'ephem_2003VS2_07nov2014'




! this section to guess the format of the ephemeris file
	open (20,file=fichier,status='old',form='formatted')
	imode= 0
	read(20,*,err=98) ith,tmn, iah,iamn,asec, idd,idmn,dsec, dis		! lecture ancien format
	goto 96
 98	imode=1
	rewind (20) 
	read(20,*,err=97) ith,tmn, a_deg, d_deg,                 dis		! lecture autre format
 	rewind (20) ; goto 96
 97	imode=2
	rewind (20) 
	read(20,'(1x,A11,1x,I2,A1,I2,5x,f11.7,1x,f11.7,1x,f16.13)',err=97) date, ith,column,imn, a_deg, d_deg, dis	! lecture autre format
	tmn= dfloat(imn)
 	goto 96
! format detected, then rewind the input file and proceed

96	continue
	rewind (20) 

	k=     0
	do i= 1, 10000
		if (imode.eq.0) read(20,*,err=99,end=99) ith,tmn, iah,iamn,asec, idd,idmn,dsec, dis
		if (imode.eq.1) read(20,*,err=99,end=99) ith,tmn, a_deg, d_deg,                 dis
		if (imode.eq.2) then
		 read(20,'(1x,A11,1x,I2,A1,I2,5x,f11.7,1x,f11.7,1x,f16.13)',err=99,end=99) date, ith,column,imn, a_deg, d_deg, dis
		 tmn= dfloat(imn)
		endif 
 		k= k+1
		t(k)=  dfloat(ith) + tmn/60.d0
		delta(k)= dis
		
		if (imode.eq.0) then
		 at(k)= dfloat(iah) + dfloat(iamn)/60.d0 + asec/3600.d0
		 if (idd.eq.0) write(*,*)  'Probleme sur le signe de delta!'
		 if (idd.lt.0) then
			idmn= -idmn ; dsec= -dsec
		 endif
		 dt(k)= dfloat(idd) + dfloat(idmn)/60.d0 + dsec/3600.d0
		endif
		
		if (imode.eq.1) then ; at(k)= a_deg/15.d0 ; dt(k)= d_deg ; endif
		if (imode.eq.2) then ; at(k)= a_deg/15.d0 ; dt(k)= d_deg ; endif
	enddo
 99	continue
	close (20)
	write(*,*)  k, ' points lus dans ', fichier
	if (imode.eq.0) write (*,*) 'format coordonnees: h:min:sec, deg:amin:asec'
	if (imode.eq.1) write (*,*) 'format coordonnees: deg.xxx, deg.xxx'

	write(*,*)  'Correction a l''ephemeride de la planete?' 
	write(*,*)  'D_alpha*cos(delta) et D_delta en **mas**'
	read(*,*)  dat_cosdt, ddt
	dat_cosdt= dat_cosdt/1.d3
	ddt      =       ddt/1.d3

	do i= 1, k
                dat= dat_cosdt/( 3600.d0*15.d0*dcos(dt(i)*radd) )
                at(i)= at(i) + dat
                dt(i)= dt(i) + ddt/3600.d0
	enddo

	return
	end
******************************** fin de ephem ***************************

      subroutine dms(heure, ihh,imn,sec)
c
c	ce sous-programme transforme l'heure decimale
c	'heure' en heure, minute, seconde: ihh:imn:sec.
c
	implicit real*8 (a-h,o-z)
	seuil= 1.d-2
c
	ihh= int(heure)
	amn= (heure - dfloat(ihh))*60.
	imn= int(amn)
	sec= (amn - dfloat(imn))*60.
c
	if(heure.lt.0.and.ihh.ne.0) then
		imn= iabs(imn)
		sec= dabs(sec)
	endif
c
	if ((60.0d0-sec).le.seuil) then
	 sec= 0.d0
	 imn= imn+1
	endif

	if (imn.eq.60) then
	 imn= 0
	 ihh= ihh+1
	endif

	return
	end
c
******************************* FIN DE DMS ********************************
