c
c	Il faut rentrer:
c
c	Le code de la date, iobs, voir liste ci-dessous
c
c	Un code special pour preciser si l'etoile du 2 aout 88 est passee
c	au sud ou au nord du centre de Neptune, iopns.
c
c	Le code qui determine le pole de Neptune selon que Triton est pris
c	en compte ou non, itriton.
c
c	L'heure, "heure", en TU decimal (ex. 10.5= 10 h 30 mn)
c
c	Le programme fournit: ....
c
c	fort.21: ksi, eta, rho, TU_sec
c       fort.22: ksig, etag geocentrique OU pour utilisation dans "map" (en faisant etag= -etag)
c	fort.24: rayon dans plan ciel, flux, heure
c	fort.25: rayon dans plan anneau, flux
c	fort.26: x, y dans plan anneau, flux
c	fort.27: hhmmss, dra*cos(del), ddec de Pluton /t a l'etoile (pour "tracking table" NACO)
c	fort.28: t (sec), 1. pour creer une courbe artificielle a utiliser dans omegat20 
c
	implicit real*8 (a-h,o-z)
	parameter (nmax=50000)
	dimension t_obs(nmax), flux_obs(nmax)
	real*8 ksi,ksig,ksip,ksir,ksif,ksigf,ksipf,ksirf     ! coordonnee ksi
	integer*4 annee 				     ! annee
 100	format(i3,1x,i2,1x,f7.3)
	pi= dacos(-1.d+00)
	rad= pi/180.d0
	arcsec= pi/(180.d0*3600.d0)
	c= 299792.458d+00				     ! NB: km/sec

	write(*,*) '1= heures choisies une a une'
	write(*,*) '2= heures choisies par pas constant entre deux limites'
	write(*,*) '3= heures choisies a partir d''un fichier d''occultation'
	read(*,*) imode

c	write(*,*)  'Echelle de hauteur? (km)'
c	read(*,*)  Hkm
c	dksi0= 0.d0
c	deta0= 0.d0
c	dt=    0.d0

	write(*,*)  'Corrections sur ksi et eta (km)'
	read(*,*)  dksi0, deta0

	if (imode.eq.2) then 
 	 write(*,*)  'Heure debut, heure fin et pas? (toutes en sec)'
	 read(*,*)  heure_0, heure_1, pas_heure
	 heure_0= heure_0/3600.d0
	 heure_1= heure_1/3600.d0
	 heure= heure_0
	 heurf= heure + 1.d0/3600.d0
	endif
	
	if (imode.eq.3) then
	 call lecture (nmax, npt,t_obs,flux_obs)
	 i_t_obs= 1	 
	 heure_0= t_obs(i_t_obs)/3600.d0
	 heure_1= t_obs(    npt)/3600.d0
	 heure= heure_0
	 heurf= heure + 1.d0/3600.d0	 
	endif
c
c .........................................................................
c
c On rentre les parametres de l'observation
c
 1	continue
c	write(*,*)  'Code de l''observation ?'
	write(*,*)  'Code observation, correction sur l''heure et ecriture?'
c	write(*,*)  '0=  15 juin      1983, HOBART'
c	write(*,*)  '1=  15 juin      1983, CFHT'
c	write(*,*)  '2=  22 juillet   1984, ESO'
c	write(*,*)  '3=  22 juillet   1984, CTIO'
c	write(*,*)  '4=  20 aout      1985, ESO'
c	write(*,*)  '5=  20 aout      1985, CTIO'
c	write(*,*)  '6=  20 aout      1985, CFHT'
c	write(*,*)  '7=  23 avril     1986, CFHT'
c	write(*,*)  '8=  23 aout      1986, CFHT'
c	write(*,*)  '9=  22 juin      1987, ESO'
c	write(*,*)  '10= 09 juillet   1987, CFHT'
c	write(*,*)  '11= 02 aout      1988, ESO'
c	write(*,*)  '12= 25 aout      1988, CFHT'
c	write(*,*)  '13= 12 septembre 1988, PIC'
c	write(*,*)  '14= 12 septembre 1988, OHP'
c	write(*,*)  '15= 07 juillet   1989, ESO'
c	write(*,*)  '16= 07 juillet   1989, TEIDE'
c	write(*,*)  '17= 07 juillet   1989, PIC'
c	write(*,*)  '18= 03 juillet   1989, CATANIA'
c	write(*,*)  '19= 03 juillet   1989, MEUDON'
c	write(*,*)  '20= 03 juillet   1989, ODESSA'
c	write(*,*)  '21= 03 juillet   1989, PIC'
c	write(*,*)  '22= 03 juillet   1989, EIN HAROD'
c	write(*,*)  '23= 03 juillet   1989, MANLEY'
c	write(*,*)  '24= 03 juillet   1989, WISE'
c	write(*,*)  '25= 03 juillet   1989, VATICAN'
c	write(*,*)  '26= 03 juillet   1989, HATFIELD'
c	write(*,*)  '27= 03 juillet   1989, NEZEL'
c	write(*,*)  '28= 03 juillet   1989, DENZAU'
c	write(*,*)  '29= 03 juillet   1989, RGO'
c	write(*,*)  '30= 11 Juillet   1990, ESO'
c	write(*,*)  '31= 18 aout      1991, CFHT'
c	write(*,*)  '32= 11 juillet   1992, ESO'
c	write(*,*)  '33= 13 decembre  1989, KITT PEAK'

c	write(*,*)  '34= 18 juillet   1997, BUNDABERG'
c	write(*,*)  '35= 18 juillet   1997, DUCABROOK'
c	write(*,*)  '36= 18 juillet   1997, LOCHINGTON'
c	write(*,*)  '37= 18 juillet   1997, BROWNSVILLE'
c	write(*,*)  '38= 18 juillet   1997, CHILLAGOE'

c	write(*,*)  '39= 13 mai	 1971, PRETORIA'
c	write(*,*)  '40= 10 octobre   1999, ESO'
c	write(*,*)  '41= 10 octobre   1999, PARANAL'
c	write(*,*)  '42= 10 octobre   1999, MONT MEGANTIC'
c	write(*,*)  '43= 10 octobre   1999, PIC DU MIDI'
c	write(*,*)  '44= 10 octobre   1999, CATALINA STATION'
c	write(*,*)  '45= 10 octobre   1999, KITT PEAK'
c	write(*,*)  '50= 03 decembre  1999, ESO'
c
c	write(*,*)  '51= 08 septembre 2001, SALINAS'
c	write(*,*)  '52= 08 septembre 2001, EL BUERAN (CUENCA)'
c	write(*,*)  '53= 08 septembre 2001, MARACAIBO'
c	write(*,*)  '54= 08 septembre 2001, BOBARES'
c	write(*,*)  '55= 08 septembre 2001, ARVAL CARACAS'
c	write(*,*)  '56= 08 septembre 2001, ARIKOK, ARUBA'
c	write(*,*)  '57= 08 septembre 2001, FORT DE FRANCE'
c	write(*,*)  '58= 08 septembre 2001, TOBAGO'
c	write(*,*)  '59= 08 septembre 2001, TRINIDAD'
c	write(*,*)  '60= 08 septembre 2001, BRIDGETOWN, BARBADOS'
c	write(*,*)  '61= 08 septembre 2001, PONTA DELGADA, ACORES'
c	write(*,*)  '62= 08 septembre 2001, OEIRAS'
c	write(*,*)  '63= 08 septembre 2001, PORTIMAO'
c	write(*,*)  '64= 08 septembre 2001, COAA ALGARVE (PORTIMAO)'
c	write(*,*)  '65= 08 septembre 2001, LINHACEIRA'
c	write(*,*)  '66= 08 septembre 2001, ALVITO'
c	write(*,*)  '67= 08 septembre 2001, GRANADA'
c	write(*,*)  '68= 08 septembre 2001, ALCUBLAS'
c	write(*,*)  '69= 08 septembre 2001, BORDEAUX'
c	write(*,*)  '70= 08 septembre 2001, PIC'
c	write(*,*)  '71= 08 septembre 2001, ST MAURICE LA CLOUERE (POITOU)'
c	write(*,*)  '72= 08 septembre 2001, PEZENAS'
c	write(*,*)  '73= 08 septembre 2001, MAUGUIO'
c	write(*,*)  '74= 08 septembre 2001, NIMES'
c	write(*,*)  '75= 08 septembre 2001, ORFEUILLES (MARSEILLE)'
c	write(*,*)  '76= 08 septembre 2001, OHP'
c	write(*,*)  '77= 08 septembre 2001, BINFIELD'
c	write(*,*)  '78= 08 septembre 2001, WORTH HILL, BOURNEMOUTH'
c	write(*,*)  '79= 08 septembre 2001, TEVERSHAM'
c	write(*,*)  '80= 08 septembre 2001, PLATEAU D''ALBION'
c	write(*,*)  '81= 08 septembre 2001, H. FABRE, MARSEILLE'
c	write(*,*)  '82= 08 septembre 2001, GUITALENS'
c	write(*,*)  '83= 08 septembre 2001, SABADELL (ARDANUY)'
c   write(*,*)  '84= 08 septembre 2001, SALON'
c	write(*,*)  '85= 08 septembre 2001, WELA, ARUBA'
c	write(*,*)  '86= 08 septembre 2001, ST MAURICE CAZEVIEILLE'
c	write(*,*)  '87= 08 septembre 2001, MARINHA GRANDE'
c	write(*,*)  '88= 08 septembre 2001, ALMEIRIM'
c	write(*,*)  '89= 08 septembre 2001, CARCAVELOS'
c	write(*,*)  '90= 08 septembre 2001, SETUBAL'
c	write(*,*)  '91= 08 septembre 2001, ALCACER DO SAL'
c	write(*,*)  '92= 08 septembre 2001, BARCELONA 1'
c	write(*,*)  '93= 08 septembre 2001, ST SAVINIEN'
c	write(*,*)  '94= 08 septembre 2001, ST MARTIN DE CRAU'
c	write(*,*)  '95= 08 septembre 2001, BARCELONA 2'
c	write(*,*)  '96= 08 septembre 2001, ST ESTEVE SESROVIRES'
c	write(*,*)  '97= 08 septembre 2001, ALELLA'
c	write(*,*)  '98= 08 septembre 2001, CASTELLON'
c	write(*,*)  '99= 08 septembre 2001, HORTONEDA'
c	write(*,*) '100= 08 septembre 2001, BARCELONA 3'
c	write(*,*) '101= 08 septembre 2001, SABADELL (CASAS)'
c	write(*,*) '102= 08 septembre 2001, ESPLUGUES DE LLOBREGAT'
c	write(*,*) '103= 08 septembre 2001, ZARAGOZA'
c	write(*,*) '104= 08 septembre 2001, CALERN'
c	write(*,*) '105= 08 septembre 2001, DAX'
c	write(*,*) '106= 08 septembre 2001, PUIMICHEL'
c	write(*,*) '107= 08 septembre 2001, CHATELLERAULT'
c
c	write(*,*) '110= 01 aout 2003, George Observatory'
c	write(*,*) '111= 01 aout 2003, Monterrey'
c
c	write(*,*) '120= 20 juillet 2002, ARICA'
c	write(*,*) '121= 20 juillet 2002, JERUSALEM'
c	write(*,*) '122= 20 juillet 2002, CUMBAYA'
c	write(*,*) '123= 20 juillet 2002, LIMA'
c	write(*,*) '124= 20 juillet 2002, MAMINA'
c	write(*,*) '125= 20 juillet 2002, PARANAL'
c	write(*,*) '126= 20 juillet 2002, LOS ARMAZONES'
c	write(*,*) '127= 20 juillet 2002, ITAJUBA'
c	write(*,*) '128= 20 juillet 2002, LA SILLA'
c	write(*,*) '129= 20 juillet 2002, EL LEONCITO'
c	write(*,*) '130= 20 juillet 2002, MERIDA'

c	write(*,*) '150= 21 aout 2002, CFHT' 
c	write(*,*) '151= 21 aout 2002, Lowell' 
c	write(*,*) '152= 21 aout 2002, Palomar' 
c	write(*,*) '153= 21 aout 2002, Lick' 

c	write(*,*) '160= 15 decembre 2002, Tomar' 
c	write(*,*) '161= 15 decembre 2002, Sabadell (Casas)' 
c	write(*,*) '162= 15 decembre 2002, Pic du Midi'
c	write(*,*) '163= 15 decembre 2002, Czarna Bialostocka'
c	write(*,*) '164= 15 decembre 2002, Vitebsk'
c	write(*,*) '165= 15 decembre 2002, Meudon' 

c	write(*,*) '166= 15 decembre 2002, Kooriyama'
c	write(*,*) '167= 15 decembre 2002, Ooe'
c	write(*,*) '168= 15 decembre 2002, Kashiwa'
c	write(*,*) '169= 15 decembre 2002, Hitachi' 
c	write(*,*) '170= 15 decembre 2002, Chichibu'
c	write(*,*) '171= 15 decembre 2002, Musashino'
c	write(*,*) '172= 15 decembre 2002, Mitaka'
c	write(*,*) '173= 15 decembre 2002, Abrera'

c	write(*,*) '200= 01 avril 2003, Pic' 
c	write(*,*) '201= 01 avril 2003, Palomar' 
c	write(*,*) '202= 01 avril 2003, La Silla' 
c	write(*,*) '203= 01 avril 2003, Tenerife (Roque Muchachos)'
c	write(*,*) '204= 01 avril 2003, New Jersey'
c	write(*,*) '205= 01 avril 2003, Max Valier Observatory' 
c	write(*,*) '206= 01 avril 2003, Nyrola Observatory '
c	write(*,*) '207= 01 avril 2003, Livermore'

c	write(*,*) '220= 14 novembre 2003, Windhoek'
c	write(*,*) '221= 14 novembre 2003, Tivoli'
c	write(*,*) '222= 14 novembre 2003, HESS'
c	write(*,*) '223= 14 novembre 2003, Hakos'
c	write(*,*) '224= 14 novembre 2003, Kleinbegin'
c	write(*,*) '225= 14 novembre 2003, Sandfontein'
c	write(*,*) '226= 14 novembre 2003, Springbok'
c	write(*,*) '227= 14 novembre 2003, Nuwerus'
c	write(*,*) '228= 14 novembre 2003, Gifberg'
c	write(*,*) '229= 14 novembre 2003, Cederberg'
c	write(*,*) '230= 14 novembre 2003, SAAO Sutherland'
c	write(*,*) '231= 14 novembre 2003, SAAO Cape Town'
c	write(*,*) '232= 14 novembre 2003, Boyden'
c	write(*,*) '240= 14 novembre 2003, Maido'
c	write(*,*) '241= 14 novembre 2003, Les Makes'
c	write(*,*) '242= 14 novembre 2003, Fournaise'
	
c	write(*,*) '243= 14 novembre 2003, Pico Veleta'
c	write(*,*) '244= 14 novembre 2003, WIRO'
c	write(*,*) '245= 14 novembre 2003, Merida (Llano del Hato)'

c	write(*,*) '260= 11 juillet 2005, La Silla'
c	write(*,*) '261= 11 juillet 2005, Paranal'
c	write(*,*) '262= 11 juillet 2005, Montevideo' 
c	write(*,*) '263= 11 juillet 2005, Bosque Alegre' 
c	write(*,*) '264= 11 juillet 2005, El Leoncito' 
c	write(*,*) '265= 11 juillet 2005, Itajuba' 
c	write(*,*) '266= 11 juillet 2005, Asuncion' 
c	write(*,*) '267= 11 juillet 2005, San Pedro de Atacama'
c	write(*,*) '268= 11 juillet 2005, Tarija'
c	write(*,*) '269= 11 juillet 2005, Marangani'
c	write(*,*) '270= 11 juillet 2005, Wykrota (B. Horizonte)'
c	write(*,*) '271= 11 juillet 2005, CEAMIG-REA (B. Horizonte)'
c	write(*,*) '272= 11 juillet 2005, Patacamaya'
c	write(*,*) '273= 11 juillet 2005, Tilomonte'
c	write(*,*) '274= 11 juillet 2005, SOAR'
c	write(*,*) ' '
c	write(*,*) '300= 10 avril 2006, VLT'
c	write(*,*) '301= 10 avril 2006, La Silla'
c	write(*,*) '310= 12 juin  2006, Mount John Obs.'
c	write(*,*) '311= 12 juin  2006, Hobart'
c	write(*,*) '312= 12 juin  2006, AAT'
c	write(*,*) '313= 12 juin  2006, Les Makes'
c	write(*,*) '314= 12 juin  2006, Gault, Blue Mountains'
c	write(*,*) '315= 12 juin  2006, Blair, Stockport Observatory'
c	write(*,*) ' '  
c	write(*,*) '320= 06 aout  2006, Pluton'
c	write(*,*) ' '
c	write(*,*) '330= 18 mars  2007, Lick (Crossley)'
c	write(*,*) '331= 18 mars  2007, Kitt Peak (Bok)'
c	write(*,*) '332= 18 mars  2007, Mt Bigelow (Kuiper)'
c	write(*,*) '333= 18 mars  2007, Tenagra'
c	write(*,*) '334= 18 mars  2007, Guanajuato'
c	write(*,*) '335= 18 mars  2007, MMT'
c	write(*,*) '336= 18 mars  2007, Pinto Valley, CA'
c	write(*,*) '337= 18 mars  2007, Hereford, AZ'
c	write(*,*) '338= 18 mars  2007, Cloudbait, CO'
c	write(*,*) '339= 18 mars  2007, Palmer Divide, CO'
c	write(*,*) '340= 18 mars  2007, Calvin-Rehoboth, NM'
c	write(*,*) '341= 18 mars  2007, Appalachian, NC'
c	write(*,*) '342= 18 mars  2007, Moore Obs, WA'
c	write(*,*) '343= 18 mars  2007, George Obs, TX'
c	write(*,*) '344= 18 mars  2007, Lowell Obs, AZ'
c	write(*,*) '345= 18 mars  2007, Palomar'
c	write(*,*) '346= 18 mars  2007, San Pedro Martir'
c	write(*,*) '347= 18 mars  2007, McDonald, TX'
c	write(*,*) '348= 18 mars  2007, Univ. of Oklahoma'
c	write(*,*) '349= 18 mars  2007, KASI, Mt Lemmon'
c
c	write(*,*) '360= 12 mai   2007, Paranal'
c	write(*,*) '361= 12 mai   2007, La Silla'
c	write(*,*) '362= 12 mai   2007, Los Armazones'
c	write(*,*) '363= 12 mai   2007, San Pedro'		
c	write(*,*) '364= 12 mai   2007, Itajuba'
c	write(*,*) '365= 12 mai   2007, Windhoek'
c	write(*,*) '366= 12 mai   2007, Hakos'
c	write(*,*) '367= 12 mai   2007, Sandfontein'
c				
c	write(*,*) '380= 14 juin  2007, Paranal'
c	write(*,*) '381= 14 juin  2007, La Silla'
c	write(*,*) '382= 14 juin  2007, SOAR'
c	write(*,*) '383= 14 juin  2007, San Pedro de Atacama'
c	write(*,*) '384= 14 juin  2007, Itajuba'	
c	write(*,*) '385= 14 juin  2007, SALT'
c	write(*,*) '386= 14 juin  2007, Hakos'
c	write(*,*) '387= 14 juin  2007, Windhoek'
c	write(*,*) '388= 14 juin  2007, Les Makes'
c	write(*,*) '389= 14 juin  2007, Mairinque'
c
c	write(*,*) '390= 14 juin  2007, Itajuba, PAR RAPPORT BARYCENTRE Pluton/Charon'	
c
c	write(*,*) '391= 09 juin  2007, SOAR'
c	write(*,*) '392= 09 juin  2007, La Silla'
c
c	write(*,*) '395= 31 juillet  2007, Mt John'
c	write(*,*) '396= 31 juillet  2007, Canopus'
c	write(*,*) '397= 31 juillet  2007, Canberra'
c
c	write(*,*) '400= 21 mai 2008, Hakos'
c	write(*,*) '401= 21 mai 2008, Windhoek'
c	write(*,*) '402= 21 mai 2008, Tivoli'
c	write(*,*) '403= 21 mai 2008, Sossusvlei'
c	write(*,*) '404= 21 mai 2008, Sprinbok'
c	write(*,*) '405= 21 mai 2008, Cederberg'
c	write(*,*) '406= 21 mai 2008, SAAO (1m)'
c	write(*,*) '407= 21 mai 2008, Cape Town'
c	write(*,*) '408= 21 mai 2008, Les Makes'
c	write(*,*) '409= 21 mai 2008, Maido'
c	write(*,*) '410= 21 mai 2008, Fournaise'
c
c	write(*,*) '420= 22 juin 2008, Pluton, Hobart'
c	write(*,*) '421= 22 juin 2008, Pluton, Blue Mountains, Sydney'
c	write(*,*) '422= 22 juin 2008, Pluton, SAAO'
c	write(*,*) '423= 22 juin 2008, Pluton, Canberra'
c	write(*,*) '424= 22 juin 2008, Pluton, Siding Spring'
c	write(*,*) '425= 22 juin 2008, Pluton, Stockport'
c	write(*,*) '426= 22 juin 2008, Pluton, Perth (Lowell Telescope)'	
c	write(*,*) '427= 22 juin 2008, Pluton, Perth, Roger Groom'
c	write(*,*) '428= 22 juin 2008, Pluton, Perth, Creg Bolt'	
c	write(*,*) '429= 22 juin 2008, Pluton, Reedy Creek, John Broughton'	
c	write(*,*) '430= 22 juin 2008, Pluton, Bankstown, Sydney, Ted Dobosz'	
c	write(*,*) '431= 22 juin 2008, Pluton, Glenlee, Steve Kerr'	
c	write(*,*) '432= 22 juin 2008, Pluton, Les Makes'	
c	write(*,*) '433= 22 juin 2008, Pluton, Sprinbok'	
c	write(*,*) '434= 22 juin 2008, Pluton, Hakos'	
c
c	write(*,*) '440= 22 juin 2008, Charon, SAAO'
c	write(*,*) '441= 22 juin 2008, Charon, Sprinbok'
c	write(*,*) '442= 22 juin 2008, Charon, Hakos'
c	write(*,*) '443= 22 juin 2008, Charon, Les Makes'
c
c	write(*,*) '450= 24 juin 2008, Pluton, CFHT'
c
c	write(*,*) '460= 30 juin 2009, 2002 MS4, VLT'
c	write(*,*) '461= 30 juin 2009, 2002 MS4, La Silla'
c	write(*,*) '462= 30 juin 2009, 2002 MS4, San Pedro'
c	write(*,*) '463= 30 juin 2009, 2002 MS4, Itajuba'
c	write(*,*) '464= 30 juin 2009, 2002 MS4, Ceamig B. Horizonte'
c
c	write(*,*) '470= 25 aout 2008, Pluton, Lick'
c	write(*,*) '471= 25 aout 2008, Pluton, Grand Rapids'
c
c	write(*,*) '480= 19 fevrier 2010, Varuna, HESS'
c	write(*,*) '481= 19 fevrier 2010, Varuna, Hakos'
c	write(*,*) '482= 19 fevrier 2010, Varuna, Natal'
c	write(*,*) '483= 19 fevrier 2010, Varuna, Fortaleza/Quixada'
c	write(*,*) '484= 19 fevrier 2010, Varuna, Sao Luis'
c	write(*,*) '485= 19 fevrier 2010, Varuna, Stellenbosch Flying Field'
c	write(*,*) '486= 19 fevrier 2010, Varuna, Recife/Camalau'
c	write(*,*) '487= 19 fevrier 2010, Varuna, Tivoli'
c	write(*,*) '488= 19 fevrier 2010, Varuna, Cecite'
c	write(*,*) '489= 19 fevrier 2010, Varuna, Florianopolis'
c
c	write(*,*) '490= 19 fevrier 2010, Varuna, Windhoek/MIT'
c	write(*,*) '491= 19 fevrier 2010, Varuna, Maceio/MIT'
c
c	write(*,*) '500= 14 fevrier 2010, Pluton, Pic'
c	write(*,*) '501= 14 fevrier 2010, Pluton, Lu'
c	write(*,*) '502= 14 fevrier 2010, Pluton, Sisteron'
c
c	write(*,*) '510= 19 mai 2010, Pluton, VLT' 
c	write(*,*) '511= 19 mai 2010, Pluton, NTT' 
c	write(*,*) '512= 19 mai 2010, Pluton, SOAR' 
c
c	write(*,*) '515= 04 juin 2010, Pluton, Mt John' 
c	write(*,*) '516= 04 juin 2010, Pluton, Hobart'
c	write(*,*) '517= 04 juin 2010, Pluton, Stu Parker, Oxford'
c	write(*,*) '518= 04 juin 2010, Pluton, Bill Allen, Vintage Lane Obs'
c	write(*,*) '519= 04 juin 2010, Pluton, Dave Gault, West Wyalong'
c	write(*,*) '520= 04 juin 2010, Pluton, Hristo Pavlov, Jugiong'
c
c	write(*,*) '530= 04 juillet 2010, Pluton, El Leoncito' 
c
c	write(*,*) '540=  17 aout 2010, Ceres, LNA' 
c	write(*,*) '541=  17 aout 2010, Ceres, CEAMIG' 
c	write(*,*) '542=  17 aout 2010, Ceres, Ponta Grossa' 
c	write(*,*) '543=  17 aout 2010, Ceres, IMPE-Sao Jose dos Campos' 
c	write(*,*) '544=  17 aout 2010, Ceres, Florianopolis' 
c
c	write(*,*) '545=  08 janvier 2011, 2003AZ84, San Pedro de Atacama' 
c	write(*,*) '546=  08 janvier 2011, 2003AZ84, TRAPPIST - La Silla' 
c
c	write(*,*) '550=  06 novembre 2010, Eris, La Silla Trappist' 
c	write(*,*) '551=  06 novembre 2010, Eris, San Pedro' 
c	write(*,*) '552=  06 novembre 2010, Eris, Leoncito' 
c	write(*,*) '553=  06 novembre 2010, Eris, Itajuba' 
c	write(*,*) '554=  06 novembre 2010, Eris, Fortaleza' 
c
c	write(*,*) '560= *05*novembre 2010, Eris, Itajuba' 
c
c	write(*,*) '570=  23 avril 2011, Makemake, NTT/Ultracam' 
c	write(*,*) '571=  23 avril 2011, Makemake, VLT/ISAAC' 
c	write(*,*) '572=  23 avril 2011, Makemake, Los Armazones' 
c	write(*,*) '573=  23 avril 2011, Makemake, San Pedro de Atacama' 
c	write(*,*) '574=  23 avril 2011, Makemake, LNA' 
c	write(*,*) '575=  23 avril 2011, Makemake, TRAPPIST' 
c
c	write(*,*) '584=  04 mai   2011, Quaoar, Las Campanas OGLE'                     
c	write(*,*) '585=  04 mai   2011, Quaoar, Cerro Tololo 4m Blanco'                     
c	write(*,*) '586=  04 mai   2011, Quaoar, Montevideo'                     
c	write(*,*) '587=  04 mai   2011, Quaoar, Brasilia'                     
c	write(*,*) '588=  04 mai   2011, Quaoar, CEAMIG'                     
c       write(*,*) '589=  04 mai   2011, Quaoar, Leoncito - AR'               
c	write(*,*) '590=  04 mai   2011, Quaoar, Armazones'                     
c	write(*,*) '591=  04 mai   2011, Quaoar, San Pedro de Atacama'  	
c	write(*,*) '592=  04 mai   2011, Quaoar, Rivera, Uruguay' 		
c	write(*,*) '593=  04 mai   2011, Quaoar, Salto, Uruguay' 		
c	write(*,*) '594=  04 mai   2011, Quaoar, San Martina (PUC Santiago)' 	
c	write(*,*) '595=  04 mai   2011, Quaoar, Itajuba' 			
c	write(*,*) '596=  04 mai   2011, Quaoar, Ponta Grossa'  		
c	write(*,*) '597=  04 mai   2011, Quaoar, Santa Catarina'		
c	write(*,*) '598=  04 mai   2011, Quaoar, Sao Jose dos Campos INPE - BR' 
c	write(*,*) '599=  04 mai   2011, Quaoar, TRAPPIST La Silla - CH' 	
c
c	write(*,*) '600= 04 juin 2011, Pluton, Trappist' 
c	write(*,*) '601= 04 juin 2011, Pluton, San Pedro' 
c	write(*,*) '602= 04 juin 2011, Pluton, Itajuba' 
c	write(*,*) '603= 04 juin 2011, Pluton, San Martina' 
c
c	write(*,*) '610= 04 juin 2011, Charon, Trappist' 
c	write(*,*) '611= 04 juin 2011, Charon, San Pedro' 
c	write(*,*) '612= 04 juin 2011, Charon, Itajuba' 
c	write(*,*) '613= 04 juin 2011, Charon, San Martina' 
c
c*---
c	write(*,*) '620= 23 juin 2011, Pluton, Kaui, Kekaha'
c	write(*,*) '621= 23 juin 2011, Pluton, Maui, Lahaina'
c	write(*,*) '622= 23 juin 2011, Pluton, CFHT'
c	write(*,*) '623= 23 juin 2011, Pluton, Kwajalein'
c	write(*,*) '624= 23 juin 2011, Pluton, Nauru'
c	write(*,*) '625= 23 juin 2011, Pluton, Majuro'
c	write(*,*) '629= 23 juin 2011, Pluton, xxx'
c
c	write(*,*) '630= 23 juin 2011, Charon, Kaui, Kekaha'
c	write(*,*) '631= 23 juin 2011, Charon, Maui, Lahaina'
c	write(*,*) '632= 23 juin 2011, Charon, CFHT'
c	write(*,*) '633= 23 juin 2011, Charon, Kwajalein'
c	write(*,*) '634= 23 juin 2011, Charon, Nauru'
c	write(*,*) '635= 23 juin 2011, Charon, Majuro'
c	write(*,*) '639= 23 juin 2011, Charon, xxx'
c
c	write(*,*) '640= 27 juin 2011, Pluton, Kaui, Kekaha'
c	write(*,*) '641= 27 juin 2011, Pluton, Maui, Lahaina'
c	write(*,*) '642= 23 juin 2011, Pluton, CFHT'
c	write(*,*) '643= 23 juin 2011, Pluton, Kwajalein'
c	write(*,*) '644= 23 juin 2011, Pluton, Nauru'
c	write(*,*) '645= 23 juin 2011, Pluton, Majuro'
c	write(*,*) '649= 27 juin 2011, Pluton, xxx'
c
c	write(*,*) '650= 27 juin 2011, Hydra, Darwin'
c	write(*,*) '659= 27 juin 2011, Hydra, xxx'
c
c	write(*,*) '660= 03 fevrier 2012, 2003 AZ84, Alicante/Spain '
c	write(*,*) '661= 03 fevrier 2012, 2003 AZ84, Mt Abu/India '
c	write(*,*) '662= 03 fevrier 2012, 2003 AZ84, Weizmann/Israel '
c	write(*,*) '663= 03 fevrier 2012, 2003 AZ84, IUCAA Girawali/India '
c	write(*,*) '664= 03 fevrier 2012, 2003 AZ84, Liverpool/Canary '
c	write(*,*) '679= 03 fevrier 2012, 2003 AZ84,  '
c
c	write(*,*) '680= 17 fevrier 2012, Quaoar, Cote d Azur/FR '
c	write(*,*) '681= 17 fevrier 2012, Quaoar, TAROT Grasse/FR '
c	write(*,*) '682= 17 fevrier 2012, Quaoar, JLX Valensole/FR '
c	write(*,*) '683= 17 fevrier 2012, Quaoar, Sposseti Belinzona/CH'
c
c	write(*,*) '685= 03 December 2015, 2002 VE95, Samford Valley, Qld, Australia'
c	write(*,*) '686= 03 December 2015, 2002 VE95, Flynn, ACT, Australia'
c	write(*,*) '687= 03 December 2015, 2002 VE95, Murrumbateman, NSW, Australia'
c
c	write(*,*) '688= July 09, 2013 - Quaoar-Mérida (CIDA), Venezuela'
c
c	write(*,*) '689= January 13, 2013 - Sedna-Cairns (Tim), Australia'
c	write(*,*) '690= January 13, 2013 - Sedna-Cairns (Joe), Australia'
c
c	write(*,*) '691= August 29, 2013 - Eris/Alice Springs, Australia'
c	write(*,*) '692= August 29, 2013 - Eris/Samford Valley, Australia'

c	write(*,*) '693= November 07, 2014 - 2003VS2/Bosque Alegre/ARG'
c	write(*,*) '694= November 07, 2014 - 2003VS2/Montevidéo/URU'
c	write(*,*) '695= November 07, 2014 - 2003VS2/Rosário/ARG'
c	write(*,*) '696= November 07, 2014 - 2003VS2/Santa Rosa/ARG'
c	write(*,*) '697= November 07, 2014 - 2003VS2/NTT/La Silla/CHI'
c	write(*,*) '698= November 07, 2014 - 2003VS2/PROMPT/Tololo/CHI'

	write(*,*) '700= November 07, 2014 - 2003VS2/ Bosque Alegre'
	write(*,*) '701= November 07, 2014 - 2003VS2/ Busoniano'
	write(*,*) '702= November 07, 2014 - 2003VS2/ Los Molinos'
	write(*,*) '703= November 07, 2014 - 2003VS2/ Santa Rosa'

*---

	write(*,*) 'code, correction Dt (sec) et iprint'

c	read(*,*)  iobs

	read(*,*)  iobs, dt, iprint
	dt= dt/3600.d0
c
c Dans le cas de l'observation du 2 aout 1988 et du 11 juillet 1990, il y a
c ambiguite sur le trajet de l'etoile. Option sud: l'etoile est passee au sud
c du centre de Neptune, option nord: elle est passee au nord.
c
	if(iobs.eq.11.or.iobs.eq.30) then
		write(*,*)  '-1= option sud'
		write(*,*)  ' 1= option nord'
		read(*,*)  iopns
	endif
c
c GMasses de Neptune et Titan pour deviation relativiste.
c
c Unite: km**3/s**2, on passe a la masse en grammes en multipliant les nombres
c        suivants par 1.499 10**22
c
	if(iobs.ge.0.and.iobs.le.17)   GM= 6.87d+06		! Neptune
	if(iobs.ge.18.and.iobs.le.29)  GM= 8.977d+03	! Titan
	if(iobs.ge.30.and.iobs.le.32)  GM= 6.87d+06		! Neptune
	if(iobs.eq.33)                 GM= 1.2669d+08 	! Jupiter
	if(iobs.ge.34.and.iobs.le.38)  GM= 1.4279d+03 	! Triton
	if(iobs.eq.39)                 GM= 1.2669d+08 	! Jupiter
	if(iobs.ge.40.and.iobs.le.45)  GM= 1.2669d+08 	! Jupiter
	if(iobs.ge.50.and.iobs.le.50)  GM= 37931272. 	! Saturne
	if(iobs.ge.51.and.iobs.le.0111)   GM= 235.3		! Titania (http://ssd.jpl.nasa.gov/?sat_phys_par)
	if(iobs.ge.0120.and.iobs.le.0153) GM= 870.7		! Pluton (memo B. Jacobson mail 17/5/05)
	if(iobs.ge.0160.and.iobs.le.0173) GM=  41.8		! Tethys 
	if(iobs.ge.0200.and.iobs.le.0207) GM=  1.2669d+08 ! Jupiter 
	if(iobs.ge.0220.and.iobs.le.0245) GM=  8.9782d+03 ! Titan (http://www.imcce.fr/solarsys/projet/chap2/2_17.html)
	if(iobs.ge.0260.and.iobs.le.0274) GM=  108.0	! Charon (JPL site)
	if(iobs.ge.0300.and.iobs.le.0315) GM=  870.7	! Pluton
	if(iobs.eq.0320)                  GM=  870.7	! Pluton
	if(iobs.ge.0330.and.iobs.le.0397) GM=  870.7	! Pluton
	if(iobs.ge.0400.and.iobs.le.0410) GM=  1427.9	! Triton
	if(iobs.ge.0420.and.iobs.le.0434) GM=  870.7	! Pluton
	if(iobs.ge.0440.and.iobs.le.0443) GM=  108.0	! Charon
	if(iobs.ge.0450.and.iobs.le.0450) GM=  870.7	! Pluton
	if(iobs.ge.0460.and.iobs.le.0464) GM=   23.0	! 2002 MS4 valeur approx. deduite de Pluton (R= 350 km vs. 1180 km)
	if(iobs.ge.0470.and.iobs.le.0471) GM=  870.7	! Pluton
	if(iobs.ge.0480.and.iobs.le.0491) GM=    4.4	! Varuna, valeur approx R= 250 km, rho= 1000 kg/m**3
	if(iobs.ge.0500.and.iobs.le.0530) GM=  870.7	! Pluton
	if(iobs.ge.0540.and.iobs.le.0544) GM=   62.9	! Ceres (Carry et al. 2008)
	if(iobs.ge.0545.and.iobs.le.0546) GM=   23.0	! 2003AZ84 (usada MS4, confirmer!!!)
	if(iobs.ge.0550.and.iobs.le.0560) GM= 1114.1	! Eris, Brown & Schaller Science 316, 1585 (2007)
	if(iobs.ge.0570.and.iobs.le.0575) GM=   230.	! Makemake, APPROXIMATIF! Suppose meme densite Pluton, R= 750 km, vs 1170 pour Pluton
	if(iobs.ge.0584.and.iobs.le.0599) GM=  106.8	! Quaoar, Fraser &  Brown ApJ 714, 1547-1550 (2010)
	if(iobs.ge.0600.and.iobs.le.0603) GM=  870.7	! Pluton
	if(iobs.ge.0610.and.iobs.le.0613) GM=  108.0	! Charon
	if(iobs.ge.0620.and.iobs.le.0629) GM=  870.7	! Pluton
	if(iobs.ge.0630.and.iobs.le.0639) GM=  108.0	! Charon
	if(iobs.ge.0640.and.iobs.le.0649) GM=  870.7	! Pluton
	if(iobs.ge.0650.and.iobs.le.0659) GM=  0.021	! Hydra (Tholen et al. AJ 2008)
	if(iobs.ge.0660.and.iobs.le.0679) GM=   23.0	! 2003AZ84 (usada MS4, confirmer!!!)
	if(iobs.ge.0680.and.iobs.le.0683) GM=  106.8	! Quaoar, Fraser &  Brown ApJ 714, 1547-1550 (2010)
	if(iobs.ge.0685.and.iobs.le.0687) GM=  2.0	! 2002 VE95, aproximado pela massa de Chariklo
	if(iobs.ge.0688.and.iobs.le.0688) GM=  106.8	! Quaoar, Fraser &  Brown ApJ 714, 1547-1550 (2010)
	if(iobs.ge.0689.and.iobs.le.0690) GM=  106.8	! Sedna, aproximado pela massa de Quaoar que tem raio parecido.
	if(iobs.ge.0691.and.iobs.le.0692) GM= 1114.1	! Eris, Brown & Schaller Science 316, 1585 (2007) (como em 0560).
	if(iobs.ge.0693.and.iobs.le.0698) GM=  4.0	! 2003 VS2, aproximado pelo dobro da massa de Chariklo.
	if(iobs.ge.0700.and.iobs.le.0703) GM=  5.0	! 2003 VS2, aprox. by R=260 km rho= 1000 kg/m**3
	write(*,*)  'Masse du corps (g): ', 1.499d+022*GM
c
c On precise si la masse de Triton est consideree ou pas pour la projection
c dans le plan des anneaux.
c
c	write(*,*)  '0= Triton n''est pas pris en compte'
c	write(*,*)  '1= Triton est pris en compte'
c	read(*,*)  itriton
	itriton = 0
c
	if (imode.eq.1) then 
	 write(*,*)  'Heure T.U. ?'
	 read(5,*) ih,im,sec
	 heure= (dfloat(ih)*3600.d+00 + dfloat(im)*60.d+00 + sec)/3600.d+00
	 sec= sec + 1.d+00
	 heurf= (dfloat(ih)*3600.d+00 + dfloat(im)*60.d+00 + sec)/3600.d+00
	endif

 2	continue
c
c Fin des parametres d'entree
c
c
c ..................................................................
c
c Position de l'etoile vue du centre de la Terre:
c
	call geocentrique (iobs,iopns,heure+dt, ksig, etag, dis )
	call geocentrique (iobs,iopns,heurf+dt, ksigf,etagf,disf)
c
c On effectue les corrections sur ksi et eta
c
	ksig =  ksig  + dksi0
	etag =  etag  + deta0
	ksigf=  ksigf + dksi0
	etagf=  etagf + deta0
c
c Effet de la parallaxe diurne
c
	call parallaxe (iobs,heure+dt,   alpha_star,delta_star,ksip, etap)
	call parallaxe (iobs,heurf+dt,   alpha_star,delta_star,ksipf,etapf)
	cos_delta_star= dcos(delta_star)
	alpha_pla= alpha_star - (ksig-dksi0)/(dis*cos_delta_star)	! NB. il ne faut pas appliquer la correction
	delta_pla= delta_star - (etag-deta0)/dis					! dksi0, deta0 pour avoir l'ephemeride de la planete

	alpha_pla_loc= alpha_pla  - ksip/(dis*cos_delta_star)		! coordonnees topocentriques planete
	delta_pla_loc= delta_pla  - etap/dis

	alpha_star    = alpha_star/(rad*15.d0)			! heure
	delta_star    = delta_star/rad					! degre
	alpha_pla     = alpha_pla/(rad*15.d0)			! heure
	delta_pla     = delta_pla/rad					! degre
	alpha_pla_loc = alpha_pla_loc/(rad*15.d0)		! heure
	delta_pla_loc = delta_pla_loc/rad				! degre

	call dms (alpha_star, iahstar,iamstar,asecstar)
	call dms (delta_star, iddstar,idmstar,dsecstar)
	write (*,'(A,I2.2,1x,I2.2,1x,f11.8,1x,I3,1x,I2.2,1x,f11.8)') 'coordonnees  de l''etoile:                 ', 
     *           iahstar,iamstar,asecstar, iddstar,idmstar,dsecstar
	call dms (alpha_pla, iahpla,iampla,asecpla)
	call dms (delta_pla, iddpla,idmpla,dsecpla)
	call dms (alpha_pla_loc, iahpla_loc,iampla_loc,asecpla_loc)
	call dms (delta_pla_loc, iddpla_loc,idmpla_loc,dsecpla_loc)
	write (*,'(A,I2.2,1x,I2.2,1x,f11.8,1x,I3,1x,I2.2,1x,f11.8)') 'coordonnees geocentriques  de la planete: ', 
     *           iahpla,iampla,asecpla, iddpla,idmpla,dsecpla
	write (*,'(A,I2.2,1x,I2.2,1x,f11.8,1x,I3,1x,I2.2,1x,f11.8)') 'coordonnees topocentriques de la planete: ', 
     *           iahpla_loc,iampla_loc,asecpla_loc, iddpla_loc,idmpla_loc,dsecpla_loc
	
c
c Calcul des coefficients aa,bb,cc,dd pour projeter dans le plan des anneaux
c
	call polpo (iobs,itriton,dis, B,P,BP,aa,bb,cc,dd,aap,bbp,ccp,ddp)
c
c Position dans l'ombre
c
	ksi=  ksig  + ksip
	eta=  etag  + etap
	ksif= ksigf + ksipf
	etaf= etagf + etapf
	vksig= ksigf-ksig
	vetag= etagf-etag
	vksi= ksif-ksi
	veta= etaf-eta
	vciele= dsqrt( vksi*vksi + veta*veta )
	vcielr= dsqrt(ksif*ksif+etaf*etaf) - dsqrt(ksi*ksi+eta*eta)
c
c Correction de la deflection relativiste
c
	relat= (4.d0*GM*dis)/c**2
	ray2=   ksi*ksi + eta*eta
	dksir=  (ksi* relat)/ray2
	detar=  (eta* relat)/ray2
	dksirf= (ksif*relat)/ray2
	detarf= (etaf*relat)/ray2
c
c Position dans le plan du ciel (correction relativiste)
c
	ksir=  ksi  + dksir
	etar=  eta  + detar
	ksirf= ksif + dksirf
	etarf= etaf + detarf
	vksir= ksirf-ksir
	vetar= etarf-etar
	vcieler= dsqrt(vksir*vksir + vetar*vetar)
	vcielrr= dsqrt(ksirf*ksirf + etarf*etarf)-dsqrt(ksir*ksir + etar*etar)
c
c Vitesse perpendiculaire dans le plan du ciel
c
	xper=   ksir*dcos(P) -  etar*dsin(P)
	yper=   ksir*dsin(P) +  etar*dcos(P)
	vxper= vksir*dcos(P) - vetar*dsin(P)
	vyper= vksir*dsin(P) + vetar*dcos(P)
	coeff= (xper*dsin(B))**2 + yper**2
	anx= xper*dsin(B)*dsin(B)/coeff
	any= yper                /coeff
	an = dsqrt( anx*anx + any*any )
	anx= anx/an
	any= any/an
	vper= vxper*anx + vyper*any
c
c elevation planetocentrique du point dans le plan du ciel
c
	sinphi= dcos(B)*( dsin(P)*ksir + dcos(P)*etar )
	phi= dasin( sinphi/dsqrt(ksir*ksir+etar*etar) )
c
c Position dans le plan des anneaux, voir le sous-programme "polpo" pour 
c la definition de x et y
c
	x=  ksir*aa  + etar*bb
	y=  ksir*cc  + etar*dd
	xf= ksirf*aa + etarf*bb
	yf= ksirf*cc + etarf*dd
c
c Position dans le plan des anneaux, par rapport au noeud ascendant sur le plan equatorial terrestre,
c voir le sous-programme "polpo" pour la definition de x et y
c
	xp= x*aap + y*bbp
	yp= x*ccp + y*ddp
	rp= dsqrt( xp*xp + yp*yp )

	if ((xp.eq.0.).and.(yp.ge.0.)) then
		along= pi/2.d0
	endif
	if ((xp.eq.0.).and.(yp.le.0.)) then
		along= -pi/2.d0
	endif
	if (xp.ne.0.) then
		along= datan(yp/xp)
		if(xp.le.0) along= along+pi
	endif

	along= along*(180./pi)
	xpf= xf*aap + yf*bbp
	ypf= xf*ccp + yf*ddp
	rpf= dsqrt( xpf*xpf + ypf*ypf )
	vxp= xpf-xp
	vyp= ypf-yp
	vequa= dsqrt( vxp*vxp + vyp*vyp )
	vrad= rpf-rp
c
c trace de l'ombre de la planete dans son plan equatorial
c
c	rp= 25250.
c	do i= 0, 72
c		angle= dfloat(i)*5.*(pi/180.)
c		xx= rp*cos(angle)
c		yy= rp*sin(angle)
c		xxx= xx*aa + yy*bb
c		yyy= xx*cc + yy*dd
c		xxxp= xxx*aap + yyy*bbp
c		yyyp= xxx*ccp + yyy*ddp
c	enddo
c
c
c
c calcul de la position des arcs
c
	tlum= dis/(3600.d+00*c)
c	heuredis= heure - tlum
c	call jjmmaa (iobs, jj,mm,annee)
c	call ajulien (jj,mm,annee, ajd,hs)
c	write(*,*)  'Jour julien de la date ', ajd
c 	call arcpos (ajd,heuredis, alib1,alib2,ega1,ega2,fra1,fra2)
c	call arcposciel (aa,bb,cc,dd,aap,bbp,ccp,ddp,alib1, aksilib1, etalib1)
c	call arcposciel (aa,bb,cc,dd,aap,bbp,ccp,ddp,alib2, aksilib2, etalib2)
c	call arcposciel (aa,bb,cc,dd,aap,bbp,ccp,ddp, ega1, aksiega1, etaega1)
c	call arcposciel (aa,bb,cc,dd,aap,bbp,ccp,ddp, ega2, aksiega2, etaega2)
c	call arcposciel (aa,bb,cc,dd,aap,bbp,ccp,ddp, fra1, aksifra1, etafra1)
c	call arcposciel (aa,bb,cc,dd,aap,bbp,ccp,ddp, fra2, aksifra2, etafra2)
c
c
	if (imode.eq.1) then
	 write(*,*) 'Position geocentrique:  ', ksig, etag
	 write(*,*) 'Correction de parallaxe:', ksip, etap
	 write(*,'(A,2(f8.4,1x))') 'Correction relativiste sur ksi et eta:', dksir, detar
	 write(*,*) 'Position dans l''ombre (ksi, eta, dist), sans et avec correc. relativiste:'
	 write(*,*)  ksi,   eta, dsqrt( ksi*ksi   + eta*eta   )
	 write(*,*)  ksir, etar, dsqrt( ksir*ksir + etar*etar )
	endif

	if (iprint.eq.1) then
		write(20,'(F13.4,3x,2(f17.6,2x))') heure*3600.d0, ksi, eta			! pour decalage
c		write(20,*) iobs, heure*3600.d0, ksi, eta			! pour decalage
		write(21,'(3(f12.3,2x),2x,f9.3)') sngl(ksi),sngl(eta), sngl(dsqrt(ksi*ksi+eta*eta)), heure*3600.d0 
c     *		sngl(vcieler), sngl(vcielrr), sngl(phi/rad)	! avec correction relativiste
		write(22,*) sngl(ksig), sngl(+etag), sngl(dsqrt(ksig**2+etag**2)), sngl(heure*3600.d0)	! geocentrique pour utilisation dans planisphere
c		write(22,*) sngl(ksig), sngl(-etag), sngl(heure)	! pour utilisation dans "map"
c		write(23,*) iobs,sngl(heure*3600.d0),sngl(ksi),sngl(eta),
c    *			    sngl(dsqrt(ksi*ksi + eta*eta))
	    if (imode.eq.3) then
		 write (24,*) sngl(dsqrt(ksi*ksi + eta*eta)), sngl(flux_obs(i_t_obs)), sngl(heure*3600.d0)
		 write (25,*) sngl(rp), sngl(flux_obs(i_t_obs))
		 write (26,*) sngl(xp), sngl(yp), sngl(flux_obs(i_t_obs))
		endif
		call dms (heure, ih_tu,im_tu,sec_tu)
		write (27,'(A,I2,A,I2,A,f5.2,1X,f6.2,1x,f6.2)')	! ephem planete /t a etoile (pour tracking table NACO)
     *  ':', ih_tu, ':', im_tu, ':', sngl(sec_tu), -sngl(ksi/(dis*arcsec)), -sngl(eta/(dis*arcsec))
		write (28,*) heure*3600.d0, sngl(1.d0)
	endif

	if (imode.eq.1) then
c	write(40,*) sngl(heure),
c     *		    sngl(ksig), sngl(vksig), sngl(etag), sngl(vetag)
c	write(40,*) sngl(heure),
c     *		    sngl(ksi ), sngl(vksi ), sngl(eta ), sngl(veta )
c	write(41,*) sngl(heure), sngl(ksip), sngl(etap)

	 write(*,*)  'Vitesses geocentrique dans le plan du ciel (km/sec):'
	 write(*,*)  'en ksi, eta, module, radiale, ss et avec corr. relat.'
	 write(*,'(3(f8.4,2x))')  vksig,  vetag, sngl(dsqrt(vksig**2 + vetag**2))

	 write(*,*)  'Vitesses dans le plan du ciel (km/sec):'
	 write(*,*)  'en ksi, eta, module, radiale, ss et avec corr. relat.'
	 write(*,'(4(f8.4,2x))')  vksi,  veta,  vciele,  vcielr
	 write(*,'(4(f8.4,2x))')  vksir, vetar, vcieler, vcielrr
	 write(*,*)  'Vitesses: module et perpend. a l''anneau equatorial local'
	 write(*,'(2(f8.4,2x))')  vcieler, vper
	 write(*,'(A,3(f9.4,2x))')  'P, B, B'':', P/rad, B/rad, BP/rad
	 write(*,*) 'distance, km/arcsec:', dis, dis*arcsec
	 
	 write(*,*)  '1= continue'

	 read(*,*)  iop
	 write(*,*)  'Elevation planetocentrique du point du plan du ciel:'
	 write(*,*)  '(planete supposee spherique)'
	 write(*,*)  phi/rad
	 write(*,*)  'Position dans le plan des anneaux:'
	 write(*,*)  x,y
	 write(*,*)  'Position dans le plan des anneaux, equateur 1950.?'
	 write(*,*)  xp, yp
	 write(*,*)  'Rayon et longitude (degre) dans le plan des anneaux'
	 write(*,*)  rp, along
	 write(*,*)  'Vitesses dans le plan equatorial (module et radiale)'
	 write(*,*)  vequa, vrad
	 write(*,*)  'Temps lumiere a la planete et distance:'
	 write(*,*)  tlum, dis 
c	write(*,*)  'Longitudes de Liberte, Egalite et Fraternite a cet instant,'
c	write(*,*)  'et positions dans le plan du ciel'
c	write(*,*)  alib1, alib2
c	write(*,*)   ega1,  ega2
c	write(*,*)   fra1,  fra2
c	write(*,*)   aksilib1, etalib1, aksilib2, etalib2
c	write(*,*)   aksiega1, etaega1, aksiega2, etaega2
c	write(*,*)   aksifra1, etafra1, aksifra2, etafra2
c
c	write(23,*) aksilib1, etalib1
c	write(23,*) aksilib2, etalib2
c	write(23,*) aksiega1, etaega1
c	write(23,*) aksiega2, etaega2
c	write(23,*) aksifra1, etafra1
c	write(23,*) aksifra2, etafra2
c
	 write(*,*)  '1= recommencer'
	 read(*,*)  iop
	 if(iop.eq.1) go to 1
	endif

	if (imode.eq.2) then
c	 call latitude (P,B,Hkm,dis,heure+dt,ksi,eta,vksi,veta)
	 heure= heure + pas_heure/3600.d0
	 heurf= heure + 1.d0/3600.d0
	 if (heure.le.heure_1) go to 2
	endif

	if (imode.eq.3) then
	 i_t_obs= i_t_obs + 1
	 heure= t_obs(i_t_obs)/3600.d0
	 heurf= heure + 1.d0/3600.d0
	 if (i_t_obs.le.npt) go to 2
	endif	
	
c
c Pour le papier anneaux Pluton: ephemeride par rapport au barycentre Pluton/Charon, voir rapport. A neutraliser
c
	if(iobs.ge.0300.and.iobs.le.0301) write (*,*) 'ATTENTION: verifier si ephemeride barycentrique utilisee pour 10 avril 2006!'
	if(iobs.ge.0390.and.iobs.le.0390) write (*,*) 'ATTENTION: verifier si ephemeride barycentrique utilisee pour 14 juin  2007!'

	stop
	end
c
c
c 				FIN DE MAIN
c
c *************************************************************************
c
c Calcule l'effet de parallaxe pour differentes observations. Les parametres
c de chaque observation a ete mis en memoire dans le sous-programme donnees
c
c	Entrees: iobs, code de l'observation.
c		 heure T.U. decimale.
c
c	Sorties: ksip, correction de parallaxe vers l'est, en kilometres.
c		 etap,     "            "             "          "
c
c
	subroutine parallaxe (iobs,heure,   alpha_star,delta_star,ksip,etap)
c
c
c
	implicit real*8 (a-h,o-z)
	real*8 ksip
c
c A partir du 8 septembre 2001 on calcule directement l'effet de
c parallaxe
c
	if (iobs.ge.51) then
		call new_parallaxe (iobs,heure,   alpha_star,delta_star,ksip,etap)
		return
	endif
c
	call donneepar (iobs,      az,bz,cz,dz,ez)
c
	ksip=      az*dsin( bz + heure*cz )
	etap= dz + ez*dcos( bz + heure*cz )

	theta= 0.d0
c
c Dans le cas des observations de Titan du 18 juillet 1997, les coordonnees
c geocentriques sont calculees en 1950.0 alors que les corrections de parallaxe
c sont toujours calculees a l'epoque. Pour etre coherent, il faut donc calculer
c les corrections de parallaxe en 1950.0 aussi, d'ou le calcul ci-dessous (voir
c le classeur Titan du 25 juin 1998).
c
	if (iobs.ge.18.and.iobs.le.29) then
		theta = -4.0742546d-3
	endif
c
c Dans le cas des observations de Triton du 18 juillet 1997, les coordonnees
c geocentriques sont calculees en J2000, alors que les corrections de parallaxe
c sont toujours calculees a l'epoque. Pour etre coherent, il faut donc calculer
c les corrections de parallaxe en J2000 aussi, d'ou le calcul ci-dessous (voir
c le classeur Triton du 26 mai 1998)
c
	if (iobs.ge.34.and.iobs.le.38) then
		theta = 2.19d-4
	endif
c
c Observation de Jupiter du 10 octobre 1999: meme chose (voir cahier Jupiter
c du 18 janvier 2000)
c
	if (iobs.ge.40.and.iobs.le.45) then
		theta = -1.129d-5
	endif

	x= ksip
	y= etap
	ksip= x*dcos(theta) - y*dsin(theta)
	etap= x*dsin(theta) + y*dcos(theta)

	return
	end
c
c
c				FIN DE PARALLAXE
c
c **************************************************************************
c
c Sous-programme qui contient les donnees relatives a l'effet de parallaxe
c
	subroutine donneepar (iobs,      az,bz,cz,dz,ez)
c
c	Entrees: iobs, code de l'observation
c
c	Sorties: az,bz,cz,dz,ez, coefficients pour le calcul de parallaxe
c
c
c
	implicit real*8 (a-h,o-z)
c
c 15 juin 1983, HOBART
c
	if(iobs.eq.0) then
		az=   4683.586745791447
		bz=  -3.798880597845694
		cz=  0.2625161707892094
		dz=  -3996.740405424362
		ez=   1766.471890391954
	endif
c
c 15 juin 1983, CFHT
c
	if(iobs.eq.1) then
		az=   6006.373323781802
		bz=  -2.802360843195703
		cz=  0.2625161707892094
		dz=   1992.163832072317
		ez=   2265.376963382829
	endif
c
c 22 juillet 1984, ESO
c
	if(iobs.eq.2) then
		az=   5571.106693606581
		bz= -0.6981503476918034
		cz=  0.2625161707892094
		dz=  -2869.345461291702
		ez=   2108.025356656989
	endif
c
c 22 juillet 1984, CTIO
c
	if(iobs.eq.3) then
		az=   5520.999756935743
		bz= -0.6996086430038968
		cz=  0.2625161707892094
		dz=  -2950.397622968390
		ez=   2089.065624083971
	endif
c
c 20 aout 1985, ESO
c
	if(iobs.eq.4) then
		az=   5571.106693606581
		bz= -0.2360265527815983
		cz=  0.2625161707892094
		dz=  -2867.947660778468
		ez=   2114.160490447605
	endif
c
c 20 aout 1985, CTIO
c
	if(iobs.eq.5) then
		az=   5520.999756935743
		bz= -0.2374848480936917
		cz=  0.2625161707892094
		dz=  -2948.960337926466
		ez=   2095.145578037398
	endif
c
c 20 aout 1985, CFHT
c
	if(iobs.eq.6) then
		az=   6006.373323781802
		bz=  -1.715029469229418
		cz=  0.2625161707892094
		dz=   1990.120360616433
		ez=   2279.338355984227
	endif
c
c 23 avril 1986, CFHT
c
	if(iobs.eq.7) then
		az=   6006.373323781802
		bz=  -3.855935732547265
		cz=  0.2625161707892094
		dz=   1991.256881551624
		ez=   2271.585784634633
	endif
c
c 23 aout 1986, CFHT
c
	if(iobs.eq.8) then
		az=   6006.373323781802
		bz=  -1.708839942982291
		cz=  0.2625161707892094
		dz=   1989.828943910858
		ez=   2281.321247361381
	endif
c
c 22 juin 1887, ESO
c
	if(iobs.eq.9) then
		az=   5571.106693606581
		bz=  -1.369004637301931
		cz=  0.2625161707892094
		dz=  -2869.685561720135
		ez=   2106.529454479794
	endif
c
c 9 juillet 1987, CFHT
c
	if(iobs.eq.10) then
		az=   6006.373323781802
		bz=  -2.546842991065560
		cz=  0.2625161707892094
		dz=   1990.965953511567
		ez=   2273.573236666826
	endif
c
c 2 aout 1988, ESO
c
	if(iobs.eq.11) then
		az=   5571.106693606581
		bz= -0.6729202272258361
		cz=  0.2625161707892094
		dz=  -2869.619655039628
		ez=   2106.819436109178
	endif
c
c 25 aout 1985, CFHT
c
	if(iobs.eq.12) then
		az=   6006.373323781802
		bz=  -1.748638080543166
		cz=  0.2625161707892094
		dz=   1990.766786209334
		ez=   2274.932663953740
	endif
c
c 12 septembre 1988, PIC
c
	if(iobs.eq.13) then
		az=   4678.867740969617
		bz=   1.279787389994279
		cz=  0.2625161707892094
		dz=   4001.553217765268
		ez=   1773.686560428301
	endif
c
c 12 septembre 1988, OHP
c
	if(iobs.eq.14) then
		az=   4601.414379116247
		bz=   1.377029863710187
		cz=  0.2625161707892094
		dz=   4074.281557135524
		ez=   1744.325186141871
	endif
c
c 7 juillet 1989, ESO
c
	if(iobs.eq.15) then
		az=   5571.106693606581
		bz=  -7.461770363133664
		cz=  0.2625161707892094
		dz=  -2872.750207790759
		ez=   2092.993600786245
	endif
c
c 7 juillet 1989, TEIDE
c
	if(iobs.eq.16) then
		az=   5620.702834158639
		bz=  -6.515347613892853
		cz=  0.2625161707892094
		dz=   2784.992548057436
		ez=   2111.626237084214
	endif
c
c 7 juillet 1989, PIC
c
	if(iobs.eq.17) then
		az=   4678.867740969617
		bz=  -6.224787163516246
		cz=  0.2625161707892094
		dz=   4007.540172218352
		ez=   1757.790826733384
	endif
c
c 3 juillet 1989, CATANIA
c
	if(iobs.eq.18) then
		az=   5054.813255330191
		bz=  0.2537239135724793
		cz=  0.2625161707892094
		dz=   3586.633478371783
		ez=   1926.641523487401
	endif
c
c 3 juillet 1989, MEUDON
c
	if(iobs.eq.19) then
		az=   4208.899448923581
		bz=  3.1245339379718740E-02
		cz=  0.2625161707892094
		dz=   4415.884080291165
		ez=   1604.221568013142
	endif
c
c 3 juillet 1989, ODESSA
c
	if(iobs.eq.20) then
		az=    4406.5482043839
		bz=   0.52083707540124
		cz=   0.26251617078921
		dz=    4248.9259022471
		ez=    1679.4831145729
	endif
c
c 3 juillet 1989, PIC
c
	if(iobs.eq.21) then
		az=   4678.867740969617
		bz= -5.2131342533995317E-03
		cz=  0.2625161707892094
		dz=   3997.882847230498
		ez=   1783.349931503832
	endif
c
c 3 juillet 1989, EIN HAROD
c
	if(iobs.eq.22) then
		az=   5380.946725810052
		bz=  0.6100619691879882
		cz=  0.2625161707892094
		dz=   3155.298661488393
		ez=   2050.947260353740
	endif
c
c 3 juillet 1989, MANLEY
c
	if(iobs.eq.23) then
		az=   3824.901176970752
		bz= -5.5475707828749554E-02
		cz=  0.2625161707892094
		dz=   4703.002520098755
		ez=   1457.860668347542
	endif
c
c 3 juillet 1989, WISE
c
	if(iobs.eq.24) then
		az=   5495.682447349881
		bz=  0.5990373060795574
		cz=  0.2625161707892094
		dz=   2984.283221205659
		ez=   2094.678768162228
	endif
c
c 3 juillet 1989, VATICAN
c
	if(iobs.eq.25) then
		az=   4766.134735342083
		bz=  0.2131159196427445
		cz=  0.2625161707892094
		dz=   3906.036071542670
		ez=   1816.611736934654
	endif
c
c 3 juillet 1989 HATFIELD
c
	if(iobs.eq.26) then
		az=   3954.766628087630
		bz= -9.2995278976585327E-03
		cz=  0.2625161707892094
		dz=   4610.897842343304
		ez=   1507.358870941745
	endif
c
c 3 juillet 1989 NEZEL
c
	if(iobs.eq.27) then
		az=   3848.754700887263
		bz=  0.1525917593634230
		cz=  0.2625161707892094
		dz=   4686.387854140953
		ez=   1466.952436398634
	endif
c
c 3 juillet 1989 DENZAU
c
	if(iobs.eq.28) then
		az=   3986.580006445167
		bz=  0.1157125159364698
		cz=  0.2625161707892094
		dz=   4587.612647251495
		ez=   1519.484536648864
	endif
c
c 3 juillet 1989 RGO
c
	if(iobs.eq.29) then
		az=   4033.163997666543
		bz= -1.7250252550293776E-03
		cz=  0.2625161707892094
		dz=   4552.839965165843
		ez=   1537.286424045690
	endif
c
c 11 juillet 1990, ESO
c
	if(iobs.eq.30) then
		az=   5571.106693606581
		bz=  -1.154739958213783
		cz=  0.2625161707892094
		dz=  -2875.444439204310
		ez=   2081.009021497123
	endif
c
c 18 aout 1991, CFHT
c
	if(iobs.eq.31) then
		az=   6006.373323781802
		bz=  -2.009446655224144
		cz=  0.2625161707892094
		dz=   1996.204285791411
		ez=   2237.472573209267
	endif
c
c 13 decembre 1989, KITT PEAK
c
	if(iobs.eq.33) then
		 az=   5418.033286716989
		 bz=  -2.239400951745198
		 cz=  0.2625161707892094
		 dz=   3089.816726590284
		 ez=  -2122.083904664141
	endif
c
c 18 juillet 1997, BUNDABERG
c
	if(iobs.eq.34) then
		 az=   5786.698039342732
		 bz=   2.574449733127754
		 cz=  0.2625161707892094
		 dz=  -2511.910382036621
		 ez=   1980.675695932511
	endif
c
c 18 juillet 1997, DUCABROOK
c
	if(iobs.eq.35) then
		 az=   5834.819414149405
		 bz=   2.488369221754767
		 cz=  0.2625161707892094
		 dz=  -2413.036306152251
		 ez=   1997.146719111260
	endif
c
c 18 juillet 1997, LOCHINGTON
c
	if(iobs.eq.36) then
		 az=   5832.684447406927
		 bz=   2.489750068081304
		 cz=  0.2625161707892094
		 dz=  -2417.442605053377
		 ez=   1996.415960963919
	endif
c
c 18 juillet 1997, BROWNSVILLE
c
	if(iobs.eq.37) then
		 az=   5737.408486303362
		 bz=  -1.787354233188314
		 cz=  0.2625161707892094
		 dz=   2609.154248137038
		 ez=   1963.804827761283
	endif
c
c 18 juillet 1997, CHILLAGOE
c
	if(iobs.eq.38) then
		 az=   6096.359978965665
		 bz=   2.437446717372623
		 cz=  0.2625161707892094
		 dz=  -1755.762578342800
		 ez=   2086.667035656214
	endif
c
c 13 mai 1971, PRETORIA
c
	if(iobs.eq.39) then
		 az=   5747.974253394440
		 bz=  0.2961094690474775
		 cz=  0.2625161707892094
		 dz=  -2595.465541466093
		 ez=   1947.570737470130
	endif
c
c 10 octobre 1999, ESO
c
	if(iobs.eq.40) then
		az=   5571.106693606581
		bz=  -1.442815571489685
		cz=  0.2625161707892094
		dz=  -3046.627456154876
		ez=  -1027.672034935083
	endif
c
c 10 octobre 1999, Paranal
c
	if(iobs.eq.41) then
		az=   5803.826501125843
		bz=  -1.437081437676362
		cz=  0.2625161707892094
		dz=  -2597.206813755162
		ez=  -1070.600603945883
	endif
c
c 10 octobre 1999, MONT MEGANTIC
c
	if(iobs.eq.42) then
		az=   4482.493468480031
		bz=  -1.450176497609971
		cz=  0.2625161707892094
		dz=   4446.137690929423
		ez=  -826.8614186877025
	endif
c
c 10 octobre 1999, PIC DU MIDI
c
	if(iobs.eq.43) then
		az=   4678.867740969617
		bz= -0.2058323718722671
		cz=  0.2625161707892094
		dz=   4250.102179859755
		ez=  -863.0855226794399
	endif
c
c 10 octobre 1999, CATALINA STATION
c
	if(iobs.eq.44) then
		az=    5391.5623189067
		bz=   -2.1409540572585
		cz=   0.26251617078921
		dz=    3342.5324532713
		ez=   -994.55245151000
	endif
c
c 10 octobre 1999, KITT PEAK
c
	if(iobs.eq.45) then
		az=    5418.0121986745
		bz=   -2.1560933340785
		cz=   0.26251617078921
		dz=    3300.4347111626
		ez=   -999.43151832017
	endif
c
c 03 decembre 1999, ESO
c
	if(iobs.eq.50) then
		az=   5571.106693606581
		bz= -0.6878888835558222
		cz=  0.2625161707892094
		dz=  -3021.611089937851
		ez=  -1243.562760978456
	endif
c
c
c
	return
	end
c
c
c				FIN DE DONNEEPAR
c
c
c ***************************************************************************
c
c Calcule la position de l'etoile par rapport au centre de la planete, vue
c depuis le centre de la Terre.
c
	subroutine geocentrique (iobs,iopns,heure, ksig,etag,dis)
c
c	Entrees: iobs, code de l'observation.
c
c		 iopns, option pour definir si l'etoile est passee au nord
c		 (1) ou au sud (-1) du centre de Neptune.
c
c		 heure: TU en heure decimale.
c
c	Sorties: ksig, etag, coordonnees geocentriques de l'etoile dans le
c		 plan du ciel (km), positif vers l'est et vers le nord.
c
c		 dis: distance de Neptune a la Terre (km)
c
	implicit real*8 (a-h,o-z)
	real*8 ksig
c
c 15 juin 1983
c
	if(iobs.eq.0.or.iobs.eq.1) then
		aksi=    1.942303474450455
		bksi=    85702.81847735686
		cksi=   -1263367.942181461
		aeta=    1.892947775312507
		beta=   -1129.914020161678
		ceta=    23655.34229969175
		dis= 4376394071.477264
	endif
c
c 22 juillet 1984
c
	if(iobs.eq.2.or.iobs.eq.3) then
		aksi=   -18.67283714592668
		bksi=    71550.57717923506
		cksi=   -469254.6150665939
		aeta=    1.480827628059842
		beta=    1158.262533592552
		ceta=    19667.95853956588
		dis= 4395137333.936507
	endif
c
c 20 aout 1985
c
	if(iobs.eq.4.or.iobs.eq.5.or.iobs.eq.6) then
		aksi=   -31.36317246977069
		bksi=    38768.58017673461
		cksi=   -234232.8115694095
		aeta=   0.4041263408104783
		beta=    2503.112985674425
		ceta=   -13460.16690614829
		dis= 4439683332.133106
	endif
c
c 23 avril 1986
c
	if(iobs.eq.7) then
		aksi=    33.55208725753437
		bksi=    27149.05526034722
		cksi=   -365676.9667433390
		aeta=    1.973074862959720
		beta=   -175.3318809121326
		ceta=    21389.78900241895
		dis= 4453679600.403858
	endif
c
c 23 aout 1986
c
	if(iobs.eq.8) then
		aksi=   -30.88380063709974
		bksi=    38031.78749215990
		cksi=   -257845.1327237141
		aeta=  -0.1338122440749885
		beta=    3063.687368666440
		ceta=   -6704.415821574944
		dis= 4439406120.422355
	endif
c
c 22 juin 1987
c
	if(iobs.eq.9) then
		aksi=    4.492478860314737
		bksi=    85165.88652529870
		cksi=   -269899.3171767063
		aeta=    1.747848190517187
		beta=    4585.489078748795
		ceta=   -36812.43174520412
		dis= 4371993502.939207
	endif
c
c 09 juillet 1987
c
	if(iobs.eq.10) then
		aksi=   -6.486260664043584
		bksi=    84218.59210283111
		cksi=   -746160.9173202395
		aeta=    1.237434074653038
		beta=    4836.442148246808
		ceta=   -40174.46126394378
		dis= 4373271584.580111
	endif
c
c 02 aout 1988, option sud
c
	if(iobs.eq.11.and.iopns.eq.-1) then
		aksi=   -19.56064509897442
		bksi=    69536.24564791580
		cksi=   -393637.7245280153
		aeta=   6.7265929632981170E-02
		beta=    5619.777111272552
		ceta=   -33707.39531111592
		dis= 4392158746.336141
	endif
c
c 02 aout 1988, option nord
c
	if(iobs.eq.11.and.iopns.eq.1) then
		aksi=   -19.56064509897442
		bksi=    69536.24564791580
		cksi=   -394249.7245280153
		aeta=   6.7265929632071675E-02
		beta=    5619.777111272562
		ceta=   -26606.39531111594
		dis= 4392158746.336141
	endif
c
c 25 aout 1988
c
	if(iobs.eq.12) then
		aksi=   -31.00300270080061
		bksi=    40931.23850190353
		cksi=   -387766.1386386958
		aeta=   -1.130498074194406
		beta=    4288.991870033216
		ceta=   -33643.21288963594
		dis= 4432503022.557647
	endif
c
c 12 septembre 1988
c
	if(iobs.eq.13.or.iobs.eq.14) then
	 	aksi=   -36.47633360781879
		bksi=    11717.77715694112
		cksi=   -209591.9119793359
		aeta=   -1.904459243852159
		beta=    2860.016088417047
		ceta=   -59480.34572117303
		dis= 4474605637.823063
	endif
c
c 07 juillet 1989
c
	if(iobs.eq.15.or.iobs.eq.16.or.iobs.eq.17) then
		aksi=  -3.257828389487372
		bksi=   85462.84527749704
		cksi=  -2051365.618035949
		aeta=   1.350815680012083
		beta=   7534.611673365305
		ceta=  -178303.6596690045
		dis=    4368534143.099570
	endif
c
c 03 juillet 1989
c
	if(iobs.ge.18.and.iobs.le.29) then
		aksi=   154.3501285232430
		bksi=   56946.23864750977
		cksi=  -1379125.481269501d+00 + 8186.95337d+00
		aeta=   18.19046992544551
		beta=   14816.36735984861
		ceta=  -352335.6030745959d+00 +  659.010413d+00
		dis=     1.34958d+9
	endif
c
c 11 juillet 1990, option sud
c
	if(iobs.eq.30.and.iopns.eq.-1) then
		aksi=  -3.863106300424306
		bksi=   85044.93269569209
		cksi=  -410135.6397301873
		aeta=   1.058657302910689
		beta=   8941.520638252781
		ceta=  -45923.34705329606
		dis=    4.367406147085697d+9
	endif
c
c 11 juillet 1990, option nord
c
	if(iobs.eq.30.and.iopns.eq.1) then
		aksi=  -3.863106300423851
		bksi=   85044.93269569208
		cksi=  -410953.4697301873
		aeta=   1.058657302911996
		beta=   8941.520638252765
		ceta=  -38717.93505329601
		dis=    4.367406147085697d+9
	endif
c
c 18 aout 1991
c
	if(iobs.eq.31) then
		aksi=  -24.40691755720457
		bksi=   60127.14789163109
		cksi=  -449951.5535721480
		aeta=  -1.750945326084633
		beta=    7931.327424852687
		ceta=  -70269.61713788468
		dis=     4.401287066d+9
	endif
c
c 11 juillet 1992
c
	if(iobs.eq.32) then
		aksi=  0.
		bksi=  84900.
		cksi=  0.
		aeta=  0.
		beta=  0.
		ceta=  0.
		dis=   4.411d+9
	endif
c
c 13 decembre 1989
c
	if(iobs.eq.33) then
		aksi=      11.00821939136290
		bksi=   56983.83056602708
		cksi= -561096.5129411838 	- 2843.
		aeta=      -4.052471612248908
		beta=   -4234.939153105840
		ceta=   48215.55184021784 	+ 2962.
		dis=   6.27098534d+8
	endif
c
c 18 juillet 1997
c le developpement polynomial de ksi et eta a ete fait a partir des 3 positions
c geocentriques donnees par astrom.for a 9:18, 10:18 et 11:18
c
c Ephemeride de Wasserman:
c
c	if(iobs.ge.34.and.iobs.le.38) then
c		write(*,*)  'Ephemeride Wasserman'
c		aksi=   353.0182477415310
c		bksi=     76216.82331460162
c		cksi=   -833528.7790223214
c		aeta=   18.63183044579000
c		beta=     29609.19774998216
c		ceta=   -308276.4216390371
c		dis=          4.358588027d+9	! distance a 10h30 UTC
c	endif
c
c Ephemeride d'Arlot:
c
	if(iobs.ge.34.and.iobs.le.38) then
		write(*,*)  'Ephemeride Arlot'
		aksi=       352.5640621443694
		bksi=     76203.31084798171
		cksi=   -798992.0179521392
		aeta=        19.01011986404762
		beta=     29585.55221196075
		ceta=   -304745.1693688219
		dis=          4.358588027d+9	! distance a 10h30 UTC
	endif
c
c 13 mai 1971
c
	if(iobs.eq.39) then
		aksi= -1.488577517268880
		bksi=  57561.83130818815
		cksi= -1111489.296174108 + 680.189
		aeta= -31.68944676065369
		beta= -9505.730498439629
		ceta=  136718.4116531998 + 1552.323
		dis=   653778296.1930600	! distance a 19h10 UTC
	endif
c
c 10 octobre 1999
c
c NB. les corrections de -9.5 et -14.8 km sur ksi et eta sont les corrections
c sur la position de l'etoile (Hipparcos + Tycho2). Voir les explications
c detaillees dans le cahier Jupiter no. 2, 23 aout 2000.
c
	if(iobs.ge.40.and.iobs.le.45) then
		aksi=       10.15180133964634      ! fit de dg 2 sur positions astrom
		bksi=    51313.6451151072	   ! a 6:25, 7:05, 7:45
		cksi=  -385946.2044279141 - 9.5d0
		aeta=        2.223185039849461
		beta=    19139.3480119430
		ceta=   -78449.38182224386 -14.8d0
		dis= 596290141.6271607             ! distance a 7h UTC
	endif
c
c 03 decembre 1999
c
	if(iobs.ge.50.and.iobs.le.50) then
		aksi=    -17.72979042000304    	! fit dg 2 sur positions astrom
		bksi=  57100.41072368165 	! de 2h a 10h par pas de 1h
		cksi= -247953.2760419623
		aeta=      -8.052092188354493
		beta=   15288.71100005913
		ceta= -111377.3256043004
		dis= 1243957273.788514	! distance a 4h30 UTC
	endif
c
c 03 decembre 1999
c
        if(iobs.ge.50.and.iobs.le.50) then
                aksi=    -17.72979042000304     ! fit dg 2 sur positions astrom
                bksi=  57100.41072368165        ! de 2h a 10h par pas de 1h
                cksi= -247953.2760419623
                aeta=      -8.052092188354493
                beta=   15288.71100005913
                ceta= -111377.3256043004
                dis= 1243957273.788514  ! distance a 4h30 UTC
        endif
c
c 03 decembre 1999
c
        if(iobs.ge.50.and.iobs.le.50) then
                aksi=    -17.72979042000304     ! fit dg 2 sur positions astrom
                bksi=  57100.41072368165        ! de 2h a 10h par pas de 1h
                cksi= -247953.2760419623
                aeta=      -8.052092188354493
                beta=   15288.71100005913
                ceta= -111377.3256043004
                dis= 1243957273.788514  ! distance a 4h30 UTC
        endif
c
c 08 septembre 2001
c
	if(iobs.ge.51.and.iobs.le.0107) then
c		aksi= -4.436821688372675  ! ephemeride Titania JPL DE405+URA027
c		bksi=  65839.09165277994  ! voir ephem_titania.for, seulement
c		cksi= -130948.1372805720  ! valide entre 01:51 et 02:09 !!!!
c		aeta=  15.46120331726638  ! attention: eta est ici l'oppose
c		beta=  35352.98011201547  ! du eta sur la carte du monde
c		ceta= -73979.79129466793  !
c		dis=   2850406650.736762  ! distance a 2h00 UTC

		aksi=  48.64235461878525  ! ephemeride Titania JPL DE405+URA027
		bksi=  65629.39695222813  ! ephem_titania.for, seulement
		cksi= -130741.6686692581  ! valide entre 01:30 et 02:30 !!!!
		aeta=  75.20103184765412  !
		beta=  35124.53064413106  !
		ceta= -73762.37329491012  !
		dis=   2850406650.736762  ! distance a 2h00 UTC
	endif
c
c 01 aout 2003
c
	if(iobs.ge.110.and.iobs.le.0111) then
		aksi= -27.8013689068921   ! ephemeride Titania JPL DE405+URA083
		bksi=  72905.5136599680   ! valide entre 03:30 et 05:30 !!!! 
		cksi= -322960.566266087   !
		aeta= -60.6401820325568   !
		beta=  15555.4391218765   !
		ceta= -72956.5579275708   !
		dis=   2855680607.33418   ! distance a 4h28 UTC
	endif
c
c 20 juillet 2002
c
    	if(iobs.ge.0120.and.iobs.le.0130) then
c		aksi=  29.6772241956678  ! ephemeride Pluton JPL PLU006.DE405
c		bksi=  58766.3004361403  ! ephem_pluton.for, seulement
c		cksi= -102002.523039657  ! valide entre 01:40 et 02:10 !!!!
c		aeta= -27.2503598956309  !
c		beta=  11446.0949212829  !
c		ceta= -20593.2438224984  !

		aksi= -53.8678038423132  ! ephemeride Pluton JPL PLU006.DE405
		bksi=  59077.2215402295  ! ephem_pluton.for, seulement
		cksi= -102290.366036936  ! valide entre 01:15 et 02:15 !!!!
		aeta=  13.0635586250019  !
		beta=  11292.2591461824  !
		ceta= -20447.3448708787  !

		dis=   4455107497.d0     ! distance a 1h44 UTC
	endif
c
c 21 aout 2002
c
        if(iobs.ge.0150.and.iobs.le.0153) then
c		aksi= -35.8833049527989    ! ephemeride Pluton JPL PLU006.DE405
c		bksi=  11532.7100960759    ! ephem_pluton.for, seulement
c		cksi= -80482.6605361820    ! valide entre 05:50 et 07:50 !!!!
c		aeta=  7.70894763097522    !
c		beta=  21219.3460739259    !
c		ceta= -148898.031823752    !
c
c Valeurs corrigees le 7 novembre 2002 (voir cahier):
c
		aksi= -35.8833046867142	   ! ephemeride Pluton JPL PLU006.DE405
		bksi=  11532.7133463513    ! ephem_pluton.for, seulement
		cksi= -80330.9746948468    ! valide entre 05:50 et 07:50 !!!!
		aeta=  7.70894735629554    !
		beta=  21219.3461622106    !
		ceta= -148898.032439279    !
		dis=  4523870353.94604d0   ! distance a UTC
        endif
c
c 15 decembre 2002
c
        if(iobs.ge.0160.and.iobs.le.0173) then
		aksi=  2368.85751276666 ! ephemeride de Tethys 
		bksi=  4293.86235460411 ! {source: SAT077 + DE-0406LE-0406}
		cksi= -936162.313487452 ! entre 18:30 et 19:30 UT
		aeta=  880.350355418138 !
		beta= -46893.8760522410 !
		ceta=  569869.950785978 !
		dis=   1204496088.78652 ! distance a 19:02 UTC
        endif
c
c Debut avril 2003
c
        if(iobs.ge.0200.and.iobs.le.0207) then
		aksi= -32.9805967613808 	! ephemeride de jupiter
		bksi=  5122.76514086033  	! source: Horizon 
		cksi= -131847.858377470		! entre 01/04/2003 00h et 09/04/2003 00h
		d3ksi= -4.153984623440521E-003 	! attention: degre 4 !!!!!!!
		d4ksi=  3.407885007083219E-006
		aeta=  9.11523648114908
		beta= -870.886758856861
		ceta= -43848.2095365587
		d3eta=  1.028052165914062E-003
		d4eta= -9.566312276628949E-007
		dis=   717996708.091873		! distance a 19:02 UTC
        endif
c
c 14 novembre 2003, 1ere occultation 
c
        if(iobs.ge.0220.and.iobs.le.0242) then
		aksi= -125.907891758718 ! ephemeride de Titan, 14 nov. 2003
		bksi=  40683.3622051077 ! {source: DE406 + TASS1.7 IMCCE}
		cksi= -6613.79469710593 ! entre 23:45 UT (13 nov.) 
		aeta=  13.8247955411625 ! et 00:45 UT (14 nov. 2003)
		beta=  5017.97969939133 ! etoile: 06:55:20.9672
		ceta=  3382.16252929752 !         22:06:07.812
		dis=   1253968185.11546 ! distance a 00:15 UT
        endif
c
c 14 novembre 2003, 2eme occultation 
c
        if(iobs.ge.0243.and.iobs.le.0245) then
		aksi=  -133.360277979437 ! ephemeride de Titan, 14 nov. 2003
		bksi=   40750.1385236864 ! {source: DE406 + TASS1.7 IMCCE}
		cksi=  -276952.762773133 ! entre 06:30 UT  
		aeta=   5.17438843970183 ! et    07:15 UT (14 nov. 2003)
		beta=   5089.47693453041 ! etoile: 06:55:17.7690
		ceta=  -36826.7926416029 !         22:06:01.226  (Manek)
		dis=   1253522913.93509  ! distance a 07:00 UT
        endif
c
c 11 juillet 2005, Charon 
c
        if(iobs.ge.0260.and.iobs.le.0274) then
                aksi=  16.0203034027472 ! ephemeride de Charon 11 Jul 2005
                bksi=  74791.8000385439 ! {source: DE413 + Pl013 JPL horizon) 
                cksi= -270871.525875259	! entre 03:11 et 04:09 UT et etoile
                aeta= -15.7760712054410	! LNA: 17:28:55.01670 -15:00:54.7260
                beta=  8381.69327259956 !
                ceta= -28789.1657240105	!
c		
c		aksi=  16.0202941795724	! ephemeride de Charon 11 Jul 2005
c		bksi=  74791.7877608504	! {source: DE413 + Pl013 JPL horizon) 
c		cksi= -272040.562475468	! entre 03:11 et 04:09 UT et etoile
c		aeta= -15.7760712614245	! Stone/JLX: 17:28:55.0130 -15:00:54.719
c 		beta=  8381.68973845561 !
c		ceta= -28636.4888139096	!
                dis=   4498462503.71029	! distance a 03:39 UT
        endif
	
c
c 10 avril 2006, Pluton p363 
c
        if(iobs.ge.0300.and.iobs.le.0301) then
!				aksi=  33.5459681996555 ! ephemeride de Pluton 10 avril 2006
!				bksi=  20368.3817852176 ! {source: DE413 JPL horizon prise en 2006)
!				cksi= -101407.403274610	! entre 03:00 et 12:00 UT et etoile
!				aeta=  1.69244235437714	! R. Behrend 17:46:06.880 -15:46:10.11
!				beta= -8040.23261734268 !
!				ceta=  55335.6250221367 !
!				dis=   4590832993.82593	! distance a 05:15 UT

!				aksi=  33.811867246305610 ! ephemeride de Pluton 10 avril 2006
!				bksi=  20364.236650643266 ! {source: DE413 JPL horizon prise en avril 2010: ~20-30 km difference avec version 2006) 
!				cksi= -101423.87372647502 ! entre 00:00 et 12:00 UT et etoile
!				aeta=  1.1969582465069379 ! R. Behrend 17:46:06.880 -15:46:10.11
!				beta= -8030.5733447102975 !
!				ceta=  55298.033328949416 !
!				dis=   4590832959.6769562 ! distance a 05:15 UT

				aksi=  34.949667348497087 ! ephemeride BARYCENTRIQUE de Pluton 10 avril 2006
				bksi=  20330.806570660599 ! {source: DE413 JPL horizon prise en avrl 2010) 
				cksi= -102613.50393876107 ! entre 00:00 et 12:00 UT et etoile
				aeta=  1.2020399248967806 ! R. Behrend 17:46:06.880 -15:46:10.11
				beta= -7949.6114426183849 !
				ceta=  54765.375185724406 !
				dis=   4590832959.6769562 ! distance a 05:14 UT
       endif
c
c 12 juin 2006, Pluton p384 
c
        if(iobs.ge.0310.and.iobs.le.0315) then
c                aksi=  8.39246409556745d0  ! ephemeride de Pluton 11 juin 2006
c                bksi=  84865.6720430911d0  ! {source: DE413 JPL horizon) 
c                cksi= -1389868.83861283d0  ! entre 15:20 et 17:20 UT et etoile
c                aeta=  4.46372095703340d0  ! R. Stone 17:41:12.0940 -15:41:34.446
c                beta=  1547.45194725201d0  !
c                ceta= -22404.9283904824d0  !
c                dis=   4506178082.97840d0  ! distance a 16:20 UT
		
                aksi=   8.39242344406739d0  ! ephemeride de Pluton 11 juin 2006
                bksi=   84865.6799008045d0  ! {source: DE413 JPL horizon) 
                cksi=  -1395263.72502086d0  ! entre 15:20 et 17:20 UT et etoile
                aeta=   4.46371147680553d0  ! R. Behrend 17:41:12.0769 -15:41:34.485
                beta=   1547.42456095162d0  !
                ceta=  -23256.4932888986d0  !
                dis=    4506177834.d0       ! distance a 16:23 UT
        endif
c
c 06 aout 2006, Pluton + Hydra 
c
        if(iobs.ge.0320.and.iobs.le.0320) then
				aksi=  -21.3201575882771d0  ! ephemeride de Pluton 06 aout 2006
                bksi=   48722.0355114296d0  ! {source: DE413 JPL horizon) 
                cksi=   12283.0228860494d0  ! entre 23:30 et 00:30 UT et etoile
                aeta=   9.19770261449203d0  ! R. Behrend 17:36:07.8427 -15:49:32.792
                beta=   14109.9818275996d0  !
                ceta=  -66742.7814413917d0  !
                dis=    4560274256.45207d0  ! distance a 00:00 UT
        endif
c
c 18 mars 2007, Pluton
c
        if(iobs.ge.0330.and.iobs.le.0349) then
				aksi=   35.6229529772918d0  ! ephemeride de Pluton 18 mars 2007
                bksi=  -25313.8448324031d0  ! {source: DE413 JPL horizon) 
                cksi=   274328.309015637d0  ! entre 10:01 et 11:59 UT et etoile
                aeta=  -1.89851853926336d0  ! R. Behrend 17:55:05.6971  -16:28:34.360 (mail 4 mars 2007)
                beta=  -6986.11411524670d0  ! ATTENTION, PENSEZ A CHANGER LE SIGNES DES COEFFS DE ETA !
                ceta=   72125.7373370588d0  !
                dis=    4677252720.10885d0  ! distance a 10:51 UT
        endif
c
c 12 mai 2007, Pluton, etoile B
c
        if(iobs.ge.0360.and.iobs.le.0367) then
				aksi=   21.5366295615331d0  ! ephemeride de Pluton 12 mai 2007
                bksi=   63087.7852356215d0  ! {source: DE413 JPL horizon) 
                cksi=  -249738.025096498d0  ! entre 03:00 et 05:40 UT et etoile (cf. cahier 28 jan 07)
                aeta=   5.13040887998829d0  ! M. Assafin 17:53:32.1006  -16:22:47.367 (mail 30 nov 2006) AVEC correction UCAC2 ---> ICRF
                beta=  -3156.25518332301d0  ! ATTENTION, PENSEZ A CHANGER LE SIGNES DES COEFFS DE ETA !
                ceta=   13674.5136952876d0  !
                dis=    4559807710.67633d0  ! distance a 04:00 UT
        endif
c
c 14 juin 2007, Pluton
c
        if(iobs.ge.0380.and.iobs.le.0389) then
!				aksi=   3.96918276864017d0  ! ephemeride de Pluton 14 juin 2007
!               bksi=   85010.2881595809d0  ! (source: DE413 JPL horizon) 
!               cksi=  -124027.662996917d0  ! entre 00:30 et 02:30 UT et etoile (cf. cahier 28 jan 07)
!				aeta=   4.64858229256220d0  ! ancienne valeur 2006 ???
!               beta=   3029.48127208737d0  ! ATTENTION, PENSEZ A CHANGER LE SIGNES DES COEFFS DE ETA !
!               ceta=  -3347.53381312977d0  !
!
!               aksi=   3.96917770138680d0  ! ephemeride de Pluton 14 juin 2007     
!               bksi=   85010.2892788936d0  ! {source: DE413 JPL horizon)  
!				cksi=  -124660.050380279d0  ! entre 00:30 et 02:30 UT et etoile (cf. cahier **28 jan 07**)
!				aeta=   4.64858687935066d0  ! etoile M. Assafin 17:50:20.7392 -16:22:42.210 cahier 04 juin 07 
!               beta=   3029.47688675422d0  ! ATTENTION, PENSEZ A CHANGER LES SIGNES DES COEFFS DE ETA !
!               ceta=  -2776.27863975939d0  !
!               dis=    4531866115.03948d0  ! distance a 01:27 UT

				aksi=   3.3069796866548131d0  ! ephemeride de Pluton 14 juin 2007     
				bksi=   85012.532709589446d0  ! {source: DE413 JPL horizon)
				cksi=  -124659.31713962542d0  ! entre 00:30 et 02:30 UT et etoile (prise en mai 2010: ~ 2 mas plus au sud que ephem 2007)
				aeta=   6.1221406119695985d0  ! etoile M. Assafin 17:50:20.7392 -16:22:42.210 cahier 04 juin 07 
				beta=   3024.6712801333706d0  ! ATTENTION, PENSEZ A CHANGER LES SIGNES DES COEFFS DE ETA !
				ceta=  -2733.1730203821617d0  !!
				dis=    4531866134.1395998d0  ! distance a 01:27 UT
        endif
		
c
c 14 juin 2007, Pluton, ephemeride BARYCENTRIQUE Pluton/Charon
c
        if(iobs.ge.0390.and.iobs.le.0390) then
				aksi=   4.4486788441063254d0  ! ephemeride BARYCENTRIQUE de Pluton 14 juin 2007     
				bksi=   85069.504399371115d0  ! {source: DE413 JPL horizon)
				cksi=  -124848.42067226613d0  ! entre 00:30 et 02:30 UT et etoile (prise en mai 2010: ~ 2 mas plus au sud que ephem 2007)
				aeta=   3.5084773926635080d0  ! etoile M. Assafin 17:50:20.7392 -16:22:42.210 cahier 04 juin 07 
				beta=   3012.5128416223347d0  ! ATTENTION, PENSEZ A CHANGER LES SIGNES DES COEFFS DE ETA !
				ceta=  -790.98742178666294d0  !
				dis=    4531866134.1395998d0  ! distance a 01:27 UT
		endif
		
c
c 09 juin 2007, Pluton
c
        if(iobs.ge.0391.and.iobs.le.0392) then
                aksi=   8.62761174944535d0  ! ephemeride de Pluton 09 juin 2007     
                bksi=   83763.6609352640d0  ! (source: DE413 JPL horizon)  
				cksi=  -806418.288115340d0  ! entre 08:30 et 10:30 UT et etoile R. Behrend 17:50:50.6535 -16:22:29.297,
				aeta=   3.81491503615933d0  ! prediction catalogue 16 avril 2002 + correction UCAC2 ---> ICRF
                beta=   2075.02949766538d0  ! ATTENTION, PENSEZ A CHANGER LES SIGNES DES COEFFS DE ETA !
                ceta=  -16840.3500093869d0  !
                dis=    4533115457.44727d0  ! distance a 09:36 UT
        endif
c
c 31 juillet 2007, Pluton
c
        if(iobs.ge.0395.and.iobs.le.0397) then
				aksi=  -22.9100147816262d0  ! ephemeride de Pluton 31 juillet 2007     
				bksi=   59614.9583410571d0  ! (source: DE413 JPL horizon)  
				cksi=  -820090.372577806d0  ! entre 12:50 et 14:50 UT et etoile R. Martins 17:45:41.98943 -16:29:31.639,
				aeta=  0.110365485090483d0  ! 2eme mail 29 juillet 07 (erratum!) + correction UCAC2 ---> ICRF
				beta=   12617.6829426895d0	! ATTENTION, PENSEZ A CHANGER LES SIGNES DES COEFFS DE ETA !
				ceta=  -168166.701438335d0	!
				dis=    4571240644.01135d0  ! distance a 13:46 UT
c
c il y a une erreur dans la solution ci-dessous (fit sur une ephemeride deja offsettee?)
c
c				aksi=  -22.9100029990122d0	! ephemeride de Pluton 31 juillet 2007     
c                bksi=  59614.9981912552d0	! (source: DE413 JPL horizon)
c				cksi=  -817874.723542103d0	! entre 12:50 et 14:50 UT et etoile R. Martins 17:45:41.98943 -16:29:31.639,
c				aeta=  0.110354108179081d0	! 2eme mail 29 juillet 07 (erratum!) + correction UCAC2 ---> ICRF
c                beta=  12617.6706653838d0  ! ATTENTION, PENSEZ A CHANGER LES SIGNES DES COEFFS DE ETA !
c                ceta= -169496.250643729d0	!
c                dis=   4571240644.01135d0  ! distance a 13:46 UT
        endif
c
c 21 mai 2008, Triton
c
        if(iobs.ge.0400.and.iobs.le.0410) then
				aksi=   122.681931482477d0  ! ephemeride de Triton 21 mai 2008   
				bksi=  -23622.8562816940d0  ! (source: JPL horizon)  
				cksi=   43794.3562521354d0  ! entre 00:00 et 04:00 UT et etoile R. Martins 21:46:11.049528 -13:46:45.9484716,
				aeta=   309.212172128612d0  ! mail 27 avril 08
				beta=  -3496.01226220196d0	! ATTENTION, PENSEZ A CHANGER LES SIGNES DES COEFFS DE ETA !
				ceta=   3874.94186234945d0	!
				dis=    4474933527.49276d0  ! distance a 01:45 UT
		endif
c
c 22 juin  2008, Pluton
c
        if(iobs.ge.0420.and.iobs.le.0434) then
				aksi=   1.25856884847781d0  ! ephemeride de Pluton 22 juin 2008   
				bksi=   85421.3100388775d0  ! (source: DE413 JPL horizon)  
				cksi=  -1639843.34686286d0  ! entre 18:10 et 20:10 UT et etoile R. Martins 17:58:33.01376909 -17:02:38.34870554,
				aeta=   2.47475168166466d0  ! mail 02 juin 08
				beta=   5695.85277015219d0	! ATTENTION, PENSEZ A CHANGER LES SIGNES DES COEFFS DE ETA !
				ceta=  -103871.343618670d0	!
				dis=    4558337860.00627d0  ! distance a 19:10 UT
		endif
c
c 22 juin  2008, Charon
c
        if(iobs.ge.0440.and.iobs.le.0443) then
				aksi=   7.67628810337737d0  ! ephemeride de Charon 22 juin 2008   
				bksi=   84856.0661675693d0  ! (source: DE413 + PLU017 JPL horizon)  
				cksi=  -1643371.13551012d0  ! entre 18:20 et 20:20 UT et etoile R. Martins 17:58:33.01376909 -17:02:38.34870554,
				aeta=   6.60709115209823d0  ! mail 02 juin 08
				beta=   6277.54499705377d0	! ATTENTION, PENSEZ A CHANGER LES SIGNES DES COEFFS DE ETA !
				ceta=  -121267.686897110d0	!
				dis=    4558324197.03955d0  ! distance a 19:20 UT
		endif
c
c 22 juin  2008, Charon, test avec DE413 + T08
c 20 oct 2008: je fais le test. En fait, donne un resultat coherent a < 0.1 km pres, ie ksi et eta de
c l'observateur par rapport a Charon differe entre les 2 methodes de la difference entre les 2 ephemerides,
c ce qui etait evident des le depart... Du moins en supposant que les 2 ephemerides donnent le meme vitesse
c pour Charon dans le plan du ciel.
c
c        if(iobs.ge.0440.and.iobs.le.0443) then
c				aksi=   1.25856884847781d0  ! ephemeride de Pluton 22 juin 2008   
c				bksi=   85421.3100388775d0  ! (source: DE413 JPL horizon)  
c				cksi=  -1639843.34686286d0  ! entre 18:10 et 20:10 UT et etoile R. Martins 17:58:33.01376909 -17:02:38.34870554,
c				aeta=   2.47475168166466d0  ! mail 02 juin 08
c				beta=   5695.85277015219d0	! ATTENTION, PENSEZ A CHANGER LES SIGNES DES COEFFS DE ETA !
c				ceta=  -103871.343618670d0	!
c				x0=  2759.8553402421035d0	! ephemeride T08 de Charon/Pluton
c				x1=  611.85665757575157d0
c				x2= -5.1860900655982221d0
c				x3= -8.65092473528612516d-2
c				y0=  15549.316745948438d0
c				y1= -297.22475821513672d0
c				y2= -18.995581740729421d0
c				y3=  0.25963191586671908d0
c				Xcharon= x0 + x1*heure + x2*(heure**2) + x3*(heure**3)
c				Ycharon= y0 + y1*heure + y2*(heure**2) + y3*(heure**3)
c				dis=    4558324197.03955d0  ! distance a 19:20 UT
c		endif
c
c 24 juin  2008, Pluton
c
        if(iobs.ge.0450.and.iobs.le.0450) then
				aksi=  -1.10011988720726d0  ! ephemeride de Pluton 24 juin 2008   
				bksi=   85273.4922374453d0  ! (source: DE413 JPL horizon)  
				cksi=  -905903.629086786d0  ! entre 09:36 et 10:36 UT et etoile R. Martins 17:58:22.392984 -17:02:49.349328,
				aeta=   6.17016734615709d0  ! mail 09 juin 08
				beta=   6027.10549612317d0	! ATTENTION, PENSEZ A CHANGER LES SIGNES DES COEFFS DE ETA !
				ceta=  -66444.4048215610d0	!
				dis=    4558640166.45430d0  ! distance a 10:36 UT
		endif
c
c 30 juin  2009, 2002 MS4
c
        if(iobs.ge.0460.and.iobs.le.0464) then
				aksi=    -4.4733154532409571d0	! ephemeride de 2002 MS4 30 juin 2009   
				bksi=    90244.834826651204d0	! (source: JPL horizon, 30 juin 2009)  
				cksi=   -2137170.9869009401d0	! entre 22:45 et 24:45 UT et etoile M. Assafin:  18:03:37.93133 -08:28:17.6971,
				aeta=   +15.065390509613124d0	! mail 25 juin 09
				beta=   +1203.8472514488851d0	! ATTENTION, PENSEZ A CHANGER LES SIGNES DES COEFFS DE ETA !
				ceta=   -35917.217592707675d0	!
				dis=     6920590048.3620806d0	! distance a 23:45 UT
		endif
c
c 25 aout 2008, Pluton
c
        if(iobs.ge.0470.and.iobs.le.0471) then
				aksi=   -35.223799185288954d0	! ephemeride de Pluton 25 aout 2008   
				bksi=    26148.120752967247d0	! (source: DE413 JPL horizon)  
				cksi=   -118441.10041256671d0	! entre 03:30 et 05:30 UT et etoile R. Martins 17:53:27.103164 -17:15:27.54630,
				aeta=    1.7819848315398303d0	! mail 23 juilet 08
				beta=    15600.931689354853d0	! ATTENTION, PENSEZ A CHANGER LES SIGNES DES COEFFS DE ETA !
				ceta=  	-72764.814667978659d0	!
				dis=     4645765739.8504114d0	! distance a 04:30 UT
		endif
c
c 19 fevrier 2010, Varuna
c
        if(iobs.ge.0480.and.iobs.le.0491) then		
				aksi=   -27.085326548563444d0	! ephemeride de Varuna 19 fevrier 2010   
				bksi=   +67762.901239410639d0 	! (source: ephem #18 JPL horizon)  
				cksi=   -1552234.8775654701d0	! entre 22:10 et 24:10 UT et etoile RM 04 fev. 2010,
				aeta=    4.6904675561950739d0	! 
				beta=   -9882.4604272275792d0 	! ATTENTION, PENSEZ A CHANGER LES SIGNES DES COEFFS DE ETA !
				ceta=  	 225457.48026839149d0	!
				dis=     6395653310.6800442d0 	! distance a 23:10 UT
		endif
c
c 14 fevrier 2010, Pluton
c
        if(iobs.ge.0500.and.iobs.le.0502) then		
				aksi=   22.079643370921985d0 	! ephemeride de Pluton 14 fevrier 2010   
				bksi=  -89561.379967677989d0  	! (source: DE413 JPL horizon) 
				cksi=   431402.07056522940d0 	! entre 03:47 et 05:47 UT et etoile  Rio 10 novembre 2009, voir table
				aeta=  -8.03334814107756756d-2  ! g3_occ_data_pluto_2010_table
				beta=  -4608.1802305182300d0  	! ATTENTION, PENSEZ A CHANGER LES SIGNES DES COEFFS DE ETA !
				ceta=  +20375.225102583456d0	!
				dis=    4847007325.6894636d0 	! distance a 04:47 UT
		endif
c
c 19 mai 2010, Pluton
c
        if(iobs.ge.0510.and.iobs.le.0512) then		
				aksi=   21.468482961861810d0    ! ephemeride de Pluton 19 mai 2010 
				bksi=   64277.838549352149d0    ! (source: DE413 JPL horizon)
				cksi=  -382721.51689993986d0    ! entre 04:56 et 06:56 UT et etoile  Rio 10 novembre 2009, voir table
				aeta=   1.9104059622432601d0    ! g3_occ_data_pluto_2010_table:
				beta=   2254.0515298318383d0    ! 18:20:16.7941 -18:11:48.070
				ceta=  -8353.9215106761403d0    !
				 dis=   4640898576.3397970d0 	! distance a 05:56 UT
		endif
c
c 04 juin 2010, Pluton
c
         if(iobs.ge.0515.and.iobs.le.0520) then		
                aksi=   13.726785943874347d0	! ephemeride de Pluton 04 juin 2010
                bksi=   78282.017758372313d0	! (source: DE413 JPL horizon)
                cksi=  -1225786.9781723302d0	! entre 14:36 et 16:36 et et etoile  Rio 10 novembre 2009, voir table
                aeta=   3.4430947319003735d0    ! g3_occ_data_pluto_2010_table:
                beta=   4925.8377412356176d0	! 18:18:47.9349 -18:12:51.794
                ceta=  -71885.836617232708d0	! 119 points, rms ksi et eta: 9.35 et 6.47 km
                 dis=   4622039468.1015949d0	! distance a 15:36
		endif
c
c 04 juillet 2010, Pluton
c
         if(iobs.ge.0530.and.iobs.le.0530) then		
                aksi=  -8.1797302350059908d0	! ephemeride de Pluton 04 juillet 2010
                bksi=   84529.155447361496d0	! (source: DE413 JPL horizon)
                cksi=  -172005.43209152107d0	! entre 01:02 et 03:02 et etoile Rio 02 juillet 2010:
                aeta=   3.3957897867394422d0	! 18:15:42.1027 -18:16:41.124
                beta=   9336.6344199989999d0	! 119 points, rms ksi et eta: 9.14 et 6.26 km
                ceta=  -12668.126339528944d0	!
                 dis=   4616032434.1436749d0	! distance a 02:02
		endif
c
c 17 aout 2010, Ceres
c
         if(iobs.ge.0540.and.iobs.le.0544) then		
                aksi=  -31.859680548147480d0	! ephemeride de Ceres, 17 aout 2010
                bksi=  -11522.259438149176d0	! (source:  JPL horizon #26)
                cksi=   277568.57026796800d0	! entre 21:41 et 23:41 et etoile UCAC2
                aeta=   1.3027782077446943d0	! 17:18:29.008  -27:26:38.89
                beta=   5762.2770568197311d0	! 119 points, rms ksi et eta: 0.66 et 0.48 km
                ceta=  -131746.03150922264d0	!
                 dis=   343144041.66042966d0	! distance a 22:41
		endif

c
c 08 janvier 2011, 2003AZ84
c
         if(iobs.ge.0545.and.iobs.le.0546) then
	        aksi=   3.7866712860850384d0	 ! ephemeride de AZ84, 08/02/2011
                bksi=   92862.070336359611d0	 ! (source:  JPL horizon #11)
                cksi=  -598409.83518198912d0	 ! entre 05:00 et 7:59 et etoile WFI
                aeta=  -8.1064245711991134d0	 ! 07 43 41.8220  +11 30 23.569
                beta=  -11328.406221697676d0	 ! 179 points, rms ksi et eta: 16.4 et 9.34 km
                ceta=   77474.762426329864d0	 !
                 dis=   6629883867.6944599d0	 ! distance km a 06:30
		endif

c
c 06 novembre 2010, Eris
c 
c ephemeride d'Eris, 06 novembre 2010 (source:  JPL horizon ***#28***)
c entre 01:20 et 03:20 et etoile Bill Owen, mail 4 nov. 2010: 01:39:09.9392, -04:21:12.140
c 119 points, rms ksi et eta: 13.72 et 19.87 km ---> 27 dec. 2010 je trouve rms ksi, eta: 29.70 km et 20.03 km?
c distance a 02:20

         if(iobs.ge.0550.and.iobs.le.0554) then
                 aksi=  -9.5271254042781948d0
                 bksi=   92641.460081157304d0
                 cksi=  -201993.01909115419d0
                 aeta=  -10.050914305598781d0
                 beta=   20455.604737461523d0
                 ceta=  -51226.932644434593d0
                  dis=   14322167213.076771d0
		 endif
c
c 06 novembre 2010, Eris
c 
c ephemeride d'Eris, 06 novembre 2010 (source:  JPL horizon ***#29***)
c entre 01:20 et 03:20 et etoile Bill Owen, mail 4 nov. 2010: 01:39:09.9392, -04:21:12.140
c 119 points, rms ksi et eta: 29.82 km et 19.96 km
c distance a 02:20
!
!        if(iobs.ge.0550.and.iobs.le.0554) then
!                 aksi=  -13.174829288458568d0
!                 bksi=   92664.649051483080d0
!                 cksi=  -205522.98371767512d0
!                 aeta=  -9.0362306132046797d0
!                 beta=   20450.211106095074d0
!                 ceta=  -52718.647603711790d0
!                  dis=   14322250622.779324d0
!		 endif
c 
c **05** novembre 2010, Eris
c 
c ephemeride d'Eris, **05** novembre 2010 (source:  JPL horizon **#29**)
c entre 01:00 et 06:00 et etoile Bill Owen, mail 4 nov. 2010: 01:39:09.9392, -04:21:12.140
c 299 points, rms ksi et eta: 30.00 km et 19.90 km
c distance a 03:30

         if(iobs.ge.0560.and.iobs.le.0560) then
                 aksi=  -9.9648885336646345d0
                 bksi=   93135.473818211933d0
                 cksi=  -2434962.6635225350d0
                 aeta=  -14.029968768212711d0
                 beta=   21117.620099796826d0
                 ceta=  -551780.64525171346d0
                  dis=   14321358526.201340d0
		 endif

c
c 23 avril 2011, Makemake
c 
c ephemeride de Makemake, 23 avril 2011 (source:  JPL horizon #43)
c entre 00:37 et 02:35 et etoile JL Ortiz, position dans la demande DDT ESO 18 avril 2011: 12:36:11.4140, +28:11:10.606
c 119 points, rms ksi et eta: 14.23 km et 10.76 km
c distance a 01:36
c
c        if(iobs.ge.0570.and.iobs.le.0575) then
c                 aksi=  -14.568764893669140d0
c                 bksi=   79372.750676897078d0
c                 cksi=  -120790.86096214066d0
c                 aeta=   19.898764792571455d0
c                 beta=  -5338.4961533930327d0
c                 ceta=   13546.204501973094d0
c                  dis=   7704196991.2267771d0
c		 endif


c
c 23 avril 2011, Makemake
c 
c ephemeride de Makemake, 23 avril 2011 (source:  JPL horizon #46)
c entre 00:37 et 02:35 et etoile WFI 01/08/2011
c fichier fort.20 (t,ksi):         119 points, rms=   14.227207     km
c fichier fort.21 (t,eta):         119 points, rms=   10.639701     km
c distance a 01:36
c
        if(iobs.ge.0570.and.iobs.le.0575) then
                 aksi=  -13.016745039421949d0
                 bksi=   79370.332632062418d0
                 cksi=  -126116.08898225593d0
                 aeta=   26.727158868591232d0
                 beta=  -5356.1413904692108d0
                 ceta=   14803.954513898290d0
                  dis=   7704217571.8428268d0
		 endif


c
c 04 mai 2011, Quaoar
c 
c ephemeride de Quaoar, 04 mai 2011 (source:  JPL horizon #17)
c entre 01:40 et 03:40 et etoile JL Ortiz, position dans la figure update 3 mai 2011:  17:28:50.8183 -15:27:42.657
c 119 points, rms ksi et eta: 12.67 km et 08.90 km
c distance a 02:40
c
c        if(iobs.ge.0584.and.iobs.le.0599) then
c				aksi=   20.936542766550701d0
c				bksi=   64292.915548761965d0
c				cksi=  -169349.10367880616d0
c				aeta=   1.5628926468725695d0
c				beta=  -13447.320977644898d0
c				ceta=   39025.625134937836d0
c				 dis=   6335560169.5390224d0
c		 endif
c
c 04 mai 2011, Quaoar
c 
c ephemeride de Quaoar, 04 mai 2011 (source:  JPL horizon **** JPL#17 ***)
c entre 02:01 et 03:59 et etoile Julio IAG 20/ago2011, obs IAG 18/ago (er RA 7mas e DEC 14mas)
c fichier fort.20 (t,ksi):         119 points, rms=   12.717681     km
c fichier fort.21 (t,eta):         119 points, rms=   8.8542242     km
c distance a 03:00
c
c        if(iobs.ge.0584.and.iobs.le.0599) then
c				aksi=	23.603559711156777d0
c				bksi=	64278.326494985005d0
c				cksi=  -177056.19566768905d0
c				aeta=	5.6470213296970542d0
c				beta=  -13470.844527931504d0
c				ceta=	35587.932683201121d0
c				 dis=	6335537408.8795786d0
c		 endif

c
c 04 mai 2011, Quaoar
c 
c ephemeride de Quaoar, 04 mai 2011 (source:  JPL horizon **** JPL#21 ***)
c entre 01:41 et 02:39 et etoile Julio IAG 20/ago2011, obs IAG 18/ago (er RA 7mas e DEC 14mas)
c fichier fort.20 (t,ksi):         119 points, rms=   12.754010     km
c fichier fort.21 (t,eta):         119 points, rms=   8.8974714     km
c distance a 02:40
c
        if(iobs.ge.0584.and.iobs.le.0599) then
				aksi=	25.224952169934113d0
				bksi=	64271.421671159311d0
				cksi=  -172122.76307958877d0
				aeta=	2.6467322079533915d0
				beta=  -13452.923614626066d0
				ceta=	36974.046549956816d0
				 dis=	6335613014.1210175d0
		 endif





c
c 04 juin 2011, Pluton
c 
c ephemeride de Pluton, 04 juin 2011 (source:  JPL horizon PLU07/DE413)
c entre 04:42 et 06:42 et etoile Rio 10 Nov. 2009 g3_occ_data_pluto_2011_data:  18 27 53.8249 -18 45  30.725 
c 119 points, rms ksi et eta: 0.0000875 km (???) et 06.48 km
c distance a 05:42: 
c
        if(iobs.ge.0600.and.iobs.le.0603) then
				aksi= -0.68253514904063195d0
				bksi=   76916.656470764021d0
				cksi=  -438146.40338949452d0
				aeta=   2.4686695530014049d0
				beta=   6248.2655140183597d0
				ceta=  -31695.961435145611d0
				dis=    4653752281.7049017d0
		 endif
c
c 04 juin 2011, Charon
c 
c ephemeride de Charon, 04 juin 2011 (source:  JPL horizon PLU07/DE413)
c entre 04:42 et 06:42 et etoile Rio 10 Nov. 2009 g3_occ_data_pluto_2011_data:  18 27 53.8249 -18 45  30.725 
c 119 points, rms ksi et eta: 9.13 km et 06.45 km
c distance a 05:42: 
c
        if(iobs.ge.0610.and.iobs.le.0613) then
				aksi=   5.2707520146668685d0
				bksi=   77080.566935426352d0
				cksi=  -424411.82476923097d0
				aeta=   2.4098484921429417d0
				beta=   5503.6940313531504d0
				ceta=  -25921.796924731709d0
				dis=    4653765029.9496012d0
		 endif
c
c 23 juin 2011, Pluton
c 
c ephemeride de Pluton, 23 juin 2011 (source:  JPL horizon barycentrique DE413 + offset Marc Buie)
c entre 10:15 et 12:15 et etoile Rio 10 Nov. 2009 g3_occ_data_pluto_2011_data:  18 25 55.4750 -18 48 07.015 
c 119 points, rms ksi et eta: 9.25 km et 06.47 km
c distance a 11:15 
c
        if(iobs.ge.0620.and.iobs.le.0629) then
                 aksi=   3.3397223015881536d0
                 bksi=   85186.385165767279d0
                 cksi=  -974443.27924155910d0
                 aeta=   3.0436657254591637d0
                 beta=   8986.8938558109730d0
                 ceta=  -101670.99971197746d0
                  dis=   4643352313.8514223d0
		 endif
c
c 23 juin 2011, Charon
c 
c ephemeride de Charon, 23 juin 2011 (source:  JPL horizon barycentrique DE413 + offset Marc Buie)
c entre 10:15 et 12:15 et etoile Rio 10 Nov. 2009 g3_occ_data_pluto_2011_data:  18 25 55.4750 -18 48 07.015 
c 119 points, rms ksi et eta: 9.30 km et 06.48 km
c distance a 11:15 
c
        if(iobs.ge.0630.and.iobs.le.0639) then
                 aksi=  -8.5254581593671901d0
                 bksi=   85643.730287262093d0
                 cksi=  -962916.88510382769d0
                 aeta=   2.9161058809503402d0
                 beta=   8241.4579473665417d0
                 ceta=  -93011.297854806733d0
                  dis=   4643352313.8514223d0
		 endif
c
c 27 juin 2011, Pluton
c 
c ephemeride de Pluton, 27 juin 2011 (source:  JPL horizon barycentrique DE413 + offset Marc Buie)
c entre 13:55 et 15:55 et etoile Rio 10 Nov. 2009 g3_occ_data_pluto_2011_data:  18 25 29.0100 -18 48 47.570  
c 119 points, rms ksi et eta: 9.22 km et 6.53 km
c distance a 14:55 
c
        if(iobs.ge.0640.and.iobs.le.0649) then				  
                 aksi= -0.87572953535618581d0
                 bksi=   85621.484178134488d0
                 cksi=  -1226969.4315801030d0
                 aeta=   2.4765295519496249d0
                 beta=   9404.2527330246048d0
                 ceta=  -133231.00506344080d0
                  dis=   4643181278.6059542d0
		 endif
c
c 27 juin 2011, Hydra
c 
c ephemeride de Hydra, 27 juin 2011 (source:  JPL horizon barycentrique DE413 + offset Marc Buie)
c entre 13:55 et 15:55 et etoile Rio 10 Nov. 2009 g3_occ_data_pluto_2011_data:  18 25 29.0100 -18 48 47.570  
c 119 points, rms ksi et eta: 9.22 km et 6.53 km
c distance a 14:55 
c
        if(iobs.ge.0650.and.iobs.le.0659) then
                 aksi=   1.1297043080276126d0
                 bksi=   85480.465547754487d0
                 cksi=  -1275145.4365299537d0
                 aeta=   1.3296748729530918d0
                 beta=   9902.8515480400602d0
                 ceta=  -143469.13941475889d0
                  dis=   4643181278.6059542d0
		 endif

c
c 03 fevrier 2012, 2003 AZ84
c 
c ephemeride de 2003 AZ84, 03 fevrier 2012 (source:  JPL horizon JPL#15)
c entre 18:47 et 20:47 et etoile JLX moyenne email 01/02/12. 
c
c  fichier fort.20 (t,ksi):	    119 points, rms=   13.432420     km
c  fichier fort.21 (t,eta):	    119 points, rms=   9.3298616     km 
c premier et dernier temps: 18 48 0.000  20 46 0.000
c distance donnee pour  t=  19 47 0.000
c


       if(iobs.ge.0660.and.iobs.le.0679) then
                 aksi=  -9.0801130611353074d0
                 bksi=   87645.659969759989d0
                 cksi=  -1737401.4624621109d0
                 aeta=  -6.4623116277441568d0
                 beta=  -18853.834217150477d0
                 ceta=   376525.77038443502d0
                  dis=   6620895942.8725510d0

		 endif


c
c 17 Fevrier 2012, Quaoar
c 
c ephemeride de Quaoar, 17 fevrier 2012 (source:  JPL horizon JPL#21)
c entre 03:30 et 5:30 et etoile WFI Rio 01 AOUT 2011: 
c
c fichier fort.20 (t,ksi):         119 points, rms=   13.551773     km
c fichier fort.21 (t,eta):         119 points, rms=   9.1429291     km
c  
c   premier et dernier temps: 03 31 0.000  05 29 0.000
c   distance donnee pour  t=  04 30 0.000
c


        if(iobs.ge.0680.and.iobs.le.0683) then
                 aksi=   17.199489037688181d0
                 bksi=  -65854.163239721660d0
                 cksi=   294399.42873631779d0
                 aeta=  -4.5830386084637666d0
                 beta=  -10445.085443234871d0
                 ceta=   39303.520299399053d0
                  dis=   6510624940.6791649d0
		 endif

c
c 03 December 2015, 2002 VE95
c 
c ephemeride de 2002 VE95, 03 December 2015 (source:  JPL horizon JPL#14)
c entre 12:00 et 13:00 et etoile GAIA DR1: 
c
c fichier fort.20 (t,ksi):          59 points, rms=   8.83168888     km
c fichier fort.21 (t,eta):          59 points, rms=   5.77530909     km
c  
c premier et dernier temps: 12 01 0.000  12 59 0.000
c distance donnee pour  t=  12 30 0.000
c


        if(iobs.ge.0685.and.iobs.le.0687) then
                 aksi=   25.683908684032758d0
                 bksi=   83151.827372091459d0
                 cksi=  -1038824.1191098636d0
                 aeta=  -27.207057301361147d0
                 beta=   13114.339518034809d0
                 ceta=  -155890.07542748141d0
                  dis=   4253776920.7331438d0
		 endif

c
c Quaoar-July 09, 2013
c 
c ephemeride de Quaoar in July 09, 2013 (source:  JPL horizon JPL#31)
c between 01:34 and 03:34 A.M for GAIA DR1 star: 
c
c fichier fort.20 (t,ksi):         119 points, rms=   12.0865288     km
c fichier fort.21 (t,eta):         119 points, rms=   8.76939201     km
c  
c premier et dernier temps: 01 35 0.000  03 33 0.000
c distance donnee pour  t=  02 34 0.000
c
  
        if(iobs.ge.0688.and.iobs.le.0688) then
                 aksi=  -10.001796395717975d0
                 bksi=   80886.469884400765d0
                 cksi=  -211905.78537631754d0
                 aeta=   2.9076350628297973d0
                 beta=   1510.7968421353185d0
                 ceta=  -3683.4688342377176d0
                  dis=   6301367743.6409912d0
		 endif

c
c Sedna-January 13, 2013
c 
c ephemeride de Sedna in January 13, 2013 (source:  JPL horizon JPL#28)
c between 11:25 and 13:25 for GAIA DR1 star: 
c
c fichier fort.20 (t,ksi):          59 points, rms=   8.83168888     km
c fichier fort.21 (t,eta):          59 points, rms=   5.77530909     km
c
c premier et dernier temps: 12 01 0.000  12 59 0.000
c distance donnee pour  t=  12 30 0.000
c
  
        if(iobs.ge.0689.and.iobs.le.0690) then
                 aksi=   25.683908684032758d0
                 bksi=   83151.827372091459d0
                 cksi=  -1038824.1191098636d0
                 aeta=  -27.207057301361147d0
                 beta=   13114.339518034809d0
                 ceta=  -155890.07542748141d0
                  dis=   4253776920.7331438d0
		 endif

c
c Eris August 29, 2013
c 
c ephemeride de Eris in August 29, 2013 (source:  JPL horizon JPL#57)
c between 14:38 and 16:38 for GAIA DR1 star: 
c
c fichier fort.20 (t,ksi):         119 points, rms=   30.3706493     km
c fichier fort.21 (t,eta):         119 points, rms=   19.9956760     km
c
c premier et dernier temps: 14 39 0.000  16 37 0.000
c distance donnee pour  t=  15 38 0.000
c
  
        if(iobs.ge.0691.and.iobs.le.0692) then
                 aksi=   20.832092892774199d0
                 bksi=   58518.673891025654d0
                 cksi=  -916668.47568179795d0
                 aeta=   5.7816161896389531d0
                 beta=   35700.282390462176d0
                 ceta=  -553670.74477819668d0
                  dis=   14326992675.357355d0
		 endif

c
c 2003 VS2 in November 07, 2014
c 
c ephemeride de 2003 VS2 in November 07, 2014 (source:  JPL horizon JPL#30)
c between 03:25 and 05:25 for GAIA DR1 star: 
c
c fichier fort.20 (t,ksi):         119 points, rms=   9.22811413     km
c fichier fort.21 (t,eta):         119 points, rms=   7.45238638     km
c
c premier et dernier temps: 03 26 0.000  05 24 0.000
c distance donnee pour  t=  04 25 0.000
c
  
        if(iobs.ge.0693.and.iobs.le.0698) then
                 aksi=   20.445461485434862d0
                 bksi=   76006.481478757400d0
                 cksi=  -335862.21452584077d0
                 aeta=   9.6040940091309963d0
                 beta=   1926.8179397260401d0
                 ceta=  -3681.9249089274977d0
                  dis=   5343005213.8857927d0
		 endif

        if(iobs.ge.0700.and.iobs.le.0703) then
                 aksi=   18.237615229178118d0
                 bksi=   76025.936731146387d0
                 cksi=  -335900.93277622899d0
                 aeta=   8.4469410132513048d0
                 beta=   1937.0656715047694d0
                 ceta=  -3716.2601611088394d0
                  dis=   5343005213.8857927d0
		 endif


*------------------------------------------------------------------------------------------------------------------------
c
c Cas particulier d'avril 2003 (Jupiter)
c
	if(iobs.ge.0200.and.iobs.le.0207) then
	 ksig= aksi*(heure**2) + bksi*heure + cksi
	 ksig= ksig + d4ksi*(heure**4) + d3ksi*(heure**3)  
	 etag= aeta*(heure**2) + beta*heure + ceta
	 etag= etag + d4eta*(heure**4) + d3eta*(heure**3)  
	 return
	endif
c
c cas particulier de Charon 22 juin 2008, ephem T08 + DE413
c
c A neutraliser!
c
c	if(iobs.ge.0440.and.iobs.le.0443) then
c	 ksig= aksi*heure*heure + bksi*heure + cksi
c	 etag= aeta*heure*heure + beta*heure + ceta
c	 ksig= ksig - Xcharon
c	 etag= etag - Ycharon
c	 return
c	endif
c
c cas general
c
	ksig= aksi*heure*heure + bksi*heure + cksi
	etag= aeta*heure*heure + beta*heure + ceta
c
c         
c
	return
	end
c
c
c				FIN DE GEOCENTRIQUE
c
c ***************************************************************************
c
	subroutine polpo (idate,itriton,dis,
     *                B,P,BP,aa,bb,cc,dd,aap,bbp,ccp,ddp)
c
c Sous-programme de calcul de l'angle de position et de l'elevation au-dessus
c du plan des anneaux. Calcul aussi les coeff. a,b,c,d de projection du plan
c du ciel dans le plan des anneaux.
c
c	Entrees: idate, code de la date.
c		 itriton, precise si Triton est pris en compte ou non dans
c		 la position du pole.
c		 dis, distance de la planete, km
c		 la position du pole de la Terre en 1950.
c
c	Sorties: B, elevation planetocentrique de la Terre (deg. decimal).
c		 P, angle de position du pole nord (deg. decimal).
c		 BP, elevation planetocentrique du Soleil (deg. decimal)
c		 aa,
c		 bb, x= aa*ksi + bb*eta			projection dans le plan des anneaux
c		 cc, y= cc*ksi + dd*eta
c		 dd.
c NB. (x,y) sont les coordonnees de la projection du point (ksi,eta), depuis le plan du ciel
c vers le plan equatorial de la planete. Ce point projete est repere par rapport aux deux vecteurs 
c unitaires "a" et "b" appartenant au plan equatorial de la planete
c Le vecteur "b" est porte par l'anse de l'anneau, i.e. appartient aussi au plan du ciel,   
c il est defini comme b= -cos(P).e + sin(P).n, ou "e" et "n" sont les vecteurs unitaires
c pointant vers l'est et le nord celestes locaux. Par exemple si P=0, "b" pointe vers l'anse ouest
c "a" est defini comme a= b x p, ou "p" est le vecteur unitaire porte par le pole de la planete.
c Donc "a", "b" et "p" forment un triedre direct. "a" est porte par le petit axe des anneaux qui
c point vers l'observateur.
c
c		 aap,
c		 bbp, xp= aap*x + bbp*x		projection dans le plan des anneaux avec comme
c		 ccp, yp= ccp*x + ddp*x		origine le noeud ascendant /t plan equatorial Terre
c		 ddp.                       ATTENTION!: utilise les (x,y) calcules plus haut 
c
c NB. (xp,yp) sont les coordonnees de la projection de (ksi,eta) dans le plan equatorial de la
c planete, avec comme axe d'origine le noeud ascendant du plan equatorial de la planete sur
c le plan equatorial terrestre. La longitude de ce point par rapport a ce noeud ascendant est donc
c arctan(yp/xp)
c
	implicit real*8 (a-h,o-z)
 100	format (1x,i3,1x,i2,1x,f7.3)
	real*8 nrp,nx,ny,nz,n1x,n1y,n1z,n1,n2x,n2y,n2z
	pi= 4.*datan(1.d+00)
	ua= 1.4959787d+08
c
c
c Coordonnees de Neptune et de son pole a differentes dates:
c
c 15 juin 1983
c
	if(idate.eq.0.or.idate.eq.1) then
		iah=   17
		iam=   51
		asec=  18.9
		idh=  -22
		idm=  -09
		dsec= -31.
c
		iahe=  12
		iame=  00
		asece= 51.41
		idhe=  89
		idme=  48
		dsece= 49.6
c
		iahs=  05
		iams=  31
		asecs= 02.2
		idhs=  23
		idms=  16
		dsecs= 29.
		diss=  1.0157697*ua
c
		if(itriton.eq.0) then
			iahp=  19
			iamp=  56
			asecp= 42.58
			idhp=  42
			idmp=  53
			dsecp= 25.4
		endif
		if(itriton.eq.1) then
			iahp=  19
			iamp=  56
			asecp= 42.58
			idhp=  42
			idmp=  53
			dsecp= 25.4
		endif
c
	endif
c
c 22 juillet 1984
c
	if(idate.eq.2.or.idate.eq.3) then
		iah=   17
		iam=   56
		asec=  49.8
		idh=  -22
		idm=  -14
		dsec= -02.
c
		iahe=  12
		iame=  00
		asece= 53.11
		idhe=  89
		idme=  48
		dsece= 27.5
c
		iahs=  08
		iams=  06
		asecs= 04.6
		idhs=  20
		idms=  17
		dsecs= 11.
		diss=  1.0160194*ua
c
		if(itriton.eq.0) then
			iahp=  19
			iamp=  56
			asecp= 44.77
			idhp=  42
			idmp=  53
			dsecp= 36.2
		endif
		if(itriton.eq.1) then
			iahp=  19
			iamp=  56
			asecp= 44.77
			idhp=  42
			idmp=  53
			dsecp= 36.2
		endif
c
	endif
c
c 20 aout 1985
c
	if(idate.eq.4.or.idate.eq.5.or.idate.eq.6) then
		iah=   18
		iam=   04
		asec=  17.2
		idh=  -22
		idm=  -18
		dsec= -06.
c
		iahe=  12
		iame=  00
		asece= 54.77
		idhe=  89
		idme=  48
		dsece= 05.9
c
		iahs=  09
		iams=  56
		asecs= 33.9
		idhs=  12
		idms=  32
		dsecs= 23.
		diss=  1.0118501*ua
c
		if(itriton.eq.0) then
			iahp=  19
			iamp=  56
			asecp= 46.92
			idhp=  42
			idmp=  53
			dsecp= 46.7
		endif
		if(itriton.eq.1) then
			iahp=  19
			iamp=  56
			asecp= 46.92
			idhp=  42
			idmp=  53
			dsecp= 46.7
		endif
c
	endif
c
c 23 avril 1986
c
	if(idate.eq.7) then
		iah=   18
		iam=   24
		asec=  49.8
		idh=  -22
		idm=  -13
		dsec= -20.
c
		iahe=  12
		iame=  00
		asece= 55.80
		idhe=  89
		idme=  47
		dsece= 52.4
c
		iahs=  02
		iams=  01
		asecs= 25.1
		idhs=  12
		idms=  21
		dsecs= 38.
		diss=  1.0053926*ua
c
		if(itriton.eq.0) then
			iahp=  19
			iamp=  56
			asecp= 48.26
			idhp=  42
			idmp=  53
			dsecp= 53.3
		endif
		if(itriton.eq.1) then
			iahp=  19
			iamp=  56
			asecp= 48.26
			idhp=  42
			idmp=  53
			dsecp= 53.3
		endif
c
	endif
c
c 23 aout 1986
c
	if(idate.eq.8) then
		iah=   18
		iam=   13
		asec=  44.4
		idh=  -22
		idm=  -19
		dsec= -20.
c
		iahe=  12
		iame=  00
		asece= 56.31
		idhe=  89
		idme=  47
		dsece= 45.7
c
		iahs=  10
		iams=  06
		asecs= 45.5
		idhs=  11
		idms=  37
		dsecs= 16.
		diss=  1.0112512*ua
c
		if(itriton.eq.0) then
			iahp=  19
			iamp=  56
			asecp= 48.92
			idhp=  42
			idmp=  53
			dsecp= 56.6
		endif
		if(itriton.eq.1) then
			iahp=  19
			iamp=  56
			asecp= 48.92
			idhp=  42
			idmp=  53
			dsecp= 56.6
		endif
c
	endif
c
c 22 juin 1987
c
	if(idate.eq.9) then
		iah=   18
		iam=   29
		asec=  25.4
		idh=  -22
		idm=  -13
		dsec= -00.
c
		iahe=  12
		iame=  00
		asece= 57.59
		idhe=  89
		idme=  47
		dsece= 29.0
c
		iahs=  06
		iams=  00
		asecs= 18.8
		idhs=  23
		idms=  26
		dsecs= 36.
		diss=  1.0163327*ua
c
		if(itriton.eq.0) then
			iahp=  19
			iamp=  56
			asecp= 50.57
			idhp=  42
			idmp=  54
			dsecp= 04.7
		endif
		if(itriton.eq.1) then
			iahp=  19
			iamp=  56
			asecp= 50.57
			idhp=  42
			idmp=  54
			dsecp= 04.7
		endif
c
	endif
c
c 09 juillet 1987
c
	if(idate.eq.10) then
		iah=   18
		iam=   27
		asec=  27.2
		idh=  -22
		idm=  -14
		dsec= -31.
c
		iahe=  12
		iame=  00
		asece= 57.66
		idhe=  89
		idme=  47
		dsece= 28.1
c
		iahs=  07
		iams=  10
		asecs= 38.8
		idhs=  22
		idms=  27
		dsecs= 00.
		diss=  1.0166832*ua
c
		if(itriton.eq.0) then
			iahp=  19
			iamp=  56
			asecp= 50.67
			idhp=  42
			idmp=  54
			dsecp= 05.2
		endif
		if(itriton.eq.1) then
			iahp=  19
			iamp=  56
			asecp= 50.67
			idhp=  42
			idmp=  54
			dsecp= 05.2
		endif
c
	endif
c
c 02 aout 1988
c
	if(idate.eq.11) then
		iah=   18
		iam=   34
		asec=  32.0
		idh=  -22
		idm=  -13
		dsec= -11.
c
		iahe=  12
		iame=  00
		asece= 59.30
		idhe=  89
		idme=  47
		dsece= 06.7
c
		iahs=  08
		iams=  49
		asecs= 27.3
		idhs=  17
		idms=  45
		dsecs= 57.
		diss=  1.0148070*ua
c
		if(itriton.eq.0) then
			iahp=  19
			iamp=  56
			asecp= 52.79
			idhp=  42
			idmp=  54
			dsecp= 15.6
		endif
		if(itriton.eq.1) then
			iahp=  19
			iamp=  56
			asecp= 52.79
			idhp=  42
			idmp=  54
			dsecp= 15.6
		endif
c
	endif
c
c 25 aout 1988
c
	if(idate.eq.12) then
		iah=   18
		iam=   32
		asec=  47.2
		idh=  -22
		idm=  -15
		dsec= -21.
c
		iahe=  12
		iame=  00
		asece= 59.40
		idhe=  89
		idme=  47
		dsece= 05.5
c
		iahs=  10
		iams=  16
		asecs= 01.7
		idhs=  10
		idms=  45
		dsecs= 37.
		diss=  1.0107324*ua
c
		if(itriton.eq.0) then
			iahp=  19
			iamp=  56
			asecp= 52.92
			idhp=  42
			idmp=  54
			dsecp= 16.2
		endif
		if(itriton.eq.1) then
			iahp=  19
			iamp=  56
			asecp= 52.92
			idhp=  42
			idmp=  54
			dsecp= 16.2
		endif
c
	endif
c
c 12 septembre 1988
c
	if(idate.eq.13.or.idate.eq.14) then
		iah=   18
		iam=   32
		asec=  08.4
		idh=  -22
		idm=  -16
		dsec= -34.
c
		iahe=  12
		iame=  00
		asece= 59.48
		idhe=  89
		idme=  47
		dsece= 04.5
c
		iahs=  11
		iams=  21
		asecs= 12.3
		idhs=  04
		idms=  10
		dsecs= 41.
		diss=  1.0064133*ua
c
		if(itriton.eq.0) then
			iahp=  19
			iamp=  56
			asecp= 53.02
			idhp=  42
			idmp=  54
			dsecp= 16.7
		endif
		if(itriton.eq.1) then
			iahp=  19
			iamp=  56
			asecp= 53.02
			idhp=  42
			idmp=  54
			dsecp= 16.7
		endif
c
	endif
c
c 07 juillet 1989
c
	if(idate.eq.15.or.idate.eq.16.or.idate.eq.17) then
		iah=   18
		iam=   46
		asec=  56.3
		idh=  -22
		idm=  -04
		dsec= -00.
c
		iahe=  12
		iame=  01
		asece= 00.73
		idhe=  89
		idme=  46
		dsece= 48.1
c
		iahs=  07
		iams=  08
		asecs= 39.3
		idhs=  22
		idms=  30
		dsecs= 17.
		diss=  1.0166867*ua
c
		if(itriton.eq.0) then
			iahp=  19
			iamp=  56
			asecp= 54.64
			idhp=  42
			idmp=  54
			dsecp= 24.7
		endif
		if(itriton.eq.1) then
			iahp=  19
			iamp=  56
			asecp= 54.64
			idhp=  42
			idmp=  54
			dsecp= 24.7
		endif
c
	endif
c
c 03 juillet 1989
c
	if( (idate.ge.18).and.(idate.le.29) ) then
		iah=   18
		iam=   45
		asec=  44.84
		idh=  -22
		idm=  -24
		dsec= -17.7
c
		iahe=  12
		iame=  01
		asece= 00.71
		idhe=  89
		idme=  46
		dsece= 48.3
c
		iahs=  06
		iams=  48
		asecs= 04.5
		idhs=  22
		idms=  58
		dsecs= 58.
		diss=  1.0167044*ua
c
		iahp=  02
		iamp=  40
		asecp= 29.08
		idhp=  83
		idmp=  29
		dsecp= 37.6
c
	endif
c
c 11 juillet 1990
c
	if(idate.eq.30) then
		iah=   18
		iam=   56
		asec=  18.7
		idh=  -21
		idm=  -55
		dsec= -59.6
c
		iahe=  12
		iame=  01
		asece= 02.28
		idhe=  89
		idme=  46
		dsece= 27.9
c
		iahs=  07
		iams=  19
		asecs= 57.1
		idhs=  22
		idms=  10
		dsecs= 14.
		diss=  1.0165690*ua
c
		if(itriton.eq.0) then
			iahp=  19
			iamp=  56
			asecp= 56.65
			idhp=  42
			idmp=  54
			dsecp= 34.6
		endif
		if(itriton.eq.1) then
			iahp=  19
			iamp=  56
			asecp= 56.65
			idhp=  42
			idmp=  54
			dsecp= 34.6
		endif
c
	endif
c
c 18 aout 1991
c
	if(idate.eq.31) then
		iah=   19
		iam=   02
		asec=  05.97
		idh=  -21
		idm=  -52
		dsec= -11.76
c
		iahe=  12
		iame=  01
		asece= 03.98
		idhe=  89
		idme=  46
		dsece= 05.8
c
		iahs=  09
		iams=  47
		asecs= 29.4
		idhs=  13
		idms=  19
		dsecs= 51.
		diss=  1.0123260*ua
c
		if(itriton.eq.0) then
			iahp=  19
			iamp=  56
			asecp= 58.85
			idhp=  42
			idmp=  54
			dsecp= 45.4
		endif
		if(itriton.eq.1) then
			iahp=  19
			iamp=  56
			asecp= 58.85
			idhp=  42
			idmp=  54
			dsecp= 45.4
		endif
c
	endif
c
c 11 juillet 1992
c
	if(idate.eq.32) then
		iah=   19
		iam=   15
		asec=  31.34
		idh=  -21
		idm=  -33
		dsec= -01.53
c
		iahe=  12
		iame=  01
		asece= 05.36
		idhe=  89
		idme=  45
		dsece= 47.8
c
		iahs=  07
		iams=  22
		asecs= 02.93
		idhs=  22
		idms=  06
		dsecs= 04.0
		diss=  1.0165949*ua
c
		if(itriton.eq.0) then
			iahp=  19
			iamp=  57
			asecp= 00.63
			idhp=  42
			idmp=  54
			dsecp= 54.2
		endif
		if(itriton.eq.1) then
			iahp=  19
			iamp=  57
			asecp= 00.63
			idhp=  42
			idmp=  54
			dsecp= 54.2
		endif
c
	endif
c
c 13 decembre 1989
c
	if(idate.eq.33 ) then
		iah=   06
		iam=   33
		asec=  27.2553
		idh=   23
		idm=   03
		dsec=  28.7457
c
		iahe=  12
		iame=  00
		asece= 00.
		idhe=  90
		idme=  00
		dsece= 00.
c
		iahs=  17
		iams=  20
		asecs= 37.03
		idhs=  -23
		idms=  -08
		dsecs= -0.8
		diss=  0.9844257*ua
c
		iahp=  17
		iamp=  52
		asecp= 09.6
		idhp=  64
		idmp=  29
		dsecp= 24.
c
	endif
c
c 18 juillet 1997 (Triton)
c
	if (idate.ge.34.and.idate.le.38) then
		iah=   20
		iam=   02
		asec=  51.2771
		idh=   -20
		idm=   -00
		dsec=  -57.157
c
		iahe=  12
		iame=  00
		asece= 00.
		idhe=  90
		idme=  00
		dsece= 00.
c
		iahs=  07
		iams=  51
		asecs= 12.
		idhs=  20
		idms=  58
		dsecs= 39.
		diss=  1.01624*ua
c
		iahp=    7
		iamp=   51
		asecp=  39.95
		idhp=  -20
		idmp=  -19
		dsecp= -13.3
c
	endif
c
c 13 mai 1971
c
	if (idate.eq.39) then
		iah=   16	! coordonnees a 18h40
		iam=   05
		asec=  27.069
		idh=   -19
		idm=   -48
		dsec=  -03.27
c
		iahe=  12
		iame=  00
		asece= 00.
		idhe=  90
		idme=  00
		dsece= 00.
c
		iahs=  03	! coordonnees APPARENTES du Soleil, pas en J2000
		iams=  19 	! ce qui doit induire une petite erreur
		asecs= 45.	! source: astronomical ephemeris
		idhs=  18
		idms=  21
		dsecs= 39.
		diss=  1.01054*ua
c
		iahp=   17 	! d'apres message de Pierre Drossart (avril 99)
		iamp=   52
		asecp=  12.
		idhp=   64
		idmp=   29
		dsecp=  24.
c
	endif
c
c 10 octobre 1999
c
	if (idate.ge.40.and.idate.le.45) then
		iah=   02	! coordonnees a 07h
		iam=   00
		asec=  23.35
		idh=   10
		idm=   37
		dsec=  28.5
c
		iahe=  12
		iame=  00
		asece= 00.
		idhe=  90
		idme=  00
		dsece= 00.
c
		iahs=  13	! coordonnees astrometriques du Soleil
		iams=  00 	! a 7h TU. Source: ephemeride JPL on line
		asecs= 58.84
		idhs=  -06
		idms=  -30
		dsecs= -12.6
		diss=  0.9987106494*ua
c
		iahp=   17 	! d'apres message de Pierre Drossart (avril 99)
		iamp=   52 	! confirme par ephemeride jpl on line
		asecp=  12.
		idhp=   64
		idmp=   29
		dsecp=  24.
c
	endif
c
c 03 decembre 1999
c
	if (idate.ge.50.and.idate.le.50) then
		iah=   02	! coordonnees a 04h30
		iam=   40
		asec=  15.61
		idh=   12
		idm=   54
		dsec=  00.4
c
		iahe=  12
		iame=  00
		asece= 00.
		idhe=  90
		idme=  00
		dsece= 00.
c
		iahs=  16	! coordonnees astrometriques du Soleil
		iams=  35 	! a 4h30 TU. Source: ephemeride JPL on line
		asecs= 49.63
		idhs=  -22
		idms=  -01
		dsecs= -47.7
		diss=  0.9858286834*ua
c
		iahp=   02	! d'apres French et al. Icarus 103, 163 (1993)
		iamp=   42	! avec taux de precession applique au pole
		asecp=  20.11	! de Saturne
		idhp=   83
		idmp=   32
		dsecp=  12.84
c
	endif
c
c 08 septembre 2001
c
        if (idate.ge.51.and.idate.le.0107) then
                iah=   21       ! coordonnees a 02h00
                iam=   38
                asec=  13.9593
                idh=  -14
                idm=  -54
                dsec= -35.675
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c
                iahs=  11		! coordonnees astrometriques du Soleil
                iams=  06		! a 2h00 TU. Source: ephemeride JPL on line
                asecs= 32.3042
                idhs=  05
                idms=  43
                dsecs= 22.322
                diss=  1.00745761587505*ua
c
                iahp=  17		! JPL Horizon
                iamp=  09		!
                asecp= 14.64	!
				idhp=  -15
                idmp=  -10
                dsecp= -30.
c
        endif
c
c 01 aout 2003
c
        if (idate.ge.0110.and.idate.le.0111) then
                iah=   22       ! coordonnees a 04h28
                iam=   15
                asec=  54.5438
                idh=  -11
                idm=  -36
                dsec= -56.020
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c
                iahs=  08 	! coordonnees astrometriques du Soleil
                iams=  43	! a 4h28 TU. Source: ephemeride JPL on line
                asecs= 41.5759
                idhs=  18
                idms=  08
                dsecs= 16.384
                diss=  1.01504892026983*ua
c
                iahp=  17	! JPL Horizon
                iamp=  09	!
                asecp= 14.64	!
				idhp=  -15
                idmp=  -10
                dsecp= -30.
c
        endif
c
c 20 juillet 2002 (Pluton)
c
        if (idate.ge.0120.and.idate.le.0130) then
                iah=   17       ! coordonnees a 01h44
                iam=   00
                asec=  18.03104
                idh=  -12
                idm=  -41
                dsec= -42.028
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c
                iahs=  07 	! coordonnees astrometriques du Soleil
                iams=  56	! a 4h30 TU. Source: ephemeride JPL on line
                asecs= 53.8203
                idhs=  20
                idms=  43
                dsecs= 29.266
                diss=  1.0161237484*ua
c
                iahp=  20	! JPL Horizon
                iamp=  52	!
                asecp= 04.8	!
		idhp=  09
                idmp=  05
                dsecp= 24.
c
        endif
c
c 21 aout 2002 (Pluton)
c
        if (idate.ge.0150.and.idate.le.0153) then
                iah=   16       ! coordonnees a 06h50 
                iam=   58
                asec=  49.43132 
                idh=  -12
                idm=  -51
                dsec= -31.869
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c
                iahs=  10 	! coordonnees astrometriques du Soleil
                iams=  00	! a 6h50 TU. Source: ephemeride JPL on line
                asecs= 48.4219 
                idhs=  12
                idms=  09
                dsecs= 30.126 
                diss=  1.0116028306*ua
c
                iahp=  20	! JPL Horizon
                iamp=  52	!
                asecp= 04.8	!
		idhp=  09
                idmp=  05
                dsecp= 24.
c
        endif
c
c 15 decembre 2003 (Tethys) 
c
        if (idate.ge.0160.and.idate.le.0173) then
                iah=   05       ! coordonnees a 19h02 
                iam=   41
                asec=  33.8412 
                idh=   22
                idm=   03
                dsec=  31.857 
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c
                iahs=  17 	! coordonnees astrometriques du Soleil
                iams=  32	! a 19h02 Source: ephemeride JPL on line
                asecs= 07.2117
                idhs= -23
                idms= -17
                dsecs=-03.175
                diss=  0.9841913083*ua
c
                iahp=  02	! JPL Horizon
                iamp=  43	!
                asecp= 29.8	!
		idhp=  83
                idmp=  32
                dsecp= 06.
c
        endif
c
c 01 Avril 2003 (Jupiter) 
c
        if (idate.ge.0200.and.idate.le.0207) then
                iah=   08       ! coordonnees de Jupiter le 1/4/2003 a 00h
                iam=   42
                asec=  45.1457 
                idh=   19
                idm=   06
                dsec=  11.474 
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c
                iahs=  00 	! coordonnees astrometriques du Soleil
                iams=  39	! a 19h02 Source: ephemeride JPL on line
                asecs= 46.6976 
                idhs=  04
                idms=  16
                dsecs= 53.989 
                diss=  0.9991281013*ua
c

		iahp=   17 	! d'apres message de Pierre Drossart (avril 99)
		iamp=   52 	! confirme par ephemeride jpl on line
		asecp=  12.
		idhp=   64
		idmp=   29
		dsecp=  24.
c
        endif
c
c 14 novembre 2003 (Titan) occultation #1
c
        if (idate.ge.0220.and.idate.le.0242) then
                iah=   06       ! coordonnees de Titan le 14/11/2003 a 00:15h
                iam=   55	! Ephem BdL DE-406 + TASS7
                asec=  20.92519 
                idh=   22
                idm=   06
                dsec=  07.0492 
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c
                iahs=  15  	! coordonnees astrometriques du Soleil
                iams=  14	! a 00h15 Source: IMCCE 
                asecs= 52.18815
                idhs= -18
                idms= -02
                dsecs=-42.0982 
                diss=  0.989486010*ua
c
		iahp=   02	! d'apres French et al. Icarus 103, 163 (1993)
		iamp=   42	! avec taux de precession applique au pole
		asecp=  19.54 	! de Saturne
		idhp=   83
		idmp=   32
		dsecp=  11.91 
c
        endif
c
c 14 novembre 2003 (Titan) occultation #2
c
        if (idate.ge.0243.and.idate.le.0245) then
                iah=   06       ! coordonnees de Titan le 14/11/2003 a 07:00h
                iam=   55	! Ephem BdL DE-406 + TASS7
                asec=  17.74810 
                idh=   22
                idm=   06
                dsec=  01.3818 
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c
                iahs=  15  	! coordonnees astrometriques du Soleil
                iams=  16	! a 07h00 Source: IMCCE 
                asecs= 01.15615
                idhs= -18
                idms= -07
                dsecs=-08.7190 
                diss=  0.989425349*ua
c
		iahp=   02	! d'apres French et al. Icarus 103, 163 (1993)
		iamp=   42	! avec taux de precession applique au pole
		asecp=  19.54 	! de Saturne
		idhp=   83
		idmp=   32
		dsecp=  11.91 
c
        endif
c
c 11 juillet 2005 (Charon) 
c
        if (idate.ge.0260.and.idate.le.0274) then
                iah=   17       ! coordonnees Charon a 03:39 (DE413 JPL)
                iam=   28       ! 
                asec=  55.0093
                idh=  -15
                idm=  -00
                dsec= -54.799
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c               
                iahs=  07       ! coordonnees astrometriques du Soleil
                iams=  21       ! a 03:39 (JPL)
                asecs= 45.9190 
                idhs=  22
                idms=  06
                dsecs= 32.336
                diss=  1.0166358653*ua
c
                iahp=  20       ! JPL Horizon
                iamp=  52       !
                asecp= 04.8     !
                idhp=  09
                idmp=  05
                dsecp= 24.
c
        endif
c
c 10 avril 2006 (Pluton) 
c
        if (idate.ge.0300.and.idate.le.0301) then
                iah=   17       ! coordonnees Pluton a 05:15 (DE413 JPL)
                iam=   46
                asec=  06.8599
                idh=  -15
                idm=  -46
                dsec= -10.702
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c               
                iahs=  01       ! coordonnees astrometriques du Soleil
                iams=  14       ! a 05:15 (JPL)
                asecs= 19.1805
                idhs=  07
                idms=  51
                dsecs= 53.301
                diss=  1.00177434505313*ua
c
c                iahp=  20       ! JPL Horizon
c                iamp=  52       !
c                asecp= 04.8     !
c                idhp=  09
c                idmp=  05
c                dsecp= 24.

                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936	!
                idhp=  06
                idmp=  10
                dsecp= 04.8

!                iahp=  20       ! Oppose de Jacobson, 11 mai 2009, pour comparaison
!                iamp=  52       ! (convention UAI)
!                asecp= 12.336   !
!                idhp=  06
!                idmp=  10
!                dsecp= 03.612
c
        endif
c
c 12 juin 2006
c
        if (idate.ge.0310.and.idate.le.0315) then	
                iah=   17       ! coordonnees Pluton a 16:20 (DE413 JPL)
                iam=   41
                asec=  12.0987
                idh=  -15
                idm=  -41
                dsec= -34.632
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c              
                iahs=  05       ! coordonnees astrometriques du Soleil
                iams=  22       ! a 16:20 (JPL)
                asecs= 53.8013
                idhs=  23
                idms=  09
                dsecs= 51.381 
                diss=  1.01546711125711*ua
c
c                iahp=  20       ! JPL Horizon
c                iamp=  52       !
c                asecp= 04.8     !
c                idhp=  09
c                idmp=  05
c                dsecp= 24.
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936    !
                idhp=  06
                idmp=  10
                dsecp= 04.8
c
        endif
c
c 06 aout 2006
c
        if (idate.ge.0320.and.idate.le.0320) then	
                iah=   17       ! coordonnees Pluton a 00:00 (DE413 JPL)
                iam=   36
                asec=  07.8042
                idh=  -15
                idm=  -49
                dsec= -29.773
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c              
                iahs=  09       ! coordonnees astrometriques du Soleil
                iams=  03       ! a 00:00 (JPL)
                asecs= 10.9051
                idhs=  16
                idms=  49
                dsecs= 06.206 
                diss=  1.01432688307509*ua
c
c                iahp=  20       ! JPL Horizon
c                iamp=  52       !
c                asecp= 04.8     !
c                idhp=  09
c                idmp=  05
c                dsecp= 24.
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936    !
                idhp=  06
                idmp=  10
                dsecp= 04.8				
c
        endif
c
c 18 mars 2007
c
        if (idate.ge.0330.and.idate.le.0349) then	
                iah=   17       ! coordonnees Pluton a 10:51 (DE413 JPL)
                iam=   55
                asec=  05.6852
                idh=  -16
                idm=  -28
                dsec= -34.188
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c              
                iahs=  23       ! coordonnees astrometriques du Soleil
                iams=  50       ! a 10:51 (JPL)
                asecs= 19.6085
                idhs= -01
                idms= -02
                dsecs=-52.994 
                diss=  0.995275421917529*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936    !
                idhp=  06
                idmp=  10
                dsecp= 04.8				
c
        endif
c
c 12 mai 2007
c
        if (idate.ge.0360.and.idate.le.0367) then	
                iah=   17       ! coordonnees Pluton a 04:00 (DE413 JPL)
                iam=   53
                asec=  32.0913
                idh=  -16
                idm=  -22
                dsec= -47.418
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c              
                iahs=  03       ! coordonnees astrometriques du Soleil
                iams=  14       ! a 04:00 (JPL)
                asecs= 09.9294
                idhs=  17
                idms=  59
                dsecs= 57.178
                diss=  1.01011703997165*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936    !
                idhp=  06
                idmp=  10
                dsecp= 04.8				
c
        endif
c
c 14 juin 2007
c
        if (idate.ge.0380.and.idate.le.0390) then	
                iah=   17       ! coordonnees Pluton a 01:27 (DE413 JPL)
                iam=   50
                asec=  20.7436
                idh=  -16
                idm=  -22
                dsec= -42.284
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c              
                iahs=  05       ! coordonnees astrometriques du Soleil
                iams=  27       ! a 01:27 (JPL)
                asecs= 33.6428 
                idhs=  23
                idms=  13
                dsecs= 44.685
                diss=  1.01565019518874*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936    !
                idhp=  06
                idmp=  10
                dsecp= 04.8				
c
        endif
c
c 09 juin 2007
c
        if (idate.ge.0391.and.idate.le.0392) then	
                iah=   17       ! coordonnees Pluton a 09:36 (DE413 JPL)
                iam=   50
                asec=  50.6582
                idh=  -16
                idm=  -22
                dsec= -29.453
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c              
                iahs=  05       ! coordonnees astrometriques du Soleil
                iams=  08       ! a 01:27 (JPL)
                asecs= 14.2849 
                idhs=  22
                idms=  54
                dsecs= 18.900
                diss=  1.01511050546129*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936    !
                idhp=  06
                idmp=  10
                dsecp= 04.8				
c
        endif
c
c 31 juillet 2007
c
        if (idate.ge.0395.and.idate.le.0397) then
                iah=   17       ! coordonnees Pluton a 13:46 (DE413 JPL)
                iam=   45
                asec=  42.0011
                idh=  -16
                idm=  -29
                dsec= -31.890
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=  08       ! coordonnees astrometriques du Soleil
                iams=  41       ! a 13:46 (JPL)
                asecs= 12.5995 
                idhs=  18
                idms=  17
                dsecs= 43.695
                diss=  1.01509370926622*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936    !
                idhp=  06
                idmp=  10
                dsecp= 04.8				
c
        endif
c
c 21 mai 2008, Triton
c
        if (idate.ge.0400.and.idate.le.0410) then
                iah=   21       ! coordonnees Triton a 01:55
                iam=   46
                asec=  11.0406
                idh=  -13
                idm=  -46
                dsec= -45.889
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=  03       ! coordonnees astrometriques du Soleil
                iams=  52       ! a 01:45 (JPL)
                asecs= 25.5594 
                idhs=  20
                idms=  12
                dsecs= 27.042
                diss=  1.01214443385841*ua
c
				iahp=   20		! pole Triton en juillet 1997 d'apres mon programme
				iamp=   01		! pole_triton
				asecp=  13.339
				idhp=   20
				idmp=   18
				dsecp=  12.154
c
        endif
c
c 22 juin 2008, Pluton
c
        if (idate.ge.0420.and.idate.le.0434) then
                iah=   17       ! coordonnees Pluton a 19:10
                iam=   58
                asec=  33.0205
                idh=  -17
                idm=  -02
                dsec= -38.630
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=  06       ! coordonnees astrometriques du Soleil
                iams=  06       ! a 19:10 (JPL)
                asecs= 58.5581 
                idhs=  23
                idms=  25
                dsecs= 42.729
                diss=  1.01636673337968*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936   !
                idhp=  06
                idmp=  10
				dsecp= 04.8
        endif
c
c 22 juin 2008, Charon
c
        if (idate.ge.0440.and.idate.le.0443) then
                iah=   17       ! coordonnees Charon a 19:20
                iam=   58
                asec=  33.0316
                idh=  -17
                idm=  -02
                dsec= -38.465
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=  06       ! coordonnees astrometriques du Soleil
                iams=  07       ! a 19:10 (JPL)
                asecs= 00.2903 
                idhs=  23
                idms=  25
                dsecs= 42.331
                diss=  1.01636715847607*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936    !
                idhp=  06
                idmp=  10
                dsecp= 04.8
        endif
c
c 24 juin 2008, Pluton
c
        if (idate.ge.0450.and.idate.le.0450) then
                iah=   17       ! coordonnees Pluton a 10:36
                iam=   58
                asec=  22.3997
                idh=  -17
                idm=  -02
                dsec= -49.265
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=  06       ! coordonnees astrometriques du Soleil
                iams=  13       ! a 10:36 (JPL)
                asecs= 48.3057 
                idhs=  23
                idms=  24
                dsecs= 01.257
                diss=  1.01646212068031*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936    !
                idhp=  06
                idmp=  10
                dsecp= 04.8
        endif
c
c 30 juin 2009, 2002 MS4
c
        if (idate.ge.0460.and.idate.le.0464) then
                iah=   18       ! coordonnees 2002 MS4 a 23:45
                iam=   03
                asec=  37.9241
                idh=  -08
                idm=  -28
                dsec= -17.732
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=  06       ! coordonnees astrometriques du Soleil
                iams=  39       ! a 10:36 (JPL)
                asecs= 55.0197 
                idhs=  23
                idms=  07
                dsecs= 14.858
                diss=  1.01664790823764d0*ua
c
                iahp=  00       ! pole inconnu
                iamp=  00
                asecp= 00.
				idhp=  00
                idmp=  00
                dsecp= 00.
        endif
c
c 25 aout 2008, Pluton
c
        if (idate.ge.0470.and.idate.le.0471) then
                iah=   17       ! coordonnees Pluton a 04:30
                iam=   53
                asec=  27.1078
                idh=  -17
                idm=  -15
                dsec= -27.434
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=  10       ! coordonnees astrometriques du Soleil
                iams=  16       ! a 04:30 (JPL)
                asecs= 51.3369 
                idhs=  10
                idms=  40
                dsecs= 49.146
                diss=  1.01076226332943*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936   !
                idhp=  06
                idmp=  10
                dsecp= 04.8
        endif
c
c 19 fevrier 2010, Varuna
c
        if (idate.ge.0480.and.idate.le.0491) then
                iah=   07       ! coordonnees Varuna a 23:10 (JPL #18)
                iam=   29
                asec=  22.4641
                idh=   26
                idm=   07
                dsec=  23.204
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=   22       ! coordonnees astrometriques du Soleil
                iams=   12       ! a 23:10 (JPL)
                asecs=  24.1657 
                idhs=  -11
                idms=  -05
                dsecs= -49.132
                diss=   0.98869525320540*ua
c
                iahp=  00       ! ARBITRAIRE !!!
                iamp=  00       !
                asecp= 00.		!
                idhp=  00
                idmp=  00
                dsecp= 00.
        endif

c
c 14 fevrier 2010, Pluton
c
        if (idate.ge.0500.and.idate.le.0502) then
                iah=   18       ! coordonnees Pluton a 04:47
                iam=   19
                asec=  14.3746
                idh=  -18
                idm=  -16
                dsec= -42.242
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=  21       ! coordonnees astrometriques du Soleil
                iams=  50       ! a 04:47 (JPL)
                asecs= 03.0532 
                idhs= -13
                idms= -06
                dsecs=-28.514
                diss=  0.98754277172922*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936   !
                idhp=  06
                idmp=  10
                dsecp= 04.8
        endif
c
c 19 mai 2010, Pluton
c
        if (idate.ge.0510.and.idate.le.0512) then				
                iah=   18       ! coordonnees Pluton a 05:56
                iam=   20
                asec=  16.7959
                idh=  -18
                idm=  -11
                dsec= -48.296
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=  03       ! coordonnees astrometriques du Soleil
                iams=  43       ! a 05:56 (JPL)
                asecs= 04.7209 
                idhs=  19
                idms=  43
                dsecs= 08.750
                diss=  1.01167805861636*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936   !
                idhp=  06
                idmp=  10
                dsecp= 04.8
        endif
c
c 04 juin 2010, Pluton
c
        if (idate.ge.0515.and.idate.le.0520) then				
                iah=   18       ! coordonnees Pluton a 15:36
                iam=   18
                asec=  47.9388
                idh=  -18
                idm=  -12
                dsec= -52.053
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=  04       ! coordonnees astrometriques du Soleil
                iams=  49       ! a 15:36 (JPL)
                asecs= 35.2112 
                idhs=  22
                idms=  27
                dsecs= 06.144
                diss=  1.01449929667757*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936   !
                idhp=  06
                idmp=  10
                dsecp= 04.8
        endif
c
c 04 juillet 2010, Pluton
c
        if (idate.ge.0530.and.idate.le.0530) then				
                iah=   18       ! coordonnees Pluton a 02:02
                iam=   15
                asec=  42.1032
                idh=  -18
                idm=  -16
                dsec= -41.407
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=  06       ! coordonnees astrometriques du Soleil
                iams=  51       ! a 02:02 (JPL)
                asecs= 37.8296 
                idhs=  22
                idms=  54
                dsecs= 27.701
                diss=  1.01668750569734*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936   !
                idhp=  06
                idmp=  10
                dsecp= 04.8
        endif
c
c 17 aout 2010, Ceres
c
        if (idate.ge.0540.and.idate.le.0544) then				
                iah=   17       ! coordonnees Ceres a 22:41
                iam=   18
                asec=  29.0165
                idh=  -27
                idm=  -26
                dsec= -38.669
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=  09       ! coordonnees astrometriques du Soleil
                iams=  48       ! a 22:41 (JPL)
                asecs= 13.8989 
                idhs=  13
                idms=  15
                dsecs= 54.692
                diss=  1.01231361900037*ua
c
                iahp=  19       ! Carry et al.  2008
                iamp=  12       ! A VERIFIER !!!
		asecp= 00.		!
                idhp=  66
                idmp=  00
                dsecp= 00.
        endif

c
c 08 jan 2011, 2003AZ84
c
        if (idate.ge.0545.and.idate.le.0546) then				
                iah=   07       ! coordonnees AZ84 = coordenada estrela
                iam=   43
                asec=  41.822
                idh=  +11
                idm=  +30
                dsec= +23.569
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=  19       ! coordonnees astrometriques du Soleil
                iams=  15       ! a 08/01/2011 06:30 (JPL)
                asecs= 42.2163 
                idhs=  -22
                idms=  17
                dsecs= 51.578
                diss=  0.98338001556303*ua
c
                iahp=  0       ! no information
                iamp=  0       ! 
		asecp= 00.	!
                idhp=  0
                idmp=  00
                dsecp= 00.
        endif

c
c 06 novembre 2010, Eris
c
        if (idate.ge.0550.and.idate.le.0560) then								
                iah=   01       ! coordonnees Eris a 02:20
                iam=   39
                asec=  09.9256
                idh=  -04
                idm=  -21
                dsec= -12.089
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=  14       ! coordonnees astrometriques du Soleil
                iams=  43       ! a 02:20 (JPL)
                asecs= 50.0058 
                idhs= -15
                idms= -51
                dsecs=-50.283
                diss=  0.99138821974729*ua
c
				iahp=  02       ! solution 1 Brown & Schaller Science 2007 Supplement
				iamp=  23		! i_E= 61.3, Omega_E= 139 ecliptique --->
				asecp= 36.		! l_E= Omega-90 = 49, beta_E= 90-i= 28.7 ecliptique --->  
				idhp=  44		! alpha= 35.9 deg, delta= 44.7 deg par coord_eclipt_to_equat.f
				idmp=  42
				dsecp= 00.
!				iahp=  00       ! solution 2 Brown & Schaller Science 2007 Supplement
!				iamp=  40.		! idem mais avec i_E= 142, Omega_E= 68
!				asecp= 24.
!				idhp= -54
!				idmp= -36
!				dsecp=-00.
        endif

c
c 23 avril 2011, Makemake
c
        if (idate.ge.0570.and.idate.le.0575) then								
                iah=   12       ! coordonnees Makemake a 01:36
                iam=   36
                asec=  11.4015
                idh=  +28
                idm=   11
                dsec=  10.471
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=  02       ! coordonnees astrometriques du Soleil
                iams=  00       ! a 01:36 (JPL)
                asecs= 54.4039 
                idhs= +12
                idms=  18
                dsecs= 40.703
                diss=  1.00526965255930*ua
c
				iahp=  00       ! INCONNU !
				iamp=  00
				asecp= 00. 
				idhp=  00
				idmp=  00
				dsecp= 00.
        endif
c
c 04 mai 2011, Quaoar
c
        if (idate.ge.0584.and.idate.le.0599) then								
                iah=   17       ! coordonnees Quaoar a 02:40
                iam=   28
                asec=  50.8132
                idh=  -15
                idm=  -27
                dsec= -42.760
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=  02       ! coordonnees astrometriques du Soleil
                iams=  42       ! a 01:36 (JPL)
                asecs= 49.3214 
                idhs= +15
                idms=  47
                dsecs= 13.676
                diss=  1.00819571985839*ua
c
				iahp=  00       ! INCONNU !
				iamp=  00
				asecp= 00. 
				idhp=  00
				idmp=  00
				dsecp= 00.
        endif
c
c 04 juin 2011, Pluton
c
        if (idate.ge.0600.and.idate.le.0603) then	
                iah=   18       ! coordonnees Pluton a 05:42
                iam=   27
                asec=  53.8241
                idh=  -18
                idm=  -45
                dsec= -30.902
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c
                iahs=  04       ! coordonnees astrometriques du Soleil
                iams=  46       ! a 05:42 (JPL)
                asecs= 51.7794 
                idhs= +22
                idms=  22
                dsecs= 24.813
                diss=  1.01444225284419*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936   !
                idhp=  06
                idmp=  10
                dsecp= 04.8
        endif
c
c 04 juin 2011, Charon
c
        if (idate.ge.0610.and.idate.le.0613) then		
                iah=   18       ! coordonnees Charon a 05:42
                iam=   27
                asec=  53.7777
                idh=  -18
                idm=  -45
                dsec= -30.970
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c
                iahs=  04       ! coordonnees astrometriques du Soleil
                iams=  46       ! a 05:42 (JPL)
                asecs= 51.7794 
                idhs= +22
                idms=  22
                dsecs= 24.813
                diss=  1.01444225284419*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936   !
                idhp=  06
                idmp=  10
                dsecp= 04.8
        endif
c
c 23 juin 2011, Pluton
c
        if (idate.ge.0620.and.idate.le.0629) then	
                iah=   18       ! coordonnees Pluton a 11:15
                iam=   25
                asec=  55.52406
                idh=  -18
                idm=  -48
                dsec= -07.0068
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c
                iahs=  06       ! coordonnees astrometriques du Soleil
                iams=  06       ! a 11:15 (JPL)
                asecs= 34.7350 
                idhs= +23
                idms=  25
                dsecs= 45.309
                diss=  1.01636831112548*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936   !
                idhp=  06
                idmp=  10
                dsecp= 04.8
        endif
c
c 23 juin 2011, Charon
c
        if (idate.ge.0630.and.idate.le.0639) then	
                iah=   18       ! coordonnees Charon a 11:15
                iam=   25
                asec=  55.47660
                idh=  -18
                idm=  -48
                dsec= -7.0182
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c
                iahs=  06       ! coordonnees astrometriques du Soleil
                iams=  06       ! a 11:15 (JPL)
                asecs= 34.7350 
                idhs= +23
                idms=  25
                dsecs= 45.309
                diss=  1.01636831112548*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936   !
                idhp=  06
                idmp=  10
                dsecp= 04.8
        endif
c
c 27 juin 2011, Pluton
c
        if (idate.ge.0640.and.idate.le.0649) then	
                iah=   18       ! coordonnees Pluton a 14:55
                iam=   25
                asec=  28.85349
                idh=  -18
                idm=  -48
                dsec= -47.9075
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c
                iahs=  06       ! coordonnees astrometriques du Soleil
                iams=  23       ! a 14:55 (JPL)
                asecs= 50.1697 
                idhs= +23
                idms=  19
                dsecs= 28.975
                diss=  1.01659073275910*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936   !
                idhp=  06
                idmp=  10
                dsecp= 04.8
        endif
c
c 27 juin 2011, Hydra
c
        if (idate.ge.0650.and.idate.le.0659) then	
                iah=   18       ! coordonnees Pluton a 14:55
                iam=   25
                asec=  29.00941
                idh=  -18
                idm=  -48
                dsec= -47.7717
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c
                iahs=  06       ! coordonnees astrometriques du Soleil
                iams=  23       ! a 14:55 (JPL)
                asecs= 50.1697 
                idhs= +23
                idms=  19
                dsecs= 28.975
                diss=  1.01659073275910*ua
c
                iahp=  20       ! Oppose de Buie et al, AJ 2008
                iamp=  52       ! (convention UAI)
                asecp= 12.936   !
                idhp=  06
                idmp=  10
                dsecp= 04.8
        endif

c
c 03 fevrier 2012, 2003 AZ84
c

        if (idate.ge.0660.and.idate.le.0679) then								
                iah=   07       ! coordonnees 20003 AZ84 19:47
                iam=   45
                asec=  54.7845
                idh=   11
                idm=   12
                dsec=  43.062
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
                iahs=  21       ! coordonnees astrometriques du Soleil
                iams=  06       ! a 19:47 (JPL)
                asecs= 43.5084 
                idhs= -16
                idms= -33
                dsecs= -48.759
                diss=  0.98562859759213*ua
c
				iahp=  00       ! INCONNU !
				iamp=  00
				asecp= 00. 
				idhp=  00
				idmp=  00
				dsecp= 00.
        endif

c
c 17 fev 2012, Quaoar
c

        if (idate.ge.0680.and.idate.le.0683) then								
                iah=   17       ! coordonnees Quaoar a 04:30
                iam=   34
                asec=  21.8488
                idh=  -15
                idm=  -42
                dsec= -10.339
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c 
c 21 59 41.1093 -12 15 30.421
                iahs=  21       ! coordonnees astrometriques du Soleil
                iams=  59       ! a 04:30 (JPL)
                asecs= 41.1093 
                idhs= -12
                idms=  -15
                dsecs= -30.421
                diss=  0.98802203390966*ua
c
				iahp=  00       ! INCONNU !
				iamp=  00
				asecp= 00. 
				idhp=  00
				idmp=  00
				dsecp= 00.
        endif

c
c 03 December 2015, 2002 VE95
c

        if (idate.ge.0685.and.idate.le.0687) then								
                iah=   5       ! coordonnees 2002 VE95 12:30
                iam=   37
                asec=  38.5129
                idh=   +8
                idm=   +5
                dsec=  +9.416
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c   
                iahs=  16       ! coordonnees astrometriques du Soleil
                iams=  36       ! a 12:30 (JPL)
                asecs= 49.84 
                idhs= -22
                idms=  -3
                dsecs= -41.6
                diss=  0.98578398180764*ua
c
				iahp=  00       ! INCONNU !
				iamp=  00
				asecp= 00. 
				idhp=  00
				idmp=  00
				dsecp= 00.
        endif

c
c Quaoar-July 09, 2013
c

        if (idate.ge.0688.and.idate.le.0688) then								
                iah=   17       ! coordonnees Quaoar et 02:40 (JPL#31). 
                iam=   34
                asec=  40.4554
                idh=   -15
                idm=   -23
                dsec=  -37.581
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c   
                iahs=  7       ! coordonnees astrometriques du Soleil et 2:40 (JPL). 
                iams=  13      
                asecs= 13.01 
                idhs=  +22
                idms=  +22
                dsecs= +15.7
                diss=  1.01667215105169*ua
c
				iahp=  00       ! INCONNU !
				iamp=  00
				asecp= 00. 
				idhp=  00
				idmp=  00
				dsecp= 00.
        endif

c
c Sedna - January 13, 2013
c

        if (idate.ge.0689.and.idate.le.0690) then								
                iah=   3       ! coordonnees Sedna et 12:30 (JPL#28). 
                iam=   32
                asec=  44.7011
                idh=   +6
                idm=   +47
                dsec=  +18.859
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c   
                iahs=  19       ! coordonnees astrometriques du Soleil et 12:30 (JPL). 
                iams=  40      
                asecs= 37.75 
                idhs=  -21
                idms=  -25
                dsecs= -28.8
                diss=  0.98360230336624*ua
c
				iahp=  00       ! INCONNU !
				iamp=  00
				asecp= 00. 
				idhp=  00
				idmp=  00
				dsecp= 00.
        endif

c
c Eris - August 29, 2013
c

        if (idate.ge.0691.and.idate.le.0692) then								
                iah=   1       ! coordonnees Eris et 15:38 (JPL#57). 
                iam=   42
                asec=  58.5021
                idh=   -3
                idm=   -21
                dsec=  -17.427
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c   
                iahs=  10       ! coordonnees astrometriques du Soleil et 15:38 (JPL). 
                iams=  32      
                asecs= 07.93 
                idhs=  +9
                idms=  +12
                dsecs= +40.3
                diss=  1.00979124847261*ua
c
				iahp=  00       ! INCONNU !
				iamp=  00
				asecp= 00. 
				idhp=  00
				idmp=  00
				dsecp= 00.
        endif

c
c 2003 VS2 in November 07, 2014
c

        if (idate.ge.0693.and.idate.le.0698) then								
                iah=   4       ! coordonnees 2003 VS2 et 04:25 (JPL#30). 
                iam=   48
                asec=  32.1383
                idh=   33
                idm=   58
                dsec=  36.214


                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c   
                iahs=  14       ! coordonnees astrometriques du Soleil et 04:25 (JPL).  
                iams=  48      
                asecs= 03.43 
                idhs=  -16
                idms=  -10
                dsecs= -49.7
                diss=  0.99106899191144*ua
c
				iahp=  00       ! INCONNU !
				iamp=  00
				asecp= 00. 
				idhp=  00
				idmp=  00
				dsecp= 00.
        endif










        if (idate.ge.0700.and.idate.le.703) then
c 2014-Nov-07 04:25     04 48 32.1383 +33 58 36.214 35.7157838684412 -14.4056320 297.039117
                iah= 04
                iam= 48
                asec= 32.1383
                idh= +33
                idm= 58
                dsec= 36.214
c
                iahe=  12
                iame=  00
                asece= 00.
                idhe=  90
                idme=  00
                dsece= 00.
c   
                iahs= 14
                iams= 48
                asecs= 03.4262
                idhs= -16
                idms= -10
                dsecs= -49.704
                diss= 0.99106899191144*ua
c
		iahp=  00       ! INCONNU !
		iamp=  00
		asecp= 00. 
		idhp=  00
		idmp=  00
		dsecp= 00.
        endif













c
c
c Angles en radian
c
	alpha= ( dfloat(iah)*3600.d0 + dfloat(iam)*60. + asec )/3600.d0
	alpha= alpha*(pi/12.)
	delta= ( dfloat(idh)*3600.d0 + dfloat(idm)*60. + dsec )/3600.d0
	delta= delta*(pi/180.)
c
	alphap= ( dfloat(iahp)*3600.d0 + dfloat(iamp)*60. + asecp )/3600.d0
	alphap= alphap*(pi/12.)
	deltap= ( dfloat(idhp)*3600.d0 + dfloat(idmp)*60. + dsecp )/3600.d0
	deltap= deltap*(pi/180.)
c
	alphae= ( dfloat(iahe)*3600.d0 + dfloat(iame)*60. + asece )/3600.d0
	alphae= alphae*(pi/12.)
	deltae= ( dfloat(idhe)*3600.d0 + dfloat(idme)*60. + dsece )/3600.d0
	deltae= deltae*(pi/180.)
c
	alphas= ( dfloat(iahs)*3600.d0 + dfloat(iams)*60. + asecs )/3600.d0
	alphas= alphas*(pi/12.)
	deltas= ( dfloat(idhs)*3600.d0 + dfloat(idms)*60. + dsecs )/3600.d0
	deltas= deltas*(pi/180.)
c
c Quantites intermediaires
c
c cos et sin des angles
c
	cosal= dcos(alpha)
	sinal= dsin(alpha)
	cosd=  dcos(delta)
	sind=  dsin(delta)
c
	cosalp= dcos(alphap)
	sinalp= dsin(alphap)
	cosdp= dcos(deltap)
	sindp= dsin(deltap)
c
	cosale= dcos(alphae)
	sinale= dsin(alphae)
	cosde=  dcos(deltae)
	sinde=  dsin(deltae)
c
	cosals= dcos(alphas)
	sinals= dsin(alphas)
	cosds=  dcos(deltas)
	sinds=  dsin(deltas)
c
	cosa=  dcos(alphap-alpha)
	sina=  dsin(alphap-alpha)
c
c
c
	sinB= -( sindp*sind + cosdp*cosd*cosa )
	cBcP= sindp*cosd - cosdp*sind*cosa
	cBsP= cosdp*sina
c
	B= dasin(sinB)
	cosB= dcos(B)
	sBcB= sinB*cosB
	P= datan(cBsP/cBcP)
	if(cBcP.lt.0.) P= P + pi
	sinP= dsin(P)
	cosP= dcos(P)
c
c coordonnees des vecteurs
c
c est
c
	ex= -sinal
	ey=  cosal
	ez=  0.
c
c nord
c
	nx= -cosal*sind
	ny= -sinal*sind
	nz=        cosd
c
	cx=  cosal*cosd
	cy=  sinal*cosd
	cz=        sind
c
	ax= -sinB*sinP*ex - sinB*cosP*nx - cosB*cx
	ay= -sinB*sinP*ey - sinB*cosP*ny - cosB*cy
	az= -sinB*sinP*ez - sinB*cosP*nz - cosB*cz
c
	bx=      -cosP*ex +      sinP*nx
	by=      -cosP*ey +      sinP*ny
	bz=      -cosP*ez +      sinP*nz
c
c pole de la planete
c
	px=  cosB*sinP*ex + cosB*cosP*nx - sinB*cx
	py=  cosB*sinP*ey + cosB*cosP*ny - sinB*cy
	pz=  cosB*sinP*ez + cosB*cosP*nz - sinB*cz
c
c pole de la Terre
c
	ptx=  cosale*cosde
	pty=  sinale*cosde
	ptz=         sinde
c
c noeud ascendant (intersection plan equatorial planete et plan equatorial Terre)
c
	n1x= pty*pz - ptz*py
	n1y= ptz*px - ptx*pz
	n1z= ptx*py - pty*px
	n1= dsqrt( n1x*n1x + n1y*n1y + n1z*n1z )
	n1x= n1x/n1
	n1y= n1y/n1
	n1z= n1z/n1
c
c p,n1,n2 triedre direct
c
	n2x= py*n1z - pz*n1y
	n2y= pz*n1x - px*n1z
	n2z= px*n1y - py*n1x
c
c Calcul des coefficients de projection le long de a et b
c
	erp= cosdp*sina
	nrp= sindp*cosd - cosdp*sind*cosa
	aa=  -erp/sBcB
	bb=  -nrp/sBcB
	cc=  -nrp/cosB
	dd=   erp/cosB
c
c Calcul des coefficients qui permettent de calculer la position du point
c par rapport au noeud ascendant (intersection avec plan equatorial Terre)
c
	aap= ax*n1x + ay*n1y + az*n1z
	bbp= bx*n1x + by*n1y + bz*n1z
	ccp= ax*n2x + ay*n2y + az*n2z
	ddp= bx*n2x + by*n2y + bz*n2z
c
c calcul de l'elevation planetocentrique du soleil
c
	sx= diss*cosals*cosds
	sy= diss*sinals*cosds
	sz= 	   diss*sinds
	plx=   dis*cosal*cosd
	ply=   dis*sinal*cosd
	plz=         dis*sind
	ux= plx - sx
	uy= ply - sy
	uz= plz - sz
	anorm= dsqrt( ux**2 + uy**2 + uz**2 )
	ux= ux/anorm
	uy= uy/anorm
	uz= uz/anorm
	sinBP= -( ux*px + uy*py + uz*pz )
	BP= dasin(sinBP)
c
c
c
	return
	end
c
c
c				FIN DE POLPO
c
c **************************************************************************
c
c
c Subroutine qui calcule la longitude des arcs L,E,F au jour julien "ajd",
c a l'heure (AU NIVEAU DE NEPTUNE) "heure".
c
c Sorties: lib1, lib2 (longitudes du debut et de la fin de Liberte)
c	   ega1, ega2 (                "		   Egalite)
c	   fra1, fra2 (                "		   Fraternite)
c
c Ces longitudes (deg) sont calculees dans le plan equatorial de Neptune, a
c partir du noeud ascendant 1950. EME (cf. memo. de Phil Nicholson du 4/4/1990)
c
	subroutine arcpos (ajd,heure, lib1,lib2,ega1,ega2,fra1,fra2)
c
	implicit real*8 (a-h,o-z)
	real*8 lib1, lib2, lib10, lib20, moy
c
c Moyen mouvement des arcs (deg/jour)
c
	moy= 820.1185
c
c Longitudes des arcs a l'epoque de reference (deg)
c
	lib10= 347.4
	lib20= 342.9
	ega10= 332.7
	ega20= 328.7
	fra10= 323.7
	fra20= 313.7
c
c L'epoque de reference t0 est le 25 aout 1989 a 4 h TU AU NIVEAU DE NEPTUNE
c
	t0= 2447763.5 + 4.d+00/24.d+00
c
c L'heure "heure" est DEJA corrigee du temps de lumiere
c
	t= (ajd-t0) + heure/24.d+00
c
	lib1= moy*t + lib10
	lib2= moy*t + lib20
	ega1= moy*t + ega10
	ega2= moy*t + ega20
	fra1= moy*t + fra10
	fra2= moy*t + fra20
c
	lib1= dmod(lib1,360.d+00)
	if (lib1.lt.0.) lib1= lib1 + 360.
	lib2= dmod(lib2,360.d+00)
	if (lib2.lt.0.) lib2= lib2 + 360.
c
	ega1= dmod(ega1,360.d+00)
	if (ega1.lt.0.) ega1= ega1 + 360.
	ega2= dmod(ega2,360.d+00)
	if (ega2.lt.0.) ega2= ega2 + 360.
c
	fra1= dmod(fra1,360.d+00)
	if (fra1.lt.0.) fra1= fra1 + 360.
	fra2= dmod(fra2,360.d+00)
	if (fra2.lt.0.) fra2= fra2 + 360.
c
c
c
	return
	end
c
c
c 				FIN DE ARCPOS
c
c **************************************************************************
c
c Sous-programme qui calcule la position des arcs dans le plan du ciel,
c connaissant leur position dans le plan des anneaux. Les formules sont les
c inverses de celles qui donnent les positions dans le plan des anneaux en
c fonction de ksi et eta
c
c entrees: aa,bb,cc,dd,aap,bbp,ccp,ddp, coefficients calcules dans polpo
c          along, longitude dans le plan des anneaux, en degre (noeud EME 1950)
c
c sorties: aksi, eta, coordonnees dans le plan du ciel
c
	subroutine arcposciel (aa,bb,cc,dd,aap,bbp,ccp,ddp,along, aksi, eta)
c
c
	implicit real*8 (a-h,o-z)
	pi= 4.d+00*datan(1.d+00)
	rad= pi/180.d+00
c
c rayon suppose des arcs
c
	ray= 62942.
c
	det = aa*dd   - bb*cc
	detp= aap*ddp - bbp*ccp

	if ((det.eq.0.).or.(detp.eq.0.)) then
		return
	endif
c
	xp= ray*dcos(rad*along)
	yp= ray*dsin(rad*along)
c
	x= ( xp*ddp - yp*bbp)/detp
	y= (-xp*ccp + yp*aap)/detp
c
	aksi= ( x*dd - y*bb)/det
	eta=  (-x*cc + y*aa)/det
c
c
c
	return
	end
c
c 				FIN DE ARCPOSCIEL
c *********************************************************************
c
c
c
c Sous-programme de calcul du Jour Julien.
c
c	Entree: heure, jj/mm/aa (date).
c	Sortie: ajd (jour julien),
c	hs0 (temps sideral a 00 T.U. du jour).
c	hs  (temps sideral courant)

	subroutine ajulien (heure,jj,mm,aa, ajd,hs0,hs)

	implicit real*8 (a-h,o-z)
	dimension nj(12)
	integer*4 aa
c
c Nombre de jours dans chaque mois
c
	nj(1)=  31
	nj(2)=  28
	if( mod(aa,4  ).eq.0 ) nj(2)=  29 ! annees bissextiles normales
        if( mod(aa,100).eq.0 ) nj(2)=  28 ! les siecles ne sont pas bissextiles
        if( mod(aa,400).eq.0 ) nj(2)=  29 ! les 4-siecles sont bissextiles
	nj(3)=  31
	nj(4)=  30
	nj(5)=  31
	nj(6)=  30
	nj(7)=  31
	nj(8)=  31
	nj(9)=  30
	nj(10)= 31
	nj(11)= 30
	nj(12)= 31
c
c Nombre de jour ecoules depuis le 1 er janvier 1950 (00 h T.U.) au 1 er 
c janvier de l'annee "aa" (00 h T.U.).
c
	jd= 0
c
c Annee 1950
c
	if (aa.eq.1950) jd= 0
c
c Annees apres 1950
c
	if (aa.gt.1950) then 
	 do i= 1950, aa-1
	  nbrj= 365
	  if( mod(i,4  ).eq.0 ) nbrj= 366 ! annees bissextiles normales
	  if( mod(i,100).eq.0 ) nbrj= 365 ! les siecles ne sont pas bissextiles
	  if( mod(i,400).eq.0 ) nbrj= 366 ! les 4-siecles sont bissextiles
	  jd= jd + nbrj
	 enddo
	endif
c
c Annees avant 1950
c
	if (aa.lt.1950) then 
	 do i= aa, 1949
	  nbrj= -365
	  if( mod(i,4  ).eq.0 ) nbrj= -366 ! annees bissextiles normales
	  if( mod(i,100).eq.0 ) nbrj= -365 ! les siecles ne sont pas bissextiles
	  if( mod(i,400).eq.0 ) nbrj= -366 ! les 4-siecles sont bissextiles
	  jd= jd + nbrj
	 enddo
	endif
c
c Nombre de jours ecoules depuis le 1 er janvier 1950 (00 h T.U.) jusqu'au
c jour "jj mm aa" (00 h T.U.).
c
	do 2 i= 1, mm-1
 2	jd= jd + nj(i)
c
	jd= jd + jj - 1
c
c Jour Julien au jour "jj mm aa" a 00 h T.U.
c
	ajd= dfloat(jd) + 2433282.5d0
c
c Heure siderale au jour "jj mm aa" a 00 h T.U.
c On prend pour reference l'heure siderale au 1 er janvier 1985 a 00 h T.U.
c ---> changee pour le 1er janvier 2000 a 00h TU le 7 aout 2006 
c ---> changee pour le 1er janvier 2010 a 00h TU le 6 novembre 2010
c
	rapp= 0.99726956633602776d0
c	ihhs0= 6		! 01 janvier 1985, from JPL Horizons
c	mms0=  42
c	secs0= 21.9674

c	ihhs0= 6		! 01 janvier 2000
c	mms0=  39
c	secs0= 51.8041d0

	ihhs0= 6		! 01 janvier 2010
	mms0=  42
	secs0= 10.1622d0

	ts0=(3600.d0*dfloat(ihhs0) + 60.d0*dfloat(mms0) + secs0)/3600.d0
c
c hs0
c
c	djd= ajd - 2446066.5d0			! origine 01 janvier 1985
c	djd= ajd - 2451544.5d0			! origine 01 janvier 2000
	djd= ajd - 2455197.5d0			! origine 01 janvier 2010

	hs= djd/rapp
	unite= 1.d0
	hs= dmod(hs,unite)
	hs= hs*24.d0 + ts0
	unite= 24.d0
	hs0= dmod(hs,unite)
	if (hs0.lt.(0.d0)) hs0= hs0 + 24.d0
c
c hs
c
	hs= hs0 + heure/rapp
	hs= dmod(hs,unite)
	if (hs.lt.(0.d0)) hs= hs + 24.d0
c
c
c
	return
	end
c
c
c			FIN DE AJULIEN
c
c
c **********************************************************************
c
c sous-programme qui donne la date jj,mm,aa en fonction de l'observation
c iobs. Sert ensuite a calculer le jour julien.
c
	subroutine jjmmaa (iobs, jj, mm, aa)
c
	implicit real*8 (a-h,o-z)
	integer*4 aa
c
c 15 juin 1983
c
	if(iobs.eq.0.or.iobs.eq.1) then
		jj= 15
		mm= 06
		aa= 1983
	endif
c
c 22 juillet 1984
c
	if(iobs.eq.2.or.iobs.eq.3) then
		jj= 22
		mm= 07
		aa= 1984
	endif
c
c 20 aout 1985
c
	if(iobs.eq.4.or.iobs.eq.5.or.iobs.eq.6) then
		jj= 20
		mm= 08
		aa= 1985
	endif
c
c 23 avril 1986
c
	if(iobs.eq.7) then
		jj= 23
		mm= 04
		aa= 1986
	endif
c
c 23 aout 1986
c
	if(iobs.eq.8) then
		jj= 23
		mm= 08
		aa= 1986
	endif
c
c 22 juin 1987
c
	if(iobs.eq.9) then
		jj= 22
		mm= 06
		aa= 1987
	endif
c
c 09 juillet 1987
c
	if(iobs.eq.10) then
		jj= 09
		mm= 07
		aa= 1987
	endif
c
c 02 aout 1988
c
	if(iobs.eq.11) then
		jj= 02
		mm= 08
		aa= 1988
	endif
c
c 25 aout 1988
c
	if(iobs.eq.12) then
		jj= 25
		mm= 08
		aa= 1988
	endif
c
c 12 septembre 1988
c
	if(iobs.eq.13.or.iobs.eq.14) then
		jj= 12
		mm= 09
		aa= 1988
	endif
c
c 07 juillet 1989
c
	if(iobs.eq.15.or.iobs.eq.16.or.iobs.eq.17) then
		jj= 07
		mm= 07
		aa= 1989
	endif
c
c 03 juillet 1989
c
	if(iobs.ge.18.and.iobs.le.29) then
		jj= 03
		mm= 07
		aa= 1989
	endif
c
c 11 juillet 1990
c
	if(iobs.eq.30) then
		jj= 11
		mm= 07
		aa= 1990
	endif
c
c 18 aout 1991
c
	if(iobs.eq.31) then
		jj= 18
		mm= 08
		aa= 1991
	endif
c
c 11 juillet 1992
c
	if(iobs.eq.32) then
		jj= 11
		mm= 07
		aa= 1992
	endif
c
c 13 decembre 1989
c
	if(iobs.eq.33) then
		jj= 13
		mm= 12
		aa= 1989
	endif
c
c 18 juillet 1997
c
	if(iobs.ge.34.and.iobs.le.38) then
		jj= 18
		mm= 07
		aa= 1997
	endif
c
c 13 mai 1971
c
	if(iobs.eq.39) then
		jj= 13
		mm= 05
		aa= 1971
	endif
c
c 10 octobre 1999
c
	if(iobs.ge.40.and.iobs.le.45) then
		jj= 10
		mm= 10
		aa= 1999
	endif
c
c 03 decembre 1999
c
	if(iobs.eq.50) then
		jj= 03
		mm= 12
		aa= 1999
	endif
c
c
c
	return
	end
c
c				FIN DE JJMMAA
c
c ***********************************************************************
c
	subroutine latitude (P,B,H,diskm,heure,ksi,eta,vksi,veta)
	implicit real*8 (a-h,o-z)
	real*8 J_2, K, ksi, mu
	pi= dacos(-1.d0)
	c= 299792.458d+00				     ! NB: km/sec
c
c Parametres de Jupiter
c
	GM     = 1.2669d+17 	! m^3/sec^2
	R_g    = 71398.d0	! rayon GRAVITATIONEL equatorial
	J_2    = 14736.d-6	! R_g et J_2 pris dans livre Murray & Dermott
	mu     = 3.68d-27	! masse moleculaire (kg) pour 90% H2 et 10% He
	K      = 4.73d-30	! refractivite moleculaire (m^3/molecule)
	dis    = diskm*1.d3	! distance (m)
c	r_eq   = 71880.d0	! rayon equatorial de MI-OCCULTATION (Hubbard et al. 72)
c	epsilon= 0.065		! ???
	r_eq   = 71698 + 30. 	! rayon equatorial de MI-OCCULTATION
				! au niveau de la planete, SANS effet relativiste
				! on doit donc mettre ici Req dans l'ombre + H
	epsilon= 0.06406	! aplatissement REEL (voir 26 mars 2001)
c
c Parametres de Saturne
c
c	GM     = 3.7931272d+16	! m^3/sec^2
c	R_g    = 60330.d0	! rayon gravitationel equatorial
c	J_2    = 16298.d-6	! R_g et J_2 pris dans livre Murray & Dermott
c	mu     = 3.54d-27    	! masse moleculaire (kg), Hubbard et al. Ic.130, 404 (9`
c	K      = 4.73d-30    	! refractivite moleculaire (m^3/molecule)
c	dis    = diskm*1.d3	! distance (m)
c	r_eq   = 60960.d0	! rayon equatorial de MI-OCCULTATION (au niveau de la planete)
c				! d'apres Hubbard et al. Icarus 130, 404 (97)
c	epsilon= 0.097962434	! d'apres site JPL
c
c
c Parametres de Titan
c
c	GM     = 8.977d+12	! m^3/sec^2
c	R_g    = 2575.d0	! rayon gravitationel equatorial (sol de Titan)
c	J_2    = 0.d0		!
c	mu     = 4.65174d-26	! masse moleculaire de l'azote
c	K      = 1.12d-29	! refractivite moleculaire azote (m^3/molecule)
c	dis    = diskm*1.d3	! distance (m)
c	r_eq   = 3000.		! rayon de miocc au niveau de Titan
c	epsilon= 0.d0		! aplatissement
c
        r_po= r_eq*( 1.d0 - epsilon )
	epsilonp = 1.d0 - dsqrt( (dsin(B))**2 + ((1.d0-epsilon)*dcos(B))**2 )
	write(*,*)  'Aplatissement et rayon polaire apparents:',
     *          epsilonp, r_eq*(1.d0 - epsilonp)
c
c On calcule la latitude planetocentrique par iteration: on estime d'abord
c le rayon minimum d'approche rmin (Baum and Code), puis on estime la latitude
c de ce point, puis le rayon de mi-hauteur a cette latitude (limbe elliptique)
c puis de nouveau rmin, etc...
c
	phi0= 60.*(pi/180.d0)
	tol= 1.d-4
 1	r_limbe= 1.d0/dsqrt( (dcos(phi0)/r_eq)**2 + (dsin(phi0)/r_po)**2 )
	z= dsqrt(ksi*ksi+eta*eta)

	call baum_code (0.d0,H,(r_limbe-H),z, flux,rmin)

	call image (P,B,GM,c,R_g,J_2,diskm,r_eq,epsilonp,rmin,ksi,eta,
     *              vksi,veta,vrel,xlimbe,ylimbe,xlimber,ylimber)

	sinphi= dcos(B)*( dsin(P)*xlimbe + dcos(P)*ylimbe )/rmin
	phi= dasin(sinphi)
	cosphi= dcos(phi)
	if (dabs(phi-phi0).ge.tol) then
		phi0= phi
		go to 1
	endif
c
c Pression a la miocc. (mubar), puis pression courante (mubar)
c Attention: H et r_limbe doivent etre exprimes dans la meme unite! (ici km)
c NB: le 10. en facteur transforme Pa ---> mubar
c
	pres2= 10.d0*(mu*GM/(K*dis))*dsqrt( H**5/( 2.d0*pi*(r_limbe**5) ))
	pres = pres2*dexp( -(rmin-r_limbe)/H )
c
c Densite moleculaire a miocc et a r_min (molecules/m^3)
c
	dens2= ( (H/diskm)*dsqrt(H/(2.d0*pi*r_limbe)) )/K
	dens = dens2*dexp( -(rmin-r_limbe)/H )
c
c Densite moleculaire surfacique (molecule/m^2)
c NB: 1.d6 transforme H**2 de km**2 ---> m**2
c
	sigma= H**2/(K*dis)*dsqrt(rmin/r_limbe)*dexp(-(rmin-r_limbe)/H)
	sigma= 1.d6*sigma
c
c Latitude geographique
c
	tanphip= ((r_eq/r_po)**2)*(sinphi/cosphi)
	phip= datan(tanphip)
c
c Longitude a la surface de Jupiter
c
	cosl=          (xlimbe*dcos(P) - ylimbe*dsin(P))/(dcos(phi)*r_limbe)
	sinl= -dsin(B)*(xlimbe*dsin(P) + ylimbe*dcos(P))/(dcos(phi)*r_limbe)
	alplan= datan(sinl/cosl)
	if (cosl.lt.0.d0) alplan= alplan + pi
	if (alplan.lt.0.d0) alplan= alplan + 2.d0*pi
c
c Correction relativiste. ATTENTION! pas encore de correction due a J2
c Noter les facteurs de conversion de GM et dis (km^3/sec^2 et km)
c
	drelat= (4.d0*(GM/1.d9)*(dis/1.d3))/((c**2)*rmin)
c
c
c
	write (30,*) sngl(ksi), sngl(eta), sngl(heure*3600.d0)
	write (31,*) sngl(xlimber), sngl(ylimber)
	write (32,*) sngl(heure*3600.d0), sngl(phi*(180./pi))
	write (33,*) sngl(heure*3600.d0), sngl(phip*(180./pi))
	xdumb= dcos(alplan)*dcos(phi)
	ydumb= dsin(alplan)*dcos(phi)
	write (34,*) sngl(xdumb), sngl(ydumb)
	xdumb= dcos(alplan)*dcos(phip)
	ydumb= dsin(alplan)*dcos(phip)
	write (35,*) sngl(xdumb), sngl(ydumb)
c
c Attention! Le rayon de plus courte approche doit tenir
c compte des effets relativistes (~ 50 km), d'ou le calcul
c de sqrt(xlimber**2 + ylimber**2) ci-dessous. Cette quantite
c est appelee rlimber un peu plus bas.
c
       write (36,*) sngl(heure*3600.d0), sngl(phi*(180./pi)),
     *               sngl(alplan*(180./pi) -90.),
     *               sngl(dsqrt(xlimber**2+ylimber**2)), sngl(flux)
	write (37,*)           sngl(heure*3600.d0), sngl(flux)
	if (pres.ne.0.d0.and.dens.ne.0.d0) then
		write (38,*)   sngl(heure*3600.d0), sngl( dlog10(pres) )
		write (39,*)   sngl(heure*3600.d0), sngl( dlog10(dens) )
	else
		write (38,*)   sngl(heure*3600.d0), sngl(-100000.d0)
		write (39,*)   sngl(heure*3600.d0), sngl(-100000.d0)
	endif
	if (sigma.ne.0.d0) then
		write (40,*)   sngl(heure*3600.d0), sngl( dlog10(sigma) )
	else
		write (40,*)   sngl(heure*3600.d0), sngl(-100000.d0)
	endif
c
c Ramene au niveau de pression de reference pres0
c
	pres0= 2.d0		! pression de reference en mubar
	if (pres.gt.1.d-30) then
	        d_z= H*dlog(pres/pres0)
	else
		d_z= 0.d0	! correction qui n'a pas vraiment de sens
	endif
	rlimber= dsqrt( xlimber**2 + ylimber**2 )
        rlimber0= rlimber + d_z
        correc= rlimber0/rlimber
        xlimber0= (xlimber*correc)
        ylimber0= (ylimber*correc)
c        write (41,*) sngl(ksi), sngl(eta), sngl(xlimber), sngl(ylimber),
c     *  sngl(pres), sngl(xlimber0), sngl(ylimber0), sngl(H), sngl(d_z)

	write (42,*) sngl(heure*3600.d0), sngl(rmin - r_limbe)

	write(*,*)  ' '
	write(*,*)  'ksi, eta, vitesse perp. au limbe',
     *  sngl(ksi), sngl(eta), sngl(vrel)
	write(*,*)  'correction relativiste radiale (pas de J2)', drelat
	write(*,*)  'xlimber, ylimber (corr. relat. et J2)',
     *  sngl(xlimber), sngl(ylimber)
	write(*,*)  'Correction radiale due refraction:',
     *  sngl(rmin-dsqrt(ksi**2+eta**2))
	write(*,*)  'Correction radiale due relativite:',
     *  sngl(dsqrt(xlimber**2+ylimber**2)-rmin)
	write(*,*) 'heure (sec), phi',sngl(heure*3600.d0),sngl(phi*(180./pi))
	write(*,*) 'heure (sec), phip',sngl(heure*3600.d0),sngl(phip*(180./pi))
	xdumb= dcos(alplan)*dcos(phi)
	ydumb= dsin(alplan)*dcos(phi)
	write(*,*) 'projection longitude, lat. planeto',sngl(xdumb),sngl(ydumb)
	xdumb= dcos(alplan)*dcos(phip)
	ydumb= dsin(alplan)*dcos(phip)
	write(*,*) 'projection longitude, lat. geogra.',sngl(xdumb),sngl(ydumb)
	write(*,*)   sngl(heure*3600.d0), sngl(phi*(180./pi)),
     *          sngl(alplan*(180./pi) -90.)
	write(*,*)  'Heure, flux:',                   sngl(heure), sngl(flux)
	write(*,*)  'Heure, pression (mubar):',       sngl(heure), sngl(pres)
	write(*,*)  'Heure, densite (molecules/m3):', sngl(heure), sngl(dens)
	write(*,*)  'heure, densite mol. surf.',      sngl(heure), sngl(sigma)
	write(*,*)  'heure, alt. au-dessus du limbe de mi-occultation local:',
     *          sngl(heure*3600.d0), sngl(rmin - r_limbe)

	return
	end

**************************** fin de latitude **********************************

	subroutine baum_code (Req,H,z_2,z, phi,r)
c
c	entrees:Req= rayon de la surface solide
c		H= echelle de hauteur
c	        z_2= rayon de mi-occultation DANS L'OMBRE
c		z=   distance de L'OBSERVATEUR au centre de l'ombre
c
c	sorties:phi= flux
c	        r=   rayon de plus courte approche dans l'atmosphere
c
c
	implicit real*8 (a-h,o-z)
	tol= 1.d-4			! tolerance
	itermax= 20			! nombre d'iterations maximal

	y= -(z-z_2)/H
c
c	point de depart
c
	if (y.lt.0.) u0= dexp(1.d0 + y)
	if (y.ge.0.) u0=      1.d0 + y
c
c	iterations
c
	err= 2.d0*tol
	iter= 0
	do while (err.ge.tol)
		if (u0.eq.0.d0) then
			go to 1
		endif
		u= u0*( (y+2.d0-dlog(u0))/(1.d0+u0) )
		err= dabs(u-u0)
		u0= u
		iter= iter + 1
		if (iter.ge.itermax) then
		   write(*,*)  'Probleme: iteration maximale atteinte!', itermax
		endif
	enddo

 1	phi= 1.d+0/(1.d0+u0)
	r=   z + H*( (1.d0/phi)-1.d0 )
	if (r.le.Req) phi= 0.		! le rayon rencontre la surface
	phi= dabs(r/z)*phi		! focalisation spherique

	return
	end
************************************ fin de baum_code ************************

	subroutine image (P,B,GM,c,R_g,J_2,diskm,r_eq,e,rmin,ksi,eta,vksi,veta,
     *                    vrel,xlimbe,ylimbe,xlimber,ylimber)
c
c attention: e= aplatissement APPARENT
c
c Calcule la position de l'image principale pour un observateur ksi, eta
c L'image se situe a la distance rmin du centre de l'ombre
c Attention: ce programme n'est valide que loin de la caustique centrale
c
	implicit real*8 (a-h,o-z)
	real*8 J_2, ksi
	pi= dacos(-1.d0)
	cosP= dcos(P)
	sinP= dsin(P)

	aa= r_eq
	bb= aa*(1.d0 - e)	! petit-axe APPARENT (bb: different de B!)

	X= ksi*cosP - eta*sinP	! coordonnees du point dans une repere aligne
	Y= ksi*sinP + eta*cosP	! avec l'ellipsoide (apparent)
	R= dsqrt(X**2 + Y**2)

	VX= vksi*cosP - veta*sinP ! vitesse du point dans une repere aligne
	VY= vksi*sinP + veta*cosP ! avec l'ellipsoide (apparent)

	alpha= datan(Y/X)	! angle du point X,Y
	if (X.le.0.d0)     alpha= alpha +      pi
	if (alpha.le.0.d0) alpha= alpha + 2.d0*pi

	denom = dsqrt( X**2 + ((1.d0-e)*Y)**2 )
	cosphi=            X/denom
	sinphi= ((1.d0-e)*Y)/denom
	coeff = aa*e*(1.d0-e)/denom

	phi= datan(sinphi/cosphi)
	if (cosphi.le.0.d0) phi= phi +      pi
	if (phi.le.0.d0)    phi= phi + 2.d0*pi

	theta0= phi
	tol= 1.d-4
	do i= 1, 20
		theta = phi + dasin(coeff*dcos(theta0)*dsin(theta0))
		if (dabs(theta-theta0).le.tol) go to 1
		theta0= theta
	enddo
 1	continue

	beta= datan((aa/bb)*dtan(theta)) ! angle de la normale au limbe avec OX
	if (dcos(theta).le.0.d0) beta= beta +      pi
	if (beta.le.0.d0)        beta= beta + 2.d0*pi

	gamma= beta - dasin( (R/rmin)*dsin(beta-alpha) )
	XX= rmin*dcos(gamma)	! coordonneees du point du limbe dans le
	YY= rmin*dsin(gamma)	! repere OXY

	xlimbe=  XX*cosP + YY*sinP ! point du limbe dans le repere O-ksi-eta
	ylimbe= -XX*sinP + YY*cosP ! du uniquement a la refraction
c
c deviation relativiste: (NB. la conversion GM/1.d9 pour km^3/sec^2)
c
	fact= (4.d0*(GM/1.d9)*diskm)/(c*rmin)**2
	DX= 1.d0 - J_2*(R_g**2)*((dcos(B))**2)*(3.d0*YY**2-XX**2)/rmin**4
	DY= 1.d0 + J_2*(R_g**2)*((dcos(B))**2)*(3.d0*XX**2-YY**2)/rmin**4
	DX= fact*DX*XX
	DY= fact*DY*YY
c
	XX= XX + DX
	YY= YY + DY
c
	xlimber=  XX*cosP + YY*sinP ! point du limbe dans le repere O-ksi-eta
	ylimber= -XX*sinP + YY*cosP ! du a la refraction et la relativite

	vrel= VX*dcos(beta) + VY*dsin(beta) ! vitesse de l'observateur (dans
					    ! l'ombre) projetee sur la normale
					    ! locale au limbe
 	return
	end
*********************************** fin de image ****************************

	subroutine new_parallaxe (iobs,heure,   alpha,delta,ksip,etap)
c
c NB. alpha, delta: coordonnees de l'etoile (radian)
c	
	implicit real*8 (a-h,o-z)
	real*8 lat, lon, nx, ny, nz, ksip
	integer*4 aa
	pi= dacos(-1.d0)
	radh= pi/12.d0
	radd= pi/180.d0

	if (iobs.ge.51.and.iobs.le.0107) then
		iah=   21		! coordonnees J2000 de l'etoile
		iam=   38
		asec=  13.96284d+0
		idd=  -14
		imd=  -54
		dsec= -35.9070d+0
		jj= 08
		mm= 09
		aa= 2001
	endif

	if (iobs.ge.0110.and.iobs.le.0111) then
		iah=   22		! coordonnees J2000 de l'etoile
		iam=   15
		asec=  54.55430+0
		idd=  -11
		imd=  -36
		dsec= -56.3584+0
		jj= 01
		mm= 08
		aa= 2003
	endif

	if (iobs.ge.0120.and.iobs.le.0130) then
		iah=   17		! coordonnees J2000 de l'etoile
		iam=   00
		asec=  18.03104d0
		idd=  -12
		imd=  -41
		dsec= -42.028d0
		jj= 20
		mm= 07
		aa= 2002
	endif

	if (iobs.ge.0150.and.iobs.le.0153) then
		iah=   16		! coordonnees J2000 de l'etoile
		iam=   58
		asec=  49.43132d0
		idd=  -12
		imd=  -51
		dsec= -31.869d0
		jj= 21
		mm= 08
		aa= 2002
	endif

	if (iobs.ge.0160.and.iobs.le.0173) then
		iah=   05		! coordonnees J2000 de l'etoile
		iam=   41
		asec=  33.8871d0
		idd=   22
		imd=   03
		dsec=  31.214 
		jj= 15
		mm= 12
		aa= 2002
	endif

	if (iobs.ge.0200.and.iobs.le.0207) then
		iah=   08		! coordonnees J2000 de l'etoile
		iam=   42
		asec=  42.4732d0 
		idd=   19
		imd=   05
		dsec=  58.880d0
		jj= 01
		mm= 04
		aa= 2003
	endif

	if (iobs.ge.0220.and.iobs.le.0242) then
		iah=   06		! coordonnees J2000 de l'etoile #1
		iam=   55
		asec=  20.9672
		idd=   22
		imd=   06
		dsec=  07.812 
		jj= 14
		mm= 11
		aa= 2003
	endif

	if (iobs.ge.0243.and.iobs.le.0245) then
		iah=   06		! coordonnees J2000 de l'etoile #2
		iam=   55		! (Manek)
		asec=  17.7690
		idd=   22
		imd=   06
		dsec=  01.226 
		jj= 14
		mm= 11
		aa= 2003
	endif

        if (iobs.ge.0260.and.iobs.le.0274) then
                iah=   17               ! coordonnees J2000 (Stone + JLX 4/9/04)
                iam=   28               
                asec=  55.0130 
                idd=  -15
                imd=  -00
                dsec= -54.719 
                jj= 11
                mm= 07
                aa= 2005
        endif

        if (iobs.ge.0300.and.iobs.le.0301) then
                iah=   17               ! coordonnees J2000 Behrend 06 mars 2006
                iam=   46               
                asec=  06.880 
                idd=  -15
                imd=  -46
                dsec= -10.11 
                jj= 10
                mm= 04
                aa= 2006
        endif

        if (iobs.ge.0310.and.iobs.le.0315) then
                iah=   17               ! coordonnees J2000 Ron Stone Aug. 2004
                iam=   41               
                asec=  12.0940
                idd=  -15
                imd=  -41
                dsec= -34.447 
                jj= 12
                mm= 06
                aa= 2006
        endif

        if (iobs.ge.0320.and.iobs.le.0320) then
                iah=   17               ! coordonnees J2000 R. Behrend 15 jul 06
                iam=   36               
                asec=  07.8427
                idd=  -15
                imd=  -49
                dsec= -32.792
                jj= 06
                mm= 08
                aa= 2006
        endif

        if (iobs.ge.0330.and.iobs.le.0349) then
                iah=   17               ! coordonnees J2000 ---> ICRF R. Behrend 04 mars 07
                iam=   55               
                asec=  05.6971
                idd=  -16
                imd=  -28
                dsec= -34.360
                jj= 18
                mm= 03
                aa= 2007
        endif

        if (iobs.ge.0360.and.iobs.le.0367) then
                iah=   17               ! coordonnees J2000 ---> ICRF M. Assafin, mail 30 nov 07 et cahier 28 jan 07
                iam=   53              
                asec=  32.1006
                idd=  -16
                imd=  -22
                dsec= -47.367
                jj= 12
                mm= 05
                aa= 2007
        endif

        if (iobs.ge.0380.and.iobs.le.0390) then			! etoile 14 juin 2007
                iah=   17               ! coordonnees J2000 ---> ICRF M. Assafin, mail 30 nov 07 et cahier 28 jan 07
                iam=   50               
                asec=  20.7412
                idd=  -16
                imd=  -22
                dsec= -42.236
                jj= 14
                mm= 06
                aa= 2007
        endif

        if (iobs.ge.0391.and.iobs.le.0392) then			! etoile 09 juin 2007
                iah=   17								! coordonnees J2000 ---> ICRF R. Behrend, catalogue 16 avril 2007
                iam=   50               
                asec=  50.6535
                idd=  -16
                imd=  -22
                dsec= -29.297
                jj= 09
                mm= 06
                aa= 2007
        endif

        if (iobs.ge.0395.and.iobs.le.0397) then			! etoile 31 juillet 2007
                iah=   17								! coordonnees J2000 ---> ICRF R. Martins, mail 30 juillet 2007
                iam=   45               
                asec=  41.98943
                idd=  -16
                imd=  -29
                dsec= -31.639
                jj= 31
                mm= 07
                aa= 2007
        endif

        if (iobs.ge.0400.and.iobs.le.0410) then			! etoile 21 mai 2008
                iah=   21								! R. Martins, mail 27 avril 2008
                iam=   46               
                asec=  11.049528d0
                idd=  -13
                imd=  -46
                dsec= -45.9484716d0
                jj= 21
                mm= 05
                aa= 2008
        endif

        if (iobs.ge.0420.and.iobs.le.0434) then			! etoile 22 juin 2008, Pluton
                iah=   17								! R. Martins, mail 02 juin 2008
                iam=   58               
                asec=  33.01376909d0
                idd=  -17
                imd=  -02
                dsec= -38.34870554d0
                jj= 22
                mm= 06
                aa= 2008
        endif
		
        if (iobs.ge.0440.and.iobs.le.0443) then			! etoile 22 juin 2008, Charon
                iah=   17								! R. Martins, mail 02 juin 2008
                iam=   58               
                asec=  33.01376909d0
                idd=  -17
                imd=  -02
                dsec= -38.34870554d0
                jj= 22
                mm= 06
                aa= 2008
        endif

        if (iobs.ge.0450.and.iobs.le.0450) then			! etoile 24 juin 2008, Pluton
                iah=   17								! R. Martins, mail 09 juin 2008
                iam=   58               
                asec=  22.392984d0
                idd=  -17
                imd=  -02
                dsec= -49.349328d0
                jj= 24
                mm= 06
                aa= 2008
        endif


        if (iobs.ge.0460.and.iobs.le.0464) then			! etoile 30 juin 2009, 2002 MS4
                iah=   18								! M. Assafin, mail 25 juin 2009
                iam=   03               
                asec=  37.93133d0
                idd=  -08
                imd=  -28
                dsec= -17.6971d0
                jj= 30
                mm= 06
                aa= 2009
        endif


        if (iobs.ge.0470.and.iobs.le.0471) then			! etoile 25 aout 2008, Pluton
                iah=   17								! R. Martins, mail 23 juillet 2008
                iam=   53               
                asec=  27.103164d0
                idd=  -17
                imd=  -15
                dsec= -27.54630d0
                jj= 25
                mm= 08
                aa= 2008
        endif

        if (iobs.ge.0480.and.iobs.le.0491) then			! etoile 19 fevrier 2010, Varuna
                iah=   07								! R. Martins, mail 04 fevrier 2010
                iam=   29               
                asec=  22.4714d0
                idd=   26
                imd=   07
                dsec=  23.173d0
                jj= 19
                mm= 02
                aa= 2010
        endif

        if (iobs.ge.0500.and.iobs.le.0502) then			! etoile 14 fevrier 2010, Pluton
                iah=   18								! Rio 10 nov. 2010, voir g3_occ_data_pluto_2010_table
                iam=   19               
                asec=  14.3851d0
                idd=  -18
                imd=  -16
                dsec= -42.313d0
                jj= 14
                mm= 02
                aa= 2010
        endif

        if (iobs.ge.0510.and.iobs.le.0512) then			! etoile 19 mai 2010, Pluton
                iah=   18								! Rio 10 nov. 2010, voir g3_occ_data_pluto_2010_table
                iam=   20               
                asec=  16.7941d0
                idd=  -18
                imd=  -11
                dsec= -48.070d0
                jj= 19
                mm= 05
                aa= 2010
        endif

        if (iobs.ge.0515.and.iobs.le.0520) then			! etoile 04 juin 2010, Pluton
                iah=   18								! Rio 10 nov. 2010, voir g3_occ_data_pluto_2010_table
                iam=   18               
                asec=  47.9349d0
                idd=  -18
                imd=  -12
                dsec= -51.794d0
                jj= 04
                mm= 06
                aa= 2010
        endif

        if (iobs.ge.0530.and.iobs.le.0530) then			! etoile 04 juillet 2010, Pluton
                iah=   18								! Rio 02 juillet 2010
                iam=   15               
                asec=  42.1027d0
                idd=  -18
                imd=  -16
                dsec= -41.124d0
                jj= 04
                mm= 07
                aa= 2010
        endif

        if (iobs.ge.0540.and.iobs.le.0544) then			! etoile Ceres 17 aout 2010, UCAC2
                iah=   17
                iam=   18               
                asec=  29.008d0
                idd=  -27
                imd=  -26
                dsec= -38.89d0
                jj= 17
                mm= 08
                aa= 2010
        endif

        if (iobs.ge.0545.and.iobs.le.0546) then			! etoile 2003AZ84 08/01/2011, WFI 01/ago/11
                iah=   07
                iam=   43               
                asec=  41.8220d0
                idd=  +11
                imd=  +30
                dsec= +23.569d0
                jj= 08
                mm= 01
                aa= 2011
        endif


        if (iobs.ge.0550.and.iobs.le.0554) then			! etoile Bill Owen, mail 4 nov. 2010, Eris **6 NOVEMBRE 2010**
				iah=   01
                iam=   39               
                asec=  09.9392d0
                idd=  -04
                imd=  -21
                dsec= -12.140
                jj= 06
                mm= 11
                aa= 2010
        endif

        if (iobs.ge.0560.and.iobs.le.0560) then			! etoile Bill Owen, mail 4 nov. 2010, Eris **5 NOVEMBRE 2010**
				iah=   01
                iam=   39               
                asec=  09.9392d0
                idd=  -04
                imd=  -21
                dsec= -12.140
                jj= 05
                mm= 11
                aa= 2010
        endif

c        if (iobs.ge.0570.and.iobs.le.0575) then			! etoile JL Ortiz demande DDT ESO 18 avril 2011, Makemake 23 avril 2011
c		iah=   12
c                iam=   36               
c                asec=  11.4140d0
c                idd=  +28
c                imd=   11
c                dsec=  10.606d0
c                jj= 23
c                mm= 04
c                aa= 2011
c        endif

        if (iobs.ge.0570.and.iobs.le.0575) then			! etoile  WFI 01/08/2011 position, Makemake 23 avril 2011
		iah=   12
                iam=   36               
                asec=  11.3938d0
                idd=  +28
                imd=   11
                dsec=  10.493d0
                jj= 23
                mm= 04
                aa= 2011
        endif

        if (iobs.ge.0584.and.iobs.le.0599) then			! star email 20/ago/2011 Julio, Quaoar 04 mai 2011
		iah=   17					!, obs IAG 18/ago (er RA 7mas e DEC 14mas)
                iam=   28               
                asec=  50.8009d0
                idd=  -15
                imd=  -27
                dsec= -42.770d0
                jj= 04
                mm= 05
                aa= 2011
        endif
 
        if (iobs.ge.0600.and.iobs.le.0603) then			! etoile 04 juin 2011, Pluton
                iah=   18								! Rio 10 nov. 2009, voir g3_occ_data_pluto_2011_table
                iam=   27               
                asec=  53.8249d0
                idd=  -18
                imd=  -45
                dsec= -30.725
                jj= 04
                mm= 06
                aa= 2011
        endif

        if (iobs.ge.0610.and.iobs.le.0613) then			! etoile 04 juin 2011, Charon
                iah=   18								! Rio 10 nov. 2009, voir g3_occ_data_pluto_2011_table
                iam=   27               
                asec=  53.8249d0
                idd=  -18
                imd=  -45
                dsec= -30.725
                jj= 04
                mm= 06
                aa= 2011
        endif

        if (iobs.ge.0620.and.iobs.le.0629) then			! etoile 23 juin 2011, Pluton
                iah=   18								! Rio 10 nov. 2009, voir g3_occ_data_pluto_2011_table
                iam=   25               
                asec=  55.4750
                idd=  -18
                imd=  -48
                dsec= -07.015
                jj= 23
                mm= 06
                aa= 2011
        endif

        if (iobs.ge.0630.and.iobs.le.0639) then			! etoile 23 juin 2011, Charon
                iah=   18								! Rio 10 nov. 2009, voir g3_occ_data_pluto_2011_table
                iam=   25               
                asec=  55.4750
                idd=  -18
                imd=  -48
                dsec= -07.015
                jj= 23
                mm= 06
                aa= 2011
        endif

        if (iobs.ge.0640.and.iobs.le.0649) then			! etoile 27 juin 2011, Pluton
                iah=   18								! Rio 10 nov. 2009, voir g3_occ_data_pluto_2011_table
                iam=   25               
                asec=  29.0100
                idd=  -18
                imd=  -48
                dsec= -47.570
                jj= 27
                mm= 06
                aa= 2011
        endif

        if (iobs.ge.0650.and.iobs.le.0659) then			! etoile 27 juin 2011, Hydra
                iah=   18								! Rio 10 nov. 2009, voir g3_occ_data_pluto_2011_table
                iam=   25               
                asec=  29.0100
                idd=  -18
                imd=  -48
                dsec= -47.570
                jj= 27
                mm= 06
                aa= 2011
        endif


        if (iobs.ge.0660.and.iobs.le.0679) then			! etoile 27 juin 2011, Hydra
                iah=   07								! Rio 10 nov. 2009, voir g3_occ_data_pluto_2011_table
                iam=   45               
                asec=  54.76965
                idd=   11
                imd=   12
                dsec=  43.0933
                jj= 03
                mm= 02
                aa= 2012
        endif
	
	
        if (iobs.ge.0680.and.iobs.le.0683) then			! etoile 27 juin 2011, Hydra
                iah=   17								! Rio 10 nov. 2009, voir g3_occ_data_pluto_2011_table
                iam=   34               
                asec=  21.8453
                idd=  -15
                imd=  -42
                dsec= -10.586
                jj= 17
                mm= 02
                aa= 2012
        endif

        if (iobs.ge.0685.and.iobs.le.0687) then			! etoile 03 December 2015, Gaia DR1, 2002 VE95
                iah=   5								
                iam=   37               
                asec=  38.52788
                idd=  +8
                imd=  +5
                dsec= +9.599577
                jj= 3
                mm= 12
                aa= 2015
        endif

        if (iobs.ge.0688.and.iobs.le.0688) then			! etoile 09 juillet 2013, Gaia DR1, Quaoar (tabela online)
                iah=   17								
                iam=   34               
                asec=  40.46381
                idd=  -15
                imd=  -23
                dsec= -37.568918
                jj= 9
                mm= 7
                aa= 2013
        endif

       if (iobs.ge.0689.and.iobs.le.0690) then			! January 13, 2013, Gaia DR1 star, Sedna (pdf na pasta)  
                iah=   3								
                iam=   32               
                asec=  44.70578
                idd=  +6
                imd=  +47
                dsec= +18.998235
                jj= 13
                mm= 1
                aa= 2013
        endif

       if (iobs.ge.0691.and.iobs.le.0692) then			! August 29, 2013, Gaia DR1 star, Eris (planilha google)    
                iah=   1								
                iam=   42               
                asec=  58.50527
                idd=  -3
                imd=  -21
                dsec= -17.342460
                jj= 29
                mm= 8
                aa= 2013
        endif

       if (iobs.ge.0693.and.iobs.le.0698) then		       ! November 07, 2014, Gaia DR1 star, 2003 VS2 (planilha google)    
                iah=   4								
                iam=   48               
                asec=  32.13899
                idd=   33
                imd=   58
                dsec=  36.407472
                jj= 7
                mm= 11
                aa= 2014
        endif


c GAIA: 04 48 32.1390	+33 58 36.407
c 4UC 620-017761  (RA 4 48 32.1457, Dec +33 58 36.365
       if (iobs.ge.0700.and.iobs.le.0703) then    
                iah=   4								
                iam=   48               
                asec=  32.13899
                idd=   33
                imd=   58
                dsec=  36.407472
                jj= 7
                mm= 11
                aa= 2014
        endif



	alpha= dfloat(iah) + dfloat(iam)/60.d0 + asec/3600.d0
	alpha= alpha*radh
	delta= dfloat(idd) + dfloat(imd)/60.d0 + dsec/3600.d0
	delta= delta*radd

	if (iobs.eq.51) isite= 31		! Salinas
	if (iobs.eq.52) isite= 32		! Cuenca
	if (iobs.eq.53) isite= 33		! Maracaibo
	if (iobs.eq.54) isite= 34		! Bobares
	if (iobs.eq.55) isite= 35		! Caracas
	if (iobs.eq.56) isite= 36		! Arikok, Aruba
	if (iobs.eq.57) isite= 37		! Fort de France
	if (iobs.eq.58) isite= 38		! Tobago
	if (iobs.eq.59) isite= 39		! Trinidad
	if (iobs.eq.60) isite= 40		! Bridgetown, Barbados
	if (iobs.eq.61) isite= 41		! Acores
	if (iobs.eq.62) isite= 42		! Oeiras
	if (iobs.eq.63) isite= 43		! Portimao
	if (iobs.eq.64) isite= 44		! COAA Algarve (Portimao)
	if (iobs.eq.65) isite= 45		! Linhaceira
	if (iobs.eq.66) isite= 46		! Alvito
	if (iobs.eq.67) isite= 47		! Granada
	if (iobs.eq.68) isite= 48		! Alcublas
	if (iobs.eq.69) isite= 49		! Bordeaux
	if (iobs.eq.70) isite=  4 		! Pic
	if (iobs.eq.71) isite= 50		! St Maurice (Poitou)
	if (iobs.eq.72) isite= 51		! Pezenas
	if (iobs.eq.73) isite= 52		! Mauguio
	if (iobs.eq.74) isite= 53		! Nimes
	if (iobs.eq.75) isite= 54		! Orfeuilles
	if (iobs.eq.76) isite=  5 		! OHP
	if (iobs.eq.77) isite= 55		! Binfield
	if (iobs.eq.78) isite= 56		! Worth Hill, Bournemouth
	if (iobs.eq.79) isite= 57		! Teversham
	if (iobs.eq.80) isite= 58		! Plateau d'Albion
	if (iobs.eq.81) isite= 59		! Bd H. Fabre, Marseille
	if (iobs.eq.82) isite= 60 		! Guitalens
	if (iobs.eq.83) isite= 61 		! Sabadell (Ardanuy)
	if (iobs.eq.84) isite= 62 		! Salon
	if (iobs.eq.85) isite= 63 		! Wela, Aruba
	if (iobs.eq.86) isite= 64 		! St Maurice Cazevieille
	if (iobs.eq.87) isite= 65		! Marinha Grande
	if (iobs.eq.88) isite= 66		! Almeirim
	if (iobs.eq.89) isite= 67		! Carcavelos
	if (iobs.eq.90) isite= 68		! Setubal
	if (iobs.eq.91) isite= 69		! Alcacer do Sal
	if (iobs.eq.92) isite= 70 		! Barcelona 1
	if (iobs.eq.93) isite= 71 		! St Savinien
	if (iobs.eq.94) isite= 72 		! St Martin de Crau
	if (iobs.eq.95) isite= 73 		! Barcelona 2
	if (iobs.eq.96) isite= 74 		! St Esteve
	if (iobs.eq.97) isite= 75 		! Alella
	if (iobs.eq.98) isite= 76 		! Castellon
	if (iobs.eq.99) isite= 77 		! Hortoneda
	if (iobs.eq.0100)isite= 78 		! Barcelona 3
	if (iobs.eq.0101)isite= 61 		! Sabadell (Casas)
	if (iobs.eq.0102)isite= 79 		! Esplugues de Llobregat
	if (iobs.eq.0103)isite= 80 		! Zaragoza
	if (iobs.eq.0104)isite= 81 		! Calern
	if (iobs.eq.0105)isite= 82 		! Dax
	if (iobs.eq.0106)isite= 182 	! Puimichel
	if (iobs.eq.0107)isite=	183		! Chatellerault

	if (iobs.eq.0110)isite=	175		! George Observatory
	if (iobs.eq.0111)isite=	184		! Monterrey (hacienda Zuazua)

	if (iobs.eq.0120)isite= 83 		! Arica
	if (iobs.eq.0121)isite= 84 		! Jerusalem
	if (iobs.eq.0122)isite= 85 		! Cumbaya
	if (iobs.eq.0123)isite= 86 		! Lima
	if (iobs.eq.0124)isite= 87 		! Mamina
	if (iobs.eq.0125)isite= 26 		! Paranal
	if (iobs.eq.0126)isite= 88 		! Los Armazones
	if (iobs.eq.0127)isite= 89 		! Itajuba
	if (iobs.eq.0128)isite=  3 		! La Silla
	if (iobs.eq.0129)isite= 90 		! El Leoncito
	if (iobs.eq.0130)isite= 91 		! Merida

	if (iobs.eq.0150)isite=   1 	! CFHT 
	if (iobs.eq.0151)isite=  92 	! Lowell
	if (iobs.eq.0152)isite=  93		! Palomar
	if (iobs.eq.0153)isite=  94		! Lick

	if (iobs.eq.0160)isite=  95 	! Tomar
	if (iobs.eq.0161)isite=  96 	! Sabadell (Casas) 
	if (iobs.eq.0162)isite=   4 	! Pic
	if (iobs.eq.0163)isite=  97 	! Czarna Bialostocka 
	if (iobs.eq.0164)isite=  98 	! Vitebsk
	if (iobs.eq.0165)isite=   8 	! Meudon 

	if (iobs.eq.0166)isite=  99		! Kooriyama
	if (iobs.eq.0167)isite= 100		! Ooe
	if (iobs.eq.0168)isite= 101		! Kashiwa
	if (iobs.eq.0169)isite= 102		! Hitachi
	if (iobs.eq.0170)isite= 103		! Chichibu
	if (iobs.eq.0171)isite= 104		! Musashino
	if (iobs.eq.0172)isite= 105		! Mitaka
	if (iobs.eq.0173)isite= 106		! Abrera 

	if (iobs.eq.0200)isite=   4		! Pic 
	if (iobs.eq.0201)isite=  93		! Palomar
	if (iobs.eq.0202)isite=   3		! La Silla
	if (iobs.eq.0203)isite= 107		! Roque de los Muchachos
	if (iobs.eq.0204)isite= 108		! New Jersey 
	if (iobs.eq.0205)isite= 109		! Max Valier Observatory 
	if (iobs.eq.0206)isite= 110		! Nyrola Observatory 
	if (iobs.eq.0207)isite= 111		! Livermore

	if (iobs.eq.0220)isite= 112		! Windhoek
	if (iobs.eq.0221)isite= 113		! Tivoli
	if (iobs.eq.0222)isite= 114		! HESS
	if (iobs.eq.0223)isite= 115		! Hakos
	if (iobs.eq.0224)isite= 116		! Kleinbegin
	if (iobs.eq.0225)isite= 117		! Sandfontein
	if (iobs.eq.0226)isite= 118		! Springbok
	if (iobs.eq.0227)isite= 119		! Nuwerus
	if (iobs.eq.0228)isite= 120		! Gifberg
	if (iobs.eq.0229)isite= 121		! Cederberg
	if (iobs.eq.0230)isite= 122		! SAAO Sutherland 
	if (iobs.eq.0231)isite= 123		! SAAO Cape Town
	if (iobs.eq.0232)isite= 124		! Boyden

	if (iobs.eq.0240)isite= 140		! Maido
	if (iobs.eq.0241)isite= 141		! Les Makes
	if (iobs.eq.0242)isite= 142		! Fournaise

	if (iobs.eq.0243)isite= 149		! Pico Veleta
	if (iobs.eq.0244)isite= 150		! WIRO
	if (iobs.eq.0245)isite=  91		! Merida

	if (iobs.eq.0260)isite=   3 		! La Silla
	if (iobs.eq.0261)isite=  26 		! Paranal 
	if (iobs.eq.0262)isite= 143 		! Montevideo 
	if (iobs.eq.0263)isite= 144 		! Bosque Alegre
	if (iobs.eq.0264)isite=  90 		! El Leoncito
	if (iobs.eq.0265)isite=  89 		! Itajuba
	if (iobs.eq.0266)isite= 145 		! Asuncion
	if (iobs.eq.0267)isite= 146 		! San Pedro de Atacama
	if (iobs.eq.0268)isite= 147 		! Tarija
	if (iobs.eq.0269)isite= 148 		! Marangani
	if (iobs.eq.0270)isite= 151 		! Wykrota, Bresil
	if (iobs.eq.0271)isite= 152 		! CEAMIG-REA, Bresil
	if (iobs.eq.0272)isite= 153 		! Patacamaya, Bolivie	
	if (iobs.eq.0273)isite= 154			! Tilomonte, Chili
	if (iobs.eq.0274)isite= 155			! SOAR
	
	if (iobs.eq.0300)isite=  26			! Paranal
	if (iobs.eq.0301)isite=   3			! La Silla

	if (iobs.eq.0310)isite= 157 		! Mount John
	if (iobs.eq.0311)isite=  14 		! Hobart
	if (iobs.eq.0312)isite= 158  		! AAT
	if (iobs.eq.0313)isite= 141			! Les Makes
	if (iobs.eq.0314)isite= 159			! Gault, Blue Mountains
	if (iobs.eq.0315)isite= 160			! Blair, Stockport Observatory
	
	if (iobs.eq.0320)isite=	 26			! Paranal
	
	if (iobs.eq.0330)isite=	 94 		! Lick
	if (iobs.eq.0331)isite=	163 		! Kitt Peak (Bok 90")
	if (iobs.eq.0332)isite=	164			! Mt Bigelow (Kuiper tel. 61")
	if (iobs.eq.0333)isite=	165			! Tenagra
	if (iobs.eq.0334)isite=	166			! Guanajuato
	if (iobs.eq.0335)isite=	167			! MMT
	if (iobs.eq.0336)isite=	168			! Pinto Valley
	if (iobs.eq.0337)isite=	169			! Hereford, AZ
	if (iobs.eq.0338)isite=	170			! Cloudbait, CO
	if (iobs.eq.0339)isite=	171			! Palmer Divide, CO
	if (iobs.eq.0340)isite=	172 		! Calvin-Rehoboth, NM
	if (iobs.eq.0341)isite=	173 		! Appalachian, NC
	if (iobs.eq.0342)isite=	174 		! Moore, WA
	if (iobs.eq.0343)isite=	175 		! George, TX
	if (iobs.eq.0344)isite=	 92 		! Lowell, AZ				
	if (iobs.eq.0345)isite=	 93 		! Palomar
	if (iobs.eq.0346)isite=	176 		! San Pedro Martir
	if (iobs.eq.0347)isite=	177 		! McDonald, TX
	if (iobs.eq.0348)isite=	178 		! Univ. Oklahoma
	if (iobs.eq.0349)isite=	179 		! KASI, Mt Lemmon
	
	if (iobs.eq.0360)isite=	 26 		! Paranal
	if (iobs.eq.0361)isite=	  3			! La Silla
	if (iobs.eq.0362)isite=	 88			! Los Armazones
	if (iobs.eq.0363)isite=	146			! San Pedro de Atacama
	if (iobs.eq.0364)isite=	 89			! Itajuba			
	if (iobs.eq.0365)isite=	112			! Windoek
	if (iobs.eq.0366)isite=	115			! Hakos
	if (iobs.eq.0367)isite=	117			! Sandfontein
			
	if (iobs.eq.0380)isite=	 26 		! Paranal
	if (iobs.eq.0381)isite=	  3  		! La Silla
	if (iobs.eq.0382)isite=	155 		! SOAR
	if (iobs.eq.0383)isite=	146 		! San Pedro de Atacama
	if (iobs.eq.0384)isite=	 89 		! Itajuba
	if (iobs.eq.0385)isite=	122  		! SALT
	if (iobs.eq.0386)isite=	115  		! Hakos
	if (iobs.eq.0387)isite=	112  		! Windhoek
	if (iobs.eq.0388)isite=	141  		! Les Makes
	if (iobs.eq.0389)isite=	180  		! Mairinque

	if (iobs.eq.0390)isite=	 89 		! Itajuba, ephemeride BARYCENTRIQUE
	
	if (iobs.eq.0391)isite=	155 		! SOAR
	if (iobs.eq.0392)isite=	  3 		! La Silla
	
	if (iobs.eq.0395)isite=	157 		! Mt John
	if (iobs.eq.0396)isite=	 14			! Hobart (Canopus)
	if (iobs.eq.0397)isite=	181			! Canberra

	if (iobs.eq.0400)isite=	115			! Hakos				Triton 21 mai 2008
	if (iobs.eq.0401)isite= 112			! Windhoek
	if (iobs.eq.0402)isite=	113			! Tivoli
	if (iobs.eq.0403)isite=	185			! Sossusvlei
	if (iobs.eq.0404)isite=	118			! Sprinbok
	if (iobs.eq.0405)isite=	121			! Cederberg
	if (iobs.eq.0406)isite=	122			! SAAO
	if (iobs.eq.0407)isite=	123			! Cape Town
	if (iobs.eq.0408)isite=	141			! Les Makes
	if (iobs.eq.0409)isite=	140			! Maido
	if (iobs.eq.0410)isite=	142			! Fournaise

	if (iobs.eq.0420)isite=	 14			! Hobart			Pluton 22 juin 2008
	if (iobs.eq.0421)isite= 159			! Blue Mountains
	if (iobs.eq.0422)isite= 122			! SAAO
	if (iobs.eq.0423)isite= 181			! Canberra
	if (iobs.eq.0424)isite= 158			! Siding Spring
	if (iobs.eq.0425)isite= 160			! Stockport
	if (iobs.eq.0426)isite= 186			! Perth, Lowell Telescope
	if (iobs.eq.0427)isite= 187			! Perth, Roger Groom,          22 juin 2008
	if (iobs.eq.0428)isite= 188			! Perth, Creg Bolt,            22 juin 2008
	if (iobs.eq.0429)isite= 189			! Reedy Creek, John Broughton, 22 juin 2008
	if (iobs.eq.0430)isite= 190			! Bankstown, Ted Dobosz,       22 juin 2008
	if (iobs.eq.0431)isite= 191			! Glenlee, Steve Kerr,         22 juin 2008
	if (iobs.eq.0432)isite= 141			! Les Makes,                   22 juin 2008
	if (iobs.eq.0433)isite= 118			! Sprinbok,                    22 juin 2008
	if (iobs.eq.0434)isite=	115			! Hakos,                       22 juin 2008

	if (iobs.eq.0440)isite=	122			! SAAO				Charon 22 juin 2008
	if (iobs.eq.0441)isite=	118			! Sprinbok
	if (iobs.eq.0442)isite=	115			! Hakos
	if (iobs.eq.0443)isite=	141			! Les Makes

	if (iobs.eq.0450)isite=	  1			! CFHT				Pluton 24 juin 2008

	if (iobs.eq.0460)isite=	 26			! Paranal			2002 MS4 30 juin 2009
	if (iobs.eq.0461)isite=	  3 		! La Silla
	if (iobs.eq.0462)isite=	146  		! San Pedro
	if (iobs.eq.0463)isite=	 89 		! Itajuba
	if (iobs.eq.0464)isite=	152 		! Ceamig B. Horizonte

	if (iobs.eq.0470)isite=	 94			! Lick				Pluton 25 aout 2008
	if (iobs.eq.0471)isite=	192 		! Grand Rapids		Pluton 25 aout 2008
	
	if (iobs.eq.0480)isite=	114 		! HESS				Varuna 19 fevrier 2010
	if (iobs.eq.0481)isite=	115 		! Hakos				Varuna 19 fevrier 2010
	if (iobs.eq.0482)isite=	193 		! Natal				Varuna 19 fevrier 2010
	if (iobs.eq.0483)isite=	194 		! Fortaleza/Quixada	Varuna 19 fevrier 2010
	if (iobs.eq.0484)isite=	195 		! Sao Luis			Varuna 19 fevrier 2010
	if (iobs.eq.0485)isite=	199			! Stellenbosch(RSA) Varuna 19 fevrier 2010
	if (iobs.eq.0486)isite=	196			! Recife/Camalau	Varuna 19 fevrier 2010
	if (iobs.eq.0487)isite=	113			! Tivoli			Varuna 19 fevrier 2010
	if (iobs.eq.0488)isite=	197			! Cecite			Varuna 19 fevrier 2010
	if (iobs.eq.0489)isite=	198			! Florianopolis		Varuna 19 fevrier 2010

	if (iobs.eq.0490)isite=	202			! Windhoek/MIT		Varuna 19 fevrier 2010
	if (iobs.eq.0491)isite=	203			! Maceio/MIT		Varuna 19 fevrier 2010

	if (iobs.eq.0500)isite=	  4			! Pic				Pluton 14 fevrier 2010
	if (iobs.eq.0501)isite=	200  		! Lu				Pluton 14 fevrier 2010
	if (iobs.eq.0502)isite=	201  		! Sisteron				Pluton 14 fevrier 2010

	if (iobs.eq.0510)isite=	 26			! VLT				Pluton 19 mai 2010
	if (iobs.eq.0511)isite=	  3 		! NTT				Pluton 19 mai 2010
	if (iobs.eq.0512)isite=	155 		! SOAR				Pluton 19 mai 2010

	if (iobs.eq.0515)isite=	157			! Mt John						 Pluton 04 juin 2010
	if (iobs.eq.0516)isite=	 14			! Hobart						 Pluton 04 juin 2010
	if (iobs.eq.0517)isite=	207 		! Stu Parker, Oxford NZ			 Pluton 04 juin 2010
	if (iobs.eq.0518)isite=	208 		! Bill Allen Vintage Lane NZ	 Pluton 04 juin 2010
	if (iobs.eq.0519)isite=	209 		! Dave Gault, West Wyalong       Pluton 04 juin 2010
	if (iobs.eq.0520)isite=	210 		! Hristo Pavlov, Jugiong		 Pluton 04 juin 2010

	if (iobs.eq.0530)isite=	 90			! Leoncito			Pluton 04 juillet 2010

	if (iobs.eq.0540)isite=	 89			! Itajuba				Ceres 17 aout 2010
	if (iobs.eq.0541)isite=	152			! CEAMIG				Ceres 17 aout 2010
	if (iobs.eq.0542)isite=	204 		! Ponta Grossa			Ceres 17 aout 2010
	if (iobs.eq.0543)isite=	205 		! INPE					Ceres 17 aout 2010
	if (iobs.eq.0544)isite=	206 		! Florianopolis, UFSC	Ceres 17 aout 2010

	if (iobs.eq.0545)isite=	146 		! 2003AZ84, 08/01/2011  Alain - San Pedro de Atacama
	if (iobs.eq.0546)isite=	3 		! 2003AZ84, 08/01/2011  TRAPPIST - La Silla

	if (iobs.eq.0550)isite=	212			! La Silla/Trappist		Eris 06 novembre 2010
	if (iobs.eq.0551)isite=	146			! San Pedro				Eris 06 novembre 2010
	if (iobs.eq.0552)isite=	 90			! El Leoncito			Eris 06 novembre 2010
	if (iobs.eq.0553)isite=	 89			! Itajuba				Eris 06 novembre 2010
	if (iobs.eq.0554)isite= 211			! Fortaleza				Eris 06 novembre 2010

	if (iobs.eq.0560)isite=	 89			! Itajuba				Eris **05** novembre 2010

	if (iobs.eq.0570)isite=	213			! NTT				Makemake, 23 avril 2011
	if (iobs.eq.0571)isite=  26			! VLT/UT4 (devrait etre UT3)	Makemake, 23 avril 2011
	if (iobs.eq.0572)isite=	 88			! Los Armazones			Makemake, 23 avril 2011
	if (iobs.eq.0573)isite=	146			! San Pedro de Atacama		Makemake, 23 avril 2011
	if (iobs.eq.0574)isite=	 89			! Itajuba			Makemake, 23 avril 2011
	if (iobs.eq.0575)isite=	212			! Trappist			Makemake, 23 avril 2011

	if (iobs.eq.0584)isite=	225			! Las Campanas OGLE          			Quaoar,   04 mai   2011
	if (iobs.eq.0585)isite=	  2			! Cerro Tololo 4m          			Quaoar,   04 mai   2011
	if (iobs.eq.0586)isite=	143			! Montevideo/UY          			Quaoar,   04 mai   2011
	if (iobs.eq.0587)isite=	224			! Brasilia - Paulo Cacella			Quaoar,   04 mai   2011
	if (iobs.eq.0588)isite=	152			! CEAMIG					Quaoar,   04 mai   2011
	if (iobs.eq.0589)isite=	 90			! Leoncito - AR					Quaoar,   04 mai   2011
	if (iobs.eq.0590)isite=	 88			! Armazones					Quaoar,   04 mai   2011
	if (iobs.eq.0591)isite=	146      		! San Pedro de Atacama		 		Quaoar,   04 mai   2011
	if (iobs.eq.0592)isite=	214			! Rivera, Uruguay				Quaoar,   04 mai   2011
	if (iobs.eq.0593)isite=	215			! Salto,  Uruguay				Quaoar,   04 mai   2011
	if (iobs.eq.0594)isite=	216			! Santa Martina					Quaoar,   04 mai   2011
	if (iobs.eq.0595)isite=	 89			! Itajuba					Quaoar,   04 mai   2011
	if (iobs.eq.0596)isite=	204			! Ponta Grossa					Quaoar,   04 mai   2011
	if (iobs.eq.0597)isite=	206			! Santa Catarina				Quaoar,   04 mai   2011
	if (iobs.eq.0598)isite=	205			! Sao Jose dos Campos - INPE			Quaoar,   04 mai   2011
	if (iobs.eq.0599)isite=	212			! TRAPPIST					Quaoar,   04 mai   2011

	if (iobs.eq.0600)isite=	212			! Trappist						Pluton 04 juin 2011
	if (iobs.eq.0601)isite=	146 		! San Pedro						Pluton 04 juin 2011
	if (iobs.eq.0602)isite=	 89			! Itajuba						Pluton 04 juin 2011
	if (iobs.eq.0603)isite=	216			! San Martina					Pluton 04 juin 2011
	if (iobs.eq.0610)isite=	212			! Trappist						Charon 04 juin 2011
	if (iobs.eq.0611)isite=	146 		! San Pedro						Charon 04 juin 2011
	if (iobs.eq.0612)isite=	 89			! Itajuba						Charon 04 juin 2011
	if (iobs.eq.0613)isite=	216			! San Martina					Pluton 04 juin 2011

*---
	if (iobs.eq.0620)isite=	218			! Kaui, Kekaha					Pluton 23 juin 2011
	if (iobs.eq.0621)isite=	220			! Maui, Lahaina					Pluton 23 juin 2011
	if (iobs.eq.0622)isite=	  1			! CFHT							Pluton 23 juin 2011
	if (iobs.eq.0623)isite=	221			! Kwajalein						Pluton 23 juin 2011
	if (iobs.eq.0624)isite=	222			! Nauru							Pluton 23 juin 2011
	if (iobs.eq.0625)isite=	223			! Majuro						Pluton 23 juin 2011

	if (iobs.eq.0630)isite=	218			! Kaui, Kekaha					Charon 23 juin 2011
	if (iobs.eq.0631)isite=	220			! Maui, Lahaina					Charon 23 juin 2011
	if (iobs.eq.0632)isite=	  1			! CFHT							Charon 23 juin 2011
	if (iobs.eq.0633)isite=	221			! Kwajalein						Charon 23 juin 2011
	if (iobs.eq.0634)isite=	222			! Nauru							Charon 23 juin 2011
	if (iobs.eq.0635)isite=	223			! Majuro						Charon 23 juin 2011

	if (iobs.eq.0640)isite=	218			! Kaui, Kekaha					Pluton 27 juin 2011
	if (iobs.eq.0641)isite=	220			! Maui, Lahaina					Pluton 27 juin 2011
	if (iobs.eq.0642)isite=	  1			! CFHT							Pluton 27 juin 2011
	if (iobs.eq.0643)isite=	221			! Kwajalein						Pluton 27 juin 2011
	if (iobs.eq.0644)isite=	222			! Nauru							Pluton 27 juin 2011
	if (iobs.eq.0645)isite=	223			! Majuro						Pluton 27 juin 2011

	if (iobs.eq.0650)isite=	219			! Darwin						Hydra  27 juin 2011

	if (iobs.eq.0660)isite=	226			! Alicante               , 2003 AZ84  03 fevrier 2012  
	if (iobs.eq.0661)isite=	227			! Mt Abu                 , 2003 AZ84  03 fevrier 2012
	if (iobs.eq.0662)isite=	228			! Weizmann/Israel        , 2003 AZ84  03 fevrier 2012
	if (iobs.eq.0663)isite=	229			! IUCAA Girawali/India   , 2003 AZ84  03 fevrier 2012
	if (iobs.eq.0664)isite=	230			! Liverpool/Canary       , 2003 AZ84  03 fevrier 2012
													 
	if (iobs.eq.0680)isite=	231			! Cote d'Azur/FR         ,  Quaoar  17 fevrier 2012
	if (iobs.eq.0681)isite=	232			! TAROT Grasse/FR        ,  Quaoar  17 fevrier 2012
	if (iobs.eq.0682)isite=	233			! JLX Valensole/FR       ,  Quaoar  17 fevrier 2012
	if (iobs.eq.0683)isite=	234			! Sposseti Belinzona/CH  ,  Quaoar  17 fevrier 2012

	if (iobs.eq.0685)isite=	235			! Samford Valley/AUS     , 2002 VE95 03 December 2015
	if (iobs.eq.0686)isite=	236			! Flynn/AUS              , 2002 VE95 03 December 2015
	if (iobs.eq.0687)isite=	237			! Murrumbateman/AUS      , 2002 VE95 03 December 2015

	if (iobs.eq.0688)isite=	91			! Mérida/Venezuela, Quaoar-July 09, 2013

	if (iobs.eq.0689)isite=	238			! Cairns/Aus (Tim), Sedna - January 13, 2013
	if (iobs.eq.0690)isite=	239			! Cairns/Aus (Joe), Sedna - January 13, 2013
						
	if (iobs.eq.0691)isite=	240			! Alice Springs/Aus - Eris August 29, 2013
	if (iobs.eq.0692)isite=	235			! Samford Valley/Aus - Eris August 29, 2013


	if (iobs.eq.0693)isite=	144			! Bosque Alegre/ARG - 2003 VS2(07/11/2014)
	if (iobs.eq.0694)isite=	143			! Montevidéo/URU    - 2003 VS2 (07/11/2014)
	if (iobs.eq.0695)isite=	241			! Rosário/ARG       - 2003 VS2(07/11/2014)
	if (iobs.eq.0696)isite=	242			! Santa Rosa/ARG    - 2003 VS2 (07/11/2014)
	if (iobs.eq.0697)isite=	213			! NTT/La Silla/CHI  - 2003 VS2 (07/11/2014)
	if (iobs.eq.0698)isite=	2			! PROMPT/Tololo/CHI - 2003 VS2 (07/11/2014)

	if (iobs.eq.0700)isite=	144			! Bosque Alegre
	if (iobs.eq.0701)isite=	241			! Busoniano
	if (iobs.eq.0702)isite=	143			! Los Molinos
	if (iobs.eq.0703)isite=	242			! Santa Rosa





c
c NB. lat en deg et lon (ouest) en heure, alt en m
c
	call sites       (isite, itopo,lat,lon,alt)
	call rayon_terre (lat,alt,itopo, rcphi,rsphi)
c	write(*,*) 'Rayon de la Terre a l''endroit considere:', sngl(dsqrt(rcphi**2+rsphi**2))

	call ajulien     (heure,jj,mm,aa, ajd,hs0,hs)
	hs_loc= hs - lon
! j'avais comment\8E la ligne ci-dessous, pourquoi? (permet d'eviter une heure siderale locale < 0 ou > 24):	
	if (hs_loc.lt.(0.d0)) hs_loc= hs_loc + 24.d0 ; if (hs_loc.gt.(24.d0)) hs_loc= hs_loc - 24.d0
	write(*,*)  'Heure siderale a Greenwich a 00h TU et a l''heure courante, heure siderale locale'
	call dms    (hs0,     ihsid0,    imsid0,    secsid0)
	call dms    (hs,      ihsid,     imsid,     secsid)
	call dms    (hs_loc,  ihsid_loc, imsid_loc, secsid_loc)
	write(*,'(3(I2.2,1x,I2.2,1x,f7.4,3x))')  ihsid0,imsid0,secsid0, ihsid,imsid,secsid, ihsid_loc,imsid_loc,secsid_loc

	call pole_terre  (ajd, alpha_p,delta_p)

	hs=           hs*radh
	alpha_p= alpha_p*radh
	delta_p= delta_p*radd

	cosap= dcos(alpha_p)
	sinap= dsin(alpha_p)
	cosdp= dcos(delta_p)
	sindp= dsin(delta_p)
c
c pole de la Terre dans repere J2000
c
	px= cosap*cosdp
	py= sinap*cosdp
	pz=       sindp
c
c X: vecteur perpendiculaire a p et dans le plan equatorial J2000
c
	xx= -sinap
	xy=  cosap
	xz=  0.d0
c
c Y: vecteur perpendiculaire a x et p, avec x, y, p direct
c
	yx= -cosap*sindp
	yy= -sinap*sindp
	yz=        cosdp
c
c Calcul du vecteur g, qui pointe du centre de la Terre a l'equateur, sur
c le meridien de Greenwich
c
	phi= hs - alpha_p
	sinphi= dsin(phi)
	cosphi= dcos(phi)

	gx= sinphi*xx - cosphi*yx
	gy= sinphi*xy - cosphi*yy
	gz= sinphi*xz - cosphi*yz
c
c Calcul du vecteur b, qui est perpendiculaire a p et g, avec g, b, p direct
c
	bx= cosphi*xx + sinphi*yx
	by= cosphi*xy + sinphi*yy
	bz= cosphi*xz + sinphi*yz
c
c Calcul du vecteur position de l'observateur, r, dans le repere J2000
c
	lon= lon*radh
	coslon= dcos(lon)
	sinlon= dsin(lon)

	rx= rcphi*coslon*gx - rcphi*sinlon*bx + rsphi*px
	ry= rcphi*coslon*gy - rcphi*sinlon*by + rsphi*py
	rz= rcphi*coslon*gz - rcphi*sinlon*bz + rsphi*pz
c
c Vecteur est, e, dans J2000
c
	cosa= dcos(alpha)
	sina= dsin(alpha)
	cosd= dcos(delta)
	sind= dsin(delta)

	ex= -sina
	ey=  cosa
	ez=  0.d0
c
c Vecteur nord celeste, n, dans J2000
c
	nx= -cosa*sind
	ny= -sina*sind
	nz=       cosd
c
c Calcul de la parallaxe
c
	ksip= rx*ex + ry*ey + rz*ez
	etap= rx*nx + ry*ny + rz*nz

	return
	end

********************************** fin de new_parallaxe ********************

c
c Sous programme de calcul de la latitude et du rayon geocentriques
c en fonction de la latitude topocentrique
c
	subroutine rayon_terre (lat,alt,itopo, rcphi,rsphi)
c
c	Entrees: lat, latitude topocentrique en degre decimal.
c		 alt, altitude en metres.
c
c	Sorties: rcphi, r geocentrique * cos(latitude geocentrique).
c		 rsphi, r geocentrique * sin(latitude geocentrique).
c
      implicit real*8 (a-h,o-z)
      real*8 lat,latp,latr
      rayon= 6378.16d+0
      pi= dacos(-1.d+00)
c
c si coordonnees topocentriques, correction a la latitude en arcsec
c
      latr= lat*(pi/180.d0)
	if (itopo.eq.1) then
 	  latp= -692.743d0*dsin(2.d0*latr) + 1.1633d0*dsin(4.d0*latr)
	  latp= latp - 0.0026d0*dsin(6.d0*latr)
	endif
	if (itopo.eq.0) then
	  latp= 0.d0
	endif
	if (itopo.ne.0.and.itopo.ne.1) then
	  write(*,*)  'itopo non defini!'
	  stop
	endif
c
	  latp= (latp/3600.d0 + lat)*(pi/180.d0)
	  rayonp= -0.000003519d0*dcos(4.d0*latr)+0.000000008d0*dcos(6.d0*latr)
	  rayonp= rayonp + 0.001676438d0*dcos(2.d0*latr) + 0.998327073d0
	  rayonp= rayon*rayonp
c
c
      rcphi= (rayonp + (alt/1000.d0))*dcos(latp)
      rsphi= (rayonp + (alt/1000.d0))*dsin(latp)
c
c
c
       write(*,*)'-----------------------------------------------'
       write(*,*) rcphi, rsphi
       write(*,*)'-----------------------------------------------'

      return
      end
c
c
c ********************** fin de rayon_terre ***********************************
c
	subroutine pole_terre  (ajd, alpha_p,delta_p)
        implicit real*8 (a-h,o-z)
c
c calcul des coordonnees J2000 du pole de la Terre. Formule obtenue sur
c le site du BdL
c
	ajd0= 2451545.d0			! JJ du 1.5 janvier 2000
	T= (ajd - ajd0)/36525.d0	! siecles juliens
	
	alpha_p= 0.d0 - 0.641d0*T
	alpha_p= alpha_p/15.d0		! transformation en heure
	alpha_p= alpha_p + 24.d0	! pour avoir alpha_p > 0

	delta_p= 90.d0 - 0.557d0*T	! degre

	call dms (alpha_p, ihap,imap,secap)
	call dms (delta_p, iddp,imdp,secdp)

	write(*,'(A,2(I2.2,1x,I2.2,1x,f7.4,3x))') 'pole Terre: ', ihap,imap,secap, iddp,imdp,secdp

	return
	end

************************* fin de pole_terre ***************************

	subroutine dms(heure, ihh,imn,sec)
c
c	ce sous-programme transforme l'heure decimale
c	'heure' en heure, minute, seconde: ihh:imn:sec.
c
	implicit real*8 (a-h,o-z)
	seuil= 1.d-3
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
c	
	return
	end
c
******************************* FIN DE DMS ********************************

	subroutine sites (isite, itopo, lat, lon, alt)
	implicit real*8 (a-h,o-z)
	real*8 lat, lon
c
c Sorties: lat en deg, lon en heure et alt en m
c
c itopo= 1, coordonnees topocentriques des observatoires.
c itopo= 0, coordonnees geocentriques  des observatoires.
c A verifier...
c
c	1= CFHT
c	2= CTIO
c	3= ESO La Silla
c	4= PIC
c	5= OHP
c	6= TEIDE
c	7= CATANIA
c	8= MEUDON
c	9= ANGERS
c	10= MANLEY
c	11= EIN HAROD
c	12= WISE
c	13= VATICAN
c	14= HOBART
c	15= HATFIELD
c	16= NEZEL
c	17= DENZAU
c	18= RGO
c	19= KITT PEAK
c	20= BUNDABERG
c	21= DUCABROOK
c	22= LOCHINGTON
c	23= BROWNSVILLE
c	24= CHILLAGOE
c	25= PRETORIA
c	26= ESO Paranal 
c	27= MONT MEGANTIC
c	28= ODESSA
c	29= CATALINA STATION
c	30= KITT PEAK
c	31= SALINAS
c	32= EL BUERAN (CUENCA)
c	33= MARACAIBO
c	34= BOBARES
c	35= ARVAL CARACAS
c	36= ARUBA
c	37= FORT DE FRANCE
c	38= TOBAGO
c	39= TRINIDAD
c	40= BRIDGETOWN, BARBADOS
c	41= ACORES
c	42= OEIRAS
c	43= PORTIMAO
c	44= COAA ALGARVE (PORTIMAO)
c	45= LINHACEIRA
c	46= ALVITO
c	47= GRANADA
c	48= ALCUBLAS
c	49= BORDEAUX
c	50= ST MAURICE (POITOU)
c	51= PEZENAS
c	52= MAUGUIO
c	53= NIMES
c	54= ORFEUILLES (MARSEILLE)
c	55= BINFIELD
c	56= WORTH HILL, BOURNEMOUTH
c	57= TEVERSHAM
c	58= PLATEAU D'ALBION
c	59= BD H. FABRE, MARSEILLE
c	60= GUITALENS
c	61= SABADELL
c   62= SALON
c	63= WELA, ARUBA
c	64= ST MAURICE CAZEVIEILLE
c	65= MARINHA GRANDE
c	66= ALMEIRIM
c	67= CARCAVELOS
c	68= SETUBAL
c	69= ALCACER DO SAL
c	70= BARCELONA 1
c	71= St SAVINIEN
c	72= St MARTIN DE CRAU
c	73= BARCELONA 2
c	74= ST ESTEVE
c	75= ALELLA
c	76= CASTELLON
c	77= HORTONEDA
c	78= BARCELONA 3
c	79= ESPLUGUES DE LLOBREGAT
c	80= ZARAGOZA
c	81= CALERN
c	82= DAX
c   83= ARICA
c	84= Jerusalem
c	85= Cumbaya
c	86= Lima
c	87= Mamina
c	88= Los Armazones
c	89= Itajuba
c	90= El Leoncito
c	91= CIDA Merida/ Venezuela
c	92= Lowell
c	93= Palomar
c	94= Lick
c	95= Tomar
c	96= Sabadell (Casas) 
c	97= Czarna Bialostocka 
c	98= Vitebsk
c	99= Kooriyama
c	100= Ooe
c	101= Kashiwa
c	102= Hitachi
c	103= Chichibu
c	104= Musashino
c	105= Mitaka
c	106= Abrera 
c	107= Roque de los Muchachos
c	108= New Jersey
c	109= Max Valier
c	110= Nyrola
c	111= Livermore
c	112= Windhoek
c	113= Tivoli
c	114= Hess
c	115= Hakos
c	116= Kleinbegin
c	117= Sandfontein
c	118= Springbok 
c	119= Nuwerus
c	120= Gifberg
c	121= Cederberg
c	122= SAAO Sutherland
c	123= SAAO Cape
c	124= Boyden
c	140= Maido
c	141= Makes
c	142= Fournaise
c	143= Montevideo
c	144= Bosque Alegre (Cordoba)
c   145= Asuncion
c   146= San Pedro de Atacama
c   147= Tarija
c   148= Marangani
c	149= Pico Veleta
c	150= WIRO
c	151= Wykrota, Bresil
c	152= CEAMIG-REA, Bresil
c	153= Patacamaya, Bolivie
c	154= Tilomonte, Chili
c	155= SOAR
c	156= IRSF
c	157= Mount John, N. Zealand
c	158= AAT
c	159= Dave Gault, Blue Mountains, occn Pluton 12 juin 06 et Pluton 22 jun 08 
c	160= Lade Blair, 12 juin 2006 (occn Pluton 12 juin 06)
c	161= Dome C, Antartique
c	162= Starhome (John Sanford)
c	163= Kitt Peak
c	164= Mt Bigelow, Kuiper 61"
c	165= Tenagra
c	166= Guanajuato
c	167= Mt Hopkins
c	168= Pinto Valley
c	169= Hereford
c	170= Cloudbait
c	171= Palmer Divide
c	172= Calvin-Rehoboth
c	173= Dark Sky Obs, Appalachian
c	174= Moore Obs.
c	175= George Obs., TX
c	176= San Pedro Martir
c	177= Mc Donald
c	178= Oklahoma Univ.
c	179= KASI (Corea tel.), Catalina Station
c	180= Mairinque
c   181= Canberra
c   182= Puimichel
c   183= Chatellerault
c   184= Monterrey, hacienda Zuazua
c   185= Sossusvlei
c   186= Perth, Lowell Telescope
c   187= Perth, Roger Groom, 22 juin 2008
c   188= Perth, Creg Bolt, 22 juin 2008
c   189= Reedy Creek, John Broughton, 22 juin 2008
c   190= Bankstown, Ted Dobosz, 22 juin 2008
c   191= Glenlee, Steve Kerr, 22 juin 2008
c   192= Grand Rapids, 25 aout 2008
c	193= Natal, 19 fevrier 2010
c	194= Fortaleza/Quixada, 19 fevrier 2010
c	195= Sao Luis, 19 fevrier 2010
c	196= Recife, 19 fevrier 2010
c	197= Cecite, 19 fevrier 2010
c	198= Florianopolis, 19 fevrier 2010
c	199= Stellenbosch, 19 fevrier 2010
c	200= Lu, 14 fevrier 2010
c	201= Sisteron, 14 fevrier 2010
c	202= Windhoek/MIT Varuna 19 fev. 2010
c	203= Maceio/MIT   Varuna 19 fev. 2010
c	204= Ponta Grossa,			Ceres 17 aout 2010
c	205= IMPE,					Ceres 17 aout 2010
c	206= Florianopolis, UFSC,	Ceres 17 aout 2010
c	207= Oxford NZ, Stu Parker,			Pluton 04 juin 2010
c	208= Vintage Lane NZ, Bill Allen,	Pluton 04 juin 2010
c	209= Dave Gault, West Wyalong,		Pluton 04 juin 2010
c	210= Hristo Pavlov, Jugiong			Pluton 04 juin 2010
c   211= Fortaleza, Dennis Weaver
c   212= La Silla, Trappist
c   213= La Silla, NTT
c   214= Rivera, Uruguay
c   215= Salto, Uruguay
c   216= Santa Martina
c   217= Santa Catarina
c   218= Kaui, Kekaha
c   219= Darwin
c	220= Lahaina
c	221= Kwajalein (iles Marshall)
c	222= Nauru     (iles Marshall)
c	223= Mjuro     (iles Marshall)
c       224= Brasilia/ Paulo Cacella
c       225= Las Campanas
c	226= Alicante		 
c	227= Mt Abu		 
c	228= Weizmann/Israel	 
c	229= IUCAA Girawali/India 
c	230= Liverpool/Canary					
c	231= Cote d'Azur/FR	 
c	232= TAROT Grasse/FR	 
c	233= JLX Valensole/FR	 
c	234= Sposseti Belinzona/CH
c	235= Samford Valley/AUS
c	236= Flynn/AUS
c	237= Murrumbateman/AUS
c	238= Cairns/AUS-Tim
c	239= Cairns/AUS-Joe
c	240= Alice Springs/AUS

c
c       Mauna Kea, CFHT:
c
	if(isite.eq.1) then
                itopo= 1
                ihho= 10
                imno= 21
                seco= 53.33
                iddo= 19
                immo= 49
                aseco= 33.9
                alt= 4205.
	endif
c
c	Cerro Tololo:
c
	if(isite.eq.2) then
		itopo= 1
		ihho=   04
		imno=   43
		seco=   15.633
		iddo=  -30
		immo=  -09
		aseco= -56.3
		alt= 2225.
	endif
c
c 	ESO La Silla: 3.6m
c
	if(isite.eq.3) then
		itopo= 1
		ihho=   04
		imno=   42
		seco=   55.61
		iddo=  -29
		immo=  -15
		aseco= -39.5
		alt= 2400.
	endif
c
c NB. 
c ESO LA Silla 3.6m: http://www.eso.org/sci/facilities/lasilla/telescopes/3p6/overview/index.html
c
c 70\A1 43' 54.1" W  -29\A1 15' 39.5" S  2400 m (WGS84)
c
c ESO La Silla NTT: http://www.eso.org/sci/facilities/lasilla/telescopes/ntt/overview/index.html
c
c 70\A1 44' 01.5" W  -29\A1 15' 32.1" S  2375 m (WGS84)
c
c ESO La Silla 2.2m: http://www.eso.org/sci/facilities/lasilla/telescopes/2p2/overview/index.html
c
c 70\A1 44' 12.0" W  -29\A1 15' 28.2" S alt ??? (WG84)
c
c Positions GPS Peter Sinclaire 17 avril 2011:
c 1m:    S29.25668 W70.73836=    29 15 24.048 S    70 44 18.096 W = 04h 42m57.2064s  W
c 2.2m:  S29.25785 W70.73665=    29 15 28.26  S    70 44 11.94  W = 04h 42m 56.796s  W
c NTT:   S29.25896 W70.73374=    29 15 32.256 S    70 44 01.464 W = 04h 42m 56.0976s W
c 3.6m:  S29.26097 W70.73169=    29 15 39.492 S    70 43 54.084 W = 04h 42m 55.6056s W
c 1.52m: S29.25591 W70.73894=    29 15 21.276 S    70 44 20.184 W = 04h 42m 57.3456s W
c
c Trappist, Mail E. Jehin 11 janvier 2011:
c                                29 15 16.59  S    70 44 21.82 W =  04h 42m 57.45466s W  alt= 2315 m

c	Pic du Midi:
c
	if(isite.eq.4) then
		itopo= 1
		ihho=   00
		imno=   00
		seco=  -34.16
		iddo=   42
		immo=   56
		aseco=  12.0
		alt= 2862.
	endif
c
c	Haute Provence:
c
	if(isite.eq.5) then
		itopo= 1
		ihho=   00
		imno=  -22
		seco=  -51.34
		iddo=   43
		immo=   55
		aseco=  45.6
		alt= 651.
	endif
c
c	Teide:
c
	if(isite.eq.6) then
		itopo= 1
		ihho=   01
		imno=   06
		seco=   01.333
		iddo=   28
		immo=   17
		aseco=  30.
		alt= 238.
	endif
c
c	Catania:
c
	if(isite.eq.7) then
		itopo= 1
		ihho=   00
		imno=  -59
		seco=  -54.8
		iddo=   37
		immo=   41
		aseco=  30.
		alt= 1725.
	endif
c
c	Meudon:
c
	if(isite.eq.8) then
		itopo= 1
		ihho=   00
		imno=  -08
		seco=  -55.5
		iddo=   48
		immo=   48
		aseco=  18.
		alt= 162.
	endif
c
c	Angers:
c
	if(isite.eq.9) then
		itopo= 1
		ihho=   00
		imno=   01
		seco=   38.87
		iddo=   47
		immo=   23
		aseco=  18.
		alt= 80.
	endif
c
c	Manley:
c
	if(isite.eq.10) then
		itopo= 1
		ihho=   00
		imno=   10
		seco=   57.
		iddo=   53
		immo=   14
		aseco=  43.1
		alt= 75.
	endif
c
c	Ein Harod:
c
	if(isite.eq.11) then
		itopo= 1
		ihho= -02
		imno= -21
		seco= -34.8
		iddo=  32
		immo=  33
		aseco= 35.
		alt= 20.
	endif
c
c	Wise:
c
	if(isite.eq.12) then
		itopo= 1
		ihho= -02
		imno= -19
		seco= -03.2
		iddo=  30
		immo=  35
		aseco= 48.
		alt= 900.
	endif
c
c	Vatican:
c
	if(isite.eq.13) then
		itopo= 1
		ihho= -00
		imno= -50
		seco= -36.4
		iddo=  41
		immo=  44
		aseco= 48.
		alt= 450.
	endif
c
c	Hobart:  source ????
c
c	if(isite.eq.14) then
c		itopo= 1
c		ihho=   14
c		imno=   10
c		seco=   16.46
c		iddo=  -42
c		immo=  -50
c		aseco= -57.3
c		alt=    287.
c	endif
c
c	Hobart: 16", email Wolfgang, 11 juillet 2006
c
	if(isite.eq.14) then
		itopo= 1
		ihho=  -09
		imno=  -49
		seco=  -43.636
		iddo=  -42
		immo=  -50
		aseco= -49.8
		alt=    272.
	endif
c
c	Hatfield:
c
	if(isite.eq.15) then
		itopo= 1
		ihho=   00
		imno=   00
		seco=   22.032
		iddo=   51
		immo=   46
		aseco=  27.8
		alt=    66.
	endif
c
c	Nezel:
c
	if(isite.eq.16) then
		itopo= 1
		ihho=  -00
		imno=  -36
		seco=  -44.133
		iddo=   52
		immo=   58
		aseco=  37.
		alt=    20.
	endif
c
c	Denzau:
c
	if(isite.eq.17) then
		itopo= 1
		ihho=  -00
		imno=  -28
		seco=  -17.007
		iddo=   51
		immo=   24
		aseco=  34.9
		alt=    94.
	endif
c
c	RGO: 
c
	if(isite.eq.18) then
		itopo= 1
		ihho=  -00 
		imno=  -01
		seco=  -20.6 
		iddo=   50
		immo=   52
		aseco=  18. 
		alt=    53.
	endif
c
c	Kitt Peak: 
c
	if(isite.eq.19) then
		itopo= 1
		ihho=   07 
		imno=   26 
		seco=   24. 
		iddo=   31
		immo=   57
		aseco=  48. 
		alt=    2120.
	endif
c
c	Bundaberg: (cf. message Wolfgang 10/7/97 + 18/2/98)
c
	if(isite.eq.20) then
		itopo= 1
		ihho=   -10
		imno=   -09
		seco=   -30.36
		iddo=   -24
		immo=   -56
		aseco=  -35.7
		alt=     10.
	endif
c
c	Ducabrook: (cf. messages Wolfgang 19/7/97, 10/2/98, Neilsen, 2/3/98)
c
	if(isite.eq.21) then
		itopo= 1
		ihho=   -09
		imno=   -49
		seco=   -46.668
		iddo=   -23
		immo=   -53
		aseco=  -55.02
		alt=    320.
	endif
c
c	Lochington: (cf. messages Wolfgang 19/7/97, 10/2/98, Neilsen, 2/3/98)
c
	if(isite.eq.22) then
		itopo= 1
		ihho=   -09
		imno=   -50
		seco=   -05.656
		iddo=   -23
		immo=   -56
		aseco=  -42.48
		alt=    270.
	endif
c
c	Brownsville: (cf. memo. http://www.lpl.arizona.edu/~rhill/planocc)
c
	if(isite.eq.23) then
		itopo= 1
		ihho=   06
		imno=   30
		seco=   08.75
		iddo=   25
		immo=   58
		aseco=  40.95
		alt=    -0.9
	endif
c
c	Chillagoe: (cf. message de Jim Elliot, 3/7/1998)
c
	if(isite.eq.24) then
		itopo= 1
		ihho=   -09
		imno=   -38
		seco=   -06.4333
		iddo=   -17
		immo=   -08
		aseco=  -57.4
		alt=      0.			! a verifier !
	endif
c
c	Pretoria (cf. message de Pierre Drossart avril 1999)
c
	if(isite.eq.25) then
		itopo= 1
		ihho=   -01
		imno=   -52
		seco=   -54.9 
		iddo=   -25
		immo=   -47
		aseco=  -18. 
		alt=   1542.			
	endif
c
c	Paranal UT4 (http://www.eso.org/paranal/site/paranal.html#GeoInfo)
c
	if(isite.eq.26) then
		itopo= 1
		ihho=   04
		imno=   41
		seco=   36.53
		iddo=  -24 
		immo=  -37
		aseco= -31.
		alt=   	2635.43	
	endif
c
c	MONT MEGANTIC 
c
	if(isite.eq.27) then
		itopo= 1
		ihho=  04 
		imno=  44
		seco=  36.8
		iddo=  45
		immo=  27
		aseco= 18.
		alt=   1114.	
	endif
c
c	ODESSA (cf. message de Horanic 22/12/1999)
c
	if(isite.eq.28) then
		itopo= 1
		ihho= -02
		imno= -01
		seco= -05.6
		iddo=  46
		immo=  23
		aseco= 48.
		alt=   50.
	endif
c
c	CATALINA STATION, MT LEMMON, 60" telescope
c   il semble qu'il y ait eu une erreur: ces sont en fait les coordonnees de Mt Bigelow, voir plus loin (11/4/07)
c   anciennes coord. mises par erreur: 07:22:55.6733= 110:34:55.1 32:25:00.7, z= 2510.
c   le site http://james.as.arizona.edu/~psmith/60inch/ donne les coordonnees corrigees ci-dessous:
c
	if(isite.eq.29) then
		itopo= 1
		ihho=  07
		imno=  23
		seco=  09.3333
		iddo=  32
		immo=  26
		aseco= 36.
		alt=   2790.
	endif
c
c	KITT PEAK
c
	if(isite.eq.30) then
		itopo= 1
		ihho=  07 
		imno=  26
		seco=  23.8533
		iddo=  31
		immo=  57
		aseco= 47.0
		alt=   2076.
	endif
c
c	SALINAS 
c
	if(isite.eq.31) then
		itopo= 1
		ihho=  05
		imno=  12 
		seco=  44.44
		iddo=  00
		immo=  28 
		aseco= 01.5
		alt=   2040.
	endif
c
c	CERRO EL BUERAN (CUENCA)
c
	if(isite.eq.32) then
		itopo= 1
		ihho=  05
		imno=  15
		seco=  53.61
		iddo=  -02
		immo=  -39
		aseco= -46.2828
		alt=   3987.
	endif
c
c	MARACAIBO
c
	if(isite.eq.33) then
		itopo= 1
		ihho=  04 
		imno=  46
		seco=  29.67
		iddo=  10
		immo=  42
		aseco= 53.
		alt=   0.
	endif
c
c	BOBARES
c
	if(isite.eq.34) then
		itopo= 1
		ihho=  04
		imno=  37
		seco=  49.37 
		iddo=  10
		immo=  16
		aseco= 19.5
		alt=   620.
	endif
c
c	ARVAL CARACAS
c
	if(isite.eq.35) then
		itopo= 1
		ihho=  04
		imno=  27
		seco=  22.63
		iddo=  10
		immo=  30
		aseco= 08.6
		alt=   915.
	endif
c
c	ARIKOK, ARUBA
c
	if(isite.eq.36) then
		itopo= 1
		ihho=  04
		imno=  39
		seco=  42.27
		iddo=  12
		immo=  29
		aseco= 55.6
		alt=   99. 
	endif
c
c	FORT DE FRANCE
c
	if(isite.eq.37) then
		itopo= 1
		ihho=  04
		imno=  04
		seco=  19.8
		iddo=  14
		immo=  37
		aseco= 04.
		alt=   0.
	endif
c
c	TOBAGO
c
	if(isite.eq.38) then
		itopo= 1
		ihho=  04
		imno=  03
		seco=  02.8
		iddo=  11
		immo=  13
		aseco= 54.
		alt=   0.
	endif
c
c	TRINIDAD	
c
	if(isite.eq.39) then
		itopo= 1
		ihho=  03 
		imno=  58 
		seco=  20.67 
		iddo=  13
		immo=  05
		aseco= 07.
		alt=   74.
	endif
c
c	BRIDGETOWN, BARBADOS
c
	if(isite.eq.40) then
		itopo= 1
		ihho=  03
		imno=  58
		seco=  20.67
		iddo=  13
		immo=  05
		aseco= 07.
		alt=   74.
	endif
c
c	PONTA DELGADA, ACORES
c
	if(isite.eq.41) then
		itopo= 1
		ihho=  01
		imno=  42
		seco=  42.53
		iddo=  37
		immo=  44
		aseco= 38.
		alt=   50.
	endif
c
c	OEIRAS
c
	if(isite.eq.42) then
		itopo= 1
		ihho=  00
		imno=  37
		seco=  17.67
		iddo=  38
		immo=  41
		aseco= 07.
		alt=   50.
	endif
c
c	PORTIMAO
c
	if(isite.eq.43) then
		itopo= 1
		ihho=  00
		imno=  34
		seco=  30.21
		iddo=  37
		immo=  08
		aseco= 28.7
		alt=   64.
	endif
c
c	COAA ALGARVE (PORTIMAO)
c
	if(isite.eq.44) then
		itopo= 1
		ihho=  00
		imno=  34
		seco=  24.12
		iddo=  37
		immo=  11
		aseco= 24.6
		alt=   65.
	endif
c
c	LINHACEIRA
c
	if(isite.eq.45) then
		itopo= 1
		ihho=  00
		imno=  33
		seco=  32.1
		iddo=  39
		immo=  31
		aseco= 22.6
		alt=   90.
	endif
c
c	ALVITO
c
	if(isite.eq.46) then
		itopo= 1
		ihho=  00
		imno=  32
		seco=  24.73
		iddo=  38
		immo=  11
		aseco= 01.8
		alt=   50. 
	endif
c
c	GRANADA
c
	if(isite.eq.47) then
		itopo= 1
		ihho=  00
		imno=  13
		seco=  32.33
		iddo=  37
		immo=  10
		aseco= 54.
		alt=   0.
	endif
c
c	ALCUBLAS
c
	if(isite.eq.48) then
		itopo= 1
		ihho=  00
		imno=  02
		seco=  45.6
		iddo=  39
		immo=  47
		aseco= 41.
		alt=   903.
	endif
c
c	BORDEAUX 
c
	if(isite.eq.49) then
		itopo= 1
		ihho=  00
		imno=  02
		seco=  06.8 
		iddo=  44
		immo=  50
		aseco= 07.
		alt=   73.
	endif
c
c 	ST MAURICE, POITOU
c
	if(isite.eq.50) then
		itopo= 1
		ihho=  -00
		imno=  -02
		seco=  -00. 
		iddo=  46
		immo=  22
		aseco= 00.
		alt=   120.
	endif
c
c	PEZENAS
c
	if(isite.eq.51) then
		itopo= 1
		ihho=  -00
		imno=  -13
		seco=  -31.07
		iddo=  43
		immo=  26
		aseco= 48.
		alt=   0.
	endif
c
c	MAUGUIO
c
	if(isite.eq.52) then
		itopo= 1
		ihho=  -00
		imno=  -15
		seco=  -48.47
		iddo=  43
		immo=  34
		aseco= 07.
		alt=   0.
	endif
c
c	NIMES
c
	if(isite.eq.53) then
		itopo= 1
		ihho=  -00
		imno=  -16
		seco=  -56.33
		iddo=  43
		immo=  48
		aseco= 05.
		alt=   170.
	endif
c
c	ORFEUILLES, MARSEILLES
c
	if(isite.eq.54) then
		itopo= 1
		ihho=  -00
		imno=  -21
		seco=  -51.67
		iddo=  43
		immo=  18
		aseco= 57.
		alt=   180.
	endif
c
c	BINFIELD
c
	if(isite.eq.55) then
		itopo= 1
		ihho=  00
		imno=  03
		seco=  09.3
		iddo=  51
		immo=  25
		aseco= 26.
		alt=   73.
	endif
c
c	WORTH MATRAVERS
c
	if(isite.eq.56) then
		itopo= 1
		ihho=  00
		imno=  08
		seco=  07.37
		iddo=  50
		immo=  35
		aseco= 52.9
		alt=   140.
	endif
c
c	TEVERSHAM
c
	if(isite.eq.57) then
		itopo= 1
		ihho=  -00
		imno=  -00
		seco=  -46.
		iddo=  52
		immo=  12
		aseco= 06.
		alt=   10.
	endif
c
c	PLATEAU D'ALBION 
c
	if(isite.eq.58) then
		itopo= 1
		ihho=  -00
		imno=  -20
		seco=  -16.36
		iddo=   44
		immo=   00
		aseco=  11.2
		alt=    1090.
	endif
c
c	BD H. FABRE, MARSEILLE 
c
	if(isite.eq.59) then
		itopo= 1
		ihho=  -00
		imno=  -21
		seco=  -39.53
		iddo=   43
		immo=   18
		aseco=  28.
		alt=    90.
	endif
c
c	GUITALENS
c
	if(isite.eq.60) then
		itopo= 1
		ihho=  -00
		imno=  -08
		seco=  -08.75 
		iddo=   43
		immo=   38 
		aseco=  34.7
		alt=    148.
	endif
c
c	SABADELL
c
	if(isite.eq.61) then
		itopo= 1
		ihho=  -00
		imno=  -08
		seco=  -21.93
		iddo=   41
		immo=   33 
		aseco=  03. 
		alt=    231.
	endif
c
c	SALON
c
	if(isite.eq.62) then
		itopo= 1
		ihho=  -00
		imno=  -20
		seco=  -24.
		iddo=   43
		immo=   36
		aseco=  00.
		alt=    40.
	endif
c
c	WELA, ARUBA
c
	if(isite.eq.63) then
		itopo= 1
		ihho=   04
		imno=   39
		seco=   51.31
		iddo=   12
		immo=   29
		aseco=  00.5
		alt=    30.
	endif
c
c	ST MAURICE CAZEVIEILLE
c
	if(isite.eq.64) then
		itopo= 1
		ihho=  -00
		imno=  -16
		seco=  -56.4
		iddo=   44
		immo=   00
		aseco=  36.
		alt=   190.
	endif
c
c	MARINHA GRANDE
c
	if(isite.eq.65) then
		itopo= 1
		ihho=  00
		imno=  36
		seco=  07.92
		iddo=  39
		immo=  45
		aseco= 00.
		alt=   30.
	endif
c
c	ALMEIRIM
c
	if(isite.eq.66) then
		itopo= 1
		ihho=  00
		imno=  34
		seco=  19.87
		iddo=  39
		immo=  10
		aseco= 50.
		alt=   10.
	endif
c
c	CARCAVELOS
c
	if(isite.eq.67) then
		itopo= 1
		ihho=  00
		imno=  37
		seco=  22.97
		iddo=  38
		immo=  41
		aseco= 11.4
		alt=   56.
	endif
c
c	SETUBAL
c
	if(isite.eq.68) then
		itopo= 1
		ihho=  00
		imno=  35
		seco=  40.
		iddo=  38
		immo=  30
		aseco= 00.
		alt=   40.
	endif
c
c	ALCACER DO SAL
c
	if(isite.eq.69) then
		itopo= 1
		ihho=  00
		imno=  33
		seco=  54.21
		iddo=  38
		immo=  21
		aseco= 38.9
		alt=   25.
	endif
c
c	BARCELONA 1
c
	if(isite.eq.70) then
		itopo= 1
		ihho=  -00
		imno=  -08
		seco=  -48.67
		iddo=  41
		immo=  25
		aseco= 18.9
		alt=   60.
	endif
c
c	ST SAVINIEN
c
	if(isite.eq.71) then
		itopo= 1
		ihho=  00
		imno=  02
		seco=  35.87
		iddo=  45
		immo=  53
		aseco= 26.5
		alt=   30.
	endif
c
c	ST MARTIN DE CRAU
c
	if(isite.eq.72) then
		itopo= 1
		ihho=  -00
		imno=  -19
		seco=  -16
		iddo=  43
		immo=  38
		aseco= 00.
		alt=   20.
	endif
c
c	BARCELONA 2
c
	if(isite.eq.73) then
		itopo= 1
		ihho=  -00
		imno=  -07
		seco=  -31.2
		iddo=  41
		immo=  23
		aseco= 32.4
		alt=   540.
	endif
c
c	ST ESTEVE
c
	if(isite.eq.74) then
		itopo= 1
		ihho=  -00
		imno=  -07
		seco=  -29.71
		iddo=  41
		immo=  29
		aseco= 41.5
		alt=   180.
	endif
c
c	ALELLA
c
	if(isite.eq.75) then
		itopo= 1
		ihho=  -00
		imno=  -09
		seco=  -12.12
		iddo=  41
		immo=  29
		aseco= 02.3
		alt=   45.
	endif
c
c	CASTELLON
c
	if(isite.eq.76) then
		itopo= 1
		ihho=  00
		imno=  00
		seco=  07.21
		iddo=  40
		immo=  00
		aseco= 31.5
		alt=   85.
	endif
c
c	HORTONEDA
c
	if(isite.eq.77) then
		itopo= 1
		ihho=  -00
		imno=  -04
		seco=  -10.33
		iddo=  42
		immo=  14
		aseco= 49.0
		alt=   1001.
	endif
c
c	BARCELONA 3
c
	if(isite.eq.78) then
		itopo= 1
		ihho=  -00
		imno=  -08
		seco=  -33.07
		iddo=  41
		immo=  23
		aseco= 06.4
		alt=   60.
	endif
c
c	ESPLUGUES DE LLOBREGAT
c
	if(isite.eq.79) then
		itopo= 1
		ihho=  -00
		imno=  -08
		seco=  -22.47
		iddo=  41
		immo=  22
		aseco= 38.0
		alt=  120.
	endif
c
c	ZARAGOZA
c
	if(isite.eq.80) then
		itopo= 1
		ihho=  00
		imno=  04
		seco=  10.03
		iddo=  41
		immo=  37
		aseco= 29.4
		alt=  330.
	endif
c
c	CALERN
c
	if(isite.eq.81) then
		itopo= 1
		ihho=  -00
		imno=  -27
		seco=  -42.4
		iddo=  43
		immo=  44
		aseco= 54.
		alt=  1270.
	endif
c
c	DAX, http://astrosurf.com/obsdax/home.htm
c
	if(isite.eq.82) then
		itopo= 1
		ihho=  00
		imno=  04
		seco=  07.32
		iddo=  43
		immo=  41
		aseco= 36.4
		alt=   11.
	endif
c
c	ARICA
c
	if(isite.eq.83) then
		itopo= 1
		ihho=  04
		imno=  39
		seco=  03.43
		iddo= -18
		immo= -26
		aseco=-53.8
		alt=   2500.
	endif
c
c	JERUSALEM, ECUADOR
c
	if(isite.eq.84) then
		itopo= 1
		ihho=  05
		imno=  13
		seco=  24.78
		iddo=  00
		immo=  00
		aseco= 03.342
		alt=   2314.
	endif
c
c	CUMBAYA
c
	if(isite.eq.85) then
		itopo= 1
		ihho=  05
		imno=  13
		seco=  43.98
		iddo= -00
		immo= -11
		aseco=-41.826
                alt=   2425.
	endif
c
c	LIMA (site M. Kretlow)
c
	if(isite.eq.86) then
		itopo= 1
		ihho=  05
		imno=  06
		seco=  07.84
		iddo= -11
		immo= -58
		aseco=-0.18
		alt=   2938.
	endif
c
c	MAMINA (D'apres message Marc Buie du 17 sept. 2002)
c
	if(isite.eq.87) then
		itopo= 1
		ihho=  04
		imno=  36
		seco=  48.0
		iddo= -20
		immo= -04
		aseco=-51.
                alt=   2870.
	endif
c
c	LOS ARMAZONES
c
	if(isite.eq.88) then
		itopo= 1
		ihho=  04
		imno=  40
		seco=  47.13
		iddo= -24
		immo= -35
		aseco=-51.
		alt=   3064.
	endif
c
c	ITAJUBA: http://www.lna.br/opd/opd_e.html
c
!	if(isite.eq.89) then
!		itopo= 1
!		ihho=  03
!		imno=  02
!		seco=  19.8
!		iddo= -22
!		immo= -32
!		aseco=-04.
!	    alt=   1864.
c
c From JPL Horizon system:
c
	if(isite.eq.89) then
		itopo= 1
		ihho=  03
		imno=  02
		seco=  19.833333
		iddo= -22
		immo= -32
		aseco=-07.8
	    alt=   1810.7
	endif
c
c	EL LEONCITO (http://www.casleo.gov.ar/Sitio/complejo.html)
c
c	if(isite.eq.90) then
c		itopo= 1
c		ihho=  04
c		imno=  37
c		seco=  12.
c		iddo= -31
c		immo= -47
c		aseco=-57.
c		alt=   2552.
c	endif
c
c
c	EL LEONCITO (position GPS Stefan Renner occn Charon 11/7/2005)
c
	if(isite.eq.90) then
		itopo= 1
		ihho=  04
		imno=  37
		seco=  10.996
		iddo= -31
		immo= -47
		aseco=-55.572
		alt=   2492.1
	endif
c
cc
c	MERIDA - In CIDA/Venezuela para Quaoar- July 09, 2013
c
	if(isite.eq.91) then
		itopo= 1
		ihho=  04
		imno=  43
		seco=  29.438
		iddo=  08
		immo=  47
		aseco= 25.80
		alt=   3600.
	endif
c
c	LOWELL
c
	if(isite.eq.92) then
		itopo= 1
		ihho=  07
		imno=  26
		seco=  39.19
		iddo=  35
		immo=  12
		aseco= 06.8
		alt=   2221.
	endif
c
c	PALOMAR 
c
	if(isite.eq.93) then
		itopo= 1
		ihho=  07
		imno=  47
		seco=  27.36
		iddo=  33
		immo=  21
		aseco= 21.6
		alt=   1713.14
	endif
c
c	LICK
c
	if(isite.eq.94) then
		itopo= 1
		ihho=  08
		imno=  06
		seco=  34.92
		iddo=  37
		immo=  20
		aseco= 24.6
		alt=   1281.39
	endif
c
c	TOMAR	
c
	if(isite.eq.95) then
		itopo= 1
		ihho=  00
		imno=  33
		seco=  32.1
		iddo=  39 
		immo=  31
		aseco= 22.6
                alt=   80.
	endif
c
c	SABADELL (CASAS) 
c
	if(isite.eq.96) then
		itopo= 1
		ihho= -00
		imno= -08
		seco= -21.93
		iddo=  41 
		immo=  33
		aseco= 04.
                alt=   225.
	endif
c
c	CZARNA BIALOSTOCKA 
c
	if(isite.eq.97) then
		itopo= 1
		ihho= -01
		imno= -34 
		seco= -22.85 
		iddo=  53 
		immo=  17  
		aseco= 53.8
                alt=   158.
	endif
c
c	VITEBSK
c
	if(isite.eq.98) then
		itopo= 1
		ihho= -02
		imno= -01 
		seco= -18.6 
		iddo=  55  
		immo=  03
		aseco= 40.
                alt=   180.
	endif
c
c       KOORIYAMA
c
        if(isite.eq.99) then
                itopo= 1
                ihho= -09
                imno= -21
                seco= -44.034
                iddo=  37 
                immo=  27
                aseco= 26.76
                alt=   260.
        endif
c
c       OOE
c
        if(isite.eq.100) then
                itopo= 1
                ihho= -09
                imno= -20
                seco= -35.493
                iddo=  38
                immo=  22
                aseco= 36.9
                alt=   182.
        endif
c
c	KASHIWA
c
        if(isite.eq.101) then
                itopo= 1
                ihho= -09
                imno= -19
                seco= -48.726
                iddo=  35
                immo=  50
                aseco= 48.3
                alt=   39.
        endif
c
c	HITACHI
c
        if(isite.eq.102) then
                itopo= 1
                ihho= -09
                imno= -22
                seco= -31.233
                iddo=  36
                immo=  32
                aseco= 19.22
                alt=   20.
        endif
c
c       CHICHIBU
c
        if(isite.eq.103) then
                itopo= 1
                ihho= -09
                imno= -16
                seco= -08.733
                iddo=  35
                immo=  57
                aseco= 53.
                alt=   355.
        endif
c
c       MUSASHINO
c
        if(isite.eq.104) then
                itopo= 1
                ihho= -09
                imno= -18
                seco= -15.433
                iddo=  35
                immo=  42
                aseco= 24.2
                alt=   66.
        endif
c
c       MITAKA 
c
        if(isite.eq.105) then
                itopo= 1
                ihho= -09
                imno= -18
                seco= -13.467
                iddo=  35
                immo=  41
                aseco= 53.
                alt=   62. 
        endif
c
c       ABRERA
c
        if(isite.eq.106) then
                itopo= 1
                ihho= -00
                imno= -07
                seco= -42.4 
                iddo=  41
                immo=  31
                aseco= 21.
                alt=   160. 
        endif
c
c       ROQUE DE LOS MUCHACHOS
c
        if(isite.eq.107) then
                itopo= 1
                ihho=  01
                imno=  11
                seco=  31.6
                iddo=  28
                immo=  45
                aseco= 36.
                alt=   2326. 
        endif
c
c       NEW JERSEY 
c
        if(isite.eq.108) then
                itopo= 1
                ihho=  04
                imno=  58
                seco=  22.8
                iddo=  39
                immo=  19
                aseco= 04.8
                alt=   0. 
        endif
c
c       MAX VALIER OBSERVATORY 
c
        if(isite.eq.109) then
                itopo= 1
                ihho= -00
                imno= -46
                seco= -00
                iddo=  46
                immo=  30
                aseco= 00.
                alt=   1350. 
        endif
c
c        NYROLA OBSERVATORY (positions prises dans JPL/Horizons)
c
        if(isite.eq.110) then
                itopo= 1
                ihho= -01
                imno= -42
                seco= -03.15 
                iddo=  62
                immo=  20
                aseco= 48.0
                alt=   205.215
        endif
c
c        LIVERMORE
c
        if(isite.eq.111) then
                itopo= 1
                ihho=  08
                imno=  06
                seco=  54.53
                iddo=  37
                immo=  43
                aseco= 08.58
                alt=    0.
        endif
c
c       WINDHOEK (SONJA'S FARM)
c	(mesures GPS B. Sicardy + laptop IOT, cf. cahier 8 nov 03)
c
        if(isite.eq.112) then
                itopo=  1
                ihho=  -01
                imno=  -08
                seco=  -26.1364 
                iddo=  -22
                immo=  -41
                aseco= -54.504 
                alt=   1920. 
        endif
c
c       TIVOLI 
c	(positionn donnee sur le site http://tivoli.ghhsNET.de/enzwo/framestart_e.html)
c
c remplace par (mail Reinhold Schreiber  24 fev. 2010, Varuna 19 fev. 2010)
c
        if(isite.eq.113) then
                itopo=  1
                ihho=  -01
                imno=  -12 
                seco=  -4.14667
                iddo=  -23
                immo=  -27
                aseco= -40.9
                alt=  1430. 
        endif
c
c	HESS 
c	(position donnee dans le site http://www.mpi-hd.mpg.de/hfm/HESS/HESS.html)
c
c        if(isite.eq.114) then
c                itopo=  1
c                ihho=  -01
c                imno=  -06 
c                seco=  -00. 
c                iddo=  -23
c                immo=  -16
c                aseco= -18.
c                alt=    1800. 
c        endif
c
c position inscrite dans headers fits de ATOM 
c mail Marcus Hauser 24 fev. 2010 (Varuna 19 fev. 2010)
c
        if(isite.eq.114) then
                itopo=  1
                ihho=  -01
                imno=  -06 
                seco=  -00.6696
                iddo=  -23
                immo=  -16
                aseco= -22.08
                alt=    1800. 
        endif
c
c	Hakos 
c 	(position donnee sur le depliant, cf. classeur)
c --->  chang'e le 27/oct/05 (email KL Bath)
c 23:14:42.0 S
c 16:21:12.0 E
c
c ---> chang'e le 2 mars 2010 en regardant la manipe Triton mai 2008:
c 16d 21.6922' E
c 23d 14.840 S
c alt. 1825. 
c
        if(isite.eq.115) then
                itopo=  1
                ihho=  -01
                imno=  -05
                seco=  -26.7688 
                iddo=  -23
                immo=  -14
                aseco= -11.
                alt=    1825. 
        endif
c
c	Kleinbegin Farm 
c 	(email Mike Kretlow 2/12/2003) 
c
        if(isite.eq.116) then
                itopo=  1
                ihho=  -01
                imno=  -16 
                seco=  -05.8532
                iddo=  -27
                immo=  -58
                aseco= -20.982
                alt=   500. 		! ????????
        endif
c
c	Sandfontein Farm 
c	(position GPS cahier 13/14 nov. 2003)
c
        if(isite.eq.117) then
                itopo=  1
                ihho=  -01
                imno=  -14
                seco=  -12.745 
                iddo=  -28
                immo=  -42
                aseco= -18.708
                alt=   515. 
        endif
c
c	Springbok 
c	(email Francois Colas 17 nov 2003)
c   confirme par Thomas Widemann 22 juin 2008:
c   17 deg 52.9820 E
c   29 deg 39.6716 S
c   alt. 951 m
c
        if(isite.eq.118) then
                itopo=  1
                ihho=  -01
                imno=  -11
                seco=  -31.928 
                iddo=  -29
                immo=  -39
                aseco= -40.296 
                alt=   951. 
        endif
c
c	Nuwerus 
c	(email Francois Colas 17 nov 2003) 
c
        if(isite.eq.119) then
                itopo=  1
                ihho=  -01
                imno=  -13
                seco=  -25.88 
                iddo=  -31
                immo=  -08
                aseco= -52.5 
                alt=   370. 
        endif
c
c	Gifberg 
c	(email Wolfgang Beisker 25 nov 2003) 
c
        if(isite.eq.120) then
                itopo=  1
                ihho=  -01
                imno=  -15 
                seco=  -07.9624
                iddo=  -31
                immo=  -48
                aseco= -33.45
                alt=   349.6 
        endif
c
c	Cederberg 
c	(email Wolfgang Beisker 25 nov 2003) 
c
        if(isite.eq.121) then
                itopo=  1
                ihho=  -01
                imno=  -17
                seco=  -00.6424
                iddo=  -32
                immo=  -29
                aseco= -58.416 
                alt=   873.1 
        endif
c
c	SAAO Sutherland 
c	(email Francois Colas 17 nov 2003) 
c
        if(isite.eq.122) then
                itopo=  1
                ihho=  -01
                imno=  -23
                seco=  -14.567 
                iddo=  -32
                immo=  -22
                aseco= -46.0 
                alt=  1760. 
        endif
c
c	SAAO Cape 
c	(email Francois Colas 17 nov 2003) 
c
        if(isite.eq.123) then
                itopo=  1
                ihho=  -01
                imno=  -13
                seco=  -54.547 
                iddo=  -33
                immo=  -56
                aseco= -04.5 
                alt=    15. 
        endif
c
c	Boyden (Bloemfontein) 
c	(position dans http://www.geocities.com/assabfn/spacetides/boyden.htm)
c
        if(isite.eq.124) then
                itopo=  1
                ihho=  -01
                imno=  -45
                seco=  -37.333 
                iddo=  -29
                immo=  -02
                aseco= -20. 
                alt=  1387. 
        endif
c
c	Maido 
c	(email JLX 01 nov. 2003)
c
c   Le 21 mai 2008, le C9 etait a (mail JLX 22/09/08)
c   +55d 23' 14".8 (vs. 55d 23' 15" ci-dessous)
c   -21d 04' 15".4, soit un vingtaine de m de difference, negligeable...
c
        if(isite.eq.140) then
                itopo=  1
                ihho=  -03
                imno=  -41
                seco=  -33. 
                iddo=  -21
                immo=  -04
                aseco= -16. 
                alt=  2205. 
        endif
c
c	Makes
c	(email JLX 01 nov. 2003)
c 21:11:57.0 S
c 55:24:35.0 E
c
c        if(isite.eq.141) then
c                itopo=  1
c                ihho=  -03
c                imno=  -41
c                seco=  -38.333 
c                iddo=  -21
c                immo=  -11
c                aseco= -57. 
c                alt=   976. 
c        endif
c
c mail E Frappa 16 oct. 2009
c
        if(isite.eq.141) then
                itopo=  1
                ihho=  -03
                imno=  -41
                seco=  -38.3 
                iddo=  -21
                immo=  -11
                aseco= -57.4
                alt=   972. 
        endif
c
c	Fournaise 
c	(email JLX 01 nov. 2003)
c
        if(isite.eq.142) then
                itopo=  1
                ihho=  -03
                imno=  -42
                seco=  -35 .733 
                iddo=  -21
                immo=  -13
                aseco= -53. 
                alt=  2350. 
        endif

c
c       Montevideo 
c       (site?)
c
        if(isite.eq.143) then
                itopo=  1
                ihho=   3
                imno=   44
                seco=   45.64 
                iddo=  -34
                immo=  -45
                aseco= -19.7
                alt=   88.8
        endif
c
c       Bosque Alegre
c	(site?)
c
        if(isite.eq.144) then
                itopo=  1
                ihho=   04
                imno=   18
                seco=   11.9
                iddo=  -31
                immo=  -35
                aseco= -54.0
                alt=  1213.1
        endif
c
c       Asuncion 
c	(email F. Doncel 6/4/05)
c
        if(isite.eq.145) then
                itopo=  1
                ihho=   03
                imno=   50
                seco=   05.08
                iddo=  -25
                immo=  -26
                aseco= -09.6
                alt=   000.
        endif
c
c       San Pedro de Atacama
c	(site?)
c
c        if(isite.eq.146) then
c                itopo=  1
c                ihho=   04
c                imno=   32
c                seco=   43.25
c                iddo=  -22
c                immo=  -57
c                aseco= -09.8
c                alt=  2347.
c        endif
c
c       San Pedro de Atacama
c	(email Erica Frappa 14/07/05))
c
        if(isite.eq.146) then
                itopo=  1
                ihho=   04
                imno=   32
                seco=   43.21
                iddo=  -22
                immo=  -57
                aseco= -08.4
                alt=  2410.
        endif
c
c       Tarija 
c	(site?)
c
        if(isite.eq.147) then
                itopo=  1
                ihho=   04
                imno=   19
                seco=   00. 
                iddo=  -21
                immo=  -31
                aseco= -00.
                alt=  2000.
        endif
c
c       Marangani
c	(www.concytec.gob.pe/space/obs.htm)
c
        if(isite.eq.148) then
                itopo=  1
                ihho=   04
                imno=   44
                seco=   42.93
                iddo=  -14
                immo=  -21
                aseco= -35.
                alt=  4060.
        endif
c
c       Pico Veleta
c	de JPL Horizon "Sierra Nevada Observatory"
c
        if(isite.eq.149) then
                itopo=  1
                ihho=   00
                imno=   13
                seco=   32.3267
                iddo=   37
                immo=   03
                aseco=  51.0
                alt=  2925.82
        endif
c
c       WIRO
c 	cf. email Bob Howell 17 nov. 2003 (BS 19 nov. 2003)
c
        if(isite.eq.150) then
                itopo=  1
                ihho=   07
                imno=   03
                seco=   54.2667
                iddo=   41
                immo=   05
                aseco=  50.
                alt=  2943.
        endif
c
c       Wykrota, Bresil
c 	cf. email Dennis Weaver 20/04/05
c
        if(isite.eq.151) then
                itopo=  1
                ihho=   02
                imno=   54
                seco=   45.67
                iddo=  -19
                immo=  -49
                aseco= -27.
                alt=  1500.
        endif
c
c       CEAMIG-REA, Bresil
c 	cf. email Dennis Weaver 20/04/05
c
        if(isite.eq.152) then
                itopo=  1
                ihho=   02
                imno=   55
                seco=   59.40
                iddo=  -19
                immo=  -49
                aseco= -49.
                alt=   825.
        endif
c
c       Patacamaya, Bolivie (100 km S. La Paz)
c 	cf. http://www.astro.edu.bo/pataca.htm
c
        if(isite.eq.153) then
                itopo=  1
                ihho=   04
                imno=   31
                seco=   48.47
                iddo=  -17
                immo=  -15
                aseco= -57.
                alt=   3789.
        endif
c
c   Tilomonte, 150 km sud S. Pedro de Atacama
c 	Position GPS Francois Colas
c
        if(isite.eq.154) then
                itopo=  1
                ihho=   04
                imno=   32
                seco=   21.55
                iddo=  -23
                immo=  -44
                aseco= -06.
                alt=   2373.
        endif
c
c       SOAR
c 	Cf. http://www.soartelescope.org/release/00home/eng_home.php
c
        if(isite.eq.155) then
                itopo=  1
                ihho=   04
                imno=   42
                seco=   56.09
                iddo=  -30
                immo=  -14
                aseco= -16.8
                alt=   2738.
        endif
c
c       IRSF 
c 	Cf. email Tetsuya Nagata 25 oct.05
c
        if(isite.eq.156) then
                itopo=  1
                ihho=  -01
                imno=  -23
                seco=  -14.5133
                iddo=  -32
                immo=  -22
                aseco= -48.0
                alt=   1760.
        endif
c
c       Mount John, New Zealand  
c 	Cf. http://www.phys.canterbury.ac.nz/research/mt_john/facilities.shtml
c
        if(isite.eq.157) then
                itopo=  1
                ihho=  -11
                imno=  -21
                seco=  -51.6
                iddo=  -43
                immo=  -59
                aseco= -12.0
                alt=   1029.
        endif
c
c       AAT Australie, Siding Spring
c 	Cf. site JPL
c
        if(isite.eq.158) then
                itopo=  1
                ihho=  -09
                imno=  -56
                seco=  -15.86
                iddo=  -31
                immo=  -16
                aseco= -37.4
                alt=   1164.5
        endif
c
c   Dave Gault, 12 juin 2006, Blue Mountains
c 	Cf. email W. Beisker 14 juin 2006
c   + mails D. Gault 22-24 jun 08
c   + mail  D. Gault 14 oct 2009 (dossier 22jun08) pour decimales
c
        if(isite.eq.159) then
                itopo=  1
                ihho=  -10
                imno=  -02
                seco=  -33.86
                iddo=  -33
                immo=  -39
                aseco= -51.9
                alt=   286.
        endif
c
c	Stockport, Lade Blair, 12 juin 2006, 
c	Cf. telephone Wolfgang 23 juin 2006
c 	puis mail Wolfgang 12 aout 2006, et mail Lade Blair 21 aout 2006
c
        if(isite.eq.160) then
c                itopo=  1
c                ihho=  -09
c                imno=  -13
c                seco=  -31.2
c                iddo=  -34
c                immo=  -25
c                aseco= -12.
c                alt=   0000.
                itopo=  1
                ihho=  -09
                imno=  -14
                seco=  -55.025333
                iddo=  -34
                immo=  -19
                aseco= -55.31
                alt=    142.
        endif
c
c	Dome C, Antartique
c	transmis par Danielle Briot, le 18 decembre 2006, cf. agenda papier
c
        if(isite.eq.161) then
                itopo=  1
                ihho=  -08
                imno=  -13
                seco=  -23.
                iddo=  -75
                immo=  -06
                aseco= -36.0
                alt=    3233.
        endif
c
c	Starhome, California 
c	transmis par John Sanford, mail 20 jan 2007 
c
        if(isite.eq.162) then
                itopo=  1
                ihho=   07
                imno=   55
                seco=   13.2
                iddo=   36
                immo=   13
                aseco=  28.
                alt=    700.
        endif
c
c	Kitt Peak 
c	GPS, 18 mars 2007 
c
        if(isite.eq.163) then
                itopo=  1
                ihho=   07
                imno=   26
                seco=   24.0568
                iddo=   31
                immo=   57
                aseco=  46.98
                alt=    2079.
        endif
c
c	Mt Bigelow, Kuiper tel. 61"
c	http://james.as.arizona.edu/~psmith/61inch/
c
        if(isite.eq.164) then
                itopo=  1
                ihho=   07
                imno=   22
                seco=   56.2867
                iddo=   32
                immo=   24
                aseco=  59.3
                alt=    2510.
        endif
c
c	Tenagra, site web
c
        if(isite.eq.165) then
                itopo=  1
                ihho=   07
                imno=   23
                seco=   30.9867
                iddo=   31
                immo=   27
                aseco=  44.4
                alt=    1312.
        endif
c
c	Guanajuato, Google Earth
c   ???
c
        if(isite.eq.166) then
                itopo=  1
                ihho=   06
                imno=   45
                seco=   01.3333
                iddo=   21
                immo=   01
                aseco=  10.
                alt=    2000.
        endif
c
c	Mount Hopkins, JPL site 696
c   MMT ???
c
        if(isite.eq.167) then
                itopo=  1
                ihho=   07
                imno=   23
                seco=   32.3067
                iddo=   31
                immo=   41
                aseco=  20.2
                alt=    2627.6
        endif
c
c	Pinto Valley, Morave Desert CA
c   Google Earth ???
c
        if(isite.eq.168) then
                itopo=  1
                ihho=   07
                imno=   45
                seco=   21
                iddo=   33
                immo=   46
                aseco=  27
                alt=    1000.
        endif
c
c	Hereford, AZ
c   http://reductionism.net.seanic.net/HD209458/ExoPlanet.html
c
        if(isite.eq.169) then
                itopo=  1
                ihho=   07
                imno=   20
                seco=   57.048
                iddo=   31
                immo=   27
                aseco=  07.9
                alt=    1420.
        endif
c
c	Cloudbait, CO
c   http://www.cloudbait.com/observatory.html
c
        if(isite.eq.170) then
                itopo=  1
                ihho=   07
                imno=   01
                seco=   56.068
                iddo=   38
                immo=   47
                aseco=  09.96
                alt=    2767.
        endif
c
c	Palmer Divide Observatory, CO, code UAI 716
c
        if(isite.eq.171) then
                itopo=  1
                ihho=   06
                imno=   59
                seco=   00.2667
                iddo=   39
                immo=   05
                aseco=  04.8
                alt=    2301.7
        endif
c
c Calvin-Rehoboth Observatory, NM
c http://www.calvin.edu/academic/phys/observatory/equipment/Rehoboth/
c et G98 JPL
c
        if(isite.eq.172) then
                itopo=  1
                ihho=   07
                imno=   14
                seco=   37.5504
                iddo=   35
                immo=   31
                aseco=  31.55
                alt=    2024.
        endif
c
c	Dark Sky Observatory, Appalachian, North Carolina, http://wikimapia.org/17823/
c
        if(isite.eq.173) then
                itopo=  1
                ihho=   05
                imno=   25
                seco=   39.3333
                iddo=   36
                immo=   15
                aseco=  09.
                alt=    1000.
        endif
c
c	Moore Obs, WA, coordonnees dans mail 23 mars (version html) 
c
        if(isite.eq.174) then
                itopo=  1
                ihho=   07
                imno=   56
                seco=   30.4333
                iddo=   46
                immo=   15
                aseco=  08.4
                alt=    125.0
        endif
c
c	George Obs, TX, coordonnees Google Earth?
c   ---> confirme par mail D. Dunham 19 fev 08, dossier Titania
c
        if(isite.eq.175) then
                itopo=  1
                ihho=   06
                imno=   22
                seco=   22.4667
                iddo=   29
                immo=   22
                aseco=  30.
                alt=    24.
        endif
c
c	San Pedro Martir, ref. ???
c
        if(isite.eq.176) then
                itopo=  1
                ihho=   07
                imno=   41
                seco=   51.2667
                iddo=   31
                immo=   02
                aseco=  39.
                alt=    2830.
        endif
c
c	McDonald, http://www.as.utexas.edu/mcdonald/mcdonald.html
c
        if(isite.eq.177) then
                itopo=  1
                ihho=   06
                imno=   56
                seco=   05.4267
                iddo=   30
                immo=   40
                aseco=  17.4
                alt=    2076.
        endif
c
c Oklahoma Univ.  http://www.observatory.ou.edu/
c Position JPL H30
c
        if(isite.eq.178) then
                itopo=  1
                ihho=   06
                imno=   29
                seco=   46.6067
                iddo=   35
                immo=   12
                aseco=  08.7
                alt=    381.6
        endif
c
c	CATALINA STATION,1m KASI (Coree)
c   cf. http://lemmon.kasi.re.kr/
c
		if(isite.eq.179) then
			itopo= 1
			ihho=  07
			imno=  23
			seco=  09.2667
			iddo=  32
			immo=  26
			aseco= 32.
			alt=   2776.
		endif
c
c	Mairinque
c   Mail Cristovao Jacques, 16 juin 2007
c
		if(isite.eq.180) then
			itopo= 1
			ihho=  03
			imno=  08
			seco=  50.6666
			iddo= -23
			immo= -34
			aseco=-10.
			alt=   933.
		endif
c
c	Canberra
c   Dave Herald, mail 02 aout 07 ---> correction 04 aout 07 
c
		if(isite.eq.181) then
			itopo= 1
			ihho= -09
			imno= -56
			seco= -15.2667
			iddo= -35
			immo= -23
			aseco=-49.3
			alt=   583.
		endif
c
c	Puimichel
c   UAI 184 (site Titania 8 septembre 2001, mail Demeautis 03 nov 2001) 
c
		if(isite.eq.182) then
			itopo= 1
			ihho= -00
			imno= -24
			seco= -08.667
			iddo=  43
			immo=  58
			aseco= 53.1
			alt=   714.
		endif
c
c	Chatellerault
c   (mail E. Bredner 29 oct 2001, dossier Titania)
c
		if(isite.eq.183) then
			itopo= 1
			ihho= -00
			imno= -02
			seco= -16.19
			iddo=  46
			immo=  50
			aseco= 21.0
			alt=   131.
		endif

c
c	Monterrey, hacienda Zuazua
c   (mail Pedro Valdes Sada, 18 fev 2008, dossier Titania)
c
		if(isite.eq.184) then
			itopo= 1
			ihho=  06
			imno=  40
			seco=  34.37333
			iddo=  25
			immo=  54
			aseco= 51.5
			alt=   385.
		endif

c
c	Sossusvlei
c   Source ???????
c
		if(isite.eq.185) then
			itopo=   1
			ihho=  -01
			imno=  -01
			seco=  -25.
			iddo=  -24
			immo=  -43
			aseco= -23.
			alt=   000.		! ???????
		endif
c
c	Perth, Lowell Telescope
c   Source JPL Horizons
c
		if(isite.eq.186) then
			itopo=   1
			ihho=  -07
			imno=  -44
			seco=  -32.4
			iddo=  -32
			immo=  -00
			aseco= -28.6
			alt=   428.
		endif
c
c	Perth, Roger Groom
c   Header images fits Pluton 22 juin 2008 et
c   mail du 30 juin 08 <--- position GPS
c
c   NB. la position dans le header des images 2006 etait probablement fausse
c   car sans doute prise dans Google Earth (???\A1
c
		if(isite.eq.187) then
			itopo=   1
			ihho=  -07
			imno=  -44
			seco=  -40.82
			iddo=  -31
			immo=  -52
			aseco= -53.0
			alt=   260.
		endif
c
c	Perth, Creg Bolt, 22 juin 2008
c   report (voir mail Dave Gault 26 juin 08 et graphique ds dossier Bolt) 
c
		if(isite.eq.188) then
			itopo=   1
			ihho=  -07
			imno=  -43
			seco=  -02.086667
			iddo=  -31
			immo=  -47
			aseco= -21.5
			alt=    45.
		endif
c
c	Reedy Creek, John Broughton, 22 juin 2008
c   report (voir dossier broughton) 
c
c		if(isite.eq.189) then
c			itopo=   1
c			ihho=  -10
c			imno=  -13
c			seco=  -35.52
c			iddo=  -28
c			immo=  -06
c			aseco= -30.3
c			alt=    66.
c		endif
c
c legerement modifie (~20 m) apres mail 15 oct 2009
c
		if(isite.eq.189) then
			itopo=   1
			ihho=  -10
			imno=  -13
			seco=  -35.4667
			iddo=  -28
			immo=  -06
			aseco= -29.9
			alt=    65.
		endif
c
c	Bankstown, Ted Dobosz, 22 juin 2008
c   mail 25 juin 2008 
c
		if(isite.eq.190) then
			itopo=   1
			ihho=  -10
			imno=  -04
			seco=  -07.
			iddo=  -33
			immo=  -55
			aseco= -56.
			alt=    24.9
		endif
c
c  Glenlee, Steve Kerr, 22 juin 2008
c  cf. mail (entre autres) de confirmation de John Talbot, 02 aug 08
c
		if(isite.eq.191) then
			itopo=   1
			ihho=  -10
			imno=  -02
			seco=  -00.05333
			iddo=  -23
			immo=  -16
			aseco= -09.6
			alt=    50.
		endif

c  Grand Rapids, Larry Molnar, 25 aout 2008
c  cf. mail 25 aout 2008, puis reference H62 prise sur site JPL
c
		if(isite.eq.192) then
			itopo=   1
			ihho=   05
			imno=   42
			seco=   21.1933333
			iddo=   42
			immo=   55
			aseco=  49.5
			alt=    252.5
		endif
c
c Natal, 19 fevrier 2010
c mail F. Colas 25 fevrier 2010
c 
		if(isite.eq.193) then
			itopo=   1
			ihho=   02
			imno=   26
			seco=   04.2
			iddo=  -06
			immo=  -15
			aseco= -46.002
			alt=   350.
		endif
c
c  Fortaleza/Quixada, 19 fevrier 2010
c  mail F. Colas 25 fev. 2010
c
		if(isite.eq.194) then
			itopo=   1
			ihho=   02
			imno=   36
			seco=   02.2396
			iddo=  -05
			immo=  -02
			aseco= -27.546
			alt=   521.
		endif
c
c  Sao Luis, 19 fevrier 2010
c  mail Francois Colas 20 fevrier 2010, puis 25 fev. 2010
c
		if(isite.eq.195) then
			itopo=   1
			ihho=   02
			imno=   56
			seco=   50.764
			iddo=  -02
			immo=  -35
			aseco= -35.88
			alt=     2.
		endif
c
c  Recife/Camalau, 19 fevrier 2010
c  rapport Audem\87rio Prazeres & Everaldo Faustino
c  confirme par mail du 26 fevrier 2010
c
		if(isite.eq.196) then
			itopo=   1
			ihho=   02
			imno=   27
			seco=   02.74
			iddo=  -07
			immo=  -54
			aseco= -29.76
			alt=   509.140
		endif
c
c  Observat\97rio Astron\99mico Genival Leite/CECITE-SEE-AL, 19 fevrier 2010
c  rapport Edmilson Souza Barreto, voir aussi http://oagll.blogspot.com/
c
		if(isite.eq.197) then
			itopo=   1
			ihho=   02
			imno=   23
			seco=   01.8
			iddo=  -09
			immo=  -38
			aseco= -29.7
			alt=    56.
		endif
c
c  Florianopolis, 19 fevrier 2010
c  rapport Felipe Braga, mail 25 fevrier 2010
c
		if(isite.eq.198) then
			itopo=   1
			ihho=   03
			imno=   14
			seco=   10.4667
			iddo=  -27
			immo=  -39
			aseco= -37.
			alt=    10.		! d'apres Google
		endif
c
c	Stellenbosch (RSA), 19 fevrier 2010
c	http://www.psychohistorian.org/astronomy/news/2010/20100220_varuna_occultation.php 
c   et mail Auke Slotegraf 24 fevrier 2010
c
        if(isite.eq.199) then
				itopo=  1
                ihho=  -01
                imno=  -15
                seco=  -17.4667 
                iddo=  -33
                immo=  -58
                aseco= -50.0 
                alt=   106. 
        endif
c
c  Lu, 14 fevrier 2010
c  Mail C. Olkin 20 mars 2010
c
		if(isite.eq.200) then
			itopo=   1
			ihho=  -00
			imno=  -41
			seco=  -28.022108
			iddo=   46
			immo=   37
			aseco=  26.3358
			alt=  1937.
		endif
c		
c  Sisteron, 14 fevrier 2010
c  Mail F Colas 15 juin 2010 + dossier
c
		if(isite.eq.201) then
			itopo=   1
			ihho=  -00
			imno=  -23
			seco=  -45.08667
			iddo=   44
			immo=   05
			aseco=  18.20
			alt=   634.
		endif
c
c  Windhoek/MIT, 19 fevrier 2010
c  mail A. Gulbis 21 juillet 2010
c
		if(isite.eq.202) then
			itopo=   1
			ihho=  -01
			imno=  -08
			seco=  -25.28
			iddo=  -22
			immo=  -26
			aseco= -09.
			alt=  1810.
		endif
c
c  Maceio/MIT, 19 fevrier 2010
c  mail A. Gulbis 21 juillet 2010
c
		if(isite.eq.203) then
			itopo=   1
			ihho=   02
			imno=   22
			seco=   57.72
			iddo=  -09
			immo=  -38
			aseco= -30.
			alt=    74.
		endif

c
c Ponta Grossa
c header image Ceres 17 aout 2010 --> plus precis avec le
c mail de Marcelo Emilio 18 mai 2011, dossier Eris
c
		if(isite.eq.204) then
			itopo=   1
			ihho=   03
			imno=   20
			seco=   23.76
			iddo=  -25
			immo=  -05
			aseco= -22.2
			alt=   909.
		endif
c
c INPE
c relatorio Andr\8E Milone Ceres 17 aout 2010
c
		if(isite.eq.205) then
			itopo=   1
			ihho=   03
			imno=   03
			seco=   26.933
			iddo=  -23
			immo=  -12
			aseco= -33.
			alt=   617.
		endif
c
c Florianopolis, UFSC
c relatorio Alexandre Mello Ceres 17 aout 2010
c
		if(isite.eq.206) then
			itopo=   1
			ihho=   03
			imno=   14
			seco=    5.369
			iddo=  -27
			immo=  -36
			aseco= -12.286
			alt=    20.
		endif

c
c Oxford NZ, Stu Parker
c mail 09 juin 2010 (occ. Pluton 04 juin 2010)
c
		if(isite.eq.207) then
			itopo=   1
			ihho=  -11
			imno=  -28
			seco=  -52.5212
			iddo=  -43
			immo=  -18
			aseco= -36.78
			alt=   221.
		endif
c
c Vintage Lane Obs. NZ, Bill Allen
c mail 04 juin 2010 (occ. Pluton 04 juin 2010)
c
		if(isite.eq.208) then
			itopo=   1
			ihho=  -11
			imno=  -35
			seco=  -21.381333
			iddo=  -41
			immo=  -29
			aseco= -36.27
			alt=    37.5
		endif

c
c Dave Gault, West Wyalong Pluton 04 juin 2010
c cf. report xls D. Gault
c
		if(isite.eq.209) then
			itopo=   1	
			ihho=  -09
			imno=  -48
			seco=  -56.973333
			iddo=  -33
			immo=  -55
			aseco= -25.5
			alt=    245.
		endif
c
c Hristo Pavlov, Jugiong, Pluton 04 juin 2010
c cf. report xlsH. Pavlov (different dans rapport Talbot!)
c
		if(isite.eq.210) then
			itopo=   1	
			ihho=  -09
			imno=  -53
			seco=  -16.1206666
			iddo=  -34
			immo=  -49
			aseco= -34.53
			alt=    335
		endif
c
c Fortaleza, Dennis Weaver
c position mail D. Weaver 18 janvier 2011 (Eris)
c
		if(isite.eq.211) then
			itopo=   1
			ihho=   02
			imno=   34
			seco=   01.8
			iddo=  -03
			immo=  -44
			aseco= -18.00
			alt=    38.
		endif
c
c ESO La Silla: Trappist
c mail Emmanuel Jehin 11 janvier 2011
c
	if(isite.eq.212) then
		itopo= 1
		ihho=   04
		imno=   42
		seco=   57.454666
		iddo=  -29
		immo=  -15
		aseco= -16.59
		alt= 2315.
	endif
c
c NTT La Silla
c Positions GPS, mail Peter Sinclaire 17 avril 2011:
c
	if(isite.eq.213) then
		itopo= 1
		ihho=   04
		imno=   42
		seco=   56.0976
		iddo=  -29
		immo=  -15
		aseco= -32.256
		alt= 2400.		! comme 3.6m ???
	endif
c
c Rivera, Uruguay
c e-mail Santiago Roland from GPS, 04/05/11
c
	if(isite.eq.214) then
		itopo= 1
		ihho=   03
		imno=   42
		seco=   23.979
		iddo=  -30
		immo=  -51
		aseco= -47.7
		alt=   200.			! APPROX Google Earth
	endif
c
c Salto, Uruguay
c Salto OLASU, Uruguai
c http://www.olasu.com.uy/ubicacion.php?idmenu=2
c Position mail Gonzalo Tancredi 6 mai 2011, confirme mail S. Bruzzone 10 mai 2011
c
	if(isite.eq.215) then
		itopo= 1
		ihho=   03
		imno=   51
		seco=   54.8620
		iddo=  -31
		immo=  -23
		aseco= -33.09
		alt=    37.				! APPROX Google Earth
	endif
c
c Santa Martina
c Position mail Rodrigo Andres Leiva 6 mai 2011
c
	if(isite.eq.216) then
		itopo= 1
		ihho=   04
		imno=   42
		seco=   08.26667
		iddo=  -33
		immo=  -16
		aseco= -09.
		alt=  1450.				
	endif
c
c Santa Catarina
c Position mail Marcelo Assafin, Quaoar 4 mai 2011
c
	if(isite.eq.217) then
		itopo= 1
		ihho=   03
		imno=   14
		seco=   05.369
		iddo=  -27
		immo=  -36
		aseco= -12.286
		alt=   20.				
	endif
c
c Kaui, Kehaha
c Pluton/Charon 23 juin 2011, mail Thomas 21 juin 2011
c
	if(isite.eq.218) then
		itopo= 1
		ldo=  +159 ; lmo=  +43 ; also=  +21.6						! longitude (W) en \A1 ' "
		iddo=  +21 ; immo= +58 ; aseco= +15.0						! latitude      en \A1 ' "
		alt=    22.													! altitude en m
		alongo= (dfloat(ldo) + dfloat(lmo)/60.d0 + also/3600.d0)/15.d0
	    call dms (alongo, ihho, imno,seco)							! longitude (W) en h m s
		if (ihho.le.0) then
		 imno= -imno ; seco= -seco
		endif
	endif
c
c Darwin
c Hydra 27 juin 2011, Google Earth
c
	if(isite.eq.219) then
		itopo= 1
		ldo=  -130 ; lmo=  -50 ; also=  -36.						! longitude (W) en \A1 ' "
		iddo=  -12 ; immo= -27 ; aseco= -38.						! latitude      en \A1 ' "
		alt=     0.													! altitude en m
		alongo= (dfloat(ldo) + dfloat(lmo)/60.d0 + also/3600.d0)/15.d0
	    call dms (alongo, ihho, imno,seco)							! longitude (W) en h m s
		if (ihho.le.0) then
		 imno= -imno ; seco= -seco
		endif
	endif
c
c Maui, Lahaina
c Pluton/Charon 23 juin 2011, Google Earth
c
	if(isite.eq.220) then
		itopo= 1
		ldo=  +156 ; lmo=  +40 ; also=  +56.82						! longitude (W) en \A1 ' "
		iddo=  +20 ; immo= +52 ; aseco= +42.03						! latitude      en \A1 ' "
		alt=    00.													! altitude en m
		alongo= (dfloat(ldo) + dfloat(lmo)/60.d0 + also/3600.d0)/15.d0
	    call dms (alongo, ihho, imno,seco)							! longitude (W) en h m s
		if (ihho.le.0) then
		 imno= -imno ; seco= -seco
		endif
	endif
c
c Kwajalein (iles Marshall)
c Pluton 23 juin 2011, Google Earth
c
	if(isite.eq.221) then
		itopo= 1
		ldo=  -167 ; lmo=  -22 ; also=  -11.12						! longitude (W) en \A1 ' "
		iddo=  +09 ; immo= +11 ; aseco= +10.34						! latitude      en \A1 ' "
		alt=     0.													! altitude en m
		alongo= (dfloat(ldo) + dfloat(lmo)/60.d0 + also/3600.d0)/15.d0
	    call dms (alongo, ihho, imno,seco)							! longitude (W) en h m s
		if (ihho.le.0) then
		 imno= -imno ; seco= -seco
		endif
	endif
c
c Nauru (iles Marshall)
c Pluton 23 juin 2011, Google Earth
c
	if(isite.eq.222) then
		itopo= 1
		ldo=  -166 ; lmo=  -55 ; also=  -50.86						! longitude (W) en \A1 ' "
		iddo=  -00 ; immo= -31 ; aseco= -21.86						! latitude      en \A1 ' "
		alt=     0.													! altitude en m
		alongo= (dfloat(ldo) + dfloat(lmo)/60.d0 + also/3600.d0)/15.d0
	    call dms (alongo, ihho, imno,seco)							! longitude (W) en h m s
		if (ihho.le.0) then
		 imno= -imno ; seco= -seco
		endif
	endif
c
c Majuro (iles Marshall)
c Pluton 23 juin 2011, Google Earth
c
	if(isite.eq.223) then
		itopo= 1
		ldo=  -171 ; lmo=  -11 ; also=  -08.22						! longitude (W) en \A1 ' "
		iddo=  +07 ; immo= +06 ; aseco= +59.41						! latitude      en \A1 ' "
		alt=     0.													! altitude en m
		alongo= (dfloat(ldo) + dfloat(lmo)/60.d0 + also/3600.d0)/15.d0
	    call dms (alongo, ihho, imno,seco)							! longitude (W) en h m s
		if (ihho.le.0) then
		 imno= -imno ; seco= -seco
		endif
	endif
	
c
c Brasilia
c Formulario Quaoar Paulo Cacella, 04/05/11
c
		if(isite.eq.224) then
			itopo=   1
			ihho=   03
			imno=   11  
			seco=   38.7067
			iddo=  -15
			immo=  -53
			aseco= -29.1
			alt=    1072.  !a confirmer !!!
		alongo= (dfloat(ldo) + dfloat(lmo)/60.d0 + also/3600.d0)/15.d0
	    call dms (alongo, ihho, imno,seco)							! longitude (W) en h m s
		if (ihho.le.0) then
		 imno= -imno ; seco= -seco
		endif		
	endif
c
c Las Campanas - OGLE JJ Kavelaars
c from wikipedia 29 00 54S 70 41 32W
c
	if(isite.eq.225) then
		itopo= 1
		ihho=   04
		imno=   42
		seco=   41.31
		iddo=  -29
		immo=  -00
		aseco= -54.00
		alt= 2380
	endif

c
c Alicante
c e-mail Ortiz data?
c
	if(isite.eq.226) then
		itopo= 1
		ihho=   -23
		imno=   -58
		seco=   -13.2
		iddo=   38
		immo=   28
		aseco=  33
		alt= 187
	endif


c	 if(isite.eq.226) then
c		 itopo= 1
c		 ldo=  -359; lmo=  -33 ; also=  -18					       ! longitude (W) en \A1 ' "
c		 iddo=  38 ; immo= 28 ; aseco= 33					       ! latitude      en \A1 ' "
c		 alt=	187.3									 ! altitude en m
c		 alongo= (dfloat(ldo) + dfloat(lmo)/60.d0 + also/3600.d0)/15.d0
c	     call dms (alongo, ihho, imno,seco) 						 ! longitude (W) en h m s
c		 if (ihho.le.0) then
c		  imno= -imno ; seco= -seco
c		 endif
c	 endif
c
c
c Mt Abu Observatory / India  
c  (e-mail Jose Ortiz 04/fev)
c
	if(isite.eq.227) then
		itopo= 1       ! longitude (W) en h,', ",  latitude en d,', ",  altitude en m
		ihho=  -04 
		imno=  -51 
		seco=  -07.133 
		iddo=  24 
		immo=  39 
		aseco= 10 
		alt= 1680
	endif

c
c Weizmann Institute of Sciences Martin S. Kraar observatory IAU/MPC code C78, Israel
c (email Ilan Manulis 5/2)
c
	if(isite.eq.228) then
		itopo= 1       ! longitude (W) en h,', ",  latitude en d,', ",  altitude en m
		ihho=  -02 
		imno=  -19 
		seco=  -15.051 
		iddo=  31 
		immo=  54 
		aseco= 29.1 
		alt= 107
	endif

c
c IUCAA Girawali Observatory -India- 
c (e-mail Jean Lecacheux 05/fev)
c
	if(isite.eq.229) then
		itopo= 1       ! longitude (W) en h,', ",  latitude en d,', ",  altitude en m
		ihho=  -04 
		imno=  -55 
		seco=  -22.733 
		iddo=  19 
		immo=  04 
		aseco= 21 
		alt= 1000
	endif

c
c http://telescope.livjm.uk/Info/TelInst/Spec/
c The Liverpool Telescope, La Palma Canary Island
c

	if(isite.eq.230) then
		itopo= 1       ! longitude (W) en h,', ",  latitude en d,', ",  altitude en m
		ihho=  01 
		imno=  11 
		seco=  31.013 
		iddo=  28 
		immo=  45 
		aseco= 44.8 
		alt= 2363
	endif

c
c Paolo Tanga    Phone: ++33(0)492003042
c Address: Laboratoire Lagrange, Observatoire de la C\F4te d'Azur
c	 BP4229 - 16304 Nice Cedex 4 - France		      
c e-mail: Paolo Tanga@oca.eu				      
c
	if(isite.eq.231) then
		itopo= 1       ! longitude (W) en h,', ",  latitude en d,', ",  altitude en m
		ihho=  -0
		imno= -29 
		seco= -03.147 
		iddo=  43 
		immo=  47 
		aseco= 22.2 
		alt= 385
	endif


c
c TAROT North automatic telescope
c email Alain Klotz / Eric Frappa 26/fev/12 rapor planoccult
c
c
	if(isite.eq.232) then
		itopo= 1       ! longitude (W) en h,', ",  latitude en d,', ",  altitude en m
		ihho=  -00 
		imno= -27 
		seco= -41.673  
		iddo=  43 
		immo=  45 
		aseco= 07.3 
		alt= 1270
	endif
	

c
c Plateau de Valensole, Haute-Provence - near VALENSOLE
c mail JLX 23/fev/2012  rapor planoccult Occ Quaoar
c

	if(isite.eq.233) then
		itopo= 1       ! longitude (W) en h,', ",  latitude en d,', ",  altitude en m
		ihho=  -0 
		imno= -24
		seco= -01.533  
		iddo=  43 
		immo=  51
		aseco= 52.4 
		alt= 622
	endif


c
c Stefano SPOSETTI - near Bellinzona
c email 26/fev/12 rapor ploanoccult Quaoar
c

	if(isite.eq.234) then
		itopo= 1       ! longitude (W) en h,', ",  latitude en d,', ",  altitude en m
		ihho=   -0
		imno=   -36
		seco=   -05.767
		iddo=   46
		immo=   13
		aseco=  53.2
		alt= 260
	endif

c
c Jonathan Bradshaw -  Samford Valley/AUS
c positive chord of occult by 2002 VE95  and Eris.
c

	if(isite.eq.235) then
		itopo= 1       ! longitude (W) en h,', ",  latitude en d,', ",  altitude en m
		ihho=  -10
		imno=  -11
		seco=  -23.53
		iddo=  -27	
		immo=  -22
		aseco= -07.00
		alt= 80
	endif

c
c John Newman - Flynn/AUS
c negative chord of occult by 2002 VE95 
c

	if(isite.eq.236) then
		itopo= 1       ! longitude (W) en h,', ",  latitude en d,', ",  altitude en m
		ihho=  -9
		imno=  -56
		seco=  -11.83
		iddo=  -35	
		immo=  -11
		aseco= -55.30
		alt= 657
	endif

c
c Dave Herald - Murrumbateman/AUS
c negative chord of occult by 2002 VE95 
c

	if(isite.eq.237) then
		itopo= 1       ! longitude (W) en h,', ",  latitude en d,', ",  altitude en m
		ihho=  -9
		imno=  -59
		seco=  -35.20
		iddo=  -34	
		immo=  -57
		aseco= -31.50
		alt= 594
	endif

c
c Tim Carruthers - Cairns/AUS
c positive chord of occultation by Sedna 
c

	if(isite.eq.238) then
		itopo= 1       ! longitude (W) en h,', ",  latitude en d,', ",  altitude en m
		ihho=  -9
		imno=  -42
		seco=  -53.47
		iddo=  -16	
		immo=  -58
		aseco= -6.8
		alt= 108
	endif
c
c Joseph Brimacombe - Cairns/AUS
c positive chord of occultation by Sedna 
c

	if(isite.eq.239) then
		itopo= 1       ! longitude (W) en h,', ",  latitude en d,', ",  altitude en m 
		ihho=  -9
		imno=  -43
		seco=  -3.31
		iddo=  -16	
		immo=  -54
		aseco= -26.3
		alt= 0
	endif

c
c William Hanna - Alice Springs/AUS
c positive chord of occultation by Eris 
c

	if(isite.eq.240) then
		itopo= 1       ! longitude (W) en h,', ",  latitude en d,', ",  altitude en m 
		ihho=  -8
		imno=  -55
		seco=  -32.212
		iddo=  -23	
		immo=  -42
		aseco= -73.79
		alt= 581
	endif

c
c Rosário/ARG
c positive chord for 2003 VS2 in November 07, 2014:
c
	if(isite.eq.241) then
		itopo= 1
		ihho=   4  
		imno=   2 
		seco=   35.28 
		iddo=  -32 
		immo=  -58
		aseco= -16
		alt=   31.		
	endif

c
c Santa Rosa/ARG
c positive chord for 2003 VS2 in November 07, 2014 (e-mail): 
c
	if(isite.eq.242) then
		itopo= 1
		ihho=  4 
		imno=  17 
		seco=  17.833 
		iddo=  -36
		immo=  -38 
		aseco= -16
		alt= 139.7		
	endif



	lon= (dfloat(ihho)*3600.d0 + dfloat(imno)*60. +  seco)/3600.d0
	lat= (dfloat(iddo)*3600.d0 + dfloat(immo)*60. + aseco)/3600.d0

	return
	end
*********************************** fin de sites ************************

	subroutine lecture (nmax, npt,t,flux)
	implicit real*8 (a-h,o-z)
	dimension t(nmax), flux(nmax)
	character*200 fichier	
1000	format(a)

 1	write(*,*) 'Fichier a lire?'
	read 1000, fichier
	open (unit=50,file=fichier,status='old',form='formatted',err=1)
	
	npt= 0
	do i= 1, nmax
	 read(50,*,end=2,err=2) t(i), flux(i)
	 npt= npt+1
	enddo
 2	continue
 	close(50)
	
	write(*,*) npt, ' points lus'
	if (npt.eq.nmax) write(*,*) 'Tous les points ne sont pas lus!'
	
	return
	end
*********************************** fin de lecture ************************	

