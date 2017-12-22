	implicit real*8 (a-h,o-z)
	real*8 K
	pi= dacos(-1.d0)
	arcsec= pi/(3600.d0*180.d0)
!
! Estime le diametre angulaire (mas) d'une etoile geante ou supergeante
! a partir de ses magnitudes V (ou B) et K. Precision typique de 10 a 20 %
!
! Source: G.R. van Belle (1999), Predicting stellar angular size,
! Publi. Astron. Soc. Pacific 111, 1515-1523.
!
	A_V= 0.669d0
	B_V= 0.223d0

	A_B= 0.648d0
	B_B= 0.220d0
	
	write (*,*) 'B, V et K?'
	read (*,*) B, V, K
	
	write (*,*) 'Distance? (km)'
	read (*,*) D
	
	diam_V= A_V + B_V*(V-K) -0.2d0*V
	diam_V= 10**(diam_V)

	diam_B= A_B + B_B*(B-K) -0.2d0*B
	diam_B= 10**(diam_B)
	
	write (*,*) 'Diametres deduits de B et V (mas):', sngl(diam_B), sngl(diam_V)
	write (*,*) 'Diametre en km pour  B et V:', sngl((diam_B*arcsec*D)/1.d3), sngl((diam_V*arcsec*D)/1.d3)
	
	
	stop
	end
