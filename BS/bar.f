c
c Ce programme genere un profil de diffraction par une barre infinie, en mode
c normal "0". 
c
c VERIFIER COMMENT EST CALCULE CHI2 (SIGMA=CSTE OU PROPORTIONNEL AU SIGNAL?)
c
c On peut aussi l'utiliser pour ajuster un profil de bord diffractant a un objet
c opaque (occultation par un gros corps). Dans ce cas, il faut donner une largeur 
c suffisamment grande a la bande (opaque) et n'ajuster que le bord gauche (immersion,
c mode "-1") ou le bord droit (emersion, mode "+1") au donnees.
c 
c L'ajustement se fait en calculant le chi2 en fonction du t0 de 1/2 occultation
c
c Si on ne veut pas ajuster avec des donnees, il suffit de ne pas le demander au debut
c et de ne pas demander non plus de reponse instrumentale. On ne demande alors qu'une valeur
c de t0 
c
c sorties:
c fort.20: ombre geometrique
c fort.21: ombre avec diffraction (convolee par la bande passante) et
c          par le diametre stellaire)
c fort.22: ombre lissee par la reponse instrumentale
c fort.23: ombre (fort.22) interpolee sur les points d'observation 
c fort.24: t0, chi2, nombre de points ajustes
c fort.25: rayon de l'etoile, chi2_min, npt fittes (NB. "append")
c fort.26: appele "chi2.dat": chi2 (eventuellement - nfit) (NB. "append"), voir par ex. donnees Hakos/Varuna 19 fev 2010
c fort.27: differences entre les points interpoles (fort.23) et les observations
c fort.28: distance perpendiculaire au centre de l'anneau, ombre lissee par diffraction et diam. stellaire,
c          i.e. equivalent a fort.21
c fort.29: ombre geometrique en distance (au lieu de temps), equivalent de fort.20
c fort.30: ombre geometrique en distance (au lieu de temps), ombre lisse par reponse instrumentale
c          i.e. equivalent a fort.22 
c fort.31: chi2's individuels (point par point)
c
c NB. si on veut integrer le flux sur un intervalle Dt, il suffit de donner
c une reponse instrumentale constante sur -Dt/2 a Dt/2. On obtient alors
c une moyenne glissante des points generes dans fort.21. Cette constante
c n'a pas a etre normalisee, car le programme s'occupe de normaliser l'aire
c de la reponse instrumentale a un. 
c
c Ex: si on a une integration du flux pendant 0.2 sec, il faut mettre dans
c le fichier de la reponse instrumentale:
c -0.1 0.5 ou -0.1 1.2
c +0.1 0.5    +0.1 1.2, etc...
c
c
        implicit real*8 (a-h,o-z)
        parameter (nmax=20000)
        real*8 i_ampli
        dimension t(-nmax:nmax), flux(-nmax:nmax)
        dimension trep(nmax), rep(nmax)
        dimension tl(nmax),   fluxl(nmax)
        dimension tobs(nmax), fobs(nmax)
        character*120 observ
 1000   format(a)


c*-----------------------------------------------------------------------------
c        write(*,*) 'Faut il supprimer le fichier chi2.dat deja present? '
c        read(*,*) suppr_chi2
c        if (suppr_chi2.eq.1) call system ('rm -f chi2.dat')
        call system ('rm -f chi2.dat')
*-----------------------------------------------------------------------------
        open (unit=26,file='chi2.dat',status='unknown',form='formatted',position='append')

*-----------------------------------------------------------------------------
        write(*,*) 'Mode:'
        write(*,*) '-1= t0 for immersion'
        write(*,*) '+1= t0 for emersion'
        write(*,*) ' 0= t0 for the center of a square well'
        write(*,*) ' 2= two square well at a time'
        read(*,*) imod
*-----------------------------------------------------------------------------
c reponse instrumentale
        write(*,*) 'Exposure time:'
        read(*,*) exptime
        trep(1)=-exptime/2.d0
        trep(2)=exptime/2.d0
        rep(1)=1.d0
        rep(2)=1.d0
        nrep=2

*-----------------------------------------------------------------------------
        write(*,*) 'Name of the file containing observations data (time x norm. flux): '
        read 1000, observ
        open(unit=19,file=observ,status='old',form='formatted')
        write(*,*)'Sigma of the observation: '
        read(*,*) sigma
        write(*,*) 'Time limits to fix the number of points to fit the model: '
        read (*,*) tobs_min_impo, tobs_max_impo
        nobs= 0
        do i= 1, nmax
          read(19,*,end=102,err=103) ttt, fff
          if (ttt.ge.tobs_min_impo.and.ttt.le.tobs_max_impo) then
            nobs= nobs + 1
            tobs(nobs)= ttt
            fobs(nobs)= fff
          endif
103       continue
        enddo
102     continue
        write(*,*) nobs, ' points read at the file'
*-----------------------------------------------------------------------------
        write(*,*) 'Wavelength? Filter (detector) wavelength width? (microns)'
        read(*,*) wvlngth, dlambda
        wvlngth= wvlngth*1.d-09        ! en km
        dlambda= dlambda*1.d-09        ! en km

        write(*,*) 'Distance geocentric-object ? Star radii ? (km)'
        read(*,*) dist, re 

        write(*,*) 'Normal velocity of the star? (at the plane of the sky, km/s)'
        read(*,*) vn

        write(*,*) 'Width of the well (put 10000 or more if infinie, km)? Transmission (0 for opaque)?'
        read(*,*) width, trans

        write(*,*) 'Duration to analyse and step? (sec) (must be bigger than time interval above)'
        read(*,*) duree, pas

        write(*,*) 'Stellar flux outside occultation, stellar flux during occultation?'
        read(*,*) phi1, phi0

        write(*,*) 'Reference instant? (sec)'
        read(*,*) t0_ref
        t_milieu= t0_ref + width/(2.d0*vn)

        write(*,*) 'Number of times to explore around reference time? (sec)'
        read(*,*) nheure, pas_heure
        t0_min= t0_ref - pas_heure*dfloat(nheure/2)
        t0_max= t0_ref + pas_heure*dfloat(nheure/2)
        t0= t0_min


        x3= 0.d0     ;  x4= 0.d0   ; opa_ampli2= 0.d0                ! les parametres de la bande 2 sont neutralises, 
        width2= 0.d0 ;  trans2= 1.d0 ; deltat= 0.d0                  ! sauf avis contraire
        if (imod.eq.2) then
          write(*,*) 'Second well width? (km), Transmission (0 for opaque)? Shift in time from center of first well? (sec)'
          read(*,*) width2, trans2, deltat
        endif

        write(*,*) 'Sigma (0= sigma constant, 1= sigma proportional to signal)'
        read(*,*) isigma

        if (imod.eq.-1) then
          x1= 0.d0                 ! bord gauche de l'ombre cale
          x2= width                ! sur l'origine (ie t0 en temps)
        endif

        if (imod.eq.0) then
          x1= -width/2.d+00         ! ombre centree sur le milieu
          x2= +width/2.d+00         ! de la bande
        endif

        if (imod.eq.1) then
          x1= -width                ! bord droit de l'ombre cale
          x2= 0.d+00                ! sur l'origine (ie t0 en temps)
        endif

        if (imod.eq.2) then
          x1= -width/2.d+00                ! ombre centree sur le milieu
          x2= +width/2.d+00                ! de la bande principale
          x3= deltat*vn - width2/2.d+00    ! ombre centree sur le milieu
          x4= deltat*vn + width2/2.d+00    ! de la bande secondaire
        endif

        npt= int( duree/(2.d+00*pas) )
        opa_ampli=  1.d+00 - dsqrt(trans )
        opa_ampli2= 1.d+00 - dsqrt(trans2)


c-------------------------- debut de l'exploration en temps ------------------
        do while (t0.le.t0_max)
c-----------------------------------------------------------------------------
c
c Trace de l'ombre geometrique
c
          call ombre_geometrique (imod,t0,duree,width,trans,width2,trans2,deltat,vn,phi1,phi0)


c-----------------------------------------------------------------------------
c convolution par le diametre stellaire
          som= 0.
          do i= -npt, npt
            x= vn*pas*dfloat(i)
            wvlngth1= wvlngth - dlambda/2.d+00
            wvlngth2= wvlngth + dlambda/2.d+00

c         etoile: sous-programme de convolution par le diametre stellaire.
c         NB.:    on appelle deux fois etoile pour convoler par la bande passante.
            call etoile (imod,re,x1,x2,x3,x4,opa_ampli,opa_ampli2,wvlngth1,dist,x, flux1)
            call etoile (imod,re,x1,x2,x3,x4,opa_ampli,opa_ampli2,wvlngth2,dist,x, flux2)
            flu=   ( flux1 + flux2 )/2.d+00
            flux(i)= flu*(phi1-phi0) + phi0
            t(i)= t0 + dfloat(i)*pas
            write(21,*)  t(i),        sngl(flux(i))
            write(28,*) (t(i)-t0)*vn, sngl(flux(i)) 
            som= som + pas*vn*(1.d+00-flux(i))
          enddo
          write (29,*) -100., 1. ; write (29,*) -width/2.d0, 1. 
          write (29,*) -width/2.d0, sngl(phi0+(phi1-phi0)*trans)
          write (29,*) +width/2.d0, sngl(phi0+(phi1-phi0)*trans)
          write (29,*) +width/2.d0, 1. ; write (29,*) +100., 1.
          close(21) ; close(28) ; close(29)

c-----------------------------------------------------------------------------
c
c Convolution par la reponse instrumentale
          call instrument(nmax,t,flux,npt,trep,rep,nrep, tl,fluxl,nptl)

c
c Ecriture du modele final (apres convolution par etoile et instrument)
          som= 0.
          tl_min=  1.d50
          tl_max= -1.d50
          do i= 1, nptl
            write(22,*)  tl(i)       ,      fluxl(i)
            write(30,*) (tl(i)-t0)*vn, sngl(fluxl(i)) 
            if (tl(i).ge.tl_max) tl_max= tl(i)
            if (tl(i).le.tl_min) tl_min= tl(i)
            som= som + pas*vn*(1.d+00-fluxl(i))
          enddo
          close(22)
c-----------------------------------------------------------------------------
c Calcul du chi2 avec les donnees
c
          nfit= 0
          chi2= 0.d0
          tobs_min=  1.d50
          tobs_max= -1.d50
          do i= 1, nobs
            fac= (tobs(i)-tl_min)*(tobs(i)-tl_max)
            if (fac.le.(0.d0)) then
              call interpol (nmax,nptl,tl,fluxl,tobs(i), fmod_inter)
              write(23,*) sngl(tobs(i)), sngl(fmod_inter)
              write(27,*) sngl(tobs(i)), sngl(fobs(i)-fmod_inter)
              if (isigma.eq.0) sigma_local= sigma                        ! ATTENTION !cas ou sigma est constant!
              if (isigma.eq.1) sigma_local= sigma*fobs(i)                ! ATTENTION: cas ou sigma est proportionnel au signal!
              chi2= chi2 + ((fmod_inter - fobs(i))/sigma_local)**2
              write(31,*) sngl(tobs(i)), sngl(((fmod_inter - fobs(i))/sigma_local)**2)
              nfit= nfit + 1
              if (tobs(i).lt.tobs_min) tobs_min= tobs(i)
              if (tobs(i).gt.tobs_max) tobs_max= tobs(i)
            endif
          enddo
          close(23) ; close(27)
          write(*,*) 't0, chi2, nfit:', t0, sngl(chi2), nfit
          write(24,*) t0, chi2, nfit
          write(26,*) t0, chi2

          t0= t0 + pas_heure       ! on incremente t0
        enddo

c-------------------------- fin de l'exploration en temps ------------------
c-----------------------------------------------------------------------------


c on cherche le temps correspondant au minimum de chi2:

        close (24)

        call error_bar (nmax,1.d0)        ! calcul 1-sigma
        call error_bar (nmax,9.d0)        ! calcul 3-sigma

        write (*,*) 'Fresnel scale (km):', sngl(dsqrt(wvlngth*dist/2.d0))
        write (*,*) 'Well width (perpendicular, km), transmission T at Earth level'
        write (*,*) 'opacity p=1-sqrt(T) at the ''ring'' plane:'
        write (*,*) sngl(width), sngl(trans), sngl(opa_ampli)
        if (isigma.eq.0) write (*,*) 'sigma constant'
        if (isigma.eq.1) write (*,*) 'sigma proportional to signal '
        
        stop
        end

c                        FIN DE MAIN
c ***********************************************************************
c
c
        subroutine bar(x1,x2,opa_ampli,wvlngth,dist,x, flux)
c
c This subroutine gives the complex amplitude resulting from diffraction of a 
c coherent planar wave, on a homogeneous semi-transparent stripe. The stripe is
c assumed to be infinite in the 0y direction, it begins at x= "x1" and ends
c at x= "x2" (x1 < x2). It removes a fraction "opa_ampli" (opacity in amplitude)
c of the amplitude of the incident wave. So, opa_ampli= 0 corresponds to a
c transparent stripe, while opa_ampli= 1 corresponds to an opaque stripe. This
c opacity is related to the fractional transmission in intensity, F, by:
c
c                        opa_ampli= 1 - sqrt(F)
c
c The subroutine call the FRESNEL subroutine, which uses dimensionless 
c quantities normalized to the Fresnel scale fr= sqrt( (wvlngth*dist)/2 ),
c where "wvlngth" is the wavelength of the observation, and "dist" is the
c distance of the observer to the diffracting object (here the stripe).
c
c The wave is assumed to have an amplitude equal to unity at the stripe,
c and the subroutine gives the complex amplitude recorded by the observer,
c whose abscissa along the 0x axis is "x". The real part of the amplitude
c is Re= 1+r_ampli, and the imaginary part is Im= i_ampli (i_ampli to be 
c declared REAL*8 !). Caution: the resulting intensity is flux= Re**2 + Im**2,
c NO sqrt !
c
c The formula which is used is:
c
c amplitude(x)= opa_ampli*( (i-1)/2 )*{ C(x-x1) - C(x-x2) + 
c                            i*[ S(x-x1) - S(x-x2) ] },
c
c where x,x1,x2 have been normalized to the Fresnel scale, i*i= -1, and
c C and S are given by the subroutine FRESNEL.
c
c
c NB. x1,x2,wvlngth,dist and x must be expressed in the same unity when entered
c in the subroutine !
c
c
c
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 i_ampli
c
c Required accuracy on the amplitude:
        eps= 1.d-16

c Normalization to the Fresnel scale
        fr= dsqrt( (wvlngth*dist)/2.d+00 )
        x1fr= x1/fr
        x2fr= x2/fr
        xfr=   x/fr

c Calculation of the amplitude
        xmx1= xfr-x1fr
        xmx2= xfr-x2fr
        call FRESNEL (xmx1,eps,c1,s1)
        call FRESNEL (xmx2,eps,c2,s2)
        cc= c1-c2
        ss= s1-s2
        r_ampli= - (cc+ss)*(opa_ampli/2.d+00)
        i_ampli=   (cc-ss)*(opa_ampli/2.d+00)

        amplir= 1.d+00 + r_ampli 
        amplii= i_ampli
        flux= amplir*amplir + amplii*amplii

        return
        end
c                                end of bar
c ***********************************************************************

        subroutine bar2(x1,x2,x3,x4,opa_ampli,opa_ampli2,wvlngth,dist,x, flux)
c
c same as subroutine bar, except that there are **two** diffracting bands
        IMPLICIT REAL*8 (A-H,O-Z)
        REAL*8 i_ampli
c
c Required accuracy on the amplitude:
        eps= 1.d-16
c
c Normalization to the Fresnel scale
        fr= dsqrt( (wvlngth*dist)/2.d+00 )

        x1fr= x1/fr
        x2fr= x2/fr
        x3fr= x3/fr
        x4fr= x4/fr

        xfr=   x/fr

c Calculation of the amplitude of the main band
        xmx1= xfr-x1fr
        xmx2= xfr-x2fr
        call FRESNEL (xmx1,eps,c1,s1)
        call FRESNEL (xmx2,eps,c2,s2)
        cc= c1-c2
        ss= s1-s2
        r_ampli= - (cc+ss)*(opa_ampli/2.d+00)
        i_ampli=   (cc-ss)*(opa_ampli/2.d+00)

c addition of the amplitude of the second band
        xmx3= xfr-x3fr
        xmx4= xfr-x4fr
        call FRESNEL (xmx3,eps,c3,s3)
        call FRESNEL (xmx4,eps,c4,s4)
        cc= c3-c4
        ss= s3-s4
        r_ampli= r_ampli - (cc+ss)*(opa_ampli2/2.d+00)
        i_ampli= i_ampli + (cc-ss)*(opa_ampli2/2.d+00)

        amplir= 1.d+00 + r_ampli 
        amplii= i_ampli
        flux= amplir*amplir + amplii*amplii

        return
        end
c                              end of bar2
C **********************************************************************

      SUBROUTINE FRESNEL(X,EPS,C,S)
C
C SOUS PROGRAMME CALCULANT LES FONCTIONS DE FRESNEL. LE RESULTAT POUR LE 
C COSINUS ET LE SINUS SONT RESPECTIVEMENT (C + OU - EPS), ET (S + OU - EPS).
C
      IMPLICIT REAL*8 (A-H,O-Z)
      N=1
      PI= 4.D+00*DATAN(1.D+00)
      A=PI*X*X
C EVITE UN OVERFLOW SUR K PLUS LOIN 
      A= DMIN1(A, (2.d+00)**(30.))
      A2=A/2.D+00
C SI X<SQRT(2/PI), CALCUL DIRECT PAR SERIE ENTIERE
      IF(A2.LT.1.D+00)GOTO7
      COSA=DCOS(A2)
      SINA=DSIN(A2)
      I=1
      K=A
      SS=1.D+00/(PI*X)
      SC=0.
      U=SS
C BOUCLE D'INTEGRATION PAR PARTIES (TANT QUE LE RESTE DECROIT)
 1    U=-U*DFLOAT(I)/A
      SC=SC+U
      IF(DABS(U).LT.EPS)GOTO5
      I=I+2
      IF(I.GT.K)GOTO4
      U=U*DFLOAT(I)/A
      SS=SS+U
      IF(DABS(U).LT.EPS)GOTO5
      I=I+2
      IF(I.LE.K)GOTO1
C CAS DU RESTE NON NEGLIGEABLE (X "PETIT")
      V=-DFLOAT(I)
      C=SC*COSA+SS*SINA
      S=SC*SINA-SS*COSA
      SC=0.
      SS=1.D+00
      K=2-I     
      UPS=EPS/DABS(U)
      GOTO2
  4    V=DFLOAT(I)
      C=SC*COSA+SS*SINA
      S=SC*SINA-SS*COSA
      K=2-I
      SC=-1.D+00
      SS=0.
      UPS=EPS/DABS(U)
      GOTO3
  7   V=X
      K=3
      C=0.
      S=0.
      SC=X
      SS=0.
      U=1.D+00
      UPS=EPS
C BOUCLE D'INTEGRATION DU RESTE PAR SERIE ENTIERE
  3   V=V*A2/DFLOAT(N)
      W=V/DFLOAT(K)
      SS=SS+W
      IF(DABS(W).LT.UPS)GOTO6
      N=N+1
      K=K+2
   2  V=-V*A2/DFLOAT(N)
      W=V/DFLOAT(K)
      SC=SC+W
      IF(DABS(W).LT.UPS)GOTO6
      N=N+1
      K=K+2
      GOTO3
   6  C=C+SC*U
      S=S+SS*U
      RETURN   
CAS DU RESTE NEGLIGEABLE (X "GRAND")
  5   C=0.5+SC*COSA+SS*SINA
      S=0.5+SC*SINA-SS*COSA
      IF(X.GT.0)RETURN
      C=C-1
      S=S-1
      RETURN
      END

C                        FIN DE FRESNEL
C *******************************************************************

        subroutine etoile (imod,re,x1,x2,x3,x4,opa_ampli,opa_ampli2,wvlngth,dist,x, flux)
c
c        Sous-programme de convolution du signal diffracte par l'etoile,
c        supposee circulaire et uniforme. 
c
        implicit real*8 (a-h,o-z)
c
c npt: echantillonnage sur le rayon stellaire pour le lissage

        npt= 12                                ! doit etre pair pour "somme" ! NB npt sert a la fois pour explorer horizontalement
                                               ! et verticalement le disque stellaire
        tranche= re/dfloat(npt)
        flux= 0.d+00
        som = 0.d+00

        if (re.eq.0) then                      ! etoile ponctuelle
          if (imod.ne.2) call bar (x1,x2,      opa_ampli,           wvlngth,dist,x, flux)
          if (imod.eq.2) call bar2(x1,x2,x3,x4,opa_ampli,opa_ampli2,wvlngth,dist,x, flux)
          return
        endif

        do  i= -npt, npt                       ! etoile taille finie
          p= dfloat(i)*tranche

**********************************************************************
! etoile uniforme:
          coeff= dsqrt( dabs(re*re - p*p) )    ! dabs: si argument tres legerement negatif
!
**********************************************************************
! etoile assombrie centre-bord:
! pour chaque valeur de p on integre ("somme") l'intensite lumineuse de 0 a dsqrt(re**2 - p**2) --->
! prend ~ 4 fois plus de temps que le disque uniforme. On peut aller plus vite en calculant une fois pour toutes
! "coeff" pour differente valeur de p.
!
!                p_norm= p/re                                        ! cause de l'appel de "somme" ---> a ne considerer que
!                fac= 1.d0 - p_norm**2                               ! dans un 2eme temps
!                ymax  = dsqrt( dabs(1.d0 - p_norm**2) )             ! dabs: si argument tres legerement negatif
!                coeff = somme(0.d0,ymax,p_norm,npt)                 ! integre l'intensite verticalement a p_norm=cste
**********************************************************************

          xx= x + p
          if (imod.ne.2) call bar (x1,x2,      opa_ampli,           wvlngth,dist,xx, fluxi)
          if (imod.eq.2) call bar2(x1,x2,x3,x4,opa_ampli,opa_ampli2,wvlngth,dist,xx, fluxi)
          flux= flux + coeff*fluxi
          som= som + coeff
        enddo

        flux= flux/som

        return
        end


c                        FIN DE ETOILE
c ***********************************************************************
c
c
        subroutine instrument (nmax, t,   flux ,npt ,
     *                         trep,rep  ,nrep,
     *                         tl,  fluxl,nptl)
c
c Sous-programme de convolution du signal par la reponse instrumentale.
c 
c t, flux,npt: tableau de npt valeurs, t= temps et flux= flux original correspondant.
c
c trep,rep,nrep: tableau de nrep valeurs, trep= temps et rep= reponse instrumentale correspondante.
c
c tl,fluxl,nptl: tableau de nptl valeurs, tl= nouveau temps et fluxl= flux lisse correspondant.
c
        implicit real*8 (a-h,o-z)
        dimension t(-nmax:nmax), flux(-nmax:nmax)
        dimension trep(nmax),   rep(nmax)
        dimension tl(nmax),    fluxl(nmax)

c
c som= somme pour construire le flux lisse
c aire= integrale de la fonction instrumentale
c dtmax= longueur en temps de la fonction instrumentale.
c
        nptl= 0
        dtmax= trep(nrep) - trep(1)

        if (dtmax.eq.0) then                ! cas ou la reponse est un delta en zero 
          do i= -npt, npt
            nptl= nptl + 1
            tl(nptl)   = t(i)
            fluxl(nptl)= flux(i)
          enddo
          goto 102
        endif
c
c On boucle d'abord sur le flux original
        do i= -npt, npt
          tmin= t(i) - dtmax
          if( tmin.lt.t(-npt) ) go to 100
          som= 0.
          aire= 0.

c On boucle ensuite sur la reponse instrumentale
          do j= i+1-npt, i+1+npt
            dt= t(i) - t(i-j+1)
            discri= (dt-trep(1))*(dt-trep(nrep))
            if(discri.le.(0.d0)) then
c On cherche le point de la rep. instru. qui est le plus proche du point
c courant (i-j+1) et on interpole.
              coeff= 0.
              do k= 1, nrep-1
                discri= ( dt-trep(k) )*( dt-trep(k+1) )
                if ( discri.le.0.d+00 ) then
                  coeff= rep(k+1) - rep(k)
                  coeff= (dt-trep(k))/(trep(k+1)-trep(k))*coeff
                  coeff= coeff + rep(k)
                  go to 101
                endif
              enddo
101           som= som + coeff*flux(i-j+1)
              aire= aire + coeff
            endif
          enddo
          nptl= nptl + 1
          tl(nptl)   = t(i)
          fluxl(nptl)= som/aire
100       continue
        enddo

102     continue
        return
        end

c                        FIN DE INSTRUMENT
c ****************************************************************************
 
        subroutine ombre_geometrique (imod,t0,duree,width,trans,width2,trans2,deltat,vn,phi1,phi0)

        implicit real*8 (a-h,o-z)
        parameter (nt=10)
        dimension t(nt), indx(nt)


        if (imod.eq.0) then
          t(1)= t0-duree/2.d0                ! milieu de la bande centre sur t0
          t(2)= t0-width/(2.d0*vn)
          t(3)= t0+width/(2.d0*vn)
          t(4)= t0+duree/2.d0
         call indexx(4,t,indx)
        endif

        if (imod.eq.-1) then 
          t(1)= t0-duree/2.d0                ! bord droit de la bande centre sur t0
          t(2)= t0
          t(3)= t0+(width/vn)
          t(4)= t0+duree/2.d0
          call indexx(4,t,indx)
        endif

        if (imod.eq.1) then 
          t(1)= t0-duree/2.d0                ! bord gauche de la bande centre sur t0
          t(2)= t0-(width/vn)
          t(3)= t0
          t(4)= t0+duree/2.d0
          call indexx(4,t,indx)
        endif

        if (imod.eq.2) then
          t(1)= t0-duree/2.d0                ! milieu de la bande centre sur t0
          t(2)= t0-width/(2.d0*vn)
          t(3)= t0+width/(2.d0*vn)
          t(4)= t0+duree/2.d0
          t(5)= t0+deltat-width2/(2.d0*vn)
          t(6)= t0+deltat+width2/(2.d0*vn)
          call indexx(6,t,indx)
        endif

        if (imod.ne.2) then                                              ! one band
          write(20,*) sngl(t(indx(1))), sngl(phi1)
          write(20,*) sngl(t(indx(2))), sngl(phi1)                        !SNGL = meme fonction que REAL (transforme un double precision en reel)
          write(20,*) sngl(t(indx(2))), sngl(phi0+(phi1-phi0)*trans)
          write(20,*) sngl(t(indx(3))), sngl(phi0+(phi1-phi0)*trans)
          write(20,*) sngl(t(indx(3))), sngl(phi1)
          write(20,*) sngl(t(indx(4))), sngl(phi1)
        endif

        if (imod.eq.2) then                                               ! two bands
          if (deltat.lt.0) then ; transg= trans2 ; transd= trans  ; endif
          if (deltat.gt.0) then ; transg= trans  ; transd= trans2 ; endif
          write(20,*) sngl(t(indx(1))), sngl(phi1)
          write(20,*) sngl(t(indx(2))), sngl(phi1)
          write(20,*) sngl(t(indx(2))), sngl(phi0+(phi1-phi0)*transg)
          write(20,*) sngl(t(indx(3))), sngl(phi0+(phi1-phi0)*transg)
          write(20,*) sngl(t(indx(3))), sngl(phi1)
          write(20,*) sngl(t(indx(4))), sngl(phi1)
          write(20,*) sngl(t(indx(4))), sngl(phi0+(phi1-phi0)*transd)
          write(20,*) sngl(t(indx(5))), sngl(phi0+(phi1-phi0)*transd)
          write(20,*) sngl(t(indx(5))), sngl(phi1)
          write(20,*) sngl(t(indx(6))), sngl(phi1)
        endif

        close(20)
        return
        end
**************************** fin de ombre_geometrique ***********************


        subroutine interpol (nmax,nmod,tmod,fmod,t, fmod_inter)
        implicit real*8 (a-h,o-z)
        real*8 tmod(nmax), fmod(nmax)

        iflag= 0
        do i= 1, nmod-1
          fac= (t-tmod(i))*(t-tmod(i+1))
          if (fac.le.0) then
            fmod_inter= (fmod(i+1)- fmod(i))/(tmod(i+1)- tmod(i))
            fmod_inter= fmod_inter*(t-tmod(i)) + fmod(i)
            flag= 1
          return
          endif
        enddo 

        if (flag.eq.0) write(*,*) 'No interpolation possible!!!'

        return
        end
*********************************** fin de interpol *********************

c
c Integrale d'une fonction par la methode de Simpson
c Integre la fonction f(x,y) sur y entre a et b, avec npt intervalles
c Attention, npt doit etre pair !
c
c
c        integrale(a,b)= (h/3)*[ f(x,a0) + 4f(x,a1) + 2f(x,a2) + 4f(x,a3) + ...
c                  ... + 2f(x,an-2) + 4f(x,an-1) + f(x,an) ]
c
c        ou h= (b-a)/npt, a0= a et an= b
c
        function somme(a,b,x,npt)
        implicit real*8 (a-h,o-z)
c
        h= (b-a)/dfloat(npt)
c
c somme des termes impairs
c
        somi= 0.d+00
        do i= 1, npt-1, 2
          y= a + dfloat(i)*h
          somi= somi + f(x,y)
        enddo
        somi= 4.d+00*somi
c
c somme des termes pairs
c
        somp= 0.d+00
        do i= 2, npt-2, 2
          y= a + dfloat(i)*h
          somp= somp + f(x,y)
        enddo
        somp= 2.d+00*somp
c
c addition des bornes
c
        somme= f(x,a) + somi + somp + f(x,b)
        somme= (h*somme)/3.d+00

        return
        end
*************************** fin de somme  **************************************
c
c Calcul de l'intensite lumineuse a (x,y) du centre de l'etoile (assombrissement centre-bord)
c
c l'intensite lumineuse I emise par un element de surface de l'etoile
c dont la normale fait un angle theta avec la ligne de visee, avec mu=cos(theta), est de la form:
c
c I(mu)= 1 - sum_1^4 a_k*[1-mu^(k/2)]
c
c NB. l'intensite vaut un au milieu du disque stellaire (ou mu=1) 
c
c Sources: Claret Astron. Astrophys. 363, 1081\D01190 (2000) et
c mail Pierre Kervella, 2 juillet 2008
c
        function f(x,y)
        implicit real*8 (a-h,o-z)
        dimension a(4)
        real*8 mu
        tol= -1.d-15

        a(1)=  0.6699d0                ! mail Pierre Kervella, 2 juillet 2008
        a(2)= -0.7671d0
        a(3)=  1.6405d0
        a(4)= -0.6607d0

        fac= 1.d0 - (x**2 + y**2)
        if (fac.lt.(0.d0)) then
          if (fac.lt.tol) write (*,*) 'error: 1-x^2-y^2=', fac
          fac= 0.d0
        endif

        mu= dsqrt(fac)

        f= 1.d0
        do k= 1, 4
          f= f - a(k)*( 1.d0-mu**(dfloat(k)/2.d0) )
        enddo
c
c NB. on peut avoir mu=0 dans la formule precedente, mais pas de pb pour
c calculer mu**(dfloat(k)/2.d0) (=0 dans ce cas)
c
        return
        end
 
****************************** fin de f ***************************

        subroutine error_bar (nmax,Dchi2)
        implicit real*8 (a-h,o-z)

        open (unit=24,file='fort.24',status='old',form='formatted',position='rewind')
        chi2_min= 1.d50
        do i= 1, nmax
          read (24,*,err=103,end=103) t0, chi2, nfit
          if (chi2.le.chi2_min) then
            chi2_min     = chi2
            t0_chi2_min  = t0
            nfit_chi2_min= nfit
          endif
        enddo
103     continue

        close (24)
        open (unit=24,file='fort.24',status='old',form='formatted',position='rewind')
        
        t_inf=  1.d50
        t_sup= -1.d50
        do i= 1, nmax
          read (24,*,err=104,end=104) t0, chi2, nfit
          if (chi2.le.(chi2_min+dchi2)) then
            if (t0.le.t_inf) t_inf= t0
            if (t0.ge.t_sup) t_sup= t0
          endif
        enddo
104     continue

        write (*,*) 
        write (*,'(A,f10.4,3x,f8.3,3x,I3)') 't0, chi2 and nfit at minimum: ', t0_chi2_min, chi2_min, nfit_chi2_min
        write (*,'(A,f5.2,2x,f9.3,2x,f9.3)')  
     *  'Dchi2, interval or chi2 < chi2_min + dchi2: ', sngl(dchi2), sngl(t_inf), sngl(t_sup)
 
    
        t0_best= (t_inf + t_sup)/2.d0 ; t0_err= (-t_inf + t_sup)/2.d0
        write (*,*) 'Result: t0 best=', t0_best, '+/-', sngl(t0_err)

        ihh= int(t0_best/3600.d0)
        amn= (t0_best/3600.d0 - dfloat(ihh))*60.
        imn= int(amn)
        sec= (amn - dfloat(imn))*60.
        write (*,'(A,I2.2,1x,I2.2,1x,f6.3)') 'or:        ', ihh, imn, sec

        close (24)

        return
        end

*********************************  fin de error_bar *****************************************
c
c Num. Rec. prend n valeurs du tableau arr, et en fait une liste croissante
c avec indx: arr(indx(1)),... arr(indx(n))
c
      SUBROUTINE indexx(n,arr,indx)
      INTEGER n,indx(n),M,NSTACK
c     REAL arr(n)
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
c     REAL a
      REAL*8 a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) then
          stop 
          write(*,*) 'NSTACK too small in indexx'
        endif 
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
************************************ fin de indexx *************************

