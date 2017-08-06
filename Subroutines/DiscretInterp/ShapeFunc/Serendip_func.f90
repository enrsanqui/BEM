SUBROUTINE Serendip_func(Ni,xsi,eta,ldim,nodes,inci)
  !--------------------------------------------------------
  ! Computes Serendipity shape functions Ni(xsi,eta)
  ! for one and two-dimensional (linear/parabolic) finite
  ! boundary elements
  !--------------------------------------------------------
  REAL,INTENT(OUT)           :: Ni(:)             ! Array with shape function
  REAL,INTENT(IN)            :: xsi,eta           ! Intrinsic coordinates
  INTEGER,INTENT(IN)         :: ldim              ! Element dimension
  INTEGER,INTENT(IN)         :: nodes             ! number of nodes
  INTEGER, INTENT(IN)        :: inci(:)           ! Element incidences
  REAL                       :: mxs,pxs,met,pet   ! temporary variables
  SELECT CASE (ldim)
  CASE(1) ! One dimensional element
     Ni(1) = 0.5*(1.0-xsi); Ni(2) = 0.5*(1.0+xsi);
     IF (nodes == 2) RETURN !linear element finished
     Ni(3) = 1.0 -xsi*xsi
     Ni(1) = Ni(1) - 0.5*Ni(3); Ni(2) = Ni(2)-0.5*Ni(3)
  CASE(2) !two dimensional element
     mxs = 1.0-xsi; pxs = 1.0+xsi; met = 1.0-eta; pet = 1.0+eta
     Ni(1) = 0.25*mxs*pet; Ni(2) = 0.25*pxs*met
     Ni(3) = 0.25*pxs*pet; Ni(4) = 0.25*mxs*pet
     IF (nodes == 4) RETURN ! Linear element finished
     IF (Inci(5) > 0) THEN ! zero node = node missing
        Ni(5) = 0.5*(1.0-xsi*xsi)*met; Ni(1) = Ni(1) -0.5*Ni(5);
        Ni(2) = Ni(2)-0.5*Ni(5)
     END IF
     IF (Inci(6) > 0) THEN
        Ni(6) = 0.5*(1.0 - eta*eta)*pxs
        Ni(2) = Ni(2) - 0.5*Ni(6) ; Ni(3) = Ni(3) - 0.5*Ni(6)
     END IF
     IF (Inci(7) > 0) THEN
        Ni(7) = 0.5*(1.0 - xsi*xsi)*pet
        Ni(3) = Ni(3)-0.5*Ni(7) ; Ni(4) = Ni(4) - 0.5*Ni(7)
     END IF
     IF (Inci(8) > 0) THEN
        Ni(8) = 0.5*(1.0 - eta*eta)*mxs
        Ni(4) = Ni(4) - 0.5*Ni(8) ; Ni(1) = Ni(1) - 0.5*Ni(8)
     END IF
  CASE DEFAULT ! Error message
     CALL Error_message('Element dimension not 1 or 2')
  END SELECT
  RETURN
  END SUBROUTINE Serendip_func
