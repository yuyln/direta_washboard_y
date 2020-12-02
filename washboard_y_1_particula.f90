module variaveis
    implicit none
    character*255 corrente, formt
    integer:: t, i_s, k, i__, j, i_, entx, enty, cx, cy, cont, nx, ny
    integer, parameter:: ni = 5000000, ns = 1, np = 7, nc = 3000, &
&                        corte = 1, relax = 100000
    doubleprecision, parameter:: fator = 1.d-6
    doubleprecision, parameter:: tf = 500.d0, hs = tf / ni, &
    &                            pi = dacos(-1.d0), lx = 30.d0, &
    &                            ly = 30.d0, xhi = 1.d0, &
    &                            temp = 0.d0, fcf = 3.d0, &
    &                            ac = 0.d0, &
    &                            u = 3.d0, rp = 0.65d0, &
    &                            hc = fcf / nc
    doubleprecision, parameter:: nac = 100.d0, omega = dble(nac) / dble(ni)
    doubleprecision, parameter:: betadamp = 5.d0, &
&                                damping = dsqrt(1.d0 / (1.d0 + betadamp ** 2.d0)), &
&                                beta = betadamp * damping
    doubleprecision, dimension(1: ns):: x, x0, y, y0
    doubleprecision, dimension(1: np):: x_pin, y_pin
    doubleprecision:: pos_(1: 2, 1: ns, 0: ni / corte)
    doubleprecision:: ssx, ssy, vx, vy, inicio, fim, fc, vmx, vmy
    doubleprecision:: k1x, k1y, k2x, k2y, fx, fy, na, p1, p2, x_temp, y_temp
    doubleprecision:: fcx, fcy, aux, hx, hy, wy
    logical:: esq, dir, cima, baixo
end module

module gera_function
    implicit none
    integer:: i, index_
    integer, parameter:: n = 1000000
    doubleprecision:: xg, hw
    doubleprecision, dimension(0: n):: y_bessel, y_wash
    doubleprecision, parameter:: aw = 0.d0, ab = 0.5d0
    doubleprecision::bb, bw, hb
    
    contains
    
    subroutine g_bessel(xhi)
        implicit none
        doubleprecision:: BESSK, xhi
        do i = 0, n
            xg = ab + i * hb
            y_bessel(i) = BESSK(1, xg / xhi)
        enddo
    end subroutine g_bessel
    
    subroutine g_wash(lx, np, u)
        implicit none
        integer:: np
        doubleprecision:: lx, aux, u
        doubleprecision, parameter:: pi = dacos(-1.d0)
        bw = lx
        hw = (bw - aw) / n
        aux = 2.d0 * pi * np / lx
        do i = 0, n
            xg = i * hw
            y_wash(i) = aux * dsin(aux * xg) * u
        enddo
    end subroutine
    
    doubleprecision function tab_bessel(x)
        implicit none
        doubleprecision:: x
        if (x > bb) then
            tab_bessel = 0.d0
        elseif (x < ab) then
            tab_bessel = y_bessel(0)
        else
            index_ = nint(x / hb)
            if (index_ > n) then
                index_ = n
            endif
            if (index_ < 0) then
                index_ = 0
            endif
            tab_bessel = y_bessel(index_)
        endif
    end function
    
    doubleprecision function tab_wash(x)
        implicit none
        integer:: index_
        doubleprecision:: x
        if (x > bw) then
            tab_wash = 0.d0
        elseif (x < aw) then
            tab_wash = y_wash(0)
        else
            index_ = nint(x / hw)
            if (index_ > n) then
                index_ = n
            endif
            if (index_ < 0) then
                index_ = 0
            endif
            tab_wash = y_wash(index_)
        endif
    end function    
end module

subroutine cortes
    use variaveis
    use gera_function
    doubleprecision:: BESSK
    do i = 0, n
        aux = dble(i) * 30.d0 / dble(n)
        if (BESSK(1, aux / xhi) <= fator) then
            bb = aux
            hb = (bb - ab) / n
            exit
        endif
    enddo
end subroutine

subroutine ler_skyr(ns, x0, y0, lx, ly)
    implicit none
    integer:: ns, i_s
    doubleprecision, dimension(1: ns):: x0, y0
    doubleprecision:: lx, ly, xr, yr
    open(1, file='skyrmion.dat', status='old')
    do i_s = 1, ns
        read(1, *) xr, yr
        if ((xr <= lx).and.(yr <= ly)) then
            x0(i_s) = xr
            y0(i_s) = yr
        endif
    enddo
    close(1)
end subroutine

subroutine skyr_wash(x0, tab_wash, wx)
    implicit none
    doubleprecision:: x0, wx
    doubleprecision, external:: tab_wash
    wx = tab_wash(x0)
end subroutine

subroutine skyr_skyr(x0, y0, x1, y1, ssx, ssy, tab_bessel, xhi)
    implicit none
    doubleprecision:: x0, y0, x1, y1
    doubleprecision:: ssx, ssy, xhi, dx, dy, d, d2
    doubleprecision, external:: tab_bessel
    call distancia(x0, y0, x1, y1, dx, dy, d, d2)
    ssx = ssx + tab_bessel(d) * dx / (d * xhi)
    ssy = ssy + tab_bessel(d) * dy / (d * xhi)
end subroutine

subroutine contorno(x0, y0, lx, ly)
    doubleprecision:: x0, y0
    doubleprecision:: lx, ly
    
    if (x0 > lx) then
        x0 = x0 - lx
    elseif (x0 < 0.d0) then
        x0 = lx + x0
    endif
    
    if (y0 > ly) then
        y0 = y0 - ly
    elseif (y0 < 0.d0) then
        y0 = ly + y0
    endif
end subroutine

subroutine forcas(desx, desy)
    use variaveis
    use gera_function
    doubleprecision:: desx, desy
    ssx = 0.d0
    ssy = 0.d0
    do i__ = 1, ns
        if (i__ /= i_s) then
            call skyr_skyr(x0(i_s) + desx, y0(i_s) + desy, x0(i__), y0(i__), ssx, ssy, tab_bessel, xhi)
        endif
    enddo
    call skyr_wash(y0(i_s) + desy, tab_wash, wy)
    call random_number(na)
    x_temp = (2.d0 * na - 1.d0) * dsqrt(temp)
    call random_number(na)
    y_temp = (2.d0 * na - 1.d0) * dsqrt(temp)
    
    fcx = fc * dcos(ac)
    fcy = fc * dsin(ac)
    p1 = ssx + x_temp + fcx
    p2 = ssy + y_temp + fcy + wy
    fx = (damping * p1 + beta * p2) / (beta ** 2.d0 + damping ** 2.d0)
    fy = (damping * p2 - beta * p1) / (beta ** 2.d0 + damping ** 2.d0)

end

subroutine integraRK2
    use variaveis
    use gera_function
    
    call forcas(0.d0, 0.d0)
    k1x = fx * hs
    k1y = fy * hs
        
    call forcas(k1x, k1y)
    k2x = fx * hs
    k2y = fy * hs
    
    fx = (k1x + k2x) / 2.d0
    fy = (k1y + k2y) / 2.d0
end subroutine

subroutine mov
    use variaveis
    use gera_function
    call integraRK2
    x(i_s) = x0(i_s) + fx
    y(i_s) = y0(i_s) + fy
    vx = (x(i_s) - x0(i_s)) / hs
    vy = (y(i_s) - y0(i_s)) / hs
end subroutine

subroutine reset(ns, x, y, x0, y0)
    integer:: ns, i_s
    doubleprecision, dimension(1: ns):: x0, y0, x, y
    do i_s = 1, ns
        x0(i_s) = x(i_s)
        y0(i_s) = y(i_s)
    enddo
end subroutine

subroutine distancia(x0, y0, x1, y1, dx, dy, d, d2)
    doubleprecision:: x0, y0, x1, y1, d, d2, dx, dy
    dx = x0 - x1
    dy = y0 - y1
    d = dsqrt(dx ** 2.d0 + dy ** 2.d0)
    d2 = d ** 2.d0
end subroutine

program principal
    use variaveis
    use gera_function
    call cortes
    call ler_skyr(ns, x0, y0, lx, ly)
    call g_bessel(xhi)
    call g_wash(ly, np, u)
    open(11, file='dados.dat', status='unknown')
    write(11, "(a20, 2f10.5)") "Range Bessel: ", ab, bb
    write(11, "(a20, 2f10.5)") "LX, LY: ", lx, ly
    write(11, "(a20, 2i10)") "NP: ", np
    write(11, "(a20, 2i10)") "NS: ", ns
    write(11, "(a20, i20, f10.5)") "NI, HS: ", ni, hs
    write(11, "(a20, i20, f10.5)") "NC, HC: ", nc, hc
    write(11, "(a20, 3f10.5)") "B/D, B, D: ", betadamp, beta, damping
    close(11)
    call system("mkdir pos")
    write(formt, 24) 2 * ns
    24 format("(", I10, "(f30.15))")
    open(21, file='theta_mov.dat', status='unknown')
    fc = 0.d0
    do t = 1, relax
        do i_s = 1, ns
            call mov
            call contorno(x(i_s), y(i_s), lx, ly)
        enddo
        call reset(ns, x, y, x0, y0)
    enddo
    
    open(22, file='inicial.dat', status='unknown')
    write(22, *) (x(i_s), y(i_s), i_s = 1, ns)
    close(22)
    
    do k = 1, nc
        fc = dble(k) * hc
        write(corrente, 11) fc
11      format("pos/pos_", f7.3, ".dat")
        open(11, file=corrente, status='unknown')
        vmx = 0.d0
        vmy = 0.d0
        do t = 1, ni
            do i_s = 1, ns
                call mov
                call contorno(x(i_s), y(i_s), lx, ly)
                vmx = vmx + vx / (ns * ni)
                vmy = vmy + vy / (ns * ni)
            enddo
            if (mod(t, corte) == 0) then
                    pos_(1, i_s, t / corte) = x(i_s)
                    pos_(2, i_s, t / corte) = y(i_s)
                end if
            call reset(ns, x, y, x0, y0)
        enddo
        write(*, 10) k * 100.d0 / nc, fc
10      format(f6.2, "%", f7.3)
        write(21, *) fc, datan(vmy / vmx) * 180.d0 / pi, vmx, vmy, vmy / vmx
        write(11, formt) pos_
        close(11)
    enddo
    close(21)
    call cpu_time(fim)
    print*, fim - inicio
end program


!função bessel (copiada)------------------------------------------------
FUNCTION BESSK(N,X)
      IMPLICIT NONE
      INTEGER N,J
      DOUBLE PRECISION X,BESSK,BESSK0,BESSK1,TOX,BK,BKM,BKP
!     ------------------------------------------------------------------------
!     CE SOUS-PROGRAMME CALCULE LA FONCTION BESSEL MODIFIFIEE 3E ESPECE
!     D'ORDRE N ENTIER POUR TOUT X REEL POSITIF > 0.  ON UTILISE ICI LA
!     FORMULE DE RECURRENCE CLASSIQUE EN PARTANT DE BESSK0 ET BESSK1.
!
!     THIS ROUTINE CALCULATES THE MODIFIED BESSEL FUNCTION OF THE THIRD
!     KIND OF INTEGER ORDER, N FOR ANY POSITIVE REAL ARGUMENT, X. THE
!     CLASSICAL RECURSION FORMULA IS USED, STARTING FROM BESSK0 AND BESSK1.
!     ------------------------------------------------------------------------ 
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
!     ------------------------------------------------------------------------
      IF (N.EQ.0) THEN
      BESSK = BESSK0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSK = BESSK1(X)
      RETURN
      ENDIF
      IF (X.EQ.0.D0) THEN
      BESSK = 1.D30
      RETURN
      ENDIF
      TOX = 2.D0/X
      BK  = BESSK1(X)
      BKM = BESSK0(X)
      DO 11 J=1,N-1
      BKP = BKM+DFLOAT(J)*TOX*BK
      BKM = BK
      BK  = BKP
   11 CONTINUE
      BESSK = BK
      RETURN
      END
!     ----------------------------------------------------------------------
      FUNCTION BESSK0(X)
!     CALCUL DE LA FONCTION BESSEL MODIFIEE DU 3EME ESPECE D'ORDRE 0
!     POUR TOUT X REEL NON NUL.
!
!     CALCULATES THE THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF 
!     ORDER ZERO FOR ANY POSITIVE REAL ARGUMENT, X.
!     ----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION X,BESSK0,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,&   
     &Q5,Q6,Q7,BESSI0
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,&
      &0.23069756D0,0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,&
      &0.2189568D-1,-0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/

      IF(X.EQ.0.D0) THEN
      BESSK0=1.D30
      RETURN
      ENDIF
      IF(X.LE.2.D0) THEN
      Y=X*X/4.D0
      AX=-DLOG(X/2.D0)*BESSI0(X)
      BESSK0=AX+(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      Y=(2.D0/X)
      AX=DEXP(-X)/DSQRT(X)
      BESSK0=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END
!     ----------------------------------------------------------------------
      FUNCTION BESSK1(X)
!     CALCUL DE LA FONCTION BESSEL MODIFIEE DE 3EME ESPECE D'ORDRE 1
!     POUR TOUT X REEL POSITF NON NUL.
!
!     CALCULATES THE THE MODIFIED BESSEL FUNCTION OF THE THIRD KIND OF 
!     ORDER ONE FOR ANY POSITIVE REAL ARGUMENT, X.
!     ----------------------------------------------------------------------
      IMPLICIT NONE
      DOUBLE PRECISION X,BESSK1,Y,AX,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,&
     &Q5,Q6,Q7,BESSI1

      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,0.15443144D0,-0.67278579D0,&  
     &-0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,&
      &-0.3655620D-1,0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/

      IF(X.EQ.0.D0) THEN
      BESSK1=1.D32
      RETURN
      ENDIF
      IF(X.LE.2.D0) THEN
      Y=X*X/4.0D0
      AX=DLOG(X/2.0D0)*BESSI1(X)
      BESSK1=AX+(1.D0/X)*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      Y=(2.D0/X)
      AX=DEXP(-X)/DSQRT(X)
      BESSK1=AX*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END
!
!     Bessel Function of the 1st kind of order zero.
!
      FUNCTION BESSI0(X)
      IMPLICIT NONE
      DOUBLE PRECISION X,BESSI0,Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,&
     &Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,&
     &1.2067429D0,& 
     &0.2659732D0,0.360768D-1,0.45813D-2/

      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,& 
     &0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,&  
     &0.2635537D-1,-0.1647633D-1,0.392377D-2/

      IF(DABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
      AX=DABS(X)
      Y=3.75D0/AX
      BX=DEXP(AX)/DSQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI0=AX*BX
      ENDIF
      RETURN
      END
!
!     Bessel Function of the 1st kind of order one.
!
      FUNCTION BESSI1(X)
      IMPLICIT NONE
      DOUBLE PRECISION X,BESSI1,Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,&
      &Q4,Q5,Q6,Q7,Q8,Q9,AX,BX

      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,& 
     &0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/

      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,& 
     &-0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,&        
     &-0.2895312D-1,0.1787654D-1,-0.420059D-2/

      IF(DABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      AX=DABS(X)
      Y=3.75D0/AX
      BX=DEXP(AX)/DSQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI1=AX*BX
      ENDIF
      RETURN
      END
!função bessel (copiada)------------------------------------------------
