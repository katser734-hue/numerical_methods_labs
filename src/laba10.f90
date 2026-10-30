subroutine LabNumber10
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB10_N") /=0 ) read (str,*) tmp_str, N
    end do
    close (2000)
end subroutine LabNumber10

real function  Funct1(x, N)
    real x
    integer N
    Funct1 = exp(2*x) + 3*x- N
end function  Funct1

real function  Funct2(x, N)
    real x, a
    integer N
    a=1.01-0.001*N
    Funct2 = x**4-(1+a)*x**3-(4-a)*x**2+4*(1+a)*x-4*a-sqrt(x-1)
end function  Funct2

real function Dixotomia(x1, x2, N, Funct, p)
    integer ::  N
    real  :: x1, x2, x, Funct, p
    x=(x1+x2)/2
  
    !write(*,*) x, 0
    do while(abs(x1 - x2) > p)
        if(Funct(x1, N)*Funct(x, N) < 0) then
            x2=x
        else
            x1=x
        end if
        x=(x1+x2)/2
        !write(*,*) x, 0
     
    end do
    Dixotomia=x
 
end function Dixotomia

real function Sekus(x1, x2, N, Funct, p)
    integer :: N
    real ::  f, f1, f2,  x1,  x2, x, t, Funct, p
    f1=Funct(x1, N)
    f2=Funct(x2, N)
    !write(*,*) x1,f1
    !write(*,*) x2,f2
    !write(*,*)
    t=1
    do while(abs(t) > p)
        f1=Funct(x1, N)
        f2=Funct(x2, N)
        
        x=x1-f1/(f1-f2)*(x1-x2)
        f=Funct(x, N)
        t=x-x1
        if(f2*f<0) then
            x1=x
        else
            x2=x1
        end if
        f1=Funct(x1, N)
        f2=Funct(x2, N)
        !write(*,*) x1,f1
        !write(*,*) x2,f2
    end do
    Sekus = x
end function Sekus

real function Prosto ( N )
    integer N
    real :: x2, x1
    x1 = 1.0
    !write(*,*) x1, log(N-3*x1)/2.0
    x2=log(N-3*x1)/2.0
    !write(*,*) x2, x2
    !write(*,*) x2, log(N-3*x2)/2.0
    do while (abs(x2-x1) > 0.001)
        x1=x2
        x2 = log(N-3*x1)/2.0
        !write(*,*) x2, x2
        !write(*,*) x2, log(N-3*x2)/2.0
    end do
    Prosto = x2
end function Prosto 

real (8) function Nuton ( N)
    implicit none
    integer N
    real (8) :: x2, x1
    x1 = 1.0
    x2=x1-(exp(2*x1) + 3*x1 - N)/(2*exp(2*x1) + 3)
    !write(*,*) x2
    do while (abs(x2-x1) > 0.001)
        x1=x2
        x2=x1-(exp(2*x1) + 3*x1- N)/(2*exp(2*x1) + 3)
        !write(*,*) x2
    end do
    Nuton = x2
end function Nuton


real function Prosto2 ( N )
    implicit none
    integer N
    real :: x2, x1, a
    x1 =3.0
    a = 1.01-0.001*N
    !write(*,*) x1, ((1+a)*x1**3+(4-a)*x1**2-4*(1+a)*x1+4*a+sqrt(x1-1))**(1/4.0)
    x2=((1+a)*x1**3+(4-a)*x1**2-4*(1+a)*x1+4*a+sqrt(x1-1))**(1/4.0)


    !write(*,*) x2, x2
    !write(*,*) x2, ((1+a)*x2**3+(4-a)*x2**2-4*(1+a)*x2+4*a+sqrt(x2-1))**(1/4.0)

    do while (abs(x2-x1) > 0.00001) !1d-5
        x1=x2
        x2=((1+a)*x1**3+(4-a)*x1**2-4*(1+a)*x1+4*a+sqrt(x1-1))**(1/4.0)

       !write(*,*) x2, x2
       !write(*,*) x2, ((1+a)*x2**3+(4-a)*x2**2-4*(1+a)*x2+4*a+sqrt(x2-1))**(1/4.0)
    end do
    Prosto2 = x2
end function Prosto2 

real  function Nuton2 ( N)
    integer N
    real  :: x2, x1
    x1 = 3.0
    a=1.01-0.001*N
    x2=x1-(x1**4-(1+a)*x1**3-(4-a)*x1**2+4*(1+a)*x1-4*a-sqrt(x1-1))/&
    (4*x1**3-3*(1+a)*x1**2-2*(4-a)*x1+4*(1+a)-1.0/(2*sqrt(x1-1)))
    !write(*,*) x2
    do while (abs(x2-x1) > 0.0001)
        x1=x2
        x2=x1-(x1**4-(1+a)*x1**3-(4-a)*x1**2+4*(1+a)*x1-4*a-sqrt(x1-1))/&
        (4*x1**3-3*(1+a)*x1**2-2*(4-a)*x1+4*(1+a)-1.0/(2*sqrt(x1-1)))
        !write(*,*) x2
    end do
    Nuton2 = x2
end function Nuton2

subroutine rezFile10
    use module_main
    implicit none
    real, external :: Funct1, Funct2
    real Dixotomia, Sekus, Prosto
    real Prosto2
    real Dix, Sek, Pro
    real Dix2, Sek2, Pro2
    real(8) :: TCH, Nut, Nuton
    real ::  Nuton2, Nut2
    real(8) :: TCH2
    real  xxx1, xxx2
    real por1, por2
    TCH=0.58743214919d0
    TCH2=2.18579
    por1=0.001 
    por2=0.0001 
    

    xxx2=0; xxx1=1
    Dix = Dixotomia(xxx1, xxx2, N, Funct1, por1)


    xxx1=0; xxx2=1
    Sek = Sekus(xxx1, xxx2, N, Funct1, por1)

    Pro = Prosto ( N )
    Nut = Nuton  ( N )
    write(*,*) "First Formula: Tch, chet, abs, otn"
    write(*,*) "Dix",  TCH, Dix, abs(TCH - Dix), abs(TCH - Dix)/TCH
    write(*,*) "Sek",  TCH, Sek, abs(TCH - Sek), abs(TCH - Sek)/TCH
    write(*,*) "Pro",  TCH, Pro, abs(TCH - Pro), abs(TCH - Pro)/TCH
    write(*,*) "Nut",  TCH, Nut, abs(TCH - Nut), abs(TCH - Nut)/TCH
    write(*,*)
    write(*,*)
    xxx2=3; xxx1=2
    Dix2 = Dixotomia(xxx1, xxx2, N, Funct2, por2)


    xxx1=2; xxx2=3
    Sek2 = Sekus(xxx1, xxx2, N, Funct2, por2)

    Pro2 = Prosto2 ( N )
    Nut2 = Nuton2  ( N )

    write(*,*) "Second Formula: Tch, chet, abs, otn"
    write(*,*) "Dix", TCH2, Dix2, abs(TCH2 - Dix2), abs(TCH2 - Dix2)/TCH2
    write(*,*) "Sek", TCH2, Sek2, abs(TCH2 - Sek2), abs(TCH2 - Sek2)/TCH2
    write(*,*) "Pro", TCH2, Pro2, abs(TCH2 - Pro2), abs(TCH2 - Pro2)/TCH2
    write(*,*) "Nut", TCH2, Nut2, abs(TCH2 - Nut2), abs(TCH2 - Nut2)/TCH2

 
end subroutine rezFile10