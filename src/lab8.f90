subroutine LabNumber8
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB8_N") /=0 ) read (str,*) tmp_str, N
        if (index(str,"LAB8_ni") /=0 ) read (str,*) tmp_str, ni
    end do
    close (2000)
end subroutine LabNumber8



real function Integ1(x, pi,  N)
    implicit none
    integer N
    real:: x, pi
    Integ1=pi
    Integ1=(x**(-0.5+0.02*N))*sin((-0.5+0.02*N+1)*log(x))/(log(x))
end function Integ1

real function Integ2(x, pi, N)
    implicit none
    integer N
    real:: x, pi
    Integ2=(x**(1.0/3))/(1+2*x*cos(pi/6+0.02*N)+x**2)
end function Integ2

real function Integ3(x, pi, N)
    implicit none
    integer N
    real:: x, pi
    Integ3=pi
    Integ3=log(1.0/x)*exp(-(3+0.02*N)*x)*sin((4-0.02*N)*x)
end function Integ3



real function integral1 (F, a, b, ni, pi, N)
    implicit none
    real F
    real:: a, b, h, pi 
    integer i, ni, N
    h=(b-a)/ni
    integral1 = F(a, pi, N)+F(b, pi, N)
    do i=1, ni-1
        if(mod(i, 2)==0) then
            integral1=integral1+2*F(a+i*h,pi, N)
        else 
            integral1=integral1+4*F(a+i*h, pi, N)
        end if
    end do
    integral1=integral1*h/3
end function integral1



subroutine rezFile8
    use module_main
    implicit none
    real integral1
    real, external :: Integ1, Integ2, Integ3
    real(8) :: Vichisl1, Vichisl2, Vichisl3
    real(8) :: TH1, TH2, TH3
 
    TH1 = 0.7842735633765
    TH2 = 1.281999424981667
    TH3 = 0.2307158998792302

    eps = 0.0001
    a = (eps/3.0*(-0.5+0.02*N+1))**(1.0/(-0.5+0.02*N+1))
    b = (1-eps/3.0)**(1.0/(-0.5+0.02*N+1))
    write(*,*)a, b,(a-b)/ni
    Vichisl1 = integral1(Integ1, a, b, ni, pi, N)

    ni=350000000
    a=0
    b=(3.0/eps)**(3.0/2)
    write(*,*) a, b, (a-b)/ni
    Vichisl2 = integral1(Integ2, a, b, ni, pi, N)

    ni=1000000
    a=sqrt(2*eps/(3*(4-0.02)*N))
    b=log(3.0/(eps*(3+0.02*N)**2))/(3+0.02*N)
    write(*,*) a, b, (a-b)/ni
    Vichisl3 = integral1(Integ3, a, b, ni, pi, N)

    write(*,*) TH1, Vichisl1, abs(TH1 - Vichisl1),abs(TH1 - Vichisl1)/TH1 
    write(*,*) TH2, Vichisl2, abs(TH2 - Vichisl2),abs(TH2 - Vichisl2)/TH2
    write(*,*) TH3, Vichisl3, abs(TH3 - Vichisl3),abs(TH3 - Vichisl3)/TH3
  


    


end subroutine rezFile8