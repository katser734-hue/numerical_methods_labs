subroutine LabNumber6
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB6_N") /=0 ) read (str,*) tmp_str, N
        if (index(str,"LAB6_ni") /=0 ) read (str,*) tmp_str, ni
    end do
    close (2000)
end subroutine LabNumber6


real function Func(x0, p)
    implicit none
    real:: x0, p
    Func=(x0**p)/(1+x0**2)+(2*x0**(-p+2))/(1+x0**2)**2/(-p+1)
end function Func

real function Func2(x0, p)
    implicit none
    real:: x0, p
    Func2=(cos(x0))**2*cos(p*x0)
end function Func2

real function Func3(x0, p)
    implicit none
    real:: x0, p
    Func3=x0**p*log(1.0/x0)**2
end function Func3

real function Func4(x0, p)
    implicit none
    real:: x0, p
    Func4=1.0/((p-1)*x0**2-2*x0+1)**(3.0/2)
end function Func4

real function integral12 (F, a, b, n, p)
    implicit none
    real F
    integer n
    real a, b, h, p
    integer :: i 
    integral12 = 0
    h=(b-a)/n
    do i=0, n-1
        integral12=integral12+F(a+h*i, p)
    end do
    integral12=integral12*h
end function integral12




real function integral13 (F, a, b, n, p)
    implicit none
    real F
    integer n
    real:: a, b, h, p
    integer :: i
    integral13 = 0
    h=(b-a)/n
    do i=1, n
        integral13=integral13+F(a+i*h,p)
    end do
    integral13=integral13*h
end function integral13

real function integral11 (F, a, b, n, p)
    implicit none
    real F
    real:: a, b, h, p
    integer  :: i , n
    h=(b-a)/n
    integral11 = F(a,p)/2.0+F(b,p)/2.0
    do i=1, n-1
        integral11=integral11+F(a+i*h,p)
    end do
    integral11=integral11*h
end function integral11



real function integral14 (F, a, b, n, p)
    implicit none
    real F
    integer n
    real:: a, b, h, p
    integer :: i
    integral14 = 0
    h=(b-a)/n
    do i=0, n-1
        integral14=integral14+F(a+h*i+h/2.0,p)
    end do
    integral14=integral14*h
end function integral14

real function integral15 (F, a, b, n, p)
    implicit none
    real F
    real:: a, b, h, p 
    integer i, n
    h=(b-a)/n
    integral15 = F(a, p)+F(b, p)
    do i=1, n-1
        if(mod(i, 2)==0) then
            integral15=integral15+2*F(a+i*h,p)
        else 
            integral15=integral15+4*F(a+i*h, p)
        end if
    end do
    integral15=integral15*h/3
end function integral15

real function integral16 (F, a, b, n, p)
    implicit none
    real F
    real:: a, b, h, p 
    integer i, n
    h=(b-a)/n
    integral16 = F(a, p)+F(b, p)
    do i=1, n-1
        if(mod(i, 3)==0) then
            integral16=integral16+2*F(a+i*h,p)
        else 
            integral16=integral16+3*F(a+i*h, p)
        end if
    end do
    integral16=integral16*h*3.0/8.0
end function integral16


subroutine rezFile6
    use module_main
    implicit none
    real , external:: Func,Func2,Func3,Func4
    real  integral11,integral12,integral13,integral14,integral15,integral16
    real Tochnoe1, Tochnoe2, Tochnoe3, Tochnoe4
    !open (unit = 2016, file = trim(resFileName))
    Tochnoe1 = 38.472232614
    Tochnoe2 = 0.0007388491
    Tochnoe3 = 0.215960
    Tochnoe4 = 1.928280
    p=0.969+N*0.001
    a=0
    b=1
    write (*,*) "N=", N
    write (*,*)  " int1Trapecia:" , 1.0/2/(1-p)+integral11(Func, a, b, ni, p)
    write (*,*)  " int1Left:" , 1.0/2/(1-p)+integral12(Func, a, b, ni, p)
    write (*,*)  " int1Right:" , 1.0/2/(1-p)+integral13(Func, a, b, ni, p)
    write (*,*)  " int1Centr:" , 1.0/2/(1-p)+integral14(Func, a, b, ni, p)
    write (*,*)  " int1Simpson:" , 1.0/2/(1-p)+integral15(Func, a, b, ni, p)
    write (*,*)  " int1Nuton:" , 1.0/2/(1-p)+integral16(Func, a, b, ni+1, p)
    write (*,*)   "Tochnoe1"
    write (*,*)  " DTrapecia:" , abs(1.0/2/(1-p)+integral11(Func, a, b, ni, p)-Tochnoe1)
    write (*,*)  " D1Left:" , abs(1.0/2/(1-p)+integral12(Func, a, b, ni, p)-Tochnoe1)
    write (*,*)  " D1Right:" , abs(1.0/2/(1-p)+integral13(Func, a, b, ni, p)-Tochnoe1)
    write (*,*)  " D1Centr:" , abs(1.0/2/(1-p)+integral14(Func, a, b, ni, p)-Tochnoe1)
    write (*,*)  " D1Simpson:" , abs(1.0/2/(1-p)+integral15(Func, a, b, ni, p)-Tochnoe1)
    write (*,*)  " D1Nuton:" , abs(1.0/2/(1-p)+integral16(Func, a, b, ni+1, p)-Tochnoe1)
    write (*,*)  "Nuton"
    write (*,*)  " DTrapecia:" , abs(1.0/2/(1-p)+integral11(Func, a, b, ni, p)-1.0/2/(1-p)-integral16(Func, a, b, ni+1, p))
    write (*,*)  " D1Left:" , abs(1.0/2/(1-p)+integral12(Func, a, b, ni, p)-1.0/2/(1-p)-integral16(Func, a, b, ni+1, p))
    write (*,*)  " D1Right:" , abs(1.0/2/(1-p)+integral13(Func, a, b, ni, p)-1.0/2/(1-p)-integral16(Func, a, b, ni+1, p))
    write (*,*)  " D1Centr:" , abs(1.0/2/(1-p)+integral14(Func, a, b, ni, p)-1.0/2/(1-p)-integral16(Func, a, b, ni+1, p))
    write (*,*)  " D1Simpson:" , abs(1.0/2/(1-p)+integral15(Func, a, b, ni, p)-1.0/2/(1-p)-integral16(Func, a, b, ni+1, p))
    write (*,*)  "Right:"
    write (*,*)  " DR1Trapecia:" , abs(1.0/2/(1-p)+integral11(Func, a, b, ni, p)-1.0/2/(1-p)-integral13(Func, a, b, ni, p))
    write (*,*)  " DR1Left:" , abs(1.0/2/(1-p)+integral12(Func, a, b, ni, p)-1.0/2/(1-p)-integral13(Func, a, b, ni, p))
    write (*,*)  " DR1Centr:" , abs(1.0/2/(1-p)+integral14(Func, a, b, ni, p)-1.0/2/(1-p)-integral13(Func, a, b, ni, p))
    write (*,*)  " DR1Simpson:" , abs(1.0/2/(1-p)+integral15(Func, a, b, ni, p)-1.0/2/(1-p)-integral13(Func, a, b, ni, p))
    write (*,*)  "Left:"
    write (*,*)  " DTrapecia:" , abs(1.0/2/(1-p)+integral11(Func, a, b, ni, p)-1.0/2/(1-p)-integral12(Func, a, b, ni, p))
    write (*,*)  " D1Centr:" , abs(1.0/2/(1-p)+integral14(Func, a, b, ni, p)-1.0/2/(1-p)-integral12(Func, a, b, ni, p))
    write (*,*)  " D1Simpson:" , abs(1.0/2/(1-p)+integral15(Func, a, b, ni, p)-1.0/2/(1-p)-integral12(Func, a, b, ni, p))
    write (*,*)  "Center:"
    write (*,*)  " DTrapecia:" , abs(1.0/2/(1-p)+integral11(Func, a, b, ni, p)-1.0/2/(1-p)-integral14(Func, a, b, ni, p))
    write (*,*)  " D1Simpson:" , abs(1.0/2/(1-p)+integral15(Func, a, b, ni, p)-1.0/2/(1-p)-integral14(Func, a, b, ni, p))
    write (*,*) "Trapecia:"
    write (*,*)  " DTrapecia:" , abs(1.0/2/(1-p)+integral15(Func, a, b, ni, p)-1.0/2/(1-p)-integral11(Func, a, b, ni, p))
    write (*,*)
    write (*,*)
    write (*,*)
    a=0
    b=pi/2
    p = 10.2 + 0.01*N

    write (*,*)  " int2Trapecia:" , integral11(Func2, a, b, ni, p)
    write (*,*)  " int2Left:" , integral12(Func2, a, b, ni, p)
    write (*,*)  " int2Right:" , integral13(Func2, a, b, ni, p)
    write (*,*)  " int2Centr:" , integral14(Func2, a, b, ni, p)
    write (*,*)  " int2Simpson:" , integral15(Func2, a, b, ni, p)
    write (*,*)  " int2Nuton:" , integral16(Func2, a, b, ni+1, p)
    write (*,*)
    write (*,*)   "Tochnoe2"
    write (*,*)  " DTrapecia:" , abs(integral11(Func2, a, b, ni, p)-Tochnoe2)
    write (*,*)  " D2Left:" , abs(integral12(Func2, a, b, ni, p)-Tochnoe2)
    write (*,*)  " D2Right:" , abs(integral13(Func2, a, b, ni, p)-Tochnoe2)
    write (*,*)  " D2Centr:" , abs(integral14(Func2, a, b, ni, p)-Tochnoe2)
    write (*,*)  " D2Simpson:" , abs(integral15(Func2, a, b, ni, p)-Tochnoe2)
    write (*,*)  " D2Nuton:" , abs(integral16(Func2, a, b, ni+1, p)-Tochnoe2)
    write (*,*)  "Nuton"
    write (*,*)  " D2Trapecia:" , abs(integral11(Func2, a, b, ni, p)-integral16(Func2, a, b, ni+1, p))
    write (*,*)  " D2Left:" , abs(integral12(Func2, a, b, ni, p)-integral16(Func2, a, b, ni+1, p))
    write (*,*)  " D2Right:" , abs(integral13(Func2, a, b, ni, p)-integral16(Func2, a, b, ni+1, p))
    write (*,*)  " D2Centr:" , abs(integral14(Func2, a, b, ni, p)-integral16(Func2, a, b, ni+1, p))
    write (*,*)  " D2Simpson:" , abs(integral15(Func2, a, b, ni, p)-integral16(Func2, a, b, ni+1, p))
    write (*,*)  "Right:"
    write (*,*)  " D2Trapecia:" , abs(+integral11(Func2, a, b, ni, p)-integral13(Func2, a, b, ni, p))
    write (*,*)  " D2Left:" , abs(integral12(Func2, a, b, ni, p)-integral13(Func2, a, b, ni, p))
    write (*,*)  " D2Centr:" , abs(integral14(Func2, a, b, ni, p)-integral13(Func2, a, b, ni, p))
    write (*,*)  " D2Simpson:" , abs(integral15(Func2, a, b, ni, p)-integral13(Func2, a, b, ni, p))
    write (*,*)  "Left:"
    write (*,*)  " D2Trapecia:" , abs(integral11(Func2, a, b, ni, p)-integral12(Func2, a, b, ni, p))
    write (*,*)  " D2Centr:" , abs(integral14(Func2, a, b, ni, p)-integral12(Func2, a, b, ni, p))
    write (*,*)  " D2Simpson:" , abs(integral15(Func2, a, b, ni, p)-integral12(Func2, a, b, ni, p))
    write (*,*)  "Center:"
    write (*,*)  " D2Trapecia:" , abs(integral11(Func2, a, b, ni, p)-integral14(Func2, a, b, ni, p))
    write (*,*)  " D2Simpson:" , abs(integral15(Func2, a, b, ni, p)-integral14(Func2, a, b, ni, p))
    write (*,*) "Trapecia:"
    write (*,*)  " D2Trapecia:" , abs(integral15(Func2, a, b, ni, p)-integral11(Func2, a, b, ni, p))
    write (*,*)
    write (*,*)
    write (*,*)
    
    p=1+0.02*N
    h=0.005
    a=0.000001
    b=1
    ni=int((b-a)/h)

    write (*,*)  " int3Trapecia:" ,Tochnoe3,  integral11(Func3, a, b, ni, p),&
    abs(integral11(Func3, a, b, ni, p)-Tochnoe3),abs(integral11(Func3, a, b, ni, p)-Tochnoe3)/Tochnoe3
    h=0.002
    ni=int((b-a)/h)
    write (*,*)  " int3Simpson:" ,Tochnoe3, integral15(Func3, a, b, ni, p),&
    abs(integral15(Func3, a, b, ni, p)-Tochnoe3),abs(integral15(Func3, a, b, ni, p)-Tochnoe3)/Tochnoe3
    write (*,*)

    p=3+0.01*N
    h=0.005
    ni=int((b-a)/h)

    write (*,*)  " int4Trapecia:" ,Tochnoe4,  integral11(Func4, a, b, ni, p),&
    abs(integral11(Func4, a, b, ni, p)-Tochnoe4),abs(integral11(Func4, a, b, ni, p)-Tochnoe4)/Tochnoe4

    h=0.001
    ni=int((b-a)/h)
    write (*,*)  " int4Simpson:" ,Tochnoe4, integral15(Func4, a, b, ni, p),&
    abs(integral15(Func4, a, b, ni, p)-Tochnoe4),abs(integral15(Func4, a, b, ni, p)-Tochnoe4)/Tochnoe4

   ! close (2016)
end subroutine rezFile6








