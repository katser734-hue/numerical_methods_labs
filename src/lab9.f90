 subroutine LabNumber9
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB9_ni") /=0 ) read (str,*) tmp_str, ni
    end do
    close (2000)
end subroutine LabNumber9

real function Fun91(y, x)
    implicit none
    real:: y, x
    Fun91=(x**2+2*y)
end function Fun91

real function Fun92(y, x)
    implicit none
    real:: y, x
    Fun92=1.0/(x+y)**2
end function Fun92

real function Fun93(y, x)
    implicit none
    real:: y, x
    Fun93=x**2/(1+y**2)
end function Fun93

real function integral91Pomog (FF, a, b, n, x)
    implicit none
    real FF
    integer n
    real:: a, b, h, x
    integer :: i
    integral91Pomog = 0
    h=(b-a)/n
    do i=0, n-1
        integral91Pomog=integral91Pomog+FF(a+h*i+h/2.0, x)
    end do
    integral91Pomog=integral91Pomog*h
end function integral91Pomog

real function integral92Pomog (FF, a, b, n, x)
    implicit none
    real,external :: FF
    real:: a, b, h, x
    integer  :: i , n
    ! interface func_
    !     real function FF(y, x)
    !         real:: y, x
    !     end function FF
    ! end interface
    h=(b-a)/n
    integral92Pomog = FF( a, x)/2.0+FF(b, x)/2.0
    do i=1, n-1
        integral92Pomog=integral92Pomog+FF(a+i*h, x)
    end do
    integral92Pomog=integral92Pomog*h
end function integral92Pomog

real function integral93Pomog ( FF, a, b, n, x)
    implicit none
    real    :: FF
    real    :: a, b, h, x
    integer :: i, n
    h=(b-a)/n
    integral93Pomog = FF(a, x)+FF(b, x)
    do i=1, n-1
        if(mod(i, 2)==0) then
            integral93Pomog=integral93Pomog+2*FF(a+i*h, x)
        else 
            integral93Pomog=integral93Pomog+4*FF( a+i*h, x)
        end if
    end do
    integral93Pomog=integral93Pomog*h/3
end function integral93Pomog


real function integral91 (F, FF, a1, b1, a2, b2, n)
    implicit none
    real F !, FF
    integer n
    real:: a1, b1, a2, b2, h
    integer :: i

    interface func_
        real function FF(y, x)
            real:: y, x
        end function FF
    end interface

    integral91 = 0
    h=(b1-a1)/n
    do i=0, n-1
        integral91=integral91+F(FF, a2, b2, n, a1+h*i+h/2.0)
    end do
    integral91=integral91*h
end function integral91

real function integral92 (F, FF, a1, b1, a2, b2, n)
    implicit none
    real F !, FF
    real:: a1, b1, a2, b2, h
    integer  :: i , n

    interface func_
        real function FF(y, x)
            real:: y, x
        end function FF
    end interface

    h=(b1-a1)/n
    integral92 = F(FF, a2, b2, n, a1)/2.0+F(FF, a2, b2, n, b1)/2.0
    do i=1, n-1
        integral92=integral92+F(FF, a2, b2, n,a1+i*h)
    end do
    integral92=integral92*h
end function integral92


real function integral93 (F, FF, a1, b1, a2, b2, n)
    implicit none
    real,external :: F , FF
    real:: a1, b1, a2, b2, h
    integer i, n

    ! interface func_
    !     real function FF(y, x)
    !         real:: y, x
    !     end function FF
    ! end interface

    h=(b1-a1)/n
    integral93 = F(FF, a2, b2, n, a1)+F(FF, a2, b2, n, b1)
    do i=1, n-1
        if(mod(i, 2)==0) then
            integral93=integral93+2*F(FF, a2, b2, n, a1+i*h)
        else 
            integral93=integral93+4*F(FF, a2, b2, n, a1+i*h)
        end if
    end do
    integral93=integral93*h/3
end function integral93

subroutine init_random_seed()
    integer :: iii, net, clock
    integer, dimension(:), allocatable :: seed
    call random_seed(size = net)
    allocate(seed(net))
    call system_clock(count=clock)
    seed = clock + 37 * (/ (iii - 1, iiii = 1, net) /)
    call random_seed(put = seed)
    deallocate(seed)
end subroutine init_random_seed




subroutine rezFile9
    use module_main
    implicit none
    real Rand1(1000000), Rand2(1000000), Rand3(1000000)
    real Rand11, Rand22, Rand33

    real , external:: integral91Pomog, integral92Pomog, integral93Pomog
    real , external:: Fun92, Fun93, Fun91 
    real  integral91 , integral92, integral93
    real  Tochnoe91, Tochnoe92, Tochnoe93
    real gr0, gr1, gr2, gr3, gr4

    real  Carlik1, Carlik2, Carlik3
    integer Chetchik1,  Chetchik2,  Chetchik3
    real  Carl1, Carl2, Carl3

    ! interface func_
    !     real function Fun91(y, x)
    !         real:: y, x
    !     end function Fun91
    ! end interface
   

    Tochnoe91 = 4.666666667
    Tochnoe92 = 0.040821994
    Tochnoe93 = 0.107250184
    gr0=0
    gr1=1
    gr2=2
    gr3=3
    gr4=4

    write(*,*) "Int1: Chet, abs, otn", Tochnoe91
    write(*,*) "Center", integral91(integral91Pomog, Fun91, gr0, gr2, gr0, gr1, ni)&
     ,abs(Tochnoe91 - integral91(integral91Pomog, Fun91, gr0, gr2, gr0, gr1, ni))&
     ,abs(Tochnoe91 - integral91(integral91Pomog, Fun91, gr0, gr2, gr0, gr1, ni))/Tochnoe91
    write(*,*) "Trap", integral92(integral92Pomog, Fun91, gr0, gr2, gr0, gr1, 1000)&
     ,abs(Tochnoe91 - integral92(integral92Pomog, Fun91, gr0, gr2, gr0, gr1, 1000))&
     ,abs(Tochnoe91 - integral92(integral92Pomog, Fun91, gr0, gr2, gr0, gr1, 1000))/Tochnoe91
    write(*,*) "Simp", integral93(integral93Pomog, Fun91, gr0, gr2, gr0, gr1, ni)&
     ,abs(Tochnoe91 - integral93(integral93Pomog, Fun91, gr0, gr2, gr0, gr1, ni))&
     ,abs(Tochnoe91 - integral93(integral93Pomog, Fun91, gr0, gr2, gr0, gr1, ni))/Tochnoe91
     write(*,*)
     write(*,*)"Int2: Chet, abs, otn", Tochnoe92
     write(*,*) "Center", integral91(integral91Pomog, Fun92, gr3, gr4, gr1, gr2, ni)&
     ,abs(Tochnoe92 - integral91(integral91Pomog, Fun92, gr3, gr4, gr1, gr2, ni))&
     ,abs(Tochnoe92 - integral91(integral91Pomog, Fun92, gr3, gr4, gr1, gr2, ni))/Tochnoe92
    write(*,*) "Trap", integral92(integral92Pomog, Fun92, gr3, gr4, gr1, gr2, ni)&
     ,abs(Tochnoe92 - integral92(integral92Pomog, Fun92, gr3, gr4, gr1, gr2, ni))&
     ,abs(Tochnoe92 - integral92(integral92Pomog, Fun92, gr3, gr4, gr1, gr2, ni))/Tochnoe92
    write(*,*) "Simp",  integral93(integral93Pomog, Fun92, gr3, gr4, gr1, gr2, ni)&
     ,abs(Tochnoe92 - integral93(integral93Pomog, Fun92, gr3, gr4, gr1, gr2, ni))&
     ,abs(Tochnoe92 - integral93(integral93Pomog, Fun92, gr3, gr4, gr1, gr2, ni))/Tochnoe92
     write(*,*)
     write(*,*)"Int3: Chet, abs, otn", Tochnoe93
     write(*,*) "Center", integral91(integral91Pomog, Fun93, gr0, gr1, gr1, gr2, ni)&
     ,abs(Tochnoe93 - integral91(integral91Pomog, Fun93, gr0, gr1, gr1, gr2, ni))&
     ,abs(Tochnoe93 - integral91(integral91Pomog, Fun93, gr0, gr1, gr1, gr2, ni))/Tochnoe93
    write(*,*) "Trap", integral92(integral92Pomog, Fun93, gr0, gr1, gr1, gr2, ni)&
     ,abs(Tochnoe93 - integral92(integral92Pomog, Fun93, gr0, gr1, gr1, gr2, ni))&
     ,abs(Tochnoe93 - integral92(integral92Pomog, Fun93, gr0, gr1, gr1, gr2, ni))/Tochnoe93
    write(*,*) "Simp", integral93(integral93Pomog, Fun93, gr0, gr1, gr1, gr2, ni)&
     ,abs(Tochnoe93 - integral93(integral93Pomog, Fun93, gr0, gr1, gr1, gr2, ni))&
     ,abs(Tochnoe93 - integral93(integral93Pomog, Fun93, gr0, gr1, gr1, gr2, ni))/Tochnoe93
     write(*,*)

    Carl1=0.261799
    Carl2=0.5
    Carl3=0.034074


     Chetchik1=0
     Chetchik2=0
     Chetchik3=0
    Carlik1=0
    Carlik2=0
    Carlik3=0


    call init_random_seed()
    do i=1, 10000
        call random_number(Rand11)
        call random_number(Rand22)
        call random_number(Rand33)
        Rand1(i)=Rand11
        Rand2(i)=Rand22
        Rand3(i)=Rand33
        if( Rand2(i) <= Rand1(i)    ) then 
            Carlik1 = Carlik1 + sqrt(Rand1(i)**2 - Rand2(i)**2)
            Chetchik1=Chetchik1 + 1
        end if
        if( Rand2(i)**2 >= Rand1(i) ) then 
            Carlik2 = Carlik2 + exp(Rand1(i)/Rand2(i))            
            Chetchik2=Chetchik2 + 1
        end if
        if( Rand1(i) <= 1 - Rand2(i)- Rand3(i) .and. Rand2(i) <=  1 - Rand3(i)) then
            Carlik3 = Carlik3 + 1.0/(Rand1(i) + Rand2(i)+Rand3(i)+1)**3 
            Chetchik3=Chetchik3 + 1
        end if
    end do
    Carlik1=Carlik1/2.0/Chetchik1
    Carlik2=Carlik2/3.0/Chetchik2
    Carlik3=Carlik3/6.0/Chetchik3

    write(*,*) "N=10 000:"
    write(*,*) "Carl1:", Carl1, Carlik1, abs(Carl1-Carlik1), abs(Carl1-Carlik1)/Carl1
    write(*,*) "Carl2:", Carl2, Carlik2, abs(Carl2-Carlik2), abs(Carl2-Carlik2)/Carl2
    write(*,*) "Carl3:", Carl3, Carlik3, abs(Carl3-Carlik3), abs(Carl3-Carlik3)/Carl3

    Chetchik1=0
    Chetchik2=0
    Chetchik3=0
    Carlik1=0
    Carlik2=0
    Carlik3=0

    call init_random_seed()
    do i=1, 100000
        call random_number(Rand11)
        call random_number(Rand22)
        call random_number(Rand33)
        Rand1(i)=Rand11
        Rand2(i)=Rand22
        Rand3(i)=Rand33
        if( Rand2(i) <= Rand1(i)    ) then 
            Carlik1 = Carlik1 + sqrt(Rand1(i)**2 - Rand2(i)**2)
            Chetchik1=Chetchik1 + 1
        end if
        if( Rand2(i)**2 >= Rand1(i) ) then 
            Carlik2 = Carlik2 + exp(Rand1(i)/Rand2(i))            
            Chetchik2=Chetchik2 + 1
        end if
        if( Rand1(i) <= 1 - Rand2(i)- Rand3(i) .and. Rand2(i) <=  1 - Rand3(i)) then
            Carlik3 = Carlik3 + 1.0/((Rand1(i) + Rand2(i)+Rand3(i)+1)**3) 
            Chetchik3=Chetchik3 + 1
        end if
    end do
    Carlik1=Carlik1/2.0/Chetchik1
    Carlik2=Carlik2/3.0/Chetchik2
    Carlik3=Carlik3/6.0/Chetchik3
    
    write(*,*) "N=100 000:"
    write(*,*) "Carl1:", Carl1, Carlik1, abs(Carl1-Carlik1), abs(Carl1-Carlik1)/Carl1
    write(*,*) "Carl2:", Carl2, Carlik2, abs(Carl2-Carlik2), abs(Carl2-Carlik2)/Carl2
    write(*,*) "Carl3:", Carl3, Carlik3, abs(Carl3-Carlik3), abs(Carl3-Carlik3)/Carl3


    Chetchik1=0
    Chetchik2=0
    Chetchik3=0
    Carlik1=0
    Carlik2=0
    Carlik3=0

    call init_random_seed()
    do i=1, 1000000
        call random_number(Rand11)
        call random_number(Rand22)
        call random_number(Rand33)
        Rand1(i)=Rand11
        Rand2(i)=Rand22
        Rand3(i)=Rand33
        if( Rand2(i) <= Rand1(i)    ) then 
            Carlik1 = Carlik1 + sqrt(Rand1(i)**2 - Rand2(i)**2)
            Chetchik1=Chetchik1 + 1
        end if
        if( Rand2(i)**2 >= Rand1(i) ) then 
            Carlik2 = Carlik2 + exp(Rand1(i)/Rand2(i))            
            Chetchik2=Chetchik2 + 1
        end if
        if( Rand1(i) <= 1 - Rand2(i)- Rand3(i) .and. Rand2(i) <=  1 - Rand3(i)) then
            Carlik3 = Carlik3 + 1.0/(Rand1(i) + Rand2(i)+Rand3(i)+1)**3 
            Chetchik3=Chetchik3 + 1
        end if
    end do
    Carlik1=Carlik1/2.0/Chetchik1
    Carlik2=Carlik2/3.0/Chetchik2
    Carlik3=Carlik3/6.0/Chetchik3

    write(*,*) "N=1 000 000:"
    write(*,*) "Carl1:", Carl1, Carlik1, abs(Carl1-Carlik1), abs(Carl1-Carlik1)/Carl1
    write(*,*) "Carl2:", Carl2, Carlik2, abs(Carl2-Carlik2), abs(Carl2-Carlik2)/Carl2
    write(*,*) "Carl3:", Carl3, Carlik3, abs(Carl3-Carlik3), abs(Carl3-Carlik3)/Carl3

end subroutine rezFile9