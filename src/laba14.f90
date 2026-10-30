subroutine LabNumber14
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB14_N") /=0 ) read (str,*) tmp_str, N

    end do
    close (2000)
end subroutine LabNumber14





real function Fu11(x0, y2)
    implicit none
    real:: x0, y2
    Fu11=-0.37
end function Fu11

real function Fu22(x0, y3)
    implicit none
    real:: x0, y3

    Fu22=0.255 
end function Fu22

real function Fu33(x0, y1)
    implicit none
    real:: x0, y1

    Fu33= -1.62
end function Fu33







real function Fu1(x0, y2)
    implicit none
    real:: x0, y2
    Fu1=(y2-1.5*sqrt(x0))/x0   

end function Fu1

real function Fu2(x0, y3)
    implicit none
    real:: x0, y3

    Fu2=(y3+2*sqrt(x0))/x0   
end function Fu2

real function Fu3(x0, y1)
    implicit none
    real:: x0, y1

    Fu3=(y1-1.5*sqrt(x0))/x0   
end function Fu3












real function integ (F, a, b, y)
    implicit none
    real F
    real:: a, b, h, y
    integer  :: i , n
 
    n=16384
   ! n=8388608
    h=(b-a)/n
    integ = 0
    do i=1, n
        integ=integ+F(a+i*h,y)
    end do
    integ=integ*h
end function integ





subroutine rezFile14
    use module_main
    implicit none

    real , external:: Fu1,Fu2,Fu3
    real , external:: Fu11,Fu22,Fu33
    real  integ


    real coef1, coef, pom1, pom2, pom3
    real pp1(17), pp2(17), pp3(17)
    real ppp1(17), ppp2(17), ppp3(17)
    real mX(17)
    real Toc1, Toc2, Toc3
    real  igrek0,  igrek1,  igrek2,  igrek3,  igrek4
    real iksik
    real TochnX(10)
    real ChitaemX(10)
    integer :: iks = 2
    real  koef 
    


    
    pom1 = sin(sqrt(14.0)*log(2.0)/10)
    pom2 = cos(sqrt(14.0)*log(2.0)/10)
    coef = -(76*pom1-2*sqrt(14.0)*pom2)/(2**0.2*sqrt(14.0)*pom1**2+2**0.2*sqrt(14.0)*pom2**2)
    coef1 = (2*sqrt(14.0)*pom1+76*pom2)/(2**0.2*sqrt(14.0)*pom1**2+2**0.2*sqrt(14.0)*pom2**2)




    do i = 1, 10
        iksik = 2+0.01*i
        ChitaemX(i) = 2+4*(iksik-2)+1.5*(iksik-2)**2+0.216666*(iksik-2)**3+0.0352*(iksik-2)**4
        TochnX(i) = iksik**(0.2)*(coef1*sin(sqrt(14.0)*log(iksik)/10.0)+coef*cos(sqrt(14.0)*log(iksik)/10.0))
    end do


    ! do i = 1, 10
    !     write(*,*) ChitaemX(i), abs(TochnX(i)-ChitaemX(i))/abs(TochnX(i))
    ! end do
    ! write(*,*)



    !ЗАДАННИЕ 2

    iks = 2
    igrek0 = 1.0/3
    igrek1 = -1.0/3
    igrek2 = -igrek0*(41.1*iks-0.1)-igrek1*(41.1*iks**2-0.2*iks+0.1)
    igrek2 = igrek2/(-6.85*iks**3-0.05*iks**2+0.05*iks-6.9)
    pom1 = igrek0*(41.1*iks-0.1)+igrek1*(41.1*iks**2-0.2*iks+0.1)
    pom2  = -6.85*iks**3-0.05*iks**2+0.05*iks-6.9
    igrek3 = -(41.1*igrek0+(41.1*iks-0.1)*igrek1+(82.2*iks-0.2)*igrek1+(41.1*iks**2-0.2*iks+0.1)*igrek2)*pom2
    igrek3 = igrek3+pom1*(-6.85*3*iks**2-0.1*iks+0.05)
    igrek3 = igrek3/pom2**2


    write(*,*)



    do i = 1, 10
        iksik = 2+0.01*i
        ChitaemX(i) = igrek0 + igrek1*(iksik-iks)+ igrek2*(iksik-iks)**2 /2+ igrek3*(iksik-iks)**3 /6
    end do


    coef = -6.85*iks**3-0.05*iks**2+0.05*iks-6.9
    coef1 = -6.85*3*iks**2-0.1*iks+0.05
    pom1 = coef1*1.0/3-coef*1.0/3
    pom2 = coef*1.0/3-2*pom1

    do i = 1, 10
        iksik = 2+0.01*i
        TochnX(i) = (iksik*pom1+pom2)/(-6.85*iksik**3-0.05*iksik**2+0.05*iksik-6.9)
    end do


    do i = 1, 10
        iksik = 2+0.01*i
        !write(*,*) iksik,'&', TochnX(i),'&', ChitaemX(i),'&', abs(TochnX(i)-ChitaemX(i)),'&',&
        !abs(TochnX(i)-ChitaemX(i))/abs(TochnX(i)),  '\\ \hline'
        !write(*,*) iksik, ChitaemX(i)
    end do

































    Ar = -1+0.001*N

    Toc1 = 0.64*Ar+0.8
    Toc2 = 0.64*Ar+1.6
    Toc3 = 0.64*Ar-0.8

    do i=1, 17
        mX(i) = 0.64+i*0.01
    end do

    do i=1, 17
        pp1(i) = integ(Fu1, 0.64, mX(i), 0.255*(mX(i)-0.64)+Toc2) + Toc1
        pp2(i) = integ(Fu2, 0.64, mX(i), -1.62*(mX(i)-0.64)+Toc3) + Toc2
        pp3(i) = integ(Fu3, 0.64, mX(i), -0.37*(mX(i)-0.64)+Toc1) + Toc3
    end do

    ! do i=1, 17
    !     write(*,*) mX(i), -sqrt(mX(i))+Ar*mX(i)
    ! end do

    do i = 1, 10
        write(*,*) mX(i), abs(-sqrt(mX(i))+Ar*mX(i)-pp3(i))/abs(-sqrt(mX(i))+Ar*mX(i))
    end do
    write(*,*)
 

    ! do i=1, 17
    !     write(*,*) mX(i),'&', sqrt(mX(i))+Ar*mX(i),'&', pp1(i),'&', abs(sqrt(mX(i))+Ar*mX(i)-pp1(i)),'&', &
    !     abs(sqrt(mX(i))+Ar*mX(i)-pp1(i))/abs(sqrt(mX(i))+Ar*mX(i)), '\\ \hline'
    ! end do
    ! write(*,*) 

    ! do i=1, 17
    !     write(*,*) mX(i),'&', 2*sqrt(mX(i))+Ar*mX(i),'&', pp2(i),'&', abs(2*sqrt(mX(i))+Ar*mX(i)-pp2(i)),'&', &
    !     abs(2*sqrt(mX(i))+Ar*mX(i)-pp2(i))/abs(2*sqrt(mX(i))+Ar*mX(i)), '\\ \hline'
    ! end do
    ! write(*,*) 

    ! do i=1, 17
    !     write(*,*) mX(i),'&', -sqrt(mX(i))+Ar*mX(i),'&', pp3(i),'&', abs(-sqrt(mX(i))+Ar*mX(i)-pp3(i)),'&', &
    !     abs(-sqrt(mX(i))+Ar*mX(i)-pp3(i))/abs(-sqrt(mX(i))+Ar*mX(i)), '\\ \hline'
    ! end do
    ! write(*,*) 

    ! write(*,*) 



    do j = 1, 1
        write(*,*) 

        write(*,*) 


        do i=1, 17
            ppp1(i) = pp1(i)
            ppp2(i) = pp2(i)
            ppp3(i) = pp3(i)


            pp1(i) = ppp2(i)*(log(mX(i))-log(0.64))-3*sqrt(mX(i))+3*sqrt(0.64)  + Toc1
            pp2(i) = ppp3(i)*(log(mX(i))-log(0.64))+4*sqrt(mX(i))-4*sqrt(0.64) + Toc2
            pp3(i) = ppp1(i)*(log(mX(i))-log(0.64))-3*sqrt(mX(i))+3*sqrt(0.64) + Toc3
        end do
        write(*,*) 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAaa:'

        do i=1, 17
            write(*,*) mX(i),'&', sqrt(mX(i))+Ar*mX(i),'&',&
             sqrt(mX(i))+Ar*mX(i),'&', abs(sqrt(mX(i))+Ar*mX(i)-pp1(i))/1000,'&', &
            abs(sqrt(mX(i))+Ar*mX(i)-pp1(i))/abs(sqrt(mX(i))+Ar*mX(i))/1000, '\\ \hline'
        end do
        write(*,*) 
    
        do i=1, 17
            write(*,*) mX(i),'&', 2*sqrt(mX(i))+Ar*mX(i),'&',&
             2*sqrt(mX(i))+Ar*mX(i),'&', abs(2*sqrt(mX(i))+Ar*mX(i)-pp2(i))/1000,'&', &
            abs(2*sqrt(mX(i))+Ar*mX(i)-pp2(i))/abs(2*sqrt(mX(i))+Ar*mX(i))/1000, '\\ \hline'
        end do
        write(*,*) 
    
        do i=1, 17
            write(*,*) mX(i),'&', -sqrt(mX(i))+Ar*mX(i),'&',&
             -sqrt(mX(i))+Ar*mX(i),'&', abs(-sqrt(mX(i))+Ar*mX(i)-pp3(i))/1000,'&', &
            abs(-sqrt(mX(i))+Ar*mX(i)-pp3(i))/abs(-sqrt(mX(i))+Ar*mX(i))/1000, '\\ \hline'
        end do
        write(*,*) 
    end do
    
 
end subroutine rezFile14