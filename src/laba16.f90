subroutine LabNumber16
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB16_N") /=0 ) read (str,*) tmp_str, N

    end do
    close (2000)
end subroutine LabNumber16

real function funk16(x, y, c)
    implicit none
    real x, y, c
    funk16 = y**2/4/c-x/4/c-5*c/4/x/x

end function funk16


subroutine rezFile16
    use module_main
    implicit none
    real funk16
    real :: tochnoeX(11), tochnoeY(11)
    real ::  chitY(11), chMi(11), chN(11)
    real cc,  counters
    real K1, K2, K3, K4
    cc = 1 + 0.001*N


    do i=1, 11
        tochnoeX(i) = 1+0.1*(i-1)
    end do

    do i=1, 11
        tochnoeY(i) = (cc + tochnoeX(i)**(1.5))/tochnoeX(i)
    end do
    chitY(1) = 1+cc
    do i=2, 11
        K1 = 0.1*funk16(tochnoeX(i-1), tochnoeY(i-1), cc)
        K2 = 0.1*funk16(tochnoeX(i-1)+0.1/2, tochnoeY(i-1)+K1/2, cc)
        K3 = 0.1*funk16(tochnoeX(i-1)+0.1/2, tochnoeY(i-1)+K2/2, cc)
        K4 = 0.1*funk16(tochnoeX(i-1)+0.1, tochnoeY(i-1)+K3, cc)
        chitY(i) = chitY(i-1) + 1/6.0*(K1+2*K2+2*K3+K4)
    end do

    do i=1, 11
        !write(*,*) tochnoeX(i),'&', tochnoeY(i),'&', chitY(i),'&', abs( tochnoeY(i)- chitY(i)),'&',&
        ! abs( tochnoeY(i)- chitY(i))/tochnoeY(i), '\\ \hline'
    end do




    write(*,*)
    write(*,*)
    write(*,*)
    write(*,*)
    chN(1) = chitY(1) 
    chN(2) = chitY(2) 
    chN(3) = chitY(3) 
    chN(4) = chitY(4) 

    
    do i=4, 11
        counters = 0 
        counters = counters + funk16(tochnoeX(i-3), chN(i-3), cc)
        counters = counters - 2*funk16(tochnoeX(i-2), chN(i-2), cc)
        counters = counters + 7*funk16(tochnoeX(i-1), chN(i-1), cc)

        chN(i) = counters*0.1/3+chN(i-2)


    end do
    write(*,*) 
    do i=1, 11
        !write(*,*) tochnoeX(i),'&', tochnoeY(i),'&', chN(i),'&', abs( tochnoeY(i)- chN(i)),'&',&
        ! abs( tochnoeY(i)- chN(i))/tochnoeY(i), '\\ \hline'
    end do

    write(*,*)
    write(*,*)
    write(*,*)
    write(*,*)
    chMi(1) = chN(1)
    chMi(2) = chN(2)
    chMi(3) = chN(3)
    chMi(4) = chN(4)
    chMi(5) = chN(5)

    
    do i=5, 11
        counters = 0 
        counters = counters + 14*funk16(tochnoeX(i-4), chN(i-4), cc)
        counters = counters + 64*funk16(tochnoeX(i-3), chN(i-3), cc)
        counters = counters + 24*funk16(tochnoeX(i-2), chN(i-2), cc)
        counters = counters +64*funk16(tochnoeX(i-1),  chN(i-1), cc)
        counters = counters + 14*funk16(tochnoeX(i),   chN(i),   cc)
        chMi(i) = counters*0.1/45+chN(i-4)


    end do
    write(*,*) 
    do i=1, 11
        write(*,*) tochnoeX(i),'&', tochnoeY(i),'&', chMi(i),'&', abs( tochnoeY(i)- chMi(i)),'&',&
         abs( tochnoeY(i)- chMi(i))/tochnoeY(i), '\\ \hline'
    end do
    write(*,*)
    write(*,*)
    write(*,*)
    do i=1, 11
        write(*,*) tochnoeX(i),chMi(i)
    end do


end subroutine rezFile16