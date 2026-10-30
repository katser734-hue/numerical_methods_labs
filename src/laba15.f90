subroutine LabNumber15
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB15_N") /=0 ) read (str,*) tmp_str, N

    end do
    close (2000)
end subroutine LabNumber15

real function funk(x, y, c)
    implicit none
    real x, y, c
    funk = y**2/4/c-x/4/c-5*c/4/x/x

end function funk


subroutine rezFile15
    use module_main
    implicit none
    real funk
    real :: tochnoeX(11), tochnoeY(11)
    real ::  chitY(11), chADAM(11), chMiln(11)
    real cc, del
    real K1, K2, K3, K4, KK
    cc = 1 + 0.001*N


    do i=1, 11
        tochnoeX(i) = 1+0.1*(i-1)
    end do

    do i=1, 11
        tochnoeY(i) = (cc + tochnoeX(i)**(1.5))/tochnoeX(i)
    end do
    chitY(1) = 1+cc
    do i=2, 11
        K1 = 0.1*funk(tochnoeX(i-1), tochnoeY(i-1), cc)
        K2 = 0.1*funk(tochnoeX(i-1)+0.1/2, tochnoeY(i-1)+K1/2, cc)
        K3 = 0.1*funk(tochnoeX(i-1)+0.1/2, tochnoeY(i-1)+K2/2, cc)
        K4 = 0.1*funk(tochnoeX(i-1)+0.1, tochnoeY(i-1)+K3, cc)
        chitY(i) = chitY(i-1) + 1/6.0*(K1+2*K2+2*K3+K4)
    end do

    do i=1, 11
        write(*,*) tochnoeX(i),'&', tochnoeY(i),'&', chitY(i),'&', abs( tochnoeY(i)- chitY(i)),'&',&
         abs( tochnoeY(i)- chitY(i))/tochnoeY(i), '\\ \hline'
    end do


    do i=1, 11
        write(*,*) tochnoeX(i), tochnoeY(i)
    end do
    write(*,*)
    
    do i=1, 11
        write(*,*) tochnoeX(i), chitY(i)
    end do
    write(*,*)
    chADAM(1) = chitY(1)
    chADAM(2) = chitY(2)
    chADAM(3) = chitY(3)
    chADAM(4) = chitY(4)

    do i=5, 11
        K1 = funk(tochnoeX(i-1), chADAM(i-1), cc)
        K2 = funk(tochnoeX(i-2), chADAM(i-2), cc)
        K3 = funk(tochnoeX(i-3), chADAM(i-3), cc)
        K4 = funk(tochnoeX(i-4), chADAM(i-4), cc)
        del = (55*K1-59*K2+37*K3-9*K4)/24.0*0.1
        chADAM(i) = chADAM(i-1) + del
        KK= funk(tochnoeX(i), chADAM(i), cc)
        chADAM(i) = chADAM(i-1) +(9*KK+19*K1-5*K2+K3)/240
    end do


    do i=1, 11
        write(*,*) tochnoeX(i),'&', tochnoeY(i),'&', chADAM(i),'&', abs( tochnoeY(i)- chADAM(i)),'&',&
         abs( tochnoeY(i)- chADAM(i))/tochnoeY(i), '\\ \hline'
    end do
    write(*,*) 

  

    write(*,*) 
    write(*,*) 
    write(*,*) 

    chMiln(1) = chitY(1)
    chMiln(2) = chitY(2)
    chMiln(3) = chitY(3)
    chMiln(4) = chitY(4)

    do i=5, 11
        K1 = funk(tochnoeX(i-1), chMiln(i-1), cc)
        K2 = funk(tochnoeX(i-2), chMiln(i-2), cc)
        K3 = funk(tochnoeX(i-3), chMiln(i-3), cc)

        del = (2*K3-K2+2*K1)*4*0.1/3.0
        chMiln(i) = chMiln(i-4) + del
        KK= funk(tochnoeX(i), chMiln(i), cc)

        chMiln(i) = chMiln(i-2)+0.1/3*(K2+4*K1+KK)
    end do


    do i=1, 11
        write(*,*) tochnoeX(i),'&', tochnoeY(i),'&', chMiln(i),'&', abs( tochnoeY(i)- chMiln(i)),'&',&
         abs( tochnoeY(i)- chMiln(i))/tochnoeY(i), '\\ \hline'
    end do
    write(*,*) 

    do i=1, 11
        write(*,*) tochnoeX(i),chMiln(i)
    end do

    write(*,*) 

    do i=1, 11
        write(*,*) tochnoeX(i), abs( tochnoeY(i)- chMiln(i))
    end do
    
end subroutine rezFile15