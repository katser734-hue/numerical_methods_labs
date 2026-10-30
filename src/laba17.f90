subroutine LabNumber17
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB17_N") /=0 ) read (str,*) tmp_str, N

    end do
    close (2000)
end subroutine LabNumber17

real function funk17(v)
    implicit none
    real v
    funk17 = v
end function funk17

real function funk172(u)
    implicit none
    real u
    funk172 = 4*u
end function funk172

subroutine rezFile17
    use module_main
    implicit none
    real funk17, funk172
    real :: tochnoeX(11), tochnoeU(11), tochnoeV(11), RazU(11), RazV(11)
    real ::  chitU(11), chitV(11)
    real K1, K2, K3, K4



    do i=1, 11
        tochnoeX(i) = 1+0.01*(i-1)
    end do

    do i=1, 11
        tochnoeU(i) = 2*exp(tochnoeX(i)*2-2) - exp(2-2*tochnoeX(i))
        tochnoeV(i) = 4*exp(tochnoeX(i)*2-2) + 2*exp(2-2*tochnoeX(i))
    end do
    chitU(1) = 1
    chitV(1) = 6
    do i=2, 11
        K1 = 0.01*funk17(tochnoeV(i-1))
        K2 = 0.01*funk17(tochnoeV(i-1)+K1/2)
        K3 = 0.01*funk17(tochnoeV(i-1)+K2/2)
        K4 = 0.01*funk17(tochnoeV(i-1)+K3)
        chitU(i) = chitU(i-1) + 1/6.0*(K1+2*K2+2*K3+K4)
    end do

    
    do i=2, 11
        K1 = 0.01*funk172(chitU(i-1))
        K2 = 0.01*funk172(chitU(i-1)+K1/2)
        K3 = 0.01*funk172(chitU(i-1)+K2/2)
        K4 = 0.01*funk172(chitU(i-1)+K3)
        chitV(i) = chitV(i-1) + 1/6.0*(K1+2*K2+2*K3+K4)
    end do

    do i=2, 11
        K1 = 0.01*funk17(chitV(i-1))
        K2 = 0.01*funk17(chitV(i-1)+K1/2)
        K3 = 0.01*funk17(chitV(i-1)+K2/2)
        K4 = 0.01*funk17(chitV(i-1)+K3)
        chitU(i) = chitU(i-1) + 1/6.0*(K1+2*K2+2*K3+K4)
    end do



    do i=1, 11
        !write(*,*) tochnoeX(i),'&', tochnoeU(i),'&', chitU(i),'&', abs( tochnoeU(i)- chitU(i)),'&',&
        ! abs( tochnoeU(i)- chitU(i))/tochnoeU(i), '\\ \hline'
    end do
    write(*,*) 
    write(*,*) 
    
    write(*,*) 
    do i=1, 11
        !write(*,*) tochnoeX(i),'&', tochnoeV(i),'&', chitV(i),'&', abs( tochnoeV(i)- chitV(i)),'&',&
        ! abs( tochnoeV(i)- chitV(i))/tochnoeV(i), '\\ \hline'
    end do
    write(*,*) 
 

    RazU(1) = 1
    RazV(1) = 6

    RazU(2) = chitU(2)
    RazV(2) = chitV(2)

    RazU(3) = chitU(3)
    RazV(3) = chitV(3)

    do i = 3, 11
        RazU(i) = RazU(i-2) + 2*0.01/N*(RazV(i) +(N-2)*RazV(i-1) + RazV(i-2))
        RazV(i) = RazV(i-2) + 8*0.01/N*(0.5*RazU(i) +(N-1)*RazU(i-1) + RazU(i-2))
        RazU(i) = RazU(i-2) + 2*0.01/N*(RazV(i) +(N-2)*RazV(i-1) + RazV(i-2))
    end do

    
    write(*,*) 
    write(*,*) 
    write(*,*) 
    do i=1, 11
        write(*,*) tochnoeX(i),'&', tochnoeU(i),'&', RazU(i),'&', abs( tochnoeU(i)- RazU(i)),'&',&
         abs( tochnoeU(i)- RazU(i))/tochnoeU(i), '\\ \hline'
    end do
    write(*,*) 
    write(*,*) 
    do i=1, 11
        write(*,*) tochnoeX(i), abs( tochnoeU(i)- RazU(i))
    end do
    
    write(*,*) 
    
    
    write(*,*) 
    do i=1, 11
        write(*,*) tochnoeX(i),'&', tochnoeV(i),'&', RazV(i),'&', abs( tochnoeV(i)- RazV(i)),'&',&
         abs( tochnoeV(i)- RazV(i))/tochnoeV(i), '\\ \hline'
    end do
    



end subroutine rezFile17