subroutine LabNumber22
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB22_N") /=0 ) read (str,*) tmp_str, N

    end do
    close (2000)
end subroutine LabNumber22



subroutine rezFile22
    use module_main
    implicit none

    real  ::  iks(11), te(11)
    real :: Utn(11, 11),  Uch(11, 11)

    real ::alfa
    real ::tay, th




    alfa = 2+0.01*N
    th = alfa/2+1
    tay = th/10
    h = 1.0/10


    do i = 1, 11
        iks(i) = (i-1)*h
        te(i) = tay*(i-1)
    end do

    if ( tay**2/alfa**2/h**2 <= 1 ) then
        write(*,*) 'SXEMA YSTOUCHIVA'
    else
        write(*,*) 'SXEMA HEYSTOUCHIVA!!!!!!'
        pause
    end if

    do i = 1, 11
        do j=1, 11
            Uch(i, j) = 0
        end do
    end do 


    do i = 1, 11
        do j=1, 11
            Utn(i, j) = 2*cos(pi*te(j)/alfa)*sin(pi*iks(i))
        end do
    end do 

    do j=1, 11
        Uch(1, j) = 0
    end do

    do j=1, 11
        Uch(11, j) = 0
    end do

    do i=1, 11
        Uch(i, 1) = 2*sin(pi*iks(i))
    end do
 
    
    do i=1, 11
        Uch(i, 2) = 2*sin(pi*iks(i))-tay**2*pi**2*sin(pi*iks(i))/alfa**2
    end do    





    do j = 3, 11
        do i=2, 10
            Uch(i, j) = tay**2/(alfa**2*h**2)*(Uch(i+1, j-1)-2.0*Uch(i, j-1)+Uch(i-1, j-1))+2*Uch(i, j-1)- Uch(i, j-2)
        end do
    end do 





    do i = 1, 11
        j = 11
        write(*,*) iks(i),'&', Utn(i, j),'&', Uch(i, j),'&',abs(Utn(i, j)-Uch(i, j)),&
        '&', abs(Utn(i, j)-Uch(i, j))/abs(Utn(i, j)), '\\ \hline'
    end do
end subroutine rezFile22