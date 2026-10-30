subroutine LabNumber20
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB20_N") /=0 ) read (str,*) tmp_str, N

    end do
    close (2000)
end subroutine LabNumber20



subroutine rezFile20
    use module_main
    implicit none

    real :: Utcn(11, 6)
    real :: U(11, 6)

    real ::alfa

    alfa = -0.5+0.001*N

    do i = 1, 11
        do j = 1, 6
            Utcn(i, j) = 3+(2*alfa+1)*((i-1)*0.1)+alfa*(((j-1)*0.2)**2)/2-alfa*(j-1)*0.2
        end do
    end do

    do i = 1, 11
        write(*,*) '0.1','&', Utcn(i, 1),'&', Utcn(i, 2),'&', Utcn(i, 3),'&', Utcn(i, 4),&
        '&', Utcn(i, 5),'&', Utcn(i, 6), '\\ \hline'
    end do

    write(*,*)
    write(*,*)

    do j = 1, 11
        U(j, 1) = 3+(2*alfa+1)*0.1*(j-1)
    end do

    do j = 1, 6
        U(1, j) = 3 - alfa*(j-1)*0.2+alfa*((j-1)*0.2)**2/2
    end do

    do i = 2, 11
        do j =2, 6
            U (i, j) = (U(i-1, j)*(alfa/0.1/(2*alfa+1)-1/0.2)+alfa*(j-1)*0.2+U(i-1, j-1)/0.2)*(2*alfa+1)*0.1/alfa
        end do
    end do

    do i = 1, 11
        j = 1
        write(*,*) '0.1','&', U(i, j),'&', U(i, j+1),'&', U(i, j+2),'&', U(i, j+3),'&',&
         U(i, j+4),'&', U(i, j+5), '\\ \hline'
    end do
    write(*,*)
    write(*,*)
    !погрешность
    do i = 1, 11
        write(*,*) '0.1','&', abs(Utcn(i, 1) - U(i, 1)),'&', abs(Utcn(i, 2)-U(i, 2)),'&', abs(Utcn(i, 3)-U(i, 3)),'&',&
         abs(Utcn(i, 4)-U(i,4)),'&', abs(Utcn(i, 5)-U(i, 5)),'&', abs(Utcn(i, 6)-U(i, 6)), '\\ \hline'
    end do

    write(*,*)
    write(*,*)


    do i = 1, 11
        write(*,*) '0.1','&', abs(Utcn(i, 1) - U(i, 1))/abs(Utcn(i, 1)),'&', abs(Utcn(i, 2)-U(i, 2))/abs(Utcn(i, 2)),'&',&
         abs(Utcn(i, 3)-U(i, 3))/abs(Utcn(i, 3)),'&', abs(Utcn(i, 4)-U(i,4))/abs(Utcn(i, 4)),'&', &
         abs(Utcn(i, 5)-U(i, 5))/abs(Utcn(i, 5)),'&', abs(Utcn(i, 6)-U(i, 6))/abs(Utcn(i, 6)), '\\ \hline'
    end do
end subroutine rezFile20