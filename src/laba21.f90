subroutine LabNumber21
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB21_N") /=0 ) read (str,*) tmp_str, N

    end do
    close (2000)
end subroutine LabNumber21



subroutine rezFile21
    use module_main
    implicit none

    real  :: Ut(11), iks(11), aa(11), bb(11)
    real :: Utn(11, 2001),  Uch(11, 2001)

    real (16)::alfa, sigma, krot, xs, ts, tema
    real :: x1 = 0
    real :: x2 = 1
    real :: t1 = 0
    real :: t2 = 20
    real :: hx = 0.1
    real :: ht = 0.01
    


    xs = int((x2-x1)/hx)+1
    ts = int((t2-t1)/ht)+1

    krot = ht/hx/hx
    alfa = -2.8+0.5*N
    sigma = 0.5




    if (sigma >= 0.5-hx**2/4/ht) then
        write(*,*) 'SXEMA YSTOUCHIVA'
    else
        write(*,*) 'SXEMA HEYSTOUCHIVA!!!!!!'
        pause
    end if

    do i = 1, 11
        iks(i) = (i-1)*0.1
    end do
   
    do i = 1, 11
        Ut(i) = alfa*(3-2*iks(i))
    end do

    do i = 1, 11
        write(*,*) iks(i), Ut(i) 
    end do


    write(*,*)
    write(*,*)

    do i = 1, 11
        do j=1, 2001
            Utn(i, j) = Ut(i)
        end do
    end do 

    do i=1, 11
        Uch(i, 1) = 2.0*alfa
    end do

    do j=1, 2001
        Uch(1, j) = 3.0*alfa
        Uch(11, j) = alfa
    end do


    do j=1, 2000
        do i=2, 10
            if (i==2) then 
                tema = 1.0+2.0*sigma*krot 
                aa(i) = sigma*krot/tema
                bb(i) = (Uch(i, j) +(1.0-sigma)*krot*(Uch(i+1, j)-2.0*Uch(i, j)+Uch(i-1, j))+sigma*krot*3.0*alfa)/tema
            else
                tema = 1.0+2.0*sigma*krot - sigma*krot*aa(i-1)
                aa(i) = sigma*krot/tema
                bb(i) = (Uch(i, j) +(1.0-sigma)*krot*(Uch(i+1, j)-2.0*Uch(i, j)+Uch(i-1, j))+sigma*krot*bb(i-1))/tema 
            end if
        end do
        
        do i = 10, 2, -1
            Uch(i, j+1) = aa(i)*Uch(i+1, j+1)+bb(i)
        end do
    end do

    write(*,*)
    write(*,*)  

    do i = 1, 11
        j = 2001
        write(*,*) iks(i),'&', Utn(i, j),'&', Uch(i, j),'&',abs(Utn(i, j)-Uch(i, j)),&
        '&', abs(Utn(i, j)-Uch(i, j))/abs(Utn(i, j)), '\\ \hline'
    end do
end subroutine rezFile21