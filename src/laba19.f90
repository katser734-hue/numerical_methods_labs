subroutine LabNumber19
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB19_N") /=0 ) read (str,*) tmp_str, N

    end do
    close (2000)
end subroutine LabNumber19



subroutine rezFile19
    use module_main
    implicit none

    real :: ik(11)
    real ::  yik(11)
    real ::  chyk(11)
    real :: bet, kt
    real :: ah = 0.1
    bet = 1+0.001*N
    ik(1) = 0.0

    do i=2, 11
        ik(i) = ik(i-1)+0.1
    end do

    yik(2) = 1/(bet-2)+(bet-1)/(bet-2)*ah+(2*bet-3)/(bet-2)*ah**2/2+(3*bet-5)/(bet-2)*ah**3/6+(4*bet-7)/(bet-2)*ah**4/24

    yik(1) = 1/(bet-2)


    do i=3, 11
        kt = ik(i)*exp(ik(i))+10*ik(i-1)*exp(ik(i-1))+ik(i-2)*exp(ik(i-2))
        yik(i) = yik(i-1)*(24+10*ah*ah*(2*bet-3))/(12-ah*ah*(2*bet-3))-yik(i-2)-(2*ah*ah*(bet-2)*kt/(12-ah*ah*(2*bet-3)))
    end do

    do i=1, 11
        chyk(i) = (ik(i) + 1/(bet-2))*exp(ik(i))
    end do


    do i=1, 11
        write(*,*) ik(i),'&', chyk(i),'&', yik(i), '&',  abs(chyk(i)-yik(i)),'&',  abs(chyk(i)-yik(i))/abs(chyk(i)), '\\ \hline'
    end do
    write(*,*) 
    write(*,*) 
    
    write(*,*) ik(2),'&', chyk(2),'&', yik(2), '&',  abs(chyk(2)-yik(2)),'&',  abs(chyk(2)-yik(2))/abs(chyk(2)), '\\ \hline'
    write(*,*) ik(6),'&', chyk(6),'&', yik(6), '&',  abs(chyk(6)-yik(6)),'&',  abs(chyk(6)-yik(6))/abs(chyk(6)), '\\ \hline'
    write(*,*) ik(11),'&', chyk(11),'&', yik(11), '&',  abs(chyk(11)-yik(11)),'&',  abs(chyk(11)-yik(11))/abs(chyk(11)), '\\ \hline'

end subroutine rezFile19