subroutine LabNumber2
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB2_n") /=0 ) read (str,*) tmp_str, ni
 
    end do
    close (2000)
end subroutine LabNumber2


subroutine rezFile2
    use module_main
    implicit none
    real(4) :: xist2=0.0, xi=0.0, xi_yi=0.0, yi=0.0, lnxi=0.0, lnxi_yi=0.0, lnxist2=0.0, xist4=0.0, xist3=0.0, yi_xist2=0.0
    real(4) :: MNK1=0.0, MNK2=0.0, MNK3=0.0
    
    integer, dimension(ni) :: vector_y
    real, dimension(ni) :: vector_x

    real, dimension(ni) :: vector_y1
    real, dimension(ni) :: vector_y2
    real, dimension(ni) :: vector_y3

    real, dimension(3) :: AA
    real, dimension(3) :: BB

    open (unit = 2001, file = 'labnumber2.txt')
    vector_x(1) =2.1
    vector_x(2) =3.9
    vector_x(3) =5.1
    vector_x(4) =6.3
    vector_x(5) =7.0
    vector_x(6) =9.5
    vector_y(1) =0
    vector_y(2) =5
    vector_y(3) =10
    vector_y(4) =15
    vector_y(5) =20
    vector_y(6) =25
   
    
    close(2001)

    do i = 1, ni  
        xi=xi+vector_x(i) 
        xist2=xist2+vector_x(i)**2
        xi_yi=xi_yi+vector_x(i)*vector_y(i)
        yi=yi+vector_y(i)
        lnxi=lnxi+log(vector_x(i))
        lnxist2=lnxist2+log(vector_x(i))**2
        lnxi_yi=lnxi_yi+log(vector_x(i))*vector_y(i)
        xist3=xist3+vector_x(i)**3
        xist4=xist4+vector_x(i)**4
        yi_xist2=yi_xist2+vector_x(i)**2*vector_y(i)
        
    end do
    
    AA(1)=(xi_yi-xi*yi/ni)/(xist2-(xi**2)/ni)
    BB(1)=(yi-xi*AA(1))/ni

    BB(2)=(lnxi_yi-yi*lnxi/ni)/(lnxist2-lnxi**2/ni)
    AA(2)=(yi-lnxi*BB(2))/ni

    BB(3) =(xi_yi-yi_xist2*xist3/xist4)/(xist2-xist3**2/xist4)
    AA(3)=(yi_xist2-xist3*BB(3))/xist4

    
    do i = 1, ni  
        vector_y1(i)=AA(1)*vector_x(i)+BB(1)
        vector_y2(i)=AA(2)+log(vector_x(i))*BB(2)
        vector_y3(i)=AA(3)*(vector_x(i))**2+BB(3)*vector_x(i)
    end do

    do i = 1, ni  
        MNK1=MNK1+(vector_y1(i)-vector_y(i))**2
        MNK2=MNK2+(vector_y2(i)-vector_y(i))**2
        MNK3=MNK3+(vector_y3(i)-vector_y(i))**2
    end do

    do i = 1, 3  
        write(*,*) "AA:", AA(i)
    end do

    do i = 1, 3  
        write(*,*) "BB:", BB(i)
    end do


    open (unit = 2002, file = trim(resFileName))
    write(*,*) "MNK1:", MNK1
    write(2002,*) "MNK1:", MNK1
    
    write(*,*) "MNK2:", MNK2
    write(2002,*) "MNK2:", MNK2

    write(*,*) "MNK3:", MNK3
    write(2002,*) "MNK3:", MNK3

    close(2002)

end subroutine rezFile2
