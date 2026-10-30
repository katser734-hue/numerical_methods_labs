subroutine LabNumber5
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB5_N") /=0 ) read (str,*) tmp_str, N
        if (index(str,"LAB5_ni") /=0 ) read (str,*) tmp_str, ni
    end do
    close (2000)
end subroutine LabNumber5


real function GAUS(matrix,vector, ni, RESH)
    integer  :: ni
    real matrix(ni, ni)
    real vector(ni)
    real RESH(ni)
    integer ::  k, i, j
    real :: upir
    GAUS=1
    do k=1, ni     
        upir=matrix(k, k)
        vector(k)=vector(k)/matrix(k,k)
        do j=k, ni 
            matrix(k,j)=matrix(k,j)/upir
        end do

        do i=k+1, ni
            vector(i)=vector(k)*matrix(i,k)-vector(i)
            upir=matrix(i,k)
            do j=k, ni 
                matrix(i,j)=matrix(k,j)*upir-matrix(i,j)
            end do
        end do
    end do
    RESH(ni)=vector(ni)
    i=ni-1
    do while(i>=1)
        upir=0
        j=ni
        do while(j>i)
            upir=upir+matrix(i,j )*RESH(j)
            j=j-1
        end do
        RESH(i) = vector(i) - upir
        i=i-1
    end do    
end function GAUS


subroutine rezFile5
    use module_main
    implicit none
    real GAUS
    real matrix1(ni, ni)
    real vector1(ni)
    real RESH(ni)
    real(4), dimension(ni) :: ZnachX
    real(4), dimension(ni) :: ZnachY
    real Pr1, Pr3, Tpr1, Tpr3
    real TRun, Run, Runyt
    h=0.25+(0.05+0.001*N)*1-0.25

    do i=1, ni
        ZnachX(i)=0.25+(0.05+0.001*N)*(i-1)
        ZnachY(i)=exp(2*sqrt(ZnachX(i))+1.0/2/sqrt(ZnachX(i)))
    end do

    write(*,*) "X", ZnachX
    write(*,*) "Y", ZnachY
    do i=1, ni
        write(*,*) ZnachX(i), ZnachY(i)
    end do


    do j=1, ni
        matrix1(1,j)=1
    end do
    do i = 2, ni    
        do j=1, ni
            matrix1(i,j)=(j-1)**(i-1)*h
        end do
    end do
    do i = 1, ni  
        vector1(i)= i-1
    end do
   x=GAUS(matrix1,vector1, ni, RESH)
   write(*,*) "CiX1", RESH
    Pr1=0
    do i=1, ni
        Pr1=Pr1+ZnachY(i)*RESH(i)
    end do


    do j=1, ni
        matrix1(1,j)=1
    end do
    do i = 2, ni    
        do j=1, ni
            matrix1(i,j)=(j-1)**(i-1)*h
        end do
    end do
    do i = 1, ni  
        vector1(i)= (i-1)*3**(i-2)
    end do
   x=GAUS(matrix1,vector1, ni, RESH)
   write(*,*) "CiX3", RESH
   Pr3=0
    do i=1, ni
        Pr3=Pr3+ZnachY(i)*RESH(i)
    end do


    Tpr1=exp(2*sqrt(ZnachX(2))+1.0/2/sqrt(ZnachX(2)))*(1.0/sqrt(ZnachX(2))-1.0/4/sqrt(ZnachX(2)**3))
    Tpr3=exp(2*sqrt(ZnachX(4))+1.0/2/sqrt(ZnachX(4)))*(1.0/sqrt(ZnachX(4))-1.0/4/sqrt(ZnachX(4)**3))
    write(*,*)"Pr1, Pr3, Tpr1, Tpr3:", Pr1, Pr3, Tpr1, Tpr3
    write(*,*)"ABSPr1, ABSPr3:", abs(Pr1-Tpr1), abs(Pr3-Tpr3)
    write(*,*)"OTNPr1, OTNPr3:", abs(Pr1-Tpr1)/Tpr1, abs(Pr3-Tpr3)/Tpr3
    


    open (unit = 2015, file = trim(resFileName))

    do j=1, 5
       write(2015,*) ZnachX(j), ZnachY(j)
    end do
    close (2015)
    

 
    do i=0, 7
        h=pi/800/(10**i)
        Trun=exp(sin(100*pi/200+N))*cos(100*pi/200+N)*100
        Run=1.0/2/h*(exp(sin(100*pi/200+N))-exp(sin(100*(pi/200-2*h)+N)))
        Runyt=(Run-(1.0/2/3.0/h*(exp(sin(100*pi/200+N))-exp(sin(100*(pi/200-2*3*h)+N)))))/8.0+Run
        
        write(*,*)"Trun, Run, Ab, Ot, RY, ab,ot:",&
            Trun, Run,abs(Run-Trun),abs(Run-Trun)/Trun, Runyt,abs(Runyt-Trun),abs(Runyt-Trun)/Trun, h
    end do
end subroutine rezFile5




