subroutine LabNumber11
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB11_N") /=0 ) read (str,*) tmp_str, N
    end do
    close (2000)
end subroutine LabNumber11


real (8) function DETER(matrix)
    real (8) :: matrix(5, 5)
    integer ::  k, i, j
    real (8):: upir
    DETER=1
    do k=1, 5      
        upir=matrix(k, k)
        DETER=DETER*upir
        do j=k, 5 
            matrix(k,j)=matrix(k,j)/upir
        end do
        do i=k+1, 5 

            upir=matrix(i,k)
            do j=k, 5 
                matrix(i,j)=matrix(i,j)-matrix(k,j)*upir
            end do
        end do
    end do
      
end function DETER

real (8) function OBR(matrix2,matrix1)
    real (8 ) ::matrix2(5, 5)
    real (8)::matrix1(5, 5)
    integer ::  k, i, j
    real (8):: upir
    OBR=1
    do k=1, 5      
        upir=matrix2(k, k)
        do j=1, k-1
            matrix1(k,j)=matrix1(k,j)/upir
        end do
        do j=k, 5
            matrix2(k,j)=matrix2(k,j)/upir
            matrix1(k,j)=matrix1(k,j)/upir
        end do
        do i=k+1, 5 
            upir=matrix2(i,k)
            do j=1, k-1
                matrix1(i,j)=matrix1(k,j)*upir-matrix1(i,j)
            end do
            do j=k, 5 
                matrix2(i,j)=matrix2(k,j)*upir-matrix2(i,j)
                matrix1(i,j)=matrix1(k,j)*upir-matrix1(i,j)
            end do
        end do
    end do
    k=5
    do while(k>1)
        i=k-1
        do while(i>=1)
            upir=matrix2(i,k)
            do j=1, 5
                matrix2(i,j)=matrix2(i,j)-matrix2(k,j)*upir
                matrix1(i,j)=matrix1(i,j)-matrix1(k,j)*upir
            end do
            i=i-1
        end do
        k=k-1
    end do
    
end function OBR


subroutine rezFile11
    use module_main
    implicit none
    real (8) :: DETER, OBR
    real (8) :: matrix1(5, 5)
    real (8) :: matrix2(5, 5)
    real (8) :: alfa
    real (8) ::normM1, normL1,normK1
    real (8) ::normM2, normL2,normK2
    real (8) ::chit

    normM1=0; normL1=0; normK1=0
    normM2=0; normL2=0; normK2=0
    alfa = 0.1-0.001*N

    matrix1(1,1) = alfa**2 - 1
    matrix1(2,2) = alfa**2 - 3.99
    matrix1(3,3) = alfa**2 + 0.29
    matrix1(4,4) = alfa**2 - 14.91
    matrix1(5,5) = alfa**2 - 6.54

    do j=2, 5 
        matrix1(1,j) = (j-1)/10.0*alfa + (j-1)/10.0
    end do

    do i=2, 5 
        matrix1(i,1) = (i-1)/10.0*alfa - (i-1)/10.0
    end do

    matrix1(2,3)=0.5*alfa+1.02; matrix1(2,4)=0.6*alfa+1.23; matrix1(2,5)=0.7*alfa+1.44
    matrix1(3,2)=0.5*alfa-0.98; matrix1(3,4)=0.8*alfa+0.36; matrix1(3,5)=0.9*alfa+0.43
    matrix1(4,2)=0.6*alfa-1.17; matrix1(4,3)=0.8*alfa+0.36; matrix1(4,5)=alfa-2.74
    matrix1(5,2)=0.7*alfa-1.36; matrix1(5,3)=0.9*alfa+0.43; matrix1(5,4)=alfa+5.26

    write(*,*) "matrix C:"
    do i = 1, 5
        write(*,*) matrix1(i,:)
    end do

   y=DETER(matrix1)
   write(*,*)" DETERM:", y



   alfa = 1-0.01*N

    do i=1, 5 
        matrix1(i,i) = 5-i+1
    end do

    matrix1(1, 2)=-0.8/alfa; matrix1(1, 3)=0.2/alfa**2; matrix1(1, 4)=-1.0/15/alfa**3; matrix1(1, 5)=1.0/30/alfa**4
    matrix1(2, 1)=-20*alfa; matrix1(2, 3)=-1/alfa; matrix1(2, 4)=1.0/3/alfa**2; matrix1(2, 5)=-1.0/6/alfa**3
    matrix1(3, 1)=45*alfa**2; matrix1(3, 2)=-9*alfa; matrix1(3, 4)=-1.0/alfa; matrix1(3, 5)=0.5/alfa**2
    matrix1(4, 1)=-60*alfa**3; matrix1(4, 2)=12*alfa**2; matrix1(4, 3)=-4*alfa; matrix1(4, 5)=-1.0/alfa
    matrix1(5, 1)=30*alfa**4; matrix1(5, 2)=-6*alfa**3; matrix1(5, 3)=2*alfa**2; matrix1(5, 4)=-alfa

    write(*,*) "matrix A:"
    do i = 1, 5
        write(*,*) matrix1(i,:)
    end do


    do i = 1, 5 
        chit=0   
        do j=1, 5 
            chit = chit + abs(matrix1(i,j))
        end do
        if(chit > normM1) normM1 = chit
    end do
    
    do j = 1, 5 
        chit=0   
        do i=1, 5 
            chit = chit + abs(matrix1(i,j))
        end do
        if(chit > normL1) normL1 = chit
    end do

    chit=0
    do i = 1, 5  
        do j=1, 5 
            chit = chit + matrix1(i,j)**2
        end do
    end do
    
    normK1 = sqrt(chit)

    do i = 1, 5    
        do j=1, 5 
            if(i==j) then
                matrix2(i,j)=1
            else
                matrix2(i,j)=0
            end if
        end do
    end do

   z=OBR( matrix1, matrix2)

    write(*,*) "OBRATNA MATRIX A:"
    do i = 1, 5
        write(*,*) matrix2(i,:)
    end do

    do i = 1, 5 
        chit=0   
        do j=1, 5 
            chit = chit + abs(matrix2(i,j))
        end do
        if(chit > normM2) normM2 = chit
    end do
    
    do j = 1, 5 
        chit=0   
        do i=1, 5 
            chit = chit + abs(matrix2(i,j))
        end do
        if(chit > normL2) normL2 = chit
    end do

    chit=0
    do i = 1, 5  
        do j=1, 5 
            chit = chit + matrix2(i,j)**2
        end do
    end do
    
    normK2 = sqrt(chit)

    write(*,*) "Norm M, L, K"
    write(*,*) "matrix A:", normM1, normL1,normK1
    write(*,*) "matrix obr(A):", normM2, normL2,normK2



end subroutine rezFile11