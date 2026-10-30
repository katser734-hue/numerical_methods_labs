subroutine LabNumber12
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB12_N") /=0 ) read (str,*) tmp_str, N
    end do
    close (2000)
end subroutine LabNumber12



subroutine rezFile12
    use module_main
    implicit none

    real (8) :: M, GAUSSS
    real (8) :: matrix1(5, 5)
    real (8) :: matrix2(5, 5)
    real (8) :: matrixM(5, 5)

    real (8) :: vectorik(5)
    real (8) :: vectorikP(5)
    real (8) :: vectorikR(5)
    real (8) :: vectorikTCH(5)

    real (8) :: alfa

    alfa = 0.01+0.002*N

    do i=1, 5 
        matrix1(i,i) = i*alfa**2
    end do
    do j=2, 5 
        matrix1(1,j) = 1+(j-1)/10.0
    end do
    do i=2, 5 
        matrix1(i,1) = ((6-i)/10.0+1)*alfa
    end do

    matrix1(2,3) = 1.1 + alfa; matrix1(2,4) = 1.2 + alfa; matrix1(2,5) = 1.3 + alfa
    matrix1(3,2) = 1.4*alfa;   matrix1(3,4) = 1.1 + alfa; matrix1(3,5) = 1.2 + 2*alfa
    matrix1(4,2) = 1.3*alfa;   matrix1(4,3) = 1.4*alfa;   matrix1(4,5) = 1.1 + 3*alfa
    matrix1(5,2) = 1.2*alfa;   matrix1(5,3) = 1.3*alfa;   matrix1(5,4) = 1.4*alfa;


    vectorikP(1) = 3.21*alfa**2 + 1.934
    vectorikP(2) = 19.58*alfa**2 - 1.326*alfa - 8.253
    vectorikP(3) = 16.74*alfa**2 - 0.631*alfa - 13.251
    vectorikP(4) = -17.16*alfa**2 + 3.061*alfa - 7.821
    vectorikP(5) = -35.55*alfa**2 + 16.5271*alfa 

    write(*,*) "matrix A1:"
    do i = 1, 5
        write(*,*) matrix1(i,:)
    end do


    write(*,*) "Vect B1:"
    write(*,*) vectorikP

    do i=1, 5      
        do j=1, 5 
            matrixM(i,j) = matrix1(i, j)
        end do
    end do

    do i=1, 4
        y = M(matrixM, i)

        do k=1, 5    
            vectorik(k) = 0  
            do j=1, 5 
                vectorik(k) = vectorik(k) + matrixM(k,j) * vectorikP(j)
            end do
        end do


        do k=1, 5 
            vectorikP(k) = vectorik(k) 
        end do

        matrixM = matmul(matrixM, matrix1)

        do k=1, 5      
            do j=1, 5 
                matrix1(k,j) = matrixM(k, j)
            end do
        end do
    end do


    z = GAUSSS(matrix1, vectorik, 5, vectorikR)


    write(*,*) "Resh X:"
    write(*,*) vectorikR




    alfa = 0.02+0.001*N

    matrix1(1,1)=10; matrix1(1,2)=alfa+0.3; matrix1(1,3)=alfa+0.2; matrix1(1,4)=alfa+0.1; matrix1(1,5)=alfa+1

    matrix1(2,1)=10*alfa-3; matrix1(2,2)=alfa**2+9.91; matrix1(2,3)=alfa**2+0.9*alfa+0.24
    matrix1(2,4)=alfa**2+0.8*alfa+0.17; matrix1(2,5)=alfa**2+1.7*alfa+1.7

    matrix1(3,1)=10*alfa-2; matrix1(3,2)=alfa**2+10.1*alfa-3.06; matrix1(3,3)=2*alfa**2+9.87
    matrix1(3,4)=2*alfa**2+0.8*alfa+0.22; matrix1(3,5)=2*alfa**2+3.5*alfa+0.22

    matrix1(4,1)=10*alfa-1; matrix1(4,2)=alfa**2+10.2*alfa-2.03; matrix1(4,3)=2*alfa**2+10.2*alfa-3.08
    matrix1(4,4)=3*alfa**2+9.86; matrix1(4,5)=3*alfa**2+6.4*alfa+2.6

    matrix1(5,1)=10*alfa-10; matrix1(5,2)=alfa**2+9.3*alfa-20.3; matrix1(5,3)=2*alfa**2+7.5*alfa-30.8
    matrix1(5,4)=3*alfa**2+4.6*alfa-41.4; matrix1(5,5)=4*alfa**2-20


    vectorik(1) = 35.07*alfa + 25.07
    vectorik(2) = 35.07*alfa**2 + 49.379*alfa + 40.122
    vectorik(3) = 69.9*alfa**2 + 88.55*alfa + 63.1691
    vectorik(4) = 103.2*alfa**2 + 155.906*alfa + 169.1416
    vectorik(5) = 124.38*alfa**2 + 71.759*alfa - 979.664

    write(*,*) "matrix A2:"
    do i = 1, 5
        write(*,*) matrix1(i,:)
    end do

    write(*,*) "Vect B2:"
    write(*,*) vectorik



    do i=1, 5      
        do j=1, 5 
            matrix2(i,j) = matrix1(i, j)
        end do
    end do

    do i=1, 5      
        vectorikP(i) = vectorik(i)
    end do
    z = GAUSSS(matrix2, vectorikP, 5, vectorikR)

    write(*,*) "Vect X TCH:"
    write(*,*) vectorikR

    do k=1, 5 
        vectorikTCH(k) = vectorikR(k) 
    end do




    do j=1, 5 
        matrix1(5,j) = matrix1(4, j)*4+ matrix1(5,j)       
    end do

    do j=1, 5 
        matrix1(5,j) = matrix1(3, j)*4 + matrix1(5,j)   
    end do
    
    do j=1, 5 
        matrix1(5,j) = matrix1(2, j)*3.8 + matrix1(5,j)  
    end do
    
    do j=1, 5 
        matrix1(5,j) = matrix1(1, j)*3 + matrix1(5,j)  
    end do

    vectorik(5) = vectorik(1)*3 +vectorik(5)
    vectorik(5) = vectorik(2)*3.8 +vectorik(5)
    vectorik(5) = vectorik(3)*4 +vectorik(5)
    vectorik(5) = vectorik(4)*4 +vectorik(5)
    
    write(*,*) "matrix A2 IZM:"
    do i = 1, 5
        write(*,*) matrix1(i,:)
    end do

    write(*,*) "Vect B2 IZM:"
    write(*,*) vectorik


    do i = 1, 5
        vectorikR(i) = 0
    end do


    Epsik = 1
    ff=0
    do while (abs(Epsik) > 0.001)

        do i=1, 5
            vectorikP(i) = vectorik(i)
            do j=1, 5
                if(j/=i) vectorikP(i) = vectorikP(i)  - matrix1(i, j)*vectorikR(j)
            end do
            vectorikP(i) = vectorikP(i)/matrix1(i, i)
        end do

        do i = 1, 5
            Pamat12(i) = vectorikR(i) 
        end do

        do i = 1, 5
            vectorikR(i) = vectorikP(i) 
        end do

        Epsik=max(abs(Pamat12(1)-vectorikR(1)), abs(Pamat12(2)-vectorikR(2)),&
        abs(Pamat12(3)-vectorikR(3)),abs(Pamat12(4)-vectorikR(4)),abs(Pamat12(5)-vectorikR(5)))

        ff=ff+1
    end do

    write(*,*) "Vect X PROST ITER:"
    write(*,*) vectorikR
    write(*,*) "Kolvo iter:", ff
    write (*,*) "Pogresh:"
    do i = 1, 5
        write (*,*) "x", i, ":  " , abs(vectorikR(i) - vectorikTCH(i))
    end do





    do i = 1, 5
        vectorikR(i) = 0
    end do


    Epsik = 1
    ff=0
    do while (abs(Epsik) > 0.001)
        do i=1, 5
            vectorikP(i) = vectorik(i)
            do j=1, 5
                if(j/=i) vectorikP(i) = vectorikP(i)  - matrix1(i, j)*vectorikR(j)
            end do
            vectorikP(i) = vectorikP(i)/matrix1(i, i)
            Pamat12(i) = vectorikR(i) 
            vectorikR(i) = vectorikP(i) 
        end do

        Epsik=max(abs(Pamat12(1)-vectorikR(1)), abs(Pamat12(2)-vectorikR(2)),&
        abs(Pamat12(3)-vectorikR(3)),abs(Pamat12(4)-vectorikR(4)),abs(Pamat12(5)-vectorikR(5)))

        ff=ff+1
    end do

    write(*,*) "Vect X metod ZEi:"
    write(*,*) vectorikR
    write(*,*) "Kolvo iter:", ff

    write (*,*) "Pogresh:"
    do i = 1, 5
        write (*,*) "x", i, ":  " , abs(vectorikR(i) - vectorikTCH(i))
    end do


end subroutine rezFile12