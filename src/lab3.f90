subroutine LabNumber3
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB3_N") /=0 ) read (str,*) tmp_str, N
        if (index(str,"LAB3_T") /=0 ) read (str,*) tmp_str, T
        if (index(str,"LAB3_ni") /=0 ) read (str,*) tmp_str, ni
    end do
    close (2000)
end subroutine LabNumber3

real function F( a, x00)
    implicit none
    real::  a, x00
    if( x00 < -a .and. x00 >= -2 ) F = 0 
    if( x00 >=  a .and. x00 <=  2 ) F = 0 
    if( x00 >=  -a .and. x00<  0 ) F =  x00 + a 
    if( x00 >=  0 .and. x00 < a ) F = -x00 + a
end function F

real function F1(x00, pi, T, k, a)
    implicit none
    real F
    real::  x00,  pi, a
    integer  :: T,k
    F1 = F(a, x00)*cos(k*pi*x00/T*2)
end function F1

real function F2(x00, pi, T, k, a)
    implicit none
    real F
    real:: x00, pi, a
    integer  :: T, k
    F2 = F(a,x00)*sin(k*pi*x00/T*2)
end function F2

real function integral_a ( a1, b1, ni, pi, T, k, a)
    implicit none
    real F1
    real:: a1, b1, h, h1, pi, a
    integer  :: i, ni, T, k
    h=(b1-a1)/ni
    h1=(b1-a1)/ni
    integral_a = F1(a1,  pi,T, k, a)+F1(b1, pi,T, k, a)
    do i= 1, ni-2
        integral_a=integral_a+2*F1(a1+h1,  pi,T, k, a)
        h1=h1+h
    end do
    integral_a=integral_a*h/2
end function integral_a

real function integral_b ( a2, b2, ni, pi, T, k, a)
    implicit none
    real F2
    real:: a2, b2, h, h1, pi, a
    integer  :: i , ni, T, k
    h=(b2-a2)/ni
    h1=(b2-a2)/ni
    integral_b = F2(a2, pi, T, k, a)+F2(b2, pi, T, k, a)
    do i= 1, ni-2
        integral_b=integral_b+2*F2(a2+h1, pi, T, k,a)
        h1=h1+h
    end do
    integral_b=integral_b*h/2
end function integral_b

real function Fure ( ni, x0, pi, T, a, nch)
    implicit none
    real integral_b, integral_a
    real::   x0, pi, a, tut
    integer  ::  T, nch, ni, k
    tut=T/2
    Fure=integral_a( -tut, tut, ni, pi, T, 0, a)/2
    do k = 1, nch
        Fure=Fure+integral_a( -tut, tut, ni, pi, T, k, a)*cos(k*pi*x0/2)/2+integral_b( -tut, tut, ni, pi, T, k, a)*sin(k*pi*x0/2)/2
    end do
end function Fure

subroutine rezFile3
    use module_main
    implicit none
    real , external::  F1, F2, integral_a, integral_b, F
    integer :: nlab3 =39
    real  Fure
    
    real, dimension(39) :: vector_xx
    real, dimension(39) :: vector_yy
    real, dimension(39) :: vector_y8
    real, dimension(39) :: vector_y9
    open (unit = 2005, file = trim(resFileName))
    a=0.2+0.002*N  
    nch=8
    do i = 1, nlab3  
        vector_xx(i) = -2+0.1*i    
    end do

    do i = 1, nlab3  
        if( vector_xx(i) < -a .and. vector_xx(i) >= -2 ) vector_yy(i) = 0 
        if( vector_xx(i) >=  a .and. vector_xx(i) <=  2 ) vector_yy(i) = 0 
        if( vector_xx(i) >=  -a .and. vector_xx(i)<  0 ) vector_yy(i) =  vector_xx(i) + a 
        if( vector_xx(i) >=  0 .and. vector_xx(i) < a ) vector_yy(i) = -vector_xx(i) + a
    end do

    do i= 1, nlab3
        vector_y8(i)=Fure( ni, vector_xx(i), pi, T, a, nch )
        vector_y9(i)=Fure( ni, vector_xx(i), pi, T, a, nch+1 )
    end do

    do i= 1, nlab3
        write(*,*) "x", vector_xx(i) , "y", vector_yy(i), "y8", vector_y8(i), "y9",vector_y9(i)
    end do

    do i= 1, nlab3 
       write(*,*) "y-y8", abs(vector_yy(i)-vector_y8(i)),&
         "y-y9", abs(vector_yy(i)-vector_y9(i)), "y8-y9", abs(vector_y8(i)-vector_y9(i))
    end do

    do i= 1, nlab3
        write(*,*) "y-y8", abs(vector_yy(i)-vector_y8(i)), "y-y9",&
         abs(vector_yy(i)-vector_y9(i)), "y8-y9", abs(vector_y8(i)-vector_y9(i))
   end do
    do i= 1, nlab3
        write(*,*) "y-dy8", abs(vector_yy(i)-vector_y8(i))/vector_yy(i) , "y-dy9", &
        abs(vector_yy(i)-vector_y9(i))/vector_yy(i)
    end do

    close (2005)
end subroutine rezFile3


real function Legandr ( n, matrix, ST)
    implicit none
    real matrix(ST, ST)

    integer    n, k, ST
    real chis1, chis2

    chis1=(2.0*(n-2)+1)/(n-1)
    chis2=-(n-2.0)/(n-1)

    do k = 1, ST
        matrix(n, k)= chis2*matrix(n-2,k)
        if(k/=1) matrix(n, k)=matrix(n, k)+chis1*matrix(n-1,k-1)
    end do
    Legandr=chis1
end function Legandr



real function integral3 (  k, ni, matrix, ST,a)
    implicit none
    real F3
    real matrix(ST, ST)
    real:: a3, b3, h, h1,a
    integer  :: i, ni, k, ST
    a3=-1.0
    b3=1.0
    h=(b3-a3)/ni
    h1=(b3-a3)/ni
    integral3 = F3( k, a3, matrix, ST,a)+F3( k,b3, matrix, ST,a)
    do i= 1, ni-2
        integral3=integral3+2*F3(k, a3+h1, matrix, ST, a)
        h1=h1+h
    end do
    integral3=integral3*h/2.0
end function integral3


real function integral5 (  k, ni, matrix, ST)
    implicit none
    real F5
    real matrix(ST, ST)
    real:: a3, b3, h, h1
    integer  :: i, ni, k, ST
    a3=-1
    b3=1
    h=(b3-a3)/ni
    h1=(b3-a3)/ni
    integral5 = F5( k, a3, matrix, ST)+F5( k,b3, matrix, ST)
    do i= 1, ni-2
        integral5=integral5+2*F5(k, a3+h1, matrix, ST)
        h1=h1+h
    end do
    integral5=integral5*h/2.0
end function integral5

real function F5( k, x00, matrix, ST)
    implicit none
    real::  x00, pomog
    integer k, i, ST
    real matrix(ST, ST)
    pomog=0
    do i=1, ST
       if(matrix(k,i)/=0) pomog=pomog+matrix(k,i)*(x00**(i-1))  
    end do

    F5=pomog*pomog
end function F5

real function F3( k, x00, matrix, ST, a)
    implicit none
    real::  x00, pomog, a
    integer k, i, ST
    real matrix(ST, ST)
    F3=abs(x00**2 + x00*(3. - 2.*a) + 3.75 - 3.*a)

    pomog=0
    do i=1, ST
       if(matrix(k,i)/=0) pomog=pomog+matrix(k,i)*(x00**(i-1))  
    end do

    F3=F3*pomog
end function F3


subroutine rezFile33
    use module_main
    implicit none
    real, dimension(ST,ST):: matrix
    real, dimension(ST1,ST1):: matrix1
    real, external :: F3, F5
    real Legandr, pom, integral3,integral5
    real, dimension(ST) :: vector_aa
    real, dimension(ST1) :: vector_aa1
    real, dimension(5) :: vecY
    real, dimension(5) :: vecX
    real, dimension(5) :: vecY1
    real, dimension(5) :: vecY2
    real, dimension(5) :: vecT

    do i = 1, ST  
        do j=1, ST 
            matrix(i,j)=0
        end do
    end do
    a=1.5-0.02*N
    matrix(1, 1)=1
    matrix(2, 1)=0; matrix(2,2)=1
    
    do i=3, ST
        pom=Legandr( i, matrix, ST)
    end do

    do k = 1, ST
        vector_aa(k)=integral3(k, ni, matrix, ST,a)/integral5(k, ni, matrix, ST)
    end do

    do j=1, 5
        vecY1(j)=0
    end do
  
    vecT(1)=-1
    do j=2, 5
        vecT(j)=vecT(j-1)+0.5
    end do

    vecX(1)=0.5
    do j=2, 5
        vecX(j)=vecX(j-1)+0.5
    end do

    do j=1, 5
        vecY(j)=abs(vecX(j)**2-2*a*vecX(j)+1.5)
    end do

    do k=1, 5
        do i = 1, ST
            do j=1, ST 
                vecY1(k)=vecY1(k)+matrix(i,j)*vector_aa(i)*vecT(k)**(j-1)
            end do
        end do
    end do
  

    do i = 1, ST1  
        do j=1, ST1 
            matrix1(i,j)=0
        end do
    end do
    a=1.5-0.02*N
    matrix1(1, 1)=1
    matrix1(2, 1)=0; matrix1(2,2)=1
    
    do i=3, ST1
        pom=Legandr( i, matrix1, ST1)
    end do

    do k = 1, ST1
        vector_aa1(k)=integral3(k, ni, matrix1, ST1,a)/integral5(k, ni, matrix1, ST1)
    end do

    do j=1, 5
        vecY2(j)=0
    end do
  
    vecT(1)=-1
    do j=2, 5
        vecT(j)=vecT(j-1)+0.5
    end do

    vecX(1)=0.5
    do j=2, 5
        vecX(j)=vecX(j-1)+0.5
    end do

    do j=1, 5
        vecY(j)=abs(vecX(j)**2-2*a*vecX(j)+1.5)
    end do

    do k=1, 5
        do i = 1, ST1
            do j=1, ST1 
                vecY2(k)=vecY2(k)+matrix1(i,j)*vector_aa1(i)*vecT(k)**(j-1)
            end do
        end do
    end do


    write(*,*) 
    write(*,*) "vecT", vecT
    write(*,*) "vecX",vecX
    write(*,*) "vecY", vecY
    write(*,*) "vecL5",vecY1
    write(*,*) "vecL10",vecY2

    do j=1, 5
        write(*,*) "x, y-L5, y-L10,L5-L10:", vecX(j), abs(vecY(j)-vecY1(j)), abs(vecY(j)-vecY2(j)), abs(vecY1(j)-vecY2(j))
    end do

    
    do j=1, 5 
        write(*,*) "x, abs(y-L5), abs(y-L10):", vecX(j), abs(vecY(j)-vecY1(j))/abs(vecY(j)), abs(vecY(j)-vecY2(j))/abs(vecY(j))
    end do

    do j=1, 5
        !write(*,*)  vecX(j), abs(vecY(j)-vecY1(j))/abs(vecY(j))
    end do
    write(*,*)
    do j=1, 5
        !write(*,*)  vecX(j), abs(vecY(j)-vecY2(j))/abs(vecY(j))
    end do
    write(*,*)
  
end subroutine rezFile33








real function ChebOrt(KNC, pi, x00)
    real  pi, x00
    integer KNC
    if(KNC==0) then
        ChebOrt=1./sqrt(pi)
    else
        ChebOrt=sqrt(2./pi)*cos(KNC*acos(x00) )
    end if
end function ChebOrt

real function ChebOrtFAN(KNC, a, pi, x00)
    real a, pi, x00
    integer KNC
    if(KNC==0) then
        ChebOrtFAN=(1./sqrt(pi))*abs(x00**2 + x00*(3. - 2.*a) + 3.75 - 3.*a)/sqrt(1-x00**2)
    else
        ChebOrtFAN=sqrt(2./pi)*cos(KNC*acos(x00) )*abs(x00**2 + x00*(3. - 2.*a) + 3.75 - 3.*a)/sqrt(1-x00**2)
    end if
end function ChebOrtFAN


real function integralA (  KNC, ni, pi, a)
    implicit none
    real ChebOrtFAN
    real:: a1, b1, h, h1, pi, a
    integer  :: i, ni, KNC
    a1=-1.0
    b1=1.0
    h=(b1-a1)/ni
    h1=(b1-a1)/ni
    integralA = ChebOrtFAN( KNC,a, pi, a1+h1)+ChebOrtFAN( KNC,a, pi, b1-h1)
    do i= 1, ni-2
        integralA=integralA+2*ChebOrtFAN(KNC, a, pi,  a1+h1)
        h1=h1+h
    end do
    integralA=integralA*h/2.0

end function integralA

subroutine rezFile333
    use module_main
    implicit none
    real, external::  ChebOrtFAN
    real  integralA, ChebOrt
    real, dimension(5) :: vecY
    real, dimension(5) :: vecX
    real, dimension(5) :: vecY1, vecY2
    real, dimension(5) :: vecT

    a = 1.5-0.02*N

    vecT(1)=-1
    do j=2, 5
        vecT(j)=vecT(j-1)+0.5
    end do

    vecX(1)=0.5
    do j=2, 5
        vecX(j)=vecX(j-1)+0.5
    end do

    do j=1, 5
        vecY(j)=abs(vecX(j)**2-2*a*vecX(j)+1.5)
    end do

    do j=1,5
        vecY1(j)=0
        do i=0, KNC
            vecY1(j)=vecY1(j)+integralA(  i, ni, pi, a)*ChebOrt(i, pi, vecT(j))
         end do
    end do

    do j=1,5
        vecY2(j)=0
        do i=0, KNC1
            vecY2(j)=vecY2(j)+integralA(  i, ni, pi, a)*ChebOrt(i, pi, vecT(j))
         end do
    end do


    write(*,*) 
    write(*,*) "vecT", vecT
    write(*,*) "vecX",vecX
    write(*,*) "vecY", vecY
    write(*,*) "vecT5",vecY1
    write(*,*) "vecT10",vecY2

    do j=1, 5
        write(*,*) "x, y-T5, y-T10,T5-T10:", vecX(j), abs(vecY(j)-vecY1(j)), abs(vecY(j)-vecY2(j)), abs(vecY1(j)-vecY2(j))
    end do

    do j=1, 5 
        write(*,*) "x, abs(y-T5), abs(y-T10):", vecX(j), abs(vecY(j)-vecY1(j))/abs(vecY(j)), abs(vecY(j)-vecY2(j))/abs(vecY(j))
    end do

    do j=1, 5
        write(*,*)  vecX(j), abs(vecY(j)-vecY1(j))/abs(vecY(j))
    end do
    write(*,*)
    do j=1, 5
        write(*,*)  vecX(j), abs(vecY(j)-vecY2(j))/abs(vecY(j))
    end do
    write(*,*)
end subroutine rezFile333

