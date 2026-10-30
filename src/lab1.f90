subroutine LabN
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    read (2000,'(a)') str
    if (index(str,"LAB_WORK") /=0 ) read (str,*) tmp_str, numLab
    close (2000)
end subroutine LabN

    
subroutine NameEnd
    use module_main
    implicit none
    resFileName(:)   = " "
    resFileName(1:8) = "_resLab_"
    write(resFileName(9:10),"(i2.2)") numLab
    resFileName(11:14) =".txt"
end subroutine NameEnd

subroutine LabNumber1
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB1_N") /=0 ) read (str,*) tmp_str, N
        if (index(str,"LAB1_n") /=0 ) read (str,*) tmp_str, ni
        if (index(str,"LAB1_X1") /=0 ) read (str,*) tmp_str, diapozonX1
        if (index(str,"LAB1_X2") /=0 ) read (str,*) tmp_str, diapozonX2
    end do
    close (2000)
end subroutine LabNumber1

integer function Factik(n)
    integer n, k 
    Factik = 1
    k=n
    do while( k > 1 )
        Factik = k*Factik
        k=k-1
    end do
end  function Factik

subroutine rezFile1
    use module_main
    implicit none
    
    integer  KK
    real(8):: chiiss

    integer, allocatable :: BINk(:)
    real(8), allocatable :: Ai(:), aai(:),chis(:), Xi(:), Ypp(:), AH(:)


    real(8), dimension(ni, ni) :: Raznosti
    integer Factik

    real(4) :: start_time, finish_time

    allocate(vectorX(ni+1))
    allocate(vectorY(ni+1))

    allocate(Ai(0:ni), aai(0:ni), BINk(ni), chis(ni+1), Xi(ni), Ypp(ni+1),AH(0:ni))


    Ai(:)=0; aai(:)=0; BINk(:)=0

    h = (diapozonX2-diapozonX1)/ni
    do i=1, ni+1
        vectorX(i) = diapozonX1+h*(i-1)
        vectorY(i) = exp(N*vectorX(i)) 
    end do
    do i=1, ni+1
        write(*,*) vectorX(i), vectorY(i)
    end do

 
    write(*,*)"VNEPx", vectorX
    write(*,*)"VNEPy", vectorY

    KK=2**ni-1

    do ff=1, ni+1
        Ai(:)=0
        BINk(:)=0
        Ai(0)=1
        chis(ff)=1
        
        do i=1, ni
            if(ff/=i) then
               Xi(i) = vectorX(i)
            else 
                Xi(i) = vectorX(ni+1)
                chiiss = vectorX(i)
            end if  
        end do
        if(ff==ni+1) chiiss = vectorX(ff)

        do i=1, ni
            chis(ff)=chis(ff)* (chiiss-Xi(i))
        end do

        do k=1, KK
            
            BINk(ni)=BINk(ni)+1
            do j=ni, 1, -1
                if (BINk(j)==2) then
                    BINk(j)=0
                    BINk(j-1)=BINk(j-1)+1
                else
                    exit
                end if
            end do
            prod=1
            do j=1, ni
                prod=prod*Xi(j)**BINk(j)
            end do
            i=sum(BINk)
            Ai(i)=Ai(i)+prod
        end do
        
        
        do i=0, ni
            Ai(i)=Ai(i)/chis(ff)*vectorY(ff)
            if(mod(i,2)/=0) Ai(i)=-Ai(i)
        end do
     
        do i=0, ni
            aai(i)=aai(i)+Ai(i)
        end do 
        
    end do
    do i=0, ni
        write(*,*) "Ai", aai(i)
    end do
    summ=0

    do j=1, ni+1 
        summ=0
        do i=1, ni+1
            summ=summ+ aai(i-1)*vectorX(j)**(ni+1-i)
        end do 
       Ypp(j)=summ
       write(*,*) vectorX(j), Ypp(j)
    end do

    !Нью
    do j=1, ni
        do i=1, ni
            Raznosti(i, j) = 0
        end do
    end do

   
    do i=1, ni
        Raznosti(i, 1) = vectorY(i+1)-vectorY(i)
    end do


    k=1
    do j=2, ni
        do i=1, ni-k
            Raznosti(i, j) = Raznosti(i+1, j-1) - Raznosti(i, j-1)
        end do
        k=k+1
    end do

    do i=1, ni
       Xi(i) = vectorX(1)+h*(i-1)
    end do
   
    do i=1, ni
        !write(*,* )"XN", Xi(i)
    end do

aai(:)=0

    do nir=1, ni
        Ai(:)=0
        BINk(:)=0
        Ai(0)=1
        KK=2**(nir)-1
        do k=1, KK
            
            BINk(nir)=BINk(nir)+1
            do j=nir, 1, -1
                if (BINk(j)==2) then
                    BINk(j)=0
                    BINk(j-1)=BINk(j-1)+1
                else
                    exit
                end if
            end do
            prod=1
            do j=1, nir
                prod=prod*Xi(j)**BINk(j)
            end do
            i=sum(BINk)
            Ai(i)=Ai(i)+prod
        end do

        do i=0, nir
            if(mod(i,2)/=0) Ai(i)=-Ai(i)
            Ai(i)=Ai(i)*Raznosti(1,nir)/Factik(nir)/(h**(nir))
        end do
  

        do i=0, nir
            aai(i) = aai(i)+Ai(nir-i)
        end do 
   
    end do
    aai(0)=aai(0)+vectorY(1)

    do i=0, ni
        write(*,*) "NUTONAi1", aai(i)
    end do
    write(*,*)

    summ=0

    do j=1, ni+1 
        summ=0
        do i=1, ni+1
            summ=summ+ aai(i-1)*vectorX(j)**(i-1)
        end do 
       Ypp(j)=summ
       write(*,*) vectorX(j), Ypp(j)
    end do

    !Нью2
    do i=1, ni
        Xi(i) = 1-h*(i-1)
    end do

    aai(:)=0

    do nir=1, ni
        Ai(:)=0
        BINk(:)=0
        Ai(0)=1
        KK=2**(nir)-1
        do k=1, KK
            
            BINk(nir)=BINk(nir)+1
            do j=nir, 1, -1
                if (BINk(j)==2) then
                    BINk(j)=0
                    BINk(j-1)=BINk(j-1)+1
                else
                    exit
                end if
            end do
            prod=1
            do j=1, nir
                prod=prod*Xi(j)**BINk(j)
            end do
            i=sum(BINk)
            Ai(i)=Ai(i)+prod
        end do

        do i=0, nir
            if(mod(i,2)/=0) Ai(i)=-Ai(i)
            Ai(i)=Ai(i)*Raznosti(ni-nir+1,nir)/Factik(nir)/(h**(nir))
        end do
  

        do i=0, nir
            aai(i) = aai(i)+Ai(nir-i)
        end do 
   
    end do
    aai(0)=aai(0)+vectorY(ni+1)

    do i=0, ni
        write(*,*) "NUTONAi2", aai(i)
    end do

    do j=1, ni+1 
        summ=0
        do i=1, ni+1
            summ=summ+ aai(i-1)*vectorX(j)**(i-1)
        end do 
       Ypp(j)=summ
       write(*,*) vectorX(j), Ypp(j)
    end do



    call cpu_time(start_time)
        do i=1, 100000
            summ=0
            do j=0, ni
                summ = summ + aai(j)*(i)**(ni-j)
            end do
        end do
    call cpu_time(finish_time)
    write(*,*)"Time Usually 100000 = ", finish_time - start_time

    call cpu_time(start_time)
        do i=1, 1000000
            summ=0
            do j=0, ni
                summ = summ + aai(j)*(i)**(ni-j)
            end do
        end do
    call cpu_time(finish_time)
    write(*,*)"Time Usually 1000000 = ", finish_time - start_time
    write(*,*)
    call cpu_time(start_time)
        do i=1, 100000
            summ = aai(0)
            do j=1, ni
                summ = summ*i + aai(j)
            end do
        end do
    call cpu_time(finish_time)
    write(*,*)"Time Gorner 100000 = ", finish_time - start_time

    call cpu_time(start_time)
        do i=1, 1000000
            summ = aai(0)
            do j=1, ni
                summ = summ*i + aai(j)
            end do
        end do
    call cpu_time(finish_time)
    write(*,*)"Time Gorner 1000000 = ", finish_time - start_time

  
    open (unit = 2004, file = trim(resFileName))
    write(2004,*) "         X:               Y:            YP:"
    do i=1, ni+1
        write(2004,*)  vectorX(i), vectorY(i), Ypp(i)
    end do
    write(2004,*) "        dY:               ddY:                      "
    do i=1, ni+1
        write(2004,*)  abs(vectorY(i)- Ypp(i)), abs(vectorY(i)- Ypp(i))/Ypp(i)
    end do 

    close(2004)


end subroutine rezFile1

