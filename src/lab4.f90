subroutine LabNumber4
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB4_N") /=0 ) read (str,*) tmp_str, N
        if (index(str,"LAB4_ni") /=0 ) read (str,*) tmp_str, ni
        if (index(str,"LAB4_X1") /=0 ) read (str,*) tmp_str, diapozonX1
        if (index(str,"LAB4_X2") /=0 ) read (str,*) tmp_str, diapozonX2
    end do
    close (2000)
end subroutine LabNumber4




subroutine rezFile4
    use module_main
    implicit none
    real(8), dimension(20) :: ZnachX
    real(8), dimension(20) :: ZnachY
    real(8), dimension(19,6) :: RazY
    real(8):: Proiz11, Proiz12, Proiz21, Proiz22 !Производные Ньютона
    real(8):: SProiz11, SProiz12, SProiz21, SProiz22 !Производные Стирлинга
    real(8):: LProiz11, LProiz12, LProiz21, LProiz22 !Производные Лагранжа
    real(8):: TCN11, TCN12, TCN21, TCN22 


     
    integer  KK
    real(8):: chiiss

    integer, allocatable :: BINk(:)
    real(8), allocatable :: Ai(:), aai(:),chis(:), Xi(:), Ypp(:), AH(:)
    allocate(vectorX(ni+1))
    allocate(vectorY(ni+1))
    allocate(Ai(0:ni), aai(0:ni), BINk(ni), chis(ni+1), Xi(ni), Ypp(ni+1),AH(0:ni))

    h=0.01

    do i=0, 19 
        ZnachX(i+1)=0.01*i
    end do
    a=0.01*N
    do i=0, 19 
        ZnachY(i+1)=sqrt(ZnachX(i+1)+a)
    end do
    do i=0, 19 
        write(*,*) ZnachX(i+1), ZnachY(i+1)
    end do
    TCN11=1.0/2.0/sqrt(a+0.07)
    TCN12=1.0/2.0/sqrt(a+0.14)
    TCN21=-1.0/4.0/(sqrt(a+0.07)**3)
    TCN22=-1.0/4.0/(sqrt(a+0.14)**3)
    
    Proiz11=0
    Proiz12=0
    Proiz21=0
    Proiz22=0

    SProiz11=0
    SProiz12=0
    SProiz21=0
    SProiz22=0

    LProiz11=0
    LProiz12=0
    LProiz21=0
    LProiz22=0
    
    do i = 1, 19   
        do j=1, 6
            RazY(i,j)=0
        end do
    end do
 
    do i=1, 19
        RazY(i,1)=ZnachY(i+1)-ZnachY(i)
    end do

   write(*,*) 
   k=0
    do j=2, 6
         do i=1, 18-k
            RazY(i,j)=RazY(i+1,j-1)-RazY(i,j-1)
        end do
        k=k+1
    end do

    do i = 1, 19   
        write(*,*) RazY(i, :)
    end do
   write(*,*) 

     do i=1, 6
        Proiz11=Proiz11+(-1)**(i-1)*RazY(8,i)/i
        Proiz12=Proiz12+(-1)**(i-1)*RazY(15,i)/i
     end do
     Proiz21=RazY(8,2)-RazY(8,3)+11.0/12.0*RazY(8,4)-5.0/6.0*RazY(8,5)
     Proiz22=RazY(15,2)-RazY(15,3)+11.0/12.0*RazY(15,4)-5.0/6.0*RazY(15,5)
     Proiz11=Proiz11/h
     Proiz12=Proiz12/h
     Proiz21=Proiz21/h**2
     Proiz22=Proiz22/h**2

    SProiz11=(RazY(7,1)+RazY(8,1))/2.0-1.0/12*(RazY(6,3)+RazY(7,3)+1.0/60*(RazY(5,5)-RazY(6,5)))
    SProiz12=(RazY(14,1)+RazY(15,1))/2.0-1.0/12*(RazY(13,3)+RazY(17,3)+1.0/60*(RazY(12,5)-RazY(13,5)))
    SProiz21=RazY(7,2)-1.0/12*RazY(6,4)+1.0/90*RazY(4,6)
    SProiz22=RazY(14,2)-1.0/12*RazY(13,4)+1.0/90*RazY(11,6)

    SProiz11=SProiz11/h
    SProiz12=SProiz12/h
    SProiz21=SProiz21/h**2
    SProiz22=SProiz22/h**2

    Ai(:)=0; aai(:)=0; BINk(:)=0
    N=5
    h = (diapozonX2-diapozonX1)/ni
    do i=1, ni+1
        vectorX(i) = diapozonX1+h*(i-1)
        vectorY(i) = sqrt(a+vectorX(i)) 
        !vectorY(i) = exp(N*vectorX(i)) 
    end do
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


    do i=1, ni
        LProiz11=LProiz11 + aai(i-1)*(ni-i+1)*0.07**(ni-i)
        LProiz12=LProiz12 + aai(i-1)*(ni-i+1)*0.14**(ni-i)
    end do
    do i=1,ni-1
        LProiz21=LProiz21 + aai(i-1)*(ni-i+1)*(ni-i)*0.07**(ni-i-1)
        LProiz22=LProiz22 + aai(i-1)*(ni-i+1)*(ni-i)*0.14**(ni-i-1)
    end do 
   
 
        

 

write(*,*)" Proiz11, Proiz12, Proiz21, Proiz22: ", Proiz11, Proiz12, Proiz21, Proiz22
write(*,*)" SProiz11, SProiz12, SProiz21, SProiz22: ", SProiz11, SProiz12, SProiz21, SProiz22
write(*,*)" LProiz11, LProiz12, LProiz21, LProiz22:", LProiz11, LProiz12, LProiz21, LProiz22
write(*,*) "TCN11, TCN12, TCN21, TCN22", TCN11, TCN12, TCN21, TCN22 
write(*,*)

write(*,*)" dAbsN ", abs(Proiz11-TCN11), abs(Proiz12-TCN12), abs(Proiz21-TCN21),abs(Proiz22-TCN22)
write(*,*)" dAbsS ",abs(SProiz11-TCN11), abs(SProiz12-TCN12), abs(SProiz21-TCN21), abs(SProiz22-TCN22)
write(*,*)" dAbsL", abs(LProiz11-TCN11), abs(LProiz12-TCN12), abs(LProiz21-TCN21), abs(LProiz22-TCN22)
write(*,*)

write(*,*)" dOTNuton:", abs(Proiz11-TCN11)/(TCN11), abs(Proiz12-TCN12)/(TCN12)
write(*,*)" dOTNStir:", abs(SProiz11-TCN11)/(TCN11), abs(SProiz12-TCN12)/(TCN12)
write(*,*)" dOTNLagr:", abs(LProiz11-TCN11)/(TCN11), abs(LProiz12-TCN12)/(TCN12)

write(*,*)" dOTNuton:", abs(Proiz21-TCN21)/(TCN21),abs(Proiz22-TCN22)/(TCN22)
write(*,*)" dOTNStir:", abs(SProiz21-TCN21)/(TCN21), abs(SProiz22-TCN22)/(TCN22)
write(*,*)" dOTNLagr:",  abs(LProiz21-TCN21)/(TCN21), abs(LProiz22-TCN22)/(TCN22)
end subroutine rezFile4

