subroutine LabNumber18
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB18_N") /=0 ) read (str,*) tmp_str, N

    end do
    close (2000)
end subroutine LabNumber18



subroutine rezFile18
    use module_main
    implicit none
    real iks(101), igr(101), igth(101)
    real (16) :: qi(101), pii(101), fii(101)
    real (16) :: mi(101), Fi(101)
    real (16) :: ci(101), di(101), nii(101)



    real  :: iks2(1001), igr2(1001), igth2(1001)
    
    real (16) :: qi2(1001), pii2(1001), fii2(1001)
    real (16) :: mi2(1001), Fi2(1001)
    real (16) :: ci2(1001), di2(1001), nii2(1001)

    real (16) :: alfa, aH, al0, al1, be0, be1, AA, BB, aH2
    alfa = 1+0.001*N

    aH = 0.01
    aH2 = 0.001

    be0 = 2
    be1 = -1
    al0 = 1
    al1 = 1
    AA = alfa*(sin(1.0)+cos(1.0))+6+cos(1.0)/2
    BB = alfa*(sin(2.0)+2*cos(2.0))+5.0/8*sin(2.0)-cos(2.0)/4-3.0/4
      

    do i=1, 101
        iks(i) = (i-1)*aH + 1
        fii(i)=-alfa*(cos(iks(i))-sin(iks(i))/iks(i)+8*cos(iks(i))/iks(i)**2)&
        -(iks(i)**2*sin(iks(i))+3*iks(i)*cos(iks(i))+5*sin(iks(i)))/2/iks(i)**3
        pii(i)=-1.0/iks(i)
        qi(i)=-8.0/iks(i)**2
        mi(i)=(-4+2*aH**2*qi(i))/(2+aH*pii(i))
        nii(i)=(2-aH*pii(i))/(2+aH*pii(i))
        Fi(i)=2*fii(i)/(2+aH*pii(i))
        igth(i) = sin(iks(i))/2/iks(i)+alfa*cos(iks(i))+(2*alfa*sin(1.0)+5)/5*iks(i)**4-1/iks(i)**2
    end do   

    ci(1) = ((mi(2)+4)*al1)/(2*aH*al0+al1*(nii(2)-3))
    di(1) = (aH*(2*AA+al1*Fi(2)*aH))/((mi(2)+4)*al1)

    do i=2, 101
        ci(i)=1/(mi(i)-nii(i)*ci(i-1))
        di(i)=aH**2*Fi(i)-nii(i)*ci(i-1)*di(i-1)
    end do



    igr(101) = (2*BB*aH-(4+ci(99))*ci(100)*di(100)+ci(99)*di(99))/(4*aH-3-4*ci(100)-ci(99)*ci(100))


    do i=100, 1, -1
        igr(i) = ci(i)*(di(i)-igr(i+1))
    end do

    do i=1, 101, 10
        write(*, *) iks(i),'&', igth(i), '&', igr(i), '&', abs(igth(i)-igr(i)), '&', abs(igth(i)-igr(i))/abs(igth(i)), '\\ \hline'
    end do
    write(*,*) 
    write(*,*) 
    



    do i=1, 1001
        iks2(i) = (i-1)*aH2 + 1
        fii2(i)=-alfa*(cos(iks2(i))-sin(iks2(i))/iks2(i)+8*cos(iks2(i))/iks2(i)**2)&
        -(iks2(i)**2*sin(iks2(i))+3*iks2(i)*cos(iks2(i))+5*sin(iks2(i)))/2/iks2(i)**3
        pii2(i)=-1.0/iks2(i)
        qi2(i)=-8.0/iks2(i)**2
        mi2(i)=(-4+2*aH2**2*qi2(i))/(2+aH2*pii2(i))
        nii2(i)=(2-aH2*pii2(i))/(2+aH2*pii2(i))
        Fi2(i)=2*fii2(i)/(2+aH2*pii2(i))
        igth2(i) = sin(iks2(i))/2/iks2(i)+alfa*cos(iks2(i))+(2*alfa*sin(1.0)+5)/5*iks2(i)**4-1/iks2(i)**2
    end do   

    ci2(1) = ((mi2(2)+4)*al1)/(2*aH2*al0+al1*(nii2(2)-3))
    di2(1) = (aH2*(2*AA+al1*Fi2(2)*aH2))/((mi2(2)+4)*al1)

    do i=2, 1001
        ci2(i)=1/(mi2(i)-nii2(i)*ci2(i-1))
        di2(i)=aH2**2*Fi2(i)-nii2(i)*ci2(i-1)*di2(i-1)
    end do



    igr2(1001) = (2*BB*aH2-(4+ci2(999))*ci2(1000)*di2(1000)+ci2(999)*di2(999))/(4*aH2-3-4*ci2(1000)-ci2(999)*ci2(1000))


    do i=1000, 1, -1
        igr2(i) = ci2(i)*(di2(i)-igr2(i+1))
    end do


    j=1
    do i=1, 1001, 100

        write(*,*) iks2(i),'&', igth2(i), '&', igr2(i), '&', abs(igth2(i)-igr2(i)),&
         '&', abs(igth2(i)-igr2(i))/abs(igth2(i)) , '&', abs(igr2(i)-igr(j)), '\\ \hline'
        j=j+10
    end do
    write(*,*) 
    write(*,*) 
    write(*,*) 
    write(*,*) 

    write(*,*) 
    write(*,*)







end subroutine rezFile18