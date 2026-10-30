subroutine LabNumber7
    use module_main
    implicit none
    open (unit = 2000, file = 'labnumbers.txt')
    do while(.true.)
        read (2000,'(a)', iostat=status) str
        if(status/=0) exit
        if (index(str,"LAB7_N") /=0 ) read (str,*) tmp_str, N
        if (index(str,"LAB7_ni") /=0 ) read (str,*) tmp_str, ni
    end do
    close (2000)
end subroutine LabNumber7


subroutine rezFile7
    use module_main
    implicit none
    real(8) :: Legan, Chebb, Lagg, Ermit
    real(8) :: THLegan, THChebb, THLagg, THErmit
    real,parameter :: nni=7

    real, dimension(ni) :: ALegandr
    real, dimension(ni) :: TLegandr
    real, dimension(ni) :: TCheb
    real, dimension(ni) :: ALagger
    real, dimension(ni) :: TLagger
    real, dimension(ni) :: AErmit
    real, dimension(ni) :: TErmit

    THLegan = 0.22860146
    THChebb = 0.53943106
    THLagg = 5.6049912
    THErmit = 0.35286800


    TLegandr(7)=0.949107912;     TLegandr(1)=-TLegandr(7)
    TLegandr(6)=0.741531186;     TLegandr(2)=-TLegandr(6)
    TLegandr(5)=0.405845151;     TLegandr(3)=-TLegandr(5)
    TLegandr(4)=0
    ALegandr(1)=0.129484966;     ALegandr(7)=ALegandr(1)
    ALegandr(2)=0.279705391;     ALegandr(6)=ALegandr(2)
    ALegandr(3)=0.381830051;     ALegandr(5)=ALegandr(3)
    ALegandr(4)=0.41795918

    Tcheb(1) = -0.8662468181;    TCheb(2) = -0.5296567753
    TCheb(3) = -0.3239118105;    TCheb(4) = 0
    TCheb(5) = 0.3239118105 ;     TCheb(6) = 0.5296567753
    Tcheb(7) = 0.8662468181

    ALagger(1)=0.4093189517;    ALagger(2)=0.4218312779
    ALagger(3)=0.1471263487;    ALagger(4)=0.0206335145
    ALagger(5)=0.0010740101;    ALagger(6)=0.0000158655
    ALagger(7)=0.0000000317
    TLagger(1)=0.1930436766;    TLagger(2)=1.0266648953
    TLagger(3)=2.5678767450;    TLagger(4)=4.9003530845
    TLagger(5)=8.1821534446;    TLagger(6)=12.7341802918
    TLagger(7)=19.3957278623

    AErmit(4)=0.8102646176;
    AErmit(3)=0.4256072526;   AErmit(2)=0.0545155828
    AErmit(1)=0.00097178125;  AErmit(7)=AErmit(1)
    AErmit(6)=AErmit(2);      AErmit(5)=AErmit(3)

    TErmit(1)=-2.6519613568;   TErmit(2)=-1.6735516288
    TErmit(3)=-0.8162878829;   TErmit(4)=0;
    TErmit(7)=-TErmit(1);      TErmit(6)=-TErmit(2)
    TErmit(5)=-TErmit(3)


    b=1+0.01*N
    a=1-0.02*N
    Legan = 0
    Chebb = 0
    Lagg  = 0
    Ermit = 0

    do i=1, ni
        Legan=Legan+(2.0/(sqrt((3+TLegandr(i))**3)*sqrt(12+4*TLegandr(i)+2*b)))*ALegandr(i)
    end do
 
 
    do i=1, ni
        Chebb=Chebb+1.0/(2*sqrt(1.5+0.5*TCheb(i)-0.01*N)*sqrt(2.5+0.5*TCheb(i)))
    end do
    Chebb=Chebb*2/ni

    
    do i=1, ni
        Lagg=Lagg+(exp(a*TLagger(i))/sqrt(TLagger(i)))*ALagger(i)
    end do

    a=-0.5+0.001*N
    do i=1, ni
        Ermit=Ermit+(1.0/sqrt(pi)*cos((2*a*TErmit(i))/(a+0.01)))*AErmit(i)
    end do


    write(*,*) "Legan",  THLegan,  Legan, abs(THLegan-Legan), abs(THLegan-Legan)/THLegan
    write(*,*) "Chebb",  THChebb,  Chebb, abs(THChebb-Chebb), abs(THChebb-Chebb)/THChebb
    write(*,*) "Lagg" ,  THLagg ,  Lagg , abs(THLagg-Lagg)  , abs(THLagg-Lagg)/THLagg   
    write(*,*) "Ermit",  THErmit,  Ermit, abs(THErmit-Ermit), abs(THErmit-Ermit)/THErmit
end subroutine rezFile7
