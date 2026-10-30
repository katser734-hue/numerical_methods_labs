module module_main
  implicit none
  integer      ::  numLab, ni, i, j,ff, N, status , T, shagX, nch, k, ST=6, ST1=11, KNC=6, KNC1=11, KR=10
  integer      ::   nir
  real(4) :: diapozonX1, diapozonX2, h,  x, p, eps, zad
  real(8):: summ, prod, y, z, Epsik
  real(4) ::  pi=3.141592653688, a, b, Ar

  real, dimension(5) :: Pamat
  real(8), dimension(5) :: Pamat12

  real(4), allocatable:: vectorX(:)
  real(4), allocatable:: vectorY(:)

  !real, allocatable:: matrix(:,:)


  character(len = 200) :: line, resFileName, str, tmp_str
end module module_main