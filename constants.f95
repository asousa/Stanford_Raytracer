module constants

real*8, parameter :: EPS0 = 8.854187817e-12_8
real*8, parameter :: PI = 3.141592653589793238462643_8
real*8, parameter :: MU0 = PI * 4e-7_8
real*8, parameter :: C = sqrt(1/EPS0/MU0)
real*8, parameter :: R_E = 6371.2e3_8

end module constants
