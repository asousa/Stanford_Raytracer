module constants
use types

real(kind=DP), parameter :: EPS0 = 8.854187817e-12_DP
real(kind=DP), parameter :: PI = 3.141592653589793238462643_DP
real(kind=DP), parameter :: MU0 = PI * 4e-7_DP
real(kind=DP), parameter :: C = sqrt(1.0_DP/EPS0/MU0)
real(kind=DP), parameter :: R_E = 6371.2e3_DP
real(kind=DP), parameter :: VERSION = 1.04_DP

end module constants
