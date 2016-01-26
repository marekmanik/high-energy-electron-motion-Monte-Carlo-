module params
implicit none

integer, parameter :: DP=kind(1.0d0)


real(DP),parameter :: c_light = 299792458.0_dp
real(DP),parameter :: amu_const     = 1.660538921e-27_dp
real(DP),parameter :: q0_const      = 1.602176565e-19_dp
real(DP),parameter :: me_const      = 9.10938291e-31_dp
!real(DP),parameter :: mGas_const    = 2.0_dp*15.9994_dp*amu_const

real(DP),parameter :: mGas_N2    = 2.0_dp*14.007_dp*amu_const
real(DP),parameter :: mGas_O2    = 2.0_dp*15.999_dp*amu_const

!OPAL (Moss Table 4)
real(DP),parameter :: eps_bar_N2  = 13.0
real(DP),parameter :: eps_bar_O2  = 17.4


real(DP),parameter :: kb_const      = 1.3806488e-23_dp
real(DP),parameter :: Townsend_const= 1e-21_dp ! V m^2
real(DP),parameter :: ElectronVolt  = q0_const ! V m^2

real(DP),parameter :: PI_MATH = 4.0_DP*atan(1.0_DP)
real(DP),parameter :: TWO_PI = 2.0_DP*PI_MATH

integer, parameter :: ndim = 3 ! dimension of space

real(DP),parameter :: Ngas  = 2.688e25_dp  ! m^{-3}  Moss
real(DP),parameter :: Te_init = 0.5_dp*q0_const/kb_const  !

real(DP),parameter :: ppN2 = 0.8_dp
real(DP),parameter :: ppO2 = 0.2_dp


real(DP), parameter :: delta_MC = 0.1_dp

real(DP), parameter :: Ek_const = 32e5_dp ! Critical feild air, V/m

!real(DP),parameter :: E_over_N  = 100.0_dp*Townsend_const ! V m^2

!real(DP),parameter :: accel = q0_const*E_over_N*Ngas/me_const ! acceleration ms^{-2}

!real(DP), parameter :: sigma_e = 6.954e-20_dp ! m^{-2}  [model C]
real(DP), parameter :: sigma_e = 10.0e-20_dp ! m^{-2}   [modely A,B,D]
real(DP), parameter :: sigma_i = 10.0e-20_dp ! m^{-2}
real(DP), parameter :: eps_i = 0.516_dp*ElectronVolt
real(DP), parameter :: k_const = 0.4e-20_dp/ElectronVolt

!real(DP), parameter :: m_red = me_const*mGas_const/(mGas_const+me_const)
!real(DP), parameter :: m_red = me_const
!real(DP), parameter :: vel_th = sqrt(2.0_dp*eps_i/m_red) ! threshold velocity

!real(DP),parameter :: vmax = 20.0_dp*vel_th
real(DP),parameter :: T_gas = 300.0_dp

end module params
