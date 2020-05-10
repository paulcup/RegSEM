module ssources


type :: Source

integer :: i_type_source, i_time_function, elem, proc
integer, dimension(0:2) :: gll
doubleprecision :: tau_b, depth,realcolat,reallong,refcolat,reflong
doubleprecision, dimension(0:3) :: fh
doubleprecision, dimension(0:2) :: refcoord, i_dir
doubleprecision, dimension(:), pointer :: timefunc
doubleprecision, dimension(0:2,0:2) :: Moment, InvGrad
doubleprecision, dimension (:,:,:,:), pointer :: coeff

end type 


contains

! ###############################################
doubleprecision function CompSource (Sour,time)

type (source) :: Sour
doubleprecision :: time

select case (Sour%i_time_function)
 case (1) 
     CompSource = Gaussian (time,Sour%tau_b)
 case (2)
     CompSource = Ricker (time,Sour%tau_b)
 case (4)
     CompSource = Error (time,Sour%tau_b)
 case (5)
     CompSource = Ricker2 (time,Sour%tau_b)
 case (6)
     CompSource = Hz_contact (time,Sour%tau_b)
end select

return
end function

! ################################################
doubleprecision function Gaussian (time,tau)

use angles

doubleprecision :: tau,time

doubleprecision :: beta_gauss
doubleprecision, parameter :: alpha_gauss = 1.628

beta_gauss = 2.d0*alpha_gauss/tau
Gaussian = beta_gauss * dexp(-(beta_gauss*(time-tau))**2) / dsqrt(pi)

return
end function
 
! ################################################
doubleprecision function Ricker (time,tau)

use angles

doubleprecision :: time,tau

doubleprecision :: beta_gauss
doubleprecision, parameter :: alpha_gauss = 1.628

beta_gauss = 2.d0*alpha_gauss/tau
Ricker = -2.d0*beta_gauss**3*(time-tau) * dexp(-(beta_gauss*(time-tau))**2) / dsqrt(pi)

return
end function

! ################################################
doubleprecision function Error (time,tau)

doubleprecision :: time,tau

doubleprecision :: cste,beta_erf
doubleprecision, parameter :: alpha_erf = 1.628

beta_erf = 2.d0*alpha_erf/tau
    cste = derf(-2.d0*alpha_erf)/2.d0
Error = derf(beta_erf*(time-tau))/2.d0 - cste

return
end function

! ################################################
doubleprecision function Ricker2 (time,tau)

use angles

doubleprecision :: time,tau

doubleprecision :: alpha,T

T = tau/pi
alpha = (time-tau)/T
Ricker2 = (1.d0 - 2.d0*alpha**2) * dexp(-alpha**2)

return
end function

! ################################################
doubleprecision function Hz_contact (time,tau)

doubleprecision :: time,tau

doubleprecision :: t,tmp

if (time<2.d0*tau) then
   t = ((time-tau)/tau)**2
   tmp = dsqrt(1.d0 - 4.d0*t/3.d0 + (t**2)/3.d0)
   Hz_contact = tmp**3
else
   Hz_contact = 0.d0
endif

return
end function

! #################################################
end module ssources
