subroutine set_attenuation_param (Tdomain)


use sdomains
use angles

implicit none

type(domain),target, intent (INOUT) :: Tdomain

integer :: n, i,j,k, n_solid, i_count, ngllx,nglly,ngllz
doubleprecision :: factor_scale_mu0, w_c_source, a_val,b_val, big_omega, factor_scale, f0_source, dt, &
                   factor_scale_mu, one_minus_sum_beta, Q_mu, T1_attenuation,T2_attenuation
doubleprecision, dimension(:), allocatable :: tau_mu,tau_sigma, beta, tauinv

doubleprecision, parameter :: f0_modele = 1.d0
!doubleprecision, parameter :: f0_modele = 1000.d0


T1_attenuation = Tdomain%T1_att;   T2_attenuation = Tdomain%T2_att
f0_source = exp((log(1.d0/T1_attenuation)+log(1.d0/T2_attenuation))/2.d0)
w_c_source = 2.d0*PI*f0_source

n_solid = Tdomain%n_sls
allocate (tau_mu(0:n_solid-1))
allocate (tau_sigma(0:n_solid-1))
allocate (beta(0:n_solid-1))

do n = 0,Tdomain%n_elem-1
 if (Tdomain%specel(n)%PML .eqv. .false.) then

    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    do i = 0,ngllx-1
     do j = 0,nglly-1
      do k = 0,ngllz-1

          Q_mu = Tdomain%specel(n)%Q(i,j,k)
          call compute_constant_Q(T1_attenuation,T2_attenuation,n_solid,Q_mu,tau_mu,tau_sigma)

          beta(:) = 1.d0 - tau_mu(:) / tau_sigma(:)

          factor_scale_mu0 = 1.d0 + 2.d0 * log(f0_source/f0_modele) / (PI*Q_mu)
          a_val = 1.d0
          b_val = 0.d0
          one_minus_sum_beta = 1.d0

          do i_count = 0,n_solid-1
              a_val = a_val - w_c_source * w_c_source * tau_mu(i_count) * &
                      (tau_mu(i_count) - tau_sigma(i_count)) / &
                      (1.d0 + w_c_source*w_c_source*tau_mu(i_count)*tau_mu(i_count))
              b_val = b_val + w_c_source * (tau_mu(i_count) - tau_sigma(i_count)) / &
                      (1.d0 + w_c_source*w_c_source*tau_mu(i_count)*tau_mu(i_count))
              one_minus_sum_beta = one_minus_sum_beta - beta(i_count)
          enddo

          Tdomain%specel(n)%onemSbeta(i,j,k) = one_minus_sum_beta

          big_omega = a_val*(dsqrt(1.d0 + b_val*b_val/(a_val*a_val))-1.d0);

          !--- quantity by which to scale mu to get mu_relaxed
          factor_scale_mu = b_val * b_val / (2.d0 * big_omega)

          !--- total factor by which to scale mu0
          factor_scale = factor_scale_mu * factor_scale_mu0

          !--- check that the correction factor is close to one
          if (factor_scale < 0.9d0 .or. factor_scale > 1.1d0) then
              stop 'incorrect correction factor in attenuation model'
          endif

          ! on corrige mu, mais avant on remplace lambda par kappa
          ! (ATTENTION a ne pas se planter dans le suite) :
          Tdomain%specel(n)%Lambda(i,j,k) = Tdomain%specel(n)%Lambda(i,j,k) + &
                                            2.d0/3.d0 * Tdomain%specel(n)%Mu(i,j,k)
          Tdomain%specel(n)%Mu(i,j,k) = factor_scale * Tdomain%specel(n)%Mu(i,j,k)

          !on initialise les coefs Runge-Kutta:
          dt = Tdomain%sTimeParam%dt
          allocate (tauinv(0:n_solid-1))
          tauinv(:) = - 1.d0 / tau_sigma(:)
          Tdomain%specel(n)%factor_common_3(:,i,j,k) = 2.d0 * beta(:) * tauinv(:)
          Tdomain%specel(n)%alphaval_3(:,i,j,k) = 1.d0 + dt*tauinv(:) + dt**2*tauinv(:)**2 / 2.d0 + &
                                                  dt**3*tauinv(:)**3 / 6.d0 + dt**4*tauinv(:)**4 / 24.d0
          Tdomain%specel(n)%betaval_3(:,i,j,k) = dt / 2.d0 + dt**2*tauinv(:) / 3.d0 + &
                                                 dt**3*tauinv(:)**2 / 8.d0 + dt**4*tauinv(:)**3 / 24.d0
          Tdomain%specel(n)%gammaval_3(:,i,j,k) = dt / 2.d0 + dt**2*tauinv(:) / 6.d0 + &
                                                  dt**3*tauinv(:)**2 / 24.d0
          deallocate (tauinv)
      enddo
     enddo
    enddo

 endif
enddo

deallocate (tau_mu, tau_sigma, beta)


return
end subroutine set_attenuation_param
