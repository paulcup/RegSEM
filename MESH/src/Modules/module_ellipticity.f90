!--------------------------------------------------------------
module module_ellipticity
!--------------------------------------------------------------
  implicit none
  public:: init_ellipticity,get_ellipticity
  private
  integer, parameter :: NBRAD=100, CST=7 !integration  degree
  doubleprecision :: rmax
  doubleprecision, parameter :: Omega=2.d0*3.141592653589793d0/24.d0/3600.d0, bigg=6.6723d-11
  doubleprecision, dimension(NBRAD) :: rad,eps,epsp
  doubleprecision, dimension(CST) :: xgll,wgll
  
contains
!-------------------------------------------------------------    
  subroutine init_ellipticity()
!using premc density for sake of simplicity.
!Dahlen & Tromp pp 600
!-------------------------------------------------------------    
    use funaro
    use module_spline
!
    implicit none
    integer :: i,j,nbi,NBINT
    doubleprecision :: bid,Inertia,PI,Masse,epsilon_a,I1,I2,tmp,dr,ypn,yp1
    doubleprecision, dimension(NBRAD) :: nu
    doubleprecision, dimension(:), allocatable ::  ri,radintgr,weight,rho
!
!    
    call def_xgll(xgll,CST)
    do i=1,CST
       xgll(i)=(xgll(i)+1.d0)/2.d0
    enddo
    call def_wgll(wgll,CST)
    wgll(:)=wgll(:)/2.d0
!
!getting number of layers
    call premc(0.5d0,bid,bid,bid,.false.,nbi=nbi)
    allocate(ri(nbi),radintgr(CST*nbi),weight(CST*nbi),rho(CST*nbi))
    call premc(0.5d0,bid,bid,bid,.false.,ri=ri)    
    ri(:)=ri(:)*1000.d0    
    rmax =ri(nbi)
    dr = rmax/(NBRAD-1)
    do i=1,NBRAD
       rad(i)=(i-1)*dr
    enddo
    call get_radintgr(rmax,ri,nbi,radintgr,weight,NBINT)
    rho(:)=0.d0
    do i=1,NBINT
       call premc(radintgr(i),rho(i),bid,bid,.false.)
    enddo
    rho(:)=rho(:)*1000.d0
    PI=4.d0*atan(1.d0)
    Masse  =4.d0     *PI*SUM(rho(:)*(radintgr(:))**2*weight(:))*rmax**2
    Inertia=8.d0/3.d0*PI*SUM(rho(:)*(radintgr(:))**4*weight(:))*rmax**4
    epsilon_a=10.d0*Omega**2*rmax**3/(bigg*Masse*(4.d0+25.d0*(1.d0  &
               -1.5d0*Inertia/Masse/rmax**2)**2))
    do i=2,NBRAD
       rho(:)=0.d0
       call get_radintgr(rad(i),ri,nbi,radintgr,weight,NBINT)
       do j=1,NBINT
          call premc(radintgr(j),rho(j),bid,bid,.false.)
       enddo
       rho(:)=rho(:)*1000.d0       
       I1=SUM(rho(:)*(radintgr(:))**4*weight(:))*rmax**4
       I2=SUM(rho(:)*(radintgr(:))**2*weight(:))*rmax**2*rad(i)**2
       nu(i) =25.d0/4.d0*(1.d0 - I1/I2)**2-1.d0
    enddo
    eps(NBRAD)=epsilon_a
    tmp=0.d0
    do i=NBRAD-1,2,-1
       tmp=tmp+nu(i+1)/rad(i+1)*dr
       eps(i)=epsilon_a*exp(-tmp)
    enddo
!    
! preparing spline:    
!  
    ypn=(eps(NBRAD)-eps(NBRAD-1))/dr
    yp1=(eps(2    )-eps(1      ))/dr
    call spline(rad,eps,yp1,ypn,epsp)    
!-------------------------------------------------------------    
  end subroutine init_ellipticity
!-------------------------------------------------------------    
!-------------------------------------------------------------    
  doubleprecision function get_ellipticity(r)
!-------------------------------------------------------------    
    use module_spline
    implicit none
    doubleprecision :: r
    get_ellipticity=splint(rad,eps,epsp,r)
!-------------------------------------------------------------    
  end function get_ellipticity
!-------------------------------------------------------------    

!-------------------------------------------------------------    
  subroutine get_radintgr(rm,ri,nbi,ra,weight,nb)
!-------------------------------------------------------------    
    implicit none
    integer, intent(in) :: nbi
    doubleprecision, intent(in) :: rm
    doubleprecision, dimension(:), intent(in) :: ri    
    integer, intent(out) :: nb
    doubleprecision, dimension(:), intent(out) :: ra,weight
    integer :: i,j
    doubleprecision :: h
    doubleprecision, dimension(nbi) :: ritmp
    
!
!    looking for the last interface below rm
    i=1
    ritmp(1)=ri(1)
    do while(rm>ri(i).and.i<nbi-1)       
       i=i+1
       ritmp(i)=ri(i)
    enddo
    i=i-1
    ritmp(i+1)=rm
    weight(:)=0.0d0
    ra    (:)=0.0d0
    nb=CST*i
    do i=1,nb/CST
       h=ritmp(i+1)-ritmp(i)-1.d-5
       do j=1,CST
          ra((i-1)*CST+j)  =ritmp(i)+xgll(j)*h
          weight((i-1)*CST+j)=h*wgll(j)
       enddo
    enddo
    ra(:)=ra(:)/rmax
!-------------------------------------------------------------    
  end subroutine get_radintgr
!-------------------------------------------------------------    
  subroutine premc(x0,ro,vp,vs,flag_g,GR,Qmu,nbi,ri)
!---------------------------------------
!
!   Given non-dimensionalized radius x0, prem returns
!   density, ro, compressional velocity, vp, and shear 
!   velocity, vs according to the PREM model (1s ref period)
!     (Anderson and Dziewonski, 1981).  
!
!   Also returns gravity gr
!
!
      integer 				:: i
      doubleprecision    		:: x0,x,ro,vp,vs
      logical				:: flag_g
      doubleprecision,optional   	:: GR,Qmu
      integer, optional :: nbi
      doubleprecision,optional,dimension(12) :: ri
      doubleprecision,parameter   	:: Rt = 6371000.d0,epsilon = 1.d-9
      doubleprecision,parameter		:: bigg = 6.6723d-11
      doubleprecision,dimension(12)	:: r,q
      doubleprecision,dimension(11,4)	:: d,p,s
      doubleprecision,dimension(11)	:: cumul
      logical				:: pastrouve
!
      doubleprecision			:: CST,r1,r2,t1,t2,t3,t4
      integer				:: j
!
   ro = 0.0d0; vp = 0.d0; vs = 0.d0;

!
   r(1) = 0.d0;     r(2) = 1221.5d0; r(3) = 3480.d0; r(4) = 3630.d0
   r(5) = 5600.d0;  r(6) = 5701.d0;  r(7) = 5771.d0; r(8) = 5971.d0
   r(9) = 6151.d0; r(10) = 6291.d0; r(11) = 6346.6d0
   r(12) = 6371.d0
!
  if (present(nbi)) then
     nbi=12
  endif
  if (present(ri)) then
     ri(:)=r(:)
  endif
!
  q(1)=84.6d0; q(2)=1.d5; q(3)=312.d0; q(4)=312.d0;q(5)=312.d0;
  q(6)=312.d0; q(7)=156.d0; q(8)=156.d0;q(9)=70.d0;q(10)=191.d0;
  q(11)=300.d0
!
  d(:,:) = 0.d0
!
  d(1,1) = 13.0885d0; 			  d(1,3) = -8.8381d0 
  d(2,1) = 12.5815d0; d(2,2) = -1.2638d0; d(2,3) = -3.6426d0; d(2,4) = -5.5281d0
  d(3,1) = 7.9565d0 ; d(3,2) = -6.4761;   d(3,3) = 5.5283d0;  d(3,4) = -3.0807d0
  d(4,1) = 7.9565d0 ; d(4,2) = -6.4761;   d(4,3) = 5.5283d0;  d(4,4) = -3.0807d0
  d(5,1) = 7.9565d0 ; d(5,2) = -6.4761;   d(5,3) = 5.5283d0;  d(5,4) = -3.0807d0
  d(6,1) = 5.3197d0 ; d(6,2) = -1.4836d0
  d(7,1) = 11.2494d0; d(7,2) = -8.0298d0
  d(8,1) = 7.1089d0 ; d(8,2) = -3.8045d0
  d(9,1) = 2.6910d0 ; d(9,2) = 0.6924d0
  d(10,1) = 2.6910d0; d(10,2) = 0.6924d0
  d(11,1) = 2.7d0  
!  d(11,1) = 2.715d0  
!
!----
!
  p(:,:) = 0.d0
!
  p(1,1) = 11.2622d0 ; p(1,3) = -6.3640d0
  p(2,1) = 11.0487d0 ; p(2,2) = -4.0362d0; p(2,3)  = 4.8023d0; p(2,4) = -13.5732d0
  p(3,1) = 15.3891d0 ; p(3,2) = -5.3181d0; p(3,3)  = 5.5242d0; p(3,4) = -2.5514d0
  p(4,1) = 24.952d0 ; p(4,2)  = -40.4673d0; p(4,3) = 51.4832d0; p(4,4) = -26.6419d0
  p(5,1) = 29.2766d0 ; p(5,2) = -23.6027d0; p(5,3) = 5.5242d0; p(5,4) = -2.5514d0
  p(6,1) = 19.0957d0 ; p(6,2)  = -9.8672d0
  p(7,1) = 39.7027d0 ; p(7,2)  = -32.6166d0
  p(8,1) = 20.3926d0 ; p(8,2)  = -12.2569d0
  p(9,1) = 4.1875d0 ; p(9,2)  = 3.9382d0
  p(10,1) = 4.1875d0 ; p(10,2) = 3.9382d0
  p(11,1) = 6.15d0 
!  p(11,1) = 6.18d0 
!
!----
!
  s(:,:) = 0.d0
!
  s(1,1) = 3.6678d0; s(1,3) = -4.4475d0

  s(3,1) = 6.9254d0; s(3,2) = 1.4672d0; s(3,3) = -2.0834d0; s(3,4) = 0.9783d0
  s(4,1) = 11.1671d0; s(4,2) = -13.7818d0; s(4,3) = 17.4575d0; s(4,4) = -9.2777d0
  s(5,1) = 22.3459d0; s(5,2) = -17.2473d0; s(5,3) = -2.0834d0; s(5,4) = 0.9783d0
  s(6,1) = 9.9839; s(6,2) = -4.9324
  s(7,1) = 22.3512d0; s(7,2) = -18.5856d0 
  s(8,1) = 8.9496d0; s(8,2) = -4.4597
  s(9,1) = 2.1519d0; s(9,2) = 2.3481d0
  s(10,1) = 2.1519d0; s(10,2) = 2.3481d0
  s(11,1) = 3.32d0 
!  s(11,1) = 3.47d0 
!
!
  r(:) = r(:) * 1000.d0
  x    = x0   * Rt
!
  pastrouve = .true.
!
 do i = 1,11
 if ( pastrouve.and.(x >= r(12)) ) then
  x0 = 1.d0
  ro = d(11,1) + x0*( d(11,2) + x0*( d(11,3) + x0*( d(11,4) )))
  vp = p(11,1) + x0*( p(11,2) + x0*( p(11,3) + x0*( p(11,4) )))
  vs = s(11,1) + x0*( s(11,2) + x0*( s(11,3) + x0*( s(11,4) )))
  if (present(Qmu))  Qmu=q(11)
  pastrouve = .false.
 endif

 if ( pastrouve.and.((x >= r(i)).and.(x < r(i+1))) ) then
  ro = d(i,1) + x0*( d(i,2) + x0*( d(i,3) + x0*( d(i,4) )))
  vp = p(i,1) + x0*( p(i,2) + x0*( p(i,3) + x0*( p(i,4) )))
  vs = s(i,1) + x0*( s(i,2) + x0*( s(i,3) + x0*( s(i,4) )))
  if (present(Qmu)) Qmu=q(i)
  pastrouve = .false.
 endif

 enddo
  if (.not. flag_g) return

!-- Gravity

	CST = 16.d0*datan(1.d0)*bigg

	do i = 1,11
		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = r(i+1)
		r1 = r(i)

		cumul(i) =  t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6) - &
			  ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
	
	enddo

	if ( x > r(12) ) x = r(12)

	do i =1,11
 	if ((x > r(i)).and.(x <= r(i+1))) then

		GR = 0.d0
		do j = 1,i-1
			GR = GR + cumul(j)
		enddo

		t1 = d(i,1)/3.d0
		t2 = d(i,2)/(Rt*4.d0)
		t3 = d(i,3)/((Rt**2)*5.d0)		
		t4 = d(i,4)/((Rt**3)*6.d0)
		r2 = x
		r1 = r(i)

		GR = GR +   t1*(r2**3) + t2*(r2**4) + t3*(r2**5) + t4*(r2**6)	&
			- ( t1*(r1**3) + t2*(r1**4) + t3*(r1**5) + t4*(r1**6) )
		GR = GR * CST
		if (x > epsilon) GR = GR/(x**2)
		return
	endif
	enddo
!
  end subroutine premc
!--------------------------------------------------------------
end module module_ellipticity 
!--------------------------------------------------------------
