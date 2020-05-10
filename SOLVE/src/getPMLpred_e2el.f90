subroutine get_PMLprediction_e2el (Tdomain, n, bega, dt, i_simu)


use sdomains

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: n, i_simu
doubleprecision, intent(IN) :: dt, bega

integer :: ne, ngllx, nglly, ngllz 


ngllx = Tdomain%specel(n)%ngllx;  nglly = Tdomain%specel(n)%nglly;  ngllz = Tdomain%specel(n)%ngllz

ne = Tdomain%specel(n)%near_edges(0)
Tdomain%specel(n)%sSimu(i_simu)%Forces(1:ngllx-2,0,0,:) = Tdomain%sEdge(ne)%sSimu(i_simu)%Veloc(:,:) + &
                                                          dt * (0.5-bega) * Tdomain%sEdge(ne)%sSimu(i_simu)%Accel(:,:)

ne = Tdomain%specel(n)%near_edges(1)
Tdomain%specel(n)%sSimu(i_simu)%Forces(ngllx-1,1:nglly-2,0,:) = Tdomain%sEdge(ne)%sSimu(i_simu)%Veloc(:,:) + &
                                                                dt * (0.5-bega) * Tdomain%sEdge(ne)%sSimu(i_simu)%Accel(:,:)

ne = Tdomain%specel(n)%near_edges(2)
Tdomain%specel(n)%sSimu(i_simu)%Forces(1:ngllx-2,nglly-1,0,:) = Tdomain%sEdge(ne)%sSimu(i_simu)%Veloc(:,:) + &
                                                                dt * (0.5-bega) * Tdomain%sEdge(ne)%sSimu(i_simu)%Accel(:,:)

ne = Tdomain%specel(n)%near_edges(3)
Tdomain%specel(n)%sSimu(i_simu)%Forces(0,1:nglly-2,0,:) = Tdomain%sEdge(ne)%sSimu(i_simu)%Veloc(:,:) + &
                                                          dt * (0.5-bega) * Tdomain%sEdge(ne)%sSimu(i_simu)%Accel(:,:)

ne = Tdomain%specel(n)%near_edges(4)
Tdomain%specel(n)%sSimu(i_simu)%Forces(ngllx-1,0,1:ngllz-2,:) = Tdomain%sEdge(ne)%sSimu(i_simu)%Veloc(:,:) + &
                                                                dt * (0.5-bega) * Tdomain%sEdge(ne)%sSimu(i_simu)%Accel(:,:)

ne = Tdomain%specel(n)%near_edges(5)
Tdomain%specel(n)%sSimu(i_simu)%Forces(1:ngllx-2,0,ngllz-1,:) = Tdomain%sEdge(ne)%sSimu(i_simu)%Veloc(:,:) + &
                                                                dt * (0.5-bega) * Tdomain%sEdge(ne)%sSimu(i_simu)%Accel(:,:)

ne = Tdomain%specel(n)%near_edges(6)
Tdomain%specel(n)%sSimu(i_simu)%Forces(0,0,1:ngllz-2,:) = Tdomain%sEdge(ne)%sSimu(i_simu)%Veloc(:,:) + &
                                                          dt * (0.5-bega) * Tdomain%sEdge(ne)%sSimu(i_simu)%Accel(:,:)

ne = Tdomain%specel(n)%near_edges(7)
Tdomain%specel(n)%sSimu(i_simu)%Forces(ngllx-1,nglly-1,1:ngllz-2,:) = Tdomain%sEdge(ne)%sSimu(i_simu)%Veloc(:,:) + &
                                                                      dt * (0.5-bega) * Tdomain%sEdge(ne)%sSimu(i_simu)%Accel(:,:)

ne = Tdomain%specel(n)%near_edges(8)
Tdomain%specel(n)%sSimu(i_simu)%Forces(ngllx-1,1:nglly-2,ngllz-1,:) = Tdomain%sEdge(ne)%sSimu(i_simu)%Veloc(:,:) + &
                                                                      dt * (0.5-bega) * Tdomain%sEdge(ne)%sSimu(i_simu)%Accel(:,:)

ne = Tdomain%specel(n)%near_edges(9)
Tdomain%specel(n)%sSimu(i_simu)%Forces(1:ngllx-2,nglly-1,ngllz-1,:) = Tdomain%sEdge(ne)%sSimu(i_simu)%Veloc(:,:) + &
                                                                      dt * (0.5-bega) * Tdomain%sEdge(ne)%sSimu(i_simu)%Accel(:,:)

ne = Tdomain%specel(n)%near_edges(10)
Tdomain%specel(n)%sSimu(i_simu)%Forces(0,nglly-1,1:ngllz-2,:) = Tdomain%sEdge(ne)%sSimu(i_simu)%Veloc(:,:) + &
                                                                dt * (0.5-bega) * Tdomain%sEdge(ne)%sSimu(i_simu)%Accel(:,:)

ne = Tdomain%specel(n)%near_edges(11)
Tdomain%specel(n)%sSimu(i_simu)%Forces(0,1:nglly-2,ngllz-1,:) = Tdomain%sEdge(ne)%sSimu(i_simu)%Veloc(:,:) + &
                                                                dt * (0.5-bega) * Tdomain%sEdge(ne)%sSimu(i_simu)%Accel(:,:)

return
end subroutine get_PMLprediction_e2el
