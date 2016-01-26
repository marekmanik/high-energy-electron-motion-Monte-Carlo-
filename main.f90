program reid
 use ziggurat
 use params
 use fgsl
 implicit none

! integer :: seed = 0 ! 1961656697  !1 400401039
 integer :: seed = 285594245  ! 1756249931

 integer  :: track_me_id = 97297 ! set the most energetic index from 'the_chosen_one' and set seed as well

 integer(8), parameter :: Ndata=10000

 integer(8), parameter :: N_time_steps = 30000
 integer(8),parameter ::  Npart = 100000

 integer :: collnr
 integer(8) :: i, i_particle, i_real
 integer :: j, istat, unit_id

 real(DP), dimension(ndim) :: vel ! electron velocity, gas_molecule velocity
 real(DP), dimension(ndim) :: pos ! electron position

 real(DP), dimension(ndim,Npart) :: coords, velocities
 real(DP) :: mean_eps, max_tot_freq !, mean_eps_clas
 
 real(DP) :: track_me_coords(ndim,N_time_steps)
 real(DP) :: track_me_velocity(ndim,N_time_steps)


! character(20) :: chout = ""

 real(DP) :: TotTime, accel, Efield
 real(DP)    ::  max_energy   = 0.0_dp
 integer(8)  :: max_energy_id = 0
 real(DP) :: rand,gam


 real(DP) :: m_elec(Npart)
 real(DP),allocatable :: collfreq(:,:),tot_freq(:),energy(:)
 real(DP),allocatable :: eps_okr(:,:) !,threshold(:)
 integer, allocatable :: coll_type(:)
 real(DP),allocatable :: sigma(:), nu(:),energy_treshold(:)
 real(DP),allocatable :: mGas(:), eps_bar(:)
!!!!!!!!!!!! !!
 type(fgsl_interp_accel)  :: gsl_acc
 type(fgsl_interp_type)   :: gsl_int_t
 type(fgsl_spline),allocatable :: gsl_spline_eps_okh(:)
 type(fgsl_spline),allocatable :: gsl_spline_collfreq(:)
!!!!!!!!!!!!!!
type css_data_tp
 character(len=100) :: descr
 integer            :: css_type
 real(DP)           :: threshold
 real(DP), dimension(Ndata) :: energy
 real(DP), dimension(Ndata) :: css
 real(DP), dimension(Ndata) :: eps_okr
end type

type(css_data_tp), allocatable :: N2_css_data(:)
type(css_data_tp), allocatable :: O2_css_data(:)
!!!!!!!!!!!!!!!!

integer ::  Ncss_N2, Ncss_O2
real(DP) :: delta_T

 call init_random_seed(seed)

 call read_css(N2_css_data, "CSN2.txt", Ndata, Ncss_N2)
 call read_css(O2_css_data, "CSO2.txt", Ndata, Ncss_O2)

 call init_coll_frequency(N2_css_data,O2_css_data,Ncss_N2,Ncss_O2,&
                               collnr,collfreq,eps_okr,tot_freq,coll_type,&
                               energy_treshold,mGas,eps_bar,energy)


 allocate(gsl_spline_eps_okh(collnr))
 allocate(gsl_spline_collfreq(collnr))
 allocate(sigma(collnr))
 allocate(nu(collnr))

 gsl_acc   = fgsl_interp_accel_alloc(); ! global
 gsl_int_t = fgsl_interp_linear;        ! global

 call interpol_init(collfreq,energy,gsl_spline_collfreq,Ndata)
 call interpol_init(eps_okr,energy,gsl_spline_eps_okh,Ndata)

 max_tot_freq =  maxval(tot_freq)
 max_energy   =  0.0_dp ! track the most energetic one

 m_elec = me_const
 delta_T = delta_MC/max_tot_freq
 
 i_real=0

! E_over_N = 750.0_dp*Townsend_const
! accel = q0_const*E_over_N*Ngas/me_const

 Efield = 20_dp*Ek_const

 call init_particle(coords,velocities,ndim,Npart)

 do i=1,N_time_steps

   do i_particle=1,Npart

     vel(:) = velocities(:,i_particle)
     pos(:) = coords(:,i_particle)

     accel = q0_const*Efield/(m_elec(i_particle))

     pos(1) = pos(1) + vel(1)*delta_T
     pos(2) = pos(2) + vel(2)*delta_T
     pos(3) = pos(3) + vel(3)*delta_T + 0.5_dp*accel*delta_T**2
     vel(3) = vel(3) + accel*delta_T

     gam  = 1.0/sqrt(1-norm2(vel)**2/c_light**2)    
     m_elec(i_particle) = gam*me_const

     mean_eps = 0.5_dp*m_elec(i_particle)*(norm2(vel))**2

     if(mean_eps > max_energy )then
       max_energy    = mean_eps
       max_energy_id = i_particle
     endif

     rand = uni()

     if (delta_MC < rand) then
      ! particles don't collide
      ! do nothing
     else 
        call lint(mean_eps,nu,energy,gsl_spline_collfreq,energy_treshold)
        call lint(mean_eps,sigma,energy,gsl_spline_eps_okh,energy_treshold)
        call null_collision_method(vel,mean_eps,m_elec(i_particle),i_real,energy_treshold,& 
                                   mGas,eps_bar,nu,max_tot_freq,sigma)
     endif

     coords(:,i_particle) = pos(:)
     velocities(:,i_particle) = vel(:)

   enddo ! i_particle

   TotTime = i*delta_T
!  write(chout,fmt='(I0.8)') i

   track_me_coords(:,i) = coords(:,track_me_id)
   track_me_velocity(:,i) = velocities(:,track_me_id)


   if (mod(i,N_time_steps/100) == 0) then
     call print_statistics(TotTime,coords,velocities,Efield,ndim,Npart)
   endif

 enddo ! dt Time steps
!!!!!!!!!!!!!!!!!!!!!!

! write all particles
 open(newunit=unit_id,file="XXXparticle.data",action='write', iostat=istat)
    do i_particle = 1, Npart
      write(unit_id,fmt=*) TotTime, (coords(j,i_particle), j=1,3),  (velocities(j,i_particle), j=1,3)
    enddo
 close(unit_id)

! write the most energetic one
 open(newunit=unit_id,file="XXXthe_chosen_one",action='write',position='append', iostat=istat)
   write(unit_id, *) seed, max_energy/ElectronVolt, max_energy_id
 close(unit_id)


! write the tracked one
 open(newunit=unit_id,file="XXXtracked_particle",action='write', iostat=istat)
 do i = 1, N_time_steps
   write(unit_id, *) i*delta_T, track_me_coords(:,i), track_me_velocity(:,i)
 enddo
 close(unit_id)


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine print_statistics(time,coords,velocities,Efield,ndim,Npart)
 integer,    intent(in) :: ndim
 integer(8), intent(in) :: Npart
 real(DP),   intent(in) :: time
 real(DP),   intent(in) :: coords(ndim,Npart)
 real(DP),   intent(in) :: velocities(ndim,Npart)
 real(DP),   intent(in) :: Efield

 integer(8) :: i
 real(DP) :: energies(Npart)

 do i=1, Npart
   gam         = 1.0_dp/sqrt(1-norm2(velocities(:,i))**2/c_light**2)
   energies(i) = me_const*c_light**2*(gam-1)
 enddo

 print *, time, sum(energies)/Npart/ElectronVolt, sum(velocities(3,:))/Npart/Efield, maxval(coords(3,:))

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_coll_frequency(N2_css_data,O2_css_data,Ncss_N2,Ncss_O2,&
                               collnr,collfreq,eps_okr,tot_freq,coll_type,&
                               threshold,mGas,eps_bar,energy)
 integer, intent(out) :: collnr
 integer, intent(in) :: Ncss_N2,Ncss_O2
 type(css_data_tp) :: N2_css_data(Ncss_N2)
 type(css_data_tp) :: O2_css_data(Ncss_O2)
 real(DP),intent(out),allocatable  :: energy(:),tot_freq(:)
 real(DP),intent(out),allocatable :: collfreq(:,:), eps_okr(:,:)
 integer,intent(out),allocatable :: coll_type(:)
 real(DP),intent(out),allocatable :: threshold(:)
 real(DP),intent(out),allocatable :: mGas(:), eps_bar(:)
 integer(8) :: i,j,j_tmp
 integer :: unit_id
 
 collnr =  Ncss_O2 + Ncss_N2

 allocate(collfreq(collnr,Ndata))
 allocate(eps_okr(collnr,Ndata))
 allocate(coll_type(collnr))
 allocate(threshold(collnr))
 allocate(mGas(collnr))
 allocate(eps_bar(collnr))
 allocate(tot_freq(Ndata))
 allocate(energy(Ndata))


 do i = 1, Ndata
   energy(i) = N2_css_data(1)%energy(i)*ElectronVolt
 enddo


! N2:
 do j = 1, Ncss_N2
   do i = 1, Ndata
     collfreq(j,i) = ppN2*Ngas*N2_css_data(j)%css(i)*velo(energy(i))*1e-4_dp ! cm2 -> m2
     eps_okr(j,i) = N2_css_data(j)%eps_okr(i)
   enddo
   coll_type(j)  = N2_css_data(j) % css_type
   threshold(j)  = N2_css_data(j) % threshold * ElectronVolt
        mGas(j)  = mGas_N2
     eps_bar(j)  = eps_bar_N2
 enddo

! O2:
 do j = 1, Ncss_O2
   j_tmp = j + Ncss_N2
   do i = 1, Ndata
     collfreq(j_tmp,i) = ppO2*Ngas*O2_css_data(j)%css(i)*velo(energy(i))*1e-4_dp ! cm2 -> m2
     eps_okr(j_tmp,i) = O2_css_data(j)%eps_okr(i)
   enddo
   coll_type(j_tmp)  = O2_css_data(j) % css_type
   threshold(j_tmp)  = O2_css_data(j) % threshold * ElectronVolt
        mGas(j_tmp)  = mGas_O2
      eps_bar(j_tmp) = eps_bar_O2
 enddo

 open(newunit=unit_id,file="tot_freq.data",action="write")

 do i=1,Ndata
   tot_freq(i) = sum(collfreq(:,i))
   write(unit_id,*) energy(i)/ElectronVolt, tot_freq(i)
 enddo

 close(unit_id)

end subroutine

real(DP) &
function velo(en) 
real(DP) :: en 
real(DP),parameter :: E0_elec = me_const*c_light**2
! velo  = c_light*sqrt(1-(me_const*c_light**2/(en*q0_const + me_const*c_light**2))**2)
! en [J]
 velo  = c_light*sqrt(1-(E0_elec/(en + E0_elec))**2)
end function velo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine collision(v1,v2,m1,m2,theta,phi,eps,v1_out,v2_out) !,i )
 real(DP), intent(in), dimension(ndim) :: v1, v2 ! initial velocities
 real(DP), intent(in) :: m1, m2 ! masses of colliding particels,
 real(DP), intent(in) :: eps !energy gain in nonelastic collision
 real(DP), intent(in) :: theta, phi ! scattering angles
 real(DP), intent(out), dimension(ndim) :: v1_out, v2_out ! resulting velocities

 real(DP), dimension(ndim) :: vt, vr, vr_init, vr_out ! center of mass velocity
 real(DP) :: m12, cr, vr_norm2, vr_norm2_prime, argument
 real(DP) :: A, cos_theta, cos_phi, sin_theta, sin_phi, red_mass
! integer(8) ::  i

 m12 = m1+m2
 red_mass  = m1*m2/m12
! red_mass  = m_red
 vt = (v1*m1+v2*m2)/m12
 vr = v1 - v2
 vr_norm2 = NORM2(vr)
 
 argument = vr_norm2**2+2.0_dp*eps/red_mass

 if( argument < 0.0_dp )then
  argument = 0.0_dp
 endif

 vr_norm2_prime = sqrt(argument)
 vr_init = vr
!! Get direction of scattered relative velocity (unit vr)
!! Byrd
 vr = vr/vr_norm2 ! normalize relative velocity vector
 A  = sqrt(vr(2)**2+vr(3)**2)
 cr = 1.0_dp

 cos_theta = cos(theta)
 sin_theta = sin(theta)
 sin_phi = sin(phi)
 cos_phi = cos(phi)

 vr_out(1) = cos_theta*vr(1) + sin_theta*sin_phi*A
 vr_out(2) = cos_theta*vr(2) + sin_theta*(cr*vr(3)*cos_phi-vr(1)*vr(2)*sin_phi)/A
 vr_out(3) = cos_theta*vr(3) - sin_theta*(cr*vr(2)*cos_phi+vr(1)*vr(3)*sin_phi)/A

 vr_out = vr_out*vr_norm2_prime
!!!!!!!!!!!!!

 v1_out = vt + m2/m12*vr_out
 v2_out = vt - m1/m12*vr_out

end subroutine collision
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine null_collision_method(vel,eps_p,m_elec,i_real,energy_treshold,mGas,eps_bar,nu,max_tot_freq,eps)
 real(DP), intent(inout), dimension(ndim)  :: vel
 integer(8), intent(inout) :: i_real
 real(DP), intent(in) :: energy_treshold(collnr), mGas(collnr),eps_bar(collnr)
 real(DP),intent(in) :: max_tot_freq, eps_p,m_elec
 real(DP) :: rand, vel_norm2, alfa,R,sum_nu
 real(DP) :: energy_gain , theta, phi,nu(collnr), eps(collnr),R_theta,eps_s
 real(DP),dimension(ndim) :: vel_out,vel_gas_out, vel_gas
 integer(8) :: j

! vel_gas(1) = rnor()
! vel_gas(2) = rnor()
! vel_gas(3) = rnor()

 vel_gas = 0.0_dp ! Gas motionless (MOSS)

 vel_norm2 = norm2(vel)
 rand = uni()
 R = uni()
 R_theta = uni()


 phi   = uni()*TWO_PI
 j=1
 alfa=0.0_dp

 sum_nu = sum(nu(:))

 if( sum_nu/max_tot_freq  < rand) then
  !null collision
 else
  do
   if ((alfa <= R) .and. (alfa+nu(j)/sum_nu > R)) then   !kolizia s prislusnou zrazkovou frekvenciou
#ifdef OKHRIM
      theta = acos(1-2*R_theta*(1-eps(j))/(1+eps(j)*(1-2*R_theta)))
#else
      theta = acos(1-2*uni())
#endif


    if ( coll_type(j) == 3 )then ! ionization
      eps_s = eps_bar(j)*tan(uni()*atan((eps_p - energy_treshold(j))/(2*eps_bar(j))))
      energy_gain = - energy_treshold(j) -  eps_s
    else
      energy_gain = - energy_treshold(j)
    endif

    call collision(vel,vel_gas,m_elec,mGas(j), theta,phi,energy_gain,vel_out,vel_gas_out)
    vel = vel_out
    i_real=i_real+1
    exit
   else         
    alfa = alfa + nu(j)/sum_nu
    j = j+1
   endif
  enddo
 endif
end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_particle(pos,vel,ndim,Npart)
 integer    :: ndim 
 integer(8) :: Npart
 real(DP), intent(out), dimension(ndim,Npart) :: pos
 real(DP), intent(out), dimension(ndim,Npart) :: vel

 pos = 0.0_dp

! Normal distribution at temperature Te_init 
 do i=1,Npart
   vel(:,i) = (/rnor(),rnor(),rnor()/)*sqrt(kb_const*Te_init/me_const)
 enddo

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_random_seed(seed)
  implicit none
  integer, intent(inout)  :: seed
  integer ::  un, istat

  if( seed == 0 )then! get random seed
    open(newunit=un, file="/dev/urandom", access="stream", &
                   form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
      read(un) seed
      close(un)
    else
      print  *, "RNG not seeded well"
      seed =124126
    endif

  endif

  seed = abs(seed)
!#ifdef DEBUG
  print *,  "# Seed: ", seed
!#endif
  call zigset(seed)
end subroutine init_random_seed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine interpol_init(y,x,gsl_spline,N)
 real(DP), intent(in)  :: y(collnr,Ndata),x(Ndata)
 type(fgsl_spline), intent(out)  :: gsl_spline(collnr)
 integer(8), intent(in) :: N
 integer(8) :: tmp

 do i = 1,  collnr
   gsl_spline(i) = fgsl_spline_alloc (gsl_int_t, Ndata);
   tmp = fgsl_spline_init(gsl_spline(i), x, y(i,:), N);
 enddo

end  subroutine interpol_init


subroutine lint(mean_eps,sigma,energy,gsl_spline,energy_treshold)
 real(DP), intent(in)  :: mean_eps
 real(DP), intent(in)  :: energy_treshold(collnr)
 real(DP), intent(out), dimension(collnr) :: sigma
 type(fgsl_spline), intent(in)  :: gsl_spline(collnr)
 real(DP), intent(in), dimension(Ndata)    :: energy
 real(DP) :: eps_tmp
 integer(8) :: j
 
 eps_tmp = mean_eps

 if( eps_tmp  < energy(1) )then
   eps_tmp = energy(1)
 elseif(eps_tmp > energy(Ndata))then
   eps_tmp = energy(Ndata)
 endif

 do j=1,collnr
    if (eps_tmp < energy_treshold(j)) then
      sigma(j) = 0.0_dp
    else
       sigma(j) = fgsl_spline_eval(gsl_spline(j), eps_tmp , gsl_acc)
    endif

 enddo

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_css(css_dt, filename, Ndata, Ncss)
 type(css_data_tp), allocatable, intent(inout) ::  css_dt(:)
 character(len=*)   :: filename
 integer(8), intent(in):: Ndata
 integer, intent(out):: Ncss

 integer :: unit_id, i
 integer(8) :: j

 open(newunit=unit_id,file=filename,action="read",status="old")

 read(unit_id,*) Ncss
! write(*,*) 'Nb of css:',  Ncss
 allocate(css_dt(Ncss))

 do i = 1, Ncss

  read(unit_id,'(a)') css_dt(i) % descr
  read(unit_id, * )   css_dt(i) % css_type
  read(unit_id, * )   css_dt(i) % threshold
!  print *, css_dt(i) % descr, css_dt(i) % css_type, css_dt(i) % threshold

  do j = 1, Ndata
   read(unit_id, * ) css_dt(i) % energy(j), &
                     css_dt(i) % css(j), &
                     css_dt(i) % eps_okr(j)
  enddo

 enddo

close(unit_id)

end subroutine read_css
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end program
