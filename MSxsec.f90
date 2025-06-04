function MSxsec(p_r_loc,q2_loc,theta_r_loc,e_i,the,phi_pq) result (sigma_eep)

  use response_functions_mod
  implicit none

  !include 'simulate.inc'

  ! event variables
  real :: p_r_loc, q2_loc, theta_r_loc
  real :: ei, q_lab, e_theta

  integer ::  i, j, read_status, int_status

  integer :: i_unit = 18  ! unit number for response function grid
	
  real:: Rl_loc, Rt_loc, Rlt_loc, Rtt_loc ! response functions

  ! for cross section calculation
  
  real :: e_i, the, q2, q2v, qv, e_f, p_f, dtr
  real :: phi_r, theta_r, phi_pq
  real :: sigma_eep
    
  real :: v_l, v_t, v_lt, v_tt, s_fact

  real, parameter :: pi = 3.141592654
  
  ! assign variables to the event variable
  phi_r = pi - phi_pq

  dtr = pi/180.
  
  !print*, ' Enter  p_r, q2, theta_r :'
  !read*, p_r_loc, q2_loc, theta_r_loc

  ! get interpolated response functions
  call resp_interp(p_r_loc, q2_loc, theta_r_loc/dtr , &
       Rl_loc, Rt_loc, Rlt_loc, Rtt_loc, int_status) 
  
  !print*, 'Rl, Rt, Rlt, Rtt = ', Rl_loc, Rt_loc, Rlt_loc, Rtt_loc
   print*, 'int_status  = ', int_status

  ! get kinematics for cross section
  !print*, 'enter Ei, theta_e, phi_r :'
  !read*, e_i, the, phi_r
  !the = the*dtr		! convert to radians
  !theta_r = theta_r_loc*dtr
  !phi_r = phi_r*dtr
  
  ! calculate final energy
  e_f = q2_loc/(4.*e_i*sin(the/2)**2)
  q2v = q2_loc + (e_i - e_f)**2
  qv = sqrt(q2v)
  
  p_f = sqrt(q2v + p_r_loc**2 - 2*qv*p_r_loc*cos(theta_r_loc))
  
  ! calculate cross sections
  !print*, '========== test_str_mod90: e_i, e_f, the, p_f, p_r_loc, phi_r'
  !print*, e_i, e_f, the, p_f, p_r_loc, phi_r
  !print*, '=========='

  ! call the function that calculates the cross section and return the response functions.
  
  sigma_eep = resp_calc_sigma(e_i, e_f, the, p_r_loc, phi_r, &  ! input
       v_l, v_t, v_lt, v_tt, s_fact)                                     ! output

  !print*, 'sigma, vl, vt, vlt, vtt, s_fact = ', sigma, v_l, v_t, v_lt, v_tt, s_fact 
  
 return  
end function MSxsec

subroutine init_MS(grid_dir, do_fsi, save_grid, wf_model)

  use response_functions_mod
  implicit none
  character*100 :: grid_dir
  integer :: do_fsi, save_grid, read_status
  integer :: wf_model ! wave function model selector 1 = Paris , 2 = AV18, 3 = CD-Bonn
     
  ! load grid an initialize interpolation
  if (wf_model .eq. 3) then          
     if (do_fsi .eq. 1) then
        read_status = resp_initialize('deut_MS/CD-Bonn/FSI/strfun_grid_CD-Bonn_FSI.data')
     else 
        read_status = resp_initialize('deut_MS/CD-Bonn/PWIA/strfun_grid_CD-Bonn_PWIA.data')
     endif
  endif

  if (read_status .ne. 0) then
     print*, 'cannot read response function file, cannot continue !'
     stop
  endif

end subroutine init_MS
