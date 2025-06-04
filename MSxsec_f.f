      real function MSxsec_f(ev)

      implicit none
      
      include 'simulate.inc'

      real*8 MSxsec
      real*8 e_i, q2_loc, p_r_loc, theta_r_loc, phi_pq
      real*8 nu, sth2, the

      type(event) ev

      e_i =  ev%EIN
      q2_loc = ev%Q2
      p_r_loc = ev%PM
      theta_r_loc = ev%theta_rq
      phi_pq = ev%phi_pq

      nu = ev%NU

      sth2 = q2_loc/(4.*e_i*(e_i-nu))
      the = 2.*asin(sqrt(sth2))
      
      MSxsec_f = MSxsec(p_r_loc,q2_loc,theta_r_loc,e_i,the,phi_pq)

      return
      end

      
