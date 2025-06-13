      real function MSxsec_f(ev)

      implicit none
      
      include 'simulate.inc'

      real*8 MSxsec
      real*8 e_i, q2_loc, p_r_loc, theta_r_loc, phi_pq
      real*8 nu, sth2, the

      real*8 toGeV, rtd

      data toGeV/1e-3/
      
      type(event) ev

      rtd = 180./pi
      
      e_i =  ev%EIN*toGeV
      q2_loc = ev%Q2*toGeV**2
      p_r_loc = ev%PM*toGeV
      theta_r_loc = ev%theta_rq
      phi_pq = ev%phi_pq

      nu = ev%NU*toGeV

      sth2 = q2_loc/(4.*e_i*(e_i-nu))
      the = 2.*asin(sqrt(sth2))

c      write(47, *) e_i, q2_loc, p_r_loc, theta_r_loc, phi_pq
      write(47,*) the*rtd
      
      MSxsec_f = MSxsec(p_r_loc,q2_loc,theta_r_loc,e_i,the,phi_pq)

      return
      end

      
