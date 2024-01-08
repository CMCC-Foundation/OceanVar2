    subroutine def_coef_bc(vect,a1,a2,a3)

       implicit none
       real*8   ::   vect(9)
       real*8   ::   a1,a2,a3,m_c
      

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!  Triggs - Sdika Backward conditions!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      m_c=1./(1.+a1-a2+a3)
      m_c=m_c/(1.-a1-a2-a3)
      m_c=m_c/(1.+a2+(a1-a3)*a3)

      vect(1) =m_c*(-a3*a1+1.-a3**2-a2)
      vect(2) =m_c*(a3+a1)*(a2+a3*a1)
      vect(3) =m_c*a3*(a1+a3*a2)
      vect(4) =m_c*(a1+a3*a2)
      vect(5) =m_c*(-(a2-1.)*(a2+a3*a1))
      vect(6) =m_c*(-a3*(a3*a1+a3**2+a2-1))
      vect(7) =m_c*(a3*a1+a2+a1**2-a2**2)
      vect(8) =m_c*((a1*a2)+(a3*a2**2)-(a1*a3**2)-(a3**3)-(a3*a2)+a3)
      vect(9) =m_c*(a3*(a1+a3*a2))

    end subroutine def_coef_bc


