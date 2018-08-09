module rayleigh
  use mathutils
  implicit none
  
  real(kind=dp), parameter  ::  h_mol=8427.3_dp ! m
  real(kind=dp), parameter  ::  h_max=45000.0_dp !m
contains
  
  pure function tau_m(wl) result(tau)
    real(kind=dp), intent(in) ::  wl
    real(kind=dp)             ::  tau
  
    tau = 0.0089_dp / (wl**4.0_dp) *&
      (1.0_dp+ &
        0.0113_dp/(wl**2.0_dp)+&
        0.00013_dp/(wl**4.0_dp))
    return
  end function tau_m
  
  pure function ext_m(wl) result(ext)
    real(kind=dp), intent(in) ::  wl
    real(kind=dp)             ::  ext
    
    ext = tau_m(wl) / h_mol
    return
  end function ext_m
  
  pure function tau_m_from_to(z0, z1, wl) result(tau)
    real(kind=dp), intent(in) :: z0, z1, wl
    real(kind=dp) :: tau, ext
    
    ext = ext_m(wl)
    tau = h_mol*(DEXP(-z0/h_mol)-DEXP(-z1/h_mol))*ext
    
    return
    
  end function tau_m_from_to
  
  pure function tau_m_h(wl, h)
    real(kind=dp), intent(in) :: wl, h
    real(kind=dp) :: tau_m_h
    
    tau_m_h = ext_m(wl)*DEXP(-h/h_mol)
  end function tau_m_h

end module rayleigh
