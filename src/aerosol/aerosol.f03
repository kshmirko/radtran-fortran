module aerosol
  use mathutils
  implicit none
contains

  pure function ext_a(cs, ndens) result(ext)
    ! cs in cm-3
    ! ext in m-1
    real(kind=dp), intent(in) ::  cs, ndens
    real(kind=dp)             ::  ext
    
    ext = cs*ndens*1.0d6 !
    return
  end function ext_a
  
  pure function tau_a_from_to(z0, z1, ext_a, hpbl) result(tau)
    real(kind=dp), intent(in) :: z0, z1, ext_a, hpbl
    real(kind=dp) :: tau
    
    
    tau = hpbl*(DEXP(-z0/hpbl)-DEXP(-z1/hpbl))*ext_a
    
    return
    
  end function tau_a_from_to
  
  
  pure function tau_a(cs, ndens, hpbl)
    real(kind=dp), intent(in) :: cs, ndens, hpbl
    real(kind=dp) :: tau_a
    tau_a = ext_a(cs, ndens)*hpbl
  end function tau_a
  
  pure function ext_a_h(ext, hpbl, h)
    real(kind=dp), intent(in) :: ext, hpbl, h
    real(kind=dp) :: ext_a_h
    ext_a_h = ext*DEXP(-h/hpbl)
    return
  end function ext_a_h
end module aerosol
