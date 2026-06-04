        real(8) function potential(x)
                implicit none
                real(8), intent(in) :: x

                ! Harmonic oscillator with force constant = 1
                potential = 0.5d0*x*x
                
                ! Morse potential
!               potential = (1d0-exp(-x/sqrt(2d0)))**2

                ! Particle in a 1D infinite box with L = pi
               ! potential = 1d8
               ! if (abs(x)<2d0*atan(1d0)) potential = 0d0
                
                return
        end function potential
