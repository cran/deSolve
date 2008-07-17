C  file satres.f
C  Initializer for parameter common block
      subroutine initmod(odeparms)
      external odeparms
      double precision parms(14)
      common /myparms/parms
      call odeparms(14, parms)
      return
      end
C  Compartments are:
C  y(1)  central compartment
C  y(2)  second compartment
C  y(3)  filtrate compartment
C  y(4)  'Gut'
C  y(5)  Total eliminated

C  Derivatives and one output variable
      subroutine derivs(neq, t, y, ydot, out, ip)
      integer neq, IP(*)
      double precision t, y(neq), ydot(neq), out(*)
      double precision Vc, Vt, kd, ka, Tm, KT, Kfil, Vfil, free, BW,
     $     Dose, Doseint, Qd, TDose
      common /myparms/Vc, Vt, kd, ka, Tm, KT, Kfil, Vfil, free, BW,
     $     Dose, Doseint, Qd, TDose

      if (ip(1) < 1) call rexit("nout should be at least 1")

      ydot(1) = (ka * y(4) - Qd * free * y(1) + Qd * y(2)) / Vc -
     $     kfil * y(1) * free + Tm * y(3) / (KT + y(3))
      ydot(2) = (free * Qd * y(1) - Qd * y(2)) / Vt
      ydot(3) = (Vc * kfil * y(1) * free - Vc * Tm * y(3) / (KT + y(3)) -
     $     Vc * kfil * y(3)) / Vfil
      ydot(4) = -ka * y(4)
      ydot(5) = Vc * kfil * y(3)

      out(1) = y(1) * Vc + y(2) * Vt + y(3) * Vfil + y(4) + y(5)
      return
      end
C  END file satres.f
