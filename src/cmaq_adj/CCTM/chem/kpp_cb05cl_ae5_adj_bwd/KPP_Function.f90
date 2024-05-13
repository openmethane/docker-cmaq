! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! The ODE Function of Chemical Model File
! 
! Generated by KPP-2.2 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : KPP_Function.f90
! Time                 : Tue Nov  9 23:09:28 2010
! Working directory    : /home/amir/carleton/amir/work/kpp/cb5_cl_ae5_aq
! Equation file        : cb5.kpp
! Output root filename : cb5
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE KPP_Function

  USE KPP_Parameters
  IMPLICIT NONE

! A - Rate for each equation
  REAL(kind=dp) :: A(KNREACT)

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Fun - time derivatives of variables - Agregate form
!   Arguments :
!      V         - Concentrations of variable species (local)
!      F         - Concentrations of fixed species (local)
!      RCT       - Rate constants (local)
!      Vdot      - Time derivative of variable species concentrations
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Fun ( V, F, RCT, Vdot )

! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)
! F - Concentrations of fixed species (local)
  REAL(kind=dp) :: F(NFIX)
! RCT - Rate constants (local)
  REAL(kind=dp) :: RCT(KNREACT)
! Vdot - Time derivative of variable species concentrations
  REAL(kind=dp) :: Vdot(NVAR)


! Computation of equation rates
  A(1) = RCT(1)*V(62)
  A(2) = RCT(2)*V(64)
  A(3) = RCT(3)*V(60)*V(63)
  A(4) = RCT(4)*V(62)*V(64)
  A(5) = RCT(5)*V(62)*V(64)
  A(6) = RCT(6)*V(63)*V(64)
  A(7) = RCT(7)*V(60)*V(62)
  A(8) = RCT(8)*V(60)
  A(9) = RCT(9)*V(60)
  A(10) = RCT(10)*V(19)
  A(11) = RCT(11)*V(19)
  A(12) = RCT(12)*V(60)*V(72)
  A(13) = RCT(13)*V(60)*V(68)
  A(14) = RCT(14)*V(66)
  A(15) = RCT(15)*V(66)
  A(16) = RCT(16)*V(63)*V(66)
  A(17) = RCT(17)*V(62)*V(66)
  A(18) = RCT(18)*V(62)*V(66)
  A(19) = RCT(19)*V(23)
  A(20) = RCT(20)*V(23)
  A(21) = RCT(21)*V(23)
  A(22) = RCT(22)*V(63)*V(63)
  A(23) = RCT(23)*V(62)*V(63)
  A(24) = RCT(24)*V(63)*V(72)
  A(25) = RCT(25)*V(26)
  A(26) = RCT(26)*V(26)*V(72)
  A(27) = RCT(27)*V(26)*V(26)
  A(28) = RCT(28)*V(62)*V(72)
  A(29) = RCT(29)*V(46)*V(72)
  A(30) = RCT(30)*V(63)*V(68)
  A(31) = RCT(31)*V(62)*V(68)
  A(32) = RCT(32)*V(32)
  A(33) = RCT(33)*V(32)*V(72)
  A(34) = RCT(34)*V(68)*V(68)
  A(35) = RCT(35)*V(68)*V(68)
  A(36) = RCT(36)*V(25)
  A(37) = RCT(37)*V(25)*V(72)
  A(38) = RCT(38)*V(19)
  A(39) = RCT(39)*V(72)
  A(40) = RCT(40)*V(64)*V(72)
  A(41) = RCT(41)*V(72)*V(72)
  A(42) = RCT(42)*V(72)*V(72)
  A(43) = RCT(43)*V(68)*V(72)
  A(44) = RCT(44)*V(64)*V(68)
  A(45) = RCT(45)*V(25)*V(64)
  A(46) = RCT(46)*V(64)*V(66)
  A(47) = RCT(47)*V(66)*V(72)
  A(48) = RCT(48)*V(66)*V(68)
  A(49) = RCT(49)*V(60)*V(66)
  A(50) = RCT(50)*V(66)*V(66)
  A(51) = RCT(51)*V(32)
  A(52) = RCT(52)*V(46)
  A(53) = RCT(53)*V(23)
  A(54) = RCT(54)*V(63)*V(67)
  A(55) = RCT(55)*V(53)*V(63)
  A(56) = RCT(56)*V(67)*V(68)
  A(57) = RCT(57)*V(53)*V(68)
  A(58) = RCT(58)*V(67)*V(67)
  A(59) = RCT(59)*V(53)*V(53)
  A(60) = RCT(60)*V(53)*V(67)
  A(61) = RCT(61)*V(59)*V(72)
  A(62) = RCT(62)*V(59)
  A(63) = RCT(63)*V(40)*V(72)
  A(64) = RCT(64)*V(40)
  A(65) = RCT(65)*V(44)*V(72)
  A(66) = RCT(66)*V(72)
  A(67) = RCT(67)*V(63)*V(69)
  A(68) = RCT(68)*V(68)*V(69)
  A(69) = RCT(69)*V(69)*V(69)
  A(70) = RCT(70)*V(43)*V(72)
  A(71) = RCT(71)*V(43)
  A(72) = RCT(72)*V(34)*V(72)
  A(73) = RCT(73)*V(61)*V(72)
  A(74) = RCT(74)*V(61)
  A(75) = RCT(75)*V(61)
  A(76) = RCT(76)*V(61)*V(64)
  A(77) = RCT(77)*V(61)*V(66)
  A(78) = RCT(78)*V(61)*V(68)
  A(79) = RCT(79)*V(37)
  A(80) = RCT(80)*V(37)*V(63)
  A(81) = RCT(81)*V(37)*V(68)
  A(82) = RCT(82)*V(27)*V(72)
  A(83) = RCT(83)*V(58)*V(64)
  A(84) = RCT(84)*V(58)*V(72)
  A(85) = RCT(85)*V(58)*V(66)
  A(86) = RCT(86)*V(58)
  A(87) = RCT(87)*V(63)*V(70)
  A(88) = RCT(88)*V(62)*V(70)
  A(89) = RCT(89)*V(20)
  A(90) = RCT(90)*V(20)
  A(91) = RCT(91)*V(68)*V(70)
  A(92) = RCT(92)*V(69)*V(70)
  A(93) = RCT(93)*V(67)*V(70)
  A(94) = RCT(94)*V(70)*V(70)
  A(95) = RCT(95)*V(28)*V(72)
  A(96) = RCT(96)*V(28)
  A(97) = RCT(97)*V(29)*V(72)
  A(98) = RCT(98)*V(57)*V(64)
  A(99) = RCT(99)*V(57)*V(72)
  A(100) = RCT(100)*V(57)*V(66)
  A(101) = RCT(101)*V(57)
  A(102) = RCT(102)*V(63)*V(71)
  A(103) = RCT(103)*V(62)*V(71)
  A(104) = RCT(104)*V(30)
  A(105) = RCT(105)*V(30)
  A(106) = RCT(106)*V(30)*V(72)
  A(107) = RCT(107)*V(68)*V(71)
  A(108) = RCT(108)*V(69)*V(71)
  A(109) = RCT(109)*V(67)*V(71)
  A(110) = RCT(110)*V(71)*V(71)
  A(111) = RCT(111)*V(70)*V(71)
  A(112) = RCT(112)*V(54)*V(72)
  A(113) = RCT(113)*V(47)
  A(114) = RCT(114)*V(47)
  A(115) = RCT(115)*V(47)*V(62)
  A(116) = RCT(116)*V(52)*V(64)
  A(117) = RCT(117)*V(52)*V(72)
  A(118) = RCT(118)*V(52)*V(60)
  A(119) = RCT(119)*V(52)*V(66)
  A(120) = RCT(120)*V(49)*V(64)
  A(121) = RCT(121)*V(49)*V(72)
  A(122) = RCT(122)*V(49)*V(60)
  A(123) = RCT(123)*V(49)*V(66)
  A(124) = RCT(124)*V(51)*V(64)
  A(125) = RCT(125)*V(51)*V(72)
  A(126) = RCT(126)*V(51)*V(60)
  A(127) = RCT(127)*V(51)*V(66)
  A(128) = RCT(128)*V(22)*V(72)
  A(129) = RCT(129)*V(31)*V(63)
  A(130) = RCT(130)*V(31)
  A(131) = RCT(131)*V(48)*V(72)
  A(132) = RCT(132)*V(48)*V(66)
  A(133) = RCT(133)*V(38)*V(62)
  A(134) = RCT(134)*V(38)*V(68)
  A(135) = RCT(135)*V(45)
  A(136) = RCT(136)*V(45)*V(72)
  A(137) = RCT(137)*V(45)*V(60)
  A(138) = RCT(138)*V(24)*V(72)
  A(139) = RCT(139)*V(39)*V(72)
  A(140) = RCT(140)*V(39)
  A(141) = RCT(141)*V(55)*V(64)
  A(142) = RCT(142)*V(55)*V(72)
  A(143) = RCT(143)*V(55)*V(60)
  A(144) = RCT(144)*V(55)*V(66)
  A(145) = RCT(145)*V(56)*V(72)
  A(146) = RCT(146)*V(56)*V(60)
  A(147) = RCT(147)*V(56)*V(66)
  A(148) = RCT(148)*V(56)
  A(149) = RCT(149)*V(50)*V(64)
  A(150) = RCT(150)*V(50)*V(72)
  A(151) = RCT(151)*V(50)*V(60)
  A(152) = RCT(152)*V(50)*V(66)
  A(153) = RCT(153)*V(18)*V(72)
  A(154) = RCT(154)*V(35)*V(72)
  A(155) = RCT(155)*V(33)*V(72)
  A(156) = RCT(156)*V(55)*V(62)
  A(157) = RCT(157)*V(17)
  A(158) = RCT(158)*V(21)
  A(159) = RCT(159)*V(60)*V(65)
  A(160) = RCT(160)*V(41)*V(41)
  A(161) = RCT(161)*V(41)*V(63)
  A(162) = RCT(162)*V(41)*V(68)
  A(163) = RCT(163)*V(42)*V(72)
  A(164) = RCT(164)*V(42)
  A(165) = RCT(165)*V(65)
  A(166) = RCT(166)*V(54)*V(65)
  A(167) = RCT(167)*V(33)*V(65)
  A(168) = RCT(168)*V(49)*V(65)
  A(169) = RCT(169)*V(52)*V(65)
  A(170) = RCT(170)*V(51)*V(65)
  A(171) = RCT(171)*V(55)*V(65)
  A(172) = RCT(172)*V(61)*V(65)
  A(173) = RCT(173)*V(58)*V(65)
  A(174) = RCT(174)*V(57)*V(65)
  A(175) = RCT(175)*V(34)*V(65)
  A(176) = RCT(176)*V(35)*V(65)
  A(177) = RCT(177)*V(36)*V(72)
  A(178) = RCT(178)*V(7)*V(63)
  A(179) = RCT(179)*V(7)*V(68)
  A(180) = RCT(180)*V(10)*V(63)
  A(181) = RCT(181)*V(10)*V(68)
  A(182) = RCT(182)*V(14)*V(72)
  A(183) = RCT(183)*V(13)*V(63)
  A(184) = RCT(184)*V(13)*V(68)
  A(185) = RCT(185)*V(16)*V(60)
  A(186) = RCT(186)*V(16)*V(72)
  A(187) = RCT(187)*V(16)*V(66)

! Aggregate function
  Vdot(1) = A(142)
  Vdot(2) = A(149)+A(150)+A(151)+A(152)
  Vdot(3) = A(153)
  Vdot(4) = A(153)
  Vdot(5) = A(178)
  Vdot(6) = A(179)
  Vdot(7) = 0.765*A(128)-A(178)-A(179)
  Vdot(8) = A(180)
  Vdot(9) = A(181)
  Vdot(10) = 0.804*A(138)-A(180)-A(181)
  Vdot(11) = A(183)
  Vdot(12) = A(184)
  Vdot(13) = 0.764*A(182)-A(183)-A(184)
  Vdot(14) = -A(182)
  Vdot(15) = A(185)+A(186)+A(187)
  Vdot(16) = -A(185)-A(186)-A(187)
  Vdot(17) = -A(157)+0.3*A(160)
  Vdot(18) = -A(153)
  Vdot(19) = A(9)-A(10)-A(11)-A(38)
  Vdot(20) = A(88)-A(89)-A(90)
  Vdot(21) = -A(158)+A(162)
  Vdot(22) = -A(128)
  Vdot(23) = A(18)-A(19)-A(20)-A(21)-A(53)
  Vdot(24) = -A(138)
  Vdot(25) = A(34)+A(35)-A(36)-A(37)+A(42)-A(45)
  Vdot(26) = 2*A(23)+A(24)-A(25)-A(26)-2*A(27)
  Vdot(27) = A(80)-A(82)+0.37*A(122)
  Vdot(28) = 0.8*A(91)-A(95)-A(96)+0.8*A(107)
  Vdot(29) = 0.2*A(91)+0.1*A(92)+0.1*A(93)-A(97)+0.2*A(107)+0.1*A(108)+0.1*A(109)
  Vdot(30) = A(103)-A(104)-A(105)-A(106)
  Vdot(31) = 0.56*A(128)-A(129)-A(130)+0.3*A(138)
  Vdot(32) = A(31)-A(32)-A(33)-A(51)
  Vdot(33) = -A(155)-A(167)
  Vdot(34) = 0.63*A(69)-A(72)-A(175)
  Vdot(35) = -A(154)-A(176)
  Vdot(36) = A(165)+A(166)+A(167)+0.3*A(170)+0.15*A(171)+A(172)+A(173)+A(174)+A(175)+A(176)-A(177)
  Vdot(37) = A(78)-A(79)-A(80)-A(81)
  Vdot(38) = 0.4*A(131)+A(132)-A(133)-A(134)
  Vdot(39) = 0.2*A(137)+0.8*A(138)-A(139)-A(140)+0.168*A(145)+0.85*A(146)
  Vdot(40) = A(56)+A(57)-A(63)-A(64)
  Vdot(41) = A(159)-2*A(160)-A(161)-A(162)
  Vdot(42) = -A(163)-A(164)+A(168)+A(169)+0.7*A(170)+0.85*A(171)
  Vdot(43) = A(68)-A(70)-A(71)+A(81)
  Vdot(44) = -A(65)+A(73)+A(74)+A(75)+A(76)+A(77)+A(86)+A(101)+0.2*A(116)+0.33*A(118)+A(120)+0.63*A(122)+0.1*A(124)+0.25&
               &*A(126)+A(135)+2*A(136)+0.69*A(137)+A(140)+0.066*A(143)+0.334*A(145)+0.225*A(146)+0.643*A(147)+0.333*A(148)&
               &+0.001*A(151)+A(163)+A(164)+A(172)
  Vdot(45) = 0.9*A(129)+0.3*A(131)-A(135)-A(136)-A(137)
  Vdot(46) = 2*A(19)+2*A(20)+A(28)-A(29)+A(48)-A(52)+A(61)+A(77)+A(85)+A(100)+A(132)+0.15*A(147)
  Vdot(47) = 0.76*A(112)-0.98*A(113)-A(114)-A(115)+0.76*A(166)
  Vdot(48) = 0.36*A(128)+A(130)-A(131)-A(132)+A(134)+0.2*A(138)
  Vdot(49) = -A(120)-A(121)-A(122)-A(123)-A(168)
  Vdot(50) = -A(149)-A(150)-A(151)-A(152)
  Vdot(51) = -A(124)-A(125)-A(126)-A(127)-A(170)
  Vdot(52) = -A(116)-A(117)-A(118)-A(119)-A(169)+0.3*A(170)
  Vdot(53) = -A(55)-A(57)-2*A(59)-A(60)+0.13*A(112)+0.04*A(113)+0.01*A(116)+0.09*A(119)+0.088*A(142)+0.25*A(150)+0.18&
               &*A(151)+0.25*A(152)+0.009*A(155)+0.13*A(166)+0.009*A(167)
  Vdot(54) = -0.66*A(61)-0.66*A(62)-1.11*A(112)-2.1*A(113)+0.2*A(116)-0.7*A(117)-A(118)-A(119)+0.1*A(124)+1.1*A(138)&
               &+0.25*A(141)+0.35*A(143)+2.4*A(144)+1.565*A(145)+0.36*A(146)+1.282*A(147)+0.832*A(148)+5.12*A(149)+1.66&
               &*A(150)+7*A(151)+2.4*A(156)-1.11*A(166)-A(169)+0.3*A(170)
  Vdot(55) = -A(141)-A(142)-A(143)-A(144)-A(156)-A(171)
  Vdot(56) = 0.75*A(141)+0.912*A(142)+0.65*A(143)+0.2*A(144)-A(145)-A(146)-A(147)-A(148)+0.2*A(156)+A(171)
  Vdot(57) = 0.33*A(61)+0.33*A(62)+0.5*A(63)+0.5*A(64)-A(98)-A(99)-A(100)-A(101)+0.05*A(112)+0.5*A(113)+0.3*A(116)+0.62&
               &*A(117)+0.32*A(118)+0.56*A(119)+0.22*A(121)+0.66*A(124)+0.7*A(125)+0.35*A(126)+0.64*A(127)+0.03*A(137)+0.15&
               &*A(143)+0.8*A(144)+0.12*A(145)+0.357*A(147)+0.15*A(149)+0.47*A(150)+0.21*A(151)+0.47*A(152)+0.05*A(154)+0.8&
               &*A(156)+0.05*A(166)+0.67*A(169)+0.55*A(170)-A(174)
  Vdot(58) = 0.33*A(61)+0.33*A(62)+0.5*A(63)+0.5*A(64)-A(83)-A(84)-A(85)-A(86)+A(102)+A(106)+0.9*A(108)+0.9*A(109)+2&
               &*A(110)+A(111)+0.06*A(112)+0.6*A(113)+0.2*A(116)+0.33*A(117)+0.18*A(118)+0.35*A(119)+1.24*A(124)+1.3*A(125)&
               &+0.65*A(126)+1.18*A(127)+0.252*A(145)+0.02*A(146)+0.067*A(148)+0.9*A(154)+0.991*A(155)+0.06*A(166)+0.991&
               &*A(167)+0.33*A(169)+0.45*A(170)-A(173)+A(176)
  Vdot(59) = A(55)-A(61)-A(62)+A(115)+0.1*A(129)+A(133)+0.8*A(144)+0.85*A(147)+0.53*A(152)+0.8*A(156)
  Vdot(60) = A(2)-A(3)-A(7)-A(8)-A(9)-A(12)-A(13)-A(49)+0.2*A(91)+0.2*A(107)-A(118)-A(122)-A(126)-A(137)-A(143)-A(146)&
               &-A(151)-A(159)
  Vdot(61) = 0.33*A(61)+0.33*A(62)+A(67)+1.37*A(69)+A(71)+A(72)-A(73)-A(74)-A(75)-A(76)-A(77)-A(78)+A(79)+A(92)+0.1&
               &*A(108)+0.2*A(116)+0.8*A(117)+0.74*A(118)+A(119)+A(120)+1.56*A(121)+A(122)+2*A(123)+0.25*A(126)+A(136)+0.7&
               &*A(137)+0.5*A(141)+0.629*A(142)+0.6*A(143)+0.167*A(145)+0.15*A(146)+0.282*A(147)+0.9*A(148)+0.28*A(150)+0.24&
               &*A(151)+0.1*A(154)+A(168)-A(172)+A(175)
  Vdot(62) = -A(1)+A(3)-A(4)-A(5)+A(6)-A(7)+A(14)+2*A(16)-A(18)+A(21)+2*A(22)-A(23)+A(26)+A(27)-A(28)+A(30)-A(31)+A(32)&
               &+A(33)+A(46)+A(47)+A(49)+2*A(50)+0.61*A(51)+A(52)+A(53)+A(54)+A(62)+A(67)+A(80)+A(87)-A(88)+A(89)+A(90)&
               &+A(102)-A(103)+A(104)+A(105)+A(106)-A(115)+A(119)+A(123)+A(127)+0.9*A(129)-A(133)+0.2*A(144)+0.47*A(152)&
               &-A(156)+A(161)
  Vdot(63) = A(1)-A(3)+A(4)-A(6)+A(15)-A(16)+A(17)-2*A(22)-A(23)-A(24)+A(25)+A(27)-A(30)-A(54)-A(55)-A(67)-A(80)-A(87)&
               &-A(102)-A(129)+0.2*A(156)-A(161)
  Vdot(64) = A(1)-A(2)-A(4)-A(5)-A(6)+A(8)+A(10)+A(14)-A(40)+A(41)-A(44)-A(45)-A(46)-A(76)-A(83)-A(98)-A(116)-A(120)&
               &-A(124)+0.5*A(126)-A(141)-A(149)
  Vdot(65) = 2*A(157)+A(158)-A(159)+1.4*A(160)+A(161)+A(163)+A(164)-A(165)-A(166)-A(167)-A(168)-A(169)-A(170)-A(171)&
               &-A(172)-A(173)-A(174)-A(175)-A(176)+A(177)
  Vdot(66) = A(5)+A(7)-A(14)-A(15)-A(16)-A(17)-A(18)+A(21)+A(29)-A(46)-A(47)-A(48)-A(49)-2*A(50)+0.39*A(51)+A(53)-A(77)&
               &-A(85)-A(100)-A(119)-A(123)-A(127)-A(132)-A(144)-A(147)-A(152)
  Vdot(67) = -A(54)-A(56)-2*A(58)-A(60)+A(63)+0.3*A(70)-A(93)+A(102)+0.9*A(108)-A(109)+2*A(110)+A(111)+0.87*A(112)+0.96&
               &*A(113)+0.2*A(116)+0.8*A(117)+0.22*A(118)+0.91*A(119)+0.7*A(120)+A(121)+A(123)+0.1*A(124)+A(125)+0.08*A(128)&
               &+0.6*A(131)+A(136)+0.03*A(137)+0.5*A(138)+A(139)+0.25*A(141)+0.991*A(142)+0.2*A(143)+A(144)+0.713*A(145)&
               &+0.064*A(146)+0.075*A(147)+0.7*A(148)+1.25*A(150)+0.76*A(151)+1.03*A(152)+0.1*A(154)+0.991*A(155)+A(156)&
               &+0.87*A(166)+0.991*A(167)+2*A(168)+2*A(169)+1.7*A(170)+A(171)
  Vdot(68) = A(12)-A(13)-A(30)-A(31)+A(32)-2*A(34)-2*A(35)+A(37)+A(38)+A(39)+A(40)-A(43)-A(44)+A(45)+A(47)-A(48)+0.61&
               &*A(51)-A(56)-A(57)+A(61)+A(62)+A(64)+A(65)+A(67)-A(68)+0.74*A(69)+0.3*A(70)+A(71)+A(72)+A(73)+2*A(74)+A(76)&
               &+A(77)-A(78)+A(79)+A(80)-A(81)+A(82)+A(86)-A(91)+0.9*A(92)+A(101)+A(102)-A(107)+A(108)+2*A(110)+A(111)+0.11&
               &*A(112)+0.94*A(113)+A(114)+0.3*A(116)+0.95*A(117)+0.44*A(118)+1.7*A(120)+A(121)+0.13*A(122)+0.1*A(124)&
               &+A(125)+0.5*A(126)+A(127)+0.44*A(128)+0.9*A(129)+A(130)+0.6*A(131)-A(134)+A(135)+2*A(136)+0.76*A(137)+0.7&
               &*A(138)+A(140)+0.25*A(141)+0.912*A(142)+0.066*A(143)+0.8*A(144)+0.503*A(145)+0.154*A(146)+0.925*A(147)+1.033&
               &*A(148)+0.75*A(150)+0.07*A(151)+0.28*A(152)+A(153)+A(154)+A(155)+0.8*A(156)-A(162)+A(164)+0.11*A(166)+A(167)&
               &+A(168)+A(169)+A(170)+A(171)+A(172)+A(175)+A(176)
  Vdot(69) = A(66)-A(67)-A(68)-2*A(69)+0.7*A(70)+A(86)+A(87)-0.1*A(92)+0.9*A(93)+2*A(94)+A(96)+A(97)+A(101)-A(108)&
               &+A(111)+A(165)
  Vdot(70) = A(83)+A(84)+A(85)-A(87)-A(88)+A(89)+A(90)-A(91)-A(92)-A(93)-2*A(94)+A(95)-A(111)+A(135)+A(136)+0.62*A(137)&
               &+A(139)+A(140)+0.21*A(145)+0.114*A(146)+0.967*A(148)+A(173)
  Vdot(71) = A(98)+A(99)+A(100)-A(102)-A(103)+A(104)+A(105)-A(107)-A(108)-A(109)-2*A(110)-A(111)+0.25*A(141)+0.2*A(143)&
               &+0.25*A(145)+0.075*A(147)+0.39*A(151)+A(174)
  Vdot(72) = 2*A(11)-A(12)+A(13)-A(24)+A(25)-A(26)-A(28)-A(29)+A(30)-A(33)+2*A(36)-A(37)+A(38)-A(39)-A(40)-2*A(41)-2&
               &*A(42)-A(43)+A(44)+A(45)-A(47)+0.39*A(51)+A(52)-A(61)-A(63)+A(64)-A(65)-A(66)-A(70)+A(71)-A(72)-A(73)+A(76)&
               &-A(82)+A(83)-A(84)-A(95)+A(96)-A(97)+A(98)-A(99)-A(106)-A(112)+0.1*A(116)-A(117)+0.1*A(118)+0.3*A(120)&
               &-A(121)+0.13*A(122)-A(125)+0.5*A(126)-A(128)-A(131)-A(136)+0.08*A(137)-A(138)-A(139)-A(142)+0.266*A(143)&
               &-A(145)+0.268*A(146)-A(150)+0.57*A(151)-A(153)-A(154)-A(155)+A(158)-A(163)-A(177)
      
END SUBROUTINE Fun

! End of Fun function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE KPP_Function

