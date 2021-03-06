!------------------------------------------------------------- 12/12/2018
! This module contains parameters needed for SerraLINE
  module parms
  implicit none

  integer, parameter :: sp= kind(0.0), dp = kind(0.d0), &
                        nlinR=9, nlinY=6    !purine and pyrimidine ring atoms

  real(dp), parameter :: pi = datan(1.0_dp)*2.0_dp*2.0_dp, &
                       & rad_to_deg = 180.0_dp/pi
  !BASES
!-----------------------------------------------------------------------

  character(len=4), dimension(6), parameter :: &

  !ADENINE
  & A_l = ["A   ", "A5  ", "A3  ", "DA  ", "DA5 ", "DA3 "], &

  !GUANINE
  & G_l = ["G   ", "G5  ", "G3  ", "DG  ", "DG5 ", "DG3 "], &

  !CYTOSINE
  & C_l = ["C   ", "C5  ", "C3  ", "DC  ", "DC5 ", "DC3 "], &

  !THYMINE
  & T_l = ["T   ", "T5  ", "T3  ", "DT  ", "DT5 ", "DT3 "], &

  !URACIL
  & U_l = ["U   ", "U5  ", "U3  ", "DU  ", "DU5 ", "DU3 "]
!-----------------------------------------------------------------------


  !RING ATOMS
!-----------------------------------------------------------------------
  character(len=4), parameter :: &

  N1 = "N1  ", C2 = "C2  ", N3 = "N3  ", C4 = "C4  ", C5 = "C5  ", &
  C6 = "C6  ", N7 = "N7  ", C8 = "C8  ", N9 = "N9  "

  private
  public sp, dp, nlinR, nlinY, &
       & pi, rad_to_deg, &
       & A_l, G_l, C_l, T_l, U_l, N1, C2, N3, C4, C5, C6, N7, C8, N9 

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

  contains

!-----------------------------------------------------------------------

  end module parms
