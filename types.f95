module types

! Declare types (single and double precision)
integer, parameter :: SP = selected_real_kind(p=6,r=37)
integer, parameter :: DP = selected_real_kind(p=13,r=200)

end module types
