! --- spc1_fix begin --- !
! Reject malformed one-line free-field SPC1 explicit lists that spill past the
! supported parent-card field range. Valid long lists must use continuation
! cards or THRU form.
         IF (JCARD(10)(1:) /= ' ') THEN
            JERR      = JERR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1129) JCARD(1), JCARD(2)
            WRITE(F06,1129) JCARD(1), JCARD(2)
         ENDIF
! --- spc1_fix end --- !

! Parent card contract for SPC1 explicit-list format:
!   parent card:       fields 4-9 -> grid IDs
!   continuation card: fields 2-9 -> grid IDs
