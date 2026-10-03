"""Give each recovered surface its own OGS header and payload."""
from pathlib import Path
import shutil
root=Path(__file__).resolve().parents[1]
p=root/'MYSTRAN/Source/LK9/L91/WRITE_ELEM_STRESSES.f90'
shutil.copy2(p,root/'test/t6_upgrade/WRITE_ELEM_STRESSES_before_ogs.f90')
s=p.read_text()
a=s.index('      SUBROUTINE WRITE_OGS1_SURFACE_STRESS')
b=s.index('      CONTAINS',a)
t=s[a:b]
t=t.replace('LOGICAL                         :: USED_RECOVERY','LOGICAL                         :: USED_RECOVERY, LEGACY_OP2')
t=t.replace("         TABLE_NAME = 'OGS1    '\n         CALL WRITE_TABLE_HEADER(TABLE_NAME)\n         OGS_ITABLE = -3\n         CALL WRITE_ITABLE(OGS_ITABLE)\n",'')
start=t.index('      IF (WRITE_OP2) THEN\n         WRITE(OP2) 146')
end=t.index('      ENDIF',start)+len('      ENDIF\n')
header=t[start:end]
t=t[:start]+t[end:]
start=t.index('      IF (WRITE_OP2 .AND. USED_RECOVERY .AND. (NOUT > 0)) THEN')
end=t.index("      ELSE IF (WRITE_OP2 .AND. (FAMILY(1:4) == 'TRIA')) THEN",start)
payload=t[start:end]
payload=payload[payload.index('         WRITE(OP2)'):]
payload=payload.replace('I=1,NOUT','I=SURF_START(SURFACE_INDEX),SURF_END(SURFACE_INDEX)')
t=t[:start]+"      IF (LEGACY_OP2 .AND. (FAMILY(1:4) == 'TRIA')) THEN"+t[end+len("      ELSE IF (WRITE_OP2 .AND. (FAMILY(1:4) == 'TRIA')) THEN"):]
t=t.replace('      ELSE IF (WRITE_OP2) THEN','      ELSE IF (LEGACY_OP2) THEN')
insert='''      LEGACY_OP2 = WRITE_OP2 .AND. .NOT.(USED_RECOVERY .AND. NOUT > 0)
      IF (WRITE_OP2 .AND. .NOT.LEGACY_OP2) THEN
         DO SURF=1,NUM_GP_SURFACE
            IF (SURF_END(SURF) < SURF_START(SURF)) CYCLE
            IF (SURF_START(SURF) <= 0) CYCLE
            CALL WRITE_RECOVERED_SURFACE_OP2(SURF)
         ENDDO
      ELSE IF (LEGACY_OP2) THEN
         CALL WRITE_TABLE_HEADER('OGS1    ')
         OGS_ITABLE = -3
         CALL WRITE_ITABLE(OGS_ITABLE)
'''+header.replace('      IF (WRITE_OP2) THEN\n','',1).rsplit('      ENDIF\n',1)[0]+'''      ENDIF

'''
at=t.index('      IF (USED_RECOVERY .AND. (NOUT > 0)) THEN')
t=t[:at]+insert+t[at:]
t=t.replace('      IF (WRITE_OP2) WRITE(OP2) NVALUES','      IF (LEGACY_OP2) WRITE(OP2) NVALUES')
t=t.replace('      IF (WRITE_OP2) THEN\n         CALL END_OP2_TABLE(OGS_ITABLE)',
            '      IF (LEGACY_OP2) THEN\n         CALL END_OP2_TABLE(OGS_ITABLE)')
helper='''
      SUBROUTINE WRITE_RECOVERED_SURFACE_OP2(SURFACE_INDEX)
      INTEGER(LONG), INTENT(IN) :: SURFACE_INDEX
      OGS_ID = GP_SURFACE_IDS(SURFACE_INDEX)
      CALL WRITE_TABLE_HEADER('OGS1    ')
      OGS_ITABLE = -3
      CALL WRITE_ITABLE(OGS_ITABLE)
'''+header.replace('      IF (WRITE_OP2) THEN\n','',1).rsplit('      ENDIF\n',1)[0]+'''
      NVALUES = NUM_WIDE*2*(SURF_END(SURFACE_INDEX)-SURF_START(SURFACE_INDEX)+1)
      WRITE(OP2) NVALUES
'''+payload+'''
      CALL END_OP2_TABLE(OGS_ITABLE)
      ITABLE = 0
      END SUBROUTINE WRITE_RECOVERED_SURFACE_OP2

'''
s=s[:a]+t+s[b:]
b=s.index('      CONTAINS',a)+len('      CONTAINS\n')
s=s[:b]+helper+s[b:]
p.write_text(s)
