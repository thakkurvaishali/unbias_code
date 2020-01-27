!---------------------------------------------------------------------!
!Written by Shalini Awasthi (ashalini@iitk.ac.in)
!---------------------------------------------------------------------!
PROGRAM WSMTD_rw_3D
IMPLICIT NONE
REAL*8 gridmin1, gridmax1, griddif1, dummy2,v, &
        gridmin2, gridmax2, griddif2,num, &
        gridmin3, gridmax3, griddif3, ss, &
        gridmin4, gridmax4, griddif4, &
        prob,den,alpha,fes,fes1,grid,Ro,prob3d, &
        cv1,cv2,ht,kt,kt0,ktb,bias_fact,ct,vbias,hill,width, &
        diff_s2,diff_s3,ds2,hh,dum,s1,s2,s3,s4,dummy11,cv3,cv4,prob2d,prob1d
ALLOCATABLE cv1(:),cv2(:),ht(:),vbias(:),ct(:),hill(:,:),width(:), &
            ss(:,:),prob(:,:,:,:),fes(:,:,:,:),fes1(:,:),grid(:,:), &
            cv3(:),cv4(:),prob1d(:),prob2d(:,:),prob3d(:,:,:)    
INTEGER mtd_steps,md_steps,dummy1,i,j,index1,index2,k,index3, &
        t_min,t_max,nbin1, nbin2,nbin3,nbin4,w_cv,w_hill,index4, &
        i_mtd,i_md,i_s2,i_s1,mtd_max,i_s3,narg,i_s4
LOGICAL pmf,inpgrid
CHARACTER*120 :: arg 
REAL*8, PARAMETER :: kb=1.9872041E-3 !kcal K-1 mol-1
REAL*8, PARAMETER :: au_to_kcal = 627.51 
REAL*8, PARAMETER :: pi = 4.d0*Atan(1.d0)
!      character(len=50) f1

OPEN(11,FILE='colvar_mtd',STATUS='unknown')
OPEN(12,FILE='parvar_mtd',STATUS='unknown')
!      OPEN(13,FILE='CONSTRAINT',STATUS='unknown')
OPEN(14,FILE='cvmdck_mtd',STATUS='unknown')
!
      CALL get_steps(11,mtd_steps)
      CALL get_steps(14,md_steps)
!
!      CALL check_files(11,w_hill)
!      CALL check_files(14,w_cv)
!      CALL check_files(13,i)
!      IF(i.ne.w_cv)THEN
!        print *, 
!     &    '!!Print Freq. of cvmdck_mtd & CONSTRAINTS are not equal!!'
!        print *, '    cvmdck_mtd Print Frq. =',w_cv
!        print *, '    CONSTRAINTS Print Frq. =',i
!        STOP 
!      END IF 

kt=1000.D0
kt0=300.D0
bias_fact=2000.D0
t_min=1
t_max=md_steps
pmf=.FALSE.
inpgrid=.false.
t_max=md_steps
narg = IARGC()
print*, narg
DO i=1,narg
  CALL GETARG(i,arg)
  IF(INDEX(arg,'-CV_T').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)kt
  ELSEIF(INDEX(arg,'-T').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)kt0
  ELSE IF(INDEX(arg,'-dT').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)bias_fact
  ELSE IF(INDEX(arg,'-tmin').NE.0)THEN
     CALL GETARG(i+1,arg)
     READ(arg,*)t_min
! ELSE IF(INDEX(arg,'-tmax').NE.0)THEN
!     CALL GETARG(i+1,arg)
!    READ(arg,*)t_max
!   IF(t_max.gt.md_steps)STOP '!!ERROR: t_max > total MD steps'
  ELSE IF(INDEX(arg,'-pmf').NE.0)THEN
     pmf=.true.
  ELSE IF(INDEX(arg,'-grid').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)gridmin1
      CALL GETARG(i+2,arg)
      READ(arg,*)gridmax1
      CALL GETARG(i+3,arg)
      READ(arg,*)griddif1
      CALL GETARG(i+4,arg)
      READ(arg,*)gridmin2
      CALL GETARG(i+5,arg)
      READ(arg,*)gridmax2
      CALL GETARG(i+6,arg)
      READ(arg,*)griddif2
      CALL GETARG(i+7,arg)
      READ(arg,*)gridmin3
      CALL GETARG(i+8,arg)
      READ(arg,*)gridmax3
      CALL GETARG(i+9,arg)
      READ(arg,*)griddif3
      CALL GETARG(i+10,arg)
      READ(arg,*)gridmin4
      CALL GETARG(i+11,arg)
      READ(arg,*)gridmax4
      CALL GETARG(i+12,arg)
      READ(arg,*)griddif4
            inpgrid=.true.
  ELSE IF(INDEX(arg,'-pfrqMD').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)w_cv
  ELSE IF(INDEX(arg,'-dtMTD').NE.0)THEN
      CALL GETARG(i+1,arg)
      READ(arg,*)w_hill
  END IF
END DO

WRITE(*,'(A,I10)')'No: of MTD steps        =',mtd_steps
WRITE(*,'(A,I10)')'No: of MD  steps        =',md_steps
WRITE(*,'(A,F9.2)')'Physical System Temperature T0 (K)=',kt0
WRITE(*,'(A,F9.2)')'CV Temperature T (K)   =',kt
WRITE(*,'(A,F9.2)')'Bias Factor (K)         =',bias_fact
WRITE(*,'(A,I10)')'Print Freq. cvmdck_mtd  =',w_cv
WRITE(*,'(A,I10)')'Freq. of Hill Update    =',w_hill
WRITE(*,'(A,I10)')'Reweigtht: Min step     =',t_min
WRITE(*,'(A,I10)')'Reweigtht: Max step     =',t_max
IF(pmf)WRITE(*,'(A)')'PMF will be written in PMF.dat'
!
bias_fact=(kt0+bias_fact)/kt0
kt=kb*kt
ktb=kt*bias_fact

!Make CV file containing CV1(t), CV2(t) from CONSTRAINT and cvmdck_mtd
ALLOCATE(cv1(md_steps),cv2(md_steps),cv3(md_steps),cv4(md_steps))
ALLOCATE( vbias(md_steps))
ALLOCATE(ht(mtd_steps),ct(mtd_steps),hill(mtd_steps,2), &
&         width(mtd_steps),ss(mtd_steps,2))

OPEN(16,file='cv.dat',status='unknown')
DO i_md=1,md_steps
!read 4 CVs
  READ(14,*)dummy1,dummy1,dummy11,dummy11,dummy11,dummy11,dummy11,dummy11, &
&            cv2(i_md),cv3(i_md),cv1(i_md),cv4(i_md)                         !V-edit
  WRITE(16,'(I10,4F16.6)')i_md,cv1(i_md),cv2(i_md),cv3(i_md),cv4(i_md)     !V-edit
END DO
WRITE(*,'(A)')'CV values written in cv.dat'
!

alpha=bias_fact/(bias_fact-1.D0)

nbin1 = NINT((gridmax1-gridmin1)/griddif1)+1
nbin2 = NINT((gridmax2-gridmin2)/griddif2)+1
nbin3 = NINT((gridmax3-gridmin3)/griddif3)+1
nbin4 = NINT((gridmax4-gridmin4)/griddif4)+1
WRITE(*,'(7X,4A10)')'GRIDMIN','GRIDMAX','GRIDBIN','GRIDSIZE'
WRITE(*,'(A10,3F8.4,I10)')'US  COORD:', &
&          gridmin1,gridmax1,griddif1,nbin1
WRITE(*,'(A10,3F8.4,I10)')'MTD1 COORD:', &
&          gridmin2,gridmax2,griddif2,nbin2
WRITE(*,'(A10,3F8.4,I10)')'MTD2 COORD:', &
&          gridmin3,gridmax3,griddif3,nbin3
WRITE(*,'(A10,3F8.4,I10)')'TEMP COORD:', &
&          gridmin4,gridmax4,griddif4,nbin4


ALLOCATE(prob(nbin1,nbin2,nbin3,nbin4),fes(nbin1,nbin2,nbin3,nbin4))
ALLOCATE(fes1(nbin2,nbin3),grid(nbin2,2))

DO i_mtd=1,mtd_steps
  READ(11,*) dummy1,hill(i_mtd,1),hill(i_mtd,2),dummy2
  READ(12,*) dummy1,dummy2,width(i_mtd),ht(i_mtd)
  ht(i_mtd)=ht(i_mtd)*au_to_kcal
END DO

CLOSE(11);CLOSE(12);CLOSE(13);CLOSE(14);CLOSE(16)


!==================================calculate c(t)===============================!
DO i_s2=1,nbin2 !grid over cv2 on which MTD is being done
  grid(i_s2,1)=gridmin2+DFLOAT(i_s2-1)*griddif2
END DO

DO i_s3=1,nbin3
  grid(i_s3,2)=gridmin3+DFLOAT(i_s3-1)*griddif3
END DO

OPEN(21,FILE='ct.dat',STATUS='unknown')
DO i_mtd=1,mtd_steps
  ds2=width(i_mtd)*width(i_mtd)
  hh=ht(i_mtd)
  num=0.D0
  den=0.D0
   
 DO i_s3=1,nbin3
     diff_s3=grid(i_s3,2)-hill(i_mtd,2)
     diff_s3=diff_s3*diff_s3*0.5D0

   DO i_s2=1,nbin2
     diff_s2=grid(i_s2,1)-hill(i_mtd,1)
     diff_s2=diff_s2*diff_s2*0.5D0

    
     fes1(i_s2,i_s3)=fes1(i_s2,i_s3)-hh*DEXP(-(diff_s3+diff_s2)/ds2)
     num=num+DEXP(-fes1(i_s2,i_s3)/kt)
     den=den+DEXP(-fes1(i_s2,i_s3)/ktb)
  END DO
 END DO
 
     ct(i_mtd)=kt*DLOG(num/den)
WRITE(21,'(I10,F16.8)')i_mtd,ct(i_mtd)
END DO
CLOSE(21)
WRITE(*,'(A)')'Ct factor written in ct.dat'

!=====================================calculate v(s,t)====================================!
DO i_md=1,md_steps
   mtd_max=(i_md*w_cv/w_hill)
!   ss(i_mtd,1)=cv2(i_md)              !changed mtd to md
!   ss(i_mtd,2)=cv3(i_md)
   !PRINT*, i_mtd, i_md
!
   dum=0.d0
 DO i_mtd=1,mtd_max

   ds2=width(i_mtd)*width(i_mtd)
   hh=ht(i_mtd)/alpha
!
   diff_s2=cv2(i_md)-hill(i_mtd,1)
   diff_s2=diff_s2*diff_s2*0.5D0
!
   diff_s3=cv3(i_md)-hill(i_mtd,2)
   diff_s3=diff_s3*diff_s3*0.5D0

   dum=dum+hh*DEXP(-(diff_s2+diff_s3)/ds2)

   END DO

  vbias(i_md)=dum
END DO
!=============================calculate prob (unbiased from MTD potential)==========================!       
      den=0.d0
      prob=0.d0
!
!
DO i_md=1,md_steps
  IF((i_md.GT.t_min).AND.(i_md.LT.t_max))THEN
    index1 = nint((cv1(i_md)-gridmin1)/griddif1) +1
    index2 = nint((cv2(i_md)-gridmin2)/griddif2) +1
    index3 = nint((cv3(i_md)-gridmin3)/griddif3) +1              !V-edit
    index4 = nint((cv4(i_md)-gridmin4)/griddif4) +1              !V-edit
    i_mtd=i_md*w_cv/w_hill+1
    dum=vbias(i_md) - ct(i_mtd)
    dum=dexp(dum/kt)
    prob(index1,index2,index3,index4)=prob(index1,index2,index3,index4) + dum
!          prob(index1,index2)=prob(index1,index2) + dum
!     den=den+1.D0
          den=den+dum
  END IF
END DO
dum=den*griddif1*griddif2*griddif3*griddif4
!      dum=den*griddif1*griddif2
den=1.D0/dum
OPEN(2,FILE='Pu_4D.dat',STATUS='unknown')
DO i_s1=1,nbin1
  s1=DFLOAT(i_s1-1)*griddif1+gridmin1
  DO i_s2=1,nbin2
    s2=DFLOAT(i_s2-1)*griddif2+gridmin2
    DO i_s3=1,nbin3
      s3=DFLOAT(i_s3-1)*griddif3+gridmin3
      DO i_s4=1,nbin4
        s4=DFLOAT(i_s4-1)*griddif4+gridmin4
        prob(i_s1,i_s2,i_s3,i_s4)=prob(i_s1,i_s2,i_s3,i_s4)*den
        WRITE(2,'(5E16.8)')s1,s2,s3,s4,prob(i_s1,i_s2,i_s3,i_s4)
      END DO
      WRITE(2,*)
    END DO
    WRITE(2,*)
  END DO
  WRITE(2,*)
END DO
WRITE(*,'(A)')'Unbiased distribution written in Pu_4D.dat'
CLOSE(2)

!!==================================2D-Probability================================!
ALLOCATE(prob2d(nbin1,nbin4))
prob2d=0.d0 ; den=0.d0
DO i_s1=1,nbin1
    DO i_s4=1,nbin4
      dum=0.d0
      DO i_s3=1,nbin3
        DO i_s2=1,nbin2
          dum=dum+prob(i_s1,i_s2,i_s3,i_s4)
        END DO
      END DO
      den=den+dum*griddif3*griddif2*griddif1*griddif4
      prob2d(i_s1,i_s4)=dum*griddif2*griddif3
    END DO 
END DO

DO i_s1=1,nbin1
    prob2d(i_s1,1:nbin4)=prob2d(i_s1,1:nbin4)/den
END DO

PRINT*,'den =',den,'dum =',dum


OPEN(2,FILE='Pu-temp.dat',STATUS='unknown')
DO i_s1=1,nbin1
  s1=DFLOAT(i_s1-1)*griddif1+gridmin1
  DO i_s4=1,nbin4
    s4=DFLOAT(i_s4-1)*griddif4+gridmin4 
   !prob2d(i_s1,i_s2)=prob2d(i_s1,i_s2)*den
    WRITE(2,'(4E16.8)')s1,s4,prob2d(i_s1,i_s4)
  END DO
  WRITE(2,*) 
END DO
CLOSE(2)

!DEALLOCATE(prob3d)
!!!!==================================3D-Probability================================!
!!ALLOCATE(prob3d(nbin1,nbin3,nbin4))
!!prob3d=0.d0 ; den=0.d0
!!DO i_s1=1,nbin1
!!  DO i_s3=1,nbin3
!!    DO i_s4=1,nbin4
!!      dum=0.d0
!!      DO i_s2=1,nbin2
!!        dum=dum+prob(i_s1,i_s2,i_s3,i_s4)
!!      END DO
!!      den=den+dum*griddif3*griddif2*griddif1*griddif4
!!      prob3d(i_s1,i_s3,i_s4)=dum*griddif2
!!    END DO
!!  END DO 
!!END DO
!!
!!DO i_s1=1,nbin1
!!  DO i_s3=1,nbin3
!!    prob3d(i_s1,i_s3,1:nbin4)=prob3d(i_s1,i_s3,1:nbin4)/den
!!  END DO
!!END DO
!!
!!PRINT*,'den =',den,'dum =',dum
!!
!!
!!OPEN(2,FILE='Pu-3D.dat',STATUS='unknown')
!!DO i_s1=1,nbin1
!!  s1=DFLOAT(i_s1-1)*griddif1+gridmin1
!!  DO i_s3=1,nbin3
!!    s3=DFLOAT(i_s3-1)*griddif3+gridmin3
!!    DO i_s4=1,nbin4
!!      s4=DFLOAT(i_s4-1)*griddif4+gridmin4 
!!     !prob2d(i_s1,i_s2)=prob2d(i_s1,i_s2)*den
!!      WRITE(2,'(4E16.8)')s1,s3,s4,prob3d(i_s1,i_s3,i_s4)
!!    END DO
!!    WRITE(2,*) 
!!  END DO
!!  WRITE(2,*)
!!END DO
!!CLOSE(2)
!!
!!!DEALLOCATE(prob3d)

!!===================================2D-Probability================================!
!
!ALLOCATE(prob2d(nbin1,nbin4))
!prob2d=0.d0;den=0.d0
!!print*,den
!DO i_s1=1,nbin1
!  DO i_s4=1,nbin4
!    dum=0.d0
!    DO i_s3=1,nbin3
!      dum=dum+prob3d(i_s1,i_s3,i_s4)
!    END DO
!    den=den+dum*griddif3*griddif1*griddif4
!    prob2d(i_s1,i_s4)=dum*griddif3
!  END DO 
!END DO
!
!DO i_s1=1,nbin1
!prob2d(nbin1,1:nbin4)=prob2d(nbin1,1:nbin4)/den
!END DO
!
!PRINT*,'den =',den,'dum =',dum
!
!OPEN(2,FILE='Pu-temp.dat',STATUS='unknown')
!DO i_s1=1,nbin1
!  s1=DFLOAT(i_s1-1)*griddif1+gridmin1
!  DO i_s4=1,nbin4
!    s4=DFLOAT(i_s4-1)*griddif4+gridmin4 
!  ! prob2d(i_s1,i_s3)=prob2d(i_s1,i_s3)*den
!  WRITE(2,'(3E16.8)')s1,s4,prob2d(i_s1,i_s4)
!  END DO
!WRITE(2,*) 
!END DO
!CLOSE(2)

DEALLOCATE(prob2d)

!=================================1D-Probability===================================!

allocate(prob1d(nbin1))
prob1d=0.d0
den=0.d0
do i_s1=1,nbin1
  dum=0.d0
  do i_s2=1,nbin2
    do i_s3=1,nbin3
      do i_s4=1,nbin4
        dum=dum+prob(i_s1,i_s2,i_s3,i_s4)
      end do
    end do 
  end do
  prob1d(i_s1)=dum*griddif2*griddif3*griddif4
  den=den+dum*griddif2*griddif3*griddif4*griddif1
end do

prob1d(1:nbin1)=prob1d(1:nbin1)/den

OPEN(2,FILE='Pu.dat',STATUS='unknown')
DO i_s1=1,nbin1
  s1=DFLOAT(i_s1-1)*griddif1+gridmin1
  if(prob1d(i_s1).ne.prob1d(i_s1)) prob1d(i_s1)=1.0D-16 !remove NaN
  WRITE(2,'(2E16.8)')s1,prob1d(i_s1)
END DO
CLOSE(2)

deallocate(prob1d)

END PROGRAM WSMTD_rw_3D
!---------------------------------------------------------------------!


!---------------------------------------------------------------------!
      SUBROUTINE get_steps(iunit,nsteps)
      IMPLICIT NONE
      INTEGER iunit, nsteps
      INTEGER ios
      nsteps=0
      REWIND(iunit)
      Read_Loop: DO
         READ(iunit,*,IOSTAT=ios)
         IF(ios.ne.0)EXIT Read_Loop
         nsteps=nsteps+1
      END DO Read_Loop 
      REWIND(iunit)
      END 
!---------------------------------------------------------------------!
      SUBROUTINE check_files(iunit,dt)
      IMPLICIT NONE
      INTEGER iunit, dt
      INTEGER ii, jj,i,ios
      dt=0
      i=2
      REWIND(iunit)
      READ(iunit,*)ii
      READ(iunit,*)jj
      dt=jj-ii
      ii=jj
      RLoop: DO 
        i=i+1
        READ(iunit,*,IOSTAT=ios)jj
        IF(ios.ne.0)EXIT RLoop
        IF(jj.ne.ii+dt)THEN
           print *, '!!ERROR: Steps are not at constant stride!!'
           print *, '!!       Unit No:',iunit,'!!'
           print *, '!!       Line No:',i,'!!'
           print *, '!! Expected stride =', dt,'!!'
           print *, '!! Actual stride =', jj-ii,'!!'
           STOP
        END IF
        ii=jj
      END DO RLoop
      REWIND(iunit)
      END 
!---------------------------------------------------------------------!
      SUBROUTINE get_gridmin_max(iunit,gridmin1,gridmax1,griddif1, &
      &                           gridmin2,gridmax2,griddif2,gridmin3,gridmax3,griddif3)
      IMPLICIT NONE 
      INTEGER :: iunit 
      REAL*8  :: gridmin1, gridmax1, griddif1, &
      &           gridmin2, gridmax2, griddif2, &
      &           gridmin3, gridmax3, griddif3
      INTEGER :: ii, ios
      REAL*8  :: cv1, cv2, cv3
      INTEGER, PARAMETER :: Def_Grid_Size=101
      REWIND(iunit)
      READ(iunit,*,IOSTAT=ios)ii,cv1,cv2,cv3
      if(ios.ne.0)stop 'ERROR reading CV.dat'
      gridmin1=cv1
      gridmax1=cv1
      gridmin2=cv2
      gridmax2=cv2
      gridmin3=cv3
      gridmax3=cv3
      RLoop: DO 
        READ(iunit,*,IOSTAT=ios)ii,cv1,cv2,cv3
        if(ios.ne.0)EXIT RLoop
        gridmin1=MIN(gridmin1,cv1)
        gridmin2=MIN(gridmin2,cv2)
        gridmin3=MAX(gridmin3,cv3)
        gridmax1=MAX(gridmax1,cv1)
        gridmax2=MAX(gridmax2,cv2)
        gridmax3=MAX(gridmax3,cv3)
      END DO RLoop
      griddif1=(gridmax1-gridmin1)/DFLOAT(Def_Grid_Size)
      griddif2=(gridmax2-gridmin2)/DFLOAT(Def_Grid_Size)
      griddif3=(gridmax3-gridmin3)/DFLOAT(Def_Grid_Size)
      END
!---------------------------------------------------------------------!
