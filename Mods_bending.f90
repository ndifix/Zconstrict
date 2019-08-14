module cdk

contains

          SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed

            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))

            CALL SYSTEM_CLOCK(COUNT=clock)

            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)

            DEALLOCATE(seed)
          END SUBROUTINE


!=========================================================

   subroutine paras(nstep,jobid,nthreads,jbend,nturnover,jtreadmil,noutput,nrforce,nftsz,zlen,nftszstart,xgap, &
              pi,delta,INVDELTA,beta,k_z,l_z,l_diag,kz_thet,kb,angle,dl_z,l_z12,l_z34,l_diag_new, &
              k_a,la_cen,l_mb_a,l_a_z,ktether,ltether,mbrad,rrad,p_tether,l_mem,k_mem,kmemb,kgly,kpep,fturgor,pturgor, &
              wthick,lsqueez,kwall,walrad,eta,rho,dxsol,dyzsol,ptread,phyd,tadd,kexcl,kpair,gap,kradial, &
              viscos,printtime,k_scale,k_scale_mb,mwid,rwid,kwid,kplane,kcir,kcircon,jgradient,jnogradient,fgrad,growstep)

   implicit none

   character chara*32

   integer(kind=8)::nstep

   integer n,jobid,nthreads,jbend,nturnover,jtreadmil,noutput,nrforce,nftsz,zlen,nftszstart,jgradient,jnogradient

   real(kind=8)::pi,delta,INVDELTA,beta,k_z,l_z,l_diag,kz_thet,kb,angle,dl_z,l_z12,l_z34,l_diag_new
   real(kind=8)::k_a,la_cen,l_mb_a,l_a_z,ktether,ltether,mbrad,rrad,p_tether,l_mem,k_mem,kmemb,kgly,kpep,fturgor,pturgor
   real(kind=8)::wthick,lsqueez,kwall,walrad,eta,rho,dxsol,dyzsol,ptread,phyd,tadd,kexcl,kpair,gap,kradial
   real(kind=8)::viscos,printtime,k_scale,k_scale_mb,mwid,rwid,kwid,kplane,kcir,kcircon,growstep,xgap,fgrad




16 format(A32)

   open(1,file='jobtype.txt')


   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='JOBID') exit
   end do

   read(1,*)jobid

   close(1)



   open(1,file='paras.txt')


   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='NTHREAD') exit
   end do

   read(1,*)nthreads

   do n=1,1000
      read(1,16)chara
      if(chara(1:2)=='PI') exit
   end do
   read(1,*)pi

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='DELTA') exit
   end do

   read(1,*)delta

   INVDELTA=1.0D0/DELTA

   do n=1,1000
      read(1,16)chara
      if(chara(1:4)=='BETA') exit
   end do

   read(1,*)beta

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='KFTSZ') exit
   end do

   read(1,*)k_z

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='LFTSZ') exit
   end do

   read(1,*)l_z

   l_diag=l_z*sqrt(2.0d0)

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='KZTHET') exit
   end do

   read(1,*)kz_thet

!  spring constant for the 4 chains of beads:

   kb=kz_thet/l_z**2

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='ZANGLE') exit
   end do

   read(1,*)angle

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='BENDIR') exit
   end do

   read(1,*)jbend

!  normalized changes in L_z when binding to GDP:

   dl_z=cos(pi*(180.0d0-angle)/180.0d0/2.0d0)

   l_z12=1.0d0+dl_z*jbend

   l_z34=1.0d0-dl_z*jbend

!  change in diagonal:

   l_diag_new=sqrt(2.0d0-dl_z**2)

!  for simplicity, FtsA = FtsZ in size. Then

   k_a=k_z

!  the distance from A center to each corner of the base:

   la_cen=l_z*sqrt(3.0d0)/2.0d0

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='LMBTOA') exit
   end do

   read(1,*)l_mb_a

!  tethering constant:

!   ktether=k_z*l_z/l_mb_a

   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='KTETHER') exit
   end do

   read(1,*)ktether

   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='LTETHER') exit
   end do

   read(1,*)ltether

   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='KRADIAL') exit
   end do

   read(1,*)kradial

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='KPLANE') exit
   end do

   read(1,*)kplane

   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='KCIRCUM') exit
   end do

   read(1,*)kcir

   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='KCIRCON') exit
   end do

   read(1,*)kcircon



!  the distance from A to Z:

   l_a_z=ltether-l_mb_a!8.0d0-l_z/2-l_z/2

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='MBRAD') exit
   end do

   read(1,*)mbrad

   do n=1,1000
      read(1,16)chara
      if(chara(1:4)=='MWID') exit
   end do

   read(1,*)mwid

   do n=1,1000
      read(1,16)chara
      if(chara(1:4)=='RWID') exit
   end do

   read(1,*)rwid

   do n=1,1000
      read(1,16)chara
      if(chara(1:4)=='KWID') exit
   end do

   read(1,*)kwid

!  Zring radius

   rrad=mbrad-ltether

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='NFTSZ') exit
   end do

   read(1,*)nftsz

   do n=1,1000
      read(1,16)chara
      if(chara(1:4)=='ZLEN') exit
   end do

   read(1,*)zlen

   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='NZSTART') exit
   end do

   read(1,*)nftszstart


   do n=1,1000
      read(1,16)chara
      if(chara(1:4)=='XGAP') exit
   end do

   read(1,*)xgap


   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='PTETHER') exit
   end do

   read(1,*)p_tether

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='LMEMB') exit
   end do

   read(1,*)l_mem

!  membrane rigidity

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='KMEMB') exit
   end do

   read(1,*)kmemb

   do n=1,1000
      read(1,16)chara
      if(chara(1:4)=='MEMB') exit
   end do

   read(1,*)k_mem
!   k_mem=k_z/10

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='WTHICK') exit
   end do

   read(1,*)wthick

   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='LSQUEEZ') exit
   end do

   read(1,*)lsqueez

!   do n=1,1000
!      read(1,16)chara
!      if(chara(1:5)=='KWALL') exit
!   end do

!   read(1,*)kwall

!  define cell wall radius = membrane radius + wall thickness

   walrad=mbrad+wthick

   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='VDWAALS') exit
   end do

   read(1,*)eta

!  smaller radius is used for lateral contact interaction
   rho=l_z/2

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='DXSOL') exit
   end do

   read(1,*)dxsol

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='DYZSOL') exit
   end do

   read(1,*)dyzsol

!  force constant mimicking glycan strands

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='GLYCAN') exit
   end do

   read(1,*)kgly

!  force constant mimicking peptide

   do n=1,1000
      read(1,16)chara
      if(chara(1:7)=='PEPTIDE') exit
   end do

   read(1,*)kpep

!  turgor pressure in atm

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='TURGOR') exit
   end do

   read(1,*)pturgor

   pturgor=pturgor/100

   fturgor=pturgor*dxsol*dyzsol

!  KWALL is calculated according to turgor pressure:

   kwall=pturgor*l_mem*l_mem/lsqueez/lsqueez

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='NTURN') exit
   end do

   read(1,*)nturnover

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='JTREAD') exit
   end do

   read(1,*)jtreadmil

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='PTREAD') exit
   end do

   read(1,*)ptread

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='PHYDRO') exit
   end do

   read(1,*)phyd

   do n=1,1000
      read(1,16)chara
      if(chara(1:4)=='TADD') exit
   end do

   read(1,*)tadd


!   do n=1,1000
!      read(1,16)chara
!      if(chara(1:4)=='ROFF') exit
!   end do

!   read(1,*)r_off

   do n=1,1000
      read(1,16)chara
      if(chara(1:4)=='KVDW') exit
   end do

   read(1,*)kexcl

!   do n=1,1000
!      read(1,16)chara
!      if(chara(1:3)=='EPS') exit
!   end do

!   read(1,*)eps

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='KPAIR') exit
   end do

   read(1,*)kpair

   do n=1,1000
      read(1,16)chara
      if(chara(1:4)=='WGAP') exit
   end do

   read(1,*)gap

   do n=1,1000
      read(1,16)chara
      if(chara(1:8)=='GROWSTEP') exit
   end do

   read(1,*)growstep

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='WGRAD') exit
   end do

   read(1,*)jgradient

   jnogradient=1-jgradient

!  gradient factor for the wall boundary

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='FGRAD') exit
   end do

   read(1,*)fgrad


!   do n=1,1000
!      read(1,16)chara
!      if(chara(1:6)=='JPRINT') exit
!   end do

!   read(1,*)jprint

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='PRINTT') exit
   end do

   read(1,*)printtime

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='NOUTPU') exit
   end do

   read(1,*)noutput


   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='VISCOS') exit
   end do

   read(1,*)viscos

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='NSTEP') exit
   end do

   read(1,*)nstep

   do n=1,1000
      read(1,16)chara
      if(chara(1:5)=='NRFOR') exit
   end do

   read(1,*)nrforce

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='KSCALE') exit
   end do

   read(1,*)k_scale

   do n=1,1000
      read(1,16)chara
      if(chara(1:6)=='KMBSCA') exit
   end do

   read(1,*)k_scale_mb

   close(1)

   end subroutine



!=========================================================

!==================================================

   subroutine makering(nmemb,nfil,nftsz,npoly,zlen,nwall,nbondwal,nxsol,nphisol,nftsa,nftsamax, &
               fil2a,fstart,flen,filid,a2mem,a2fil,bondwal,dxsol,dphisol,xboundmin,xboundmax, &
                rwall,lwall,xwall,ywall,zwall,xnorwall,ynorwall,znorwall,xsurf,ysurf,zsurf, &
                 xmemb,ymemb,zmemb,xcen,ycen,zcen,xa,ya,za,xb,yb,zb)

   implicit none

   integer,value:: nmemb,nfil,nftsz,npoly,zlen,nwall,nbondwal,nxsol,nphisol,nftsa,nftsamax
   integer nbond,n,natom,izero,ires,jatom,jx,jp,jstart,j,na,jm,ja,n0
   integer natom_z1,natom_z2,natom_z3,natom_z4,natom_zc,natom_a,natom_ab,natom_z2a
   integer nbond_z1,nbond_z2,nbond_z3,nbond_z4,nbond_zc,nbond_z2a,nbond_a2m,nbond0
   integer j1,j2,j3,j4,j5,j6,j7,j8,nline,i,jobid

   integer,allocatable,dimension(:),intent(in)::fil2a,fstart,flen,filid,a2mem,a2fil
   integer,allocatable,dimension(:,:),intent(in)::bondwal
   integer,allocatable,dimension(:,:)::bond
   integer,allocatable,dimension(:)::gtp

   real(kind=8)::charg,mass,w1,w2,zero
   real(kind=8),value::dxsol,dphisol,xboundmin,xboundmax

   real(kind=8),allocatable,dimension(:,:),intent(in)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall,lwall
   real(kind=8),allocatable,dimension(:,:),intent(in)::xsurf,ysurf,zsurf
   real(kind=8),allocatable,dimension(:),intent(in)::xmemb,ymemb,zmemb,rwall
   real(kind=8),allocatable,dimension(:),intent(in)::xcen,ycen,zcen,xa,ya,za!,x_a_old,y_a_old,z_a_old
   real(kind=8),allocatable,dimension(:,:),intent(in)::xb,yb,zb!,x_a_old,y_a_old,z_a_old

   character(8) tex,typ,res,segid,resno




20  FORMAT(I8,1X,A4,1X,A4,1X,A3,2X,A3,2X,A4,2X,F9.6,6X,F8.4,11X,I1)
21  FORMAT(I8,1X,A)

42 FORMAT(A6,I5,2X,A3,1X,A3,1X,I5,4X,3F8.3,2F6.2,6X,A4)
52 FORMAT(A50)

   open(1,file='zring.psf')

   open(2,file='startzring.pdb')

!  start with cell wall:
   nbond=nbondwal

!  FtsZ bonds:
   nbond=nbond+5*nftsz

!  to visualize tethering of FtsZ to membrane via FtsA:

   nbond=nbond+2*nftsamax

   allocate(bond(2,nbond))

!  calculate total number of atoms:

   natom=nwall+nmemb

!  marking filament junctions at 4 beads and their center:

   natom_z1=natom+nftsz*2
   nbond_z1=nbondwal+nftsz

   natom_z2=natom_z1+nftsz*2
   nbond_z2=nbond_z1+nftsz

   natom_z3=natom_z2+nftsz*2
   nbond_z3=nbond_z2+nftsz

   natom_z4=natom_z3+nftsz*2
   nbond_z4=nbond_z3+nftsz

   natom_zc=natom_z4+nftsz*2
   nbond_zc=nbond_z4+nftsz

!  adding FtsA:
   natom=natom_zc+nftsamax

   natom_a=natom

!  to visualize FtsA base surface:

!   natom=natom_a+4*nftsamax

!   natom_ab=natom


!  to visualize Z-to-A and A-to-membrane tethering:

   natom_z2a=natom+2*nftsamax
   nbond_z2a=nbond_zc+nftsamax
   nbond_a2m=nbond_z2a+nftsamax

   natom=natom+2*nftsamax+2*nftsamax


   write(1,21)natom,'!NATOM'

   charg=0.0d0
   mass=1.0d0
   izero=0
   tex='ATOM'
   w1=1.0d0
   w2=0.0d0

   zero=0.0

!  for cell wall

   res='WAL'
   typ='WAL'
   segid='WAL'

   ires=1
   write(resno,'(i1)')ires

   jatom=0

   do jx=1,nxsol

      do jp=1,nphisol

         jatom=jatom+1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xwall(jp,jx)/10,ywall(jp,jx)/10,zwall(jp,jx)/10,w1,w2,segid

      end do

   end do

   bond(1:2,1:nbondwal)=bondwal(1:2,1:nbondwal)

   nbond=nbondwal

!  for membrane:

   res='MBR'
   typ='MBR'
   segid='MBR'

   ires=2
   write(resno,'(i1)')ires

   do n=1,nmemb

      jatom=jatom+1

      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

      write(2,42)tex,jatom,typ,res,ires,xmemb(n)/10,ymemb(n)/10,zmemb(n)/10,w1,w2,segid

   end do


!  for FtsZ:

   segid='FTZ'

   typ='FTZ'

!  Z bead type 1:

   res='ZB1'

   ires=3
   write(resno,'(i1)')ires

   do n=1,nfil

      jstart=fstart(n)

      do j=1,flen(n)-1

         jatom=jatom+1

         na=jstart+j-1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xb(1,na)/10,yb(1,na)/10,zb(1,na)/10,w1,w2,segid

         jatom=jatom+1

         na=jstart+j

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xb(1,na)/10,yb(1,na)/10,zb(1,na)/10,w1,w2,segid

         nbond=nbond+1

         bond(1,nbond)=jatom-1

         bond(2,nbond)=jatom

      end do

   end do

   if(nbond<nbond_z1)then

      nbond0=nbond

      do n=nbond0+1,nbond_z1

         jatom=jatom+1
         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero
         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         jatom=jatom+1
         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero
         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         nbond=nbond+1

         bond(1,nbond)=jatom-1

         bond(2,nbond)=jatom

      end do

   end if



!  Z bead type 2:

   res='ZB2'

   ires=4
   write(resno,'(i1)')ires

   do n=1,nfil

      jstart=fstart(n)

      do j=1,flen(n)-1

         jatom=jatom+1

         na=jstart+j-1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xb(2,na)/10,yb(2,na)/10,zb(2,na)/10,w1,w2,segid

         jatom=jatom+1

         na=jstart+j

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xb(2,na)/10,yb(2,na)/10,zb(2,na)/10,w1,w2,segid

         nbond=nbond+1

         bond(1,nbond)=jatom-1

         bond(2,nbond)=jatom

      end do

   end do


   if(nbond<nbond_z2)then

      nbond0=nbond

      do n=nbond0+1,nbond_z2

         jatom=jatom+1
         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero
         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         jatom=jatom+1
         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero
         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         nbond=nbond+1

         bond(1,nbond)=jatom-1

         bond(2,nbond)=jatom
   
      end do

   end if


!  Z bead type 3:

   res='ZB3'

   ires=5
   write(resno,'(i1)')ires

   do n=1,nfil

      jstart=fstart(n)

      do j=1,flen(n)-1

         jatom=jatom+1

         na=jstart+j-1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xb(3,na)/10,yb(3,na)/10,zb(3,na)/10,w1,w2,segid

         jatom=jatom+1

         na=jstart+j

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xb(3,na)/10,yb(3,na)/10,zb(3,na)/10,w1,w2,segid

         nbond=nbond+1

         bond(1,nbond)=jatom-1

         bond(2,nbond)=jatom

      end do

   end do

   if(nbond<nbond_z3)then

      nbond0=nbond

      do n=nbond0+1,nbond_z3

         jatom=jatom+1
         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero
         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         jatom=jatom+1
         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero
         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         nbond=nbond+1

         bond(1,nbond)=jatom-1

         bond(2,nbond)=jatom
   
      end do

   end if


!  Z bead type 4:

   res='ZB4'

   ires=6
   write(resno,'(i1)')ires

   do n=1,nfil

      jstart=fstart(n)

      do j=1,flen(n)-1

         jatom=jatom+1

         na=jstart+j-1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xb(4,na)/10,yb(4,na)/10,zb(4,na)/10,w1,w2,segid

         jatom=jatom+1

         na=jstart+j

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xb(4,na)/10,yb(4,na)/10,zb(4,na)/10,w1,w2,segid

         nbond=nbond+1

         bond(1,nbond)=jatom-1

         bond(2,nbond)=jatom

      end do

   end do

   if(nbond<nbond_z4)then

      nbond0=nbond

      do n=nbond0+1,nbond_z4

         jatom=jatom+1
         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero
         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         jatom=jatom+1
         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero
         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         nbond=nbond+1

         bond(1,nbond)=jatom-1

         bond(2,nbond)=jatom
   
      end do

   end if


!  Z bead type center:

   res='ZBC'

   ires=7
   write(resno,'(i1)')ires

   do n=1,nfil

      jstart=fstart(n)

      do j=1,flen(n)-1

         jatom=jatom+1

         na=jstart+j-1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xcen(na)/10,ycen(na)/10,zcen(na)/10,w1,w2,segid

         jatom=jatom+1

         na=jstart+j

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,xcen(na)/10,ycen(na)/10,zcen(na)/10,w1,w2,segid

         nbond=nbond+1

         bond(1,nbond)=jatom-1

         bond(2,nbond)=jatom

      end do

   end do


   if(nbond<nbond_zc)then

      nbond0=nbond

      do n=nbond0+1,nbond_zc

         jatom=jatom+1
         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero
         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         jatom=jatom+1
         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero
         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         nbond=nbond+1

         bond(1,nbond)=jatom-1

         bond(2,nbond)=jatom
   
      end do

   end if


!  for FtsA:

   segid='FTA'

   res='AB1'
   typ='FTA'

   ires=8
   write(resno,'(i1)')ires


   do n=1,nftsa

      jatom=jatom+1

      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

      write(2,42)tex,jatom,typ,res,ires,xa(n)/10,ya(n)/10,za(n)/10,w1,w2,segid

   end do

   if(jatom<natom_a)then

      n0=jatom

      do n=n0+1,natom_a

         jatom=jatom+1

         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

      end do

   end if

!  for FtsA base:

!   res='AB2'

!   ires=9
!   write(resno,'(i1)')ires

!   do n=1,nftsa

!      do j=1,4

!         jatom=jatom+1

!         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

!         write(2,42)tex,jatom,typ,res,ires,x_a_old(j,n)/10,y_a_old(j,n)/10,z_a_old(j,n)/10,w1,w2,segid

!      end do

!   end do

!   if(jatom<natom_ab)then

!      n0=jatom

!      do n=n0+1,natom_ab

!         jatom=jatom+1

!         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

!         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

!      end do

!   end if


!  to visualize Z-A tethering:

   res='Z2A'
   typ='Z2A'
   segid='Z2A'

   ires=9
   write(resno,'(i1)')ires

   do ja=1,nftsa

      n=a2fil(ja)


      jatom=jatom+1

      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

!      write(2,42)tex,jatom,typ,res,ires,(xb(1,n)+xb(2,n))/20,(yb(1,n)+yb(2,n))/20,(zb(1,n)+zb(2,n))/20,w1,w2,segid
      write(2,42)tex,jatom,typ,res,ires,xcen(n)/10,ycen(n)/10,zcen(n)/10,w1,w2,segid

      jatom=jatom+1

      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

!      ja=fil2a(n)

!      write(2,42)tex,jatom,typ,res,ires,sum(x_a_old(1:4,ja))/40,sum(y_a_old(1:4,ja))/40,sum(z_a_old(1:4,ja))/40,w1,w2,segid

      write(2,42)tex,jatom,typ,res,ires,xa(ja)/10,ya(ja)/10,za(ja)/10,w1,w2,segid


      nbond=nbond+1

      bond(1,nbond)=jatom-1

      bond(2,nbond)=jatom

   end do

   if(nbond<nbond_z2a)then

      nbond0=nbond

      do n=nbond0+1,nbond_z2a

         jatom=jatom+1
         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero
         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         jatom=jatom+1
         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero
         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         nbond=nbond+1

         bond(1,nbond)=jatom-1

         bond(2,nbond)=jatom
   
      end do

   end if


!  visualize A-membrane tethering:

   res='A2M'
   typ='A2M'
   segid='A2M'

   ires=10
   write(resno,'(i2)')ires

   do n=1,nftsa

      jatom=jatom+1

      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

      write(2,42)tex,jatom,typ,res,ires,xa(n)/10,ya(n)/10,za(n)/10,w1,w2,segid

      jatom=jatom+1

      write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

      jm=a2mem(n)

      write(2,42)tex,jatom,typ,res,ires,xmemb(jm)/10,ymemb(jm)/10,zmemb(jm)/10,w1,w2,segid

      nbond=nbond+1

      bond(1,nbond)=jatom-1

      bond(2,nbond)=jatom

   end do

   if(nbond<nbond_a2m)then

      nbond0=nbond

      do n=nbond0+1,nbond_a2m

         jatom=jatom+1
         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero
         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         jatom=jatom+1
         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero
         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

         nbond=nbond+1

         bond(1,nbond)=jatom-1

         bond(2,nbond)=jatom

      end do

   end if

!   if(jatom<natom)then

!      n0=jatom

!      do n=n0+1,natom

!         jatom=jatom+1

!         write(1,20)jatom,segid,resno,res,typ,tex,charg,mass,izero

!         write(2,42)tex,jatom,typ,res,ires,zero,zero,zero,w1,w2,segid

!      end do

!   end if


   tex='END'

   write(2,52)tex

   close(2)

!-------------------------------------------

! -- WRITE LIST OF BONDS:

   NLINE=NBOND/4

   WRITE(1,*)
   WRITE(1,21)NBOND,'!NBOND: bonds'

   DO N=1,NLINE
      I=1+4*(N-1)
      J1=BOND(1,I); J2=BOND(2,I); J3=BOND(1,I+1); J4=BOND(2,I+1)
      J5=BOND(1,I+2); J6=BOND(2,I+2); J7=BOND(1,I+3); J8=BOND(2,I+3)

      WRITE(1,'(8I8)')J1,J2,J3,J4,J5,J6,J7,J8
   END DO

   IF(MOD(NBOND,4)==1)THEN
      J1=BOND(1,NBOND); J2=BOND(2,NBOND)
      WRITE(1,'(2I8)')J1,J2

   ELSE IF(MOD(NBOND,4)==2)THEN
      J1=BOND(1,NBOND-1); J2=BOND(2,NBOND-1)
      J3=BOND(1,NBOND); J4=BOND(2,NBOND)
      WRITE(1,'(4I8)')J1,J2,J3,J4

   ELSE IF(MOD(NBOND,4)==3)THEN
      J1=BOND(1,NBOND-2); J2=BOND(2,NBOND-2)
      J3=BOND(1,NBOND-1); J4=BOND(2,NBOND-1)
      J5=BOND(1,NBOND); J6=BOND(2,NBOND)
      WRITE(1,'(6I8)')J1,J2,J3,J4,J5,J6

   END IF

   CLOSE(1)

!---------------------------------------------

!  write coordinate

   open(1,file='coor000.inp',form='unformatted')

   write(1)dxsol,dphisol

   write(1)rwall(1:nxsol),lwall(1:nphisol,1:nxsol)
   write(1)xwall(1:nphisol,1:nxsol)
   write(1)ywall(1:nphisol,1:nxsol)
   write(1)zwall(1:nphisol,1:nxsol)
   write(1)xnorwall(1:nphisol,1:nxsol)
   write(1)ynorwall(1:nphisol,1:nxsol)
   write(1)znorwall(1:nphisol,1:nxsol)

   write(1)nmemb
   write(1)xmemb(1:nmemb)
   write(1)ymemb(1:nmemb)
   write(1)zmemb(1:nmemb)
   write(1)xboundmin,xboundmax
   write(1)xsurf(1:nphisol,1:nxsol)
   write(1)ysurf(1:nphisol,1:nxsol)
   write(1)zsurf(1:nphisol,1:nxsol)
   write(1)xnorwall(1:nphisol,1:nxsol)
   write(1)ynorwall(1:nphisol,1:nxsol)
   write(1)znorwall(1:nphisol,1:nxsol)

   write(1)npoly

   do n=1,npoly

      write(1)xb(1:4,n),xcen(n)
      write(1)yb(1:4,n),ycen(n)
      write(1)zb(1:4,n),zcen(n)

   end do

   write(1)nftsa

   write(1)xa(1:nftsa)
   write(1)ya(1:nftsa)
   write(1)za(1:nftsa)

!   do n=1,nftsa

!      write(1)x_a_old(1:4,n),xa(n)
!      write(1)y_a_old(1:4,n),ya(n)
!      write(1)z_a_old(1:4,n),za(n)
!   end do

   close(1)

!---------------------------------------------

!  write configuration

   open(1,file='conf000.inp')

   write(1,*)'VISUAL numbers'
   write(1,*)natom,natom_z1,natom_z2,natom_z3,natom_z4,natom_zc,natom_a,natom_z2a
   write(1,*)

   write(1,*)'CELLWALL'
   write(1,*)nwall,nxsol,nphisol
   write(1,*)

   write(1,*)'MEMBRANE'
   write(1,*)nmemb
   write(1,*)

   write(1,*)'FTSZ'

   write(1,*)npoly,nfil,nftsz,zlen
   write(1,*)filid(1:npoly)
   write(1,*)fil2a(1:npoly)

   allocate(gtp(npoly))
   gtp=1
   write(1,*)gtp(1:npoly)


   write(1,*)flen(1:nfil)
   write(1,*)fstart(1:nfil)

   write(1,*)'FTSA'

   write(1,*)nftsa,nftsamax
   write(1,*)a2mem(1:nftsa)
   write(1,*)a2fil(1:nftsa)


   close(1)

!----------------------------------
!----------------------------------
   open(1,file='restart.inp')

   write(1,'(a65)')'      nstart       jfile       jstart        time0        runtime'

   write(1,*)0,0,1,0,0.0
!   write(1,*)0,0.0
   close(1)


   jobid=2

   open(1,file='jobtype.txt')
   write(1,'(a79)')'to decide type of job: start from the beginning if jobid=1, continue if jobid=2'
   write(1,*)
   write(1,'(a5)')'JOBID'
   write(1,'(i1)')jobid
   close(1)

   end subroutine

!=========================================================

   subroutine getinfo(nstart,jfile,time0,jstart,runtime,natom,natom_z1,natom_z2,natom_z3,natom_z4,natom_zc, &
         natom_a,natom_z2a,nwall,nxsol,nphisol,nmemb,nftsz,npoly,nfil,zlen,nftsa,nftsamax)

   implicit none

   integer nstart,jfile
   integer(kind=8)::time0,jstart
   double precision runtime

   integer natom,natom_z1,natom_z2,natom_z3,natom_z4,natom_zc
   integer natom_z2a,natom_a,nwall,nxsol,nphisol,nmemb,nftsz,npoly,nfil,zlen,nftsa,nftsamax

   character zero*1,charid1*1,charid2*2,charid3*3
   character (len=128) fileconf,chara



   open(1,file='restart.inp')
   read(1,'(a65)')chara
   read(1,*)nstart,jfile,jstart,time0,runtime
   close(1)

   write(zero,'(i1)')0

   if(nstart<10)then
      write(charid1,'(i1)')nstart
      fileconf='conf'//zero//zero//charid1//'.inp'
   elseif(nstart<100)then
      write(charid2,'(i2)')nstart
      fileconf='conf'//zero//charid2//'.inp'
   else
      write(charid3,'(i3)')nstart
      fileconf='conf'//charid3//'.inp'
   end if

   open(1,file=fileconf)

22 read(1,*)chara

   if(chara(1:4)/='VISU')then
      goto 22
   end if

   read(1,*)natom,natom_z1,natom_z2,natom_z3,natom_z4,natom_zc,natom_a,natom_z2a

2  read(1,*)chara

   if(chara(1:4)/='CELL')then
      goto 2
   end if

   read(1,*)nwall,nxsol,nphisol

12 read(1,*)chara

   if(chara(1:4)/='MEMB')then
      goto 12
   end if

   read(1,*)nmemb

3  read(1,*)chara

   if(chara(1:4)/='FTSZ')then
      goto 3
   end if

   read(1,*)npoly,nfil,nftsz,zlen

4  read(1,*)chara

   if(chara(1:4)/='FTSA')then
      goto 4
   end if

   read(1,*)nftsa,nftsamax

   close(1)

   end subroutine

!=========================================================

   subroutine ringin(nstart,nmemb,nfil,npoly,nftsa,nphisol,nxsol, &
               filid,fil2a,flen,fstart,a2mem,a2fil,endfil,gtp, &
                dxsol,dphisol,xboundmin,xboundmax,rwall,lwall,xwall,ywall,zwall,xnorwall,ynorwall,znorwall, &
                 xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf,xmemb,ymemb,zmemb, &
                  xcen,ycen,zcen,xa,ya,za,xb,yb,zb)

   implicit none

   integer,value::nstart,nmemb,npoly,nftsa,nphisol,nxsol,nfil
   integer n1,n2,n,jread(4)

   integer,allocatable,dimension(:)::filid,fil2a,flen,fstart,a2mem,a2fil,endfil,gtp

   real(kind=8)::dxsol,dphisol,xboundmin,xboundmax

   real(kind=8),allocatable,dimension(:,:)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall,lwall
   real(kind=8),allocatable,dimension(:,:)::xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf
   real(kind=8),allocatable,dimension(:)::xmemb,ymemb,zmemb,rwall
   real(kind=8),allocatable,dimension(:,:)::xb,yb,zb!,x_a_old,y_a_old,z_a_old
   real(kind=8),allocatable,dimension(:)::xcen,ycen,zcen,xa,ya,za!,x_a_old,y_a_old,z_a_old

   character zero*1,charid1*1,charid2*2,charid3*3
   character (len=64) filecoor,fileconf,chara

   write(zero,'(i1)')0

   if(nstart<10)then
      write(charid1,'(i1)')nstart
      fileconf='conf'//zero//zero//charid1//'.inp'
      filecoor='coor'//zero//zero//charid1//'.inp'
   elseif(nstart<100)then
      write(charid2,'(i2)')nstart
      fileconf='conf'//zero//charid2//'.inp'
      filecoor='coor'//zero//charid2//'.inp'
   else
      write(charid3,'(i3)')nstart
      fileconf='conf'//charid3//'.inp'
      filecoor='coor'//charid3//'.inp'
   end if

!  read coordinates

   open(1,file=filecoor,form='unformatted')

   read(1)dxsol,dphisol
   read(1)rwall(1:nxsol),lwall(1:nphisol,1:nxsol)
   read(1)xwall(1:nphisol,1:nxsol)
   read(1)ywall(1:nphisol,1:nxsol)
   read(1)zwall(1:nphisol,1:nxsol)
   read(1)xnorwall(1:nphisol,1:nxsol)
   read(1)ynorwall(1:nphisol,1:nxsol)
   read(1)znorwall(1:nphisol,1:nxsol)

   read(1)jread(1)
   read(1)xmemb(1:nmemb)
   read(1)ymemb(1:nmemb)
   read(1)zmemb(1:nmemb)
   read(1)xboundmin,xboundmax
   read(1)xsurf(1:nphisol,1:nxsol)
   read(1)ysurf(1:nphisol,1:nxsol)
   read(1)zsurf(1:nphisol,1:nxsol)
   read(1)xnorsurf(1:nphisol,1:nxsol)
   read(1)ynorsurf(1:nphisol,1:nxsol)
   read(1)znorsurf(1:nphisol,1:nxsol)

   read(1)jread(1)

   do n=1,npoly

      read(1)xb(1:4,n),xcen(n)
      read(1)yb(1:4,n),ycen(n)
      read(1)zb(1:4,n),zcen(n)

   end do

   read(1)jread(1)

   read(1)xa(1:nftsa)
   read(1)ya(1:nftsa)
   read(1)za(1:nftsa)

!   do n=1,nftsa

!      read(1)x_a_old(1:4,n),xa(n)
!      read(1)y_a_old(1:4,n),ya(n)
!      read(1)z_a_old(1:4,n),za(n)
!   end do

   close(1)


!  read configuration

   open(1,file=fileconf)

3  read(1,*)chara

   if(chara(1:4)/='FTSZ')then
      goto 3
   end if

   read(1,*)jread(1:4)
   read(1,*)filid(1:npoly)
   read(1,*)fil2a(1:npoly)
   read(1,*)gtp(1:npoly)
   read(1,*)flen(1:nfil)
   read(1,*)fstart(1:nfil)

   n=0

!   startfil=0

   endfil=0

   do n1=1,nfil

!      startfil(fstart(n1))=1

      do n2=1,flen(n1)

         n=n+1


         if(n2==flen(n1)) endfil(n)=1

      end do

   end do

4  read(1,*)chara

   if(chara(1:4)/='FTSA')then
      goto 4
   end if

   read(1,*)jread(2)
   read(1,*)a2mem(1:nftsa)
   read(1,*)a2fil(1:nftsa)


   close(1)

   end subroutine

!=========================================================

   subroutine dcdheader(junit,jfile,ntotal)

   implicit none

   character coor*4,zero*1,charid1*1,charid2*2,charid3*3
   character (len=80) string1,string2,filedcd
   integer ifirst,nframe,nfreq,zeros5(5),peroff,zeros7(7),two,twentyfour,ntot
   integer,value:: junit,jfile,ntotal
   integer*8 jdelta

   write(zero,'(i1)')0

   if(jfile<10)then
      write(charid1,'(i1)')jfile
      filedcd='traj'//zero//zero//charid1//'.dcd'
   elseif(jfile<100)then
      write(charid2,'(i2)')jfile
      filedcd='traj'//zero//charid2//'.dcd'
   elseif(jfile<1000)then
      write(charid3,'(i3)')jfile
      filedcd='traj'//charid3//'.dcd'
   else
      print*,'too many dcd file already, stopping now...'
      stop
   end if

   open(junit,file=filedcd,form='unformatted')

   COOR='CORD'; NFRAME=10000; IFIRST=0; NFREQ=1; NTOT=100
   ZEROS5=0; JDELTA=1; PEROFF=0; ZEROS7=0; TWENTYFOUR=24; TWO=2
   STRING1='HELLOOOOOOO'; STRING2='WHAT THE HELL!'

   WRITE(JUNIT)COOR,NFRAME,IFIRST,NFREQ,NTOT,ZEROS5,JDELTA,PEROFF,ZEROS7,TWENTYFOUR
   WRITE(JUNIT)TWO,STRING1,STRING2
   WRITE(JUNIT)NTOTAL

   end subroutine

!=========================================================

   subroutine writedcd(junit,nframe,natom,natom_z1,natom_z2,natom_z3,natom_z4,natom_zc, &
               natom_a,natom_z2a,nxsol,nphisol,nmemb,nftsz,nfil,nftsa, &
                flen,a2mem,a2fil,fstart,xb,yb,zb,xcen,ycen,zcen,xa,ya,za, &
                 xmemb,ymemb,zmemb,xwall,ywall,zwall)

   implicit none

   integer,value::junit,natom,natom_z1,natom_z2,natom_z3,natom_z4,natom_zc,natom_a,natom_z2a
   integer,value::nxsol,nphisol,nmemb,nftsz,nfil,nftsa
   integer nframe
   integer jw,jx,jp,n,j,jm,jstart,ja
   integer,allocatable,intent(in),dimension(:)::flen,a2mem,a2fil,fstart

   real,allocatable,dimension(:)::xw,yw,zw
   real(kind=8),allocatable,intent(in),dimension(:,:)::xb,yb,zb!,x_a_old,y_a_old,z_a_old
   real(kind=8),allocatable,intent(in),dimension(:)::xcen,ycen,zcen,xa,ya,za!,x_a_old,y_a_old,z_a_old
   real(kind=8),allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb
   real(kind=8),allocatable,intent(in),dimension(:,:)::xwall,ywall,zwall

!print*,'natom',natom_z1,natom_z2,natom_z3,natom_z4,natom_zc,natom_a,natom_z2a,natom

   nframe=nframe+1

   allocate(xw(natom),yw(natom),zw(natom))

   jw=0

   do jx=1,nxsol

      do jp=1,nphisol

         jw=jw+1

         xw(jw)=0.1*xwall(jp,jx)
         yw(jw)=0.1*ywall(jp,jx)
         zw(jw)=0.1*zwall(jp,jx)

      end do

   end do

   xw(jw+1:jw+nmemb)=0.1*xmemb(1:nmemb)
   yw(jw+1:jw+nmemb)=0.1*ymemb(1:nmemb)
   zw(jw+1:jw+nmemb)=0.1*zmemb(1:nmemb)

   jw=jw+nmemb

!  visualize FtsZ:


   do n=1,nfil

      jstart=fstart(n)

      do j=1,flen(n)-1

         jw=jw+1

         xw(jw)=0.1*xb(1,jstart+j-1)
         yw(jw)=0.1*yb(1,jstart+j-1)
         zw(jw)=0.1*zb(1,jstart+j-1)

         jw=jw+1

         xw(jw)=0.1*xb(1,jstart+j)
         yw(jw)=0.1*yb(1,jstart+j)
         zw(jw)=0.1*zb(1,jstart+j)

      end do

   end do

   if(jw<natom_z1)then

      xw(jw+1:natom_z1)=xw(jw)
      yw(jw+1:natom_z1)=yw(jw)
      zw(jw+1:natom_z1)=zw(jw)

      jw=natom_z1

   end if

   do n=1,nfil

      jstart=fstart(n)

      do j=1,flen(n)-1

         jw=jw+1

         xw(jw)=0.1*xb(2,jstart+j-1)
         yw(jw)=0.1*yb(2,jstart+j-1)
         zw(jw)=0.1*zb(2,jstart+j-1)

         jw=jw+1

         xw(jw)=0.1*xb(2,jstart+j)
         yw(jw)=0.1*yb(2,jstart+j)
         zw(jw)=0.1*zb(2,jstart+j)

      end do

   end do

   if(jw<natom_z2)then

      xw(jw+1:natom_z2)=xw(jw)
      yw(jw+1:natom_z2)=yw(jw)
      zw(jw+1:natom_z2)=zw(jw)

      jw=natom_z2

   end if

   do n=1,nfil

      jstart=fstart(n)

      do j=1,flen(n)-1

         jw=jw+1

         xw(jw)=0.1*xb(3,jstart+j-1)
         yw(jw)=0.1*yb(3,jstart+j-1)
         zw(jw)=0.1*zb(3,jstart+j-1)

         jw=jw+1

         xw(jw)=0.1*xb(3,jstart+j)
         yw(jw)=0.1*yb(3,jstart+j)
         zw(jw)=0.1*zb(3,jstart+j)

      end do

   end do

   if(jw<natom_z3)then

      xw(jw+1:natom_z3)=xw(jw)
      yw(jw+1:natom_z3)=yw(jw)
      zw(jw+1:natom_z3)=zw(jw)

      jw=natom_z3

   end if

   do n=1,nfil

      jstart=fstart(n)

      do j=1,flen(n)-1

         jw=jw+1

         xw(jw)=0.1*xb(4,jstart+j-1)
         yw(jw)=0.1*yb(4,jstart+j-1)
         zw(jw)=0.1*zb(4,jstart+j-1)

         jw=jw+1

         xw(jw)=0.1*xb(4,jstart+j)
         yw(jw)=0.1*yb(4,jstart+j)
         zw(jw)=0.1*zb(4,jstart+j)

      end do

   end do

   if(jw<natom_z4)then

      xw(jw+1:natom_z4)=xw(jw)
      yw(jw+1:natom_z4)=yw(jw)
      zw(jw+1:natom_z4)=zw(jw)

      jw=natom_z4

   end if

   do n=1,nfil

      jstart=fstart(n)

      do j=1,flen(n)-1

         jw=jw+1

         xw(jw)=0.1*xcen(jstart+j-1)
         yw(jw)=0.1*ycen(jstart+j-1)
         zw(jw)=0.1*zcen(jstart+j-1)

         jw=jw+1

         xw(jw)=0.1*xcen(jstart+j)
         yw(jw)=0.1*ycen(jstart+j)
         zw(jw)=0.1*zcen(jstart+j)

      end do

   end do

   if(jw<natom_zc)then

      xw(jw+1:natom_zc)=xw(jw)
      yw(jw+1:natom_zc)=yw(jw)
      zw(jw+1:natom_zc)=zw(jw)

      jw=natom_zc

   end if

!  for FtsA:

   xw(jw+1:jw+nftsa)=0.1*xa(1:nftsa)
   yw(jw+1:jw+nftsa)=0.1*ya(1:nftsa)
   zw(jw+1:jw+nftsa)=0.1*za(1:nftsa)

   jw=jw+nftsa

   if(jw<natom_a)then

      xw(jw+1:natom_a)=xw(jw)
      yw(jw+1:natom_a)=yw(jw)
      zw(jw+1:natom_a)=zw(jw)

      jw=natom_a

   end if




!  Z-to-A tethering

   do ja=1,nftsa

      n=a2fil(ja)


         jw=jw+1


         xw(jw)=0.1*xcen(n)
         yw(jw)=0.1*ycen(n)
         zw(jw)=0.1*zcen(n)


         jw=jw+1


         xw(jw)=0.1*xa(ja)
         yw(jw)=0.1*ya(ja)
         zw(jw)=0.1*za(ja)


   end do

   if(jw<natom_z2a)then

      xw(jw+1:natom_z2a)=xw(jw)
      yw(jw+1:natom_z2a)=yw(jw)
      zw(jw+1:natom_z2a)=zw(jw)

      jw=natom_z2a

   end if

!  for A-2-membrane tethering:

   do n=1,nftsa

      jw=jw+1

      xw(jw)=0.1*xa(n)
      yw(jw)=0.1*ya(n)
      zw(jw)=0.1*za(n)

      jm=a2mem(n)

      jw=jw+1

      xw(jw)=0.1*xmemb(jm)
      yw(jw)=0.1*ymemb(jm)
      zw(jw)=0.1*zmemb(jm)

   end do

   if(jw<natom)then

      xw(jw+1:natom)=xw(jw)
      yw(jw+1:natom)=yw(jw)
      zw(jw+1:natom)=zw(jw)

      jw=natom

   end if

   write(junit)xw(1:natom)
   write(junit)yw(1:natom)
   write(junit)zw(1:natom)

   deallocate(xw,yw,zw)

   end subroutine

!=========================================================

   subroutine solidset(nxsol,nphisol,nmemb,npoly,jsursol,jmbsol,jzsol,pi,delta, &
              dxsol,dphisol,xmemb,ymemb,zmemb,xcen,ycen,zcen,xwall,ywall,zwall, &
              xnorwall,ynorwall,znorwall,xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf)

   implicit none

   integer,value:: nxsol,nphisol,nmemb,npoly
   integer n,jx,jp,jx0,jp0,j1,j2,jxget,jpget

   integer,allocatable,dimension(:,:)::jsursol,jmbsol,jzsol

   real(kind=8),value::pi,delta,dxsol,dphisol
   real(kind=8)::xmin,dxsolby2,piby2,dphisolby2,twopi,phi,arg
   real(kind=8)::dx,dy,dz,d2,dist2

   real(kind=8),allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb
   real(kind=8),allocatable,intent(in),dimension(:)::xcen,ycen,zcen
   real(kind=8),allocatable,intent(in),dimension(:,:)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall
   real(kind=8),allocatable,intent(in),dimension(:,:)::xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf


   xmin=-(nxsol-1)/2*dxsol

   dxsolby2=0.5d0*dxsol

   piby2=0.5d0*pi

   dphisolby2=0.5d0*dphisol

   twopi=2*pi


!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,phi,arg,jp,jx0,jp0,j1,j2,jxget,jpget,dx,dy,dz,d2,dist2) &
!$omp shared(nmemb,xmemb,xmin,dxsol,dxsolby2,nxsol,jmbsol,ymemb,delta,zmemb,piby2) &
!$omp shared(dphisolby2,nphisol,pi,twopi,dphisol) &
!$omp shared(xwall,ywall,zwall,xnorwall,ynorwall,znorwall,jsursol) &
!$omp shared(xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf) &
!$omp shared(npoly,xcen,ycen,zcen,jzsol)

!$omp do schedule(guided,64)

   do n=1,nmemb

      jx=1+(xmemb(n)-xmin)/dxsol

      if(xmemb(n)-xmin-(jx-1)*dxsol>dxsolby2)then
         jx=jx+1
      end if

      if(jx<1)then
         jx=1
      end if

      if(jx>nxsol)then
         jx=nxsol
      end if

      jsursol(1,n)=jx

      jx0=jx


      if(abs(ymemb(n))<delta)then
         if(zmemb(n)>0.0d0)then
            phi=piby2
         else
            phi=-piby2
         end if

      else

         arg=zmemb(n)/ymemb(n)

         phi=atan(arg)

         if(ymemb(n)<0.0d0)then
            phi=phi+pi
         end if

         if(arg<0.0d0.and.ymemb(n)>0.0d0)then
            phi=phi+twopi
         end if


      end if

      if(phi<dphisolby2)then
         jp=nphisol
      else

         jp=phi/dphisol

         if(phi-jp*dphisol>dphisolby2)then
            jp=jp+1
         end if

      end if

      jsursol(2,n)=jp


      jp0=jp

      dist2=1000000.0d0

      jxget=0

      do j1=1,20

         jx=jx0+j1-10

         if(jx<1) cycle

         if(jx>nxsol) cycle

         do j2=1,20

            jp=jp0+j2-10

            if(jp<1) jp=jp+nphisol

            if(jp>nphisol) jp=jp-nphisol

            if(ymemb(n)*ywall(jp,jx)+zmemb(n)*zwall(jp,jx)<0.0) cycle

            dx=xmemb(n)-xwall(jp,jx)
            dy=ymemb(n)-ywall(jp,jx)
            dz=zmemb(n)-zwall(jp,jx)

            arg=dx*xnorwall(jp,jx)+dy*ynorwall(jp,jx)+dz*znorwall(jp,jx)

            if(arg<0.0d0)cycle

            dx=dx-arg*xnorwall(jp,jx)
            dy=dy-arg*ynorwall(jp,jx)
            dz=dz-arg*znorwall(jp,jx)

            d2=dx*dx+dy*dy+dz*dz

            if(dist2>d2)then

               dist2=d2

               jxget=jx

               jpget=jp

            end if

         end do

      end do

      if(jxget==0)then
         print*,'error: could not get indices',n,jx0,jp0
         stop
      end if

      jmbsol(1,n)=jxget

      jmbsol(2,n)=jpget


   end do


!$omp end do nowait

!$omp do

   do n=1,npoly

      jx=1+(xcen(n)-xmin)/dxsol

      if(xcen(n)-xmin-(jx-1)*dxsol>dxsolby2)then
         jx=jx+1
      end if

      if(jx<1)then
         jx=1
      end if

      if(jx>nxsol)then
         jx=nxsol
      end if

      jx0=jx

      if(abs(ycen(n))<delta)then
         if(zcen(n)>0.0d0)then
            phi=piby2
         else
            phi=-piby2
         end if

      else

         arg=zcen(n)/ycen(n)

         phi=atan(arg)

         if(ycen(n)<0.0d0)then
            phi=phi+pi
         end if

         if(arg<0.0d0.and.ycen(n)>0.0d0)then
            phi=phi+twopi
         end if

      end if

      if(phi<dphisolby2)then
         jp=nphisol
      else

         jp=phi/dphisol

         if(phi-jp*dphisol>dphisolby2)then
            jp=jp+1
         end if

      end if

      jp0=jp

      dist2=1000000.0d0

      jxget=0

      do j1=1,20

         jx=jx0+j1-10

         if(jx<1) cycle

         if(jx>nxsol) cycle

         do j2=1,20

            jp=jp0+j2-10

            if(jp<1) jp=jp+nphisol

            if(jp>nphisol) jp=jp-nphisol

            if(ycen(n)*ysurf(jp,jx)+zcen(n)*zsurf(jp,jx)<0.0) cycle

            dx=xcen(n)-xsurf(jp,jx)
            dy=ycen(n)-ysurf(jp,jx)
            dz=zcen(n)-zsurf(jp,jx)

            arg=dx*xnorsurf(jp,jx)+dy*ynorsurf(jp,jx)+dz*znorsurf(jp,jx)

            if(arg<0.0d0)cycle

            dx=dx-arg*xnorsurf(jp,jx)
            dy=dy-arg*ynorsurf(jp,jx)
            dz=dz-arg*znorsurf(jp,jx)

            d2=dx*dx+dy*dy+dz*dz

            if(dist2>d2)then

               dist2=d2

               jxget=jx

               jpget=jp

            end if

         end do

      end do

      if(jxget==0)then
         print*,'error: could not get indices',n,jx0,jp0
         stop
      end if

      jzsol(1,n)=jxget

      jzsol(2,n)=jpget


   end do

!$omp end do nowait

!$omp end parallel



   end subroutine


!=========================================================
   subroutine rforceset(jrforce,nrforce,nmemb,nftsz,nftsamax,pi,rxmemb,rymemb,rzmemb, &
              rxbead,rybead,rzbead,rx_a,ry_a,rz_a,k_scale,k_scale_mb)

   implicit none

   integer,value::nrforce,nmemb,nftsz,nftsamax
   integer jrforce,j,nmax,n,nh1,nh2

   real(kind=8),value::pi,k_scale,k_scale_mb
   real(kind=8),allocatable,dimension(:,:)::rxbead,rybead,rzbead,rx_a,ry_a,rz_a
   real(kind=8),allocatable,dimension(:,:)::rxmemb,rymemb,rzmemb
   real(kind=8),allocatable,dimension(:)::rh1,rh2,rx,ry,rz
   integer seed,iseed,omp_get_thread_num

   nmax=nmemb+nftsz+nftsamax
   nh1=(nmax+1)/2
   nh2=nmax-nh1

   jrforce=0

   allocate(rh1(nh1),rh2(nh1),rx(nmax),ry(nmax),rz(nmax))

   CALL SYSTEM_CLOCK(COUNT=seed)

!$omp parallel &
!$omp default(none) &
!$omp private(j,rh1,rh2,rx,ry,rz,n,iseed) &
!$omp shared(nrforce,nmax,nh1,nh2,nmemb,nftsz,nftsamax,pi,k_scale,k_scale_mb) &
!$omp shared(rxbead,rybead,rzbead,rxmemb,rymemb,rzmemb,rx_a,ry_a,rz_a,seed) 

   iseed=seed+omp_get_thread_num( )

!$omp do schedule(guided,32)

   do j=1,nrforce

      call r4vec_uniform_01 ( nh1, iseed, rh1 )
      call r4vec_uniform_01 ( nh1, iseed, rh2 )

      rx(1:nh1)=sqrt(-2*log(rh1(1:nh1)))*sin(2*pi*rh2(1:nh1))
      rx(nh1+1:nmax)=sqrt(-2*log(rh1(1:nh2)))*cos(2*pi*rh2(1:nh2))


      call r4vec_uniform_01 ( nh1, iseed, rh1 )
      call r4vec_uniform_01 ( nh1, iseed, rh2 )

      ry(1:nh1)=sqrt(-2*log(rh1(1:nh1)))*sin(2*pi*rh2(1:nh1))
      ry(nh1+1:nmax)=sqrt(-2*log(rh1(1:nh2)))*cos(2*pi*rh2(1:nh2))

      call r4vec_uniform_01 ( nh1, iseed, rh1 )
      call r4vec_uniform_01 ( nh1, iseed, rh2 )

      rz(1:nh1)=sqrt(-2*log(rh1(1:nh1)))*sin(2*pi*rh2(1:nh1))
      rz(nh1+1:nmax)=sqrt(-2*log(rh1(1:nh2)))*cos(2*pi*rh2(1:nh2))

      rxmemb(1:nmemb,j)=rx(1:nmemb)*k_scale_mb
      rymemb(1:nmemb,j)=ry(1:nmemb)*k_scale_mb
      rzmemb(1:nmemb,j)=rz(1:nmemb)*k_scale_mb

      n=nmemb

      rxbead(1:nftsz,j)=k_scale*rx(n+1:n+nftsz)
      rybead(1:nftsz,j)=k_scale*ry(n+1:n+nftsz)
      rzbead(1:nftsz,j)=k_scale*rz(n+1:n+nftsz)

      n=n+nftsz

      rx_a(1:nftsamax,j)=k_scale*rx(n+1:n+nftsamax)
      ry_a(1:nftsamax,j)=k_scale*ry(n+1:n+nftsamax)
      rz_a(1:nftsamax,j)=k_scale*rz(n+1:n+nftsamax)

   end do

!$omp enddo nowait
!$omp end parallel


   deallocate(rh1,rh2)


   end subroutine

!=========================================================

subroutine r4vec_uniform_01 ( n, seed, r )

!*****************************************************************************80
!
!! R4VEC_UNIFORM_01 returns a unit pseudorandom R4VEC.
!
!  Discussion:
!
!    An R4VEC is an array of real ( kind = 4 ) values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, Number 2, 1969, pages 136-143.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of entries in the vector.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value,
!    which should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 4 ) R(N), the vector of pseudorandom values.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ), parameter :: i4_huge = 2147483647
  integer ( kind = 4 ) k
  integer ( kind = 4 ) seed
  real ( kind = 8 ) r(n)

  if ( seed == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R4VEC_UNIFORM_01 - Fatal error!'
    write ( *, '(a)' ) '  Input value of SEED = 0.'
    stop 1
  end if

  do i = 1, n

    k = seed / 127773

    seed = 16807 * ( seed - k * 127773 ) - k * 2836

    if ( seed < 0 ) then
      seed = seed + i4_huge
    end if

    r(i) = real ( seed, kind = 8 ) * 4.656612875E-10

  end do

  return

  end subroutine

!=========================================================

   subroutine treadmill(jtreadmil,nfil,nftsa,nmemb,fstart,flen,fil2a,a2fil,a2mem,gtp,jzsol, &
                        ptread,phyd,dttread,xcen,ycen,zcen,xa,ya,za,xb,yb,zb,xmemb,ymemb,zmemb)

   implicit none

   integer,value::nfil,nftsa,nmemb,jtreadmil

   integer n,jstart,length,ja_end,j,jz,ja,jpick,jm,jend,jcheck

   integer,allocatable,dimension(:),intent(in)::fstart,flen
   integer,allocatable,dimension(:)::fil2a,a2fil,a2mem,gtp
   integer,allocatable,dimension(:,:)::jzsol
   integer,allocatable,dimension(:)::mark

   real(kind=8),value::ptread,phyd,dttread

   real(kind=8)::r,dx,dy,dz,d2,dist2

   real(kind=8),allocatable,dimension(:)::xcen,ycen,zcen,xa,ya,za!,x_a_old,y_a_old,z_a_old
   real(kind=8),allocatable,dimension(:,:)::xb,yb,zb!,x_a_old,y_a_old,z_a_old
   real(kind=8),allocatable,dimension(:),intent(in)::xmemb,ymemb,zmemb

   allocate(mark(nmemb))

   mark=0

   mark(a2mem(1:nftsa))=1

   do n=1,nfil

      jstart=fstart(n)

      length=flen(n)

!      jcheck=0

!      do j=1,length

!         jz=jstart+j-1

!         ja=fil2a(jz)

!         if(ja>0) jcheck=jcheck+1

!         if(jcheck==2)then
!            jend=j+1
!            exit
!         end if

!      end do

!      jend=max(jend,length/3)

      do j=length,1,-1

         jz=jstart+j-1

         if(gtp(jz)==0) cycle

         call random_number(r)

!         phyd=ptread*j/(length-j+1)

         if(phyd*dttread>r) gtp(jz)=0

         exit

      end do

      if(jtreadmil==0) cycle

      call random_number(r)

      if(ptread*dttread<r) cycle



      ja_end=fil2a(jstart+length-1)

      do j=length,2,-1

         jz=jstart+j-1

         xb(1:4,jz)=xb(1:4,jz-1)
         yb(1:4,jz)=yb(1:4,jz-1)
         zb(1:4,jz)=zb(1:4,jz-1)

         xcen(jz)=xcen(jz-1)
         ycen(jz)=ycen(jz-1)
         zcen(jz)=zcen(jz-1)

         ja=fil2a(jz-1)

         fil2a(jz)=ja

         if(ja>0) a2fil(ja)=jz

         jzsol(1:2,jz)=jzsol(1:2,jz-1)

         gtp(jz)=gtp(jz-1)

      end do


!     add a new unit at the filament tip:
!     Mar 12, 2018: add new units to the circumferential direction.
!     Mar 14, 2018: reverse back to the previous alignment of the tip

      xb(1:4,jstart)=2*xb(1:4,jstart+1)-xb(1:4,jstart+2)
!      xb(1:4,jstart)=xb(1:4,jstart+1)
      yb(1:4,jstart)=2*yb(1:4,jstart+1)-yb(1:4,jstart+2)
      zb(1:4,jstart)=2*zb(1:4,jstart+1)-zb(1:4,jstart+2)

      xcen(jstart)=2*xcen(jstart+1)-xcen(jstart+2)
!      xcen(jstart)=xcen(jstart+1)
      ycen(jstart)=2*ycen(jstart+1)-ycen(jstart+2)
      zcen(jstart)=2*zcen(jstart+1)-zcen(jstart+2)

      if(ja_end>0)then

         fil2a(jstart)=ja_end

         a2fil(ja_end)=jstart

         mark(a2mem(ja_end))=0

         dist2=100000.0d0

         jpick=0

         do jm=1,nmemb

            if(mark(jm)==1) cycle

            dx=xcen(jstart)-xmemb(jm)
            dy=ycen(jstart)-ymemb(jm)
            dz=zcen(jstart)-zmemb(jm)

            d2=dx*dx+dy*dy+dz*dz

            if(d2<dist2)then

               dist2=d2

               jpick=jm

            end if

         end do

         if(jpick==0) print*,'ERROR: jpick for a2mem not assigned'

         a2mem(ja_end)=jpick

         mark(jpick)=1

         xa(ja_end)=0.5d0*(xcen(jstart)+xmemb(jpick))
         ya(ja_end)=0.5d0*(ycen(jstart)+ymemb(jpick))
         za(ja_end)=0.5d0*(zcen(jstart)+zmemb(jpick))


      else

         fil2a(jstart)=0

      end if

!     note: jzsol for the tip is temporarily unchanged


   end do


   end subroutine

!=========================================================

   subroutine depoly(nfil,npoly,nftsa,flen,fstart,filid,fil2a,a2mem,a2fil,startfil,endfil, &
              jzsol,ptread,dt,xcen,ycen,zcen,xa,ya,za,xb,yb,zb)

   implicit none

   integer,value::nfil
   integer npoly,nftsa
   integer check,n,length,jstart,jend,npar,j,jz,ja

   integer,allocatable,dimension(:)::flen,fstart,filid,fil2a,a2mem,a2fil,startfil,endfil
   integer,allocatable,dimension(:,:)::jzsol

   real(kind=8)::r,dx,dy,dz

   real(kind=8),value::ptread,dt

   real(kind=8),allocatable,dimension(:)::xcen,ycen,zcen,xa,ya,za!,x_a_old,y_a_old,z_a_old
   real(kind=8),allocatable,dimension(:,:)::xb,yb,zb!,x_a_old,y_a_old,z_a_old


   check=0

   do n=1,nfil

      length=flen(n)


      call random_number(r)

      if(ptread*dt>r)then

         check=check+1

         flen(n)=length-1

      end if

   end do

   if(check==0)then
      return
   end if

!  cleaning up the system:

   npoly=0; nftsa=0

   startfil=0

   endfil=0

   do n=1,nfil


      jstart=fstart(n)

      fstart(n)=npoly+1

      startfil(npoly+1)=1

      do j=1,flen(n)

         jz=jstart+j-1

         npoly=npoly+1

!         ringid(nftsz)=ringid(jz)

         filid(npoly)=n

         ja=fil2a(jz)

         if(ja>0)then

            nftsa=nftsa+1

            a2mem(nftsa)=a2mem(ja)

!            jasol(1:2,nftsa)=jasol(1:2,ja)

            xa(nftsa)=xa(ja)
            ya(nftsa)=ya(ja)
            za(nftsa)=za(ja)

!            x_a_old(1:4,nftsa)=x_a_old(1:4,ja)

!            y_a_old(1:4,nftsa)=y_a_old(1:4,ja)

!            z_a_old(1:4,nftsa)=z_a_old(1:4,ja)

            fil2a(npoly)=nftsa

            a2fil(nftsa)=npoly

         else

            fil2a(npoly)=0


         end if

         jzsol(1:2,npoly)=jzsol(1:2,jz)

         xcen(npoly)=xcen(jz)

         ycen(npoly)=ycen(jz)

         zcen(npoly)=zcen(jz)

         xb(1:4,npoly)=xb(1:4,jz)
         yb(1:4,npoly)=yb(1:4,jz)
         zb(1:4,npoly)=zb(1:4,jz)


         if(j==flen(n)) endfil(npoly)=1


      end do

   end do

   end subroutine

!=========================================================

   subroutine addfil(nmemb,nxsol,nphisol,zlen,nmono,nftsa,nfil,npoly,fstart,flen,fil2a,a2mem,a2fil, &
              filid,gtp,endfil,jzsol,padd,dtadd,xzmin,xzmax,ltether,l_z,pi,eta,p_tether, &
              xcen,ycen,zcen,xa,ya,za,xb,yb,zb,xmemb,ymemb,zmemb,xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf)

   implicit none

   integer, value::nmemb,nxsol,nphisol,zlen,nmono
   integer nftsa,nfil,npoly

   integer jx,jp,length,jdir,j,jstart,j1,j2,ja1,ja2,joverlap,n,jx0,jp0,jxget,jpget,jskip,nskip,npick,nm

   integer,allocatable,dimension(:)::fstart,flen,fil2a,a2mem,a2fil,filid,gtp,endfil

   integer,allocatable,dimension(:,:)::jzsol

   integer,allocatable,dimension(:)::mark

   real(kind=8),value::padd,dtadd,xzmin,xzmax,ltether,l_z,pi,eta,p_tether

   real(kind=8)::r,x0,y0,z0,xnor,ynor,znor,rad,dphi,tan0,phi0,phi,dist2,dx,dy,dz,proj,d2
   real(kind=8)::ax,ay,az,bx,by,bz,l_z_by2,dist,invdist,twoltether

   real(kind=8),allocatable,dimension(:)::xcen,ycen,zcen,xa,ya,za

   real(kind=8),allocatable,dimension(:,:)::xb,yb,zb

   real(kind=8),allocatable,dimension(:),intent(in)::xmemb,ymemb,zmemb
   real(kind=8),allocatable,dimension(:,:),intent(in)::xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf

   call random_number(r)

   if(padd*dtadd<r) return

   print*,'add a new filament'

   l_z_by2=0.5d0*l_z

   twoltether=2*ltether

   allocate(mark(nmemb))

   mark(a2mem(1:nftsa))=1

71 call random_number(r)


   jx=nxsol*r+1

   if(xsurf(1,jx)<xzmin.or.xsurf(1,jx)>xzmax) goto 71

   call random_number(r)

   jp=nphisol*r+1

   xnor=xnorsurf(jp,jx)
   ynor=ynorsurf(jp,jx)
   znor=znorsurf(jp,jx)


!  center of the new filament:

   x0=xsurf(jp,jx)+xnor*ltether
   y0=ysurf(jp,jx)+ynor*ltether
   z0=zsurf(jp,jx)+znor*ltether

   rad=sqrt(y0*y0+z0*z0)

   dphi=l_z/rad

   tan0=z0/y0

   phi0=atan(tan0)

   if(y0<0.0d0) phi0=phi0+pi


!  pick the length for filament

   call random_number(r)

   length=zlen+(r-0.5)*10

   if(length>nmono) length=nmono

!  direction of filament

   call random_number(r)

   if(r>0.5)then
      jdir=1
   else
      jdir=-1
   end if

!  coordinates of filament

   xcen(npoly+1:npoly+length)=x0

   do j=1,length/2

      phi=phi0-(j-1)*jdir*dphi

      ycen(npoly+length/2-j+1)=rad*cos(phi)

      zcen(npoly+length/2-j+1)=rad*sin(phi)

   end do

   do j=length/2+1,length

      phi=phi0+(j-length/2)*jdir*dphi

      ycen(npoly+j)=rad*cos(phi)

      zcen(npoly+j)=rad*sin(phi)

   end do

!  check if overlaping with existing filaments occurs

   joverlap=0

   do n=1,nfil

      jstart=fstart(n)

      if(abs(x0-xcen(jstart))>eta) cycle

      do j1=1,flen(n)

         ja1=jstart+j1-1

         do j2=1,length

            ja2=npoly+j2

            if(abs(ycen(ja1)-ycen(ja2))<eta)then
               joverlap=1
               exit
            end if

            if(abs(zcen(ja1)-zcen(ja2))<eta)then
               joverlap=1
               exit
            end if

         end do

         if(joverlap==1) exit

      end do

      if(joverlap==1) exit

   end do

   if(joverlap==1) goto 71

!  assign solid angle of first half

   jx0=jx
   jp0=jp

   jzsol(1,npoly+length/2)=jx0
   jzsol(2,npoly+length/2)=jp0

   do j=2,length/2

      dist2=1e10

      jxget=0

      y0=ycen(npoly+length/2-j+1)

      z0=zcen(npoly+length/2-j+1)

      do j1=1,3

         jx=jx0+j1-2

         if(jx<1)then
            cycle
         end if

         if(jx>nxsol)then
            cycle
         end if

         do j2=1,3

            jp=jp0+j2-2

            if(jp<1)then
               jp=jp+nphisol
            end if

            if(jp>nphisol)then
               jp=jp-nphisol
            end if

            dx=x0-xsurf(jp,jx)
            dy=y0-ysurf(jp,jx)
            dz=z0-zsurf(jp,jx)

            xnor=xnorsurf(jp,jx)
            ynor=ynorsurf(jp,jx)
            znor=znorsurf(jp,jx)


            proj=dx*xnor+dy*ynor+dz*znor

            d2=dx*dx+dy*dy+dz*dz-proj*proj

            if(dist2>d2)then

               dist2=d2

               jxget=jx

               jpget=jp

            end if

         end do

      end do

      if(jxget==0)then
         print*,'error in jxget',jx0,jp0
         stop
      end if

      jx0=jxget

      jp0=jpget

      jzsol(1,npoly+length/2-j+1)=jx0

      jzsol(2,npoly+length/2-j+1)=jp0

   end do

!  assign solid angle of 2nd half

   jx0=jzsol(1,npoly+length/2)
   jp0=jzsol(2,npoly+length/2)

   do j=length/2+1,length

      dist2=1e10

      jxget=0

      y0=ycen(npoly+j)

      z0=zcen(npoly+j)

      do j1=1,3

         jx=jx0+j1-2

         if(jx<1)then
            cycle
         end if

         if(jx>nxsol)then
            cycle
         end if

         do j2=1,3

            jp=jp0+j2-2

            if(jp<1)then
               jp=jp+nphisol
            end if

            if(jp>nphisol)then
               jp=jp-nphisol
            end if

            dx=x0-xsurf(jp,jx)
            dy=y0-ysurf(jp,jx)
            dz=z0-zsurf(jp,jx)

            xnor=xnorsurf(jp,jx)
            ynor=ynorsurf(jp,jx)
            znor=znorsurf(jp,jx)


            proj=dx*xnor+dy*ynor+dz*znor

            d2=dx*dx+dy*dy+dz*dz-proj*proj

            if(dist2>d2)then

               dist2=d2

               jxget=jx

               jpget=jp

            end if

         end do

      end do

      if(jxget==0)then
         print*,'error in jxget',jx0,jp0
         stop
      end if

      jx0=jxget

      jp0=jpget

      jzsol(1,npoly+j)=jx0

      jzsol(2,npoly+j)=jp0

   end do

!  assign FtsA and coordinates for 4 Z-subunits

   jskip=0

   do j=1,length

      y0=ycen(npoly+j)

      z0=zcen(npoly+j)

      jx=jzsol(1,npoly+j)

      jp=jzsol(2,npoly+j)

      xnor=xnorsurf(jp,jx)

      ynor=ynorsurf(jp,jx)

      znor=znorsurf(jp,jx)

      ax=xnor*l_z_by2
      ay=ynor*l_z_by2
      az=znor*l_z_by2


      dist=sqrt(1.0d0-xnor*xnor)

      invdist=1.0d0/dist

      bx=dist*l_z_by2

      by=-xnor*ynor*invdist*l_z_by2

      bz=-xnor*znor*invdist*l_z_by2

      xb(1,npoly+j)=x0-ax-bx
      yb(1,npoly+j)=y0-ay-by
      zb(1,npoly+j)=z0-az-bz

      xb(2,npoly+j)=x0-ax+bx
      yb(2,npoly+j)=y0-ay+by
      zb(2,npoly+j)=z0-az+bz

      xb(3,npoly+j)=x0+ax+bx
      yb(3,npoly+j)=y0+ay+by
      zb(3,npoly+j)=z0+az+bz

      xb(4,npoly+j)=x0+ax-bx
      yb(4,npoly+j)=y0+ay-by
      zb(4,npoly+j)=z0+az-bz

!     now FtsA

      call random_number(r)

      nskip=p_tether*jskip

      if(p_tether>r.or.nskip>=2)then

         npick=0


         dist2=4*ltether*ltether

         do nm=1,nmemb

            if(mark(nm)==1) cycle

            dx=x0-xmemb(nm)
            dy=y0-ymemb(nm)
            dz=z0-zmemb(nm)

            if(abs(dx)>twoltether) cycle

            if(abs(dy)>twoltether) cycle

            if(abs(dz)>twoltether) cycle


            d2=dx*dx+dy*dy+dz*dz

            if(d2<dist2)then

               dist2=d2

               npick=nm

            end if

         end do

         if(npick==0)then

            jskip=jskip+1

            print*,'could not find npick'
            cycle
         end if

         nftsa=nftsa+1

         fil2a(npoly+j)=nftsa

         a2mem(nftsa)=npick

         a2fil(nftsa)=npoly+j

         mark(npick)=1

         xa(nftsa)=0.5d0*(x0+xmemb(npick))
         ya(nftsa)=0.5d0*(y0+ymemb(npick))
         za(nftsa)=0.5d0*(z0+zmemb(npick))

         jskip=0

      else

         jskip=jskip+1

      end if

   end do

   deallocate(mark)

!  filament identities

   nfil=nfil+1

   fstart(nfil)=npoly+1

   flen(nfil)=length

   filid(npoly+1:npoly+length)=nfil

   gtp(npoly+1:npoly+length)=1

   npoly=npoly+length

   endfil(npoly)=1

   end subroutine

!=========================================================

   subroutine allpairs(npoly,nmemb,nxclude,npair5_z,npair5_mb,filid, &
              pair5_z,pair5_mb,cutoff,cos_t_2,xcen,ycen,zcen,xmemb,ymemb,zmemb)

   implicit none

   integer,value::npoly,nmemb,nxclude
   integer npair5_z,npair5_mb

   integer nnei,n1,n,nf,n2,jmb,jo,j

   integer,allocatable,dimension(:),intent(in)::filid
   integer,allocatable,dimension(:,:)::pair5_z,pair5_mb

   integer,allocatable,dimension(:,:)::pairtem,pairtemmb
   integer,allocatable,dimension(:)::pairnum,pairnummb

   real(kind=8),value::cutoff,cos_t_2
   real(kind=8)::cutoff2,d2max,dx,dy,dz,d2,dx2,rad1_2,rad2_2,cos_2

   real(kind=8),allocatable,dimension(:),intent(in)::xcen,ycen,zcen
   real(kind=8),allocatable,dimension(:),intent(in)::xmemb,ymemb,zmemb

!  for FtsZ pairs for Lennard-Jones

   cutoff2=2*cutoff

   allocate(pairtem(npoly,npoly),pairnum(npoly))

!  for membrane pairs:

   d2max=100.0d0!5*d2max

   nnei=500

   allocate(pairtemmb(nnei,nmemb),pairnummb(nmemb))

!$omp parallel &

!$omp default(none) &

!$omp private(n1,n,nf,n2,jmb,jo,dx,dy,dz,d2,dx2,rad1_2,rad2_2,cos_2) &
!$omp shared(npoly,nmemb,nxclude,nnei,filid,pairtem,pairtemmb,pairnum,pairnummb) &
!$omp shared(cutoff2,d2max,cos_t_2,xcen,ycen,zcen,xmemb,ymemb,zmemb)

!$omp do schedule(guided,32)

   do n1=1,npoly-1


      n=0

      nf=filid(n1)

      do n2=n1+1,npoly

         if(filid(n2)==nf.and.n2<n1+nxclude) cycle

         dx=xcen(n1)-xcen(n2)
         dy=ycen(n1)-ycen(n2)
         dz=zcen(n1)-zcen(n2)

         d2=abs(dx)+abs(dy)+abs(dz)!dx*dx+dy*dy+dz*dz

         if(d2<cutoff2)then

            n=n+1

            pairtem(n,n1)=n2

         end if

      end do

      pairnum(n1)=n

   end do

!$omp end do nowait



!$omp do schedule(guided,32)


   do jmb=1,nmemb-1

      n=0

      do jo=jmb+1,nmemb

         dx=0.5d0*(xmemb(jmb)-xmemb(jo))

         d2=abs(2*dx)+abs(ymemb(jmb)-ymemb(jo))+abs(zmemb(jmb)-zmemb(jo))

         if(d2>d2max)then
            cycle
         end if

         dx2=dx*dx

         rad1_2=dx2+ymemb(jmb)*ymemb(jmb)+zmemb(jmb)*zmemb(jmb)

         rad2_2=dx2+ymemb(jo)*ymemb(jo)+zmemb(jo)*zmemb(jo)

         cos_2=(dx2-ymemb(jmb)*ymemb(jo)-zmemb(jmb)*zmemb(jo))**2/rad1_2/rad2_2

         if(cos_2>cos_t_2)then

            if(n==nnei)then
               exit
            end if

            n=n+1
            pairtemmb(n,jmb)=jo


         end if

      end do

      pairnummb(jmb)=n

   end do

!$omp end do nowait
!$omp end parallel



   npair5_z=0

   do n=1,npoly-1

      if(pairnum(n)==0)then
         cycle
      end if

      do j=1,pairnum(n)

         npair5_z=npair5_z+1

         pair5_z(1,npair5_z)=n

         pair5_z(2,npair5_z)=pairtem(j,n)

      end do

   end do



   npair5_mb=0

   do n=1,nmemb-1

      if(pairnummb(n)==0)then
         cycle
      end if

      do j=1,pairnummb(n)!nnei

         npair5_mb=npair5_mb+1

         pair5_mb(1,npair5_mb)=n

         pair5_mb(2,npair5_mb)=pairtemmb(j,n)


      end do

   end do

   deallocate(pairtem,pairnum,pairtemmb,pairnummb)

   end subroutine

!==========================================

   subroutine solidupdate(nxsol,nphisol,nmemb,npoly,nftsa,jmbsol,jsursol,jzsol,a2mem, &
              pi,delta,dxsol,dphisol,xmemb,ymemb,zmemb,xcen,ycen,zcen,xwall,ywall,zwall, &
              xnorwall,ynorwall,znorwall,xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf)


   implicit none

   integer,value:: nxsol,nphisol,nmemb,npoly,nftsa
   integer n,jx,jp,jx0,jp0,j1,j2,jxget,jpget,jm

   integer,allocatable,dimension(:,:)::jmbsol,jsursol,jzsol
   integer,allocatable,dimension(:),intent(in)::a2mem

   real(kind=8),value::pi,delta,dxsol,dphisol
   real(kind=8)::xmin,dxsolby2,piby2,dphisolby2,twopi,phi,arg
   real(kind=8)::dx,dy,dz,d2,dist2

   real(kind=8),allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb
   real(kind=8),allocatable,intent(in),dimension(:)::xcen,ycen,zcen
   real(kind=8),allocatable,intent(in),dimension(:,:)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall
   real(kind=8),allocatable,intent(in),dimension(:,:)::xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf


   xmin=-(nxsol-1)/2*dxsol

   dxsolby2=0.5d0*dxsol

   piby2=0.5d0*pi

   dphisolby2=0.5d0*dphisol

   twopi=2*pi


!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,phi,arg,jp,jx0,jp0,j1,j2,jxget,jpget,dx,dy,dz,d2,dist2) &
!$omp shared(nmemb,xmemb,xmin,dxsol,dxsolby2,nxsol,jsursol,ymemb,delta,zmemb,piby2) &
!$omp shared(dphisolby2,nphisol,pi,twopi,dphisol) &
!$omp shared(xwall,ywall,zwall,xnorwall,ynorwall,znorwall,jmbsol) &
!$omp shared(xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf) &

!$omp shared(npoly,xcen,ycen,zcen,jzsol)

!$omp do schedule(guided,64)

   do n=1,nmemb

!     for membrane to interact with the ring:

      jx=1+(xmemb(n)-xmin)/dxsol

      if(xmemb(n)-xmin-(jx-1)*dxsol>dxsolby2)then
         jx=jx+1
      end if

      if(jx<1)then
         jx=1
      end if

      if(jx>nxsol)then
         jx=nxsol
      end if

      jsursol(1,n)=jx

      if(abs(ymemb(n))<delta)then
         if(zmemb(n)>0.0d0)then
            phi=piby2
         else
            phi=-piby2
         end if

      else

         arg=zmemb(n)/ymemb(n)

         phi=atan(arg)

         if(ymemb(n)<0.0d0)then
            phi=phi+pi
         end if

         if(arg<0.0d0.and.ymemb(n)>0.0d0)then
            phi=phi+twopi
         end if


      end if

      if(phi<dphisolby2)then
         jp=nphisol
      else

         jp=phi/dphisol

         if(phi-jp*dphisol>dphisolby2)then
            jp=jp+1
         end if

      end if

      jsursol(2,n)=jp

!------------------
!     for membrane to interact with the wall:

      jx0=jmbsol(1,n)

      jp0=jmbsol(2,n)

      jxget=jx0
      jpget=jp0

      dist2=10000000.0d0

      do j1=1,3

         jx=jx0-2+j1

         if(jx<1) cycle

         if(jx>nxsol) cycle

         do j2=1,3

            jp=jp0-2+j2

            if(jp<1) jp=jp+nphisol

            if(jp>nphisol) jp=jp-nphisol


            dx=xmemb(n)-xwall(jp,jx)
            dy=ymemb(n)-ywall(jp,jx)
            dz=zmemb(n)-zwall(jp,jx)

            arg=dx*xnorwall(jp,jx)+dy*ynorwall(jp,jx)+dz*znorwall(jp,jx)

            if(arg<0.0d0)cycle

            dx=dx-arg*xnorwall(jp,jx)
            dy=dy-arg*ynorwall(jp,jx)
            dz=dz-arg*znorwall(jp,jx)

            d2=dx*dx+dy*dy+dz*dz

            if(dist2>d2)then

               dist2=d2

               jxget=jx

               jpget=jp

            end if

         end do

      end do

      jmbsol(1,n)=jxget

      jmbsol(2,n)=jpget



   end do

!$omp end do nowait

!----------------------------


!$omp do schedule (guided,32)

   do n=1,npoly


      jx0=jzsol(1,n)

      jp0=jzsol(2,n)

      jxget=jx0
      jpget=jp0

      dist2=10000000.0d0


      do j1=1,3

         jx=jx0-2+j1

         if(jx<1) cycle

         if(jx>nxsol) cycle


         do j2=1,3

            jp=jp0+j2-2

            if(jp<1) jp=jp+nphisol

            if(jp>nphisol) jp=jp-nphisol


            dx=xcen(n)-xsurf(jp,jx)
            dy=ycen(n)-ysurf(jp,jx)
            dz=zcen(n)-zsurf(jp,jx)


            arg=dx*xnorsurf(jp,jx)+dy*ynorsurf(jp,jx)+dz*znorsurf(jp,jx)

            if(arg<0.0d0)cycle

            dx=dx-arg*xnorsurf(jp,jx)
            dy=dy-arg*ynorsurf(jp,jx)
            dz=dz-arg*znorsurf(jp,jx)

            d2=dx*dx+dy*dy+dz*dz

            if(dist2>d2)then

               dist2=d2

               jxget=jx

               jpget=jp

            end if

         end do

      end do

      jzsol(1,n)=jxget

      jzsol(2,n)=jpget

   end do

!$omp end do nowait



!$omp end parallel

!   do n=1,nftsa

!      jm=a2mem(n)

!      jasol(1:2,n)=jsursol(1:2,jm)

!   end do


   end subroutine


!=========================================================

   subroutine setpair(nmemb,npair5_z,npair5_mb,npair_z,npair_mb,pair5_z,pair5_mb, &
              pair_z,pair_mb,cutoff,l_pair,thet2by2,xcen,ycen,zcen,xmemb,ymemb,zmemb, &
              pairpart,boundtyp,pairtyp,l_mem,xboundmin,xboundmax,shift)

   implicit none

   integer,value::npair5_z,npair5_mb,nmemb
   integer npair_z,npair_mb

   integer n,jmb,jo,ja,jb,jc,jd,np1,np2,nnew

   integer n1,n2,n3,n4,jexit,j1,j2,j3,j4,m1,m2,jp

   integer,allocatable,intent(in),dimension(:,:)::pair5_z,pair5_mb!,jmbsol
   integer,allocatable,dimension(:,:)::pair_z,pair_mb,pairpart
   integer,allocatable,dimension(:)::boundtyp,pairtyp

   integer,allocatable,dimension(:)::mark,markmb,partlen

   integer,allocatable,dimension(:,:)::partner


   real(kind=8),value::cutoff,l_pair,thet2by2,l_mem,xboundmin,xboundmax,shift

   real(kind=8):: cutoff2,dx,dy,dz,d2,dx2,rad1,rad2,arg,dist,xc,yc,zc,xpsum,ypsum,zpsum

   real(kind=8):: xl1,yl1,zl1,xl2,yl2,zl2,xp1,yp1,zp1,xp2,yp2,zp2,xp3,yp3,zp3
   real(kind=8)::EX1,EY1,EZ1,EX2,EY2,EZ2,TX,TY,TZ,PX,PY,PZ,QX,QY,QZ,BX,BY,BZ,CX,CY,CZ
   real(kind=8)::DET,INVDET,T,U,V,delta,low,up,up1
   real(kind=8)::xmin,xmax,d2max,xmino,xmaxo

   real(kind=8),allocatable,intent(in),dimension(:)::xcen,ycen,zcen
   real(kind=8),allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb!,xboundmin,xboundmax

   real(kind=8),allocatable,dimension(:)::pairlen,lentem

   allocate(mark(npair5_z))

   mark=0

   allocate(markmb(npair5_mb))

   markmb=0

   allocate(lentem(npair5_mb))

   cutoff2=cutoff*cutoff

!$omp parallel &
!$omp default(none) &
!$omp private(n,n1,n2,jmb,jo,dx,dy,dz,d2,dx2,rad1,rad2,arg) &
!$omp shared(npair5_z,npair5_mb,pair5_z,pair5_mb,xcen,ycen,zcen) &
!$omp shared(cutoff2,mark,xmemb,ymemb,zmemb,thet2by2,markmb,lentem)

!$omp do schedule (guided,32)

  do n=1,npair5_z

      n1=pair5_z(1,n)

      n2=pair5_z(2,n)

      dx=xcen(n1)-xcen(n2)
      dy=ycen(n1)-ycen(n2)
      dz=zcen(n1)-zcen(n2)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<cutoff2)then
         mark(n)=1
      end if

   end do

!$omp end do nowait


!$omp do schedule (guided,32)


   do n=1,npair5_mb

      jmb=pair5_mb(1,n)

      jo=pair5_mb(2,n)

      dx=0.5d0*(xmemb(jmb)-xmemb(jo))

      dx2=dx*dx

      rad1=sqrt(dx2+ymemb(jmb)*ymemb(jmb)+zmemb(jmb)*zmemb(jmb))

      rad2=sqrt(dx2+ymemb(jo)*ymemb(jo)+zmemb(jo)*zmemb(jo))

      arg=1.0d0+(dx2-ymemb(jmb)*ymemb(jo)-zmemb(jmb)*zmemb(jo))/rad1/rad2

      if(arg<thet2by2)then

         markmb(n)=1

         lentem(n)=arg

      end if

   end do

!$omp end do nowait
!$omp end parallel

   npair_z=0

   do n=1,npair5_z

      if(mark(n)==1)then

         npair_z=npair_z+1

         pair_z(1:2,npair_z)=pair5_z(1:2,n)

      end if

   end do


   allocate(pairlen(npair5_mb))

   npair_mb=0


   do n=1,npair5_mb

      if(markmb(n)==1)then

         npair_mb=npair_mb+1

         pair_mb(1:2,npair_mb)=pair5_mb(1:2,n)

         pairlen(npair_mb)=lentem(n)

      end if

   end do

   deallocate(lentem,mark,markmb)

!--------------------
!  boundary problem:

   boundtyp=0

   xmin=xboundmin+l_mem

   xmax=xboundmax-l_mem

   do n=1,nmemb


      if(xmemb(n)<xmin)then

         boundtyp(n)=1


      end if

      if(xmemb(n)>xmax)then

         boundtyp(n)=2


      end if

   end do

   d2max=5*l_mem


   do n1=1,nmemb

      if(boundtyp(n1)/=1) cycle

      do n2=1,nmemb

         if(boundtyp(n2)/=2) cycle

         dx=0.5d0*(xmemb(n1)+shift-xmemb(n2))

         d2=abs(2*dx)+abs(ymemb(n1)-ymemb(n2))+abs(zmemb(n1)-zmemb(n2))

         if(d2>d2max) cycle

         dx2=dx*dx

         rad1=sqrt(dx2+ymemb(n1)*ymemb(n1)+zmemb(n1)*zmemb(n1))

         rad2=sqrt(dx2+ymemb(n2)*ymemb(n2)+zmemb(n2)*zmemb(n2))

         arg=1.0d0+(dx2-ymemb(n1)*ymemb(n2)-zmemb(n1)*zmemb(n2))/rad1/rad2

         if(arg<thet2by2)then

            npair_mb=npair_mb+1

            pair_mb(1,npair_mb)=n1

            pair_mb(2,npair_mb)=n2

            pairlen(npair_mb)=arg


         end if

      end do

   end do


!-------------------------------------

   call sorting(pairlen(:npair_mb),pair_mb(:,:npair_mb))

!  remove cross-over pairs:

!   allocate(pairtyp(npair_mb))
   pairtyp=1

   yp3=0.0d0; zp3=0.0d0


   dist=4*l_pair

   delta=0.000000000001d0

   low=-delta
   up=1.0d0+delta
   up1=1.02d0

   do np1=1,npair_mb-1

      if(pairtyp(np1)>1)then
         cycle
      end if

      ja=pair_mb(1,np1)

      jb=pair_mb(2,np1)


      xp1=xmemb(ja)
      yp1=ymemb(ja)
      zp1=zmemb(ja)


      xp2=xmemb(jb)
      yp2=ymemb(jb)
      zp2=zmemb(jb)

      if(boundtyp(ja)==1.and.boundtyp(jb)==2) xp1=xp1+shift

      xc=0.5d0*(xp1+xp2)
      yc=0.5d0*(yp1+yp2)
      zc=0.5d0*(zp1+zp2)

      xp1=xc+(xp1-xc)*up1
      yp1=yc+(yp1-yc)*up1
      zp1=zc+(zp1-zc)*up1

      xp2=xc+(xp2-xc)*up1
      yp2=yc+(yp2-yc)*up1
      zp2=zc+(zp2-zc)*up1

      xp3=0.5d0*(xp1+xp2)


!     original sum:

      xpsum=xp1+xp2
      ypsum=yp1+yp2
      zpsum=zp1+zp2

!     original coordinates

      xmino=min(xp1,xp2)
      xmaxo=max(xp1,xp2)

!     then scaled:

      xp1=xp1+xp1-xp3
      yp1=yp1+yp1
      zp1=zp1+zp1

      xp2=xp2+xp2-xp3
      yp2=yp2+yp2
      zp2=zp2+zp2

!     math

      EX1=XP2-XP1
      EY1=YP2-YP1
      EZ1=ZP2-ZP1

      EX2=XP3-XP1
      EY2=YP3-YP1
      EZ2=ZP3-ZP1

      CX=EY2*EZ1-EY1*EZ2
      CY=EZ2*EX1-EZ1*EX2
      CZ=EX2*EY1-EX1*EY2



!$omp parallel &
!$omp default(none) &
!$omp private(np2,jc,jd) &
!$omp private(dx,dy,dz,xl1,yl1,zl1,xl2,yl2,zl2,xc,yc,zc) &
!$omp shared(EX1,EY1,EZ1,EX2,EY2,EZ2,CX,CY,CZ) &
!$omp private(TX,TY,TZ,DET,INVDET,T,U,V,BX,BY,BZ)&

!$omp shared(np1,npair_mb,pairtyp,ja,jb,xp1,yp1,zp1,xmino,xmaxo) &
!$omp shared(xmemb,ymemb,zmemb,dist,pair_mb,delta,low,up,up1,l_mem,xpsum,ypsum,zpsum) &
!$omp shared(shift,boundtyp)

!$omp do schedule(guided,64)


      do np2=np1+1,npair_mb

         if(pairtyp(np2)>2)cycle


         jc=pair_mb(1,np2)

         jd=pair_mb(2,np2)

         if(jc==ja.or.jc==jb.or.jd==ja.or.jd==jb) cycle

         if(boundtyp(ja)==1.and.boundtyp(jb)==2.and.boundtyp(jc)==0.and.boundtyp(jd)==0) cycle

         if(boundtyp(ja)==0.and.boundtyp(jb)==0.and.boundtyp(jc)==1.and.boundtyp(jd)==2) cycle

         xl1=xmemb(jc)
         yl1=ymemb(jc)
         zl1=zmemb(jc)


         xl2=xmemb(jd)
         yl2=ymemb(jd)
         zl2=zmemb(jd)

         if(boundtyp(ja)==1.and.boundtyp(jb)==2)then

            if(boundtyp(jc)==1) xl1=xl1+shift

            if(boundtyp(jd)==1) xl2=xl2+shift

         else

            if(boundtyp(ja)==1.and.boundtyp(jc)==1.and.boundtyp(jd)==2) xl2=xl2-shift

            if(boundtyp(ja)==2.and.boundtyp(jc)==1.and.boundtyp(jd)==2) xl1=xl1+shift

         end if



         if(abs(ypsum-yl1-yl2)>dist) cycle


         if(abs(zpsum-zl1-zl2)>dist) cycle


         if(abs(xpsum-xl1-xl2)>dist) cycle

         if(max(xl1,xl2)<xmino.or.min(xl1,xl2)>xmaxo) cycle

         xc=0.5d0*(xl1+xl2)
         yc=0.5d0*(yl1+yl2)
         zc=0.5d0*(zl1+zl2)

         xl1=xc+(xl1-xc)*up1
         yl1=yc+(yl1-yc)*up1
         zl1=zc+(zl1-zc)*up1

         xl2=xc+(xl2-xc)*up1
         yl2=yc+(yl2-yc)*up1
         zl2=zc+(zl2-zc)*up1

!------------------------------------------


         DX=XL2-XL1
         DY=YL2-YL1
         DZ=ZL2-ZL1

!         EX1=XP2-XP1
!         EY1=YP2-YP1
!         EZ1=ZP2-ZP1

!         EX2=XP3-XP1
!         EY2=YP3-YP1
!         EZ2=ZP3-ZP1

         TX=XL1-XP1
         TY=YL1-YP1
         TZ=ZL1-ZP1

!         PX=DY*EZ2-EY2*DZ
!         PY=DZ*EX2-EZ2*DX
!         PZ=DX*EY2-EX2*DY

!         QX=TY*EZ1-EY1*TZ
!         QY=TZ*EX1-EZ1*TX
!         QZ=TX*EY1-EX1*TY


!         DET=PX*EX1+PY*EY1+PZ*EZ1

         DET=DX*CX+DY*CY+DZ*CZ


         if(abs(det)<delta) cycle

         INVDET=1.0D0/DET

!         T=(QX*EX2+QY*EY2+QZ*EZ2)*INVDET

         T=-(TX*CX+TY*CY+TZ*CZ)*INVDET


         IF(T<low.OR.T>up) cycle

         BX=DY*TZ-DZ*TY
         BY=DZ*TX-DX*TZ
         BZ=DX*TY-DY*TX


!         U=(PX*TX+PY*TY+PZ*TZ)*INVDET

         U=-(EX2*BX+EY2*BY+EZ2*BZ)*INVDET


         IF(U<low.OR.U>up) cycle

!         V=(QX*DX+QY*DY+QZ*DZ)*INVDET

         V=(EX1*BX+EY1*BY+EZ1*BZ)*INVDET


         if(u+v<low.or.u+v>up) cycle

         pairtyp(np2)=pairtyp(np2)+1


      end do

!$omp end do nowait
!$omp end parallel

   end do

   nnew=0

   do n=1,npair_mb

      if(pairtyp(n)>2)then
         cycle
      end if

      n1=pair_mb(1,n)
      n2=pair_mb(2,n)

      if(n1==n2)then
         cycle
      end if

      nnew=nnew+1

      pair_mb(1,nnew)=n1!pair_mb2(1:2,n)

      pair_mb(2,nnew)=n2

      pairtyp(nnew)=pairtyp(n)

   end do

   npair_mb=nnew

   deallocate(pairlen)

!print*,'pairs',npair_mb

!-------------------------------

!  list of tetrahedrons:

!  partners of the pairs in tetrahedrons:

   allocate(partner(20,nmemb),partlen(nmemb))

   partlen=0

   do n=1,npair_mb

      if(pairtyp(n)==2) cycle

      n1=pair_mb(1,n)

      n2=pair_mb(2,n)

      partlen(n1)=partlen(n1)+1

      partner(partlen(n1),n1)=n2

      partlen(n2)=partlen(n2)+1

      partner(partlen(n2),n2)=n1


   end do

   pairpart=0

   do n=1,npair_mb

      if(pairtyp(n)==2) cycle

      n1=pair_mb(1,n)

      n2=pair_mb(2,n)

      n3=0

      n4=0

      jexit=0

      do j1=1,partlen(n1)

         m1=partner(j1,n1)

         do j2=1,partlen(n2)

            m2=partner(j2,n2)

            if(m1==m2)then

               if(n3==0)then

                  n3=m1

               else

                  n4=m1

                  jexit=1

                  exit

               end if

            end if

         end do

         if(jexit==1)then

            exit

         end if

      end do

      if(n4>0)then

         pairpart(1,n)=n3
         pairpart(2,n)=n4




      end if

   end do


   deallocate(partner,partlen)

   end subroutine

!=========================================================

recursive   subroutine sorting(pairlen,pair)

   implicit none
   integer iq
   integer, intent(in out),dimension(:,:)::pair
   double precision,intent(in out),dimension(:)::pairlen

   if(size(pairlen)>1)then

      call partition(pairlen,pair,iq)

      call sorting(pairlen(:iq-1),pair(:,:iq-1))

      call sorting(pairlen(iq:),pair(:,iq:))

   end if

   end subroutine sorting

!-----------------------------------------------------

   subroutine partition(pairlen,pair,marker)

   integer marker,i,j,n1,n2
   integer, intent(in out),dimension(:,:)::pair
   double precision,intent(in out),dimension(:)::pairlen
   double precision x,temp

   x=pairlen(1)

   i=0

   j=size(pairlen)+1

   do

      j=j-1

      do

         if(pairlen(j)<=x) exit
         j=j-1

      end do

      i=i+1

      do

        if(pairlen(i)>=x) exit
        i=i+1

      end do

      if(i<j)then

         temp=pairlen(i)

         n1=pair(1,i)
         n2=pair(2,i)

         pairlen(i)=pairlen(j)

         pair(1:2,i)=pair(1:2,j)

         pairlen(j)=temp

         pair(1,j)=n1

         pair(2,j)=n2

      elseif(i==j)then
         marker=i+1
         return
      else
         marker=i
         return
      end if

   end do

   end subroutine partition

!==========================================



   subroutine checkdistance(nskip1,nskip2,nmemb,nphisol,nxsol,jmbsol,ndist,jdist,wthick,gap, &
                            xmemb,ymemb,zmemb,xwall,ywall,zwall,xnorwall,ynorwall,znorwall)

   implicit none

   integer,value::nskip1,nskip2,nmemb,nphisol,nxsol

   integer ndist

   integer n,jx,jp

   integer,allocatable,intent(in),dimension(:,:)::jmbsol
   integer,allocatable,dimension(:,:)::jdist
   integer,allocatable,dimension(:,:)::check
   double precision,value::wthick,gap

   double precision dx,dy,dz,arg

   double precision,allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb
   double precision,allocatable,dimension(:,:),intent(in)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall
   double precision,allocatable,dimension(:,:)::dist

   ndist=ndist+1

   if(ndist<nskip2/nskip1/2) return

   allocate(dist(nphisol,nxsol),check(nphisol,nxsol))

   dist=1000.0d0

   check=0

   do n=1,nmemb

      jx=jmbsol(1,n)

      jp=jmbsol(2,n)

      check(jp,jx)=1

      dx=xmemb(n)-xwall(jp,jx)
      dy=ymemb(n)-ywall(jp,jx)
      dz=zmemb(n)-zwall(jp,jx)

      arg=dx*xnorwall(jp,jx)+dy*ynorwall(jp,jx)+dz*znorwall(jp,jx)

      if(arg<dist(jp,jx)) dist(jp,jx)=arg

   end do

   do jx=1,nxsol-1

      do jp=1,nphisol

         if(check(jp,jx)==0) cycle

         if(dist(jp,jx)>wthick+gap) jdist(jp,jx)=jdist(jp,jx)+1

      end do

   end do


   deallocate(dist,check)

   end subroutine

!==========================================

   subroutine surfremod(nmemb,nthreads,nphisol,nxsol,nwall,jgradient,jnogradient,ndist,jsursol,jmbsol,nsurf,jdist, &
               wthick,gap,growstep,dphisol,wallwidth,pi,wallcons,xmemb,ymemb,zmemb,xsurf,ysurf,zsurf,xnorsurf, &
                ynorsurf,znorsurf,xwall,ywall,zwall,xnorwall,ynorwall,znorwall,rwall,lwall,gradient)

   implicit none

   integer,value::nmemb,nthreads,nphisol,nxsol,nwall,jgradient,jnogradient
   integer ndist
   integer n,jx,jp,jx1,jx2,jp1,jp2,tid,omp_get_thread_num,nsum

   integer,allocatable,intent(in),dimension(:,:)::jsursol,jmbsol
   integer,allocatable,dimension(:,:)::nsurf,jdist

   integer,allocatable,dimension(:,:,:)::nsurfcount!,nwallcount
   integer,allocatable,dimension(:,:)::surfmark!,wallmark

   double precision,value::wthick,gap,growstep,dphisol,wallwidth,pi
   double precision wallcons

   double precision dx1,dy1,dz1,dx2,dy2,dz2,xn,yn,zn,invdist,cir,y0,z0,sinphisolby2
   double precision dx,dy,dz,arg,d2,rad,dwallcons,rad2surf,radsurf,dlwall
!   double precision c0,c1,y1,z1,ratio

   double precision,allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb,gradient

   double precision,allocatable,dimension(:,:)::xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf
   double precision,allocatable,dimension(:,:)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall,lwall
   double precision,allocatable,dimension(:,:,:)::rsurftemp!,dist,radmemb
   double precision,allocatable,dimension(:)::rwall
   double precision,allocatable,dimension(:)::radsurfmax!,radwallmax

   allocate(rsurftemp(nthreads,nphisol,nxsol))
   allocate(nsurfcount(nthreads,nphisol,nxsol))!,nwallcount(nthreads,nphisol,nxsol))
!   allocate(dist(nthreads,nphisol,nxsol),radmemb(nthreads,nphisol,nxsol))

   nsurfcount=0

   rsurftemp=0.0d0

!   nwallcount=0

!   dist=10000.0d0

!   radmemb=0.0d0


!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,jp,tid) &
!$omp shared(nmemb,jsursol,nsurfcount,ymemb,zmemb,rsurftemp) 
!$omp do schedule(guided,64)

   do n=1,nmemb

      tid=omp_get_thread_num()+1
!tid=1

!     for membrane surface:

      jx=jsursol(1,n)

      jp=jsursol(2,n)


      nsurfcount(tid,jp,jx)=nsurfcount(tid,jp,jx)+1

!      radmemb(tid,jp,jx)=sqrt(ymemb(n)*ymemb(n)+zmemb(n)*zmemb(n))

      rsurftemp(tid,jp,jx)=rsurftemp(tid,jp,jx)+sqrt(ymemb(n)*ymemb(n)+zmemb(n)*zmemb(n))!radmemb(tid,jp,jx)

!     for cell wall:

!      jx=jmbsol(1,n)

!      jp=jmbsol(2,n)

!      nwallcount(tid,jp,jx)=nwallcount(tid,jp,jx)+1

!      dx=xmemb(n)-xwall(jp,jx)
!      dy=ymemb(n)-ywall(jp,jx)
!      dz=zmemb(n)-zwall(jp,jx)

!      arg=dx*xnorwall(jp,jx)+dy*ynorwall(jp,jx)+dz*znorwall(jp,jx)


!      if(arg<dist(tid,jp,jx))then
!         dist(tid,jp,jx)=arg
!      end if


   end do

!$omp enddo nowait
!$omp end parallel


!--------------------------------------------------

   allocate(surfmark(nphisol,nxsol))

   surfmark=0



   radsurf=0.0d0


!  remodel cell wall if the distance maintains for 10% of duration:
   ndist=ndist/10


   do jx=1,nxsol-1



      do jp=1,nphisol

!        for membrane surface:

         nsum=sum(nsurfcount(1:nthreads,jp,jx))

         if(jx==1) nsum=nsum+sum(nsurfcount(1:nthreads,jp,nxsol))

         nsurf(jp,jx)=nsum

         if(nsum>0)then

            rad=sum(rsurftemp(1:nthreads,jp,jx))/nsum

            if(jx==1) rad=rad+sum(rsurftemp(1:nthreads,jp,nxsol))/nsum

            ysurf(jp,jx)=rad*cos(jp*dphisol)
            zsurf(jp,jx)=rad*sin(jp*dphisol)

            if(rad>radsurf)then

               radsurf=rad

            end if

         else

            surfmark(jp,jx)=1

         end if

!        for cell wall:

!         nsum=sum(nwallcount(1:nthreads,jp,jx))

!         if(jx==1) nsum=nsum+sum(nwallcount(1:nthreads,jp,nxsol))

!         if(nsum>0)then

!            d2=minval(dist(1:nthreads,jp,jx))

!            if(jx==1) d2=min(d2,minval(dist(1:nthreads,jp,nxsol)))

!            radmaxmemb=maxval(radmemb(1:nthreads,jp,jx))

!            if(jx==1) radmaxmemb=max(radmaxmemb,maxval(radmemb(1:nthreads,jp,nxsol)))

!            rad=sqrt(ywall(jp,jx)*ywall(jp,jx)+zwall(jp,jx)*zwall(jp,jx))

!            if(d2>wthick+gap.and.rad>radmaxmemb)then

         if(jdist(jp,jx)>ndist)then


            dwallcons=growstep*gap

            dwallcons=dwallcons*jnogradient+dwallcons*jgradient*gradient(jx)


!            dlwall=2*pi*dwallcons/nphisol

!            y0=ywall(jp,jx)

!            z0=zwall(jp,jx)

!            rad=sqrt(y0**2+z0**2)

!            ratio=(rad-dwallcons)/rad

!            y1=y0*ratio

!            z1=z0*ratio

!            c0=0.0d0

!            c1=0.0d0

!            if(jp==1)then

!               c0=c0+sqrt((ywall(nphisol,jx)-y0)**2+(zwall(nphisol,jx)-z0)**2)

!               c1=c1+sqrt((ywall(nphisol,jx)-y1)**2+(zwall(nphisol,jx)-z1)**2)

!            else

!               c0=c0+sqrt((ywall(jp-1,jx)-y0)**2+(zwall(jp-1,jx)-z0)**2)

!               c1=c1+sqrt((ywall(jp-1,jx)-y1)**2+(zwall(jp-1,jx)-z1)**2)

!            end if

!            if(jp==nphisol)then

!               c0=c0+sqrt((ywall(1,jx)-y0)**2+(zwall(1,jx)-z0)**2)

!               c1=c1+sqrt((ywall(1,jx)-y1)**2+(zwall(1,jx)-z1)**2)

!            else

!               c0=c0+sqrt((ywall(jp+1,jx)-y0)**2+(zwall(jp+1,jx)-z0)**2)

!               c1=c1+sqrt((ywall(jp+1,jx)-y1)**2+(zwall(jp+1,jx)-z1)**2)

!            end if

!            if(c0>c1)then

!               dlwall=0.5d0*(c0-c1)

!               if(jp==1.or.jp==nphisol) lwall(1,jx)=lwall(1,jx)-dlwall

!               if(jp>1) lwall(jp,jx)=lwall(jp,jx)-dlwall

!               if(jp<nphisol) lwall(jp+1,jx)=lwall(jp+1,jx)-dlwall

!            end if



            dlwall=dwallcons*pi/nphisol


            if((jp==1.or.jp==nphisol).and.lwall(1,jx)>dlwall) lwall(1,jx)=lwall(1,jx)-dlwall

            if(jp>1.and.lwall(jp,jx)>dlwall) lwall(jp,jx)=lwall(jp,jx)-dlwall

            if(jp<nphisol.and.lwall(jp+1,jx)>dlwall) lwall(jp+1,jx)=lwall(jp+1,jx)-dlwall

               


!            end if


         end if


      end do

   end do

!  periodic boundary:

   ysurf(1:nphisol,nxsol)=ysurf(1:nphisol,1)
   zsurf(1:nphisol,nxsol)=zsurf(1:nphisol,1)

   surfmark(1:nphisol,nxsol)=surfmark(1:nphisol,1)

   lwall(1:nphisol,nxsol)=lwall(1:nphisol,1)

!  reset distance checking

   ndist=0

   jdist=0


!-------------------------------------------------

!  calculate cell wall constriction

   sinphisolby2=sin(dphisol/2)


   do jx=1,nxsol

      cir=0.0d0

      y0=ywall(nphisol,jx)

      z0=zwall(nphisol,jx)

      do jp=1,nphisol

         dy=ywall(jp,jx)-y0

         dz=zwall(jp,jx)-z0

         cir=cir+sqrt(dy*dy+dz*dz)

         y0=ywall(jp,jx)

         z0=zwall(jp,jx)

      end do


      rad=cir/nphisol/2/sinphisolby2

      wallcons=wallcons+(rwall(jx)-rad)/nxsol

      rwall(jx)=rad



   end do


!------------------------------------

!  update "empty" membrane surface


   rad2surf=radsurf*radsurf





   do jx=1,nxsol

      do jp=1,nphisol

!        for membrane surface:

         if(surfmark(jp,jx)==1)then

            if(rad2surf<ysurf(jp,jx)*ysurf(jp,jx)+zsurf(jp,jx)*zsurf(jp,jx))then

               ysurf(jp,jx)=radsurf*cos(jp*dphisol)
               zsurf(jp,jx)=radsurf*sin(jp*dphisol)

            end if

         end if


      end do

   end do


!-------------------------------------

!  update normal vectors



   do jx=1,nxsol-1

      if(jx==1)then
         jx1=nxsol-1
      else
         jx1=jx-1
      end if

      jx2=jx+1

      do jp=1,nphisol



         if(jp==1)then
            jp1=nphisol
            jp2=jp+1
         elseif(jp==nphisol)then
            jp2=1
            jp1=jp-1
         else
            jp1=jp-1
            jp2=jp+1
         end if


!        for membrane surface:

         dx1=xsurf(jp,jx2)-xsurf(jp,jx1)
         dy1=ysurf(jp,jx2)-ysurf(jp,jx1)
         dz1=zsurf(jp,jx2)-zsurf(jp,jx1)

         dx2=xsurf(jp2,jx)-xsurf(jp1,jx)
         dy2=ysurf(jp2,jx)-ysurf(jp1,jx)
         dz2=zsurf(jp2,jx)-zsurf(jp1,jx)

         if(jx==1) dx1=dx1+wallwidth


         xn=dy1*dz2-dy2*dz1
         yn=dz1*dx2-dz2*dx1
         zn=dx1*dy2-dx2*dy1

         invdist=1.0d0/sqrt(xn*xn+yn*yn+zn*zn)

         xnorsurf(jp,jx)=xn*invdist
         ynorsurf(jp,jx)=yn*invdist
         znorsurf(jp,jx)=zn*invdist

!        for cell wall:

         dx1=xwall(jp,jx2)-xwall(jp,jx1)
         dy1=ywall(jp,jx2)-ywall(jp,jx1)
         dz1=zwall(jp,jx2)-zwall(jp,jx1)

         dx2=xwall(jp2,jx)-xwall(jp1,jx)
         dy2=ywall(jp2,jx)-ywall(jp1,jx)
         dz2=zwall(jp2,jx)-zwall(jp1,jx)

         if(jx==1) dx1=dx1+wallwidth

         xn=dy1*dz2-dy2*dz1
         yn=dz1*dx2-dz2*dx1
         zn=dx1*dy2-dx2*dy1

         invdist=1.0d0/sqrt(xn*xn+yn*yn+zn*zn)

         xnorwall(jp,jx)=xn*invdist
         ynorwall(jp,jx)=yn*invdist
         znorwall(jp,jx)=zn*invdist




      end do

   end do

!  periodic boundary:

   xnorsurf(1:nphisol,nxsol)=xnorsurf(1:nphisol,1)
   ynorsurf(1:nphisol,nxsol)=ynorsurf(1:nphisol,1)
   znorsurf(1:nphisol,nxsol)=znorsurf(1:nphisol,1)

   xnorwall(1:nphisol,nxsol)=xnorwall(1:nphisol,1)
   ynorwall(1:nphisol,nxsol)=ynorwall(1:nphisol,1)
   znorwall(1:nphisol,nxsol)=znorwall(1:nphisol,1)


   deallocate(nsurfcount,rsurftemp,surfmark)


   end subroutine

!=========================================================

   subroutine constraints(nthreads,nmemb,npoly,nftsa,nxsol,nphisol,jmbsol,jsursol,jzsol,a2mem,a2fil,fil2a,endfil,nsurf, &
               kwall,wthick,lsqueez,l_mem,k_mem,l_a_z,l_mb_a,ktether,kradial,kplane,kcircon, &
                xwall,ywall,zwall,xnorwall,ynorwall,znorwall,xsurf,ysurf,zsurf, &
                 xnorsurf,ynorsurf,znorsurf,xmemb,ymemb,zmemb,xboundmin,xboundmax, &
                  xcen,ycen,zcen,xa,ya,za,fxmemb,fymemb,fzmemb, &
                   fxrep,fyrep,fzrep,fxa,fya,fza,nskip1,nskip2,ndist,jdist,gap)



   implicit none

   integer,value::nthreads,nmemb,npoly,nftsa,nxsol,nphisol,nskip1,nskip2

   integer ndist

   integer n,jx,jp,jm,ja,tid,omp_get_thread_num,jm1,jm2

   integer,allocatable,intent(in),dimension(:,:)::jmbsol,jzsol,nsurf,jsursol
   integer,allocatable,intent(in),dimension(:)::a2mem,a2fil,fil2a,endfil
   integer,allocatable,dimension(:,:)::jdist

   real(kind=8),value::kwall,wthick,lsqueez,l_mem,k_mem,l_a_z,l_mb_a,ktether,kradial,xboundmin,xboundmax,kplane,kcircon,gap

   real(kind=8)::xn,yn,zn,dist,invdist,dx,dy,dz,d2,f,dfx,dfy,dfz,invrad,proj


   real(kind=8),allocatable,intent(in),dimension(:,:)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall
   real(kind=8),allocatable,intent(in),dimension(:,:)::xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf

   real(kind=8),allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb!,xboundmin,xboundmax
   real(kind=8),allocatable,intent(in),dimension(:)::xcen,ycen,zcen,xa,ya,za

   real(kind=8),allocatable,dimension(:)::fxmemb,fymemb,fzmemb
   real(kind=8),allocatable,dimension(:)::fxrep,fyrep,fzrep,fxa,fya,fza
   real(kind=8),allocatable,dimension(:,:)::fxsurf,fysurf,fzsurf

   double precision,allocatable,dimension(:,:)::ftem,checkdist2

   real(kind=8),allocatable,dimension(:)::checkdist1

   integer,allocatable,dimension(:,:)::check


   allocate(checkdist1(nmemb))

!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,jp) &
!$omp private(dist,xn,yn,zn,dx,dy,dz,f) &
!$omp shared(nmemb,jmbsol) &
!$omp shared(kwall,wthick,lsqueez) &
!$omp shared(xboundmin,xboundmax)&
!$omp shared(xnorwall,ynorwall,znorwall) &
!$omp shared(xwall,ywall,zwall,xmemb,ymemb,zmemb) &
!$omp shared(fxmemb,fymemb,fzmemb,checkdist1) 

!$omp do schedule(guided,64)


   do n=1,nmemb


      jx=jmbsol(1,n)

      jp=jmbsol(2,n)

!     normal vector of the wall:

      xn=xnorwall(jp,jx)
      yn=ynorwall(jp,jx)
      zn=znorwall(jp,jx)

      dx=xmemb(n)-xwall(jp,jx)
      dy=ymemb(n)-ywall(jp,jx)
      dz=zmemb(n)-zwall(jp,jx)


      dist=dx*xn+dy*yn+dz*zn-wthick

      if(dist<0.0d0)then


         f=kwall*dist*dist

         fxmemb(n)=fxmemb(n)+f*xn
         fymemb(n)=fymemb(n)+f*yn
         fzmemb(n)=fzmemb(n)+f*zn

      elseif(dist<lsqueez)then

         f=-kwall*dist*dist

         fxmemb(n)=fxmemb(n)+f*xn
         fymemb(n)=fymemb(n)+f*yn
         fzmemb(n)=fzmemb(n)+f*zn

      else

         f=-kwall*lsqueez*lsqueez

         fxmemb(n)=fxmemb(n)+f*xn
         fymemb(n)=fymemb(n)+f*yn
         fzmemb(n)=fzmemb(n)+f*zn

      end if


!     blocking from the x boundaries:

      if(xmemb(n)<xboundmin)then
         fxmemb(n)=fxmemb(n)+xboundmin-xmemb(n)
      else if(xmemb(n)>xboundmax)then
         fxmemb(n)=fxmemb(n)+xboundmax-xmemb(n)
      end if

!     record distance between the membrane and cell wall:

      checkdist1(n)=dist+wthick

   end do

!$omp enddo nowait

!$omp end parallel



!  check if cell wall remodeling is signaled:

   ndist=ndist+1

   if(ndist>nskip2/nskip1/2)then

      allocate(checkdist2(nphisol,nxsol),check(nphisol,nxsol))

      checkdist2=1000.0d0

      check=0

      do n=1,nmemb

         jx=jmbsol(1,n)

         jp=jmbsol(2,n)

         check(jp,jx)=1

         if(checkdist1(n)<checkdist2(jp,jx)) checkdist2(jp,jx)=checkdist1(n)

      end do

      do jx=1,nxsol-1

         do jp=1,nphisol

            if(check(jp,jx)==0) cycle

            if(checkdist2(jp,jx)>wthick+gap) jdist(jp,jx)=jdist(jp,jx)+1

         end do

      end do

      deallocate(checkdist2,check)

   end if

   deallocate(checkdist1)


!  tethering FtsZ-FtsA-membrane:



   allocate(ftem(nphisol,nxsol))
   ftem=0.0d0


   jm1=0


   do n=1,npoly


!     first, membrane serves as a boundary for FtsZ:

      jx=jzsol(1,n)

      jp=jzsol(2,n)


      dx=xcen(n)-xsurf(jp,jx)
      dy=ycen(n)-ysurf(jp,jx)
      dz=zcen(n)-zsurf(jp,jx)

      xn=xnorsurf(jp,jx)
      yn=ynorsurf(jp,jx)
      zn=znorsurf(jp,jx)

      dist=dx*xn+dy*yn+dz*zn

      if(dist<l_mem)then

         if(dist<0.0d0) dist=0.0d0

         f=0.25d0*k_mem*(l_mem-dist)/(dist+1.0d0)

!        for each of 4 beads on the FtsZ subunit feels 1/4 of the force

         fxrep(n)=fxrep(n)+f*xn
         fyrep(n)=fyrep(n)+f*yn
         fzrep(n)=fzrep(n)+f*zn

!         tid=omp_get_thread_num()+1

         ftem(jp,jx)=ftem(jp,jx)-f*4


      end if


!     to constraints A2MEM to the circumferential alignment:

!      if(endfil(n)==1)then
!         jm1=0
!         cycle
!      end if

      ja=fil2a(n)

      if(ja==0)then
         if(endfil(n)==1) jm1=0
         cycle
      end if

      jm2=a2mem(ja)

      if(jm1>0)then

         dx=xmemb(jm1)-xmemb(jm2)

         fxmemb(jm1)=fxmemb(jm1)-kcircon*dx

         fxmemb(jm2)=fxmemb(jm2)+kcircon*dx

      end if

      jm1=jm2

      if(endfil(n)==1) jm1=0


   end do







!  tethering FtsA to membrane and FtsZ




   do ja=1,nftsa

!      n=a2fil(ja)

!      jx=jzsol(1,n)

!      jp=jzsol(2,n)


!      xn=xnorsurf(jp,jx)
!      yn=ynorsurf(jp,jx)
!      zn=znorsurf(jp,jx)

!     to membrane

      jm=a2mem(ja)

      dx=xmemb(jm)-xa(ja)
      dy=ymemb(jm)-ya(ja)
      dz=zmemb(jm)-za(ja)

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      invdist=1.0d0/dist

      f=ktether*(dist-l_mb_a)*invdist

      dfx=f*dx
      dfy=f*dy
      dfz=f*dz

!      fxa(ja)=fxa(ja)+dfx
!      fya(ja)=fya(ja)+dfy
!      fza(ja)=fza(ja)+dfz

!      fxmemb(jm)=fxmemb(jm)-dfx
!      fymemb(jm)=fymemb(jm)-dfy
!      fzmemb(jm)=fzmemb(jm)-dfz

!     constrain the connector parallel to the division plane

      dfx=dfx+kplane*dx

!      fxa(ja)=fxa(ja)+dfx

!      fxmemb(jm)=fxmemb(jm)-dfx

!     align MB-A connection rigidly along the radial direction

      invrad=1.0d0/sqrt(ymemb(jm)**2+zmemb(jm)**2)

      yn=invrad*ymemb(jm)

      zn=invrad*zmemb(jm)

      proj=dy*yn+dz*zn

      invdist=1.0d0/sqrt(dy*dy+dz*dz)

      dfy=dfy+kradial*(dy-proj*yn)*invdist
      dfz=dfz+kradial*(dz-proj*zn)*invdist

      fxmemb(jm)=fxmemb(jm)-dfx
      fymemb(jm)=fymemb(jm)-dfy
      fzmemb(jm)=fzmemb(jm)-dfz

      fxa(ja)=fxa(ja)+dfx
      fya(ja)=fya(ja)+dfy
      fza(ja)=fza(ja)+dfz


!     align MB-A connection rigidly along the normal vector of the membrane

!      proj=dx*xn+dy*yn+dz*zn

!      dx=(dx-proj*xn)*invdist
!      dy=(dy-proj*yn)*invdist
!      dz=(dz-proj*zn)*invdist



!      fxmemb(jm)=fxmemb(jm)-knormal*dx
!      fymemb(jm)=fymemb(jm)-knormal*dy
!      fzmemb(jm)=fzmemb(jm)-knormal*dz

!      fxa(ja)=fxa(ja)+knormal*dx
!      fya(ja)=fya(ja)+knormal*dy
!      fza(ja)=fza(ja)+knormal*dz


!     to FtsZ

      n=a2fil(ja)

      dx=xcen(n)-xa(ja)
      dy=ycen(n)-ya(ja)
      dz=zcen(n)-za(ja)

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      invdist=1.0d0/dist

      f=ktether*(dist-l_a_z)*invdist

      dfx=f*dx
      dfy=f*dy
      dfz=f*dz

!     constrain the connector parallel to the division plane

      dfx=dfx+kplane*dx

!     align Z-A connection rigidly along the radial direction

      proj=dy*yn+dz*zn

      invdist=1.0d0/sqrt(dy*dy+dz*dz)

      dfy=dfy+kradial*(dy-proj*yn)*invdist
      dfz=dfz+kradial*(dz-proj*zn)*invdist


      fxa(ja)=fxa(ja)+dfx
      fya(ja)=fya(ja)+dfy
      fza(ja)=fza(ja)+dfz

!     each of 4 beads on FtsZ subunit feel 1/4 of the force:

      fxrep(n)=fxrep(n)-0.25d0*dfx
      fyrep(n)=fyrep(n)-0.25d0*dfy
      fzrep(n)=fzrep(n)-0.25d0*dfz

!     align Z-A connection rigidly along the normal vector of the membrane

!      proj=dx*xn+dy*yn+dz*zn

!      dx=(dx-proj*xn)*invdist
!      dy=(dy-proj*yn)*invdist
!      dz=(dz-proj*zn)*invdist



!      dfx=knormal*dx
!      dfy=knormal*dy
!      dfz=knormal*dz

!      fxa(ja)=fxa(ja)+dfx!knormal*dx
!      fya(ja)=fya(ja)+dfy!knormal*dy
!      fza(ja)=fza(ja)+dfz!knormal*dz

!      fxrep(n)=fxrep(n)-0.25d0*dfx!knormal*dx
!      fyrep(n)=fyrep(n)-0.25d0*dfy!knormal*dy
!      fzrep(n)=fzrep(n)-0.25d0*dfz!knormal*dz


!     to make the full-length tether rigid

!      dx=xcen(n)-xmemb(jm)
!      dy=ycen(n)-ymemb(jm)
!      dz=zcen(n)-zmemb(jm)

!      dist=sqrt(dx*dx+dy*dy+dz*dz)

!      invdist=1.0d0/dist

!      f=ktether*(dist-ltether)*invdist

!      dfx=f*dx
!      dfy=f*dy
!      dfz=f*dz

!      fxmemb(jm)=fxmemb(jm)+dfx
!      fymemb(jm)=fymemb(jm)+dfy
!      fzmemb(jm)=fzmemb(jm)+dfz

!     each of 4 beads on FtsZ subunit feel 1/4 of the force:

!      fxrep(n)=fxrep(n)-0.25d0*dfx
!      fyrep(n)=fyrep(n)-0.25d0*dfy
!      fzrep(n)=fzrep(n)-0.25d0*dfz


!     constrain the tether in parallel to the division plane

!      dfx=kplane*dx

!      fxmemb(jm)=fxmemb(jm)+dfx

!      fxrep(n)=fxrep(n)-0.25d0*dfx


   end do




!  update force on membrane

   allocate(fxsurf(nphisol,nxsol),fysurf(nphisol,nxsol),fzsurf(nphisol,nxsol))



   do jx=1,nxsol-1

      do jp=1,nphisol

         f=ftem(jp,jx)/nsurf(jp,jx)

         if(jx==1) f=f+ftem(jp,nxsol)/nsurf(jp,jx)

         fxsurf(jp,jx)=f*xnorsurf(jp,jx)
         fysurf(jp,jx)=f*ynorsurf(jp,jx)
         fzsurf(jp,jx)=f*znorsurf(jp,jx)

      end do

   end do

!  periodic boundary

   fxsurf(1:nphisol,nxsol)=fxsurf(1:nphisol,1)
   fysurf(1:nphisol,nxsol)=fysurf(1:nphisol,1)
   fzsurf(1:nphisol,nxsol)=fzsurf(1:nphisol,1)



!$omp parallel &
!$omp default(none) &
!$omp private(n,jx,jp) &
!$omp shared(nmemb,jsursol,fxmemb,fymemb,fzmemb,fxsurf,fysurf,fzsurf)

!$omp do

   do n=1,nmemb

      jx=jsursol(1,n)
      jp=jsursol(2,n)

      fxmemb(n)=fxmemb(n)+fxsurf(jp,jx)
      fymemb(n)=fymemb(n)+fysurf(jp,jx)
      fzmemb(n)=fzmemb(n)+fzsurf(jp,jx)

   end do

!$omp enddo nowait
!$omp end parallel


   deallocate(fxsurf,fysurf,fzsurf,ftem)


   end subroutine

!============================================

   subroutine membrane(npair_mb,pair_mb,pairpart,pairtyp,boundtyp,kpair,l_pair,l_mem, &
              kmemb,shift,xmemb,ymemb,zmemb,fxmemb,fymemb,fzmemb)

   implicit none

   integer,value::npair_mb
   integer n,n1,n2,n3,n4

   integer,allocatable,intent(in),dimension(:,:)::pair_mb,pairpart

   integer,allocatable,intent(in),dimension(:)::boundtyp,pairtyp

   real(kind=8),value::kpair,l_pair
   real(kind=8),value::l_mem,kmemb,shift

   real(kind=8)::dx,dy,dz,d2,dist,invdist,f,dfx,dfy,dfz
   real(kind=8):: invd2,dx34,dy34,dz34,proj,xn,yn,zn

   real(kind=8),allocatable,intent(in),dimension(:)::xmemb,ymemb,zmemb
   real(kind=8),allocatable,dimension(:)::fxmemb,fymemb,fzmemb

!$omp parallel &
!$omp default(none) &
!$omp private(n,dx,dy,dz,d2,dist,invdist,f,dfx,dfy,dfz) &
!$omp shared(npair_mb,pair_mb) &
!$omp shared(kpair,l_pair,l_mem,xmemb,ymemb,zmemb,fxmemb,fymemb,fzmemb) &
!$omp private(n1,n2,n3,n4,invd2,dx34,dy34,dz34,proj,xn,yn,zn) &
!$omp shared(pairpart,boundtyp,pairtyp,kmemb,shift)

!$omp do schedule(guided,64)


   do n=1,npair_mb

      n1=pair_mb(1,n)
      n2=pair_mb(2,n)

!if(n1==4.or.n2==4)then
!print*,'check pair',n,n1,n2

!end if

      dx=xmemb(n1)-xmemb(n2)
      dy=ymemb(n1)-ymemb(n2)
      dz=zmemb(n1)-zmemb(n2)

      if(boundtyp(n1)==1.and.boundtyp(n2)==2) dx=dx+shift

      d2=dx*dx+dy*dy+dz*dz

      dist=sqrt(d2)

!     tethering:

      if(dist>l_pair)then

         f=kpair*(dist-l_pair)*(dist-l_pair)/dist

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxmemb(n1)=fxmemb(n1)-dfx
!$omp atomic
         fymemb(n1)=fymemb(n1)-dfy
!$omp atomic
         fzmemb(n1)=fzmemb(n1)-dfz

!$omp atomic
         fxmemb(n2)=fxmemb(n2)+dfx
!$omp atomic
         fymemb(n2)=fymemb(n2)+dfy
!$omp atomic
         fzmemb(n2)=fzmemb(n2)+dfz

      elseif(dist<l_mem)then

         f=-kpair*(dist-l_mem)*(dist-l_mem)/dist

         dfx=f*dx!*invdist
         dfy=f*dy!*invdist
         dfz=f*dz!*invdist

!$omp atomic
         fxmemb(n1)=fxmemb(n1)-dfx
!$omp atomic
         fymemb(n1)=fymemb(n1)-dfy
!$omp atomic
         fzmemb(n1)=fzmemb(n1)-dfz

!$omp atomic
         fxmemb(n2)=fxmemb(n2)+dfx
!$omp atomic
         fymemb(n2)=fymemb(n2)+dfy
!$omp atomic
         fzmemb(n2)=fzmemb(n2)+dfz


!if(n1==4.or.n2==4) print*,fxmemb(n1),fxmemb(n2)


      end if


!  to keep layer structure of the membrane:


      if(pairtyp(n)==2) cycle

      n3=pairpart(1,n)
      n4=pairpart(2,n)

      if(n3==0)then
         cycle
      end if

!if(n1==n3.or.n1==n4.or.n3==n4.or.n2==n3.or.n2==n4)then
!print*,'pair',n,n1,n2,n3,n4
!stop
!end if

!if(n3==4.or.n4==4)print*,'pairpart',n1,n2,n3,n4

      dx34=xmemb(n3)-xmemb(n4)
      dy34=ymemb(n3)-ymemb(n4)
      dz34=zmemb(n3)-zmemb(n4)

      if(boundtyp(n3)==1.and.boundtyp(n4)==2) dx34=dx34+shift

      if(boundtyp(n3)==2.and.boundtyp(n4)==1) dx34=dx34-shift


      xn=dy*dz34-dz*dy34
      yn=dz*dx34-dx*dz34
      zn=dx*dy34-dy*dx34

      invdist=1.0d0/sqrt(xn*xn+yn*yn+zn*zn)

!if(n3==4.or.n4==4)then
! print*,xn,yn,zn,invdist
! print*,n3,n4
!end if


      xn=xn*invdist
      yn=yn*invdist
      zn=zn*invdist

      dx=xmemb(n1)-xmemb(n3)
      dy=ymemb(n1)-ymemb(n3)
      dz=zmemb(n1)-zmemb(n3)

      if(boundtyp(n1)==1.and.boundtyp(n3)==2) dx=dx+shift

      if(boundtyp(n1)==2.and.boundtyp(n3)==1) dx=dx-shift

      proj=dx*xn+dy*yn+dz*zn

      dx=proj*xn
      dy=proj*yn
      dz=proj*zn


      dfx=kmemb*dx
      dfy=kmemb*dy
      dfz=kmemb*dz

!$omp atomic
      fxmemb(n1)=fxmemb(n1)-dfx
!$omp atomic
      fymemb(n1)=fymemb(n1)-dfy
!$omp atomic
      fzmemb(n1)=fzmemb(n1)-dfz

!$omp atomic
      fxmemb(n2)=fxmemb(n2)-dfx
!$omp atomic
      fymemb(n2)=fymemb(n2)-dfy
!$omp atomic
      fzmemb(n2)=fzmemb(n2)-dfz

!$omp atomic
      fxmemb(n3)=fxmemb(n3)+dfx
!$omp atomic
      fymemb(n3)=fymemb(n3)+dfy
!$omp atomic
      fzmemb(n3)=fzmemb(n3)+dfz

!$omp atomic
      fxmemb(n4)=fxmemb(n4)+dfx
!$omp atomic
      fymemb(n4)=fymemb(n4)+dfy
!$omp atomic
      fzmemb(n4)=fzmemb(n4)+dfz

   end do


!$omp enddo nowait

!$omp end parallel


   end subroutine

!============================================

   subroutine wallforce(nxsol,nphisol,kgly,kpep,dxsol,pturgor,delta,invdelta,lwall,xwall,ywall,zwall,fywall,fzwall)

   implicit none

   integer,value:: nxsol,nphisol

   integer jx,jp,jp1,psign(2),n,j


   real(kind=8),value::kgly,kpep,dxsol,pturgor,delta,invdelta

   real(kind=8)::y0,z0,dx,dy,dz,dist,f,fy,fz,turgor,v0,dv

   real(kind=8)::xi,yi,zi,xj,yj,zj,xk,yk,zk,xc,yc,zc

   real(kind=8),allocatable,dimension(:,:),intent(in)::xwall,ywall,zwall,lwall


   real(kind=8),allocatable,dimension(:,:)::fywall,fzwall



!  peptide crosslinks

   fywall(1:nphisol,1)=0.0d0
   fzwall(1:nphisol,1)=0.0d0

   do jx=1,nxsol-1

      do jp=1,nphisol

         dx=xwall(jp,jx)-xwall(jp,jx+1)
         dy=ywall(jp,jx)-ywall(jp,jx+1)
         dz=zwall(jp,jx)-zwall(jp,jx+1)

         dist=sqrt(dx*dx+dy*dy+dz*dz)

         f=kpep*(dist-dxsol)/dist

         fy=f*dy
         fz=f*dz

         fywall(jp,jx)=fywall(jp,jx)-fy
         fzwall(jp,jx)=fzwall(jp,jx)-fz

         fywall(jp,jx+1)=fy
         fzwall(jp,jx+1)=fz

      end do

   end do


!  periodic boundary:

!   fywall(1:nphisol,1)=fywall(1:nphisol,1)+fywall(1:nphisol,nxsol)
!   fzwall(1:nphisol,1)=fzwall(1:nphisol,1)+fzwall(1:nphisol,nxsol)




!  glycan rigidity



   do jx=1,nxsol-1

      y0=ywall(nphisol,jx)

      z0=zwall(nphisol,jx)



      do jp=1,nphisol

         dy=ywall(jp,jx)-y0

         dz=zwall(jp,jx)-z0

         dist=sqrt(dy*dy+dz*dz)


         f=kgly*(dist-lwall(jp,jx))/dist

         fy=f*dy

         fz=f*dz

         if(jp==1)then

            fywall(nphisol,jx)=fywall(nphisol,jx)+fy

            fzwall(nphisol,jx)=fzwall(nphisol,jx)+fz

         else

            fywall(jp-1,jx)=fywall(jp-1,jx)+fy

            fzwall(jp-1,jx)=fzwall(jp-1,jx)+fz

         end if


         fywall(jp,jx)=fywall(jp,jx)-fy

         fzwall(jp,jx)=fzwall(jp,jx)-fz



         y0=ywall(jp,jx)

         z0=zwall(jp,jx)



      end do



   end do


!------------------------------------

!  turgor pressure


   xc=0.0d0
   yc=0.0d0
   zc=0.0d0

   turgor=pturgor*invdelta

   do jx=1,nxsol-1

      jp1=nphisol

      xi=xwall(jp1,jx)
      yi=ywall(jp1,jx)
      zi=zwall(jp1,jx)

      xj=xwall(jp1,jx+1)
      yj=ywall(jp1,jx+1)
      zj=zwall(jp1,jx+1)

      do j=1,nphisol

         jp=jp1

         jp1=j

         xk=xj
         yk=yj
         zk=zj

         xj=xwall(jp1,jx)
         yj=ywall(jp1,jx)
         zj=zwall(jp1,jx)

!        Volume of the tetrahedron formed by (jp,jx),(jp1,jx),(jp,jx+1) and the center:

         v0=vol(xi,yi,zi,xj,yj,zj,xk,yk,zk,xc,yc,zc)

!        force on (jp,jx) along y-axis

         dv=vol(xi,yi+delta,zi,xj,yj,zj,xk,yk,zk,xc,yc,zc)-v0

         fywall(jp,jx)=fywall(jp,jx)+turgor*dv

!        force on (jp,jx) along z-axis

         dv=vol(xi,yi,zi+delta,xj,yj,zj,xk,yk,zk,xc,yc,zc)-v0

         fzwall(jp,jx)=fzwall(jp,jx)+turgor*dv

!        force on (jp1,jx) along y-axis

         dv=vol(xi,yi,zi,xj,yj+delta,zj,xk,yk,zk,xc,yc,zc)-v0

         fywall(jp1,jx)=fywall(jp1,jx)+turgor*dv

!        force on (jp1,jx) along z-axis

         dv=vol(xi,yi,zi,xj,yj,zj+delta,xk,yk,zk,xc,yc,zc)-v0

         fzwall(jp1,jx)=fzwall(jp1,jx)+turgor*dv

!        force on (jp,jx+1) along y-axis

         dv=vol(xi,yi,zi,xj,yj,zj,xk,yk+delta,zk,xc,yc,zc)-v0

         fywall(jp,jx+1)=fywall(jp,jx+1)+turgor*dv

!        force on (jp,jx+1) along z-axis

         dv=vol(xi,yi,zi,xj,yj,zj,xk,yk,zk+delta,xc,yc,zc)-v0

         fzwall(jp,jx+1)=fzwall(jp,jx+1)+turgor*dv

!        tetrahedron formed by (jp1,jx),(jp1,jx+1),(jp,jx+1) and the center:

         xi=xj
         yi=yj
         zi=zj

         xj=xwall(jp1,jx+1)
         yj=ywall(jp1,jx+1)
         zj=zwall(jp1,jx+1)

!        Volume of the tetrahedron formed by (jp1,jx),(jp1,jx+1),(jp,jx+1) and the center:

         v0=vol(xi,yi,zi,xj,yj,zj,xk,yk,zk,xc,yc,zc)

!        force on (jp1,jx) along y-axis

         dv=vol(xi,yi+delta,zi,xj,yj,zj,xk,yk,zk,xc,yc,zc)-v0

         fywall(jp1,jx)=fywall(jp1,jx)+turgor*dv

!        force on (jp1,jx) along z-axis

         dv=vol(xi,yi,zi+delta,xj,yj,zj,xk,yk,zk,xc,yc,zc)-v0

         fzwall(jp1,jx)=fzwall(jp1,jx)+turgor*dv

!        force on (jp1,jx+1) along y-axis

         dv=vol(xi,yi,zi,xj,yj+delta,zj,xk,yk,zk,xc,yc,zc)-v0

         fywall(jp1,jx+1)=fywall(jp1,jx+1)+turgor*dv

!        force on (jp1,jx+1) along z-axis

         dv=vol(xi,yi,zi,xj,yj,zj+delta,xk,yk,zk,xc,yc,zc)-v0

         fzwall(jp1,jx+1)=fzwall(jp1,jx+1)+turgor*dv

!        force on (jp,jx+1) along y-axis

         dv=vol(xi,yi,zi,xj,yj,zj,xk,yk+delta,zk,xc,yc,zc)-v0

         fywall(jp,jx+1)=fywall(jp,jx+1)+turgor*dv

!        force on (jp,jx+1) along z-axis

         dv=vol(xi,yi,zi,xj,yj,zj,xk,yk,zk+delta,xc,yc,zc)-v0

         fzwall(jp,jx+1)=fzwall(jp,jx+1)+turgor*dv

      end do

   end do

!--------------------

!  periodic boundary

   fywall(1:nphisol,1)=fywall(1:nphisol,1)+fywall(1:nphisol,nxsol)
   fzwall(1:nphisol,1)=fzwall(1:nphisol,1)+fzwall(1:nphisol,nxsol)

!--------------------

!  at the boundary

   psign(1)=-1
   psign(2)=1

   yk=0.0d0
   zk=0.0d0

   do n=1,2

      if(n==1)then
         xk=xwall(1,1)
      else
         xk=xwall(1,nxsol)
      end if

      xi=xk
      xj=xk

      jp1=nphisol

!     due to periodic boundary, the following two lines are true for both caps

      yj=ywall(nphisol,1)
      zj=zwall(nphisol,1)

      do j=1,nphisol

         jp=jp1

         yi=yj
         zi=zj

         jp1=j

         yj=ywall(jp1,1)
         zj=zwall(jp1,1)

!        Volume of the tetrahedron

         v0=vol(xi,yi,zi,xj,yj,zj,xk,yk,zk,xc,yc,zc)

!        force on jp along y-axis

         dv=vol(xi,yi+delta,zi,xj,yj,zj,xk,yk,zk,xc,yc,zc)-v0

         fywall(jp,1)=fywall(jp,1)+turgor*dv*psign(n)

!        force on jp along z-axis

         dv=vol(xi,yi,zi+delta,xj,yj,zj,xk,yk,zk,xc,yc,zc)-v0

         fzwall(jp,1)=fzwall(jp,1)+turgor*dv*psign(n)

!        force on jp1 along y-axis

         dv=vol(xi,yi,zi,xj,yj+delta,zj,xk,yk,zk,xc,yc,zc)-v0

         fywall(jp1,1)=fywall(jp1,1)+turgor*dv*psign(n)

!        force on jp1 along z-axis

         dv=vol(xi,yi,zi,xj,yj,zj+delta,xk,yk,zk,xc,yc,zc)-v0

         fzwall(jp1,1)=fzwall(jp1,1)+turgor*dv*psign(n)

      end do

   end do

   end subroutine

!============================================

      DOUBLE PRECISION function vol(xA,yA,zA,xB,yB,zB,xC,yC,zC,xD,yD,zD)
      DOUBLE PRECISION,value:: xA,yA,zA,xB,yB,zB,xC,yC,zC,xD,yD,zD
      DOUBLE PRECISION x1,y1,z1,x2,y2,z2,x3,y3,z3,x23,y23,z23,a,b,c

      x1=xa-xd
      y1=ya-yd
      z1=za-zd

      x2=xb-xd
      y2=yb-yd
      z2=zb-zd

      x3=xc-xd
      y3=yc-yd
      z3=zc-zd

      x23=y2*z3-y3*z2
      y23=z2*x3-z3*x2
      z23=x2*y3-x3*y2

      vol=abs(x1*x23+y1*y23+z1*z23)/6

      a=(y1-y2)*(z2-z3)-(z1-z2)*(y2-y3)
      b=(z1-z2)*(x2-x3)-(x1-x2)*(z2-z3)
      c=(x1-x2)*(y2-y3)-(y1-y2)*(x2-x3)

      if(a*x1+b*y1+c*z1<0.0d0)then
         vol=-vol
      end if

      end function vol

!============================================

   subroutine ftsz(npoly,npair_z,endfil,gtp,pair_z,k_z,kb,l_z,l_diag,l_diag_new,l_z12,l_z34, &
              eta,kexcl,kplane,kcir,xzmin,xzmax,kwid,xcen,ycen,zcen,xb,yb,zb,fx_z,fy_z,fz_z)

   implicit none

   integer,value::npoly,npair_z

   integer n,n1,n2,ja,jb!,omp_get_thread_num

   integer na,nz,jz,j,j1,jx,jp

   integer nbeg(2),nend(2),jrepeat

   integer,allocatable,intent(in),dimension(:)::endfil,gtp
   integer,allocatable,intent(in),dimension(:,:)::pair_z!,jzsol

   real(kind=8),value::k_z,kb,l_z,l_diag,l_diag_new,l_z12,l_z34,eta,kexcl,kplane,kcir,xzmin,xzmax,kwid

   real(kind=8)::dx,dy,dz,d2,dist,f,fx,fy,fz,l_12,l_34,l_d,dx1,dy1,dz1,dx2,dy2,dz2,dist1,dist2
   real(kind=8)::l0,l2min,invdist,xnor,ynor,znor,proj,xnor1,ynor1,znor1,invdnor


   real(kind=8),allocatable,intent(in),dimension(:)::xcen,ycen,zcen!,fxrep,fyrep,fzrep
   real(kind=8),allocatable,intent(in),dimension(:,:)::xb,yb,zb!,xnorsurf,ynorsurf,znorsurf
   real(kind=8),allocatable,dimension(:,:)::fx_z,fy_z,fz_z!,fx_a,fy_a,fz_a





   l2min=eta*eta ! at this distance Z beads start to push one another.

   nbeg(1)=1
   nbeg(2)=2

   if(mod(npoly,2)==1)then
      nend(1)=npoly
      nend(2)=npoly-1
   else
      nend(1)=npoly-1
      nend(2)=npoly
   end if

   jrepeat=0

1  jrepeat=jrepeat+1


!$omp parallel &
!$omp default(none) &
!$omp private(ja,jb,dx,dy,dz,d2,dist,invdist,f,fx,fy,fz) &
!$omp private(dx1,dy1,dz1,dist1,dx2,dy2,dz2,dist2,l0,l_12,l_34,l_d) &
!$omp shared(jrepeat,nbeg,nend,k_z,l_z,l_diag,l_diag_new,l_z12,l_z34,kplane,endfil,kcir,gtp) &
!$omp shared(xcen,ycen,zcen,xb,yb,zb,fx_z,fy_z,fz_z,kwid,xzmin,xzmax,kb)
!omp shared(xnorsurf,ynorsurf,znorsurf,jzsol) 

!$omp do schedule(guided,32)


!  linear stiffness

   do ja=nbeg(jrepeat),nend(jrepeat),2

!     confining FtsZ within RWID

      if(xcen(ja)<xzmin) fx_z(1:4,ja)=fx_z(1:4,ja)+kwid*(xzmin-xcen(ja))

      if(xcen(ja)>xzmax) fx_z(1:4,ja)=fx_z(1:4,ja)+kwid*(xzmax-xcen(ja))


!      tid=omp_get_thread_num()+1
!tid=1

      jb=ja+1


!      invrad=1.0d0/sqrt(ycen(ja)*ycen(ja)+zcen(ja)*zcen(ja))

!      ynor=ycen(ja)*invrad
!      znor=zcen(ja)*invrad



!     make 4 beads of ja rigid on a square:

      dx=xb(1,ja)-xb(2,ja)
      dy=yb(1,ja)-yb(2,ja)
      dz=zb(1,ja)-zb(2,ja)

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      f=k_z*(dist-l_z)/dist

      fx=f*dx
      fy=f*dy
      fz=f*dz


      fx_z(1,ja)=fx_z(1,ja)-fx
      fy_z(1,ja)=fy_z(1,ja)-fy
      fz_z(1,ja)=fz_z(1,ja)-fz

      fx_z(2,ja)=fx_z(2,ja)+fx
      fy_z(2,ja)=fy_z(2,ja)+fy
      fz_z(2,ja)=fz_z(2,ja)+fz

!---------------------

      dx=xb(2,ja)-xb(3,ja)
      dy=yb(2,ja)-yb(3,ja)
      dz=zb(2,ja)-zb(3,ja)

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      invdist=1.0d0/dist

      f=k_z*(dist-l_z)*invdist

      fx=f*dx
      fy=f*dy
      fz=f*dz


      fx_z(2,ja)=fx_z(2,ja)-fx
      fy_z(2,ja)=fy_z(2,ja)-fy
      fz_z(2,ja)=fz_z(2,ja)-fz

      fx_z(3,ja)=fx_z(3,ja)+fx
      fy_z(3,ja)=fy_z(3,ja)+fy
      fz_z(3,ja)=fz_z(3,ja)+fz

!     assuming FtsZ-membrane connector is rigid as a pillar then prevent 1-4 and 2-3 sides
!     from deviating from the plane formed by the normal vector and circumferential


!     normal vector (nor) of membrane:

!      jx=jzsol(1,ja)

!      jp=jzsol(2,ja)


!      xnor=xnorsurf(jp,jx)
!      ynor=ynorsurf(jp,jx)
!      znor=znorsurf(jp,jx)

!     assuming perfect circular membrane that nor_memb points to the axis
!     normal vector (nor1) of the plane form by nor_memb and the circumferential direction:

!      invdnor=1.0d0/sqrt(1.0d0-xnor*xnor)

!      xnor1=(1.0d0-xnor*xnor)*invdnor

!      ynor1=-xnor*ynor*invdnor

!      znor1=-xnor*znor*invdnor

!     projection of ftsz_23 on nor1:

!      proj=(dx*xnor1+dy*ynor1+dz*znor1)*invdist

!      dx=proj*xnor1
!      dy=proj*ynor1
!      dz=proj*znor1

!      fx_z(2,ja)=fx_z(2,ja)-knormal*dx
!      fy_z(2,ja)=fy_z(2,ja)-knormal*dy
!      fz_z(2,ja)=fz_z(2,ja)-knormal*dz

!      fx_z(3,ja)=fx_z(3,ja)+knormal*dx
!      fy_z(3,ja)=fy_z(3,ja)+knormal*dy
!      fz_z(3,ja)=fz_z(3,ja)+knormal*dz


!     assuming the FtsZ-membrane connector is rigid and constrained to the division plane
!     then the 2-3 side is also constrained to the division plane

      fx=kplane*dx

      fx_z(2,ja)=fx_z(2,ja)-fx

      fx_z(3,ja)=fx_z(3,ja)+fx


!----------------------

      dx=xb(3,ja)-xb(4,ja)
      dy=yb(3,ja)-yb(4,ja)
      dz=zb(3,ja)-zb(4,ja)

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      f=k_z*(dist-l_z)/dist

      fx=f*dx
      fy=f*dy
      fz=f*dz


      fx_z(3,ja)=fx_z(3,ja)-fx
      fy_z(3,ja)=fy_z(3,ja)-fy
      fz_z(3,ja)=fz_z(3,ja)-fz

      fx_z(4,ja)=fx_z(4,ja)+fx
      fy_z(4,ja)=fy_z(4,ja)+fy
      fz_z(4,ja)=fz_z(4,ja)+fz

!----------------------

      dx=xb(4,ja)-xb(1,ja)
      dy=yb(4,ja)-yb(1,ja)
      dz=zb(4,ja)-zb(1,ja)

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      invdist=1.0d0/dist

      f=k_z*(dist-l_z)*invdist

      fx=f*dx
      fy=f*dy
      fz=f*dz


      fx_z(4,ja)=fx_z(4,ja)-fx
      fy_z(4,ja)=fy_z(4,ja)-fy
      fz_z(4,ja)=fz_z(4,ja)-fz

      fx_z(1,ja)=fx_z(1,ja)+fx
      fy_z(1,ja)=fy_z(1,ja)+fy
      fz_z(1,ja)=fz_z(1,ja)+fz

!     assuming FtsZ-membrane connector is rigid as a pillar
!     then constrain 1-4 sides in the plane formed by the normal vector and circumferential

!      proj=(dx*xnor1+dy*ynor1+dz*znor1)*invdist

!      dx=proj*xnor1
!      dy=proj*ynor1
!      dz=proj*znor1


!      fx_z(4,ja)=fx_z(4,ja)-knormal*dx
!      fy_z(4,ja)=fy_z(4,ja)-knormal*dy
!      fz_z(4,ja)=fz_z(4,ja)-knormal*dz

!      fx_z(1,ja)=fx_z(1,ja)+knormal*dx
!      fy_z(1,ja)=fy_z(1,ja)+knormal*dy
!      fz_z(1,ja)=fz_z(1,ja)+knormal*dz


!     assuming the FtsZ-membrane connector is rigid and constrained to the division plane
!     then the 4-1 side is also constrained to the division plane

      fx=kplane*dx

      fx_z(4,ja)=fx_z(4,ja)-fx

      fx_z(1,ja)=fx_z(1,ja)+fx


!     diagonals of the square:

      dx=xb(3,ja)-xb(1,ja)
      dy=yb(3,ja)-yb(1,ja)
      dz=zb(3,ja)-zb(1,ja)

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      f=k_z*(dist-l_diag)/dist

      fx=f*dx
      fy=f*dy
      fz=f*dz


      fx_z(3,ja)=fx_z(3,ja)-fx
      fy_z(3,ja)=fy_z(3,ja)-fy
      fz_z(3,ja)=fz_z(3,ja)-fz

      fx_z(1,ja)=fx_z(1,ja)+fx
      fy_z(1,ja)=fy_z(1,ja)+fy
      fz_z(1,ja)=fz_z(1,ja)+fz

      dx=xb(4,ja)-xb(2,ja)
      dy=yb(4,ja)-yb(2,ja)
      dz=zb(4,ja)-zb(2,ja)
      
      dist=sqrt(dx*dx+dy*dy+dz*dz)

      f=k_z*(dist-l_diag)/dist

      fx=f*dx
      fy=f*dy
      fz=f*dz


      fx_z(4,ja)=fx_z(4,ja)-fx
      fy_z(4,ja)=fy_z(4,ja)-fy
      fz_z(4,ja)=fz_z(4,ja)-fz

      fx_z(2,ja)=fx_z(2,ja)+fx
      fy_z(2,ja)=fy_z(2,ja)+fy
      fz_z(2,ja)=fz_z(2,ja)+fz




      if(endfil(ja)==1) cycle


!     make mid plane (center_ja(1,4),center_ja(2,3),center_jb(1,4),center_jb(2,3)) rigid

!     side ja_14--jb_14

      dx=0.5d0*(xb(1,ja)+xb(4,ja)-xb(1,jb)-xb(4,jb))
      dy=0.5d0*(yb(1,ja)+yb(4,ja)-yb(1,jb)-yb(4,jb))
      dz=0.5d0*(zb(1,ja)+zb(4,ja)-zb(1,jb)-zb(4,jb))

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      f=k_z*(dist-l_z)/dist

      fx=f*dx
      fy=f*dy
      fz=f*dz

      fx_z(1,ja)=fx_z(1,ja)-fx
      fy_z(1,ja)=fy_z(1,ja)-fy
      fz_z(1,ja)=fz_z(1,ja)-fz

      fx_z(4,ja)=fx_z(4,ja)-fx
      fy_z(4,ja)=fy_z(4,ja)-fy
      fz_z(4,ja)=fz_z(4,ja)-fz

      fx_z(1,jb)=fx_z(1,jb)+fx
      fy_z(1,jb)=fy_z(1,jb)+fy
      fz_z(1,jb)=fz_z(1,jb)+fz

      fx_z(4,jb)=fx_z(4,jb)+fx
      fy_z(4,jb)=fy_z(4,jb)+fy
      fz_z(4,jb)=fz_z(4,jb)+fz

!     side ja_23--jb_23

      dx=0.5d0*(xb(2,ja)+xb(3,ja)-xb(2,jb)-xb(3,jb))
      dy=0.5d0*(yb(2,ja)+yb(3,ja)-yb(2,jb)-yb(3,jb))
      dz=0.5d0*(zb(2,ja)+zb(3,ja)-zb(2,jb)-zb(3,jb))

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      f=k_z*(dist-l_z)/dist

      fx=f*dx
      fy=f*dy
      fz=f*dz

      fx_z(2,ja)=fx_z(2,ja)-fx
      fy_z(2,ja)=fy_z(2,ja)-fy
      fz_z(2,ja)=fz_z(2,ja)-fz

      fx_z(3,ja)=fx_z(3,ja)-fx
      fy_z(3,ja)=fy_z(3,ja)-fy
      fz_z(3,ja)=fz_z(3,ja)-fz

      fx_z(2,jb)=fx_z(2,jb)+fx
      fy_z(2,jb)=fy_z(2,jb)+fy
      fz_z(2,jb)=fz_z(2,jb)+fz

      fx_z(3,jb)=fx_z(3,jb)+fx
      fy_z(3,jb)=fy_z(3,jb)+fy
      fz_z(3,jb)=fz_z(3,jb)+fz



!     prevent twisting: ja_1--jb_2 = ja_2--jb_1

      dx1=xb(1,ja)-xb(2,jb)
      dy1=yb(1,ja)-yb(2,jb)
      dz1=zb(1,ja)-zb(2,jb)

      dist1=sqrt(dx1*dx1+dy1*dy1+dz1*dz1)

      dx2=xb(2,ja)-xb(1,jb)
      dy2=yb(2,ja)-yb(1,jb)
      dz2=zb(2,ja)-zb(1,jb)

      dist2=sqrt(dx2*dx2+dy2*dy2+dz2*dz2)

      f=k_z*(dist1-dist2)/dist1

      fx=f*dx1
      fy=f*dy1
      fz=f*dz1

      fx_z(1,ja)=fx_z(1,ja)-fx
      fy_z(1,ja)=fy_z(1,ja)-fy
      fz_z(1,ja)=fz_z(1,ja)-fz

      fx_z(2,jb)=fx_z(2,jb)+fx
      fy_z(2,jb)=fy_z(2,jb)+fy
      fz_z(2,jb)=fz_z(2,jb)+fz

      fx=f*dx2
      fy=f*dy2
      fz=f*dz2

      fx_z(2,ja)=fx_z(2,ja)+fx
      fy_z(2,ja)=fy_z(2,ja)+fy
      fz_z(2,ja)=fz_z(2,ja)+fz

      fx_z(1,jb)=fx_z(1,jb)-fx
      fy_z(1,jb)=fy_z(1,jb)-fy
      fz_z(1,jb)=fz_z(1,jb)-fz


!     prevent twisting: ja_3--jb_4 = ja_4--jb_3

      dx1=xb(3,ja)-xb(4,jb)
      dy1=yb(3,ja)-yb(4,jb)
      dz1=zb(3,ja)-zb(4,jb)

      dist1=sqrt(dx1*dx1+dy1*dy1+dz1*dz1)

      dx2=xb(4,ja)-xb(3,jb)
      dy2=yb(4,ja)-yb(3,jb)
      dz2=zb(4,ja)-zb(3,jb)

      dist2=sqrt(dx2*dx2+dy2*dy2+dz2*dz2)

      f=k_z*(dist1-dist2)/dist1

      fx=f*dx1
      fy=f*dy1
      fz=f*dz1

      fx_z(3,ja)=fx_z(3,ja)-fx
      fy_z(3,ja)=fy_z(3,ja)-fy
      fz_z(3,ja)=fz_z(3,ja)-fz

      fx_z(4,jb)=fx_z(4,jb)+fx
      fy_z(4,jb)=fy_z(4,jb)+fy
      fz_z(4,jb)=fz_z(4,jb)+fz

      fx=f*dx2
      fy=f*dy2
      fz=f*dz2

      fx_z(4,ja)=fx_z(4,ja)+fx
      fy_z(4,ja)=fy_z(4,ja)+fy
      fz_z(4,ja)=fz_z(4,ja)+fz

      fx_z(3,jb)=fx_z(3,jb)-fx
      fy_z(3,jb)=fy_z(3,jb)-fy
      fz_z(3,jb)=fz_z(3,jb)-fz


!     contrain filament to the circumferential direction

      dx=xcen(ja)-xcen(jb)

      fx=-kcir*dx

      fx_z(1:4,ja)=fx_z(1:4,ja)+fx

      fx_z(1:4,jb)=fx_z(1:4,jb)-fx


!     bending stiffness

!     central length

!      dx=xcen(ja)-xcen(jb)
      dy=ycen(ja)-ycen(jb)
      dz=zcen(ja)-zcen(jb)

      l0=sqrt(dx*dx+dy*dy+dz*dz)

      if(gtp(ja)==1)then

         l_12=l0
         l_34=l0

         l_d=l0*l_diag/l_z

      else

         l_12=l0*l_z12
         l_34=l0*l_z34

         l_d=l0*l_diag_new
      end if


!     constraint of ja_1--jb_1

      dx=xb(1,ja)-xb(1,jb)
      dy=yb(1,ja)-yb(1,jb)
      dz=zb(1,ja)-zb(1,jb)

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      f=kb*(dist-l_12)/dist

      fx=f*dx
      fy=f*dy
      fz=f*dz

!print*,'before'
!print*,fx_z(1,ja),fy_z(1,ja),fz_z(1,ja)

      fx_z(1,ja)=fx_z(1,ja)-fx
      fy_z(1,ja)=fy_z(1,ja)-fy
      fz_z(1,ja)=fz_z(1,ja)-fz

      fx_z(1,jb)=fx_z(1,jb)+fx
      fy_z(1,jb)=fy_z(1,jb)+fy
      fz_z(1,jb)=fz_z(1,jb)+fz

!print*,'after'
!print*,fx_z(1,ja),fy_z(1,ja),fz_z(1,ja)
!stop

!     constraint of ja_2--jb_2

      dx=xb(2,ja)-xb(2,jb)
      dy=yb(2,ja)-yb(2,jb)
      dz=zb(2,ja)-zb(2,jb)

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      f=kb*(dist-l_12)/dist

      fx=f*dx
      fy=f*dy
      fz=f*dz

      fx_z(2,ja)=fx_z(2,ja)-fx
      fy_z(2,ja)=fy_z(2,ja)-fy
      fz_z(2,ja)=fz_z(2,ja)-fz

      fx_z(2,jb)=fx_z(2,jb)+fx
      fy_z(2,jb)=fy_z(2,jb)+fy
      fz_z(2,jb)=fz_z(2,jb)+fz

!     constraint of ja_3--jb_3

      dx=xb(3,ja)-xb(3,jb)
      dy=yb(3,ja)-yb(3,jb)
      dz=zb(3,ja)-zb(3,jb)

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      f=kb*(dist-l_34)/dist

      fx=f*dx
      fy=f*dy
      fz=f*dz

      fx_z(3,ja)=fx_z(3,ja)-fx
      fy_z(3,ja)=fy_z(3,ja)-fy
      fz_z(3,ja)=fz_z(3,ja)-fz

      fx_z(3,jb)=fx_z(3,jb)+fx
      fy_z(3,jb)=fy_z(3,jb)+fy
      fz_z(3,jb)=fz_z(3,jb)+fz

!     constraint of ja_4--jb_4

      dx=xb(4,ja)-xb(4,jb)
      dy=yb(4,ja)-yb(4,jb)
      dz=zb(4,ja)-zb(4,jb)

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      f=kb*(dist-l_34)/dist

      fx=f*dx
      fy=f*dy
      fz=f*dz

      fx_z(4,ja)=fx_z(4,ja)-fx
      fy_z(4,ja)=fy_z(4,ja)-fy
      fz_z(4,ja)=fz_z(4,ja)-fz

      fx_z(4,jb)=fx_z(4,jb)+fx
      fy_z(4,jb)=fy_z(4,jb)+fy
      fz_z(4,jb)=fz_z(4,jb)+fz


!     constraint of diagonal ja_1--jb_4

      dx=xb(1,ja)-xb(4,jb)
      dy=yb(1,ja)-yb(4,jb)
      dz=zb(1,ja)-zb(4,jb)

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      f=kb*(dist-l_d)/dist

      fx=f*dx
      fy=f*dy
      fz=f*dz

      fx_z(1,ja)=fx_z(1,ja)-fx
      fy_z(1,ja)=fy_z(1,ja)-fy
      fz_z(1,ja)=fz_z(1,ja)-fz

      fx_z(4,jb)=fx_z(4,jb)+fx
      fy_z(4,jb)=fy_z(4,jb)+fy
      fz_z(4,jb)=fz_z(4,jb)+fz

!     constraint of diagonal ja_4--jb_1

      dx=xb(4,ja)-xb(1,jb)
      dy=yb(4,ja)-yb(1,jb)
      dz=zb(4,ja)-zb(1,jb)

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      f=kb*(dist-l_d)/dist

      fx=f*dx
      fy=f*dy
      fz=f*dz

      fx_z(4,ja)=fx_z(4,ja)-fx
      fy_z(4,ja)=fy_z(4,ja)-fy
      fz_z(4,ja)=fz_z(4,ja)-fz

      fx_z(1,jb)=fx_z(1,jb)+fx
      fy_z(1,jb)=fy_z(1,jb)+fy
      fz_z(1,jb)=fz_z(1,jb)+fz

!     constraint of diagonal ja_2--jb_3

      dx=xb(2,ja)-xb(3,jb)
      dy=yb(2,ja)-yb(3,jb)
      dz=zb(2,ja)-zb(3,jb)

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      f=kb*(dist-l_d)/dist

      fx=f*dx
      fy=f*dy
      fz=f*dz

      fx_z(2,ja)=fx_z(2,ja)-fx
      fy_z(2,ja)=fy_z(2,ja)-fy
      fz_z(2,ja)=fz_z(2,ja)-fz

      fx_z(3,jb)=fx_z(3,jb)+fx
      fy_z(3,jb)=fy_z(3,jb)+fy
      fz_z(3,jb)=fz_z(3,jb)+fz

!     constraint of diagonal ja_3--jb_2

      dx=xb(3,ja)-xb(2,jb)
      dy=yb(3,ja)-yb(2,jb)
      dz=zb(3,ja)-zb(2,jb)

      dist=sqrt(dx*dx+dy*dy+dz*dz)

      f=kb*(dist-l_d)/dist

      fx=f*dx
      fy=f*dy
      fz=f*dz

      fx_z(3,ja)=fx_z(3,ja)-fx
      fy_z(3,ja)=fy_z(3,ja)-fy
      fz_z(3,ja)=fz_z(3,ja)-fz

      fx_z(2,jb)=fx_z(2,jb)+fx
      fy_z(2,jb)=fy_z(2,jb)+fy
      fz_z(2,jb)=fz_z(2,jb)+fz




   end do


!$omp enddo nowait
!$omp end parallel



   if(jrepeat==1) goto 1








!  Volume exclusion among Z beads:

!omp do

   do n=1,npair_z

!      tid=omp_get_thread_num()+1
!tid=1
      n1=pair_z(1,n)

      n2=pair_z(2,n)



      dx=xcen(n1)-xcen(n2)
      dy=ycen(n1)-ycen(n2)
      dz=zcen(n1)-zcen(n2)

      d2=dx*dx+dy*dy+dz*dz

      if(d2<l2min)then

         dist=sqrt(d2)

         f=-kexcl*(dist-eta)/dist

         fx=f*dx
         fy=f*dy
         fz=f*dz

         do j=1,4

            fx_z(j,n1)=fx_z(j,n1)+fx
            fy_z(j,n1)=fy_z(j,n1)+fy
            fz_z(j,n1)=fz_z(j,n1)+fz

            fx_z(j,n2)=fx_z(j,n2)-fx
            fy_z(j,n2)=fy_z(j,n2)-fy
            fz_z(j,n2)=fz_z(j,n2)-fz

         end do


      end if




   end do

!omp end do nowait





!omp end parallel



!omp parallel &
!omp default(none) &
!omp private(n,j) &
!omp shared(npoly,nthreads) &

!omp shared(fxrep,fyrep,fzrep,fxtem,fytem,fztem,fx_z,fy_z,fz_z)

!omp do schedule(guided,32)

!   do n=1,npoly

!      do j=1,4

!         fx_z(j,n)=fxrep(n)+sum(fxtem(1:nthreads,j,n))
!         fy_z(j,n)=fyrep(n)+sum(fytem(1:nthreads,j,n))
!         fz_z(j,n)=fzrep(n)+sum(fztem(1:nthreads,j,n))

!      end do



!   end do


!omp end do nowait
!omp end parallel


   end subroutine


!============================================

   subroutine newcoor(npoly,nftsa,nmemb,fmaxmb2,viscos,dt,dttread,dtadd,runtime,fxmemb,fymemb,fzmemb, &
               fxa,fya,fza,fx_z,fy_z,fz_z,xmemb,ymemb,zmemb,xa,ya,za,xcen,ycen,zcen,xb,yb,zb, &
                jforce,nxsol,nphisol,wallwidth,fywall,fzwall,xwall,ywall,zwall,xnorwall,ynorwall,znorwall)

   implicit none

   integer,value::npoly,nftsa,nmemb,nxsol,nphisol,jforce
   integer n,j,jx,jp,jx1,jx2,jp1,jp2

   integer omp_get_thread_num,tid

   real(kind=8),value::fmaxmb2,viscos,wallwidth

   real(kind=8)::dt,runtime,dttread,dtadd!,fmax

   real(kind=8)::gam,dx1,dy1,dz1,dx2,dy2,dz2,xn,yn,zn,invdist

   real(kind=8),allocatable,dimension(:),intent(in)::fxmemb,fymemb,fzmemb,fxa,fya,fza
   real(kind=8),allocatable,dimension(:,:),intent(in)::fx_z,fy_z,fz_z,fywall,fzwall,xwall

   real(kind=8),allocatable,dimension(:)::xmemb,ymemb,zmemb,xa,ya,za,xcen,ycen,zcen
   real(kind=8),allocatable,dimension(:,:)::xb,yb,zb,ywall,zwall,xnorwall,ynorwall,znorwall

   real(kind=8)::fmax2

!   real(kind=8),allocatable,dimension(:)::fmax2

!   allocate(fmax2(nthreads))



   fmax2=max(fmaxmb2,maxval(fxa(1:nftsa)*fxa(1:nftsa)+fya(1:nftsa)*fya(1:nftsa)+fza(1:nftsa)*fza(1:nftsa)))




   do n=1,npoly


      do j=1,4

         fmax2=max(fmax2,fx_z(j,n)*fx_z(j,n)+fy_z(j,n)*fy_z(j,n)+fz_z(j,n)*fz_z(j,n))

      end do


   end do


   do jx=1,nxsol-1

      do jp=1,nphisol

         fmax2=max(fmax2,fywall(jp,jx)*fywall(jp,jx)+fzwall(jp,jx)*fzwall(jp,jx))

      end do


   end do







   gam=0.01d0/sqrt(fmax2)

   dt=viscos*gam

   runtime=runtime+dt

   dttread=dttread+dt

   dtadd=dtadd+dt
!--------------------------------

!  update cell wall

   do jx=1,nxsol-1

      do jp=1,nphisol

         ywall(jp,jx)=ywall(jp,jx)+gam*fywall(jp,jx)

         zwall(jp,jx)=zwall(jp,jx)+gam*fzwall(jp,jx)

      end do

   end do

!  periodic boundary:

   ywall(1:nphisol,nxsol)=ywall(1:nphisol,1)
   zwall(1:nphisol,nxsol)=zwall(1:nphisol,1)

!  normal vector


   if(jforce==1)then

      do jx=1,nxsol-1

         if(jx==1)then
            jx1=nxsol-1
         else
            jx1=jx-1
         end if

         jx2=jx+1

         do jp=1,nphisol

            if(jp==1)then
               jp1=nphisol
               jp2=jp+1
            elseif(jp==nphisol)then
               jp2=1
               jp1=jp-1
            else
               jp1=jp-1
               jp2=jp+1
            end if

            dx1=xwall(jp,jx2)-xwall(jp,jx1)
            dy1=ywall(jp,jx2)-ywall(jp,jx1)
            dz1=zwall(jp,jx2)-zwall(jp,jx1)

            dx2=xwall(jp2,jx)-xwall(jp1,jx)
            dy2=ywall(jp2,jx)-ywall(jp1,jx)
            dz2=zwall(jp2,jx)-zwall(jp1,jx)

            if(jx==1) dx1=dx1+wallwidth

            xn=dy1*dz2-dy2*dz1
            yn=dz1*dx2-dz2*dx1
            zn=dx1*dy2-dx2*dy1

            invdist=1.0d0/sqrt(xn*xn+yn*yn+zn*zn)

            xnorwall(jp,jx)=xn*invdist
            ynorwall(jp,jx)=yn*invdist
            znorwall(jp,jx)=zn*invdist


         end do

      end do

!     periodic boundary

      xnorwall(1:nphisol,nxsol)=xnorwall(1:nphisol,1)
      ynorwall(1:nphisol,nxsol)=ynorwall(1:nphisol,1)
      znorwall(1:nphisol,nxsol)=znorwall(1:nphisol,1)


   end if

!------------------------------

!$omp parallel &
!$omp default(none) &
!$omp private(n,j) &
!$omp shared(npoly,nftsa,nmemb,gam) &
!$omp shared(xmemb,ymemb,zmemb,fxmemb,fymemb,fzmemb) &
!$omp shared(xa,ya,za,fxa,fya,fza,xb,yb,zb,xcen,ycen,zcen,fx_z,fy_z,fz_z)

!$omp do schedule(guided,32)

   do n=1,nmemb

      xmemb(n)=xmemb(n)+gam*fxmemb(n)
      ymemb(n)=ymemb(n)+gam*fymemb(n)
      zmemb(n)=zmemb(n)+gam*fzmemb(n)

   end do

!$omp end do nowait


!$omp do 

   do n=1,nftsa


      xa(n)=xa(n)+gam*fxa(n)
      ya(n)=ya(n)+gam*fya(n)
      za(n)=za(n)+gam*fza(n)


   end do

!$omp end do nowait



!$omp do schedule(guided,32)

   do n=1,npoly

      do j=1,4

         xb(j,n)=xb(j,n)+gam*fx_z(j,n)
         yb(j,n)=yb(j,n)+gam*fy_z(j,n)
         zb(j,n)=zb(j,n)+gam*fz_z(j,n)

      end do

      xcen(n)=0.25d0*(xb(1,n)+xb(2,n)+xb(3,n)+xb(4,n))
      ycen(n)=0.25d0*(yb(1,n)+yb(2,n)+yb(3,n)+yb(4,n))
      zcen(n)=0.25d0*(zb(1,n)+zb(2,n)+zb(3,n)+zb(4,n))

   end do

!$omp end do nowait





!$omp end parallel

   end subroutine

!============================================

   subroutine ringout(nstart,natom,natom_z1,natom_z2,natom_z3,natom_z4,natom_zc,natom_a,natom_z2a, &
               nwall,nxsol,nphisol,nmemb,nftsz,npoly,nfil,zlen,nftsa,nftsamax,fstart,flen,filid,fil2a,gtp,a2mem,a2fil, &
                dxsol,dphisol,xboundmin,xboundmax,rwall,lwall,xwall,ywall,zwall,xnorwall,ynorwall,znorwall, &
                 xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf,xmemb,ymemb,zmemb, &
                  xcen,ycen,zcen,xb,yb,zb,xa,ya,za)

   implicit none

   integer,value::natom,natom_z1,natom_z2,natom_z3,natom_z4,natom_zc,natom_a,natom_z2a
   integer,value::nwall,nxsol,nphisol,nmemb,nftsz,npoly,nfil,zlen,nftsa,nftsamax
   integer nstart
   integer n

   integer,allocatable,dimension(:),intent(in)::fstart,flen,filid,fil2a,gtp,a2mem,a2fil

   real(kind=8),value::dxsol,dphisol,xboundmin,xboundmax

   real(kind=8),allocatable,dimension(:,:),intent(in)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall,lwall
   real(kind=8),allocatable,dimension(:,:),intent(in)::xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf
   real(kind=8),allocatable,dimension(:),intent(in)::xmemb,ymemb,zmemb,xcen,ycen,zcen,rwall
   real(kind=8),allocatable,dimension(:,:),intent(in)::xb,yb,zb
   real(kind=8),allocatable,dimension(:),intent(in)::xa,ya,za

   character zero*1,charid1*1,charid2*2,charid3*3
   character (len=64) filecoor,fileconf


   nstart=nstart+1

   write(zero,'(i1)')0

   if(nstart<10)then
      write(charid1,'(i1)')nstart
      fileconf='conf'//zero//zero//charid1//'.inp'
      filecoor='coor'//zero//zero//charid1//'.inp'
   elseif(nstart<100)then
      write(charid2,'(i2)')nstart
      fileconf='conf'//zero//charid2//'.inp'
      filecoor='coor'//zero//charid2//'.inp'
   else
      write(charid3,'(i3)')nstart
      fileconf='conf'//charid3//'.inp'
      filecoor='coor'//charid3//'.inp'
   end if

!  coordinates

   open(1,file=filecoor,form='unformatted')

   write(1)dxsol,dphisol
   write(1)rwall(1:nxsol),lwall(1:nphisol,1:nxsol)
   write(1)xwall(1:nphisol,1:nxsol)
   write(1)ywall(1:nphisol,1:nxsol)
   write(1)zwall(1:nphisol,1:nxsol)
   write(1)xnorwall(1:nphisol,1:nxsol)
   write(1)ynorwall(1:nphisol,1:nxsol)
   write(1)znorwall(1:nphisol,1:nxsol)


   write(1)nmemb
   write(1)xmemb(1:nmemb)
   write(1)ymemb(1:nmemb)
   write(1)zmemb(1:nmemb)
   write(1)xboundmin,xboundmax
   write(1)xsurf(1:nphisol,1:nxsol)
   write(1)ysurf(1:nphisol,1:nxsol)
   write(1)zsurf(1:nphisol,1:nxsol)
   write(1)xnorsurf(1:nphisol,1:nxsol)
   write(1)ynorsurf(1:nphisol,1:nxsol)
   write(1)znorsurf(1:nphisol,1:nxsol)

   write(1)npoly

   do n=1,npoly

      write(1)xb(1:4,n),xcen(n)
      write(1)yb(1:4,n),ycen(n)
      write(1)zb(1:4,n),zcen(n)

   end do


   write(1)nftsa

   write(1)xa(1:nftsa)
   write(1)ya(1:nftsa)
   write(1)za(1:nftsa)

!   do n=1,nftsa

!      write(1)x_a_old(1:4,n),xa(n)
!      write(1)y_a_old(1:4,n),ya(n)
!      write(1)z_a_old(1:4,n),za(n)

!   end do


   close(1)

!  write configuration

   open(1,file=fileconf)


   write(1,*)'VISUAL numbers'
   write(1,*)natom,natom_z1,natom_z2,natom_z3,natom_z4,natom_zc,natom_a,natom_z2a
   write(1,*)

   write(1,*)'CELLWALL'
   write(1,*)nwall,nxsol,nphisol
   write(1,*)

   write(1,*)'MEMBRANE'
   write(1,*)nmemb
   write(1,*)


   write(1,*)'FTSZ'

   write(1,*)npoly,nfil,nftsz,zlen
   write(1,*)filid(1:npoly)
   write(1,*)fil2a(1:npoly)
   write(1,*)gtp(1:npoly)
   write(1,*)flen(1:nfil)
   write(1,*)fstart(1:nfil)



   write(1,*)'FTSA'

   write(1,*)nftsa,nftsamax
   write(1,*)a2mem(1:nftsa)
   write(1,*)a2fil(1:nftsa)



   close(1)

   end subroutine


!=========================================================


end module
