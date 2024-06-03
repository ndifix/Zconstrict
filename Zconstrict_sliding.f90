program cytokinesis

    use ieee_arithmetic

use declare
use cdk
USE IEEE_ARITHMETIC
implicit none
integer omp_get_thread_num,omp_get_num_threads,tid

   CALL system_clock(count_rate=cr)
   rate = REAL(cr)
   call system_clock(timestart)

!-------------------------------------------------

   call paras(nstep,jobid,nthreads,jbend,nturnover,jtreadmil,jdepo,noutput,nrforce,zlen,nring,intering,overlap, &
              xgap,pi,delta,INVDELTA,beta,k_z,l_z,kz_thet,thet0_z,rateremod, &
              k_a,l_mb_a,l_a_z,ktether,ltether,mbrad,rrad,p_tether,l_mem,k_mem,kmemb,kgly,kpep,fturgor,pturgor, &
              wthick,lsqueez,kwall,walrad,rho,dxsol,dyzsol,ptread,phyd,pdepo,tadd,kpair,gap,eps,rcutoff, &
              viscos,printtime,k_scale,k_scale_mb,mwid,rwid,kwid,jgradient,jnogradient,fgrad,growstep)

!  Generating FtsZ within limit:

   xzmax=rwid/2
   xzmin=-xzmax


   call omp_set_num_threads(nthreads)

!  Generate random seed:
   CALL INIT_RANDOM_SEED()


!-------------------------------------------------


   write(*,*)'================================================'
   write(*,*)'FtsZ constriction: the sliding model'

!-----------------------------------------------------

!======================================================
!  Direction of the job:

   if(jobid==1)then
      print*,'to generate the fellowship of the ring ...'
   elseif(jobid==2)then
      goto 29
   else
      print*,'suitable jobid was not found'
      stop
   end if

!--------------------------------------------------

!--------------------------------------------------

!  This section is to create the ring.


!  First create a cell wall

!  number of cell wall bead seprated at 10 nm on one hoop of radius =1.0*walrad, also number of solid angle:
   nphisol=2*pi*walrad/dyzsol

!  the angle between these beads is also the small solid angle component:
   dphisol=2*pi/nphisol

!  number of hoops separated xdel= 10 nm within membrane width, also component of solid angle:

!OK   dxsol=10.0d0

!  number of cell wall hoops:
   nxsol=(mwid+1.0d0)/dxsol

   if(mod(nxsol,2)==0)then
      nxsol=nxsol+1
   end if

   allocate(xwall(nphisol,nxsol),ywall(nphisol,nxsol),zwall(nphisol,nxsol))

   allocate(xnorwall(nphisol,nxsol),ynorwall(nphisol,nxsol),znorwall(nphisol,nxsol))

!   rwall=walrad

!  the total number of wall beads:
   nwall=nxsol*nphisol

!  number of bonds on the hoops is the same as nwall
   nbondwal=nwall

!  adding bonds connecting hoops:
   nbondwal=nbondwal+(nxsol-1)*nphisol

   allocate(bondwal(2,nbondwal))

   nw=0
   nb=0
   xmin=-(nxsol-1)/2*dxsol
   xmax=-xmin
   xboundmin=xmin
   xboundmax=xmax

!   wallwidth=xmax-xmin

!print*,'check'

   do jx=1,nxsol
      x0=xmin+(jx-1)*dxsol
      xwall(1:nphisol,jx)=x0
      do jp=1,nphisol
         nw=nw+1
         ywall(jp,jx)=walrad*cos(jp*dphisol)
         zwall(jp,jx)=walrad*sin(jp*dphisol)
         nb=nb+1
         bondwal(1,nb)=nw

         if(jp<nphisol)then
            bondwal(2,nb)=nw+1
         else
            bondwal(2,nb)=nw+1-nphisol
         end if

         ynorwall(jp,jx)=-ywall(jp,jx)/walrad
         znorwall(jp,jx)=-zwall(jp,jx)/walrad

      end do

!     adding bonds between hoops:

      if(jx==1)then
         cycle
      end if

      do jp=nw-nphisol+1,nw
         nb=nb+1
         bondwal(1,nb)=jp-nphisol
         bondwal(2,nb)=jp
      end do

   end do

   allocate(rwall(nxsol),lwall(nphisol,nxsol))

   rwall=walrad
   l0=sqrt((ywall(1,1)-ywall(2,1))**2+(zwall(1,1)-zwall(2,1))**2)

!  to make glycan force balancing turgor at the beginning

   lwall=l0-fturgor*walrad/kgly/l0
   lpep=dxsol

!print*,'check'

!------------------------------------------------------

!------------------------------------------------------

!  Create a membrane composed of a single-layer

!  if the bead size is l_mem then number of bead per hoop:

   nmeho=2*pi*mbrad/l_mem
   dphi=2*pi/nmeho

!  number of hoops:

   nhoop=(xmax-xmin+1.0d0)/l_mem

   if(mod(nhoop,2)==0)then
      nhoop=nhoop+1
   end if

!  total number of beads:

   nmemb=nhoop*nmeho

   allocate(xmemb(nmemb),ymemb(nmemb),zmemb(nmemb))

   nm=0
   do jh=1,nhoop
      x0=-(nhoop-1)/2*l_mem+(jh-1)*l_mem !xmin+(jh-1)*l_mem
      xmemb(nm+1:nm+nmeho)=x0
      do jm=1,nmeho
         nm=nm+1
         ymemb(nm)=mbrad*cos(dphi*jm)
         zmemb(nm)=mbrad*sin(dphi*jm)
      end do
   end do

!  coarse representation of membrane surface:

   allocate(xsurf(nphisol,nxsol),ysurf(nphisol,nxsol),zsurf(nphisol,nxsol))

   xsurf=xwall

!  but:

   ysurf=ywall*mbrad/walrad
   zsurf=zwall*mbrad/walrad

!------------------------------------------------------


!  Generating FtsZ within limit:

!   xzmax=rwid/2
!   xzmin=-xzmax

!  the "origin" shifted:

   xorigin=-(nring-1)*xgap/2!-xgap/2

   if((nring-1)*xgap>rwid)then
      print*,'rings out of bound!'
      print*,'NRING =',nring
      print*,'gap between rings XGAP =',xgap
      print*,'calculated width = (nring-1)*xgap =',(nring-1)*xgap
      print*,'ring width RWID =',rwid
      stop
   end if

!  number of subuints per a complete ring:
   nsunit=2*pi*rrad/l_z

!  over-estimate total number of subunits:
   nftsz=nring*(nsunit+2*overlap*(nsunit/zlen))

!  each subunit is represented with 4 beads:
!   allocate(xb(4,nftsz),yb(4,nftsz),zb(4,nftsz))

!  center of the subunits:
   allocate(xcen(nftsz),ycen(nftsz),zcen(nftsz))

!  over-estimated # of filaments:
   nfil=4*nftsz/zlen

!  filament id:
   allocate(filid(nftsz))

!  ring id
   allocate(ringid(nftsz))

!  filament configuration:
   allocate(fstart(nfil),flen(nfil))

!  membrane tethering via FtsA:
   allocate(fil2a(nftsz))
   fil2a=0

!  limit # of FtsA:
   nftsamax=nftsz!*p_tether*2

   allocate(xa(nftsamax),ya(nftsamax),za(nftsamax))
   allocate(a2mem(nftsamax),a2fil(nftsamax))
   a2mem=0
   a2fil=0

   allocate(mark(nmemb))
   mark=0

!  filling space with subunits:

!   allocate(fill(nsunit,nring+1))
!   fill=0
   nftsa=0
!   nmono=nftsz
   npoly=0
   dphi=2*pi/nsunit
!   rad12=rrad+l_z/2
!   rad34=rad12-l_z
!   fstart0=1
!   length0=0

   flen=0
   nfil=0

!-----------------------

   do nr=1,nring
!      x0=xorigin+(nr-1)*xgap
!     pick the angle of first bead:
      call random_number(r)
      phi0=2*pi*r
      ncover=0
      do n=1,1000
         if(ncover==nsunit) exit
         nfil=nfil+1
         call random_number(r)
         length=zlen+(r-0.5)*10
         if(nsunit-ncover-length<zlen/2) length=nsunit-ncover
         filid(npoly+1:npoly+length+overlap)=nfil
         ringid(npoly+1:npoly+length+overlap)=nr
         flen(nfil)=length+overlap
         fstart(nfil)=npoly+1
         dxshift=rho/(length+overlap)
         x0=xorigin+(nr-1)*xgap
         call random_number(r)
         if(r>0.5) then
            jstart=ncover+1
            jdir=1
         else
            dxshift=-dxshift
            jstart=ncover+length+overlap
            x0=x0+rho
            jdir=-1
         end if

!        now assign coordinates:

         do jf=1,length+overlap
            npoly=npoly+1
            xcen(npoly)=x0
            phi=phi0+(jstart+(jf-1)*jdir)*dphi
            ycen(npoly)=rrad*cos(phi)
            zcen(npoly)=rrad*sin(phi)

!           assign tethering on the beads:
            if(mod(jf,2)==1)then
               jh=(x0-xmin)/l_mem+1
               jm0=(jh-1)*nmeho
               dist2=9*ltether**2
               npick=0

               do jm=jm0+1,jm0+2*nmeho
                  if(mark(jm)==1)then
                     cycle
                  end if

                  dx=xcen(npoly)-xmemb(jm)
                  dy=ycen(npoly)-ymemb(jm)
                  dz=zcen(npoly)-zmemb(jm)
                  d2=dx*dx+dy*dy+dz*dz
                  if(d2<dist2)then
                     dist2=d2
                     npick=jm
                  end if
               end do

               if(npick==0)then
                  print*,'could not find npick'
                  cycle
               end if

               nftsa=nftsa+1
               fil2a(npoly)=nftsa
               a2mem(nftsa)=npick
               a2fil(nftsa)=npoly
               mark(npick)=1
               xa(nftsa)=0.5d0*(xcen(npoly)+xmemb(npick))
               ya(nftsa)=0.5d0*(ycen(npoly)+ymemb(npick))
               za(nftsa)=0.5d0*(zcen(npoly)+zmemb(npick))
            end if

            x0=x0+dxshift
         end do

         ncover=ncover+length
      end do
   end do

   nftsz=npoly
   nftsamax=nftsa
!   nfilmax=nfil

!------------------------------
!  write out the ring:

   call makering(nmemb,nfil,nftsz,npoly,zlen,nwall,nbondwal,nxsol,nphisol,nftsa,nftsamax, &
               fil2a,fstart,flen,filid,ringid,a2mem,a2fil,bondwal,dxsol,dphisol,xboundmin,xboundmax, &
                rwall,lwall,xwall,ywall,zwall,xnorwall,ynorwall,znorwall,xsurf,ysurf,zsurf, &
                 xmemb,ymemb,zmemb,xcen,ycen,zcen,xa,ya,za)

!  clean up

   deallocate(bondwal,xcen,ycen,zcen,filid,ringid,fstart,flen,fil2a,mark)
   deallocate(rwall,lwall,xwall,ywall,zwall,xnorwall,ynorwall,znorwall,xmemb,ymemb,zmemb)
   deallocate(xsurf,ysurf,zsurf,xa,ya,za,a2mem,a2fil)
!stop

!=================================================


!  getting the initial ring:

29 call getinfo(nstart,jfile,time0,jstart,runtime,natom,natom_zc, &
         natom_a,natom_z2a,nwall,nxsol,nphisol,nmemb,nftsz,npoly,nfil,zlen,nftsa,nftsamax)

   if(jstart>nstep) stop

   print*,'=================================='
   print*,'ring dynamics'

!  setting parameters:
   ptread=ptread/1000000 ! unit = inverse micro second
   phyd=phyd/1000000 ! unit = inverse micro second
   pdepo=pdepo/1000000 ! unit = inverse micro second

!  addition rate of FtsZ filament
   padd=1.0d0/tadd/1000000 ! unit = inverse micro second

!  limit of wall remodeling rate
   rateremod=rateremod/1000000 ! unit = nm/micro second

!  cut-off distance:
   cutoff=rcutoff*rho
   rswitch=cutoff-0.5d0*rho
   fswitch=6*eps*(2*(rho/rswitch)**12-(rho/rswitch)**6)/rswitch/(cutoff-rswitch)**2

!  number of adjacent beads to be exclude from interaction:
   nxclude=cutoff/l_z+1

!  parameters for attraction between membrane beads:
   l_pair=2.0d0*l_mem
   thetpair=4*l_mem/rrad
   cos_t_2=(1.0d0-(2*thetpair)**2/2)**2
   thet2by2=thetpair*thetpair/2
   oldtime=runtime
   dtime=0.0d0
!-----------------------------------

!  system configuration:
   allocate(rwall(nxsol),lwall(nphisol,nxsol),xwall(nphisol,nxsol),ywall(nphisol,nxsol),zwall(nphisol,nxsol))
   allocate(xnorwall(nphisol,nxsol),ynorwall(nphisol,nxsol),znorwall(nphisol,nxsol))
   allocate(xmemb(nmemb),ymemb(nmemb),zmemb(nmemb),jmbsol(2,nmemb))
   allocate(xsurf(nphisol,nxsol),ysurf(nphisol,nxsol),zsurf(nphisol,nxsol))
   allocate(jsursol(2,nmemb),nsurf(nphisol,nxsol))
   allocate(xnorsurf(nphisol,nxsol),ynorsurf(nphisol,nxsol),znorsurf(nphisol,nxsol))

!  check distance between membrane and cell wall
   allocate(jdist(nphisol,nxsol))

!  center of the subunits:
   allocate(xcen(nftsz),ycen(nftsz),zcen(nftsz))
   allocate(fstart(nfil),flen(nfil),jzsol(2,nftsz),filid(nftsz),ringid(nftsz),gtp(nftsz),endfil(nftsz))

!  FtsA:
   allocate(xa(nftsamax),ya(nftsamax),za(nftsamax))

!  membrane tethering:
   allocate(fil2a(nftsz),a2mem(nftsamax),a2fil(nftsamax))

   call ringin(nstart,nmemb,nfil,npoly,nftsa,nphisol,nxsol, &
               filid,ringid,fil2a,flen,fstart,a2mem,a2fil,endfil,gtp, &
                dxsol,dphisol,xboundmin,xboundmax,rwall,lwall,xwall,ywall,zwall,xnorwall,ynorwall,znorwall, &
                 xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf,xmemb,ymemb,zmemb, &
                  xcen,ycen,zcen,xa,ya,za)

!------------------------------------------------------------

!  forces:
   allocate(fywall(nphisol,nxsol),fzwall(nphisol,nxsol))
   allocate(fxmemb(nmemb),fymemb(nmemb),fzmemb(nmemb))
!   allocate(fx_z(4,nftsz),fy_z(4,nftsz),fz_z(4,nftsz))
   allocate(fxfil(nftsz),fyfil(nftsz),fzfil(nftsz))
   allocate(fxrep(nftsz),fyrep(nftsz),fzrep(nftsz))
!   allocate(fxtem(nthreads,4,nftsz),fytem(nthreads,4,nftsz),fztem(nthreads,4,nftsz))
   allocate(fxa(nftsamax),fya(nftsamax),fza(nftsamax))

!  to add random forces:
   allocate(rxmemb(nmemb,nrforce),rymemb(nmemb,nrforce),rzmemb(nmemb,nrforce))
   allocate(rxbead(nftsz,nrforce),rybead(nftsz,nrforce),rzbead(nftsz,nrforce))
   allocate(rx_a(nftsamax,nrforce),ry_a(nftsamax,nrforce),rz_a(nftsamax,nrforce))
   jrforce=nrforce
!------------------------------

!  preserving membrane stiffness and integrity

!  pairs of membrane beads:
   allocate(pair5_mb(2,nmemb*500),pair_mb(2,nmemb*100))
   allocate(pairpart(2,nmemb*100),pairtyp(nmemb*100))

!  boundary problem:
   allocate(boundtyp(nmemb))

!  for connecting two boundaries:
   wallwidth=xboundmax-xboundmin
   shift=wallwidth+l_mem

!  pairs of FtsZ beads:
   allocate(pair5_z(2,nftsz*500),pair_z(2,nftsz*100))
!------------------------------------------------------------

!  To write out trajectory
   call random_number(r)
   junit=80*r+11
   junit1=junit+1
   junit2=junit+2
   jfile=jfile+1
   call dcdheader(junit,jfile,natom)
   nframe=0

   call writedcd(junit,nframe,natom,natom_zc,natom_a,natom_z2a,nxsol,nphisol,nmemb,nftsz,nfil,nftsa, &
                flen,a2mem,a2fil,fstart,xcen,ycen,zcen,xa,ya,za, &
                 xmemb,ymemb,zmemb,xwall,ywall,zwall)

!  initial setup solid angle indices:
   call solidset(nxsol,nphisol,nmemb,npoly,jsursol,jmbsol,jzsol,pi,delta, &
              dxsol,dphisol,xmemb,ymemb,zmemb,xcen,ycen,zcen,xwall,ywall,zwall, &
              xnorwall,ynorwall,znorwall,xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf)

!  setup wall growth gradient:
   allocate(gradient(nxsol))

   jxwid=(nxsol-1)/2
   do jx=1,nxsol
      dx=1.0d0*(-jxwid+jx-1)
      dx=dx/jxwid
      gradient(jx)=exp(-dx*dx*fgrad)
   end do

!  measuring constriction rate
   timecheck0=runtime*0.000001

12 format(f10.7,10x,f10.7,3x,f10.7)

   open(junit2,file='constrict_rate.txt')

   if(runtime<1.0)then
      write(junit2,'(a43)')'   time(s)  wall_constrict(nm)  rate (nm/s)'
      wallcons=0.0d0
   else
      length=1
      read(junit2,'(a43)')a43
      do n=1,10000000
         read(junit2,*)timecheck,rwallcons,cons_rate
         wallcons=rwallcons
         if(timecheck>=timecheck0)then
            exit
         end if
         length=length+1
      end do

      close(junit2)

      open(junit2,file='constrict_rate.txt')
      do n=1,length
         read(junit2,*)
      end do
   end if

   oldwallcons=wallcons
   ndist=1000 ! just to prevent cell wall remodeling from occurring too early
   jdist=0
   nskip1=10
   nskip2=10000
!=============================================================

!  Dynamics of the ring:
   print*,'start dynamics'
   write(*,*)'timestep     runtime     '

!  start time step:
!   dt=0.0d0
   dtmodz=0.0d0
   dtremod=0.0d0
   do jstep=jstart,nstep
      if(mod(jstep-1,nskip1)==0.or.jstep==jstart)then
         jforce=1
      else
         jforce=0
      end if

!     set random force
      if(jrforce==nrforce) &
         call rforceset(jrforce,nrforce,nmemb,nftsz,nftsamax,pi,rxmemb,rymemb,rzmemb, &
              rxbead,rybead,rzbead,rx_a,ry_a,rz_a,k_scale,k_scale_mb)

      if(mod(jstep-1,nskip2)==0.or.jstep==jstart)then
         if(mod(jstep-1,nturnover)==0.or.jstep==jstart)then
!           treadmilling of FtsZ
            if(jtreadmil==1) then
               call treadmill(nfil,nftsa,nmemb,fstart,flen,fil2a,a2fil,a2mem,gtp,jzsol, &
                     ptread,phyd,dtmodz,xcen,ycen,zcen,xa,ya,za,xmemb,ymemb,zmemb)

               dtmodz=0.0d0
            elseif(jdepo==1)then
               call depoly(nfil,npoly,nftsa,flen,fstart,filid,fil2a,a2mem,a2fil,endfil,ringid, &
                    jzsol,pdepo,dtmodz,xcen,ycen,zcen,xa,ya,za)

               dtmodz=0.0d0
            end if

!              setup pairs for long run:
            call allpairs(npoly,nmemb,nxclude,intering,npair5_z,npair5_mb,filid,ringid, &
               pair5_z,pair5_mb,cutoff,cos_t_2,xcen,ycen,zcen,xmemb,ymemb,zmemb)
         end if

!           categorize beads into solid angles:
         call solidupdate(nxsol,nphisol,nmemb,npoly,nftsa,jmbsol,jsursol,jzsol,a2mem, &
              pi,delta,dxsol,dphisol,xmemb,ymemb,zmemb,xcen,ycen,zcen,xwall,ywall,zwall, &
              xnorwall,ynorwall,znorwall,xsurf,ysurf,zsurf,xnorsurf,ynorsurf,znorsurf)

!           pair setup for exclusion effect:
         call setpair(nmemb,npair5_z,npair5_mb,npair_z,npair_mb,pair5_z,pair5_mb, &
              pair_z,pair_mb,cutoff,l_pair,thet2by2,xcen,ycen,zcen,xmemb,ymemb,zmemb, &
              pairpart,boundtyp,pairtyp,l_mem,xboundmin,xboundmax,shift)

         call surfremod(nmemb,nthreads,nphisol,nxsol,nwall,jgradient,jnogradient,ndist,jsursol,jmbsol,nsurf,jdist, &
               wthick,gap,growstep,dphisol,wallwidth,pi,wallcons,xmemb,ymemb,zmemb,xsurf,ysurf,zsurf,xnorsurf, &
                ynorsurf,znorsurf,xwall,ywall,zwall,xnorwall,ynorwall,znorwall,rwall,lwall,gradient,dtremod,rateremod)
      end if

!----------------------------------------------------------------

!     calculate forces:
      if(jforce==1)then
         jrforce=jrforce+1
         fxmemb(1:nmemb)=rxmemb(1:nmemb,jrforce)
         fymemb(1:nmemb)=rymemb(1:nmemb,jrforce)
         fzmemb(1:nmemb)=rzmemb(1:nmemb,jrforce)

         fxa(1:nftsa)=rx_a(1:nftsa,jrforce)
         fya(1:nftsa)=ry_a(1:nftsa,jrforce)
         fza(1:nftsa)=rz_a(1:nftsa,jrforce)

         fxrep(1:npoly)=rxbead(1:npoly,jrforce)
         fyrep(1:npoly)=rybead(1:npoly,jrforce)
         fzrep(1:npoly)=rzbead(1:npoly,jrforce)

!        constraint forces: actin-membrane tether, wall blocking, boundaries
         call constraints(nthreads,nmemb,npoly,nftsa,nxsol,nphisol,jmbsol,jsursol,jzsol,a2mem,a2fil,fil2a,endfil,nsurf, &
               kwall,wthick,lsqueez,l_mem,k_mem,l_a_z,l_mb_a,ktether, &
                xwall,ywall,zwall,xnorwall,ynorwall,znorwall,xsurf,ysurf,zsurf, &
                 xnorsurf,ynorsurf,znorsurf,xmemb,ymemb,zmemb,xboundmin,xboundmax, &
                  xcen,ycen,zcen,xa,ya,za,fxmemb,fymemb,fzmemb, &
                   fxrep,fyrep,fzrep,fxa,fya,fza)

!        interaction among membrane beads
         call membrane(npair_mb,pair_mb,pairpart,pairtyp,boundtyp,kpair,l_pair,l_mem, &
              kmemb,shift,xmemb,ymemb,zmemb,fxmemb,fymemb,fzmemb)

!        cell wall forces
         call wallforce(nxsol,nphisol,kgly,kpep,dxsol,pturgor,delta,invdelta,lwall,xwall,ywall,zwall,fywall,fzwall)

!        FtsZ bending
         call zangle(nfil,fstart,flen,beta,delta,invdelta,kz_thet,thet0_z,xcen,ycen,zcen,fxrep,fyrep,fzrep)

!        calculate distance between membrane and cell wall
         call checkdistance(nskip1,nskip2,nmemb,nphisol,nxsol,jmbsol,ndist,jdist,wthick,gap, &
                            xmemb,ymemb,zmemb,xwall,ywall,zwall,xnorwall,ynorwall,znorwall)
      end if

!     rigidity of FtsZ and FtsA:
      fxfil(1:npoly)=fxrep(1:npoly)
      fyfil(1:npoly)=fyrep(1:npoly)
      fzfil(1:npoly)=fzrep(1:npoly)

      call zrigid(nfil,npair_z,fstart,flen,pair_z,k_z,l_z,cutoff,rho,eps,rswitch,fswitch, &
           xcen,ycen,zcen,fxfil,fyfil,fzfil)
!----------------------------------------------------

!     update coordinates:
      if(jforce==1)then
         fmaxmb2=maxval(fxmemb*fxmemb+fymemb*fymemb+fzmemb*fzmemb)
      end if

      call newcoor(npoly,nftsa,nmemb,fmaxmb2,viscos,dtmodz,dtremod,runtime,fxmemb,fymemb,fzmemb, &
               fxa,fya,fza,fxfil,fyfil,fzfil,xmemb,ymemb,zmemb,xa,ya,za,xcen,ycen,zcen, &
                jforce,nxsol,nphisol,wallwidth,fywall,fzwall,xwall,ywall,zwall,xnorwall,ynorwall,znorwall)
!----------------------------------------------------

      if(runtime-oldtime>printtime)then
         dtime=dtime+runtime-oldtime
         oldtime=runtime

         write(*,16)jstep,runtime*0.000001,'sec'
         call writedcd(junit,nframe,natom,natom_zc,natom_a,natom_z2a,nxsol,nphisol,nmemb,nftsz,nfil,nftsa, &
                flen,a2mem,a2fil,fstart,xcen,ycen,zcen,xa,ya,za, &
                 xmemb,ymemb,zmemb,xwall,ywall,zwall)

         if(mod(nframe,5)==0)then
            dtime=dtime*0.000001
            cons_rate=(wallcons-oldwallcons)/dtime
            rwallcons=wallcons
            rtime=runtime*0.000001
            write(junit2,*)rtime,rwallcons,cons_rate
            oldwallcons=wallcons
            dtime=0.0d0
         end if

         if(nframe>=noutput)then
            close(junit)
            call ringout(nstart,natom,natom_zc,natom_a,natom_z2a,nwall,nxsol,nphisol,nmemb,nftsz,npoly,nfil, &
               zlen,nftsa,nftsamax,fstart,flen,filid,ringid,fil2a,gtp,a2mem,a2fil, &
                dxsol,dphisol,xboundmin,xboundmax,rwall,lwall,xwall,ywall,zwall,xnorwall,ynorwall,znorwall, &
                 xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf,xmemb,ymemb,zmemb,xcen,ycen,zcen,xa,ya,za)
            call system_clock(timerun)

            OPEN(junit1,FILE='restart.inp')
            write(junit1,'(a65)')'      nstart       jfile       jstart        time0        runtime'

            WRITE(junit1,*)NSTART,JFILE,jstep+1,(TIMERUN-TIMESTART)/rate+time0,runtime
            CLOSE(junit1)

               jfile=jfile+1
               call dcdheader(junit,jfile,natom)
               nframe=0
         end if
      end if
   end do
!------------------------------

   call system_clock(timerun)

   TIMERUN=(TIMERUN-TIMESTART)/rate+time0

   days=timerun/86400
   timerun=timerun-86400*days

   hours=timerun/3600
   timerun=timerun-3600*hours

   mins=timerun/60
   secs=timerun-60*mins

   write(*,1)'RUNNING TIME =',DAYS,'days : ',HOURS,'hours : ',MINS,'mins : ',secs,'secs'

1  FORMAT(A14,2X,I2,1X,A7,I2,1X,A8,1X,I2,1X,A7,1x,I2,1x,A4)

16 format(i16,2x,f10.5,1x,a3)
100 print*,'end of simulation'
end

