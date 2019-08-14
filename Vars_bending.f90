module declare

   implicit none

!  for setup:

   integer(kind=8)::count_rate,cr,timestart,time0,timerun,jstart
   integer days,hours,mins,secs,nthreads,jobid,junit1,nstart,jfile,noutput

   integer zlen,nftszstart,nfil,jtreadmil,intering,nphisol,nxsol,nwall,nbondwal,nw,nb,jx,jp
   integer nmeho,nhoop,nmemb,nm,jh,jm,nftsz,fstart0,length,length0,nzbond,n,jf,jm0
   integer npick,nbox,jb,nftsa,nftsamax,tension,nring
   integer nfilmax,nsunit,jgradient,jnogradient,jxwid
   integer npoly,nmono,nx,nx0,nphi,nphi0,jocc,jdir,jbend,jskip,nskip1,nskip2

   integer,allocatable,dimension(:,:)::bondwal,fill
   integer,allocatable,dimension(:)::fil2a,a2mem,a2fil,fstart,flen,filid,gtp,mark

   real(kind=8)::rate,runtime
   real(kind=8)::mbrad,mwid,rrad,rwid,pi,delta,invdelta,beta,k_z,l_z,l_diag,k_a,la_cen,l_mb_a,ktether,l_a_z
   real(kind=8)::kz_thet,kb,ltether,p_tether,l_mem,k_mem,l_z12,l_z34,dl_z,l_diag_new
   real(kind=8)::wthick,lsqueez,kwall,walrad,dphisol,dxsol,dyzsol,xmin,xmax,xzmin,xzmax
   real(kind=8)::xboundmin,xboundmax,shift,rad12,rad34,x14,x23,cosp,sinp
   real(kind=8)::x0,y0,z0,dphi,rho,eta,r,phi0,phi,dist,dist2,dx,dy,dz,d2,xbox
   real(kind=8)::xorigin,angle,wallwidth,viscos,fgrad,kwid,kgly,lpep,kpep,fturgor,pturgor,l0

   real(kind=8),allocatable,dimension(:,:)::xwall,ywall,zwall,xnorwall,ynorwall,znorwall,lwall
   real(kind=8),allocatable,dimension(:,:)::xnorsurf,ynorsurf,znorsurf,xsurf,ysurf,zsurf
   real(kind=8),allocatable,dimension(:)::xmemb,ymemb,zmemb,rwall
   real(kind=8),allocatable,dimension(:)::xcen,ycen,zcen,xa,ya,za!,x_a_old,y_a_old,z_a_old
   real(kind=8),allocatable,dimension(:,:)::xb,yb,zb,x_a_old,y_a_old,z_a_old
   real(kind=8),allocatable,dimension(:,:)::phi_start,phi_end


!-------------------------------------------

!  for dynamics

   integer(kind=8)::nstep,jstep
   integer natom,natom_z1,natom_z2,natom_z3,natom_z4,natom_zc,natom_a,natom_ab,natom_z2a

   integer nrforce,jrforce,nxclude,junit,junit2,nframe,jprint,jforce
   integer npair5_mb,npair_mb,npair5_z,npair_z
   integer check,j,nturnover,jcutoff,ndist

   integer,allocatable,dimension(:)::startfil,endfil
   integer,allocatable,dimension(:,:)::jmbsol,jsursol,jasol,jzsol,pair5_mb,pair_mb,pair5_z,pair_z
   integer,allocatable,dimension(:,:)::pairpart,nsurf,jdist
   integer,allocatable,dimension(:)::boundtyp,pairtyp

   real(kind=8)::ptread,phyd,padd,tadd,k_scale,k_scale_mb,kmemb!_thet,thet_memb
   real(kind=8)::r_off,rvdw,r_on,r_on2,r_off2,kexcl!,fvdwmax
   real(kind=8)::eps,cutoff,kpair,l_pair,thetpair,thet2by2,cos_t_2,gap,xgap,rmin
   real(kind=8)::printtime,oldtime,dt,dttread,dtadd,fmaxmb2,wallcons,dtime,oldwallcons
!   real(kind=8)::fa1max2,fa2max2,fa3max2,fa4max2,fb1max2,fb2max2,fb3max2,fb4max2,facenmax2
   real(kind=8)::fmax,kplane,growstep,kcir,kcircon,kradial

   real(kind=8),allocatable,dimension(:)::fxmemb,fymemb,fzmemb
   real(kind=8),allocatable,dimension(:,:)::fx_z,fy_z,fz_z,fx_a,fy_a,fz_a
   real(kind=8),allocatable,dimension(:)::fxrep_a,fyrep_a,fzrep_a,fxrep,fyrep,fzrep
   real(kind=8),allocatable,dimension(:,:)::rx_a,ry_a,rz_a,fywall,fzwall
   real(kind=8),allocatable,dimension(:,:)::rxmemb,rymemb,rzmemb,rxbead,rybead,rzbead
   real(kind=8),allocatable,dimension(:,:,:)::fxtem,fytem,fztem
   real(kind=8),allocatable,dimension(:)::fxa,fya,fza,gradient

   character(len=43)::a43

   real timecheck,timecheck0,rtime,rwallcons,cons_rate

end module
