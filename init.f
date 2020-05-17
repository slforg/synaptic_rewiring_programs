      program init
c
c*********************************************************************
c # make new initial geniculo-cortical input 
c     cortex field (48*48)
c     retinae (24*24)c  
c
c                copy right         NTC: M.Miyashita,  ECU: S.Tanaka
c*********************************************************************
      implicit real(a-h,k-l,o-s,u-z)
      implicit real*8(t)
c
      real                qi,qe,dx,beta,qir,ler,qe1,qi1,qer,
     &                    qe0,le0,qec,qic,lec,lic
      integer             nvcx0,nvcy0,ncff,
     &                    nrfx0,nrfy0,nrf0,nvc0
c
c cortex size
      parameter          (nvcx0=48,nvcy0=nvcx0,nvc0=nvcx0*nvcy0)
c LGN size & expand LGN area for stimulation
      parameter          (nrfx0=24,nrfy0=nrfx0,nrf0=nrfx0*nrfy0,
     &                    lgnx=nrfx0*3,lgny=lgnx,lgox=25,lgoy=25)
c cortical interaction size  
      parameter          (ncff=9,ncffw=2*ncff+1,nrfw=2*nrfx0+1,
     &                    ncffws=ncffw**2,nrfwsq=nrfw**2)
c dendritic field 
      parameter          (ndnd=9,icdt=3,
     &                    icxd=(nvcx0-1)*icdt+1+ndnd*2,
     &                    ndndw=(2*ndnd)+1,nd0=ndndw*ndndw)
c number of LGN type
c on-nonlagged, on-lagged, off-nonlagged, off-lagged
      parameter          (itc=4, nkc=4)
c define value of pi
      parameter          (pi=3.141592653)
c # LGN RF temporal parameter
c rlmt: Decay constant of temporal RF
c rtu0: Latency of lagged response
c ft0:  Temporal frequency of RF
c rtf0: Time window of temporal RF
c ftstm: Stimulus frequency
c itsmp & dt: time step 
      parameter          (rlmt0=0.1, rtu0=0.0569725, ft0=5.0, 
     &                    rtf0=1.0/ft0, itsmp=24, scl=1.0,
     &                    ftstm=4.0,
     &                    dt=(1.0/ftstm)/real(itsmp))
c # LGN RF spatial parameter
c ler: Extent of center subfield
c lir: Extent of surround sub field
c dxl: 1 pixcel 0.5[degree]
      parameter          (ler=0.225, lir=0.9, dxl=0.5)
c # Cortical interaction
c lec: Extent of excitatory interaction function
c lic: Extent of inhibitory interaction function
c dx:  Distance between nearest neighbor neurons
      parameter          (lec=82, lic=382, dx=50.0)
c
c rctx:  Strength of cortical interaction
c a00:   Efficiency ratio of LTD to LTP
c rmus0: Maximum amount of p75 NTR expression
c rms:   Time window that enables synapses to be rewired
      parameter          (raff=1.0,rctx=4.0, a00=2.0, 
     &                    rmus0=0.8, rms=4000)
c beta: Fictitious inverse temperature
      parameter          (beta=250.0)
c
c rsf1: Spatial frequency of gratings
c Spatial frequency = 1.0/(rsf1*dxl)
      parameter          (rsf1=5.0)
c stimulus directions (every 15[deg], 24 directions)
      parameter          (idrc=24, idrh=idrc/2)
c
c rsm: number of synapse site
c v00: Contribution of afferent inputs to membrane potential
c rsnm: Number of synapse site
      parameter          (rsnm=254, v00=0.045)
c
      integer             i,icn,iexch,ik1,ik1x,ik1y,ik2,ik2x,ik2y,
     &                    ikd1,ikd2,
     &                    istpt,itime,j,j1,i1x,i1y,jc,jd,jdx,
     &                    jdy,jkd,jrf,jx,jy,
     &                    k1,k1a,k2,
     &                    lvx(1-ncff:nvcx0+ncff),
     &                    lvy(1-ncff:nvcy0+ncff),
     &                    ncffw, nrfw, ncffws,
     &                    nrfwsq,nrfxh,s(nvc0*nd0),
     &                    dnf(nvcx0,ndndw),     
     &                    spo(nvc0*nd0),spt(nvc0*nd0)
      real                ae,aer,ai,
     &                    cr1,cr2,cr3,dis2,rkd2,le,li,
     &                    rv,sn,rcrx,rcry,
     &                    stuo0, sttl0,xyvc,xyrf
      double precision    etl(lgnx,lgny,itc,itsmp,idrc),
     &                    dpch(nvcx0,nvcy0,idrc,itsmp+1),
     &                    dp0(idrc,itsmp),dp1(idrc,itsmp),
     &                    dp01(idrc,itsmp),dp11(idrc,itsmp),
     &                    dp02(idrc,itsmp),dp12(idrc,itsmp),
     &                    rsh(nvcx0,nvcy0),
     &                    rsh2(nvcx0,nvcy0),
     &                    de,de0,de1,w,ws,de11,de12,de2,p2,de3,
     &                    de01,de111,de02,de112,de0s,de11s,des,
     &                    vd(ndndw,ndndw),vdi(ndndw,ndndw),
     &                    vd0(ndndw,ndndw),vc((2*ncff+1)**2),
     &                    rtcl,rtcl2,cpt(4,4),
     &                    z0,z1,z2,z3,z4,z5,z6,z7,z8,
     &                    p0,p3,rw0,rw1,p1,rw2,
     &                    znm,znm2,rid,ridt,ridt0,ridt1,
     &                    dcx(icxd,icxd), vcdn(icxd,icxd,2),
     &                    rdih(idrc),t1,bin(idrc),bin2(idrc),
     &                    de0n,de1n,rw0n,rw1n,p0n,p1n,
     &                    z11,z01,z02,z12,
     &                    vce((2*ncff+1)**2),vci((2*ncff+1)**2),
     &                    stb(ncffw,ncffw,nvcx0,nvcy0),
     &                    rdh(lgnx,lgny,itc,idrc),rave
c
c for random seed
      integer :: seedsize, ii
      integer, allocatable :: seed(:)
c 
      real*8              engy,et, tstrt,tstrt1
      character           inpf*20,inpf2*20,inpf3*20,inpf4*20,inpf5*20,
     &                    outf*20,outf2*20,outf3*20
      real*4              st1
c--------------------------------------------------------------------
c
c for random number call
      call random_seed(size=seedsize)
      allocate(seed(seedsize))
      do ii = 1, seedsize
         call system_clock(count=seed(ii))
      end do
      call random_seed(put=seed(:))
c
c average factor
      rave=real(idrc*itsmp)
c
c ### Input Files ###
      write(6,*)'# Input  old initial connection file =?'
      read(5,*) inpf2
      open(57,file=inpf2,status='old')
c
c ### Output Files ###
      write(6,*)'# Output new initial connection file =?'
      read(5,*) outf
      open(59,file=outf,status='unknown')
c
c*********************************************************************
c set common values
      rcdt=real(icdt)
      rdnd=real(ndnd)
      rdnd2=rdnd**2+1.0**2
      idndc=ndnd+1
      rvcx=real(nvcx0)
      rvcy=real(nvcx0)
c
      nvc=nvcx0*nvcy0
      xyvc=real(nvc0)
      nrf=nrfx0*nrfy0
      xyrf=real(nrf)
      rcrx=real(nvcx0/nrfx0)
      rcry=real(nvcy0/nrfy0)
c
      nrfx=nrfx0
      nrfy=nrfy0
      nvcx=nvcx0
      nvcy=nvcy0
c
c*********************************************************************
c     old spin configuration
c*********************************************************************
	do 50 j=1, nvc
	do 50 i=1, nd0
        read(57,*) s((j-1)*nd0+i)
 50     continue
        close(57)
c
c-----------------------------------------------------------------
c
      idstp=100*int(rsnm)*nvc0
c
      do 170 i=1,idstp
c
 70   call random_number(rv)
      if(rv.ge.1.0) goto 70
      i1=int(rv*xyvc)+1
      i1x=mod(i1-1,nvcx)+1
      i1y=(i1-1)/nvcx+1
c
 81   call random_number(rv)
      if(rv.ge.1.0) goto 81
      j1=int(rv*real(nd0))+1
      jx1=mod(j1-1,ndndw)+1
      jy1=(j1-1)/ndndw+1
      j1=(i1-1)*nd0+j1
      k1=s(j1)
      if(k1.le.0.0) goto 81
      rjx=real(i1x)+real(jx1-jx0)/rcdt
      rjy=real(i1y)+real(jy1-jy0)/rcdt
      ik1=mod(k1-1,nrf)+1
      ik1x=mod(ik1-1,nrfx)+1
      ik1y=(ik1-1)/nrfx+1
      isp1t=(k1-1)/nrf+1
c
 80   call random_number(rv)
      if(rv.ge.1.0) goto 80
      ik2=int(rv*nrf0)+1
      ik2x=mod(ik2-1,nrfx)+1
      ik2y=(ik2-1)/nrfx+1
      rk2x=real(ik2x)*rcrx+0.5
      rk2y=real(ik2y)*rcry+0.5
      rax=abs(rk2x-rjx)
      raxd=abs(rk2x-real(nvcx)-rjx)
      if(rax.gt.raxd) rax=raxd
      raxd=abs(rk2x+real(nvcx)-rjx)
      if(rax.gt.raxd) rax=raxd
      ray=abs(rk2y-rjy)
      rayd=abs(rk2y-real(nvcy)-rjy)
      if(ray.gt.rayd) ray=rayd
      rayd=abs(rk2y+real(nvcy)-rjy)
      if(ray.gt.rayd) ray=rayd
      rkd2=rax**2 + ray**2
c
      arbf=rkd2
      rv=9.0**2
      if(arbf.gt.rv) goto 80
c
 83   call random_number(rv)
      if(rv.ge.1.0) goto 83
      sp2t=aint(rv*real(itc))+1.0
      isp2t=int(sp2t)
      k2 =ik2+nrf*(isp2t-1)
c
c**** retry random number *****
      if(s(j1).eq.k2) goto 81
c
      s(j1) = k2
c
 170  continue
c
c-----------------------------------------------------------------------
c I/O
c-----------------------------------------------------------------------
c
      do 51 j=1, nvc
      do 51 i=1, nd0
      ix=mod(i-1,ndndw)+1
      iy=(i-1)/ndndw+1
      rr=real((ix-idndc)**2+(iy-idndc)**2)
      if(rr.gt.rdnd2) s((j-1)*nd0+i)=0 
      write(59,*) s((j-1)*nd0+i)
 51   continue
c
 900  close(59)
c
c-----------------------------------------------------------------------
c
       stop
       end
