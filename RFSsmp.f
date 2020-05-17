      program RFSsmp
c
c*********************************************************************
c #Selforganization of spatiotemporal cortical RF
c     cortex field (48*48)
c     retinae (24*24)
c     moving grating input
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
      parameter          (itc=4,nkc=4)
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
c dxl: 1 pixel 0.5[degree] 
      parameter          (ler=0.225, lir=0.9, dxl=0.5)
c # Cortial interaction
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
c For output of RF data
c every 10 msec
      parameter          (dt0=0.01,itc0=27)
      real                cp(nrfx0,nrfy0),cpt(nkc,itc0),
     &                    sprf(itc0,nrfx0,nrfy0) 
c--------------------------------------------------------------
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
     &                    rtcl,rtcl2,
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
c 
      real*8              engy,et, tstrt,tstrt1
      character           inpf*20,inpf2*20,inpf3*20,inpf4*20,inpf5*20,
     &                    outf*20,outf2*20,outf3*20
      real*4              st1
c--------------------------------------------------------------------
c
c ### Input Files ###
      write(6,*)'# Input synaptic connection file =?'
      read(5,*) inpf
      open(56,file=inpf,status='old')
c
      write(6,*)'# Cortex cell location =(X,Y)?'
      read(5,*) jcx, jcy
c
c ### Output Files ###
      write(6,*)'# Output spatio-tempral RF file =?'
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
c## Cortical interaction ##
      aec=(lec/dx)**2
      aic=(lic/dx)**2
      qec=0.5*qec/aec/pi
      qic=0.5*qic/aic/pi
c
c## Spatial part of LGN RF ## 
      qer=1.0
      qir=1.0
      aer=(ler/dxl)**2
      air=(lir/dxl)**2
      cr1=0.5*qer/pi/aer
      cr2=0.5*qir/pi/air
c
c*********************************************************************
c      Spatial component
c*********************************************************************
      jx0=1
      jy0=1
      do 40 jy1=1, nrfx0
      do 40 jx1=1, nrfx0
	dis2=real((jx0-jx1)**2+(jy0-jy1)**2)
	cp(jx1,jy1)=cr1*exp(-0.5*dis2/aer)-cr2*exp(-0.5*dis2/air) 
 40   continue
c
c*********************************************************************         
c     temporal component
c*********************************************************************
c
c non-lagged cell
      do 42 jt=1, itc0
      atm=real(jt-1)*dt0
      if(atm.gt.rtf0) goto 43
      cpt(1,jt)= exp(-atm/rlmt0)*sin(2.0*pi*ft0*atm)
      cpt(3,jt)=-exp(-atm/rlmt0)*sin(2.0*pi*ft0*atm)
      goto 42
 43   cpt(1,jt)=0.0
      cpt(3,jt)=0.0
 42   continue
c lagged cell
      do 44 jt=1, itc0
      atm=real(jt-1)*dt0
      if(atm.lt.rtu0) goto 45
      if(atm.gt.rtu0+rtf0) goto 45
      cpt(2,jt)=-exp(-(rtf0+rtu0-atm)/rlmt0)*sin(2.0*pi*ft0*(atm-rtu0))
      cpt(4,jt)= exp(-(rtf0+rtu0-atm)/rlmt0)*sin(2.0*pi*ft0*(atm-rtu0))
      goto 44
 45   cpt(2,jt)=0.0
      cpt(4,jt)=0.0
 44   continue
c
c*********************************************************************
c     Input connection file
c*********************************************************************
	do 50 j=1, nvc0
	do 50 i=1, nd0
           read(56,*) s((j-1)*nd0+i)
 50     continue
        close(56)
c
c*********************************************************************
c     Receptive Field configuration
c*********************************************************************
      do 60 j=1,nvc0
      do 60 i=1,nd0
       spo((j-1)*nd0+i)=mod((s((j-1)*nd0+i)-1),nrf)+1
       spt((j-1)*nd0+i)=(s((j-1)*nd0+i)-1)/nrf+1
 60   continue
c
      do 65 jt=1,itc0
      do 65 j1=1,nrfy0
      do 65 j2=1,nrfx0
	sprf(jt,j1,j2)=0.0
 65   continue
c
      jc=jcx+(jcy-1)*nvcx0
      j0=(jc-1)*nd0
c
      do 70 i=1,itc0
      do 80 iy=1,nrfy0 
      do 80 ix=1,nrfx0
      do 100 jdy=1, ndndw
      do 110 jdx=1, ndndw
	j =j0+(jdy-1)*ndndw+jdx
	if(s(j).eq.0) goto 110
        ias = int(abs(s(j)))
        jk0 = mod((ias-1),nrf0)
        jkx = mod(jk0,nrfx0)+1
        jky = int(jk0/nrfx0)+1
        ix0=abs(ix-jkx)+1
        iy0=abs(iy-jky)+1
        jt=int((ias-1)/nrf0)+1
        sprf(i,ix,iy)=sprf(i,ix,iy)+cp(ix0,iy0)*cpt(jt,i)/rsnm
 110  continue
 100  continue
 90   continue
 80   continue
 70   continue
c
c---------------------------------------------------------------------
c
      do 130 i=1, itc0
      do 130 iy=1,nrfy0
      do 130 ix=1,nrfx0
      write(59,*) sprf(i,ix,iy)
  130 continue
      close(59)
c
c-----------------------------------------------------------------------
c
       stop
       end
