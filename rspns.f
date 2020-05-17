      program rspns
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
c stimulus directions (360 directions)
      parameter          (idrc=360, idrh=idrc/2)
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
     &                    drsp(nvcx0,nvcy0,idrc),
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
c for randomnumber call
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
      write(6,*)'# Input synaptic connection file =?'
      read(5,*) inpf2
      open(57,file=inpf2,status='old')
c ### Output Files ###
      write(6,*)'# Output direction response data  (drsp#) =?'
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
c
c*****************************************************************
c Responses of LGN cells
c*****************************************************************
c
      rxc=real(lgnx)-0.5
      ryc=real(lgny)-0.5
c
c###### Optimal tmporal frequency stimuli #####
c
      ft=ftstm
c
c** common part
      rlm2=rlmt0**2
      a1=rlm2/(1.0+(2.0*pi*(ft+ft0)*rlmt0)**2)
      a2=rlm2/(1.0+(2.0*pi*(ft0-ft)*rlmt0)**2)
      b0=exp(-rtf0/rlmt0)
      s0p=sin(2.0*pi*(ft0+ft)*rtf0)
      c0p=cos(2.0*pi*(ft0+ft)*rtf0)
      s0n=sin(2.0*pi*(ft0-ft)*rtf0)
      c0n=cos(2.0*pi*(ft0-ft)*rtf0)
      b11=b0*(-s0p/rlmt0-2.0*pi*(ft0+ft)*c0p)+2.0*pi*(ft0+ft)
      b12=b0*(-s0n/rlmt0-2.0*pi*(ft0-ft)*c0n)+2.0*pi*(ft0-ft)
      b13=b0*(-c0p/rlmt0+2.0*pi*(ft0+ft)*s0p)+1.0/rlmt0
      b14=b0*(-c0n/rlmt0+2.0*pi*(ft0-ft)*s0n)+1.0/rlmt0
      b21=   ( s0p/rlmt0-2.0*pi*(ft0+ft)*c0p)+b0*(2.0*pi*(ft0+ft))
      b22=   ( s0n/rlmt0-2.0*pi*(ft0-ft)*c0n)+b0*(2.0*pi*(ft0-ft))
      b23=   ( c0p/rlmt0+2.0*pi*(ft0+ft)*s0p)-b0*(1.0/rlmt0)
      b24=   ( c0n/rlmt0+2.0*pi*(ft0-ft)*s0n)-b0*(1.0/rlmt0)
      c11=a1*b11+a2*b12
      c12=a1*b13-a2*b14
      c21=a1*b21+a2*b22
      c22=a1*b23-a2*b24
c temporal part
      do 410 it0=1, itsmp
      rtn=2.0*pi*ft*real(it0-1)*dt
      rtl=rtn-2.0*pi*ft*rtu0
c spatial part
      do 400 iy0=1, lgny
      do 400 ix0=1, lgnx
      rx=real(ix0-1)-rxc
      ry=real(iy0-1)-ryc
      do 400 ith=1, idrc
c
      d01=1.0
      rth=real(ith-1)*2.0*pi/real(idrc)
      rss=rx*cos(rth)+ry*sin(rth)
c ****** LGN respnse
      swm01=2.0*pi/rsf1
      cs01n= 0.5*cos(swm01*rss-rtn)
      ss01n=-0.5*sin(swm01*rss-rtn)
      cs01l= 0.5*cos(swm01*rss-rtl)
      ss01l=-0.5*sin(swm01*rss-rtl)
      ss= d01*(cs01n*c11-ss01n*c12)
      etl(ix0,iy0,1,it0,ith)=ss*0.5*(sign(1.0,ss)+1.0)
      ss=-d01*(cs01l*c21-ss01l*c22)
      etl(ix0,iy0,2,it0,ith)=ss*0.5*(sign(1.0,ss)+1.0)
      ss=-d01*(cs01n*c11-ss01n*c12)
      etl(ix0,iy0,3,it0,ith)=ss*0.5*(sign(1.0,ss)+1.0)
      ss= d01*(cs01l*c21-ss01l*c22)
      etl(ix0,iy0,4,it0,ith)=ss*0.5*(sign(1.0,ss)+1.0)
 400  continue
 410  continue
c
c  orientation response magnitude correction
      do 311 id0=1,idrc
      do 311 i=1,nkc
      do 311 iyd=1,lgny
      do 311 ixd=1,lgnx
         rdh(ixd,iyd,i,id0)=0.0
 311  continue
c
      do 316 id0=1,idrc
      do 316 i=1,nkc
      do 316 iyd=1,lgny
      do 316 ixd=1,lgnx
         z0=0.0
      do 315 it0=1,itsmp
      z0=z0+etl(ixd,iyd,i,it0,id0)
 315  continue
      rdh(ixd,iyd,i,id0)=z0
 316  continue
c
      do 317 id0=1,idrc
      do 317 it0=1,itsmp
      do 317 i=1,nkc
      do 317 iyd=1,lgny
      do 317 ixd=1,lgnx
      etl(ixd,iyd,i,it0,id0)=etl(ixd,iyd,i,it0,id0)
     &    /rdh(ixd,iyd,i,id0)
 317  continue
c 
c*********************************************************************
c     initial spin configuration
c*********************************************************************
c
      do 50 j=1, nvc0
      do 50 i=1, nd0
        read(57,*) s((j-1)*nd0+i)
 50   continue
      close(57)
c
c*********************************************************************
c     responses of cortical cells
c*********************************************************************
c set LGN type table
      do 60 j=1, nvc0
      do 60 i=1, nd0
      is0=s((j-1)*nd0+i)
      if(is0.eq.0) goto 61
      spo((j-1)*nd0+i)=mod((is0-1),nrf0)+1
      spt((j-1)*nd0+i)=(is0-1)/nrf0+1
      goto 60
 61   spo((j-1)*nd0+i)=0
      spt((j-1)*nd0+i)=0
 60   continue
c
      do 250 iys=1, nvcy0
      do 250 ixs=1, nvcx0
      i0=(iys-1)*nvcx0+ixs
      is0=(i0-1)*nd0
      iclx=int((ixs-1)/2)+1
      icly=int((iys-1)/2)+1
c
      do 200 i=1,itsmp
      do 202 idr=1,idrc
      buf=0.0
      do 204 iyd=1,ndndw
      do 204 ixd=1,ndndw
      i0=(iyd-1)*ndndw+ixd
      imd=spt(is0+i0)
      if(imd.eq.0) goto 204
      ix=abs(ndnd+1-ixd)+1
      iy=abs(ndnd+1-iyd)+1
      ik0=spo(is0+i0)
      ikx=mod(ik0-1,nrfx0)+1
      iky=(ik0-1)/nrfx0+1
c      
      idx=iclx-ikx
      ikx0=ikx
      if(idx.lt.-nrfx0/2) ikx0=ikx-nrfx0
      if(idx.gt. nrfx0/2) ikx0=ikx+nrfx0
      ikx=ikx0+nrfx0
      idy=icly-iky
      iky0=iky
      if(idy.lt.-nrfy0/2) iky0=iky-nrfy0
      if(idy.gt. nrfy0/2) iky0=iky+nrfy0
      iky=iky0+nrfy0
c
      buf=buf+etl(ikx,iky,imd,i,idr)
 204  continue
      dpch(ixs,iys,idr,i)=buf
 202  continue
 200  continue
c
 250  continue
c
c *** threshold *** 
      do 136 iyd=1,nvcy0
      do 136 ixd=1,nvcx0
      rth1=0.0
      do 135 id=1,idrc  
      do 135 it=1,itsmp
      rth1=rth1+dpch(ixd,iyd,id,it)
 135  continue
      rsh(ixd,iyd)=rth1/real(idrc*itsmp)
 136  continue
c
      j00=ndnd+1
c
c*********************************************************************
c     direction rsponses
c*********************************************************************
c direction preference
      do 70 iy=1, nvcy0
      do 70 ix=1, nvcx0
      do 80 id=1, idrc
         z1=0.0
      do 90 it=1, itsmp
         z2=dpch(ix,iy,id,it)-rsh(ix,iy)
         z1=z1+z2*0.5*(sign(1.0,real(z2))+1.0)
 90   continue
      drsp(ix,iy,id)=z1/real(itsmp)
 80   continue
 70   continue
c
c-----------------------------------------------------------------------
c I/O
c-----------------------------------------------------------------------
c
      do 131 iys0=1,nvcy0
      do 131 ixs0=1,nvcx0
      do 131 i=1,idrc
      write(59,*) drsp(ixs0,iys0,i)
  131 continue
c
  900 close(59)
c
c-----------------------------------------------------------------------
c
       stop
       end

