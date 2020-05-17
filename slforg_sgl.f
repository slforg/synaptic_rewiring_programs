      program slforg_sgl
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
      parameter          (lec=82, lic=382, dx=55.0)
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
      write(6,*)'# Input synaptic connection file =?'
      read(5,*) inpf2
      open(57,file=inpf2,status='old')
c
      write(6,*)'Monte Caro Step =?'
      read(5,*) rdstp
      idstp=int(rdstp*real(nvc0)*rsnm)
c
      write(6,*)'Enter the number of executed steps =?'
      read(5,*) rnum
      rnn0=rnum*rsnm
      rmus00=rmus0*(rnn0/(rms+rnn0))**4
c
c ### Output Files ###
      write(6,*)'# Output synaptic connection file =?'
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
c set cut off area
      call setlv(lvx,1-ncff,nvcx+ncff,nvcx)
      call setlv(lvy,1-ncff,nvcy+ncff,nvcy)
c
c*******************************************************************
c     Cortical cell interaction matrix
c*******************************************************************
      z=0.0
      rcff=real(ncff)
      jx = ncff+1
      jy = jx
      jc = ncff*ncffw+ncff+1
      z0=0.0
      z1=0.0
      do 30 jd=1,ncffws
      dis2 = (jx-(mod(jd-1,ncffw)+1))**2 +
     &       (jy-((jd-1)/ncffw+1))**2

      vci(jd)=exp(-.5*dis2/aic)
      z1=z1+vci(jd)
      vce(jd)=exp(-.5*dis2/aec)
      z0=z0+vce(jd)
 30   continue
c
      do 35 jd=1,ncffws
      vce(jd)=vce(jd)/z0
      vci(jd)=vci(jd)/z1
 35   continue
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
c  orientation response magnitute correction
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
c     synapse connetion
c*********************************************************************
c
      jx0=ndnd+1
      jy0=jx0
c
c CPU time
c      st1=secnds(0.0)
c      write(6,*) 'start'
c
      do 170 i=1,idstp
c
      rnn0=rnum*rsnm+real(i)/real(nvc0)
      rmus=rmus0*(rnn0/(rms+rnn0))**4
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
c**** direction *****
c$omp parallel shared (dp0,dp1,i1x,i1y,dpch)
c$omp&  private(it0,idr)
c$omp do
      do 205 idr=1,idrc  
      do 205 it0=1, itsmp
      dp0(idr,it0)=dpch(i1x,i1y,idr,it0)
      dp1(idr,it0)=dpch(i1x,i1y,idr,it0)
 205  continue
c$omp end do nowait
c$omp end parallel
c$omp barrier
c
c save original location  
      ik1x0=ik1x
      ik1y0=ik1y
      ik2x0=ik2x
      ik2y0=ik2y
c location transform (ik1x,ik1y) & (ik2x,ik2y)       
      iclx=int((i1x-1)/2)+1
      icly=int((i1y-1)/2)+1     
      idx=iclx-ik1x
      ikx01=ik1x
      if(idx.lt.-nrfx0/2) ikx01=ik1x-nrfx0
      if(idx.gt. nrfx0/2) ikx01=ik1x+nrfx0
      ik1x=ikx01+nrfx0
      idy=icly-ik1y
      iky01=ik1y
      if(idy.lt.-nrfy0/2) iky01=ik1y-nrfy0
      if(idy.gt. nrfy0/2) iky01=ik1y+nrfy0
      ik1y=iky01+nrfy0
c      
      idx=iclx-ik2x
      ikx02=ik2x
      if(idx.lt.-nrfx0/2) ikx02=ik2x-nrfx0
      if(idx.gt. nrfx0/2) ikx02=ik2x+nrfx0
      ik2x=ikx02+nrfx0
      idy=icly-ik2y
      iky02=ik2y
      if(idy.lt.-nrfy0/2) iky02=ik2y-nrfy0
      if(idy.gt. nrfy0/2) iky02=ik2y+nrfy0
      ik2y=iky02+nrfy0
c
c response of cortical neuron
      ix=abs(jx0-jx1)+1
      iy=abs(jy0-jy1)+1
      rth1=0.0
      rth0=0.0
c$omp parallel shared (dp1,etl,ik1x,ik1y,ik2x,ik2y,
c$omp& isp1t,isp2t,vd0,ix,iy,idrc)
c$omp& private (it0,idr)
c$omp do
      do 206 idr=1,idrc  
      do 208 it0=1,itsmp
      dp1(idr,it0)=dp1(idr,it0)
     &+( etl(ik2x,ik2y,isp2t,it0,idr)
     &  -etl(ik1x,ik1y,isp1t,it0,idr))
 208  continue
 206  continue
c$omp end do nowait
c$omp end parallel
c$omp barrier
c
      do 209 idr=1,idrc  
      do 209 it0=1,itsmp
      rth1=rth1+dp1(idr,it0)
 209  continue
      rth1=rth1/real(idrc*itsmp)
c
c$omp parallel shared (j0,s,jx1,jy1,nrf,nrfx,
c$omp& iclx,icly,dmt,de0,vd,etl,z0,z1,ik1x,ik1y,ik2x,ik2y,zr0,
c$omp& dp1,dp0,dp11,dp01,ik1,ik2,isp1t,isp2t,de0n,
c$omp& rsh,i1x,i1y,rmax,rth1,a00,idrc)
c$omp& private (idx,idy,it0,idrss,j,jx,jy,jkd,ik0x,ik0y,ikx01,iky01,
c$omp& ist,ids1,ids2,idt1,idt2,jdx,jdy,idr,buf0,de01,z11,z01,idr0,z0)
c
c$omp do
      do 220 idr=1,idrc
      do 222 it0=1,itsmp
      buf0=dp0(idr,it0)-rsh(i1x,i1y)
      dp01(idr,it0)=0.5*(sign(1.0,buf0)+1.0)
      buf0=dp1(idr,it0)-rth1
      dp11(idr,it0)=0.5*(sign(1.0,buf0)+1.0)
 222  continue
 220  continue
c$omp end do nowait
c$omp barrier
c
c*** fixed values
      j0 = (i1-1)*nd0
c
c## energy of a cortical cell ##
c$omp single
      de0= 0.0
c$omp end single nowait
c$omp barrier
      de01=0.0
c
c$omp do
      do 90 jdy=1,ndndw
      do 90 jdx=1,ndndw
      j =j0+(jdy-1)*ndndw+jdx
      if(s(j).eq.0) goto 90
      jx=abs(jdx-jx1)+1
      jy=abs(jdy-jy1)+1
      jkd = mod(s(j)-1,nrf)
      ik0x= mod(jkd,nrfx0)+1
      ik0y= jkd/nrfx0+1
      idx=iclx-ik0x
      ikx01=ik0x
      if(idx.lt.-nrfx0/2) ikx01=ik0x-nrfx0
      if(idx.gt. nrfx0/2) ikx01=ik0x+nrfx0
      ik0x=ikx01+nrfx0
      idy=icly-ik0y
      iky01=ik0y
      if(idy.lt.-nrfy0/2) iky01=ik0y-nrfy0
      if(idy.gt. nrfy0/2) iky01=ik0y+nrfy0
      ik0y=iky01+nrfy0   
      ist=int(real(s(j)-1)/real(nrf))+1
c  
c ### Hebbian learning term
      do 92 idr0=1, 2
      idrss=(idr0-1)*(idrc/2)+1   
      do 92 it0=1, itsmp
        z01=etl(ik0x,ik0y,ist,it0,idrss)*(
     &  etl(ik2x,ik2y,isp2t,it0,idrss)*dp11(idrss,it0)
     & -etl(ik1x,ik1y,isp1t,it0,idrss)*dp01(idrss,it0))
        z11=etl(ik0x,ik0y,ist,it0,idrss)*(
     &  etl(ik2x,ik2y,isp2t,it0,idrss)
     & -etl(ik1x,ik1y,isp1t,it0,idrss))
      de01=de01-((1+a00)*z01-a00*z11)
c
 92   continue
c
 90   continue
c
c$omp end do nowait
c$omp critical
      de0=de0+de01
c$omp end critical
c$omp end parallel
c
c## Energy of cortical interaction ##
c$omp single
      de11=0.0
c
c$omp end single nowait
c$omp barrier
c
c$omp parallel shared (vs00,rdih,etl,ik1x,ik1y,isp1t,dp11,dp01,
c$omp& ik2x,ik2y,isp2t,lvx,lvy,i1x,vc,dpch,rsh,zr0,de11,idrc,
c$omp& stb,vce,vci)
c$omp& private (idrss,it0,vs0,z0,idx,idy,idx0,idy0,vc0,rid,ridt,
c$omp& de111,idr0,z0,z01,z11)
      de111=0.0
c$omp do
      do 112 idr0=1,2
      idrss=(idr0-1)*(idrc/2)+1
      do 112 it0=1,itsmp
      z01=(etl(ik2x,ik2y,isp2t,it0,idrss)*dp11(idrss,it0)
     &    -etl(ik1x,ik1y,isp1t,it0,idrss)*dp01(idrss,it0))
      z11=(etl(ik2x,ik2y,isp2t,it0,idrss)
     &    -etl(ik1x,ik1y,isp1t,it0,idrss))
c
c*** cortical interaction 
      do 112 idy=1,ncffw
      idy0=lvy(i1y-ncff-1+idy)
      do 112 idx=1,ncffw
      idx0=lvx(i1x-ncff-1+idx)
      vc0=( vce((idy-1)*ncffw+idx)-vci((idy-1)*ncffw+idx) )
c
      rid=dpch(idx0,idy0,idrss,it0)-rsh2(idx0,idy0)
      ridt=rid*0.5*(sign(1.0,real(rid))+1.0)
      de111=de111-((1+a00)*z01-a00*z11)*ridt*vc0
c
 112  continue
c$omp end do nowait
c$omp critical
      de11=de11+de111
c$omp end critical
c$omp barrier
c$omp end parallel
c
      de=(raff*v00*de0+rctx*de11)/rave
c
      des=rmus
c
c--------------------------------------------------------------------
c
 120    w =1.0/(1.0+dexp(dble(beta)*de))
        ws=2.0/(1.0+dexp(dble(beta)*des))
 131    call random_number(rv0)
        if(rv0.ge.1.0) goto 131
	if (ws .lt. rv0) goto 170
 130    call random_number(rv0)
        if(rv0.ge.1.0) goto 130
	if (w .lt. rv0) goto 170
          s(j1) = k2
c
c*** membrane potential   
      do 216 idr=1,idrc
c      idr=(idr0-1)*(idrc/2)+1 
      do 218 it0=1,itsmp
      dpch(i1x,i1y,idr,it0)=dpch(i1x,i1y,idr,it0)
     &+( etl(ik2x,ik2y,isp2t,it0,idr)
     &  -etl(ik1x,ik1y,isp1t,it0,idr))
 218  continue
 216  continue
c
c*** threshold
      rsh(i1x,i1y)=rth1
c
 170  continue
c
c CPU time
c      write(6,*) 'stop'
c      eltm=secnds(st1)
c      write(*,*) eltm
c
c-----------------------------------------------------------------------
c I/O
c-----------------------------------------------------------------------
c
 500  do 51 j=1, nvc
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
      write(6,*) 'value of mu ',rmus00,rmus,ws
      write(6,*) '--------------------------------'
c
c-----------------------------------------------------------------------
c
       stop
       end
c
c
c*********************************************************************
c     cutt off area set up
c*********************************************************************
      subroutine setlv(lv,ls,le,n)
      integer   lv(ls:le)
      do 1 i = ls,0
	 lv(i) = i + (1-int(i/n))*n
    1 continue
      do 2 i = 1,le
	 lv(i) = i - int((i-1)/n)*n
    2 continue
      return
      end
c
c*********************************************************************
c     randomiztion routine
c*********************************************************************
c *option* -NO
      function random(ix0)
      common/rand1/m(521),j
      dimension ia(521)
c
      if (ix0 .ne. 0) then
      ix = ix0
      do i=1,521
      ix=ix*69069
      ia(i)=isign(1,ix)
      end do
      do j=1,521
      ih=mod((j-1)*32,521)+1
      mj=0
      do i=1,31
      ii=mod(ih+i-2,521)+1
      mj=2*mj+(ia(ii)-1)/(-2)
      ij=mod(ii+488,521)+1
      ia(ii)=ia(ii)*ia(ij)
      end do
      m(j)=mj
      ii=mod(ih+30,521)+1
      ij=mod(ii+488,521)+1
      ia(ii)=ia(ii)*ia(ij)
      end do
      j=0
      endif
      j=j+1
      if (j.gt.521) j=1
      k=j-32
      if(k.le.0) k=k+521
      m(j)=ieor(m(j),m(k))
      random=float(m(j))*0.4656613e-9
      return
      end
