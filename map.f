      program map
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
c
c stimulus directions (360 directions)
      parameter          (idrc=360, iori=180) 
c
c rsm: number of synapse site
c v00: Contribution of afferent inputs to membrane potential
c rsnm: Number of synapse site
      parameter          (rsnm=254, v00=0.045)
c
      integer             bin(6)
      real                drsp(nvcx0,nvcy0,idrc),
     &                    drsp2(nvcx0,nvcy0,iori),
     &                    cs(nvcx0,nvcy0),sn(nvcx0,nvcy0),
     &                    rori(nvcx0,nvcy0),rmg(nvcx0,nvcy0),
     &                    rdir(nvcx0,nvcy0)
c
      real*8              engy,et, tstrt,tstrt1
      character           inpf*20,inpf2*20,inpf3*20,inpf4*20,inpf5*20,
     &                    outf*20,outf2*20,outf3*20
      real*4              st1
c
c--------------------------------------------------------------------
c
c ### Input Files ###
      write(6,*)'# response file (drsp#) =?'
      read(5,*) inpf2
      open(57,file=inpf2,status='old')
c
c ### Output Files ###
      write(6,*)'# Output dir. ori. magnitude data (dmap#) =?'
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
c     input response file (for direction)
c*********************************************************************
      do 100 iys0=1,nvcy0
      do 100 ixs0=1,nvcx0
      do 100 i=1,idrc
      read(57,*) drsp(ixs0,iys0,i)
 100  continue
c
c direction
      do 50 iys0=1,nvcy0
      do 50 ixs0=1,nvcx0
         z0=0.0
         id0=0.0
      do 60 i=1,idrc
         if(z0.gt.drsp(ixs0,iys0,i)) goto 60
         z0=drsp(ixs0,iys0,i)
         id0=i
 60   continue
      rdir(ixs0,iys0)=real(id0)
 50   continue
c
c orientaion
      do 110 iys0=1,nvcy0
      do 110 ixs0=1,nvcx0
      do 120 i=1,90
         drsp2(ixs0,iys0,90+i)=drsp(ixs0,iys0,i)
 120  continue
      do 130 i=91,180
         drsp2(ixs0,iys0,i-90)=drsp(ixs0,iys0,i)
 130  continue
      do 140 i=181,270
         drsp2(ixs0,iys0,i-90)=drsp2(ixs0,iys0,i-90)+drsp(ixs0,iys0,i)
 140  continue
      do 150 i=271,idrc
         drsp2(ixs0,iys0,i-270)=drsp2(ixs0,iys0,i-270)+drsp(ixs0,iys0,i)
 150  continue
c
 110  continue
c
c*********************************************************************
c     orientaion information from vector sum method
c*********************************************************************
      do 200 iys0=1,nvcy0
      do 200 ixs0=1,nvcx0
         cs(ixs0,iys0)=0.0
         sn(ixs0,iys0)=0.0
 200  continue
c
      do 210 iys0=1,nvcy0
      do 210 ixs0=1,nvcx0
      do 220 i=1,iori     
         rdg=real(i-1)*pi/180.0
         cs(ixs0,iys0)=cs(ixs0,iys0)+drsp2(ixs0,iys0,i)*cos(2.0*rdg)
         sn(ixs0,iys0)=sn(ixs0,iys0)+drsp2(ixs0,iys0,i)*sin(2.0*rdg)
 220  continue
      cs(ixs0,iys0)=cs(ixs0,iys0)/real(iori)
      sn(ixs0,iys0)=sn(ixs0,iys0)/real(iori)
      rmg(ixs0,iys0)=sqrt(cs(ixs0,iys0)**2 + sn(ixs0,iys0)**2)
 210  continue
c
      do 230 iys0=1,nvcy0
      do 230 ixs0=1,nvcx0
         z0=cs(ixs0,iys0)
         z1=sn(ixs0,iys0)
         z2=atan(z1/z0)/2.0
         if(z0.lt.0.0) goto 250
         if(z1.lt.0.0) goto 240        
         rori(ixs0,iys0)=z2*180.0/pi
         goto 230
 240     rori(ixs0,iys0)=z2*180.0/pi+180.0
         goto 230
 250     rori(ixs0,iys0)=z2*180.0/pi+90.0
         goto 230
 260     rori(ixs0,iys0)=90.0
 230     continue
c
         z0=0.0
         do 300 iys0=1,nvcy0
         do 300 ixs0=1,nvcx0      
            if(rmg(ixs0,iys0).gt.z0) z0=rmg(ixs0,iys0)
 300     continue
         write(6,*) 'max magnitude', z0
c
         do 1000 i=1, 6
            bin(i)=0
 1000    continue
         do 1010 iys0=1,nvcy0
         do 1010 ixs0=1,nvcx0
            ir=rori(ixs0,iys0)
            if(ir.ge.165) bin(1)=bin(1)+1
            if(ir.ge.15) goto 1020
            bin(1)=bin(1)+1
            goto 1010
 1020       if(ir.ge.45) goto 1030
            bin(2)=bin(2)+1
            goto 1010
 1030       if(ir.ge.75) goto 1040
            bin(3)=bin(3)+1
            goto 1010
 1040       if(ir.ge.105) goto 1050
            bin(4)=bin(4)+1
            goto 1010
 1050       if(ir.ge.135) goto 1060
            bin(5)=bin(5)+1
            goto 1010
 1060       if(ir.ge.165) goto 1010
            bin(6)=bin(6)+1
 1010     continue
c
          write(6,*) 'ori. ratio of cells'
          write(6,*) '  0  ', real(bin(1))/real(nvc0)
          write(6,*) ' 30  ', real(bin(2))/real(nvc0)
          write(6,*) ' 60  ', real(bin(3))/real(nvc0)
          write(6,*) ' 90  ', real(bin(4))/real(nvc0)
          write(6,*) '120  ', real(bin(5))/real(nvc0)
          write(6,*) '150  ', real(bin(6))/real(nvc0)
c
c-----------------------------------------------------------------------
c I/O
c-----------------------------------------------------------------------
c
      do 131 iys0=1,nvcy0
      do 131 ixs0=1,nvcx0
      write(59,*)rdir(ixs0,iys0),0,0,0,0,
     &           rori(ixs0,iys0),rmg(ixs0,iys0),0
 131  continue
      close(59)
c
c-----------------------------------------------------------------------
c
       stop
       end
