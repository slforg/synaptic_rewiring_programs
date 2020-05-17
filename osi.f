      program osi_idx
c-----------------------------------------------------------------
c     OSI index                                          2020.4.4
c-----------------------------------------------------------------
c

      parameter          (nvcx0=48,nvcy0=nvcx0,nvc0=nvcx0*nvcy0)
      parameter          (idrc=360, iori=180, dr=0.25)
c
      real                drsp(nvcx0,nvcy0,idrc),
     &                    dori0(nvcx0,nvcy0,iori),
     &                    dori(nvcx0,nvcy0,iori*3),
     &                    zmax,zmin,z0,zw,zs,ze,
     &                    rosi(nvcx0,nvcy0)
c
      character           inpf*30,inpf2*30,inpf3*30,inpf4*30,inpf5*30,
     &                    outf*30,outf2*30,outf3*30
c
c-------------------------------------------------
c
      write(6,*)'# Input response file (drsp#) =?'
      read(5,*) inpf
      open(57,file=inpf,status='old')
c
      write(6,*)'# Output OSI file     (dosi#) =?'
      read(5,*) outf
      open(59,file=outf,status='unknown')
c
c input data
      do 100 iys0=1,nvcy0
      do 100 ixs0=1,nvcx0
      do 100 i=1,idrc
      read(57,*) drsp(ixs0,iys0,i)
 100  continue
      close(57)
c
c direction -> orientation
      do 200 iys0=1,nvcy0
      do 200 ixs0=1,nvcx0
      do 200 i=1, iori/2
         dori0(ixs0,iys0,i)=drsp(ixs0,iys0,i+90)+drsp(ixs0,iys0,i+270)
 200  continue
c
      do 210 iys0=1,nvcy0
      do 210 ixs0=1,nvcx0
      do 210 i=iori/2+1, iori
         dori0(ixs0,iys0,i)=drsp(ixs0,iys0,i-90)+drsp(ixs0,iys0,i+90)   
 210  continue
c
c OSI
      do 300 iys0=1,nvcy0
      do 300 ixs0=1,nvcx0
      do 300 i=1,iori
         dori(ixs0,iys0,i)=dori0(ixs0,iys0,i)
         dori(ixs0,iys0,i+iori)=dori0(ixs0,iys0,i)
         dori(ixs0,iys0,i+iori+iori)=dori0(ixs0,iys0,i)
 300  continue
c
      do 400 iys0=1,nvcy0
      do 400 ixs0=1,nvcx0
         zmax=0.0
         zmin=1000.0
         ir0=0
      do 410 i=iori+1,2*iori
         z0=dori(ixs0,iys0,i)
         if(zmax.gt.z0) goto 420
         zmax=z0
         ir0=i
 420     if(zmin.lt.z0) goto 410
         zmin=z0
 410  continue
c
c half width
      zh=zmax/2
         icnt=0
         icnts=0
         zs=0.0
         ifs=0
      do 430 i=1,iori/2
         i0=ir0-i
         z0=dori(ixs0,iys0,i0)
         if(z0.gt.zh) goto 435
         if(ifs.eq.1) goto 430
         ifs=1
         zs=real(i)
         goto 430
 435     if(ifs.eq.0) goto 430
         icnts=icnts+1
 430  continue
c
      z1=0.0
         icnte=0
         ze=0.0
         ife=0
      do 440 i=1,iori/2
         i0=ir0+i
         z0=dori(ixs0,iys0,i0)
         if(z0.gt.zh) goto 445
         if(ife.eq.1) goto 440
         ife=1
         ze=real(i)
         goto 440
 445     if(ife.eq.0) goto 440
         icnte=icnte+1
 440  continue
c
      icnt=icnts+icnte
c
c OSI
      rosi(ixs0,iys0)=0.0
      if(icnt.gt.2) goto 400
      zw=zs+ze
      rosi(ixs0,iys0)=(1.0-zmin/zmax)*(1.0-zw/180.0)
c
 400  continue     
c   
c output data
      do 500 iys0=1,nvcy0
      do 500 ixs0=1,nvcx0
      write(59,*) rosi(ixs0,iys0)
 500  continue
      close(59)
c
      stop
      end

