      program osi_dst
c   
c*********************************************************************
c OSI distribution
c input   OSI data
c output  OSI distribution on screen
c
c                                                     NTC: M.Miyashita
c*********************************************************************
      implicit real(a-h,k-l,o-s,u-z)
      implicit real*8(t)
c
      real                qi,qe,dx,beta,qir,ler,qe1,qi1,qer,
     &                    qe0,le0,qec,qic,lec,lic
      common /pa1/        qi,qe,dx,beta,qir,ler,
     &                    qe0,le0,qec,qic,lec,lic
      common /pa2/        strd,stuo,stoi,sttl
      integer             nvc,nvcx,nvcy,nrf,nrfx,nrfy
      common /pa3/        nvc,nvcx,nvcy,nrf,nrfx,nrfy
      real*8              timl
      common /time1/      timl
c
      integer             nvcx0,nvcy0,ncff,
     &                    nrfx0,nrfy0,nrf0,nvc0
      integer             mgstep,iwr
      parameter          (nvcx0=48,nvcy0=nvcx0,nvc0=nvcx0*nvcy0,
     &                    nrfx0=24,nrfy0=nrfx0,nrf0=nrfx0*nrfy0)
c
c OSI bin
      parameter          (nbin=12, rstp=0.08) 
      real                dosi(nvcx0,nvcy0)
      integer             osibin(nbin)
c
      character           inpf*30,inpf2*30,inpf3*30,inpf4*30,inpf5*30,
     &                    outf*30,outf2*30,outf3*30
      real*4              st1
c--------------------------------------------------------------------
c
c ### Input Files ###
      write(6,*)'# input OSI data (dosi#) =?'
      read(5,*) inpf2
      open(57,file=inpf2,status='old')
c
c*********************************************************************
c     input response file (for direction)
c*********************************************************************
      do 100 iys0=1,nvcy0
      do 100 ixs0=1,nvcx0
      read(57,*) dosi(ixs0,iys0)
 100  continue
      close(57)
c
c clear bin
      do 110 i=1,nbin
        osibin(i)=0
 110  continue
c
      z2=0.0
      z3=0.0
      z4=0.0
      do 120 iys0=1,nvcy0
      do 120 ixs0=1,nvcx0
         z0=dosi(ixs0,iys0)
         z2=z2+z0
         ir=int(z0/rstp)+1
         osibin(ir)=osibin(ir)+1
         if(z0.eq.0) goto 120
         z3=z3+z0
         z4=z4+1
 120  continue
c
      ir=0
      z0=0
      do 130 i=1,nbin
         if(osibin(i).lt.z0) goto 130
         ir=i
         z0=osibin(i)
 130  continue
c
      z1=0.0
      do 140 i=1,nbin
         z1=z1+osibin(i)*real(i-1)*rstp
 140  continue
c
      write(6,*) '                                '
      write(6,*) '   ** OSI distribution **'
      write(6,*) '    OSI         ratio of cells'
      write(6,*) '--------------------------------'
      do 150 i=1,nbin
         write(6,*) real(i-1)*rstp,  real(osibin(i))/real(nvc0)
 150  continue
      write(6,*) '--------------------------------'
c
c      write(6,*) z2/nvc0
c      write(6,*) 'max cell number', z0
c      write(6,*) 'max OSI in bin ', real(ir-1)*rstp
c      write(6,*) 'OSI average    ', z2/nvc0
c      write(6,*) 'OSI average without OSI=0   ', z3/z4
c
c-----------------------------------------------------------------------
c
       stop
       end
