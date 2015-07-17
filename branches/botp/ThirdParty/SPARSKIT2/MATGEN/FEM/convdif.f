      program convdif 
c-----------------------------------------------------------------------
c this driver will generate a finite element matrix for the
c convection-diffusion problem 
c
c                  - Div ( K(x,y) Grad u ) + C grad u = f
c                    u = 0 on boundary 
c 
c (Dirichlet boundary conditions). 
c----------------------------------------------------------------------- 
c this code will prompt for desired mesh (from 0 to 9, with one being
c a user input one) and will general the matrix in harewell-Boeing format
c in the file mat.hb. It also generates two post script files, one 
c showing the pattern of the matrix (mat.ps) and the other showing the
c corresponding mesh (msh.ps).
c-----------------------------------------------------------------------
c the structure and organization of the fem codes follow very closely 
c that of the book by Noborou Kikuchi (Finite element methods in 
c mechanics,   Cambridge Univ. press, 1986). 
c-----------------------------------------------------------------------
c coded Y. Saad and S. Ma -- this version dated August 11, 1993 
c----------------------------------------------------------------------- 
      implicit none 
      integer maxnx, maxnel
      parameter (maxnx = 8000, maxnel = 15000)
      real*8 a(7*maxnx),x(maxnx),y(maxnx),f(3*maxnx)
      integer ijk(3,maxnel),ichild(12,maxnel),iparnts(2,maxnx),
     *     ia(maxnx),ja(7*maxnx),iwk(maxnx),jwk(maxnx),nodcode(maxnx),
     *     iperm(maxnx)
c 
      character matfile*20, title*72, munt*2,key*8, type*3
      real size 
c     
      integer iin,node,nx,nelx,iout,ndeg,na,nmesh,nref,nxmax,nelmax,nb,
     *    ii,nxnew, nelxnew,ierr,job,n,ncol,mode,ptitle,ifmt 
      external xyk, funb, func, fung
      data iin/7/,node/3/,nx/0/,nelx/0/,iout/8/,ndeg/12/, na/3000/
c--------------------------------------------------------------
c choose starting mesh   --- 
c--------------------------------------------------------------
c     files for output 
c 
      open(unit=10,file='mat.hb') 
      open(unit=11,file='msh.ps') 
      open(unit=12,file='mat.ps') 
c-----------------------------------------------------------------------
      print *, ' enter chosen mesh '
      read (*,*) nmesh
      if (nmesh .eq. 0) then 
         print *, 'enter input file for initial mesh '
         read(*,'(a20)') matfile
         open (unit=7,file=matfile)
      endif
c----------------------------------------------------------------------- 
      print *, ' Enter the number of refinements desired '
      read (*,*) nref
      call inmesh (nmesh,iin,nx,nelx,node,x,y,nodcode,ijk,iperm) 
c
c     ...REFINE THE GRID
c 
      nxmax = maxnx
      nelmax= maxnel
      nb = 0 
      do 10 ii = 1, nref
c     
c     estimate the number nx and nelx at next refinement level.
c
         call checkref(nx,nelx,ijk,node,nodcode,nb,nxnew,nelxnew)
         if (nxnew .gt. nxmax .or. nelxnew .gt. nelmax) then
	    print *, ' Was able to do only ', ii-1 ,'  refinements'
	    goto 11
         endif
c
c     ...if OK refine  all elements 
c     
         call refall(nx,nelx,ijk,node,ndeg,x,y,ichild,iparnts,nodcode,
     *        nxmax, nelmax, ierr)
         if (ierr .ne. 0) print *, '** ERROR IN REFALL : ierr =',ierr 
 10   continue
 11   continue
c-----------------------------------------------------------------------
      job = 0
c-----------------------------------------------------------------------
c     assemble the matrix in CSR format 
c----------------------------------------------------------------------- 
      call assmbo (nx,nelx,node,ijk,nodcode,x,y,
     *     a,ja,ia,f,iwk,jwk,ierr,xyk,funb,func,fung) 
      n = nx 
c----------------------Harewell-Boeing matrix---------------------------
c---------------------1---------2---------3---------5---------6
c            12345678901234567890123456789012345678901234567890
      title='1Sample matrix from SPARSKIT  ' 
      key = 'SPARSKIT'
      type = 'rua'
      ifmt = 6
      job = 2 
      call prtmt (n,n,a,ja,ia,f,'NN',title,key,type,ifmt,job,10)       
c----------------------Plot of mesh-------------------------------------
c----------------------------------------------------------------------- 
      size = 6.0 
      munt = 'in' 
      mode = 0 
      title ='Finite element mesh ' 
      ptitle = 1
      call psgrid (nx,ja,ia,x,y,title,ptitle,size,munt,11) 
c      hsize  = 5.6 
c      vsize = 5.6 
c      xleft = 0.0
c      bot = 0.0
c      job = 30 
c      call texgrd(nx,ja,ia,x,y,munt,size,vsize,hsize,
c     *     xleft,bot,job,title,ptitle,ijk,node,nelx,11) 
c----------------------Plot of matrix-pattern---------------------------
c----------------------------------------------------------------------- 
      size = 5.5
      mode = 0
      title = 'Assembled Matrix'
      ptitle = 1
      ncol = 0
      iout = 12
      call pspltm(n,n,mode,ja,ia,title,ptitle,size,munt,ncol,iwk,12) 
c      xleft = 0.00 
c      bot = 0.70 
c      job =  0
c      call texplt(nx,nx,mode,ja,ia,munt,size,vsize,hsize,xleft,bot,
c     *     job,title,ptitle,ncol,iwk,12) 
      stop
c-----end-of-program-convdif--------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------
