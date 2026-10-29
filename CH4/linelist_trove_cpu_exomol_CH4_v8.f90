program linelist_trove
 implicit none

 integer(4),parameter ::  jmax_ = 100,ipolyad_max_= 50
 real(8),parameter   :: error_max = 2.0, unc_max = 0.2
 integer(4) j,ilevel,ilevelf_,ileveli_,info,i,iparity,icontr,iterm,nterms
 integer(4) nfiles, jmax, ilines,Ncontr,igamma,igammai_,id_t,ID_m,jmax_inp,ipolyad_max
 integer(4) nlevels,v_(30),sym_t,J_t,imis,v_norm(30),npolyad(0:100),tau_,ipolyad,isym,v0,n0
 integer(8) :: idmax,id,idf,idi,i_t,idmin,nCD
 integer,allocatable :: sym(:),Jktau(:,:),NN(:,:,:),numjgp(:,:,:),numj(:,:)

 real(8),allocatable :: energies(:),coeff_(:)
 integer(4),allocatable :: g(:)
 integer(8)  :: nswap,iswap,nswap_,iselect,nselect,kswaps,jswap
 real(8),allocatable :: acoef(:),tranfreq(:)
 integer(4),allocatable :: ilevelf(:),ileveli(:)
 integer(4),allocatable :: ilevel_(:),v_l(:,:),v_norm_(:,:),igammav_(:),igammar_(:)
 integer(8),allocatable :: id_(:),Nlines(:)
 character(4),allocatable :: gamma_(:),gammar_(:),gammav_(:)
 integer(4),allocatable :: icontr_local(:,:),icontr_normal(:,:),lvib(:),icontr_vib(:)
 integer(4),allocatable :: QN_MARVEL(:,:),ncount_marvel(:,:,:)
 !real(8),allocatable :: Energy_MARVEL(:),unc_marvel(:)

 real(8) :: acoef_,abscoef,energy,energyf,energyi,tranfreq_,linestr,ZPE,coeff,nu1,nu2,unc,w(0:100)
 integer ,allocatable     :: gns(:)
 character(len=4),allocatable   :: symm(:)
  
 logical :: energeyfile_do = .false. , ifopen, ifexist
 logical :: MARVEL_do = .false.,marvel_yes = .false.
 logical :: complete = .false.
 character(4) gamma,gammav,gammar
 character(100) intfilename(4000),enrfilename,enrfilenameout,intfilenameout,dscrfilename,filename
 character(100) enrsubsfilename
 character(500) fname
 character(2)  :: ch_t2,sym2
 character(1)  :: ef
 character(4)  :: fchar
 character(7)  :: fchar1,fchar2
 character(400) :: form
 integer(4)    :: unit_f,nmodes,nsym
 integer(8)    :: Nlevels_marvel,ilevel8

 real(8),allocatable :: Energy_MARVEL(:),unc_marvel(:),tau(:),coeff_contr(:)
 integer(4),allocatable :: ID_marvel(:)
 character(2),allocatable  :: gammav_marvel(:)
 character(2)  :: flag

 real(8),parameter :: planck=6.6260693d-27,avogno=6.0221415d+23,vellgt=2.99792458d+10,boltz=1.380658d-16

 type numT
      !
      integer   :: nmax
      integer,pointer   :: nn(:)
      integer,pointer   :: sym(:)
      integer,pointer   :: J(:)
      real(8),pointer   :: energy(:)
      real(8),pointer   :: unc(:)
      !
 end type numT

 type(numT),allocatable :: level(:,:)
 type(numT),allocatable :: MARVEL(:,:,:)

    !read number of intensity files
    read*,nfiles
    !
    if (nfiles>4000) then 
      print('("Too many files (>400):",i0)'), nfiles
      stop 'Too many files'
    endif
    !
    !read intensities filenames
    do i = 1,nfiles
     read*,intfilename(i)
    enddo
    !
    !read energies filename
    read*,enrfilename
    !
    !read descr-filename - in
    read*,dscrfilename
    !
    !read intensity filename - out
    read*,intfilenameout
    !
    !read energies filename - out 
    read*,enrfilenameout
    !
    !read file with energies to substitute 
    read*,enrsubsfilename

    !read nmodes, Nsym
    read*,nmodes,Nsym


    allocate(symm(nsym),gns(nsym),stat=info); 
    if (info/=0) stop 'error: symm,gns  out of memory'

    !read nuclear spin statistics weights
    read*,(symm(isym),gns(isym),isym=1,nsym)

    !read harmonic freq errors 
    read*,w(0:nmodes)

    !read the frequency window
    read*,nu1,nu2
    !
    print*,'range =  ',nu1,' - ',nu2

    !read the number of lines to keep in memory
    read*,nswap
    !
    if (nswap<0) nswap = huge(nswap)/100
    !
    print*,'nswap = ',nswap
    !
    print("('Read in the contr-descr-files and re-assign the numbering of the levels ..')")
    filename = trim(dscrfilename)
    !
    print("( 'pre-process ',a)"), filename
    !
    open(unit=14,file=trim(filename),action='read',status='old')
    !
    do i = 1,13+nmodes
      read(14,"(a400)") form
    enddo
    !
    read(14,*) Ncontr
    !
    allocate(icontr_local(Ncontr,0:nmodes),icontr_normal(Ncontr,0:nmodes),coeff_contr(Ncontr),stat=info); 
    if (info/=0) stop 'error: icontr is out of memory'
    !
    i = 0
    ipolyad_max = 0
    !
    do 
       !
       read(14,"(a400)",end=14) form
       !
       if (form(1:3)=='End') exit
       !
       read(form,*) ilevel,igamma,ileveli_,igammai_,energy,v0,v_(1:nmodes),n0,v_norm(1:nmodes),coeff_contr(ilevel)
       !
       if (ilevel==1) then
         ZPE = energy
         print*,"ZPE=",ZPE
       endif
       !
       i = i + 1
       !
       icontr_local (ileveli_,0) = v0
       icontr_local (ileveli_,1:nmodes) = v_(1:nmodes)
       icontr_normal(ileveli_,0) = n0
       icontr_normal(ileveli_,1:nmodes) = v_norm(1:nmodes)
       !
       ipolyad = 2.0d0*sum( v_(1:4) ) + sum(v_(5:9))
       ipolyad_max = max(ipolyad,ipolyad_max)
       !
       cycle
    14 exit
    end do
    !
    close(14)
    !
    if ( all(trim(enrsubsfilename)/=(/'none','NONE','null','NULL'/)) ) then  
      !
      marvel_do = .true.
      !
      print("('Read MARVEL/HITRAN again... ')")
      filename = trim(enrsubsfilename)
      !
      print("( 'process ',a)"), filename
      !
      open(unit=24,file=trim(filename),action='read',status='old')
      !
      !
    endif
    !
    if (marvel_do) then 
       !
       !
       print("('Read in the file containing MARVEL/HITRAN energies and QNs to substitute TROVE ... ')")
       filename = trim(enrsubsfilename)
       !
       print("( 'process ',a)"), filename
       !
       open(unit=24,file=trim(filename),action='read',status='old')
       !
       i = 0
       Nlevels_MARVEL = 0
       !
       do 
          !
          read(24,*,end=24) ilevel8
          !
          i = i + 1
          Nlevels_MARVEL = max(Nlevels_MARVEL,ilevel8)
          !
          cycle
       24 exit
       end do
       !
       nlevels = i
       !
       rewind(24)
       !
       allocate(Energy_MARVEL(nlevels),ID_marvel(Nlevels_MARVEL),unc_marvel(nlevels),stat=info); 
       if (info/=0) then 
         print*,'error: MARVEL is out of memory'
         stop 'error: MARVEL is out of memory'
       endif
  
       allocate(QN_MARVEL(nlevels,8),gammav_marvel(nlevels),stat=info); 
       !
       !QN_MARVEL(nlevels,1:5),
       if (info/=0) then 
         print*,'error: QN MARVEL is out of memory'
         stop 'error: QN MARVEL is out of memory'
       endif
       !
       unc_marvel = -1
       energy_Marvel = 0
       QN_MARVEL = -1
       ID_marvel = 0
       gammav_marvel = 'XX'
       !
       i = 1
       !
       do 
          !
          read(24,*,end=25) ilevel8,energy,unc !,QN_MARVEL(i,1:6),gammav_marvel(i),QN_MARVEL(i,7:8)
          !
          ID_marvel(ilevel8) = i
          !
          energy_Marvel(i) = energy
          unc_marvel(i) = abs(unc)
          !
          i = i + 1
          !
          cycle
       25 exit
       end do
       !
       close(24)
       !
    endif
    !
    !read the energy file 
    open(unit=11,file=trim(enrfilename))
    !
    i = 0
    jmax = 0
    !
    print("('Counting states...')")
    !
    idmax = 0
    idmin = 0
    !
    do
       !
       !read(11,"(2x,a4,i7,f14.6,9x,i3)",end=11) gamma,ilevel,energy,j
       !
       read(11,*,end=11) id,energy,i_t,j,gamma,v_(1:9) ! ,gammav_,j_t,tau_,coeff,ch_t2
       !
       ! 957020  1562.509565     75       7  A1     0  0  0  0  0  0  0  0  1  F2    2   0 F2    0.88 ::         1   0  0  0  0  0  0  1  1  1

       !
       i = i + 1
       !
       jmax = max(jmax,j)
       !
       idmax = max(id,idmax)
       idmin = min(id,idmin)
       !
       cycle
11     exit
    end do
    !
    print("(' jmax = ',i4)"),jmax
    !
    nlevels = i
    print("(' nlevels = ',i0)"), nlevels
    !
    rewind(11)
    !
    allocate(numjgp(0:Jmax,nsym,0:100),stat=info)
    if (info/=0) stop 'error: numjgp is out of memory'
    !
    allocate(sym(nlevels),energies(nlevels),Jktau(nlevels,3),g(nlevels),icontr_vib(nlevels),stat=info); 
    if (info/=0) stop 'error: sym,energies,Jktau,g,icontr is out of memory'
    !
    allocate(ilevel_(idmin:idmax),stat=info); 
    if (info/=0) stop 'error: ilevel is  of memory'
    !
    if ( all(trim(enrfilenameout)/=(/'none','NONE','null','NULL'/)) ) then  
      !
      energeyfile_do = .true.
      !
      open(unit=12,file=trim(enrfilenameout),action='write',status='replace')
      !
    endif
    !
    sym = 0
    !
    sym_t = 1
    J_t =  0
    !
    print("('Reading energies...')")
    !
    npolyad = 0
    ilevel = 0 
    !
    allocate(gamma_(nlevels),gammar_(nlevels),gammav_(nlevels),v_l(nmodes,nlevels),coeff_(nlevels),v_norm_(nmodes,nlevels),&
             igammav_(nlevels),igammar_(nlevels),lvib(nlevels),stat=info)
    if (info/=0) stop 'error: swap-objects  out of memory'
    
    allocate(id_(nlevels),stat=info)
    if (info/=0) stop 'error: id_  out of memory'
    !
    numjgp = 0 
    !
    if (jmax>400) stop 'Jmax>400'
    !
!          1     0.000000      1       0  A1g    0  0  0  0  0  0  0  A1g   0   0 A1g   1.00 ::         1   0  0  0  0  0  0  0        1
!          2  1230.381731      1       0  A1g    0  0  0  1  0  1  0  A1g   0   0 A1g   1.00 ::         2   0  0  0  1  0  1  0        4
!          3  1449.114864      1       0  A1g    0  0  0  0  0  0  2  A1g   0   0 A1g   1.00 ::         3   0  0  0  0  0  0  2        9

    do i = 1,nlevels
       !
       read(11,*) id,energies(i),g(i),J,gamma,v_(1:nmodes),gammav,Jktau(i,2:3),gammar,coeff,ch_t2,iterm,v_norm(1:nmodes),icontr
       !
       !read(ch_t200,"(2x,a4,i7,f14.6,4x,a4,1x,3i3,4x,a4,1x,<nmodes>i4,7x,f5.2,2x,<nmodes>i4)") 
       !gamma,ilevel,energies(i),gammar,Jktau(i,1:3),gammav,v_(1:nmodes),coeff,v_norm(1:nmodes)
       !
       id_(i)          = id
       gamma_(i)       = gamma
       ilevel_(id)     = i
       gammar_(i)      = gammar
       gammav_(i)      = gammav
       v_l(1:nmodes,i) = v_(1:nmodes)
       coeff_(i)       = coeff
       icontr_vib(i) = icontr
       !
       Jktau(i,1) = J
       !
       do isym = 1,nsym
         if (trim( adjustl(gamma)) == symm(isym)) then 
           igamma = isym
           exit
         endif
       enddo
       !
       sym(i) = igamma
       !
       ipolyad = 2.0d0*sum( v_(1:4) ) + sum(v_(5:9))
       !
       numjgp(j,igamma,ipolyad) = numjgp(j,igamma,ipolyad) + 1
       !
       lvib(i) = numjgp(j,igamma,ipolyad)
       !
       v_norm_(1:nmodes,i) = icontr_normal(icontr,1:nmodes)
       !
       if (J==0) then
         !
         coeff_(i) = coeff_contr(icontr)**2
         !
       endif
       !
       !v_norm_(1:nmodes,i) = v_norm(1:nmodes)
       !
    enddo
    !
    print("('Printing states file...')")
    !
    if (energeyfile_do) then
       !
       do i = 1,nlevels
         !
         energy = energies(i)
         !
         v_(1:nmodes) = -1
         J_t = Jktau(i,1)
         sym_t = sym(i)
         ID_t = id_(i)
         !
         v_(1:nmodes) = v_l(1:nmodes,i)
         !
         ipolyad = 2.0d0*sum( v_(1:4) ) + sum(v_(5:9))
         !
         unc = sum(w(1:nmodes)*v_norm_(1:nmodes,i))+w(0)*J_t
         iterm = lvib(i)
         !
         flag = "Ca"
         !
         energy = energies(i)
         !
         if (id_t<=Nlevels_MARVEL) then
            if (ID_marvel(id_t)/=0) then
               !
               ID_m = ID_marvel(id_t)
               !
               energy = Energy_MARVEL(ID_m)
               !
               unc = unc_marvel(ID_m)
               !v_norm_(1:nmodes,i) = QN_MARVEL(id_m,1:nmodes)
               !jktau(i,2:3) = QN_MARVEL(id_m,7:8)
               !gammav_(i) = gammav_MARVEL(id_m)
               !
               flag = "Ma"
               !
            endif
         endif
         !
         unc = max(unc,0.000001)
         !
         !if (energies(i)<=99999.99) then 
         !
         write(12,"(i12,1x,f12.6,1x,i6,1x,i7,1x,1x,f12.6,2x,a3,1x,i4,1x,i8,2x,<nmodes>i4,1x,a3,1x,1i4,1x,a3,1x,i7,2x,f5.2,3x,<nmodes>i4,1x,a2,2x,f12.6)") & 
                    i,energy,g(i),jktau(i,1),unc,gamma_(i),ipolyad,iterm,v_norm_(1:nmodes,i),gammav_(i),jktau(i,2),gammar_(i),icontr_vib(i),&
                    coeff_(i),v_l(1:nmodes,i),flag,energies(i)


                       !i,energies(i),g(i),jktau(i,1),gamma_(i),v_l(1:nmodes,i),gammav_(i),jktau(i,1:3),gammar_(i),coeff_(i),v_norm_(1:nmodes,i),,energies(i)
                       !
         !if (v_(1)==-1) then ! MARVEL 
         !  ! CO2
         !  write(12,"(i12,1x,f12.6,1x,i6,1x,i7,1x,f12.2,1x,a3,1x,a1,1x,<nmodes>i4,i4,2x,f5.2,3x,5i4,3x,<nmodes>i4,1x,a3,i4,1x,a3)") & 
         !            i,energy,g(i),jktau(i,1),unc,gamma_(i),ef,v_norm_(1:2,i),lvib(i),v_norm_(3,i),coeff_(i),v_(1:5),v_l(1:nmodes,i)! ,gammav_(i),jktau(i,2),gammar_(i)
         !else
         !   write(12,"(i12,1x,f12.6,1x,i6,1x,i7,1x,f12.4,1x,a3,1x,a1,1x,<nmodes>i4,i4,2x,f5.2,3x,5i4,3x,<nmodes>i4,1x,a3,i4,1x,a3)") & 
         !            i,energy,g(i),jktau(i,1),unc,gamma_(i),ef,v_norm_(1:2,i),lvib(i),v_norm_(3,i),coeff_(i),v_(1:5),v_l(1:nmodes,i)! ,gammav_(i),jktau(i,2),gammar_(i)
         !
         !endif
         !
         !else
         !   write(12,"(i12,1x,e12.6,1x,i6,1x,i7,1x,i4,2x,<nmodes>i4,1x,i4,3x,4i4,2x,i8,2x,f5.2,3x,<nmodes>i4)") & 
         !              !i,energies(i),gns(sym(i))*(2*jktau(i,1)+1),jktau(i,1),sym(i),v_norm(1:4),igammav,jktau(i,1:3),igammar,ilevel,coeff,v_(1:4)
         !              !
         !              i,energies(i),gns(sym(i))*(2*jktau(i,1)+1),jktau(i,1),sym(i),v_norm_(1:nmodes,i),igammav_(i),jktau(i,1:3),igammar_(i),0,coeff_(i),v_l(1:nmodes,i)
         !       
         !endif
         ! 
       end do
       !
    endif
    !
    deallocate(gamma_,gammar_,gammav_,v_l,coeff_,v_norm_,igammav_,igammar_,lvib,icontr_normal,icontr_local,stat=info)
    if (info/=0) stop 'error: swap-objects  out of memory'
    !
    if (marvel_do) then 
       deallocate(MARVEL,stat=info)

       deallocate(Energy_MARVEL,ID_marvel,unc_marvel,stat=info); 
       if (info/=0) then 
         print*,'error: Energy_MARVEL is out of memory'
         stop 'error: Energy_MARVEL is out of memory'
       endif
  
       deallocate(QN_MARVEL,gammav_marvel,stat=info); 
       !
       if (info/=0) then 
         print*,'error: QN_MARVEL is out of memory'
         stop 'error: QN_MARVEL is out of memory'
       endif
    endif
    !
    print*,"jmax= ",jmax
    !
    close(11)
    !
    if (energeyfile_do) close(12)
    !
    if (nu1<-0.00001) then
      print*,"Negative frequency range, finish"
      stop "Negative frequency range, finish"
    endif
    !
    !deallocate(nvib,lvib,tauvib,vvib,symvib)
    !
    print("(' N of levels =  ',i0)"), nlevels
    !
    ! finish if the intensity output is not requested
    !
    if (trim(intfilenameout)=="none".or.trim(intfilenameout)=="NONE") stop 
    !
    !
    print("('Generate the Transition file ...')")
    !
    !
    !start loop over all transitions
    !
    allocate(Nlines(0:jmax))
    !
    allocate(tranfreq(nswap),ilevelf(nswap),ileveli(nswap),acoef(nswap),stat=info)
    if (info/=0) stop 'error: swap-objects  out of memory'
    !
    ilines = 0
    Nlines = 0
    !
    write(fchar1, '(i7)') nint(nu1)
    write(fchar2, '(i7)') nint(nu2)
    !write(fchar, '(i4)') old_swaps
    !
    fname = trim(intfilenameout)//'_'//trim(adjustl(fchar1))//'-'//trim(adjustl(fchar2))//'_'//trim(adjustl(fchar))//'.tra'
    !
    INQUIRE (FILE = fname, opened = ifopen)
    INQUIRE (FILE = fname, exist = ifexist)
    !
    unit_f = 20
    !
    !if (.not.ifopen)  open(unit=unit_f,file=fname,action='write',buffered='yes',status='replace')
    !
    !open(unit=unit_f,file=fname,action='write',buffered='yes',status='replace')
    !
    do i = 1,nfiles
      !
      open(unit=1,file=trim(intfilename(i)))
      !
      print("( 'process ',a,3x,i0)"), intfilename(i),ilines
      !
      fname = trim(intfilename(i))//'_'//trim(adjustl(fchar1))//'-'//trim(adjustl(fchar2))//'.tra'
      !
      unit_f = 20
      !
      open(unit=unit_f,file=fname,action='write',buffered='yes',status='replace')
      !
      ilines = 0
      imis = 0
      kswaps = 0 
      !
      loop_file : do
        !
        nswap_ = nswap
        iselect = 0
        nselect = 0
        complete  = .false.
        !
        kswaps = kswaps + 1
        !
        print*,"Reading transitions ...  ",kswaps,nswap
        !
        do iswap = 1,nswap
          !
          !read(1,*,end=20) ilevelf(iswap),ileveli(iswap),acoef(iswap)
     
          !   read new line from intensities file
          !
          iselect = iselect + 1
          !
          !read(1,*,end=20) tranfreq(iselect),ilevelf(iselect),jf(iselect),igammaf(iselect),ch_t2,ileveli(iselect),ji(iselect),igammai(iselect),acoef(iselect)
          !
          read(1,*,end=20) idf,idi,acoef(iselect)
          !
          ilevelf(iselect) = ilevel_(idf)
          ileveli(iselect) = ilevel_(idi)
          tranfreq(iselect) = energies(ilevelf(iselect))-energies(ileveli(iselect))
          !
          !j = max(Jktau(ilevelf(iselect),1),Jktau(ileveli(iselect),1))
          !
          j = Jktau(ileveli(iselect),1)
          !
          if ( tranfreq(iselect)<nu1-1e-6 .or. tranfreq(iselect)>=(nu2-1e-6) .or. j>jmax) iselect = iselect - 1
          !
          Nlines(j)= Nlines(j) + 1
          !
          cycle 
          !
        20 continue 
          !
          nswap_ = iswap-1
          iselect = iselect - 1
          !
          !complete  = .true.
          !
          exit
          !
        enddo
        !
        !print('((i10,2i10,e17.9,f15.6))'),nswap_,ilevelf(iselect),ileveli(iselect),acoef(iselect),tranfreq(iselect)
        !
        print('(i10,a)'),iselect," transitions"
        !
        nselect = iselect
        !
        if (nswap_<1.and.kswaps==1) then 
           write(6,"('No transitions for this range ')") 
           !stop 'No transitions for this range, nswap_=0'
        endif
        !
        if (nswap_<nswap) print('(a,i14,a)'),"last swap with ... ",nswap_," lines"
        !
        !
        if (complete) then 
          !
          ! sort 
          !
          print('(a)'),"Sort ...."
          !
          do iswap = 1,nselect
            !
            tranfreq_ = tranfreq(iswap)
            !
            do jswap = iswap+1,nselect
              !
              if ( tranfreq_>tranfreq(jswap) ) then
                !
                ileveli_ = ileveli(jswap)
                ilevelf_ = ilevelf(jswap)
                acoef_ = acoef(jswap) 
                tranfreq_ = tranfreq(jswap) 
                !
                ileveli(jswap) = ileveli(iswap)
                ilevelf(jswap) = ilevelf(iswap)
                acoef(jswap) = acoef(iswap) 
                tranfreq(jswap) = tranfreq(iswap)
                !
                ileveli(iswap) = ileveli_
                ilevelf(iswap) = ilevelf_
                acoef(iswap) = acoef_
                tranfreq(iswap) = tranfreq_
                !
              end if
              !
            end do
            !
          enddo
          !
          print('(a)'),"...done!"
          !
          print('(a)'),"Check duplicates ..."
          !
          do iswap = 2,nselect
            !
            ileveli_ = ileveli(iswap)
            ilevelf_ = ilevelf(iswap)
            acoef_ = acoef(iswap) 
            acoef_ = tranfreq(iswap) 
            !
            if ( ileveli_== ileveli(iswap-1) .and. ilevelf_== ilevelf(iswap-1) ) then
                !
                print('(a,2(2i10,e12.4,f15.6))'),"D: ",iswap,ileveli_,ilevelf_,acoef_,acoef_,ileveli(iswap-1), &
                      ilevelf(iswap-1),acoef(iswap),tranfreq(iswap)
                !
            end if
            !
          enddo
          !
          print('(a)'),"...done!"        
          !
        endif
        !
        print*,"Printing to trans-file...",kswaps,nselect
        !
        if (nselect>0) then
          !
          write(unit_f,"(i12,1x,i12,1x,es10.4,1x,f15.6)") (ilevelf(iswap),ileveli(iswap),acoef(iswap),tranfreq(iswap),iswap=1,nselect)
          !
        endif
        ilines = ilines + nselect
        !
        if (nswap_<nswap) exit loop_file
        !
      enddo loop_file
      !
      close(unit_f)
      !
      !if (imis /= 0) print("( 'missassigned idescr = ',i)"), imis
      !
    enddo 
    !
    do i = 0,jmax
     print("(i4,2x,i14)"), i,Nlines(i)
    enddo
    !
    deallocate(sym,energies,Jktau)
    !
    print("( 'Done!')")
    !
    close(1)
    close(13)


end program linelist_trove
