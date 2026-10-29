program linelist_trove
 use, intrinsic :: ieee_arithmetic
 use, intrinsic :: iso_fortran_env, only: real32

 implicit none

 integer(4) j,ilevel,ilevelf_,ileveli_,info,i,iparity,icontr,iterm
 integer(4) nfiles, jmax, jmax_inp, ilines,Ncontr,igamma,igammai_,id_t,ID_m
 integer(4) nlevels,v_(30),sym_t,J_t,imis,v_norm(30),npolyad(0:100)
 integer(8) :: idmax,id,idf,idi,i_t,idmin

 real(8),allocatable :: energies(:),coeff_(:)
 integer(4),allocatable :: sym(:),Jktau(:,:),g(:)
 integer(8)  :: nswap,iswap,nswap_,iselect,nselect,kswaps,jswap
 real(8),allocatable :: acoef(:),tranfreq(:)
 integer(4),allocatable :: ilevelf(:),ileveli(:)
 integer(4),allocatable :: ilevel_(:),v_l(:,:),v_norm_(:,:),igammav_(:),igammar_(:),v_spec(:,:)
 integer(8),allocatable :: id_(:),Nlines(:)
 character(4),allocatable :: gamma_(:),gammar_(:),gammav_(:)
 integer(4),allocatable :: icontr_local(:,:),icontr_normal(:,:),lvib(:),icontr_MARVEL(:)
 integer(4),allocatable :: QN_MARVEL(:,:),ID_marvel(:)
 real(8),allocatable :: Energy_MARVEL(:),unc_marvel(:),coeff_contr(:)
 character(2),allocatable :: label_MARVEL(:)

 logical :: energeyfile_do = .false. , ifopen, ifexist
 logical :: complete = .false.

 real(8) :: acoef_,energy,tranfreq_,coeff,nu1,nu2,zpe,unc,w1,w2,w3,brot,enermaxlow=1000000.0,enermaxup=1000000.0

 character(4) gamma,gammav,gammar
 character(100) intfilename(4000),enrfilename,enrfilenameout,intfilenameout,dscrfilename,filename
 character(100) enrsubsfilename
 character(500) fname
 character(2)  :: ch_t2,flag
 character(1)  :: ef
 character(4)  :: fchar
 character(7)  :: fchar1,fchar2
 character(400) :: form
 integer(4)    :: unit_f,nmodes,nsym
 character(3)  :: iv(5)
 integer(8)    :: Nlevels_marvel,ilevel8
 real(real32) :: nan

 real(8),parameter :: planck=6.6260693d-27,avogno=6.0221415d+23,vellgt=2.99792458d+10,boltz=1.380658d-16

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

    !read Nlevels 
    !read*,Nlevels_MARVEL

    !read harmonic freq errors 
    read*,w1,w2,w3,brot


    !read the frequency window
    read*,nu1,nu2
    !
    print*,'range =  ',nu1,' - ',nu2

    !read the number of lines to keep in memory
    read*,nswap
    !
    ! read enermaxlow,enermaxup
    !
    read*,enermaxlow,enermaxup
    !
    read*,jmax_inp
    !
    if (nswap<0) nswap = huge(nswap)/100
    !
    print*,'nswap = ',nswap
    !
    print("('Read in the contr-descr-files and re-assign the numbering of the levels ..')")
    filename = trim(dscrfilename)
    !
    print("( 'process ',a)"), filename
    !
    open(unit=14,file=trim(filename),action='read',status='old')
    !
    do i = 1,13+nmodes
      read(14,"(a400)") form
    enddo
    !
    read(14,*) Ncontr
    !
    allocate(icontr_local(Ncontr,0:3),icontr_normal(Ncontr,0:3),coeff_contr(Ncontr),stat=info); 
    if (info/=0) stop 'error: icontr is out of memory'
    !
    i = 0
    !
    do 
       !
       read(14,"(a400)",end=14) form
       !
       if (form(1:3)=='End') exit
       !
       read(form,*) ilevel,igamma,ileveli_,igammai_,energy,icontr_local(ilevel,0:3),icontr_normal(ilevel,0:3),coeff_contr(ilevel)
       !
       if (ilevel==1) then
         ZPE = energy
         print*,"ZPE=",ZPE
       endif
       !
       i = i + 1
       !
       cycle
    14 exit
    end do
    !
    close(14)

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
    allocate(Energy_MARVEL(nlevels),QN_MARVEL(nlevels,1:5),ID_marvel(Nlevels_MARVEL),unc_marvel(nlevels),icontr_MARVEL(Ncontr),&
             label_marvel(nlevels),stat=info); 
    if (info/=0) then 
      print*,'error: MARVEL is out of memory'
      stop 'error: MARVEL is out of memory'
    endif
    !
    unc_marvel = -1
    energy_Marvel = 0
    QN_MARVEL = -1
    ID_marvel = 0
    icontr_MARVEL = 0 
    !
    i = 1
    !
    do 
       !
       read(24,*,end=25) ilevel8,energy,v_(1:5),unc,icontr,flag
       !
       ID_marvel(ilevel8) = i
       !
       energy_Marvel(i) = energy
       QN_MARVEL(i,1:5) = v_(1:5)
       unc_marvel(i) = abs(unc)
       label_marvel(i) = flag
       !
       icontr_MARVEL(icontr) = i
       !
       i = i + 1
       !
       cycle
    25 exit
    end do
    !
    close(24)
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
       ! A1        1      0.000000   ( A1 ;  0  0  0) ( A1 ;   0   0   0)      1.00 (   0   0   0   0) (     1)
       !
       read(11,*,end=11) id,energy,i_t,j
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
    allocate(sym(nlevels),energies(nlevels),Jktau(nlevels,3),g(nlevels),stat=info); 
    if (info/=0) stop 'error: sym,v,n,l,Jktau is out of memory'
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
    allocate(gamma_(nlevels),gammar_(nlevels),gammav_(nlevels),v_l(nmodes,nlevels),v_spec(5,nlevels),&
             coeff_(nlevels),v_norm_(nmodes,nlevels),&
             igammav_(nlevels),igammar_(nlevels),lvib(nlevels),stat=info)
    if (info/=0) stop 'error: swap-objects  out of memory'
    
    allocate(id_(nlevels),stat=info)
    if (info/=0) stop 'error: id_  out of memory'
    !
    if (jmax>400) stop 'Jmax>400'
    !
!          1     0.000000      1       0  A1g    0  0  0  0  0  0  0  A1g   0   0 A1g   1.00 ::         1   0  0  0  0  0  0  0        1
!          2  1230.381731      1       0  A1g    0  0  0  1  0  1  0  A1g   0   0 A1g   1.00 ::         2   0  0  0  1  0  1  0        4
!          3  1449.114864      1       0  A1g    0  0  0  0  0  0  2  A1g   0   0 A1g   1.00 ::         3   0  0  0  0  0  0  2        9

    do i = 1,nlevels
       !
       read(11,*) id,energies(i),g(i),Jktau(i,1),gamma,v_norm(1:nmodes),gammav,Jktau(i,2:3),gammar,coeff,ch_t2,iterm,v_(1:nmodes),icontr
       !
       !read(ch_t200,"(2x,a4,i7,f14.6,4x,a4,1x,3i3,4x,a4,1x,<nmodes>i4,7x,f5.2,2x,<nmodes>i4)") 
       !gamma,ilevel,energies(i),gammar,Jktau(i,1:3),gammav,v_(1:nmodes),coeff,v_norm(1:nmodes)
       !
       id_(i)          = id
       gamma_(i)       = gamma
       ilevel_(id)     = i
       gammar_(i)      = gammar
       gammav_(i)      = gammav
       v_l(1:nmodes,i) = icontr_local(icontr,1:nmodes)
       coeff_(i)       = coeff**2
       !
       if (jktau(i,1)==0) coeff_(i) = coeff_contr(icontr)**2
       !
       lvib(i) = icontr_normal(icontr,0)
       v_norm_(1:nmodes,i) = icontr_normal(icontr,1:3)
       !
       v_spec(1,i) = icontr_normal(icontr,1)
       v_spec(2,i) = lvib(i)
       v_spec(3,i) = lvib(i)
       v_spec(4,i) = icontr_normal(icontr,3)
       v_spec(5,i) = 1
       !
       if (icontr_MARVEL(icontr)/=0) then
          ID_m = icontr_MARVEL(icontr)
          v_(1:5) = QN_MARVEL(ID_m,1:5)
          v_spec(1:5,i) = v_(1:5)
       endif
       !
    enddo
    !
    print("('Printing states file...')")
    !
    if (energeyfile_do) then
       !
       do i = 1,nlevels
         !
         select case (trim(gamma_(i)))
           !
         case ("A1","B2")
           !
           iparity = 1
           !
         case ("A2","B1")
           !
           iparity = -1
           !
         end select          
         !
         ef = "f" ; if( (-1)**jktau(i,1)==iparity ) ef = "e"
         !
         flag = "Ca"
         !
         energy = energies(i)
         !
         nan = IEEE_VALUE(nan, IEEE_QUIET_NAN)
         v_(1:5) = nan
         iv(1:5) = "nan"
         !
         v_(1:5) = v_spec(1:5,i)
         !
         ID_t = id_(i)
         !
         unc = w1*v_norm_(1,i)+w2*v_norm_(2,i)+w3*v_norm_(3,i) + jktau(i,1)*(jktau(i,1)+1)*brot
         !
         if (id_t<=Nlevels_MARVEL) then
            if (ID_marvel(id_t)/=0) then
               !
               ID_m = ID_marvel(id_t)
               !
               energy = Energy_MARVEL(ID_m)
               !
               if( abs(energy-energies(i))>1) print*,"error",id_m,energy,energies(i)
               !
               v_(1:5) = QN_MARVEL(ID_m,1:5)
               !
               write(iv(1:5),"(i3)") v_(1:5)  
               !
               unc = max(unc_marvel(ID_m),1e-6)
               !
               flag = label_marvel(id_m)
               !
            endif
         endif
         !
         write(12,"(i12,1x,f12.6,1x,i6,1x,i7,1x,f12.6,1x,a3,1x,a1,1x,<nmodes>i4,i4,2x,f5.2,3x,5i4,3x,<nmodes>i4,1x,a2,1x,f12.6)") & 
                     i,energy,g(i),jktau(i,1),unc,gamma_(i),ef,v_norm_(1:2,i),lvib(i),v_norm_(3,i),coeff_(i),v_(1:5),v_l(1:nmodes,i),flag,energies(i) !gammav_(i),jktau(i,2),gammar_(i)

         !write(12,"(i12,1x,f12.6,1x,i6,1x,i7,1x,f12.6,1x,a3,1x,a1,1x,<nmodes>i4,i4,2x,f5.2,2x,5(1x,a3),3x,<nmodes>i4,1x,a2,1x,f12.6)") & 
         !            i,energy,g(i),jktau(i,1),unc,gamma_(i),ef,v_norm_(1:2,i),lvib(i),v_norm_(3,i),coeff_(i),iv(1:5),v_l(1:nmodes,i),flag,energies(i) !gammav_(i),jktau(i,2),gammar_(i)

         ! 
       end do
       !
    endif
    !
    deallocate(gamma_,gammar_,gammav_,v_l,coeff_,v_norm_,igammav_,igammar_,lvib,icontr_normal,icontr_local,stat=info)
    if (info/=0) stop 'error: swap-objects  out of memory'
    !
    deallocate(QN_marvel,Energy_MARVEL,ID_marvel,stat=info)
    if (info/=0) stop 'error: deall MARVEL'

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
          !print*,idf,idi,ilevelf(iselect),ileveli(iselect),energies(ilevelf(iselect)),energies(ileveli(iselect))
          !
          j = max(Jktau(ilevelf(iselect),1),Jktau(ileveli(iselect),1))
          !
          !print*,"J=",J
          !print*,tranfreq(iselect),nu1,nu2
          !print*,enermaxlow,enermaxup,energies(ileveli(iselect)),energies(ilevelf(iselect))
          !
          if ( tranfreq(iselect)<nu1-1e-6 .or. tranfreq(iselect)>=(nu2-1e-6) .or. j>jmax.or.j>jmax_inp.or. &
                energies(ileveli(iselect))<enermaxlow.or.energies(ilevelf(iselect))>enermaxup) iselect = iselect - 1
          !
          Nlines(j)= Nlines(j) + 1
          !
          !print*,idf,idi,ilevelf(iselect),ileveli(iselect),iselect
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
        print('((i10,2i10,e17.9,f15.6))'),nswap_,ilevelf(iselect),ileveli(iselect),acoef(iselect),tranfreq(iselect)
        !
        print('(i10,a)'),iselect," transitions"
        !
        nselect = iselect
        !
        if (nswap_<1) then 
           write(6,"('nswap = ')") nswap_
           stop 'nswap_=0'
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
