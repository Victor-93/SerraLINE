!------------------------------------------------------------- 12/12/2018
! This module contains functions and subroutines needed for reading
! and writing.
  module IO_mod

  use parms
  implicit none

  private arrdims
  public topology_amber, coordinates_amber_WrLINE, SerraLINE_inputs

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

  contains

!-----------------------------------------------------------------------
!  Reads topology file from amber topology file
!    Gets number bp (nbp), checks if there's a "box" and the sequence on
!    both strands

  subroutine topology_amber(nbp,box,seq,ring_index,top,str)
  implicit none
  integer ::  nbp, box, str                !if box == 1 then there is a box
  integer, allocatable :: res_index(:), seq_index(:), ring_index(:)
  integer :: i, j, ierror, n_atoms, n_res, iaux(10), nlines, l, f_lines, &
           & pur, pyr  !purines, pyrimidines and type of structure
  character(360) :: top, caux                 !caux helps in reading
  character(1), allocatable, intent(out) :: seq(:)
  character(4), allocatable :: atom_names(:), res_names(:)

  !Get dimension of file top. Here, "l" is a dummy variable
  call arrdims(top,f_lines,l)

! READING TOPOLOGY----------------------------------------------------------------
  !Open topology file
  open(unit=10, file=trim(top), status="old",action="read",iostat=ierror)
    if (ierror/=0) stop 'Error in opening topology file'

    do j=1,f_lines
      read(10,"(A)") caux

      !Escape statement: if three quarters of file have been readed, then it
      !  it would automatically exit the loop. Hopefully, this will never happend,
      !  hence, it only is an emergency escape before collapsing.
      if (j > 3*f_lines / 4) exit

      !FLAG POINTERS
      if ( trim(caux) == "%FLAG POINTERS" ) then

        do                                       !In case there are any comments
          read(10,"(A)") caux                              
          if (caux(1:7) == "%FORMAT" ) exit                 
        end do

        read(10, 1001) iaux                                 !Has dimension 10 
        n_atoms=iaux(1)                                     !n_atoms
        read(10, 1001) iaux  
        n_res=iaux(2)                                       !n_res
        read(10, 1001) iaux 
        box=iaux(8)                                         !box

        !If there are no atoms stop.
        if (n_atoms <= 0) stop "No atoms in topology"
        !In case of no residues the program will try to
        !identify them (it could go wrong)

        !Print info
        write(6,"(1A20,1I10)") "Number of atoms = ", n_atoms
        write(6,"(1A20,1I10)") "Number of residues = ", n_res
        write(6,"(1A20,1I10)") "BOX = ", box

        !allocate atom_names
        allocate(atom_names(n_atoms), stat=ierror)
        if(ierror/=0) stop 'Error in allocating atom_names'
        !Not point in allocating residues since we can try identify them 
      end if
      !End of POINTERS section 


      !FLAG ATOM_NAME
      if ( trim(caux) == "%FLAG ATOM_NAME" ) then

         do                                       !In case there are any comments
          read(10,"(A)") caux                              
          if (caux(1:7) == "%FORMAT" ) exit                 
        end do

        !calculate lines in ATOM_NAME section
        if ( mod(n_atoms,20) == 0 ) then
          nlines = n_atoms / 20                       !20 atoms per line
        else
          nlines = n_atoms / 20 + 1        
        end if

        !READ atoms names
        l=1
        if (nlines > 1) then                          !in case (n_atoms < 20) which is unlikely
          do i=1,nlines-1
            read(10, 1002) atom_names(l:19+l)
            l=l+20
          end do
        end if
        read(10, 1002) atom_names(l:n_atoms)
      end if
      !End of ATOM_NAME SECTION


      !FLAG RESIDUE_LABEL
      if ( trim(caux) == "%FLAG RESIDUE_LABEL") then

         do                                       !In case there are any comments
          read(10,"(A)") caux                              
          if (caux(1:7) == "%FORMAT" ) exit                 
        end do

        !calculate lines in RESIDUE_LABEL section
        if ( mod(n_res,20) == 0 ) then
          nlines = n_res / 20                        !20 residues per line
        else
          nlines = n_res / 20 + 1                 
        end if
        
        !allocate res_names (is a temporary variable)
        allocate(res_names(n_res), stat=ierror)
        if (ierror /= 0 ) stop "Error in allocating res_names"

        !Read residues names
        l=1
        if (nlines > 1) then                         !in case (nlines < 20) which is unlikely
          do i=1,nlines-1
            read(10, 1002) res_names(l:19+l)
            l=l+20
          end do
        end if
        read(10, 1002) res_names(l:n_res)
      end if
      !End of RESIDUE_LABEL SECTION
      
      !FLAG RESIDUE_POINTER
       if ( trim(caux) == "%FLAG RESIDUE_POINTER") then
 
         do                                       !In case there are any comments
          read(10,"(A)") caux                              
          if (caux(1:7) == "%FORMAT" ) exit                 
        end do

        !calculate lines in RESIDUE_POINTER section
        if ( mod(n_res,10) == 0 ) then
          nlines = n_res / 10                         !10 residues per line
        else
          nlines = n_res / 10 + 1                 
        end if
        
        !allocate res_index (is a temporary variable)
        allocate(res_index(n_res), stat=ierror)
        if (ierror /= 0 ) stop "Error in allocating res_names"

        !Read residues indices (start)
        l=1
        if (nlines > 1) then                          !in case (nlines < 10) which is unlikely
          do i=1,nlines-1
            read(10, 1001) res_index(l:9+l)
            l=l+10
          end do
        end if
        read(10, 1001) res_index(l:n_res)

        !!!!!!!!!!!!!!!!!We have all we need!!!!!!!!!!!!!!!
        exit        

      end if
      !End of RESIDUE POINTER SECTION 
      
    end do
  close(10, iostat=ierror)
  if (ierror /= 0) stop "Error in closing topology file"
! END OF READING -------------------------------------------------------


  !Check if data obtained from topology file is enough to go on
  if ( .not. allocated(atom_names) ) then
    stop "Atoms names could not be identified"
  end if 

  !Check if double stranded structure is complete
  if (str == 2 .and. mod(n_res,2) /=0 ) then
    stop "Mismatch number of residues in double stranded structure"
  end if 
   
! COUNT NBP -------------------------------------------------------
  !The structure will be treated as a double stranded structure or a single strandedi one.
  !In case of dsDNA, if the number of residues in both strands is not the same, then
  !this will be treated as an error. 

  nbp=0
  do i=1,n_res
    !Identify Adenine
    if (any( A_l(:) .eq. res_names(i) ) ) then
      nbp=nbp+1
    end if
    !Identify Guanine
    if (any( G_l(:) .eq. res_names(i) ) ) then
      nbp=nbp+1
    end if
    !Identify Cytosine
    if (any( C_l(:) .eq. res_names(i) ) ) then
      nbp=nbp+1
    end if
    !Identify Thymine
    if (any( T_l(:) .eq. res_names(i) ) ) then
      nbp=nbp+1
    end if
    !Identify Uracil
    if (any( U_l(:) .eq. res_names(i) ) ) then
      nbp=nbp+1
    end if
  end do   !close i
 
       
! GET SEQUENCE -------------------------------------------------------

  !Allocate sequence and sequence index
  allocate(seq(nbp),seq_index(nbp), stat=ierror) 
  if(ierror/=0) stop 'Error in allocating sequence'

  !Get the correct number of nbp and check if structure is incomplete in case of 
  !a double stranded structure
  if ( str == 2 ) then
    if ( mod(nbp,2) /= 0 ) stop "Incomplete double stranded structure"
    nbp = nbp/2                 !We counted bases from both strands
  end if  
  !if str=1 then its a single stranded structure. If its 2 then is double stranded
   
  !Check if molecule is longer than 4bp, if not then stop
  if ( nbp/2 .le. 4 ) stop "Invalid DNA fragment, the fragment must be larger than 4 bp"

  !Similar loop as before but now identifying residues
  l=0 !will help us count bps
  pur = 0 ; pyr = 0 !purines and pyrimidines counters
  do i=1,n_res
    !Identify Adenine
    if (any( A_l(:) .eq. res_names(i) ) ) then
      l=l+1
      seq(l) = "A" 
      seq_index(l) = res_index(i)
      pur = pur +1
    end if
    !Identify Guanine
    if (any( G_l(:) .eq. res_names(i) ) ) then
      l=l+1
      seq(l) = "G" 
      seq_index(l) = res_index(i) 
      pur = pur +1
    end if
    !Identify Cytosine
    if (any( C_l(:) .eq. res_names(i) ) ) then
      l=l+1
      seq(l) = "C" 
      seq_index(l) = res_index(i) 
      pyr = pyr +1
    end if
    !Identify Thymine
    if (any( T_l(:) .eq. res_names(i) ) ) then
      l=l+1
      seq(l) = "T" 
      seq_index(l) = res_index(i) 
      pyr = pyr +1
    end if
    !Identify Uracil
    if (any( U_l(:) .eq. res_names(i) ) ) then
      l=l+1
      seq(l) = "U" 
      seq_index(l) = res_index(i) 
      pyr = pyr +1
    end if
  end do   !close i
     
  !Check we read all correctly the sequence and indices
  if (l /= nbp*str) stop 'Error in extracting sequence and res pointers'

  write(6,*) "Base pairs read = ", nbp

! INDECES OF RING ATOMS -------------------------------------------------------

  !Allocate ring index
  allocate(ring_index(pyr*6+pur*9), stat=ierror)
  if(ierror/=0) stop 'Error in allocating ring_index'

  ring_index = 0 !initializing, this will be true for SerraLINE
  l = 0          !This will be an auxiliar

  do i=1,nbp*str-1 !l is number of nucleotides or pyr+pur
    !For purines
    if (seq(i) .eq. "A" .or. seq(i) .eq. "G" ) then
      do j=seq_index(i),seq_index(i+1)
        if (atom_names(j) .eq. N9) then
          ring_index(l+1) = j
        end if
        if (atom_names(j) .eq. C8) then
          ring_index(l+2) = j
        end if
        if (atom_names(j) .eq. N7) then
          ring_index(l+3) = j
        end if
        if (atom_names(j) .eq. C5) then
          ring_index(l+4) = j
        end if
        if (atom_names(j) .eq. C6) then
          ring_index(l+5) = j
        end if
        if (atom_names(j) .eq. N1) then
          ring_index(l+6) = j
        end if
        if (atom_names(j) .eq. C2) then
          ring_index(l+7) = j
        end if
        if (atom_names(j) .eq. N3) then
          ring_index(l+8) = j
        end if
        if (atom_names(j) .eq. C4) then
          ring_index(l+9) = j
        end if
     end do
     l=l+9
    else 
      !If its pyrimidine
      do j=seq_index(i),seq_index(i+1)
        if (atom_names(j) .eq. N1) then
          ring_index(l+1) = j
        end if
        if (atom_names(j) .eq. C2) then
          ring_index(l+2) = j
        end if
        if (atom_names(j) .eq. N3) then
          ring_index(l+3) = j
        end if
        if (atom_names(j) .eq. C4) then
          ring_index(l+4) = j
        end if
        if (atom_names(j) .eq. C5) then
          ring_index(l+5) = j
        end if
        if (atom_names(j) .eq. C6) then
          ring_index(l+6) = j
        end if
      end do
      l=l+6 
    end if
  end do

 
 1001 format (10I8)   !Format used in FLAG pointers and residue pointers
 1002 format (20A4)   !Format used in atom and residue names

  end subroutine topology_amber
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
 !Reads the atom coordinates from WrLINE coordinate file (amber style)

 subroutine coordinates_amber_WrLINE(coords,N,frames,traj,box)
 real(dp), allocatable :: coords(:,:,:), b(:)
 integer :: nlin, ncol, frames, N, box, ierror, i, k, l, rl, &
          & rn
 character(360) :: traj, aux

  !get dimensions of trajectory file
  call arrdims(traj,nlin,ncol)

  ncol=3*N/10 !Let's recycle ncol
  if (3*N-10*ncol > 0) then
    ncol = ncol+1   !the extra line
  end if

  !Get number of frames in file

  if (box /= 0) then           !Just in case let's put it here
    frames = nlin/(ncol+1)
  else
    frames = nlin/ncol
  end if

  write(6,"(1A20,1I10)") "Number of frames =", frames
  !Lets recycle nlin, it will be the number of lines that can be readed completly by 3
  nlin = N/10  

  !rl is the number of remaining lines
  rl=3*mod(N,10)
  if ( mod(rl,10) == 0 ) then
    rl = rl/10
  else
    rl = 1+ rl/10
  end if

  !rn is number of atoms remaining
  rn = N - 10*nlin    !Its 10 atoms per nlin (3 lines)
 
  allocate(coords(3,N,frames), b(rn*3), stat=ierror)
  if (ierror /= 0) stop "Error in allocating coordinates array"

  open(unit=10, file=trim(traj), status="old",action="read",iostat=ierror)
    if (ierror/=0) stop 'Error in opening trajectory file'

  !Begin reading
  read(10,*) aux !skips first line

  do k=1,frames

  l = 0

    do i=1,nlin
      read(10, 1003) coords(1:3,l+1,k), coords(1:3,l+2,k), &
                   & coords(1:3,l+3,k), coords(1,l+4,k)

      read(10, 1003) coords(2:3,l+4,k), coords(1:3,l+5,k), &
                   & coords(1:3,l+6,k), coords(1:2,l+7,k)

      read(10, 1003) coords(3,l+7,k), coords(1:3,l+8,k), &
                   & coords(1:3,l+9,k), coords(1:3,l+10,k)
      l = l+10
    end do !close i

    ! If there are more lines to read
    if (rl > 0) then
      !If there's one more
      if (rl >= 1) then
        if (rl < 2) then
          read(10, 1003) b(1:rn*3)
        else
          read(10, 1003) b(1:10)
        end if
      end if 
    
      !If there's more than one
      if (rl >= 2) then
        if (rl < 3) then
          read(10, 1003) b(11:rn*3)
        else
          read(10, 1003) b(11:20)
        end if
      end if 
    
      !If there's 3 left
      if (rl == 3) then
        read(10, 1003) b(21:rn*3)
      end if 
      
      !Now, collect the remaining coordinates
      do i=1,rn
        coords(1,l+i,k) = b((i-1)*3+1)
        coords(2,l+i,k) = b((i-1)*3+2)
        coords(3,l+i,k) = b((i-1)*3+3)
      end do

    end if !close rl>0 
    
    !Just in case that there is a box, skip this line
    if (box /= 0) then
      read(10,*) aux
    end if
 
  end do !close k

  close(10) !close trajectory

 1003 format (10F8.3) !Format used in amber trajectory file 

 end subroutine coordinates_amber_WrLINE
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
 !This subroutine is used in SerraLINE
 ! It writes the output file with the bending angles
 subroutine write_av_bends(circle_str,nbp,ndim,mdim,frames,seq,avstd_bends)
 implicit none

 logical, intent(in) :: circle_str
 integer, intent(in) :: nbp, ndim, mdim, frames
 character(1), intent(in) :: seq(nbp)
 real(dp), intent(in) :: avstd_bends(2,mdim)

 character(360) :: output_file
 integer :: i, j, l, ierror

 output_file = "bending_parameters.out"

 open(unit=10, file=trim(output_file), status="replace", action="write", iostat=ierror)
 if (ierror/=0) stop 'Error in opening output bending file'

  write(10,*) "BENDING PARAMETERS"
  write(10,*) ""
  if (circle_str) then
    write(10,*) "CIRCULAR STRUCTURE"
  else
    write(10,*) "NON CIRCULAR STRUCTURE"
  end if
  write(10,*) "BASE PAIRS ", ndim
  write(10,*) "SEQUENCE"
  write(10,*) "", seq(1:nbp)
  write(10,*) "SNAPSHOTS ANALYSED", frames
  write(10,*) ""
  write(10,*) "First column averages, second column standard deviations"
  l=0
  do j=1,ndim-1
    write(10,1007) j, "bp"
    write(10,1005) "base-step", "Bending angle"
    write(10,*) "--------------------------------------"

    do i=1,ndim-j
      l=l+1
      write(10,1006) i, "-", i+j, " "," ",seq(i),seq(i+j), " "," ", &
               & avstd_bends(1,l), avstd_bends(2,l)
    end do
    write(10,*) ""
  end do

 close(10) 

 1005 format(A15,A20)
 1006 format(I4,A1,I4,6A1,2F10.3)
 1007 format(I10,A10)
! 1008 format(A25)
! 1009 format(A25,I10)

 end subroutine write_av_bends

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Obtains the lines and columns of the file "filen"
 subroutine arrdims(filen,nlin,ncol)
 character(260), intent(in) :: filen
 integer, intent(out) :: nlin,ncol
 character(260) :: fwc
 integer :: ierr

   fwc='.functionsmod_wc-ou.tmp'   ! auxilliary file, where we store number
                                   ! of lines and columns separated by newline

   call system('head -n 1 '//filen//' | wc -w  > '//trim(fwc),ierr)
   if (ierr == 0) then
     call system('cat '//filen//' | wc -l  >> '//trim(fwc),ierr)
     if (ierr == 0) then
       open(unit=10,file=trim(fwc),status='old',action='read')
       read(10,*) ncol
       read(10,*) nlin        ! We obtain the number of lines and number of
       close(10)              ! columns of the data file
     else                
       stop                
     end if
   else
     stop
   end if

   if (max(nlin,ncol) == 0) then
     write(6,'(a)') 'Error in arrdims. Did not get valid dimensions'
     write(6,'(a,i10,a,i10)') '  nlin=',nlin,'   ncol=',ncol
     stop
   end if
 end subroutine arrdims

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
  !Subroutine used by SerraLINE.
  ! Reads directory of trajectory (traj) and topology (top) files.
  ! Also, it returns a logical variable circle_str which indicates
  ! If the structure being analysed is a closed structure.
  ! All this information is indicated in the input file "s_line.in
  subroutine SerraLINE_inputs(traj,top,circle_str)
  implicit none

  character(360), intent(out) :: top, traj
  logical, intent(out) :: circle_str
  integer :: str

  !Closed structure?
  call nextpar
  read(5,*) str

  !Topology
  call nextpar
  read(5,"(A)") top
  top = adjustl(top)
 
  !Trajectory
  call nextpar
  read(5,"(A)") traj
  traj = adjustl(traj)

  if (str == 1) then
    circle_str = .true.
    write(6,"(A)") "CLOSED STRUCTURE"
  else if (str == 0) then
    circle_str = .false.
    write(6,"(A)") "OPEN STRUCTURE"
  else
    stop "Tell me if its a closed or open structure (1 or 0)"
  end if           

  write(6,"(2A)") "Topology file =", trim(top) 
  write(6,"(2A)") "Trajectory file =", trim(traj) 

  end subroutine SerraLINE_inputs 
!----------------------------------------------------------------------

!-----------------------------------------------------------------------
  !Reads chr, if first character is different from 'c' (comment line), then 
  !it exits and then next line should be readed
  subroutine nextpar
  implicit none
    character(len=1) :: chr
    do 
      read(5,'(A1)') chr
      if (chr(1:1) /= 'c') exit
    end do

  end subroutine nextpar
!----------------------------------------------------------------------

  end module IO_mod
