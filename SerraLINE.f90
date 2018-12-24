!						12/12/2018
! SerraLINE: MAIN
! Calculates bending angles

  program SerraLINE

  use IO_mod
  use parms
  use functions_mod

  implicit none

  logical :: circle_str
  integer :: nbp, box, str, frames, ierror, &
           & ndim, mdim, & !mdim is the extended dimension
           & i
  integer, allocatable :: ring_index(:)
  character(1), allocatable :: seq(:)
  character(360) :: top, traj

  real(dp), allocatable :: coords(:,:,:), & !coordinates
          & tangents(:,:,:), bends(:,:), &  !tangent vectors, bending angles,
          & avstd_bends(:,:)             !average and standard dev of bending

  !INPUT SECTION
  !-----------------------------------------------------------------------------------
  !Get the directories of trajectory and topology, and 
  !read if its a closed structure or not
  call SerraLINE_inputs(traj,top,circle_str)
  str=1 !Reading will be performed as if it were a single stranded structure

  !READING SECTION
  !-----------------------------------------------------------------------------------
  write(6,*) "Reading topology file"
  call topology_amber(nbp,box,seq,ring_index,top,str)

  !We do not care about ring atoms, we suppose that the topology file only 
  !contains the atoms needed to draw the curvature

  !Reads all the coordinates from the amber trajectory file
  !It assumes that the number of bp is the same number of atoms
  ! e.g. there's not more atoms in the trajectory/topology
  write(6,*) "Reading trajectory file"
  call coordinates_amber_WrLINE(coords,nbp,frames,traj,box)

  !-----------------------------------------------------------------------------------

  !TANGENT VECTORS SECTION
  !-----------------------------------------------------------------------------------
 
  !Get tangent vectors to the curve (rj-ri)
 
  !If the structure is not a circle then the last bp does not has tangent vector
  write(6,*) "Calculating tangent vectors"
  if (circle_str) then
    allocate(tangents(3,nbp,frames), stat=ierror)
  else
    allocate(tangents(3,nbp-1,frames), stat=ierror)
  end if
  if(ierror/=0) stop "Error in allocating tangent vectors"
 
  call get_tangent_vectors(coords,nbp,frames,circle_str,tangents)
 
  !Coordinates are not longer needed
  deallocate(coords, stat=ierror)
  if(ierror/=0) stop "Error in deallocating coordinates"

  !-----------------------------------------------------------------------------------
 
  !BEDINGS SECTION 
  !-----------------------------------------------------------------------------------

  !Preparing for calculating bends
  if (circle_str) then
    ndim = nbp
  else
    ndim = nbp-1
  end if

  mdim = ndim*(ndim-1)/2
  allocate(bends(frames,mdim), avstd_bends(2,mdim), stat=ierror)
  if(ierror/=0) stop "Error in allocating bendings"

  !Calculate bending angles
  write(6,*) "Calculating bending angles"
  call get_bendings(ndim,frames,tangents,bends)
 
  deallocate(tangents, stat=ierror) !No longer needed
  if (ierror/=0) stop "Error in deallocating tangent vectors"

  !Averages, standard deviations SECTION 
  !-----------------------------------------------------------------------------------

  !Now calculate averages and standard deviations using function
  !average_std()
  write(6,*) "Calculating averages and standard deviations of bendings"
  do i=1,mdim
    avstd_bends(:,i) = average_std(bends(1:frames,i),frames)
  end do

  deallocate(bends, stat=ierror) !No longer needed
  if (ierror/=0) stop "Error in deallocating bendings"

  !-----------------------------------------------------------------------------------

  !WRITING SECTION
  !-----------------------------------------------------------------------------------
  call write_av_bends(circle_str,nbp,ndim,mdim,frames,seq,avstd_bends)
 
  end program SerraLINE
