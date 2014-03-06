module push_array
  ! push a val onto the allocatable array
  ! Subsequent val must match the dimension of the first val pushed,
  ! otherwise a runtime SIGSEV - invalid memory reference will
  ! be throwan at runtime.
  ! usage:
  ! array = push(array,val)
  interface push
     module procedure push_scalar_int_onto_rank1_int2
     module procedure push_rank1_int_onto_rank2_int
     module procedure push_rank1_real8_onto_rank2_real8
  end interface push
contains
  function push_scalar_int_onto_rank1_int2 (array,val) result (new_array)
    integer(2),intent(in),allocatable :: array(:)
    integer,intent(in) :: val
    integer(2),allocatable :: new_array(:)
    integer :: length
    if (allocated(array)) then
       length = size(array) + 1
    else
       length = 1
    end if
    allocate(new_array(size(array) + 1))
    if (allocated(array)) new_array(:) = array(:)
    new_array(length) = val
    return
  end function push_scalar_int_onto_rank1_int2
  function push_rank1_int_onto_rank2_int (array,val) result (new_array)
    integer,intent(in),allocatable :: array(:,:)
    integer,intent(in) :: val(:)
    integer,allocatable :: new_array(:,:)
    integer :: length
    if (allocated(array)) then
       length = size(array,2) + 1
    else
       length = 1
    end if
    allocate(new_array(1:size(val),length))
    if (allocated(array)) new_array(1:size(val),:) = array(1:size(val),:)
    new_array(1:size(val),length) = val
    return
  end function push_rank1_int_onto_rank2_int
  function push_rank1_real8_onto_rank2_real8 (array,val) result (new_array)
    real(8),intent(in),allocatable :: array(:,:)
    real(8),intent(in) :: val(:)
    real(8),allocatable :: new_array(:,:)
    integer :: length
    if (allocated(array)) then
       length = size(array,2) + 1
    else
       length = 1
    end if
    allocate(new_array(1:size(val),length))
    if (allocated(array)) new_array(1:size(val),:) = array(1:size(val),:)
    new_array(1:size(val),length) = val
    return
  end function push_rank1_real8_onto_rank2_real8
end module push_array
