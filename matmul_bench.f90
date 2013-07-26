! Janne Blomqvist 2005, 2011
! Test matrix multiplication performance with the intrinsic matmul and
! blas for different array sizes.
program matmul_bench
  implicit none

  integer, parameter :: sp = selected_real_kind(4), &
       dp = selected_real_kind(15), &
       i64 = selected_int_kind(18)
  integer :: cnt, ii, nmax
  character(len=20) :: arg
  real(dp) :: mult
  
  cnt = command_argument_count()
  if (cnt >= 1) then
     call get_command_argument(1, arg)
     read(arg, *) nmax
  else
     nmax = 2500
  end if
  if (cnt >= 2) then
     call get_command_argument(2, arg)
     read(arg, *) mult
  else
     mult = 2
  end if

  call runsbench (nmax, mult)
  call rundbench (nmax, mult)
  call runlbench (nmax, mult)

contains

  ! Run single precision matrix mult benchmark with different sized arrays.
  subroutine runsbench (nmax, mult)
    integer, intent(in) :: nmax
    real(sp), allocatable, dimension(:,:) :: a, b, res
    integer :: n, loop
    real(dp) :: time, flops, time2, mult

    print *, ' Single precision matrix multiplication test'
    print *, ' Matrix size    Matmul (Gflops/s)     sgemm (Gflops/s)'
    print *, ' ====================================================='
    n = 2
    allocate (a(nmax,nmax), b(nmax,nmax), res(nmax,nmax))
    call random_number (a)
    call random_number (b)
    do
       res(1:n,1:n) = 0.0_sp
       call loops (n, loop, flops)
       call smatmul_timing (a(1:n,1:n), b(1:n,1:n), res(1:n,1:n), loop, time)
       res(1:n,1:n) = 0.0_sp
       call sgemm_timing (n, a, b, res, loop, time2)
       print '(1X,I12,12X,F9.3,12X,F9.3)', n, &
            flops * real(loop, dp) / time / 1.0e9_dp, &
            flops * real (loop, dp) / time2 / 1.0e9_dp
       n = ceiling(n * mult)
       if (n > nmax) exit
    end do
    deallocate (a, b, res)
  end subroutine runsbench

  ! Run double precision matrix mult benchmark with different sized arrays.
  subroutine rundbench (nmax, mult)
    integer, intent(in) :: nmax
    real(dp), allocatable, dimension(:,:) :: a, b, res
    integer :: n, loop
    real(dp) :: time, flops, time2, mult

    print *, ' Double precision matrix multiplication test'
    print *, ' Matrix size    Matmul (Gflops/s)     dgemm (Gflops/s)'
    print *, ' ====================================================='
    n = 2
    allocate (a(nmax,nmax), b(nmax,nmax), res(nmax,nmax))
    call random_number (a)
    call random_number (b)
    do
       res(1:n,1:n) = 0.0_dp
       call loops (n, loop, flops)
       call dmatmul_timing (a(1:n,1:n), b(1:n,1:n), res(1:n,1:n), loop, time)
       res(1:n,1:n) = 0.0_dp
       call dgemm_timing (n, a, b, res, loop, time2)
       print '(1X,I12,12X,F9.3,12X,F9.3)', n, &
            flops * real(loop, dp) / time / 1.0e9_dp, &
            flops * real (loop, dp) / time2 / 1.0e9_dp
       n = ceiling(n * mult)
       if (n > nmax) exit
    end do
    deallocate (a, b, res)
  end subroutine rundbench

  ! Run logical matrix mult benchmark with different sized arrays.
  subroutine runlbench (nmax, mult)
    integer, intent(in) :: nmax
    logical, allocatable, dimension(:,:) :: a, b, res
    real(dp), allocatable, dimension(:,:) :: rtmp
    integer :: n, loop
    real(dp) :: time, ops, mult

    print *, ' Default kind logical matrix multiplication test'
    print *, ' Matrix size    Matmul (Gops/s)'
    print *, ' =============================='
    n = 2
    allocate (a(nmax,nmax), b(nmax,nmax), res(nmax,nmax), rtmp(nmax,nmax))
    call random_number (rtmp)
    a = .false.
    where (rtmp > 0.5)
       a = .true.
    end where
    call random_number (rtmp)
    b = .false.
    where (rtmp > 0.5)
       b = .true.
    end where
    do
       res(1:n,1:n) = .false.
       call loops (n, loop, ops)
       call lmatmul_timing (a(1:n,1:n), b(1:n,1:n), res(1:n,1:n), loop, time)
       print '(1X,I12,10X,F9.3)', n, &
            ops * real(loop, dp) / time / 1.0e9_dp
       n = ceiling(n * mult)
       if (n > nmax) exit
    end do
    deallocate (a, b, res, rtmp)
  end subroutine runlbench  

  ! How many loops to do, to reduce clock jitter.
  subroutine loops (n, loop, ops)
    integer, intent(in) :: n
    integer, intent(out) :: loop
    real(dp), intent(out) :: ops

    ! matmul for square matrix is (2n-1)*n**2 ops.
    ops = (2.0_dp * real (n, dp) - 1.0_dp) * real (n, dp)**2
    ! Assuming an on average 1 gop/s cpu, 1e8 gops takes about 0.1 second and
    ! should be enough. We also do a maximum of 1e4 loops, since
    ! for small arrays the overhead is large.
    loop = max (min (int (1.0e8_dp / ops), 10**5), 1)
  end subroutine loops

  ! Actual routine, and timing.
  subroutine smatmul_timing (a, b, res, loop, time)
    real(sp), intent(in), dimension(:, :) :: a, b
    real(sp), intent(inout) :: res(:,:)
    integer, intent(in) :: loop
    real(dp), intent(out) :: time
    integer(i64) :: t1, t2, rate
    integer :: i

    if (size (a, 1) < 300) then
       ! Do a dry run.
       res = matmul (a, b)
    end if
    call system_clock (t1, rate)
    do i = 1, loop
       res = matmul (a, b)
    end do
    call system_clock (t2)
    time = real((t2 - t1), dp) / rate
  end subroutine smatmul_timing

  ! Actual routine, and timing.
  subroutine dmatmul_timing (a, b, res, loop, time)
    real(dp), intent(in), dimension(:, :) :: a, b
    real(dp), intent(inout) :: res(:,:)
    integer, intent(in) :: loop
    real(dp), intent(out) :: time
    integer(i64) :: t1, t2, rate
    integer :: i

    if (size (a,1) < 300) then
       res = matmul (a, b)
    end if
    call system_clock (t1, rate)
    do i = 1, loop
       res = matmul (a, b)
    end do
    call system_clock (t2)
    time = real((t2 - t1), dp) / rate
  end subroutine dmatmul_timing

  ! Actual routine, and timing.
  subroutine lmatmul_timing (a, b, res, loop, time)
    logical, intent(in), dimension(:, :) :: a, b
    logical, intent(inout) :: res(:,:)
    integer, intent(in) :: loop
    real(dp), intent(out) :: time
    integer(i64) :: t1, t2, rate
    integer :: i

    if (size (a, 1) < 300) then
       res = matmul (a, b)
    end if
    call system_clock (t1, rate)
    do i = 1, loop
       res = matmul (a, b)
    end do
    call system_clock (t2)
    time = real((t2 - t1), dp) / rate
  end subroutine lmatmul_timing

  subroutine sgemm_timing (n, a, b, res, loop, time)
    real(sp), intent(in), dimension(:, :) :: a, b
    real(sp), intent(inout) :: res(:,:)
    integer, intent(in) :: n, loop
    real(dp), intent(out) :: time
    integer(i64) :: t1, t2, rate
    integer :: i, nmax

    nmax = size (a, 1)
    if (n < 300) then
       call sgemm('n','n',n, n, n, 1.0_sp, a(1,1), nmax, b(1,1), nmax, &
            0.0_sp, res(1,1), nmax)
    end if
    call system_clock (t1, rate)
    do i = 1, loop
       call sgemm('n','n',n, n, n, 1.0_sp, a(1,1), nmax, b(1,1), nmax, &
            0.0_sp, res(1,1), nmax)
    end do
    call system_clock (t2)
    time = real((t2 - t1), dp) / rate
  end subroutine sgemm_timing

  subroutine dgemm_timing (n, a, b, res, loop, time)
    real(dp), intent(in), dimension(:, :) :: a, b
    real(dp), intent(inout) :: res(:,:)
    integer, intent(in) :: n, loop
    real(dp), intent(out) :: time
    integer(i64) :: t1, t2, rate
    integer :: i, nmax

    nmax = size (a, 1)
    if (n < 300) then
       call dgemm('n','n',n, n, n, 1.0_dp, a(1,1), nmax, b(1,1), nmax, &
            0.0_dp, res(1,1), nmax)
    end if
    call system_clock(t1, rate)
    do i = 1, loop
       call dgemm('n','n',n, n, n, 1.0_dp, a(1,1), nmax, b(1,1), nmax, &
            0.0_dp, res(1,1), nmax)
    end do
    call system_clock(t2)
    time = real((t2 - t1), dp) / rate
  end subroutine dgemm_timing

    
end program matmul_bench
