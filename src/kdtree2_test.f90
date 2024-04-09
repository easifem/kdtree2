
!
! MAIN PROGRAM HERE.  This is just an example so you know
! how to use it.   This file is in the public domain.
!

MODULE time_kdtree
USE kdtree2_precision_module
USE kdtree2_module
CONTAINS

REAL FUNCTION time_search(tree, nsearch, mode, nn, r2)
  !
  !  Return CPU time, in seconds, for searching 'nsearch' reference points
  !  using any specific search mode.
  !

  TYPE(kdtree2), POINTER :: tree
  INTEGER, INTENT(in) :: nsearch ! how many reference points
  INTEGER, INTENT(in) :: mode ! what kind of search
  INTEGER, INTENT(in) :: nn ! number of neighbors
  REAL(kdkind), INTENT(in) :: r2 ! radius^2
  !
  REAL(kdkind) :: qv(tree%dimen), rv ! query vector, random variate
  INTEGER :: i, random_loc, nf
  REAL :: t0, t1
  TYPE(kdtree2_result), ALLOCATABLE :: results(:)

  REAL(KIND(0.0D0)) :: nftotal
  CALL CPU_TIME(t0) ! get initial time.
  nftotal = 0.0

  ALLOCATE (results(nn))
  DO i = 1, nsearch

    SELECT CASE (mode)
    CASE (1)
      !
      !  Fixed NN search around randomly chosen point
      !
      CALL RANDOM_NUMBER(qv)
      CALL kdtree2_n_nearest(tp=tree, qv=qv, nn=nn, results=results)
    CASE (2)
      !
      ! Fixed NN seasrch around randomly chosen point in dataset
      ! with 100 correlation time
      !
      CALL RANDOM_NUMBER(rv)
      random_loc = FLOOR(rv * tree%n) + 1
      CALL kdtree2_n_nearest_around_point(tp=tree, idxin=random_loc, &
                                       correltime=100, nn=nn, results=results)
    CASE (3)
      !
      ! fixed r2 search
      !
      CALL RANDOM_NUMBER(qv)
      CALL kdtree2_r_nearest(tp=tree, qv=qv, r2=r2, nfound=nf, &
                             nalloc=nn, results=results)
      nftotal = nftotal + nf
    CASE default
      WRITE (*, *) 'Search type ', mode, ' not implemented.'
      time_search = -1.0 ! invalid
      RETURN
    END SELECT

  END DO
  CALL CPU_TIME(t1)

  time_search = t1 - t0
  IF (nftotal .GT. 0.0) THEN
!       write (*,*) 'Average number of neighbors found = ', nftotal / real(nsearch)
  END IF
  DEALLOCATE (results)
  RETURN
END FUNCTION time_search

REAL FUNCTION searches_per_second(tree, mode, nn, r2) RESULT(res)
  !
  !
  ! return estimated number of searches per second.
  ! Will call "time_search" with increasing numbers of reference points
  ! until CPU time taken is at least 1 second.
  !
  TYPE(kdtree2), POINTER :: tree
  INTEGER, INTENT(in) :: mode ! what kind of search
  INTEGER, INTENT(in) :: nn ! number of neighbors
  REAL(kdkind), INTENT(in) :: r2 ! radius^2
  !
  INTEGER :: nsearch
  REAL :: time_taken

  nsearch = 50 ! start with 50 reference points
  DO
    time_taken = time_search(tree, nsearch, mode, nn, r2)
    IF (time_taken .LT. 1.0) THEN
!          write (*,*) 'Intermediate result : ', time_taken, &
!           'sec, nref=',nsearch
      nsearch = nsearch * 5
      CYCLE
    ELSE
      res = REAL(nsearch) / time_taken
      RETURN
    END IF
  END DO
  RETURN
END FUNCTION searches_per_second

REAL FUNCTION average_number_within_ball(tree, navg, r2) RESULT(res)
  ! return the arithmetical average number of points within ball of size
  ! 'r2'
  !
  !  log(avg) = 1/navg log(N_within_ball(i))
  TYPE(kdtree2), POINTER :: tree
  INTEGER, INTENT(in) :: navg
  REAL(kdkind), INTENT(in) :: r2
  !
  INTEGER :: i, cnt
  REAL(kdkind) :: sum, qv(tree%dimen)
  sum = 0.0
  DO i = 1, navg
    CALL RANDOM_NUMBER(qv)
    cnt = kdtree2_r_count(tree, qv, r2)
    IF (cnt .GT. 0) sum = sum + REAL(cnt)
  END DO

  res = sum / REAL(navg)
  RETURN
END FUNCTION average_number_within_ball

END MODULE time_kdtree

PROGRAM kd_tree_test
USE kdtree2_module
USE time_kdtree

INTEGER :: n, d
REAL(kdkind), DIMENSION(:, :), ALLOCATABLE :: my_array
REAL(kdkind), ALLOCATABLE :: query_vec(:)

TYPE(kdtree2), POINTER :: tree, tree2, tree3
! this is how you declare a tree in your main program

INTEGER :: k

TYPE(kdtree2_result), ALLOCATABLE :: results(:), resultsb(:)
INTEGER :: nnbrute, rind
REAL :: t0, t1, sps, avgnum, maxdeviation
REAL(kdkind) :: rv
INTEGER, PARAMETER :: nnn = 5
INTEGER, PARAMETER :: nr2 = 5
REAL r2array(5)
DATA r2array/1.0E-4, 1.0E-3, 5.0E-3, 1.0E-2, 2.0E-2/
INTEGER :: nnarray(5)
DATA nnarray/1, 5, 10, 25, 500/

PRINT *, "Type in N and d"
READ *, n, d

ALLOCATE (my_array(d, n))
ALLOCATE (query_vec(d))
CALL RANDOM_NUMBER(my_array) !fills entire array with built-in-randoms

CALL CPU_TIME(t0)
tree => kdtree2_create(my_array, sort=.FALSE., rearrange=.FALSE.) ! this is how you create a tree.
CALL CPU_TIME(t1)
  write (*,*) real(n)/real(t1-t0), ' points per second built for non-rearranged tree.'

CALL CPU_TIME(t0)
tree2 => kdtree2_create(my_array, sort=.FALSE., rearrange=.TRUE.) ! this is how you create a tree.
CALL CPU_TIME(t1)
  write (*,*) real(n)/real(t1-t0), ' points per second built for rearranged tree.'

CALL CPU_TIME(t0)
tree3 => kdtree2_create(my_array, sort=.TRUE., rearrange=.TRUE.) ! this is how you create a tree.
CALL CPU_TIME(t1)
  write (*,*) real(n)/real(t1-t0), ' points per second built for 2nd rearranged tree.'

nnbrute = 50
ALLOCATE (results(nnbrute), resultsb(nnbrute))

WRITE (*, *) 'Comparing search of ', nnbrute, ' neighbors to brute force.'
DO k = 1, 50
  !
  !
  CALL RANDOM_NUMBER(rv)
  rind = FLOOR(rv * tree3%n) + 1
  query_vec = my_array(:, rind)

  ! find five nearest neighbors to.
  results(:)%idx = -666
  resultsb(:)%idx = -777

     call kdtree2_n_nearest_brute_force(tp=tree3, qv=query_vec, nn=nnbrute, results=resultsb)

     call kdtree2_n_nearest_around_point(tp=tree3,idxin=rind,correltime=-1,nn=nnbrute, results=results)
  ! negative 1 correlation time will get all points.

  maxdeviation = MAXVAL(ABS(results(1:nnbrute)%dis - resultsb(1:nnbrute)%dis))
  IF (ANY(results(1:nnbrute)%idx .NE. resultsb(1:nnbrute)%idx) .OR. &
      (maxdeviation .GT. 1.0E-8)) THEN
    WRITE (*, *) 'MISMATCH! @ k=', k

    PRINT *, "Tree indexes    = ", results(1:nnbrute)%idx
    PRINT *, "Brute indexes   = ", resultsb(1:nnbrute)%idx
        print *, "Tree-brute distances  = ", results(1:nnbrute)%dis- resultsb(1:nnbrute)%dis
  END IF
END DO

10 FORMAT('R^2 search, r2/d=', G10.2, ':', F13.0, A)
20 FORMAT(A, ' NN=', I7, ':', F10.0, ' searches/s in ', A)

IF (.TRUE.) THEN
  ! simple testing to mimic kdtree_test_old
  !
  DO k = 1, nnn
    sps = searches_per_second(tree2, 1, nnarray(k), 1.0_KDKIND)
    WRITE (*, 20) 'Random pts', nnarray(k), sps, 'rearr. tree.'

  END DO
  DO k = 1, nnn
    sps = searches_per_second(tree2, 2, nnarray(k), 1.0_KDKIND)
    WRITE (*, 20) 'in-data pts', nnarray(k), sps, 'reearr. tree.'
  END DO
ELSE
  DO k = 1, nr2
    avgnum = average_number_within_ball(tree2, 500, r2array(k) * &
                                        REAL(d, kdkind))
  WRITE (*, *) 'Avg number within ', r2array(k) * REAL(d, kdkind), '=', avgnum
  END DO

  DO k = 1, nnn
    sps = searches_per_second(tree, 1, nnarray(k), 1.0_KDKIND)
    WRITE (*, 20) 'Random pts', nnarray(k), sps, 'regular tree.'
    sps = searches_per_second(tree2, 1, nnarray(k), 1.0_KDKIND)
    WRITE (*, 20) 'Random pts', nnarray(k), sps, 'rearr. tree.'
    sps = searches_per_second(tree3, 1, nnarray(k), 1.0_KDKIND)
    WRITE (*, 20) 'Random pts', nnarray(k), sps, 'rearr./sorted tree.'
  END DO

  DO k = 1, nnn
    sps = searches_per_second(tree, 2, nnarray(k), 1.0_KDKIND)
    WRITE (*, 20) 'in-data pts', nnarray(k), sps, 'regular tree.'
    sps = searches_per_second(tree2, 2, nnarray(k), 1.0_KDKIND)
    WRITE (*, 20) 'in-data pts', nnarray(k), sps, 'reearr. tree.'
    sps = searches_per_second(tree3, 2, nnarray(k), 1.0_KDKIND)
    WRITE (*, 20) 'in-data pts', nnarray(k), sps, 'rearr./sorted tree.'
  END DO

  DO k = 1, nr2
    sps = searches_per_second(tree, 3, 20000, r2array(k) * REAL(d, kdkind))
    WRITE (*, 10) r2array(k), sps, ' searches/s'
    sps = searches_per_second(tree2, 3, 20000, r2array(k) * REAL(d, kdkind))
    WRITE (*, 10) r2array(k), sps, ' searches/s in rearranged tree'
    sps = searches_per_second(tree3, 3, 20000, r2array(k) * REAL(d, kdkind))
    WRITE (*, 10) r2array(k), sps, ' searches/s in rearranged/sorted tree'
  END DO

END IF
CALL kdtree2_destroy(tree)
CALL kdtree2_destroy(tree2)
CALL kdtree2_destroy(tree3)
! this releases memory for the tree BUT NOT THE ARRAY OF DATA YOU PASSED
! TO MAKE THE TREE.

DEALLOCATE (my_array)
! deallocate the memory for the data.

END PROGRAM kd_tree_test
