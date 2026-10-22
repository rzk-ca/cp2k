!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Types for time-independent DFT methods
!> \par History
!>       2020.04 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE tidft_types
   USE almo_scf_types,                  ONLY: optimizer_options_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE message_passing,                 ONLY: mp_para_env_type
   USE cp_dbcsr_api,                    ONLY: dbcsr_type
   USE kinds,                           ONLY: default_path_length, default_string_length, dp
   USE input_section_types,             ONLY: section_vals_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tidft_types'

   PUBLIC :: tidft_env_type, tidft_env_release

   ! Note that closed-shell DFT (RKS) is not explicitly considered here
   ! because excitations unpair electrons. In the (rare) cases of
   ! excitations of paired electrons, use ROKS with max unpairings set to 0

   ! Key matrices: mat_var - matrices of indpendently varied coefficients,
   !                         there are nvars of them, occupied by:
   !                         pairs, alpha-electrons, beta-electrons
   !               mat_t   - matrices of alpha and beta coefficients,
   !                         there are nspins of them, nspins is 2 with
   !                         the exception multip = 1 and max_upairings = 0

   ! Several methods are implemented:
   !  1.   UKS-mixed - spin-contamination is accepted as unavoidable
   !                   nspins = 2, nvars = 2
   !                   mat_var = mat_t
   !  2.  ROKS-mixed - single Slater determinants are used.
   !                   This leads to mixed spin states,
   !                   unless electrons pair up naturally during energy minimization
   !                   In addition to the desired multiplicity,
   !                   user specifies maximum number of allowed unpairings.
   !                   This number is used to compute how many singly occupied
   !                   orbitals are used in the optimization.
   ! The table lists sizes of mat_var matrices for different settings.
   ! The order of entries singly-occ alpha, singly-occ beta, doubly-occ
   !               max_unpair=0      max_unpair=1      max_upair=2
   !    singlet:    0,0,N/2           1,1,N/2-1         2,2,N/2-2
   !    doublet:    1,0,(N-1)/2       2,1,(N-1)/2-1     3,2,(N-1)/2-2
   !    triplet:    2,0,(N-2)/2       3,1,(N-2)/2-1     4,2,(N-2)/2-2
   !  Special cases: IF (max_unpairings=0) THEN
   !                   IF (singlet) nspins = 1, nvars = 1, mat_var = mat_t
   !                   ELSE nspins = 2, nvars = 2, mat_var != mat_t
   !                 ELSE
   !                   nspins = 2, nvars = 3, mat_var != mat_t
   !                 ENDIF
   !  2B. ROKS-pure  - enforce spin pure states by tracking additional determinants.
   !                   This approach works only for closed-shell singlet ground
   !                   states and single-excitations (i.e. max_unpairings=1).
   !                   NOTE: It automatically produces spin pure states for max_unpairings=0
   !                   The spin-pure flag does not require a separate storage of mat_var
   !                   because the single-beta mat can be interpreted as alpha on-the-fly.
   !                   However, the spin-flipped determinant requires its own alpha and beta
   !                   densities and matrices of proper dimensions have to be created

   ! This structure describes a single Slater determinant
   TYPE slater_det_type

      ! For a single determinant, multiplicity only means
      ! the difference between the number of alpha and beta electrons:
      ! desired_multip = nalpha - nbeta + 1
      ! The desired_multip is contained here (not only in tidft_env) because
      ! ROKS-pure approach tracks states of different mutiplicity
      INTEGER                                            :: desired_multip
      ! Number of alpha and beta electrons (no. of cols in mat_t)
      INTEGER, DIMENSION(2)                              :: nelectrons_ab
      ! Number of alpha and beta virtual states (no. of cols in mat_v)
      INTEGER, DIMENSION(2)                              :: nvirtuals_ab
      ! Number of doubly and singly occupied states (no. of cols in mat_var)
      INTEGER, DIMENSION(3)                              :: ncols_var

      ! Matrix of independently varied MO coefficients
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: mat_var

      ! Occupied MO coefficients NxOCC(spin)
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: mat_t

      ! Virtual MO coefficients NxVIRT(spin)
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: mat_v
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: mat_vo

      ! Projector onto occupied subspace NxN - for each spin
      ! Do not store - it is easy to compute from MOs
      !TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: mat_proj

      ! Kohn-Sham Hamiltonian NxN - for each spin
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: mat_ks
      REAL(KIND=dp)                                      :: slater_energy

   END TYPE slater_det_type

   ! This structure decribes an electronic state
   TYPE e_state_type

      ! Number of Slater determinants describing the electronic state
      INTEGER                                            :: nslater
      TYPE(slater_det_type), DIMENSION(:), ALLOCATABLE   :: slater_det
      REAL(KIND=dp)                                      :: state_energy
      REAL(KIND=dp)                                      :: excitation_energy
      REAL(kind=dp)                                      :: osc_strength

   END TYPE e_state_type

   ! the structure contains general info about the system
   TYPE tidft_env_type

      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env

      LOGICAL                                            :: calc_forces, &
                                                            orthogonal_basis
      LOGICAL                                            :: smear = .FALSE.

      ! Electronic states
      LOGICAL                                            :: restricted ! ROKS or UKS
      ! Spin pure states requested?
      LOGICAL                                            :: spin_pure
      ! Compute triplet correction to obtain spin-pure states?
      LOGICAL                                            :: triplet_correction
      ! Number of slater determinants per state
      INTEGER                                            :: ndets
      ! Number of spins treated indepedently (1 or 2)
      INTEGER                                            :: nspins
      ! Spin factors are 2 for doubly occupied and 1 for signly occupied orbtals
      INTEGER, DIMENSION(3)                              :: spin_factor_var
      ! Maximum number of electron pairs broken beyond allowed by multiplicity
      INTEGER                                            :: unpairings
      ! The multiplicity of states
      INTEGER                                            :: desired_multip
      ! Total number of electronic states
      INTEGER                                            :: nstates
      ! Number of electrons
      INTEGER                                            :: nelectrons
      ! Descriptors of electronic states
      TYPE(e_state_type), DIMENSION(:), ALLOCATABLE      :: e_states

      ! TIDFT control
      LOGICAL, DIMENSION(100)                            :: loss_function_components
      LOGICAL                                            :: state_by_state
      LOGICAL                                            :: pre_orthogonalize_states
      LOGICAL                                            :: numerical_gradient
      REAL(KIND=dp)                                      :: numerical_gradient_dx
      INTEGER                                            :: guess_max_excitation_order
      REAL(KIND=dp)                                      :: guess_min_excitation_energy
      REAL(KIND=dp)                                      :: guess_max_excitation_energy
      LOGICAL                                            :: guess_excite_alpha_only

      ! penalty related variables
      REAL(KIND=dp)                                      :: intrastate_penalty_coeff
      REAL(KIND=dp),  DIMENSION(:), ALLOCATABLE          :: intrastate_nonorthogonality_penalty_coeff
      REAL(KIND=dp)                                      :: interstate_nonorthogonality_penalty_coeff
      ! aeph: augmented Lagrangian quadratic penalty coefficient (mu)
      REAL(KIND=dp)                                      :: augmented_penalty_coeff
      ! aeph: input value of mu, used to reset each penalty-strength pass
      REAL(KIND=dp)                                      :: augmented_penalty_coeff0
      ! aeph: threshold for normalized det(phi), read from PENALTY%EPS_DET
      REAL(KIND=dp)                                      :: eps_det
      ! aeph: safeguard ratio tau 
      REAL(KIND=dp)                                      :: alm_ratio
      REAL(KIND=dp)                                      :: penalty_strength_inc_factor
      REAL(KIND=dp)                                      :: allowed_max_norm_phi_mI
      REAL(KIND=dp)                                      :: initial_guess_random_fraction
      INTEGER                                            :: initial_guess_type
      INTEGER                                            :: initial_guess_random_seed
      ! CINDY
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: grad_wrt_penalty_coeff
      REAL(KIND=dp)                                      :: grad_wrt_inter_penalty_coeff
      INTEGER                                            :: nsteps_penalty_strength
      REAL(KIND=dp)                                      :: loose_eps_scf

      ! Signed-overlap Lagrange multipliers for the interstate orthogonality
      ! constraints of the current state (loss function component 5).
      ! One multiplier per previously optimized state J:
      !   term = SUM_J mu_J * g_IJ,  g_IJ = Phi_IJ / SQRT(Phi_II * Phi_JJ)
      LOGICAL                                            :: multiplier_component
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: overlap_multipliers
      ! grow lambda only if the residual failed to drop by this ratio
      REAL(KIND=dp)                                      :: multiplier_safeguard_ratio
      ! estimated excitation energy used in the curvature-corrected dual step
      REAL(KIND=dp)                                      :: multiplier_curvature_shift
      ! inner tolerance for the next solve = this fraction of the current residual
      ! (never looser than the previous solve, never tighter than EPS_ERROR)
      REAL(KIND=dp)                                      :: multiplier_inner_tol_factor
      INTEGER                                            :: multiplier_max_iter
      ! residual of the current row, max_J |g_IJ|, and the Schur-complement
      ! part of the constraint that belongs to the current state,
      ! -ln(1 - s^T A^-1 s) = ln(Phi_II * (Phi^-1)_II)
      REAL(KIND=dp)                                      :: multiplier_residual
      REAL(KIND=dp)                                      :: schur_constraint
      ! reference copy of the current state's MOs (one per spin) taken at the
      ! last dual update; used to detect a sign (gauge) flip of the orbitals
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: mat_t_gauge_ref
      LOGICAL                                            :: gauge_ref_set = .FALSE.

      ! dipole evaluation
      INTEGER                                            :: dipole_reference_type
      !REAL(kind=dp), DIMENSION(:), POINTER               :: dipole_reference_point

      ! Keep some numerical thresholds
      REAL(KIND=dp)                                      :: eps_filter
      INTEGER                                            :: order_lanczos = 3
      REAL(KIND=dp)                                      :: eps_lanczos = 1.E-3_dp
      INTEGER                                            :: max_iter_lanczos = 128
      INTEGER                                            :: sigma_inv_algorithm = 1 ! Cholesky
      INTEGER                                            :: suggested_mo_block_size
      ! Options for various optimizers
      TYPE(optimizer_options_type)                       :: tidft_optimizer_pcg

      ! CINDY add weight alpha for modified QHQ_S prec
      REAL(KIND=dp)                                      :: weight_alpha_prec

      !!!!!!!!!!!!!!!!!!!!!!!
      !!!!!! MATRICES !!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!

      ! Basis fucntions overlap NxN
      TYPE(dbcsr_type)                                   :: mat_s
      TYPE(dbcsr_type)                                   :: mat_s_inv

      ! MO overlap OxO (stateI, stateJ, spinS)
      TYPE(dbcsr_type), DIMENSION(:,:,:,:), ALLOCATABLE  :: mat_sigma
      TYPE(dbcsr_type), DIMENSION(:,:,:,:), ALLOCATABLE  :: mat_adj_sigma
      ! inverse MO overlap OxO - for each state (stateI, spinS)
      TYPE(dbcsr_type), DIMENSION(:,:,:), ALLOCATABLE    :: mat_sigma_inv
      TYPE(dbcsr_type), DIMENSION(:,:,:), ALLOCATABLE    :: mat_sigma_vv
      TYPE(dbcsr_type), DIMENSION(:), ALLOCATABLE        :: mat_var_sigma
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE         :: det_sig_diagonal
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE         :: det_sig_normalized
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE         :: omega_intra

      ! Overlap of electronic states (stateI, stateJ)
      REAL(KIND=dp), DIMENSION(:,:,:), ALLOCATABLE       :: mat_phi_spin
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE         :: mat_phi
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE         :: mat_phi_inv
      REAL(KIND=dp)                                      :: det_mat_phi
      REAL(KIND=dp)                                      :: max_norm_mat_phi_mI

      ! Pointers to input sections
      TYPE(section_vals_type), POINTER                   :: tidft_print_section

   END TYPE tidft_env_type

CONTAINS

! **************************************************************************************************
!> \brief release tidft envirnoment
!> \brief It is kept in this nonstandard place to avoid circular dependencies
!> \param tidft_env ...
!> \par History
!>       2020.04 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_env_release(tidft_env)

      TYPE(tidft_env_type), POINTER                      :: tidft_env

      !DEALLOCATE (tidft_env%dipole_reference_point)
      DEALLOCATE (tidft_env)

   END SUBROUTINE tidft_env_release

END MODULE tidft_types

