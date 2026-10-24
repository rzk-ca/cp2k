!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Environment for TIDFT methods
!> \par History
!>       2020.04 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE tidft_env_methods

   USE cp_units,                        ONLY: cp_unit_to_cp2k
   USE almo_scf_types,                  ONLY: almo_mat_dim_aobasis,&
                                              almo_mat_dim_occ, &
                                              almo_mat_dim_virt
   USE tidft_types,                     ONLY: tidft_env_type,&
                                              tidft_mult_determinant,&
                                              tidft_mult_spin_channel,&
                                              tidft_mult_orbital_pair
   USE string_utilities,                ONLY: uppercase
   USE tidft_qs,                        ONLY: matrix_layout,&
                                              tidft_init_ks_matrices
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_release,&
                                              mp_para_env_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_release
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_complete_redistribute, dbcsr_copy, dbcsr_create, &
        dbcsr_desymmetrize, dbcsr_distribution_get, dbcsr_distribution_new, &
        dbcsr_distribution_release, dbcsr_distribution_type, dbcsr_filter, dbcsr_finalize, &
        dbcsr_get_block_p, dbcsr_get_info, dbcsr_get_num_blocks, dbcsr_get_stored_coordinates, &
        dbcsr_multiply, dbcsr_nblkcols_total, dbcsr_nblkrows_total, dbcsr_p_type, dbcsr_release, &
        dbcsr_reserve_block2d, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_work_create, &
        dbcsr_type_symmetric, dbcsr_print, dbcsr_add_on_diag
   USE input_constants,                 ONLY: &
        use_mom_ref_coac, use_mom_ref_com, use_mom_ref_user, use_mom_ref_zero, &
        almo_constraint_distance, almo_deloc_none, almo_deloc_xalmo_1diag, &
        almo_domain_layout_atomic, almo_domain_layout_molecular, almo_frz_crystal, &
        almo_mat_distr_molecular, almo_scf_diag, almo_scf_skip, almo_scf_trustr, cg_hager_zhang, &
        do_bondparm_vdw, molecular_guess, tensor_orthogonal, virt_full, virt_minimal, virt_number, &
        xalmo_trial_r0_out, optimizer_pcg, optimizer_trustr
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,set_qs_env
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tidft_env_methods'

   PUBLIC :: tidft_env_create, tidft_init, tidft_env_deallocate
   !, tidft_env_safe_deallocate_state, tidft_env_deallocate_remaining
CONTAINS

! **************************************************************************************************
!> \brief Creation and basic initialization of the tidft environment
!> \param qs_env ...
!> \par History
!>       2020.04 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_env_create(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_env_create', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, nallocate
      TYPE(tidft_env_type), POINTER                      :: tidft_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: input

      CALL timeset(routineN, handle)

      ALLOCATE (tidft_env)
      !ALLOCATE (tidft_env%dipole_reference_point(3))

      ! get basic quantities from the qs_env
      CALL get_qs_env(qs_env, input=input, dft_control=dft_control)

      ! parse the input file and set appropriate quantities
      CALL tidft_read_write_input(input, tidft_env)

      ! Cannot call init here because qs env does not contain necessary data yet
      !CALL tidft_init(qs_env, tidft_env)

      ! let qs_env know about the newly allocalted TIDFT environment
      CALL set_qs_env(qs_env, tidft_env=tidft_env)

      CALL timestop(handle)

   END SUBROUTINE tidft_env_create

! **************************************************************************************************
!> \brief Parses the TIDFT input section
!> \param input ...
!> \param tidft_env ...
!> \par History
!>       2021.09 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_read_write_input(input, tidft_env)
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_read_write_input', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, n_rep, suggested_blk_size, &
                                                            jj, GS_multips, nmultips 
      INTEGER, DIMENSION(:), POINTER                     :: tmplist
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: nstates
      CHARACTER(LEN=80)                                  :: constraint_name
      LOGICAL                                            :: explicit
      TYPE(section_vals_type), POINTER                   :: tidft_pcg_optimizer_section, &
                                                            tidft_section, &
                                                            penalty_section, WF_guess_section

      CALL timeset(routineN, handle)

      tidft_section => section_vals_get_subs_vals(input, "DFT%TIDFT")

      CALL section_vals_val_get(tidft_section, &
              "EPS_FILTER", &
              r_val=tidft_env%eps_filter)
      CALL section_vals_val_get(tidft_section, &
              "SPIN_PURE", &
              l_val=tidft_env%spin_pure)
      CALL section_vals_val_get(tidft_section, &
              "MAX_UNPAIR", &
              i_val=tidft_env%unpairings)
      CALL section_vals_val_get(tidft_section, &
              "SUGGESTED_MO_BLOCK_SIZE", &
              i_val=tidft_env%suggested_mo_block_size)
      CALL section_vals_val_get(tidft_section, &
              "MO_GUESS", &
              i_val=tidft_env%initial_guess_type)
      CALL section_vals_val_get(tidft_section, &
              "MO_GUESS_RANDOM_FRAC", &
              r_val=tidft_env%initial_guess_random_fraction)
      CALL section_vals_val_get(tidft_section, &
              "MO_GUESS_RANDOM_SEED", &
              i_val=tidft_env%initial_guess_random_seed)
      CALL section_vals_val_get(tidft_section, &
              "ENERGY_COMPONENT", &
              l_val=tidft_env%loss_function_components(1))
      CALL section_vals_val_get(tidft_section, &
              "INTRASTATE_COMPONENT", &
              l_val=tidft_env%loss_function_components(2))
      CALL section_vals_val_get(tidft_section, &
              "INTERSTATE_COMPONENT", &
              l_val=tidft_env%loss_function_components(3))
      ! aeph: augmented Lagrangian interstate term toggle (component 4)
      CALL section_vals_val_get(tidft_section, &
              "AUGMENTED_COMPONENT", &
              l_val=tidft_env%loss_function_components(4))
      ! aeph: ALM (component 4) requires the interstate component (3); abort otherwise
      IF (tidft_env%loss_function_components(4) .AND. &
          .NOT. tidft_env%loss_function_components(3)) THEN
         CPABORT("AUGMENTED_COMPONENT requires INTERSTATE_COMPONENT=.TRUE.")
      ENDIF
      ! signed-overlap multiplier term (component 5)
      CALL section_vals_val_get(tidft_section, &
              "MULTIPLIER_COMPONENT", &
              l_val=tidft_env%multiplier_component)
      IF (tidft_env%multiplier_component) THEN
         IF (.NOT. tidft_env%loss_function_components(3)) THEN
            CPABORT("MULTIPLIER_COMPONENT requires INTERSTATE_COMPONENT=.TRUE.")
         ENDIF
         IF (tidft_env%loss_function_components(4)) THEN
            CPABORT("MULTIPLIER_COMPONENT cannot be combined with AUGMENTED_COMPONENT")
         ENDIF
      ENDIF
      CALL section_vals_val_get(tidft_section, &
              "STATE_BY_STATE", &
              l_val=tidft_env%state_by_state)
      IF (tidft_env%multiplier_component .AND. .NOT. tidft_env%state_by_state) THEN
         CPABORT("MULTIPLIER_COMPONENT requires STATE_BY_STATE=.TRUE.")
      ENDIF
      CALL section_vals_val_get(tidft_section, &
              "MULTIPLIER_CONSTRAINT", &
              c_val=constraint_name)
      CALL uppercase(constraint_name)
      SELECT CASE (TRIM(constraint_name))
      CASE ("DETERMINANT")
         tidft_env%multiplier_constraint = tidft_mult_determinant
      CASE ("SPIN_CHANNEL")
         tidft_env%multiplier_constraint = tidft_mult_spin_channel
      CASE ("ORBITAL_PAIR")
         tidft_env%multiplier_constraint = tidft_mult_orbital_pair
      CASE DEFAULT
         CPABORT("MULTIPLIER_CONSTRAINT must be DETERMINANT, SPIN_CHANNEL or ORBITAL_PAIR")
      END SELECT
      tidft_env%multiplier_spin_resolved = &
         (tidft_env%multiplier_constraint .EQ. tidft_mult_spin_channel)
      CALL section_vals_val_get(tidft_section, &
              "PRE_ORTHOGONALIZE_STATES", &
              l_val=tidft_env%pre_orthogonalize_states)
      CALL section_vals_val_get(tidft_section, &
              "NUMERICAL_GRADIENT", &
              l_val=tidft_env%numerical_gradient)
      CALL section_vals_val_get(tidft_section, &
              "NUMERICAL_GRADIENT_DX", &
              r_val=tidft_env%numerical_gradient_dx)
      CALL section_vals_val_get(tidft_section, &
              "N_EXCITED_STATES", &
              i_val=tidft_env%nstates)
      ! add the ground state
      tidft_env%nstates = tidft_env%nstates + 1
      ! Read in the number of excited states of different multiplicities
      !CALL section_vals_val_get(tidft_section,"N_EXCITED_STATES",i_rep_val=1, i_vals=tmplist)
      !nmultips = SIZE(tmplist)
      CALL section_vals_val_get(tidft_section, &
              "DIPOLE_REFERENCE_TYPE", &
              i_val=tidft_env%dipole_reference_type)
      !CALL section_vals_val_get(tidft_section, &
      !        "DIPOLE_REFERENCE_POINT", &
      !        explicit=explicit)
      !IF (explicit) THEN
      !   CALL section_vals_val_get(tidft_section, &
      !           "REFERENCE_POINT", r_vals=tidft_env%dipole_reference_point)
      !   DO jj = 1, 3 
      !      tidft_env%dipole_reference_point = cp_unit_to_cp2k(&
      !         tidft_env%dipole_reference_point(jj), "angstrom")
      !   ENDDO
      !ELSE
      !   NULLIFY (tidft_env%dipole_reference_point)
      !   IF (tidft_env%dipole_reference_type == use_mom_ref_user) THEN
      !      CPABORT("User-defined reference point should be given explicitly")
      !   END IF
      !END IF

      tidft_pcg_optimizer_section => section_vals_get_subs_vals(tidft_section, &
                                                         "TIDFT_OPTIMIZER_PCG")

      ! optimizers
      CALL section_vals_val_get(tidft_pcg_optimizer_section, "EPS_ERROR", &
                                r_val=tidft_env%tidft_optimizer_pcg%eps_error)
      CALL section_vals_val_get(tidft_pcg_optimizer_section, "WHICH_GRADIENT_NORM", &
                                i_val=tidft_env%tidft_optimizer_pcg%which_gnorm)
      CALL section_vals_val_get(tidft_pcg_optimizer_section, "MAX_ITER", &
                                i_val=tidft_env%tidft_optimizer_pcg%max_iter)
      CALL section_vals_val_get(tidft_pcg_optimizer_section, "MAX_ITER_OUTER_LOOP", &
                                i_val=tidft_env%tidft_optimizer_pcg%max_iter_outer_loop)
      CALL section_vals_val_get(tidft_pcg_optimizer_section, "LIN_SEARCH_EPS_ERROR", &
                                r_val=tidft_env%tidft_optimizer_pcg%lin_search_eps_error)
      CALL section_vals_val_get(tidft_pcg_optimizer_section, "LIN_SEARCH_STEP_SIZE_GUESS", &
                                r_val=tidft_env%tidft_optimizer_pcg%lin_search_step_size_guess)
      CALL section_vals_val_get(tidft_pcg_optimizer_section, "CONJUGATOR", &
                                i_val=tidft_env%tidft_optimizer_pcg%conjugator)
      CALL section_vals_val_get(tidft_pcg_optimizer_section, "PRECONDITIONER", &
                                i_val=tidft_env%tidft_optimizer_pcg%preconditioner)
      ! CINDY 
      CALL section_vals_val_get(tidft_pcg_optimizer_section, "PREC_WEIGHT_OVERLAP", &
                                r_val=tidft_env%weight_alpha_prec)

      penalty_section => section_vals_get_subs_vals(tidft_section, "PENALTY")

      CALL section_vals_val_get(penalty_section, &
                                "INTERSTATE_NON_PENALTY_STRENGTH", &
                                r_val=tidft_env%interstate_nonorthogonality_penalty_coeff)
      tidft_env%interstate_nonorthogonality_penalty_coeff = &
         ABS(tidft_env%interstate_nonorthogonality_penalty_coeff)
      !aeph store the input tau for per-pass reset in penalty_strength_variation
      CALL section_vals_val_get(penalty_section, "AUGMENTED_RATIO", &
                          r_val=tidft_env%alm_ratio)
      tidft_env%alm_ratio = ABS(tidft_env%alm_ratio)
      ! aeph: augmented Lagrangian quadratic penalty coefficient (mu)
      CALL section_vals_val_get(penalty_section, &
                                "AUGMENTED_PENALTY_STRENGTH", &
                                r_val=tidft_env%augmented_penalty_coeff)
      tidft_env%augmented_penalty_coeff = &
         ABS(tidft_env%augmented_penalty_coeff)
      ! aeph: store the input mu for per-pass reset in penalty_strength_variation
      tidft_env%augmented_penalty_coeff0 = tidft_env%augmented_penalty_coeff
      CALL section_vals_val_get(penalty_section, &
                                "INTRASTATE_NON_PENALTY_STRENGTH", &
                                r_val=tidft_env%intrastate_penalty_coeff)
      CALL section_vals_val_get(penalty_section, &
                                "PENALTY_STRENGTH_INCREASE_FACTOR", &
                                r_val=tidft_env%penalty_strength_inc_factor)
      CALL section_vals_val_get(penalty_section, &
                                "MAX_NORM_PHI_OFF_DIAG", &
                                r_val=tidft_env%allowed_max_norm_phi_mI)
      ! aeph: read "EPS_DET", the threshold for normalized det(phi)
      CALL section_vals_val_get(penalty_section, &
                                "EPS_DET", &
                                r_val=tidft_env%eps_det)
      ! aeph: same abs as previous penalties                          
      tidft_env%eps_det = ABS(tidft_env%eps_det)

      ! CINDY 
      CALL section_vals_val_get(penalty_section, "NSTEPS_STRENGTH", &
                                i_val=tidft_env%nsteps_penalty_strength)
      CALL section_vals_val_get(penalty_section, &
                                "LOOSE_EPS_SCF", &
                                r_val=tidft_env%loose_eps_scf)
      ! signed-overlap multiplier term settings
      CALL section_vals_val_get(penalty_section, &
                                "MULTIPLIER_SAFEGUARD_RATIO", &
                                r_val=tidft_env%multiplier_safeguard_ratio)
      tidft_env%multiplier_safeguard_ratio = ABS(tidft_env%multiplier_safeguard_ratio)
      CALL section_vals_val_get(penalty_section, &
                                "MULTIPLIER_CURVATURE_SHIFT", &
                                r_val=tidft_env%multiplier_curvature_shift)
      CALL section_vals_val_get(penalty_section, &
                                "MULTIPLIER_INNER_TOL_FACTOR", &
                                r_val=tidft_env%multiplier_inner_tol_factor)
      tidft_env%multiplier_inner_tol_factor = ABS(tidft_env%multiplier_inner_tol_factor)
      CALL section_vals_val_get(penalty_section, &
                                "MULTIPLIER_MAX_ITER", &
                                i_val=tidft_env%multiplier_max_iter)

      ! CINDY: Pass input from WF_GUESS section
      WF_guess_section => section_vals_get_subs_vals(tidft_section, "WF_GUESS")

      CALL section_vals_val_get(WF_guess_section, "EXCITE_ALPHA_ONLY", &
                              l_val=tidft_env%guess_excite_alpha_only)
      CALL section_vals_val_get(WF_guess_section, "MAX_EXC_ORDER", &
                              i_val=tidft_env%guess_max_excitation_order)
      CALL section_vals_val_get(WF_guess_section, "MAX_EXC_ENERGY", &
                              r_val=tidft_env%guess_max_excitation_energy)
      tidft_env%guess_max_excitation_energy = cp_unit_to_cp2k( &
                              tidft_env%guess_max_excitation_energy, "eV")
      CALL section_vals_val_get(WF_guess_section, "MIN_EXC_ENERGY", &
                              r_val=tidft_env%guess_min_excitation_energy)
      tidft_env%guess_min_excitation_energy = cp_unit_to_cp2k( &
                              tidft_env%guess_min_excitation_energy, "eV")

      tidft_env%tidft_print_section => section_vals_get_subs_vals(tidft_section, "PRINT")

      ! check for conflicts between options

      CALL timestop(handle)

   END SUBROUTINE tidft_read_write_input

! **************************************************************************************************
!> \brief Initialization of tidft_env
!> \param qs_env ...
!> \param tidft_env ...
!> \param calc_forces ...
!> \par History
!>       2020.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_init(qs_env, tidft_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_init', routineP = moduleN//':'//routineN

      INTEGER                                            :: ao, handle, i, iao, idomain, ispin, &
                                                            multip, naos, natoms, ndomains, nelec, &
                                                            nelec_a, nelec_b, nmols, nspins, &
                                                            unit_nr, nstates, ivar, nmultips, &
                                                            nslater, slaterI, stateI, unpair
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: input

      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: grad_wrt_penalty_coeff

      CALL timeset(routineN, handle)

      ! define the output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! Set all components to true
      !tidft_env%loss_function_components(:) = .TRUE.

      ! set optimizers' types
      tidft_env%tidft_optimizer_pcg%optimizer_type = optimizer_pcg

      ! get info from the qs_env
      CALL get_qs_env(qs_env, &
                      nelectron_total=tidft_env%nelectrons, &
                      matrix_s=matrix_s, &
                      dft_control=dft_control, &
                      input=input, &
                      has_unit_metric=tidft_env%orthogonal_basis, &
                      para_env=tidft_env%para_env, &
                      blacs_env=tidft_env%blacs_env)
      CALL tidft_env%para_env%retain()
      CALL tidft_env%blacs_env%retain()

      ! copy basic quantities
      IF (dft_control%uks) THEN
         tidft_env%restricted = .FALSE.
         IF (tidft_env%spin_pure) THEN
            CPABORT("Impossible to obtain spin-pure determinants with unrestricted methods")
         ENDIF
      ELSE
         tidft_env%restricted = .TRUE.
      ENDIF

      tidft_env%desired_multip = dft_control%multiplicity

      tidft_env%spin_factor_var = (/ 2, 1, 1/)
      IF (tidft_env%restricted) THEN
         IF (tidft_env%unpairings .EQ. 0) THEN
            IF (tidft_env%desired_multip.EQ.1) THEN
               tidft_env%nspins = 1
            ELSE
               tidft_env%nspins = 2
            ENDIF
         ELSE
               tidft_env%nspins = 2
         ENDIF
      ELSE
         tidft_env%nspins = 2
      ENDIF

      tidft_env%triplet_correction = .FALSE.
      IF (tidft_env%spin_pure) THEN
         IF ( tidft_env%desired_multip.EQ.1 .AND. tidft_env%unpairings.EQ.1 ) THEN
            tidft_env%triplet_correction = .TRUE.
         ELSE IF ( tidft_env%unpairings.EQ.0 ) THEN
            tidft_env%triplet_correction = .FALSE.
         ELSE
            CPABORT("Spin pure states cannot be enforced")
         ENDIF
      ENDIF

      IF (tidft_env%nspins .NE. dft_control%nspins) THEN
         IF (tidft_env%nspins .EQ. 1 .AND. dft_control%nspins .EQ. 2) THEN
            CPABORT("Use RKS, not ROKS, for singlets with 0 unpairings")
         ENDIF
         IF (tidft_env%nspins .EQ. 2 .AND. dft_control%nspins .EQ. 1) THEN
            CPABORT("Use ROKS, not RKS, with unpaired singlets")
         ENDIF
         CPABORT("Number of spins in DFT and TIDFT are different")
      ENDIF
      nspins=tidft_env%nspins

      tidft_env%ndets = 1
      IF (tidft_env%triplet_correction) tidft_env%ndets = 2

      nstates = tidft_env%nstates
      nslater = tidft_env%ndets

      ! Allocate an array that stores the gradient of the objective function
      ! wrt the intrastate penalty strength Cp
      ALLOCATE(tidft_env%grad_wrt_penalty_coeff(nstates))
      ALLOCATE(tidft_env%intrastate_nonorthogonality_penalty_coeff(nstates))
      tidft_env%intrastate_nonorthogonality_penalty_coeff(:) = &
         ABS(tidft_env%intrastate_penalty_coeff)

      ! Number of basis set functions
      CALL dbcsr_get_info(matrix_s(1)%matrix, nfullrows_total=naos)

      ALLOCATE(tidft_env%e_states(nstates))
      ALLOCATE(tidft_env%mat_phi_spin(nstates,nstates,nspins))
      ALLOCATE(tidft_env%mat_phi(nstates,nstates))
      ALLOCATE(tidft_env%mat_phi_inv(nstates,nstates))
      ! signed-overlap multipliers (one per state, used for states J < current I)
      ALLOCATE(tidft_env%overlap_multipliers(nstates))
      tidft_env%overlap_multipliers(:) = 0.0_dp
      ALLOCATE(tidft_env%multiplier_channel(nstates))
      tidft_env%multiplier_channel(:) = 1
      ALLOCATE(tidft_env%pair_r(nstates), tidft_env%pair_P(nstates), &
               tidft_env%pair_c(nstates), tidft_env%pair_q(nstates))
      tidft_env%pair_r(:) = 0.0_dp
      tidft_env%pair_P(:) = 0.0_dp
      tidft_env%pair_c(:) = 0.0_dp
      tidft_env%pair_q(:) = 1.0_dp
      tidft_env%multiplier_residual = 0.0_dp
      tidft_env%schur_constraint = 0.0_dp
      tidft_env%gauge_ref_set = .FALSE.
      IF (tidft_env%multiplier_component .AND. nspins .NE. 2) THEN
         CPABORT("MULTIPLIER_COMPONENT requires two spin channels (nspins=1 is not supported)")
      ENDIF
      ALLOCATE(tidft_env%omega_intra(nstates,nspins))
      ALLOCATE(tidft_env%det_sig_diagonal(nstates,nspins))
      ALLOCATE(tidft_env%det_sig_normalized(nstates,nspins))

      DO stateI = 1, nstates

         tidft_env%e_states(stateI)%nslater = nslater
         ALLOCATE(tidft_env%e_states(stateI)%slater_det(nslater))
         tidft_env%e_states(stateI)%slater_det(1)%desired_multip = tidft_env%desired_multip

         IF (tidft_env%triplet_correction) tidft_env%e_states(stateI)%slater_det(2)%desired_multip = &
            tidft_env%desired_multip + 2

         ! calculate the number of alpha and beta occupied orbitals from
         ! the number of electrons and multiplicity
         ! Na + Nb = Ne
         ! Na - Nb = Mult - 1 (assume Na > Nb as we do not have more info from get_molecule_set_info)
         ! Na = (Ne + Mult - 1)/2
         ! Nb = (Ne - Mult + 1)/2
         DO slaterI = 1, nslater

            ! set the number of alpha/beta electrons
            tidft_env%e_states(stateI)%slater_det(slaterI)%nelectrons_ab(1) = &
               ( tidft_env%nelectrons + tidft_env%e_states(stateI)%slater_det(slaterI)%desired_multip - 1 ) / 2
            tidft_env%e_states(stateI)%slater_det(slaterI)%nelectrons_ab(2) = &
               ( tidft_env%nelectrons - tidft_env%e_states(stateI)%slater_det(slaterI)%desired_multip + 1 ) / 2

            ! set the number of alpha/beta virtual orbitals
            tidft_env%e_states(stateI)%slater_det(slaterI)%nvirtuals_ab(1) = &
               naos - tidft_env%e_states(stateI)%slater_det(slaterI)%nelectrons_ab(1)
            tidft_env%e_states(stateI)%slater_det(slaterI)%nvirtuals_ab(2) = &
               naos - tidft_env%e_states(stateI)%slater_det(slaterI)%nelectrons_ab(2)
            IF ((tidft_env%e_states(stateI)%slater_det(slaterI)%nvirtuals_ab(1) + &
                 tidft_env%e_states(stateI)%slater_det(slaterI)%nvirtuals_ab(2)) .EQ. 0) THEN
               CPABORT("There are no virtual states")
            ENDIF

            ! set the number of independent variable matrices
            tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(:)=0
            IF (tidft_env%restricted) THEN
               ! 1st entry: doubly occupied orbitals: (ne-(mult-1))/2-unpair
               IF (tidft_env%triplet_correction .AND. slaterI.EQ.2) THEN
                  unpair = 0
               ELSE
                  unpair = tidft_env%unpairings
               ENDIF
               ivar = (tidft_env%nelectrons - tidft_env%e_states(stateI)%slater_det(slaterI)%desired_multip + 1) / 2 - &
                  unpair
               IF (ivar .LT. 0) ivar = 0
               tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(1) = ivar
               ! 2nd entry: singly occupied alpha orbitals: mult-1+unpair
               ivar = tidft_env%e_states(stateI)%slater_det(slaterI)%desired_multip - 1 + unpair
               IF (ivar .LT. 0) ivar = 0
               tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(2) = ivar
               ! 3rd entry: singly occupied beta orbitals: unpair
               ivar = unpair
               IF (ivar .LT. 0) ivar = 0
               tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(3) = ivar
               !IF (tidft_env%triplet_correction) tidft_env%e_states(stateI)%slater_det(2)%desired_multip = &
            ELSE
               ! There are no electron pairs in UKS
               tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(1) = 0
               tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(2) = &
                  tidft_env%e_states(stateI)%slater_det(slaterI)%nelectrons_ab(1)
               tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(3) = &
                  tidft_env%e_states(stateI)%slater_det(slaterI)%nelectrons_ab(2)
            ENDIF

            IF ( (tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(1) + &
                  tidft_env%e_states(stateI)%slater_det(slaterI)%ncols_var(2)) .EQ. 0) THEN
               CPABORT("There are no alpha electrons")
            ENDIF

         ENDDO ! slaterI
      ENDDO ! stateI

      CALL tidft_env_create_matrices(tidft_env, matrix_s(1)%matrix)

      CALL tidft_init_ao_overlap(tidft_env, matrix_s(1)%matrix)
      CALL tidft_init_ks_matrices(qs_env, tidft_env)

      !tiidft      CALL construct_qs_mos(qs_env, tidft_env)

      ! FINISH setting job parameters here, print out job info
      CALL tidft_print_job_info(tidft_env, unit_nr)

      CALL timestop(handle)

   END SUBROUTINE tidft_init

! **************************************************************************************************
!> \brief create various matrices
!> \param tidft_env ...
!> \param matrix_s0 ...
!> \par History
!>       2021.09 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_env_create_matrices(tidft_env, matrix_s0)

      TYPE(tidft_env_type), INTENT(INOUT)             :: tidft_env
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s0

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_env_create_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, &
                                                            nstates, suggested_blk_size, &
                                                            idet, naos, &
                                                            nspins, ivar, &
                                                            stateI, ndets, &
                                                            nmos_stateI, &
                                                            nmos_stateI_flip, &
                                                            stateJ, n_left, &
                                                            nmos_stateJ

      CALL timeset(routineN, handle)

      ! RZK-WARNING: It seems that making S symmetric does not break the code
      ! but could save time in some of the operations
      suggested_blk_size = tidft_env%suggested_mo_block_size
      CALL matrix_layout( &
              matrix_ao_template = matrix_s0, &
              matrix_new = tidft_env%mat_s, &
              name_new = 'm_s', &
              suggested_mo_blk_size = suggested_blk_size, & ! does not affect AOxAO matrices
              orbital_type_new = (/almo_mat_dim_aobasis, almo_mat_dim_aobasis/), &
              symmetry_new = dbcsr_type_no_symmetry)
      CALL dbcsr_get_info(tidft_env%mat_s, nfullrows_total=naos)

      CALL dbcsr_create(tidft_env%mat_s_inv, template=tidft_env%mat_s)

      nstates = tidft_env%nstates
      nspins = tidft_env%nspins
      ndets = tidft_env%ndets

      ALLOCATE(tidft_env%mat_sigma(nstates,nstates,ndets,nspins))
      ALLOCATE(tidft_env%mat_adj_sigma(nstates,nstates,ndets,nspins))
      ALLOCATE(tidft_env%mat_sigma_inv(nstates,ndets,nspins))
      ALLOCATE(tidft_env%mat_var_sigma(3))
      ALLOCATE(tidft_env%mat_sigma_vv(nstates,ndets,nspins))
      DO idet = 1, ndets
         DO stateI = 1, nstates
            ALLOCATE(tidft_env%e_states(stateI)%slater_det(idet)%mat_t(nspins))
            ALLOCATE(tidft_env%e_states(stateI)%slater_det(idet)%mat_ks(nspins))
            ALLOCATE(tidft_env%e_states(stateI)%slater_det(idet)%mat_v(nspins))
            ALLOCATE(tidft_env%e_states(stateI)%slater_det(idet)%mat_vo(nspins))
            !ALLOCATE(tidft_env%e_states(stateI)%slater_det(idet)%mat_proj(nspins))
            ! It is not necessary to allocate mat_var for the other determinant(s)
            IF (idet.EQ.1) ALLOCATE(tidft_env%e_states(stateI)%slater_det(1)%mat_var(3))
         ENDDO
      ENDDO

      DO idet = 1, ndets

         DO ispin = 1, nspins

            DO stateI = 1, nstates

               nmos_stateI = tidft_env%e_states(stateI)%slater_det(idet)%nelectrons_ab(ispin)

               IF ( (naos - nmos_stateI) .GT. 0) THEN ! non-zero virtuals

                  CALL matrix_layout( &
                          matrix_ao_template = matrix_s0, &
                          matrix_new = tidft_env%e_states(stateI)%slater_det(idet)%mat_v(ispin), &
                          name_new = 'm_V_Ia', &
                          orbital_type_new = (/almo_mat_dim_aobasis, almo_mat_dim_virt/), &
                          suggested_mo_blk_size = suggested_blk_size, &
                          symmetry_new = dbcsr_type_no_symmetry, &
                          nmos = (/-1, naos - nmos_stateI/) )

                  CALL matrix_layout( &
                          matrix_ao_template = matrix_s0, &
                          matrix_new = tidft_env%e_states(stateI)%slater_det(idet)%mat_vo(ispin), &
                          name_new = 'm_VO_Ia', &
                          orbital_type_new = (/almo_mat_dim_virt, almo_mat_dim_occ/), &
                          suggested_mo_blk_size = suggested_blk_size, &
                          symmetry_new = dbcsr_type_no_symmetry, &
                          nmos = (/naos - nmos_stateI, nmos_stateI/) )

                  CALL matrix_layout( &
                          matrix_ao_template = matrix_s0, &
                          matrix_new = tidft_env%mat_sigma_vv(stateI,idet,ispin), &
                          name_new = 'm_VV_Ia', &
                          orbital_type_new = (/almo_mat_dim_virt, almo_mat_dim_virt/), &
                          suggested_mo_blk_size = suggested_blk_size, &
                          symmetry_new = dbcsr_type_no_symmetry, &
                          nmos = (/naos - nmos_stateI, naos - nmos_stateI/) )

               ENDIF ! non-zero virtuals

               IF ( nmos_stateI .GT. 0) THEN ! non-zero occupied

                  CALL matrix_layout( &
                          matrix_ao_template = matrix_s0, &
                          matrix_new = tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin), &
                          name_new = 'm_T_Ia', &
                          orbital_type_new = (/almo_mat_dim_aobasis, almo_mat_dim_occ/), &
                          suggested_mo_blk_size = suggested_blk_size, &
                          symmetry_new = dbcsr_type_no_symmetry, &
                          nmos = (/-1, nmos_stateI/) )

                  CALL dbcsr_create(tidft_env%e_states(stateI)%slater_det(idet)%mat_ks(ispin), template=tidft_env%mat_s)

                  DO stateJ = stateI, nstates

                     nmos_stateJ = tidft_env%e_states(stateJ)%slater_det(idet)%nelectrons_ab(ispin)

                     IF ( nmos_stateJ .GT. 0) THEN

                        CALL matrix_layout( &
                                matrix_ao_template = matrix_s0, &
                                matrix_new = tidft_env%mat_sigma(stateI, stateJ, idet, ispin) , &
                                name_new = 'm_sig_IJab', &
                                orbital_type_new = (/almo_mat_dim_occ, almo_mat_dim_occ/), &
                                suggested_mo_blk_size = suggested_blk_size, &
                                symmetry_new = dbcsr_type_no_symmetry, &
                                nmos = (/nmos_stateI, nmos_stateJ/) )

                        CALL dbcsr_create(tidft_env%mat_adj_sigma(stateI, stateJ, idet, ispin), &
                                template=tidft_env%mat_sigma(stateI, stateJ, idet, ispin))

                     ENDIF ! zero MOs for state J

                  ENDDO ! stateJ

                  CALL dbcsr_create(tidft_env%mat_sigma_inv(stateI, idet, ispin), &
                          template=tidft_env%mat_sigma(stateI, stateI, idet, ispin))

               ENDIF ! zero MOs

            ENDDO ! stateI

         ENDDO ! ispin

      ENDDO ! idet

      ! create non-zero matrices of independent variables
      ! NOTE: only the first determinant contains the variables
      DO stateI = 1, nstates

         DO ivar = 1, 3

            nmos_stateI = tidft_env%e_states(stateI)%slater_det(1)%ncols_var(ivar)

            IF ( nmos_stateI .GT. 0) THEN

               CALL matrix_layout( &
                       matrix_ao_template = matrix_s0, &
                       matrix_new = tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar) , &
                       name_new = 'm_X_Ia', &
                       orbital_type_new = (/almo_mat_dim_aobasis, almo_mat_dim_occ/), &
                       suggested_mo_blk_size = suggested_blk_size, &
                       symmetry_new = dbcsr_type_no_symmetry, &
                       nmos = (/-1, nmos_stateI/) )

               ! Create template, one state is sufficient
               IF (stateI .EQ. 1) THEN

                  CALL matrix_layout( &
                          matrix_ao_template = matrix_s0, &
                          matrix_new = tidft_env%mat_var_sigma(ivar) , &
                          name_new = 'm_var_sigma', &
                          orbital_type_new = (/almo_mat_dim_occ, almo_mat_dim_occ/), &
                          suggested_mo_blk_size = suggested_blk_size, &
                          symmetry_new = dbcsr_type_no_symmetry, &
                          nmos = (/nmos_stateI, nmos_stateI/) ) ! these are not MOs but independent vars

               ENDIF

            ENDIF ! non-zero ncols_var

         ENDDO ! ivar

      ENDDO ! stateI

      CALL timestop(handle)

   END SUBROUTINE tidft_env_create_matrices

! **************************************************************************************************
!> \brief deallocate tidft envirnoment data
!> \param tidft_env ...
!> \par History
!>       2020.04 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_env_deallocate(tidft_env)

      TYPE(tidft_env_type), INTENT(INOUT)             :: tidft_env

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_env_deallocate', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: handle, ispin, &
                                                            nstates, idet, &
                                                            nspins, ivar, &
                                                            stateI, ndets, &
                                                            nmos_stateI, &
                                                            nmos_stateI_flip, &
                                                            stateJ, &
                                                            nmos_stateJ

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ndets = tidft_env%ndets
      nstates = tidft_env%nstates
      nspins = tidft_env%nspins

      CALL dbcsr_release(tidft_env%mat_s)
      CALL dbcsr_release(tidft_env%mat_s_inv)

      DO idet = 1, ndets

         DO ispin = 1, nspins

            DO stateI = 1, nstates

               nmos_stateI = tidft_env%e_states(stateI)%slater_det(idet)%nelectrons_ab(ispin)

               IF ( nmos_stateI .GT. 0) THEN

                  CALL dbcsr_release(tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin))
                  CALL dbcsr_release(tidft_env%e_states(stateI)%slater_det(idet)%mat_ks(ispin))
                  !CALL dbcsr_release(tidft_env%e_states(stateI)%slater_det(idet)%mat_proj(ispin))

                  DO stateJ = stateI, nstates

                     nmos_stateJ = tidft_env%e_states(stateJ)%slater_det(idet)%nelectrons_ab(ispin)

                     IF ( nmos_stateJ .GT. 0) THEN

                         CALL dbcsr_release(tidft_env%mat_sigma(stateI, stateJ, idet, ispin))
                         CALL dbcsr_release(tidft_env%mat_adj_sigma(stateI, stateJ, idet, ispin))

                     ENDIF ! zero MOs for J

                  ENDDO ! stateJ

                  CALL dbcsr_release(tidft_env%mat_sigma_inv(stateI, idet, ispin))

               ENDIF ! zero MOs

            ENDDO ! stateI

         ENDDO ! ispin

      ENDDO ! idet

      DO stateI = 1, nstates

         DO ivar = 1, 3

            nmos_stateI = tidft_env%e_states(stateI)%slater_det(1)%ncols_var(ivar)

            IF ( nmos_stateI .GT. 0) THEN

               CALL dbcsr_release(tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar))

               IF (stateI .EQ. 1) THEN
                  CALL dbcsr_release(tidft_env%mat_var_sigma(ivar))
               ENDIF

            ENDIF

         ENDDO ! ivar

      ENDDO ! stateI

      DEALLOCATE(tidft_env%mat_sigma)
      DEALLOCATE(tidft_env%mat_adj_sigma)
      DEALLOCATE(tidft_env%mat_sigma_inv)
      DEALLOCATE(tidft_env%mat_var_sigma)
      DEALLOCATE(tidft_env%mat_sigma_vv)
      DO idet = 1, ndets
         DO stateI = 1, nstates
            DEALLOCATE(tidft_env%e_states(stateI)%slater_det(idet)%mat_t)
            DEALLOCATE(tidft_env%e_states(stateI)%slater_det(idet)%mat_ks)
            DEALLOCATE(tidft_env%e_states(stateI)%slater_det(idet)%mat_v)
            DEALLOCATE(tidft_env%e_states(stateI)%slater_det(idet)%mat_vo)
            !DEALLOCATE(tidft_env%e_states(stateI)%slater_det(idet)%mat_proj)
            IF (idet.EQ.1) DEALLOCATE(tidft_env%e_states(stateI)%slater_det(1)%mat_var)
         ENDDO
      ENDDO

      DEALLOCATE(tidft_env%intrastate_nonorthogonality_penalty_coeff)
      DEALLOCATE(tidft_env%grad_wrt_penalty_coeff)
      DEALLOCATE(tidft_env%mat_phi_inv)
      DEALLOCATE(tidft_env%mat_phi)
      IF (ALLOCATED(tidft_env%overlap_multipliers)) DEALLOCATE(tidft_env%overlap_multipliers)
      IF (ALLOCATED(tidft_env%multiplier_channel)) DEALLOCATE(tidft_env%multiplier_channel)
      IF (ALLOCATED(tidft_env%pair_r)) DEALLOCATE(tidft_env%pair_r)
      IF (ALLOCATED(tidft_env%pair_P)) DEALLOCATE(tidft_env%pair_P)
      IF (ALLOCATED(tidft_env%pair_c)) DEALLOCATE(tidft_env%pair_c)
      IF (ALLOCATED(tidft_env%pair_q)) DEALLOCATE(tidft_env%pair_q)
      IF (ALLOCATED(tidft_env%pair_u)) DEALLOCATE(tidft_env%pair_u)
      IF (ALLOCATED(tidft_env%pair_v)) DEALLOCATE(tidft_env%pair_v)
      DEALLOCATE(tidft_env%mat_phi_spin)
      DEALLOCATE(tidft_env%e_states)
      DEALLOCATE(tidft_env%omega_intra)
      DEALLOCATE(tidft_env%det_sig_diagonal)
      DEALLOCATE(tidft_env%det_sig_normalized)

      CALL mp_para_env_release(tidft_env%para_env)
      CALL cp_blacs_env_release(tidft_env%blacs_env)

      CALL timestop(handle)

   END SUBROUTINE tidft_env_deallocate

! **************************************************************************************************
!> \brief deallocate memory for a selective state I after optimization of that state is finished
!> \param tidft_env ...
!> \par History
!>       2023.08 created [Cindy Pham]
!> \author Cindy Pham
! **************************************************************************************************
   SUBROUTINE tidft_env_safe_deallocate_state(tidft_env, state_I)

      TYPE(tidft_env_type), INTENT(INOUT)             :: tidft_env
      INTEGER, INTENT(IN)                             :: state_I

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_env_safe_deallocate_state', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: handle, ispin, &
                                                            nstates, idet, &
                                                            nspins, ivar, &
                                                            state, ndets, &
                                                            nmos_stateI, &
                                                            nmos_stateI_flip, &
                                                            stateI, stateJ, &
                                                            nmos_stateJ

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ndets = tidft_env%ndets
      nstates = tidft_env%nstates  ! nstates=tidft_env%nstates=stateI (from tidft_variatonal_optimization)
      nspins = tidft_env%nspins

      DO idet = 1, ndets

         DO ispin = 1, nspins

            DO stateI = state_I, nstates

               nmos_stateI = tidft_env%e_states(stateI)%slater_det(idet)%nelectrons_ab(ispin)

               IF ( nmos_stateI .GT. 0) THEN

                  CALL dbcsr_release(tidft_env%e_states(stateI)%slater_det(idet)%mat_ks(ispin))

                  DO stateJ = stateI, nstates

                     nmos_stateJ = tidft_env%e_states(stateJ)%slater_det(idet)%nelectrons_ab(ispin)

                     IF ( nmos_stateJ .GT. 0) THEN

                     !   CALL dbcsr_release(tidft_env%mat_sigma(stateI, stateJ, idet, ispin))
                        CALL dbcsr_release(tidft_env%mat_adj_sigma(stateI, stateJ, idet, ispin))

                     ENDIF ! zero MOs for J

                  ENDDO ! stateJ

                  !CALL dbcsr_release(tidft_env%mat_sigma_inv(stateI, idet, ispin))

               ENDIF ! zero MOs

            ENDDO ! stateI

         ENDDO ! ispin

      ENDDO ! idet

      DO stateI = state_I, nstates

         DO ivar = 1, 3

            nmos_stateI = tidft_env%e_states(stateI)%slater_det(1)%ncols_var(ivar)

            IF ( nmos_stateI .GT. 0) THEN
               IF (stateI .EQ. 1) THEN
                  CALL dbcsr_release(tidft_env%mat_var_sigma(ivar))

               ENDIF

            ENDIF

         ENDDO ! ivar

      ENDDO ! stateI

      CALL timestop(handle)

   END SUBROUTINE tidft_env_safe_deallocate_state

! **************************************************************************************************
!> \brief deallocate the remaining after calling the selective state cleanup
!> \param tidft_env ...
!> \par History
!>       2023.08 created [Cindy Pham]
!> \author Cindy Pham
! **************************************************************************************************
   SUBROUTINE tidft_env_deallocate_remaining(tidft_env)

      TYPE(tidft_env_type), INTENT(INOUT)             :: tidft_env

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_env_deallocate', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: handle, ispin, &
                                                            nstates, idet, &
                                                            nspins, ivar, &
                                                            stateI, ndets, &
                                                            nmos_stateI, &
                                                            nmos_stateI_flip, &
                                                            stateJ, &
                                                            nmos_stateJ

      CALL timeset(routineN, handle)

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ndets = tidft_env%ndets
      nstates = tidft_env%nstates
      nspins = tidft_env%nspins

      CALL dbcsr_release(tidft_env%mat_s)
      CALL dbcsr_release(tidft_env%mat_s_inv)

      DO idet = 1, ndets

         DO ispin = 1, nspins

            DO stateI = 1, nstates

               nmos_stateI = tidft_env%e_states(stateI)%slater_det(idet)%nelectrons_ab(ispin)

               IF ( nmos_stateI .GT. 0) THEN

                  CALL dbcsr_release(tidft_env%e_states(stateI)%slater_det(idet)%mat_t(ispin))

                  CALL dbcsr_release(tidft_env%mat_sigma(stateI, stateJ, idet, ispin))

                  CALL dbcsr_release(tidft_env%mat_sigma_inv(stateI, idet, ispin))

               ENDIF ! zero MOs

            ENDDO ! stateI

         ENDDO ! ispin

      ENDDO ! idet

      DO stateI = 1, nstates

         DO ivar = 1, 3

            nmos_stateI = tidft_env%e_states(stateI)%slater_det(1)%ncols_var(ivar)

            IF ( nmos_stateI .GT. 0) THEN

               CALL dbcsr_release(tidft_env%e_states(stateI)%slater_det(1)%mat_var(ivar))

               IF(unit_nr > 0) THEN
                  WRITE (unit_nr,'(T2,A49,I11)') "Released 3 mainvar of state", stateI
               ENDIF

            ENDIF

         ENDDO ! ivar

      ENDDO ! stateI

      !DEALLOCATE(tidft_env%mat_sigma) ! comment out because the job stopped at this step with mem error
      !DEALLOCATE(tidft_env%mat_adj_sigma)
      !DEALLOCATE(tidft_env%mat_sigma_inv)
      DEALLOCATE(tidft_env%mat_var_sigma)
      DEALLOCATE(tidft_env%mat_sigma_vv)

      DO idet = 1, ndets
         DO stateI = 1, nstates
            DEALLOCATE(tidft_env%e_states(stateI)%slater_det(idet)%mat_t)
            DEALLOCATE(tidft_env%e_states(stateI)%slater_det(idet)%mat_ks)
            DEALLOCATE(tidft_env%e_states(stateI)%slater_det(idet)%mat_v)
            DEALLOCATE(tidft_env%e_states(stateI)%slater_det(idet)%mat_vo)
            IF (idet.EQ.1) DEALLOCATE(tidft_env%e_states(stateI)%slater_det(1)%mat_var)
            IF(unit_nr > 0) THEN
               WRITE (unit_nr,'(T2,A49,I11)') "Released mat_KS, mat_vm, mat_vo of state", stateI
            ENDIF
         ENDDO
      ENDDO

      IF(unit_nr > 0) WRITE (*,*) "Start deallocating other matrices"
      DEALLOCATE(tidft_env%intrastate_nonorthogonality_penalty_coeff)
      DEALLOCATE(tidft_env%grad_wrt_penalty_coeff)
      DEALLOCATE(tidft_env%mat_phi_inv)
      DEALLOCATE(tidft_env%mat_phi)
      IF (ALLOCATED(tidft_env%overlap_multipliers)) DEALLOCATE(tidft_env%overlap_multipliers)
      IF (ALLOCATED(tidft_env%multiplier_channel)) DEALLOCATE(tidft_env%multiplier_channel)
      IF (ALLOCATED(tidft_env%pair_r)) DEALLOCATE(tidft_env%pair_r)
      IF (ALLOCATED(tidft_env%pair_P)) DEALLOCATE(tidft_env%pair_P)
      IF (ALLOCATED(tidft_env%pair_c)) DEALLOCATE(tidft_env%pair_c)
      IF (ALLOCATED(tidft_env%pair_q)) DEALLOCATE(tidft_env%pair_q)
      IF (ALLOCATED(tidft_env%pair_u)) DEALLOCATE(tidft_env%pair_u)
      IF (ALLOCATED(tidft_env%pair_v)) DEALLOCATE(tidft_env%pair_v)
      DEALLOCATE(tidft_env%mat_phi_spin)
      DEALLOCATE(tidft_env%e_states)

      CALL cp_para_env_release(tidft_env%para_env)
      CALL cp_blacs_env_release(tidft_env%blacs_env)

      CALL timestop(handle)

   END SUBROUTINE tidft_env_deallocate_remaining

! **************************************************************************************************
!> \brief Initializes the TIDFT SCF copy of the AO overlap matrix
!> \param matrix_s0 ...
!> \param tidft_env ...
!> \par History
!>       2021.10 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE tidft_init_ao_overlap(tidft_env, matrix_s0)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_s0
      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_init_ao_overlap', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, unit_nr
      TYPE(cp_logger_type), POINTER                      :: logger

      CALL timeset(routineN, handle)

      IF (tidft_env%orthogonal_basis) THEN
         CALL dbcsr_set(tidft_env%mat_s, 0.0_dp)
         CALL dbcsr_add_on_diag(tidft_env%mat_s, 1.0_dp)
      ELSE
         CALL dbcsr_copy(tidft_env%mat_s, matrix_s0)
      ENDIF

      CALL dbcsr_filter(tidft_env%mat_s, tidft_env%eps_filter)

      ! invert S using cholesky for dense matrices
      CALL dbcsr_copy(tidft_env%mat_s_inv, tidft_env%mat_s)
      CALL cp_dbcsr_cholesky_decompose(tidft_env%mat_s_inv, &
              para_env=tidft_env%para_env, &
              blacs_env=tidft_env%blacs_env)
      CALL cp_dbcsr_cholesky_invert(tidft_env%mat_s_inv, &
              para_env=tidft_env%para_env, &
              blacs_env=tidft_env%blacs_env, &
              upper_to_full=.TRUE.)
      CALL dbcsr_filter(tidft_env%mat_s_inv, &
              tidft_env%eps_filter)

      CALL timestop(handle)

   END SUBROUTINE tidft_init_ao_overlap

! **************************************************************************************************
!> \brief Prints out a short summary about the TIDFT job
!> \param tidft_env ...
!> \param unit_nr ...
!> \par History
!>       2021.10 created [Cindy Pham, Rustam Z. Khaliullin]
!> \author Cindy Pham
! **************************************************************************************************
   SUBROUTINE tidft_print_job_info(tidft_env, unit_nr)

      TYPE(tidft_env_type), INTENT(INOUT)             :: tidft_env
      INTEGER, INTENT(IN)                             :: unit_nr

      CHARACTER(len=*), PARAMETER :: routineN = 'tidft_print_job_info', &
         routineP = moduleN//':'//routineN

      INTEGER                                         :: handle, ispin, stateI, stateJ, ivar, slaterI, nslater, nstates

      CALL timeset(routineN, handle)

      IF (unit_nr > 0) THEN

         WRITE (unit_nr, '()')
         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
         WRITE (unit_nr, '(T2,A,T19,A,T78,A)') REPEAT("-", 3), &
            "Welcome to time-independent DFT (TI-DFT) -", REPEAT("-", 3)
         WRITE (unit_nr, '(T2,A,T16,A,T78,A)') REPEAT("-", 3), &
            "an awesome variational method for excited states", REPEAT("-", 3)
         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)

         ! RKS, ROKS or UKS
         IF (tidft_env%restricted) THEN
            IF (tidft_env%nspins.EQ.2) THEN
               WRITE (unit_nr, '(T2,A)') "Restricted open-shell KS method is used"
            ELSE
               WRITE (unit_nr, '(T2,A)') "Restricted closed-shell KS method is used"
            ENDIF
         ELSE
            WRITE (unit_nr, '(T2,A)') "Unrestricted KS method is used"
         ENDIF

         ! Spin pure
         IF (tidft_env%spin_pure) THEN
            WRITE (unit_nr, '(T2,A)') "States are enforced to be spin pure"
            IF (tidft_env%triplet_correction) THEN
               WRITE (unit_nr, '(T4,A)') "Triplet states are optimized together with spin-flipped singlet states"
            ENDIF
         ELSE
            WRITE (unit_nr, '(T2,A)') "Spin purity is not enforced"
         ENDIF

         ! nspins
         SELECT CASE (tidft_env%nspins)
         CASE (1)
            WRITE (unit_nr, '(T2,A)') "Alpha and beta electrons are described using one set of MOs"
         CASE (2)
            WRITE (unit_nr, '(T2,A)') "Alpha and beta electrons are described using two sets of MOs"
         END SELECT

         ! nstates = GS + ES
         WRITE (unit_nr, '(T2,A,T70,I11)') "No. of electronic states", tidft_env%nstates

         ! deisred_multip: Na - Nb + 1
         WRITE (unit_nr, '(T2,A,T70,I11)') "Desired multiplicity", tidft_env%desired_multip

         ! unpairings
         WRITE (unit_nr, '(T2,A,T70,I11)') "Max no. of allowed unpairings", tidft_env%unpairings

         WRITE (unit_nr, '()')
         WRITE (unit_nr, '(T2,A,T70,I11)') "No. of electrons", tidft_env%nelectrons
         DO stateI = 1, 1 ! tidft_env%nstates ! no need to run over all states - these vars are the same
            DO slaterI = 1, tidft_env%e_states(stateI)%nslater
               IF (slaterI .EQ. 2) WRITE (unit_nr, '(T2,A,T70,I11)') "Triplet state"
               WRITE (unit_nr, '(T4,A,T70,I11)') "No. of alpha electrons", &
                  tidft_env%e_states(stateI)%slater_det(slaterI)%nelectrons_ab(1)
               WRITE (unit_nr, '(T4,A,T70,I11)') "No. of beta electrons", &
                  tidft_env%e_states(stateI)%slater_det(slaterI)%nelectrons_ab(2)
               WRITE (unit_nr, '(T4,A,T70,I11)') "No. of doubly occupied MOs", &
                  tidft_env%e_states(stateI)%slater_det(1)%ncols_var(1)
               WRITE (unit_nr, '(T4,A,T70,I11)') "No. of singly occupied alpha MOs", &
                  tidft_env%e_states(stateI)%slater_det(1)%ncols_var(2)+(slaterI-1)
               WRITE (unit_nr, '(T4,A,T70,I11)') "No. of singly occupied beta MOs", &
                  tidft_env%e_states(stateI)%slater_det(1)%ncols_var(3)-(slaterI-1)
            ENDDO
         ENDDO

         ! suggested_blk_size
         WRITE (unit_nr, '()')
         WRITE (unit_nr, '(T2,A,T70,I11)') "Suggested block size for MOs", tidft_env%suggested_mo_block_size

         WRITE (unit_nr, '(T2,A)') REPEAT("-", 79)
         WRITE (unit_nr, '()')

      ENDIF ! unit_nr > 0

      CALL timestop(handle)

   END SUBROUTINE tidft_print_job_info

END MODULE tidft_env_methods
