!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Variable-metric TIDFT: methods to generate initial guess
!> \par History
!>       2022.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
MODULE tidft_guess_methods
   USE tidft_guess_types,               ONLY: tidft_excitation_list_type, &
                                              tidft_excitation_type, &
                                              orbital_set_type
   USE almo_scf_methods,                ONLY: almo_scf_t_to_proj,&
                                              get_overlap,&
                                              orthogonalize_mos
   USE almo_scf_types,                  ONLY: optimizer_options_type
   USE cp_dbcsr_cholesky,               ONLY: cp_dbcsr_cholesky_decompose,&
                                              cp_dbcsr_cholesky_invert
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_syevd
   USE cp_dbcsr_api,                    ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_dot, dbcsr_filter, dbcsr_multiply, dbcsr_norm, &
        dbcsr_norm_maxabsnorm, dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, &
        dbcsr_type_no_symmetry, dbcsr_distribution_type, dbcsr_get_info, dbcsr_get_diag, &
        dbcsr_transposed, dbcsr_print, dbcsr_set_diag, dbcsr_dot, dbcsr_desymmetrize, &
        dbcsr_type_symmetric, dbcsr_init_random
   USE input_constants,                 ONLY: &
        cg_dai_yuan, cg_fletcher, cg_fletcher_reeves, &
        cg_hager_zhang, cg_hestenes_stiefel, cg_liu_storey, cg_polak_ribiere, cg_zero
   USE iterate_matrix,                  ONLY: determinant
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE machine,                         ONLY: m_flush,&
                                              m_walltime
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE tidft_types,                     ONLY: tidft_env_type
   USE tidft_optimizer_types,           ONLY: tidft_optimizer_storage_type
   USE tidft_methods,                   ONLY: tidft_optimizer_mainvar_dot, &
                                              tidft_optimizer_mainvar_add, &
                                              tidft_optimizer_mainvar_multiply, &
                                              tidft_optimizer_mainvar_multiply_diff, &
                                              tidft_optimizer_mainvar_copy, &
                                              tidft_optimizer_gradient_norm, &
                                              tidft_optimizer_mainvar_print, &
                                              tidft_optimizer_mainvar_set, &
                                              tidft_mainvar_orthonormalize, &
                                              dbcsr_set_element, &
                                              copy_columns, &
                                              stitch_columns, &
                                              split_columns, &
                                              determinant_local, &
                                              matrix_invert_local, &
                                              gradVAR_scale_split_increment, &
                                              precMO_to_precVAR_increment, &
                                              swap_columns
   USE physcon,                         ONLY: evolt
   USE tidft_qs,                        ONLY: tidft_dm_to_tidft_ks
   USE util,                            ONLY: sort
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'tidft_guess_methods'

   PUBLIC :: construct_canonical_orbitals, &
             guess_excited_states_orbitals

   LOGICAL, PARAMETER :: debug_mode = .TRUE.
   LOGICAL, PARAMETER :: safe_mode = .FALSE.

CONTAINS

! **************************************************************************************************
!> \brief Constructs initial guess for molecular oribatals of excited states
!> \param tidft_env ...
!> \par History
!>       2022.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE guess_excited_states_orbitals(tidft_env)

      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env

      CHARACTER(len=*), PARAMETER :: routineN = 'guess_excited_states_orbitals', &
         routineP = moduleN//':'//routineN

      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(tidft_excitation_list_type)                   :: excitation_list
      CHARACTER(LEN=30)                                  :: spin_string
      INTEGER                                            :: unit_nr, nspins, &
                                                            ispin, n, ivirtual, &
                                                            idet, norbs, iorb, &
                                                            count_unpairings, &
                                                            nalpha, nbeta
      INTEGER                                            :: nexcitations, &
                                                            iexcitation, &
                                                            generated_nstates, &
                                                            istate, state_index, &
                                                            current_order, estate, &
                                                            current_spin, &
                                                            nmos_stateI, &
                                                            nmos_stateV
      INTEGER                                            :: max_order, spin_excitation, iorder
      INTEGER, DIMENSION(:,:), ALLOCATABLE               :: start_orbital_indices
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: index_true, orbital_counters
      LOGICAL                                            :: check_unpairings
      REAL(KIND=dp)                                      :: eps_diff_orb_energy
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: delta, energies_temp
      TYPE(orbital_set_type)                             :: orbital_set
      TYPE(tidft_excitation_type)                        :: current_excitation
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: alph_indx_occ
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: alph_indx_vir
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: beta_indx_occ
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: beta_indx_vir

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ! RZK TODO:
      !   * Create an input section within &TIDFT called &WF_GUESS
      !     It should contain controls for this procedure:
      !       MAX_EXC_ENERGY - maximum allowed excitation energy for the guess (pre-SCF)
      !       MIN_EXC_ENERGY - minimum allowed excitation energy for the guess (pre-SCF)
      !       MAX_EXC_ORDER  - maximum allowed number of excitations from the ground state
      !     The corresponding variables for these input keywords are:
      !       tidft_env%guess_max_excitation_energy (to be created)
      !       tidft_env%guess_min_excitation_energy (already implemented)
      !       tidft_env%guess_max_excitation_order (to be created)
      !     More controls will be added later
      !   * Test this procedure for UKS and RKS for He crystal (fast, login node) to see whether
      !     excitation pairs are correctly generated. That is, the test run should only print
      !     occupied-virtual pairs generated for each excited state guess. It is not necessary to
      !     create the guess orbitals yet.
      !   * NExt test is to do the same for ROKS singlet and ROKS triplet (with triplet correction FALSE).
      !   * Next text: same system ROKS singlet with triplet correction TRUE.
      !   * Repeat the test for another simple system (fast, login node). F3B..NH3 with small cell can be used.

      ! MAJOR STEP: Obtain canonical orbital energies and info on orbital location in matrices
      nspins = tidft_env%nspins
      idet = 1 ! Use the "mixed" state to produce canonical orbitals

      ! The first index is spin, the second is occ/virt
      ALLOCATE(orbital_set%orbital_list(nspins, 2))
      DO ispin = 1, nspins

         DO ivirtual = 1, 2

            IF (ivirtual .EQ. 1) THEN
               CALL dbcsr_get_info(tidft_env%mat_sigma_inv(1,idet,ispin), nfullrows_total=norbs)
            ELSE
               CALL dbcsr_get_info(tidft_env%mat_sigma_vv(1,idet,ispin), nfullrows_total=norbs)
            ENDIF

            orbital_set%orbital_list(ispin, ivirtual)%norbs = norbs
            orbital_set%orbital_list(ispin, ivirtual)%ispin = ispin
            orbital_set%orbital_list(ispin, ivirtual)%ivirtual = ivirtual - 1

            IF (norbs .GT. 0) THEN
               ALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%orb_energy(norbs))
               ALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%index_within_matrix(norbs))
               ALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%block_within_matrix(norbs))
               ALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%index_within_block(norbs))
               ALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%orb_has_spin_double(norbs))
               IF (nspins .EQ. 1) THEN ! RKS for all states
                  orbital_set%orbital_list(ispin, ivirtual)%orb_has_spin_double(:) = .TRUE.
               ELSE ! UKS and ROKS
                  orbital_set%orbital_list(ispin, ivirtual)%orb_has_spin_double(:) = .FALSE.
               ENDIF
            ENDIF

         ENDDO

      ENDDO

      estate = 1
      CALL construct_canonical_orbitals(tidft_env, estate = estate, &
              orbital_set = orbital_set)

      ! Obtain info on the position of orbitals within the matrices
      DO ispin = 1, nspins

         nmos_stateI = tidft_env%e_states(estate)%slater_det(idet)%nelectrons_ab(ispin)
         nmos_stateV = tidft_env%e_states(estate)%slater_det(idet)%nvirtuals_ab(ispin)

         IF ( nmos_stateI .GT. 0) THEN

            ! Obtain info on the position of orbitals within the
            ! occupied-orbital matrix
            CALL get_orbital_indices(tidft_env%e_states(estate)%slater_det(idet)%mat_t(ispin), &
                    orbital_set%orbital_list(ispin, 1)%block_within_matrix, &
                    orbital_set%orbital_list(ispin, 1)%index_within_block)

         ENDIF ! non-zero occupied

         IF ( nmos_stateV .GT. 0) THEN

            ! Obtain info on the position of orbitals within the matrix
            CALL get_orbital_indices(tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), &
                    orbital_set%orbital_list(ispin, 2)%block_within_matrix, &
                    orbital_set%orbital_list(ispin, 2)%index_within_block)

         ENDIF ! non-zero virtuals

      ENDDO

      ! Update main variables so they are consistent with occupied canonical orbitals
      CALL MO2Var( &
              mat_t = tidft_env%e_states(estate)%slater_det(idet)%mat_t, &
              mat_var = tidft_env%e_states(estate)%slater_det(idet)%mat_var, &
              ncols_var = tidft_env%e_states(estate)%slater_det(idet)%ncols_var, &
              mat_s = tidft_env%mat_s, &
              mat_oo_template = tidft_env%mat_sigma(estate,estate,idet,:), &
              eps_filter = tidft_env%eps_filter)

      ! MAJOR STEP: Find out what orbitals are doubly and singly occupied
      ! For virtual orbitals, find out whether they CAN contain two or one electron(s)
      ! RZK-WARNING: This primitive algorithm only compares energies of alpha-beta orbitals
      ! In the future, norms of Delta T = T_alpha - T_beta should be analyzed

      ! This procedure is irrelevant for RKS and UKS because enforcing
      ! the correct number of unpairings is not required.
      ! Perform only for the case of ROKS
      count_unpairings = -1
      check_unpairings = tidft_env%restricted .AND. tidft_env%unpairings.GT.0
      IF (check_unpairings) THEN

         eps_diff_orb_energy = SQRT(tidft_env%eps_filter)

         count_unpairings = 0
         DO ivirtual = 1, 2

            nalpha = orbital_set%orbital_list(1, ivirtual)%norbs
            nbeta = orbital_set%orbital_list(2, ivirtual)%norbs
            ! RZK-WARNING: Implement cases with unequal number of electrons later
            CPASSERT(nalpha.EQ.nbeta)

            ALLOCATE(delta(nalpha))
            delta(:) = ABS(orbital_set%orbital_list(1, ivirtual)%orb_energy(:) - &
                       orbital_set%orbital_list(2, ivirtual)%orb_energy(:))
            DO iorb = 1, nalpha
               IF ( delta(iorb) .LT. eps_diff_orb_energy ) THEN
                  orbital_set%orbital_list(1, ivirtual)%orb_has_spin_double(iorb) = .TRUE.
                  orbital_set%orbital_list(2, ivirtual)%orb_has_spin_double(iorb) = .TRUE.
               ELSE
                  orbital_set%orbital_list(1, ivirtual)%orb_has_spin_double(iorb) = .FALSE.
                  orbital_set%orbital_list(2, ivirtual)%orb_has_spin_double(iorb) = .FALSE.
                  IF (ivirtual .EQ. 1) count_unpairings = count_unpairings + 1
               ENDIF
            ENDDO ! iorb
            DEALLOCATE(delta)

         ENDDO ! ivirtual

         IF (count_unpairings .GT. tidft_env%unpairings) THEN
            CPABORT("WF analysis found too many unpaired electrons")
         ENDIF

         ! RZK-WARNING: for triple corrction number of unpairings must be exctly equal to 1
         ! Here and in excitation generation

      ENDIF ! ROKS

      ! MAJOR STEP: Select occupied and virtual orbital from the set using various conditions
      ! RZK-WARNING: NYI

      ! MAJOR STEP: Loop over select occupied and virtual orbitals, keep desired in the list
      ALLOCATE(start_orbital_indices(nspins, 2))
      DO ispin = 1, nspins
         ! occupied are looped from top to bottom
         start_orbital_indices(ispin, 1) = &
            orbital_set%orbital_list(ispin, 1)%norbs
         ! virtuals are looped from bottom to top
         start_orbital_indices(ispin, 2) = 1
      ENDDO ! ispin

      ! Init storage structures
      nexcitations = tidft_env%nstates - 1
      max_order = tidft_env%guess_max_excitation_order
      CPASSERT(nexcitations.GT.0)
      CPASSERT(max_order.GT.0)
      excitation_list%max_nexcitations = nexcitations
      excitation_list%max_order = max_order
      excitation_list%cur_nexcitations = 0
      ALLOCATE(excitation_list%excitations(nexcitations))

      DO iexcitation = 1, nexcitations
         excitation_list%excitations(iexcitation)%order = -1
         excitation_list%excitations(iexcitation)%excitation_energy = 0.0_dp
         ALLOCATE(excitation_list%excitations(iexcitation)%delta_epsilon(max_order))
         ALLOCATE(excitation_list%excitations(iexcitation)%occupied(max_order))
         ALLOCATE(excitation_list%excitations(iexcitation)%virtual(max_order))
         ALLOCATE(excitation_list%excitations(iexcitation)%spin(max_order))
         excitation_list%excitations(iexcitation)%delta_epsilon(:) = 0.0_dp
         excitation_list%excitations(iexcitation)%occupied(:) = -1
         excitation_list%excitations(iexcitation)%virtual(:) = -1
         excitation_list%excitations(iexcitation)%spin(:) = -1
      ENDDO

      ALLOCATE(current_excitation%delta_epsilon(max_order))
      ALLOCATE(current_excitation%occupied(max_order))
      ALLOCATE(current_excitation%virtual(max_order))
      ALLOCATE(current_excitation%spin(max_order))
      current_excitation%order = 0
      current_excitation%excitation_energy = 0.0_dp
      current_excitation%unpairings_count = count_unpairings
      current_excitation%delta_epsilon(:) = 0.0_dp
      current_excitation%occupied(:) = -1
      current_excitation%virtual(:) = -1
      current_excitation%spin(:) = -1

      ! recursive subroutine for generaing multiple single excitations
      CALL generate_single_excitations( &
              current_order = 1, &
              excitation_list = excitation_list, &
              current_excitation = current_excitation, &
              orbital_set = orbital_set, &
              start_orbital_indices = start_orbital_indices, &
              min_excitation_energy = tidft_env%guess_min_excitation_energy, &
              check_unpairings = check_unpairings, &
              max_allowed_unpairings = tidft_env%unpairings, &
              excite_alpha_only = tidft_env%guess_excite_alpha_only)

      DEALLOCATE(current_excitation%delta_epsilon)
      DEALLOCATE(current_excitation%occupied)
      DEALLOCATE(current_excitation%virtual)
      DEALLOCATE(current_excitation%spin)

      ! MAJOR STEP: Order the list of excitations
      generated_nstates = excitation_list%cur_nexcitations
      ALLOCATE(index_true(generated_nstates))
      ALLOCATE(energies_temp(generated_nstates))

      DO istate = 1, generated_nstates
         energies_temp(istate) = excitation_list%excitations(istate)%excitation_energy
      ENDDO
      CALL sort(energies_temp, generated_nstates, index_true)

      ! RZK-WARNING: Change tidft_env%nstates
      !              if fewer excitations generated
      !              than requested

      ! CINDY: printout info for each excitation
      IF (unit_nr > 0) &
         WRITE (unit_nr,'(T2,A)') REPEAT("-", 78)

      DO istate = 1, generated_nstates

         state_index = index_true(istate)

         IF (unit_nr > 0) &
            WRITE(unit_nr,'(T3, A10, I5, A13, F8.3, A5)') &
               "EXCITATION", istate, "Guess energy", &
               evolt*excitation_list%excitations(state_index)%excitation_energy, &
               "(eV)"

         ! Arrays of appropriate lengths to feed into copy_columns
         current_order = excitation_list%excitations(state_index)%order

         ALLOCATE(orbital_counters(nspins))
         orbital_counters(:) = 0

         DO iorder=1, current_order

            ! determines if the spin is alpha or beta
            current_spin = excitation_list%excitations(state_index)%spin(iorder)

            IF (current_spin .EQ. 1) THEN ! alpha e-
               spin_string = "(alpha)"
               orbital_counters(1) = orbital_counters(1) + 1
            ELSE
               spin_string = "(beta)"
               orbital_counters(2) = orbital_counters(2) + 1
            ENDIF

            IF (unit_nr > 0) &
               WRITE(unit_nr,'(T2, A24, A, I3.3, A8, A12, A, I3.3, A8)') "HOMO", "-", &
                  orbital_set%orbital_list(current_spin, 1)%norbs - &
                     excitation_list%excitations(state_index)%occupied(iorder), &
                  spin_string, "--->  LUMO", "+", &
                  excitation_list%excitations(state_index)%virtual(iorder) - 1, &
                  spin_string

         ENDDO ! iorder

         IF (unit_nr > 0) &
            WRITE (unit_nr,'(T2,A)') REPEAT("-", 78)

         ! Now we know the number of alpha and beta transitions
         IF ( orbital_counters(1) .GT. 0) THEN
            ALLOCATE(alph_indx_occ(orbital_counters(1)))
            ALLOCATE(alph_indx_vir(orbital_counters(1)))
         ENDIF
         IF ( orbital_counters(2) .GT. 0) THEN
            ALLOCATE(beta_indx_occ(orbital_counters(2)))
            ALLOCATE(beta_indx_vir(orbital_counters(2)))
         ENDIF

         orbital_counters(:) = 0

         DO iorder=1, current_order

            ! determines if the spin is alpha or beta
            current_spin = excitation_list%excitations(state_index)%spin(iorder)

            IF (current_spin .EQ. 1) THEN ! alpha e-

               orbital_counters(1) = orbital_counters(1) + 1

               alph_indx_occ( orbital_counters(1) ) = &
                  excitation_list%excitations(state_index)%occupied(iorder)
               alph_indx_vir( orbital_counters(1) ) = &
                  excitation_list%excitations(state_index)%virtual(iorder)

            ELSE

               orbital_counters(2) = orbital_counters(2) + 1

               beta_indx_occ( orbital_counters(2) ) = &
                  excitation_list%excitations(state_index)%occupied(iorder)
               beta_indx_vir( orbital_counters(2) ) = &
                  excitation_list%excitations(state_index)%virtual(iorder)

            ENDIF

         ENDDO ! iorder

         DO ispin = 1, nspins

            CALL dbcsr_copy( &
                    tidft_env%e_states(istate+1)%slater_det(idet)%mat_t(ispin), &
                    tidft_env%e_states(1)%slater_det(idet)%mat_t(ispin))

            IF ( orbital_counters(ispin) .GT. 0) THEN

               !WRITE(*,*) "MAT_V spin:", ispin
               !CALL dbcsr_print(tidft_env%e_states(1)%slater_det(idet)%mat_v(ispin))

               !WRITE(*,*) " Old MAT_T:", ispin
               !CALL dbcsr_print(tidft_env%e_states(istate+1)%slater_det(idet)%mat_t(ispin))

               IF (ispin .EQ. 1) THEN
                  CALL copy_columns( &
                          matrix_dest = tidft_env%e_states(istate+1)%slater_det(idet)%mat_t(ispin), &
                          matrix_source = tidft_env%e_states(1)%slater_det(idet)%mat_v(ispin), &
                          cols_dest = alph_indx_occ, &
                          cols_source = alph_indx_vir)
               ELSE ! beta
                  CALL copy_columns( &
                          matrix_dest = tidft_env%e_states(istate+1)%slater_det(idet)%mat_t(ispin), &
                          matrix_source = tidft_env%e_states(1)%slater_det(idet)%mat_v(ispin), &
                          cols_dest = beta_indx_occ, &
                          cols_source = beta_indx_vir)
               ENDIF

               !WRITE(*,*) " New MAT_T:", ispin
               !CALL dbcsr_print(tidft_env%e_states(istate+1)%slater_det(idet)%mat_t(ispin))

            ENDIF

         ENDDO ! ispin

         DEALLOCATE(orbital_counters)

         IF ( ALLOCATED(alph_indx_occ) ) THEN
            DEALLOCATE(alph_indx_occ)
            DEALLOCATE(alph_indx_vir)
         ENDIF
         IF ( ALLOCATED(beta_indx_occ) ) THEN
            DEALLOCATE(beta_indx_occ)
            DEALLOCATE(beta_indx_vir)
         ENDIF

         CALL MO2Var( &
                 mat_t = tidft_env%e_states(istate+1)%slater_det(idet)%mat_t, &
                 mat_var = tidft_env%e_states(istate+1)%slater_det(idet)%mat_var, &
                 ncols_var = tidft_env%e_states(istate+1)%slater_det(idet)%ncols_var, &
                 mat_s = tidft_env%mat_s, &
                 mat_oo_template = tidft_env%mat_sigma(istate+1,istate+1,idet,:), &
                 eps_filter = tidft_env%eps_filter)

         !DO iorder=1, 3
         !   IF (tidft_env%e_states(istate+1)%slater_det(idet)%ncols_var(iorder).GT.0) THEN
         !      WRITE(*,*) "VARIABLE", iorder
         !      CALL dbcsr_print(tidft_env%e_states(istate+1)%slater_det(idet)%mat_var(iorder))
         !   ELSE
         !      WRITE(*,*) "NO INDEPENDENT VARIABLE: ", iorder
         !   ENDIF
         !ENDDO

      ENDDO ! istate

      CALL store_guess_occupations(tidft_env, excitation_list, index_true, generated_nstates)

      ! MAJOR STEP: Clean up
      DEALLOCATE(index_true)
      DEALLOCATE(energies_temp)

      DO iexcitation = 1, nexcitations
         DEALLOCATE(excitation_list%excitations(iexcitation)%delta_epsilon)
         DEALLOCATE(excitation_list%excitations(iexcitation)%occupied)
         DEALLOCATE(excitation_list%excitations(iexcitation)%virtual)
         DEALLOCATE(excitation_list%excitations(iexcitation)%spin)
      ENDDO

      DEALLOCATE(excitation_list%excitations)

      DO ispin = 1, nspins

         DO ivirtual = 1, 2

            !IF (unit_nr > 0) THEN
            !   WRITE(unit_nr, *) "Orbital energies for spin, virt: ", ispin, ivirtual
            !   WRITE(unit_nr, *) orbital_set%orbital_list(ispin, ivirtual)%orb_energy(:)
            !ENDIF

            IF ( orbital_set%orbital_list(ispin, ivirtual)%norbs .GT. 0 ) THEN
               DEALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%orb_energy)
               DEALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%index_within_matrix)
               DEALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%block_within_matrix)
               DEALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%index_within_block)
               DEALLOCATE(orbital_set%orbital_list(ispin, ivirtual)%orb_has_spin_double)
            ENDIF

         ENDDO

      ENDDO
      DEALLOCATE(start_orbital_indices)
      DEALLOCATE(orbital_set%orbital_list)

   END SUBROUTINE guess_excited_states_orbitals

   RECURSIVE SUBROUTINE generate_single_excitations(current_order, excitation_list, &
      current_excitation, orbital_set, start_orbital_indices, min_excitation_energy, &
      check_unpairings, max_allowed_unpairings, excite_alpha_only)

      INTEGER, INTENT(IN)                                :: current_order
      TYPE(tidft_excitation_list_type), INTENT(INOUT)    :: excitation_list
      TYPE(tidft_excitation_type), INTENT(IN)            :: current_excitation
      TYPE(orbital_set_type), INTENT(IN)                 :: orbital_set
      INTEGER, DIMENSION(:,:), INTENT(IN)                :: start_orbital_indices
      REAL(KIND=dp), INTENT(IN)                          :: min_excitation_energy
      LOGICAL, INTENT(IN)                                :: check_unpairings
      INTEGER, INTENT(IN)                                :: max_allowed_unpairings

      INTEGER                                            :: nspins, max_order, iocc, &
                                                            ispin, opposite_spin, &
                                                            nvirtuals, avirt, &
                                                            incr_unpairings, &
                                                            iexcitation, &
                                                            insertion_index
      INTEGER, DIMENSION(:,:), ALLOCATABLE               :: next_level_start_indices
      TYPE(tidft_excitation_type)                        :: next_level_excitation
      REAL(KIND=dp)                                      :: spin_factor, max_energy
      LOGICAL                                            :: add_to_excitation_list, found, &
                                                            originates_from_doubly_occupied, &
                                                            goes_to_doubly_occupied_virtual, &
                                                            max_order_condition, &
                                                            deeper_level
      ! CINDY: add keyword to excite alpha electron only
      LOGICAL, INTENT(IN)                                ::  excite_alpha_only

      nspins = SIZE(orbital_set%orbital_list, 1)
      max_order = excitation_list%max_order

      ALLOCATE(next_level_start_indices(nspins, 2))

      ALLOCATE(next_level_excitation%delta_epsilon(max_order))
      ALLOCATE(next_level_excitation%occupied(max_order))
      ALLOCATE(next_level_excitation%virtual(max_order))
      ALLOCATE(next_level_excitation%spin(max_order))

      ! copy indices
      next_level_start_indices(:,:) = start_orbital_indices(:,:)

      IF (nspins .EQ. 2) THEN
         spin_factor = 1.0_dp
      ELSE
         spin_factor = 2.0_dp
      ENDIF

      DO ispin = 1, nspins

         IF (excite_alpha_only .AND. ispin .EQ. 2) THEN
            CYCLE
         ENDIF
         opposite_spin = 1
         IF (nspins .EQ. 2) THEN
            IF (ispin .EQ. 1) THEN
               opposite_spin = 2
            ELSE
               opposite_spin = 1
            ENDIF
         ENDIF

         next_level_start_indices(opposite_spin, :) = &
            start_orbital_indices(opposite_spin, :)

         DO iocc = start_orbital_indices(ispin, 1), 1, -1

            nvirtuals = orbital_set%orbital_list(ispin, 2)%norbs

            DO avirt = start_orbital_indices(ispin, 2), nvirtuals, +1

               ! copy current excitation into the next level excitation
               CALL copy_excitation(dest = next_level_excitation, &
                       source = current_excitation)

               ! update the next level excitation
               next_level_excitation%order = current_order
               next_level_excitation%spin(current_order) = ispin
               next_level_excitation%occupied(current_order) = iocc
               next_level_excitation%virtual(current_order) = avirt
               next_level_excitation%delta_epsilon(current_order) = &
                  spin_factor * ( &
                  orbital_set%orbital_list(ispin, 2)%orb_energy(avirt) - &
                  orbital_set%orbital_list(ispin, 1)%orb_energy(iocc) )
               next_level_excitation%excitation_energy = &
                  next_level_excitation%excitation_energy + &
                  next_level_excitation%delta_epsilon(current_order)

               ! copy the newly created excitation to the list if it satisfies:
               !   energy filter requirememnts (delE should be more than the min excitation energy, no max)
               !   unpairings no. requirement (ROKS)
               !   RZK_WARNING: NYI - spin symmetry requirement (in some cases alpha excitations are equivalent to beta)
               add_to_excitation_list = .TRUE.

               ! energy filter requirement
               IF (SUM(next_level_excitation%delta_epsilon) .LT. min_excitation_energy) THEN
                  add_to_excitation_list = .FALSE.
               ENDIF

               ! No need to check if the previous ADD criteria failed
               IF (check_unpairings .AND. add_to_excitation_list) THEN

                  incr_unpairings = 0

                  ! Check if the excitation is from a doubly-occupied orbital
                  originates_from_doubly_occupied = .FALSE.
                  IF (orbital_set%orbital_list(ispin, 1)%orb_has_spin_double(iocc)) THEN
                     ! This is an occupied orbital that can hold two electrons
                     ! Is there another electron - of the opposite spin - in the orbital?
                     ! Check if the other electron is already excited
                     found = .FALSE.
                     DO iexcitation = 1, current_order - 1
                        IF (current_excitation%occupied(iexcitation) .EQ. iocc .AND. &
                            current_excitation%spin(iexcitation) .EQ. opposite_spin ) &
                            found = .TRUE.
                     ENDDO
                     IF (.NOT. found) originates_from_doubly_occupied = .TRUE.
                  ENDIF
                  IF (originates_from_doubly_occupied) &
                     incr_unpairings = incr_unpairings + 1

                  ! Check if the excitation is to a doubly-occupied orbital
                  goes_to_doubly_occupied_virtual = .FALSE.
                  IF (orbital_set%orbital_list(ispin, 2)%orb_has_spin_double(avirt)) THEN
                     ! This is a virtual orbital that can hold two electrons
                     ! Is there another electron - of the opposite spin - already in the orbital?
                     found = .FALSE.
                     DO iexcitation = 1, current_order - 1
                        IF (current_excitation%virtual(iexcitation) .EQ. avirt .AND. &
                            current_excitation%spin(iexcitation) .EQ. opposite_spin ) &
                            found = .TRUE.
                     ENDDO
                     IF (found) goes_to_doubly_occupied_virtual = .TRUE.
                  ENDIF
                  IF (goes_to_doubly_occupied_virtual) &
                     incr_unpairings = incr_unpairings - 1

                  next_level_excitation%unpairings_count = &
                     next_level_excitation%unpairings_count + incr_unpairings

                  IF (next_level_excitation%unpairings_count .GT. max_allowed_unpairings) &
                     add_to_excitation_list = .FALSE.

               ENDIF ! check unpairings

               IF (add_to_excitation_list) THEN

                  insertion_index = -1

                  ! If there are empty slots add to the next empty slot
                  IF (excitation_list%cur_nexcitations .LT. excitation_list%max_nexcitations) THEN
                     insertion_index = excitation_list%cur_nexcitations + 1
                     excitation_list%cur_nexcitations = excitation_list%cur_nexcitations + 1
                  ELSE
                  ! Find the location with the maximum excitation energy
                     max_energy = 0.0_dp
                     DO iexcitation = 1, excitation_list%cur_nexcitations
                        IF (excitation_list%excitations(iexcitation)%excitation_energy .GT. &
                           max_energy) THEN
                           insertion_index = iexcitation
                           max_energy = excitation_list%excitations(iexcitation)%excitation_energy
                        ENDIF
                     ENDDO ! excitation loop
                     IF (next_level_excitation%excitation_energy .GT. max_energy) &
                        insertion_index = -1
                  ENDIF

                  ! Copy excitation to the list
                  IF (insertion_index .GT. 0) THEN
                     CALL copy_excitation(dest = excitation_list%excitations(insertion_index), &
                             source = next_level_excitation)
                  ENDIF

               ENDIF

               ! prepare for the next level excitation
               next_level_start_indices(ispin, 1) = iocc - 1
               next_level_start_indices(ispin, 2) = avirt + 1

               max_order_condition = current_order .LT. max_order
               deeper_level = max_order_condition

               IF (deeper_level) THEN

                  CALL generate_single_excitations( &
                          current_order = current_order + 1, &
                          excitation_list = excitation_list, &
                          current_excitation = next_level_excitation, &
                          orbital_set = orbital_set, &
                          start_orbital_indices = next_level_start_indices, &
                          min_excitation_energy = min_excitation_energy, &
                          check_unpairings = check_unpairings, &
                          max_allowed_unpairings = max_allowed_unpairings,&
                          excite_alpha_only = excite_alpha_only)

               ENDIF

            ENDDO ! a virt

         ENDDO ! iocc

      ENDDO ! ispin

      DEALLOCATE(next_level_excitation%delta_epsilon)
      DEALLOCATE(next_level_excitation%occupied)
      DEALLOCATE(next_level_excitation%virtual)
      DEALLOCATE(next_level_excitation%spin)

      DEALLOCATE(next_level_start_indices)

   END SUBROUTINE generate_single_excitations

! **************************************************************************************************
!> \brief Constructs canonical orbitals. They are stored in mat_t of the environment
!> \param tidft_env ...
!> \param estate - electronic state for which the orbitals are constructed
!> \param orbital_set - storage for orbital energies
!> \par History
!>       2022.06 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE construct_canonical_orbitals(tidft_env, estate, orbital_set)

      TYPE(tidft_env_type), INTENT(INOUT)             :: tidft_env
      INTEGER, INTENT(IN)                             :: estate
      TYPE(orbital_set_type), INTENT(INOUT)           :: orbital_set

      CHARACTER(len=*), PARAMETER :: routineN = 'construct_canonical_orbitals', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: unit_nr, ispin, &
                                                            idet, ivirtual, &
                                                            ndets, nspins, &
                                                            nmos_stateI, &
                                                            nmos_stateV, &
                                                            ndiag
      TYPE(dbcsr_type)                                   :: tempNV1, tempVOcc1, &
                                                            tempVOcc2, tempVV1, &
                                                            tempVV2, tempNO1, &
                                                            tempOO1, tempOO2
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: fake_e_per_block
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: diagonal
      TYPE(cp_logger_type), POINTER                      :: logger

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      ALLOCATE(fake_e_per_block(1))
      fake_e_per_block(:) = -1

      nspins = tidft_env%nspins
      ndets=1
      IF (tidft_env%triplet_correction) ndets = 2
      idet = 1 ! Use the singlet state to produce canonical orbitals

      DO ispin = 1, nspins

         nmos_stateI = tidft_env%e_states(estate)%slater_det(idet)%nelectrons_ab(ispin)

         IF ( nmos_stateI .GT. 0) THEN

            CALL dbcsr_create(tempNO1, &
                              template=tidft_env%e_states(estate)%slater_det(idet)%mat_t(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_create(tempOO1, &
                              template=tidft_env%mat_sigma(estate,estate,idet,ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_create(tempOO2, &
                              template=tidft_env%mat_sigma(estate,estate,idet,ispin), &
                              matrix_type=dbcsr_type_no_symmetry)

            CALL orthogonalize_mos(ket=tidft_env%e_states(estate)%slater_det(idet)%mat_t(ispin), &
                                   overlap=tidft_env%mat_sigma(estate,estate,idet,ispin), &
                                   metric=tidft_env%mat_s, &
                                   retain_locality=.FALSE., &
                                   only_normalize=.FALSE., &
                                   nocc_of_domain=fake_e_per_block, &
                                   eps_filter=tidft_env%eps_filter, &
                                   order_lanczos=tidft_env%order_lanczos, &
                                   eps_lanczos=tidft_env%eps_lanczos, &
                                   max_iter_lanczos=tidft_env%max_iter_lanczos)

            ! compute the OxO block of the KS matrix
            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_ks(ispin), &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_t(ispin), &
                                0.0_dp, tempNO1, &
                                filter_eps=tidft_env%eps_filter)

            CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_t(ispin), &
                                tempNO1, &
                                0.0_dp, tempOO1, &
                                filter_eps=tidft_env%eps_filter)

            CALL cp_dbcsr_syevd(tempOO1, tempOO2, &
                                orbital_set%orbital_list(ispin, 1)%orb_energy, &
                                !o_eigenvalues(:,ispin), &
                                para_env=tidft_env%para_env, &
                                blacs_env=tidft_env%blacs_env)

            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_t(ispin), &
                                tempOO2, &
                                0.0_dp, tempNO1, &
                                filter_eps=tidft_env%eps_filter)

            CALL dbcsr_copy(tidft_env%e_states(estate)%slater_det(idet)%mat_t(ispin), tempNO1)

            ! update occ overlap
            CALL dbcsr_get_info(tidft_env%mat_sigma(estate,estate,idet,ispin), nfullrows_total=ndiag)
            ALLOCATE (diagonal(ndiag))
            diagonal(:) = 1.0_dp
            CALL dbcsr_set(tidft_env%mat_sigma(estate,estate,idet,ispin), 0.0_dp)
            CALL dbcsr_set_diag(tidft_env%mat_sigma(estate,estate,idet,ispin), diagonal)
            CALL dbcsr_filter(tidft_env%mat_sigma(estate,estate,idet,ispin), tidft_env%eps_filter)
            ! update occ inverse overlap
            CALL dbcsr_set(tidft_env%mat_sigma_inv(estate,idet,ispin), 0.0_dp)
            CALL dbcsr_set_diag(tidft_env%mat_sigma_inv(estate,idet,ispin), diagonal)
            CALL dbcsr_filter(tidft_env%mat_sigma_inv(estate,idet,ispin), tidft_env%eps_filter)
            DEALLOCATE (diagonal)

            CALL dbcsr_release(tempNO1)
            CALL dbcsr_release(tempOO1)
            CALL dbcsr_release(tempOO2)

         ENDIF ! non-zero occupied

         nmos_stateV = tidft_env%e_states(estate)%slater_det(idet)%nvirtuals_ab(ispin)

         IF ( nmos_stateV .GT. 0) THEN

            ! Virtual orbitals
            CALL dbcsr_create(tempNV1, &
                              template=tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_create(tempVV1, &
                              template=tidft_env%mat_sigma_vv(estate,idet,ispin), &
                              matrix_type=dbcsr_type_no_symmetry)
            CALL dbcsr_create(tempVV2, &
                              template=tidft_env%mat_sigma_vv(estate,idet,ispin), &
                              matrix_type=dbcsr_type_no_symmetry)

            ! If there are virtual orbitals, generate random virtual matrix
            CALL dbcsr_init_random(tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), &
                    mini_seed=estate+ispin*2+tidft_env%initial_guess_random_seed)

            ! Project the occupied subspace out if there are any
            ! occupied orbitals to project out
            IF ( nmos_stateI .GT. 0) THEN

               CALL dbcsr_create(tempVOcc1, &
                                 template=tidft_env%e_states(estate)%slater_det(idet)%mat_vo(ispin), &
                                 matrix_type=dbcsr_type_no_symmetry)
               CALL dbcsr_create(tempVOcc2, &
                                 template=tidft_env%e_states(estate)%slater_det(idet)%mat_vo(ispin), &
                                 matrix_type=dbcsr_type_no_symmetry)

               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   tidft_env%mat_s, &
                                   tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), &
                                   0.0_dp, tempNV1, &
                                   filter_eps=tidft_env%eps_filter)

               CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                   tempNV1, &
                                   tidft_env%e_states(estate)%slater_det(idet)%mat_t(ispin), &
                                   0.0_dp, tempVOcc1, &
                                   filter_eps=tidft_env%eps_filter)

               CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                   tempVOcc1, &
                                   tidft_env%mat_sigma_inv(estate,idet,ispin), &
                                   0.0_dp, tempVOcc2, &
                                   filter_eps=tidft_env%eps_filter)

               CALL dbcsr_multiply("N", "T", 1.0_dp, &
                                   tidft_env%e_states(estate)%slater_det(idet)%mat_t(ispin), &
                                   tempVOcc2, &
                                   0.0_dp, tempNV1, &
                                   filter_eps=tidft_env%eps_filter)

               CALL dbcsr_add(tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), tempNV1, 1.0_dp, -1.0_dp)

               CALL dbcsr_release(tempVOcc1)
               CALL dbcsr_release(tempVOcc2)

            ENDIF ! project out occupied subspace

            CALL orthogonalize_mos(ket=tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), &
                                   overlap=tempVV1, &
                                   metric=tidft_env%mat_s, &
                                   retain_locality=.FALSE., &
                                   only_normalize=.FALSE., &
                                   nocc_of_domain=fake_e_per_block, &
                                   eps_filter=tidft_env%eps_filter, &
                                   order_lanczos=tidft_env%order_lanczos, &
                                   eps_lanczos=tidft_env%eps_lanczos, &
                                   max_iter_lanczos=tidft_env%max_iter_lanczos)

            ! compute VxV block of the KS matrix
            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_ks(ispin), &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), &
                                0.0_dp, tempNV1, &
                                filter_eps=tidft_env%eps_filter)

            CALL dbcsr_multiply("T", "N", 1.0_dp, &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), &
                                tempNV1, &
                                0.0_dp, tempVV1, &
                                filter_eps=tidft_env%eps_filter)

            CALL cp_dbcsr_syevd(tempVV1, tempVV2, &
                                orbital_set%orbital_list(ispin, 2)%orb_energy, &
                                para_env=tidft_env%para_env, &
                                blacs_env=tidft_env%blacs_env)

            CALL dbcsr_multiply("N", "N", 1.0_dp, &
                                tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), &
                                tempVV2, &
                                0.0_dp, tempNV1, &
                                filter_eps=tidft_env%eps_filter)

            CALL dbcsr_copy(tidft_env%e_states(estate)%slater_det(idet)%mat_v(ispin), tempNV1)

            CALL dbcsr_release(tempNV1)
            CALL dbcsr_release(tempVV1)
            CALL dbcsr_release(tempVV2)

         ENDIF ! non-zero virtuals

      ENDDO

      DEALLOCATE(fake_e_per_block)

   END SUBROUTINE construct_canonical_orbitals

   SUBROUTINE get_orbital_indices(matrix, block_within_matrix, index_within_block)

      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      INTEGER, DIMENSION(:), INTENT(INOUT)               :: block_within_matrix, &
                                                            index_within_block

      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes
      INTEGER                                            :: mo_start, iblock, imo

      CALL dbcsr_get_info(matrix, col_blk_size=blk_sizes)
      mo_start = 1
      DO iblock = 1, SIZE(blk_sizes)
         block_within_matrix(mo_start:mo_start+blk_sizes(iblock)-1) = iblock
         index_within_block(mo_start:mo_start+blk_sizes(iblock)-1) = (/ (imo, imo = 1, blk_sizes(iblock)) /)
         mo_start = mo_start + blk_sizes(iblock)
      ENDDO

   END SUBROUTINE get_orbital_indices

   SUBROUTINE copy_excitation(dest, source)

      TYPE(tidft_excitation_type), INTENT(IN)             :: source
      TYPE(tidft_excitation_type), INTENT(INOUT)            :: dest

      CPASSERT(SIZE(source%delta_epsilon).EQ.SIZE(dest%delta_epsilon))

      ! copy data
      dest%order = source%order
      dest%excitation_energy = source%excitation_energy
      dest%unpairings_count = source%unpairings_count
      dest%delta_epsilon(:) = source%delta_epsilon(:)
      dest%occupied(:) = source%occupied(:)
      dest%virtual(:) = source%virtual(:)
      dest%spin(:) = source%spin(:)

   END SUBROUTINE copy_excitation

! **************************************************************************************************
!> \brief Convert MOs to independent variables
!> \param
!> \par History
!>       2024.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE MO2Var(mat_t, mat_var, ncols_var, mat_s, mat_oo_template, &
      eps_filter)

      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)      :: mat_t
      TYPE(dbcsr_type), DIMENSION(3), INTENT(INOUT)   :: mat_var
      INTEGER, DIMENSION(3), INTENT(IN)               :: ncols_var
      TYPE(dbcsr_type), INTENT(IN)                    :: mat_s
      TYPE(dbcsr_type), DIMENSION(:), INTENT(IN)      :: mat_oo_template
      REAL(KIND=dp), INTENT(IN)                       :: eps_filter

      CHARACTER(len=*), PARAMETER :: routineN = 'MO2Var', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: unit_nr, &
                                                            nspins, &
                                                            NPaired, &
                                                            ispin, &
                                                            iorb, &
                                                            copied_counter, &
                                                            empty_counter, &
                                                            orbital_index, &
                                                            isearch
      TYPE(dbcsr_type)                                   :: temp_t, temp_varB
      INTEGER, DIMENSION(2)                              :: NUnpaired
      INTEGER, DIMENSION(:,:), ALLOCATABLE               :: unpaired
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: orbitals_to_copy, &
                                                            where_to_copy, &
                                                            orbitals_to_copy_short, &
                                                            where_to_copy_short
      LOGICAL                                            :: paired_orb
      TYPE(cp_logger_type), POINTER                      :: logger

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      nspins = SIZE(mat_t)
      NPaired = ncols_var(1)
      DO ispin = 1, 2
         NUnpaired(ispin) = ncols_var(ispin+1)
      ENDDO

      IF ( nspins .EQ. 1 ) THEN
      ! This is fully restricted case
      ! Simply copy MOs to VAR(1)

         CALL dbcsr_copy(mat_var(1), mat_t(1))

      ELSE ! two spins

         IF ( NPaired .EQ. 0 ) THEN
         ! No electron pairs are required (could be UKS)
         ! Simply copy MOs to VAR(2) and VAR(3)

            IF (NUnpaired(1) .GT. 0) &
               CALL dbcsr_copy(mat_var(2), mat_t(1))
            IF (NUnpaired(2) .GT. 0) &
               CALL dbcsr_copy(mat_var(3), mat_t(2))

         ELSE
         ! This is restricted open-shell case

            ! First determine which orbitals are paired
            ! and which are unpaired
            ALLOCATE(unpaired( MAX(NUnpaired(1),NUnpaired(2)), 2 ))
            unpaired(:,:) = 0
            CALL PairedUnpaired( mat_t = mat_t, &
                    mat_var = mat_var, &
                    mat_s = mat_s, &
                    unpaired = unpaired, &
                    mat_oo_template = mat_oo_template, &
                    eps_filter = eps_filter)

            ! Collect mat_t collumns to appropriate mat_var
            DO ispin = 1, 2

               IF (NUnpaired(ispin).GT.0) THEN

                  ALLOCATE(orbitals_to_copy(NUnpaired(ispin)))
                  ALLOCATE(where_to_copy(NUnpaired(ispin)))

                  ! Copy unpaired columns directly to mat_var
                  DO iorb = 1, NUnpaired(ispin)

                     orbitals_to_copy(iorb) = unpaired(iorb, ispin)
                     where_to_copy(iorb) = iorb

                  ENDDO ! iorb

                  CALL copy_columns(matrix_dest = mat_var(ispin+1), &
                          matrix_source = mat_t(ispin), &
                          cols_dest = where_to_copy, &
                          cols_source = orbitals_to_copy)

                  DEALLOCATE(orbitals_to_copy)
                  DEALLOCATE(where_to_copy)

               ENDIF

            ENDDO ! ispin

            IF (NPaired.GT.0) THEN

               ! Use paired alpha orbtials and hope beta are the same
               ispin = 1
               IF (NUnpaired(ispin) .EQ. 0) ispin = 2

               ALLOCATE(orbitals_to_copy(NUnpaired(ispin)))
               ALLOCATE(where_to_copy(NUnpaired(ispin)))

               ! Do not copy paired columns directly to mat_var - too
               ! many columns, will be slow.
               ! Instead, copy all paired orbitals to the leftmost
               ! columns of mat_t and then split mat_t.
               ! This is much faster than copying all paired orbitals
               ! into the corresponding mat_var(1)

               ! find and record all paired orbitals in the rightmost
               ! columns of mat_t
               copied_counter = 0
               empty_counter = 0
               DO iorb = 1, NUnpaired(ispin)

                  orbital_index = iorb + NPaired

                  paired_orb = .TRUE.
                  DO isearch = 1, NUnpaired(ispin)
                     IF (orbital_index .EQ. unpaired(isearch, ispin)) &
                        paired_orb = .FALSE.
                  ENDDO

                  IF (paired_orb) THEN

                     copied_counter = copied_counter + 1
                     orbitals_to_copy(copied_counter) = orbital_index

                  ENDIF

                  ! determine orbitals to overwrite in the leftmost columns
                  IF ( unpaired(iorb,ispin) .LE. NPaired ) THEN

                     empty_counter = empty_counter + 1
                     where_to_copy(empty_counter) = unpaired(iorb,ispin)

                  ENDIF

               ENDDO ! iorb

               CPASSERT(empty_counter.EQ.copied_counter)

               ALLOCATE(orbitals_to_copy_short(copied_counter))
               ALLOCATE(where_to_copy_short(copied_counter))

               DO iorb = 1, copied_counter
                  orbitals_to_copy_short(iorb) = orbitals_to_copy(iorb)
                  where_to_copy_short(iorb) = where_to_copy(iorb)
               ENDDO

               DEALLOCATE(orbitals_to_copy)
               DEALLOCATE(where_to_copy)

               CALL dbcsr_create(temp_t, template = mat_t(ispin))
               CALL dbcsr_copy(temp_t, mat_t(ispin))
               CALL copy_columns(matrix_dest = temp_t, &
                       matrix_source = mat_t(ispin), &
                       cols_dest = where_to_copy_short, &
                       cols_source = orbitals_to_copy_short)
               CALL dbcsr_create(temp_varB, template = mat_var(ispin+1))
               CALL split_columns( &
                       temp_t, &
                       mat_var(1), &
                       temp_varB )

               CALL dbcsr_release(temp_t)
               CALL dbcsr_release(temp_varB)

               DEALLOCATE(orbitals_to_copy_short)
               DEALLOCATE(where_to_copy_short)

            ENDIF ! NPaired == 0

            DEALLOCATE(unpaired)

         ENDIF ! UKS or ROKS

      ENDIF ! RKS or (UKS, ROKS)

   END SUBROUTINE MO2Var

! **************************************************************************************************
!> \brief Given two sets of orbitals (mat_t) determine which orbitals are paired.
!>        The required number of paired and unpaired orbitals
!>        is given by the number of columns of mat_var
!>        ASSUMPTION: orthonormalized orbitals in mat_t
!> \param
!> \par History
!>       2024.02 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE PairedUnpaired(mat_t, mat_var, mat_s, unpaired, &
      mat_oo_template, eps_filter)

      TYPE(dbcsr_type), DIMENSION(2), INTENT(IN)      :: mat_t
      TYPE(dbcsr_type), DIMENSION(3), INTENT(IN)      :: mat_var
      TYPE(dbcsr_type), INTENT(IN)                    :: mat_s
      INTEGER, DIMENSION(:,:), INTENT(INOUT)          :: unpaired
      TYPE(dbcsr_type), DIMENSION(2), INTENT(IN)      :: mat_oo_template
      REAL(KIND=dp), INTENT(IN)                       :: eps_filter

      CHARACTER(len=*), PARAMETER :: routineN = 'PairedUnpaired', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: unit_nr, ispin, &
                                                            NOrbs, nspins, &
                                                            jspin, &
                                                            NUnpaired, &
                                                            ival, orb_index
      TYPE(dbcsr_type)                                   :: tempNN_symm,&
                                                            tempNO_ispin, &
                                                            tempNO_ispin_2, &
                                                            tempOO_ispin
      INTEGER, DIMENSION(:), ALLOCATABLE                 :: true_position
      REAL(KIND=dp)                                      :: warning_threshold, &
                                                            norm_of_next
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: norms
      TYPE(cp_logger_type), POINTER                      :: logger

      ! get a useful output_unit
      logger => cp_get_default_logger()
      IF (logger%para_env%is_source()) THEN
         unit_nr = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         unit_nr = -1
      ENDIF

      nspins = 2 ! This procedure should not be called for closed-shell
      CPASSERT(SIZE(mat_t).EQ.2)

      DO ispin = 1, nspins

         ! Spin of the projector
         jspin = 2
         IF (ispin .EQ. 2) jspin = 1

         CALL dbcsr_get_info(mat_var(ispin+1), nfullcols_total=NUnpaired)
         CALL dbcsr_get_info(mat_t(ispin), nfullcols_total=NOrbs)

         ! Compute norms of the projected orbitals
         ! norm_proj = diag [ Ti^t S Tj sig_inv_j Tj^t S Ti ]
         ! norm = [ Ti^t S Ti ]
         CALL dbcsr_create(tempNN_symm, &
                           template=mat_s, &
                           matrix_type=dbcsr_type_symmetric)
         CALL dbcsr_create(tempNO_ispin, &
                           template=mat_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(tempNO_ispin_2, &
                           template=mat_t(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)
         CALL dbcsr_create(tempOO_ispin, &
                           template=mat_oo_template(ispin), &
                           matrix_type=dbcsr_type_no_symmetry)

         ! R = Tj.(sigma_inv_j).Tj^t
         ! For orthogonal orbitals R = Tj.Tj^t
         ! It is safe to assume orthonormalized orbitals
         CALL dbcsr_multiply("N", "T", 1.0_dp, &
                             mat_t(jspin), &
                             mat_t(jspin), &
                             0.0_dp, tempNN_symm, &
                             filter_eps=eps_filter)

         ! norms = diag(Ti^t.S.R.S.Ti)
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             mat_s, &
                             mat_t(ispin), &
                             0.0_dp, tempNO_ispin, &
                             filter_eps=eps_filter)
         CALL dbcsr_multiply("N", "N", 1.0_dp, &
                             tempNN_symm, &
                             tempNO_ispin, &
                             0.0_dp, tempNO_ispin_2, &
                             filter_eps=eps_filter)
         CALL dbcsr_multiply("T", "N", 1.0_dp, &
                             tempNO_ispin, &
                             tempNO_ispin_2, &
                             0.0_dp, tempOO_ispin, &
                             filter_eps=eps_filter)

         ALLOCATE(norms(NOrbs))
         CALL dbcsr_get_diag(tempOO_ispin, norms)

         CALL dbcsr_release(tempNN_symm)
         CALL dbcsr_release(tempNO_ispin)
         CALL dbcsr_release(tempNO_ispin_2)
         CALL dbcsr_release(tempOO_ispin)

         ALLOCATE(true_position(NOrbs))
         CALL sort(norms, NOrbs, true_position)

         DO ival = 1, NOrbs

            orb_index = true_position(ival)

            IF ( ival .LE. NUnpaired ) THEN

               unpaired(ival, ispin) = orb_index

            !ELSE

            !   paired(ival - NUnpaired, ispin) = orb_index

            ENDIF

         ENDDO

         ! Perform simple checks before selecting orbitals
         ! Allowed deviation from the perfect 1.0 norm for paired orbitals
         warning_threshold = 0.05_dp
         norm_of_next = norms(true_position(NUnpaired+1))
         IF ( norm_of_next .LT. (1.0_dp - warning_threshold) ) THEN

            IF (unit_nr > 0) THEN
!               WRITE(unit_nr, '(T2,F10.5,T50, 2I6)') norm_of_next, &
               WRITE(*, *) norm_of_next, &
                  " norm for orbital ", NUnpaired+1, &
                  "(spin ", ispin , &
                  ") is too low to be considered spin paired."
               CPWARN("Unpaired orbitals will be forced to pair up!")
            ENDIF

         ENDIF

         ! It is OK for formally upaired orbitals to be paired

         DEALLOCATE(norms)
         DEALLOCATE(true_position)

         ENDDO ! ispin

   END SUBROUTINE PairedUnpaired

! **************************************************************************************************
!> \brief Store canonical-guess occupation IDs. GS occupied columns are 1..nocc;
!>        a promotion replaces occupied(k) by nocc+virtual(k). Cosine rank r
!>        is the number of IDs in one state that are not in the other.
! **************************************************************************************************
   SUBROUTINE store_guess_occupations(tidft_env, excitation_list, index_true, generated_nstates)

      TYPE(tidft_env_type), INTENT(INOUT)                :: tidft_env
      TYPE(tidft_excitation_list_type), INTENT(IN)       :: excitation_list
      INTEGER, DIMENSION(:), INTENT(IN)                  :: index_true
      INTEGER, INTENT(IN)                                :: generated_nstates

      INTEGER                                            :: estate, i, idet, iorder, ispin, istate, &
                                                            max_occ, n, n_full, nocc, nspins, &
                                                            occ_col, state_index, vir_col
      LOGICAL                                            :: found

      n_full = SIZE(tidft_env%e_states)
      nspins = tidft_env%nspins
      idet = 1
      IF (tidft_env%triplet_correction) idet = 2
      max_occ = 0
      DO estate = 1, n_full
         DO ispin = 1, nspins
            max_occ = MAX(max_occ, &
               tidft_env%e_states(estate)%slater_det(idet)%nelectrons_ab(ispin))
         ENDDO
      ENDDO
      IF (max_occ .LT. 1) max_occ = 1

      IF (ALLOCATED(tidft_env%guess_nocc)) DEALLOCATE(tidft_env%guess_nocc)
      IF (ALLOCATED(tidft_env%guess_occ)) DEALLOCATE(tidft_env%guess_occ)
      ALLOCATE(tidft_env%guess_nocc(n_full, nspins))
      ALLOCATE(tidft_env%guess_occ(n_full, nspins, max_occ))
      tidft_env%guess_nocc(:,:) = 0
      tidft_env%guess_occ(:,:,:) = -1

      DO ispin = 1, nspins
         n = tidft_env%e_states(1)%slater_det(idet)%nelectrons_ab(ispin)
         tidft_env%guess_nocc(1, ispin) = n
         DO i = 1, n
            tidft_env%guess_occ(1, ispin, i) = i
         ENDDO
      ENDDO

      DO istate = 1, generated_nstates
         estate = istate + 1
         IF (estate .GT. n_full) CYCLE
         state_index = index_true(istate)
         DO ispin = 1, nspins
            n = tidft_env%guess_nocc(1, ispin)
            tidft_env%guess_nocc(estate, ispin) = n
            tidft_env%guess_occ(estate, ispin, 1:n) = tidft_env%guess_occ(1, ispin, 1:n)
         ENDDO
         DO iorder = 1, excitation_list%excitations(state_index)%order
            ispin = excitation_list%excitations(state_index)%spin(iorder)
            occ_col = excitation_list%excitations(state_index)%occupied(iorder)
            vir_col = excitation_list%excitations(state_index)%virtual(iorder)
            nocc = tidft_env%guess_nocc(1, ispin)
            found = .FALSE.
            DO i = 1, nocc
               IF (tidft_env%guess_occ(estate, ispin, i) .EQ. occ_col) THEN
                  tidft_env%guess_occ(estate, ispin, i) = nocc + vir_col
                  found = .TRUE.
                  EXIT
               ENDIF
            ENDDO
            IF (.NOT. found) CPABORT("Guess occupation replace missed occupied column")
         ENDDO
      ENDDO

   END SUBROUTINE store_guess_occupations

END MODULE tidft_guess_methods

