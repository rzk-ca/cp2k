!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief input for the TIDFT section
!> \author Rustam Khaliullin
! **************************************************************************************************
MODULE input_cp2k_tidft
   USE bibliography,                    ONLY: Khaliullin2007,&
                                              Khaliullin2008,&
                                              Khaliullin2013,&
                                              Scheiber2018
   USE cp_output_handling,              ONLY: cp_print_key_section_create,&
                                              low_print_level,&
                                              debug_print_level,&
                                              add_last_numeric
   USE input_constants,                 ONLY: &
        use_mom_ref_coac, use_mom_ref_com, use_mom_ref_user, use_mom_ref_zero, &
        almo_deloc_none, almo_deloc_scf, almo_deloc_x, almo_deloc_x_then_scf, &
        almo_deloc_xalmo_1diag, almo_deloc_xalmo_scf, almo_deloc_xalmo_x, almo_frz_crystal, &
        almo_frz_none, almo_scf_diag, almo_scf_pcg, almo_scf_skip, almo_scf_trustr, atomic_guess, &
        cg_dai_yuan, cg_fletcher, cg_fletcher_reeves, cg_hager_zhang, cg_hestenes_stiefel, &
        cg_liu_storey, cg_polak_ribiere, cg_zero, restart_guess, op_loc_berry, op_loc_pipek, &
        optimizer_diis, optimizer_lin_eq_pcg, optimizer_pcg, optimizer_trustr, &
        spd_inversion_dense_cholesky, spd_inversion_ls_hotelling, spd_inversion_ls_taylor, &
        trustr_cauchy, trustr_dogleg, trustr_steihaug, tidft_prec_zero, &
        tidft_prec_s, xalmo_trial_r0_out, xalmo_trial_simplex, tidft_prec_h_s, tidft_prec_qhq_s, &
        tidft_prec_qhq_qsq
   USE input_val_types,                 ONLY: char_t,&
                                              integer_t,&
                                              lchar_t,&
                                              logical_t,&
                                              real_t
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_cp2k_excitation_analysis,  ONLY: create_nto_analysis_section,&
                                              create_mos_molden_section
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: s2a
   USE cp_dbcsr_api,                    ONLY: dbcsr_norm_frobenius, &
                                              dbcsr_norm_maxabsnorm
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_tidft'

   INTEGER, PARAMETER, PRIVATE :: optimizer_tidft_pcg = 1
   INTEGER, PARAMETER, PRIVATE :: optimizer_tidft_trustr = 2
   INTEGER, PARAMETER, PRIVATE :: optimizer_newton_pcg_solver = 3
   INTEGER, PARAMETER, PRIVATE :: optimizer_nlmo_pcg = 4

   PUBLIC :: create_tidft_section

CONTAINS

! **************************************************************************************************
!> \brief create the TIDFT section
!> \param section ...
!> \par History
!>       2020.04 created [Rustam Z Khaliullin]
!> \author Cindy Pham, Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE create_tidft_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_tidft_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="TIDFT", &
                          description="Time-independent density functional theory.", &
                          n_keywords=30, n_subsections=10, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_FILTER", &
                          description="Matrix sparsity filter.", &
                          usage="EPS_FILTER 1.e-6", default_r_val=1.e-9_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="N_EXCITED_STATES", &
                          description="Specifies the number of excited states (same multiplicity as the ground state).", &
                          usage="N_EXCITED_STATES 3", default_i_val=2)
      !CALL keyword_create(keyword, __LOCATION__, name="N_EXCITED_STATES", &
      !                    description="Specifies the number of excited states: singlet, doublet, triplet, ...", &
      !                    usage="N_EXCITED_STATES 3 1 4 0", &
      !                    n_var=-1, type_of_var=integer_t)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_UNPAIR", &
                          description="Specifies the number of allowed electron unpairings in ROKS", &
                          usage="MAX_UNPAIR 0", default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SPIN_PURE", &
                          description="Obtain spin pure electronic states.", &
                          usage="SPIN_PURE .TRUE.", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="SUGGESTED_MO_BLOCK_SIZE", &
                          description="User suggested block size for molecular orbital DBCSR matrices", &
                          usage="SUGGESTED_MO_BLOCK_SIZE 50", default_i_val=20)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MO_GUESS", &
                          description="Method that generates initial MOs.", &
                          usage="MO_GUESS RESTART", &
                          default_i_val=atomic_guess, &
                          enum_c_vals=s2a("RESTART", "ATOMIC"), &
                          enum_desc=s2a("Read from file.", &
                                        "Superposition of atomic densities."), &
                          enum_i_vals=(/restart_guess, atomic_guess/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MO_GUESS_RANDOM_SEED", &
                          description="Seed for generating random MOs for the excited states.", &
                          usage="MO_GUESS_RANDOM_SEED 10302", default_i_val=1000)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MO_GUESS_RANDOM_FRAC", &
                          description="Fraction of random MOs in the excited states.", &
                          usage="MO_GUESS_RANDOM_FRAC 0.4", default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="ENERGY_COMPONENT", &
                          description="Turn on the energy component of the loss function.", &
                          usage="ENERGY_COMPONENT .FALSE.", default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INTERSTATE_COMPONENT", &
                          description="Turn on the interstate nonorthogonality peanlty.", &
                          usage="INTERSTATE_COMPONENT .FALSE.", default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! aeph: toggle the second-order augmented Lagrangian interstate term (component 4)
      CALL keyword_create(keyword, __LOCATION__, name="AUGMENTED_COMPONENT", &
                          description="Turn on the second-order augmented Lagrangian term for the "// &
                          "interstate nonorthogonality penalty (loss function component 4). "// &
                          "Requires INTERSTATE_COMPONENT=.TRUE., else will abort.", &
                          usage="AUGMENTED_COMPONENT .TRUE.", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIPLIER_COMPONENT", &
                          description="Turn on the signed-overlap Lagrange multiplier term for the "// &
                          "interstate orthogonality constraints (loss function component 5). "// &
                          "One multiplier per previously optimized state is updated after each "// &
                          "inner optimization (method of multipliers), so that the states become "// &
                          "orthogonal at a finite INTERSTATE_NON_PENALTY_STRENGTH. "// &
                          "Requires INTERSTATE_COMPONENT=.TRUE. and STATE_BY_STATE=.TRUE.; "// &
                          "cannot be combined with AUGMENTED_COMPONENT.", &
                          usage="MULTIPLIER_COMPONENT .TRUE.", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIPLIER_CONSTRAINT", &
                          description="Used with MULTIPLIER_COMPONENT. Constraint function of each pair of "// &
                          "states to which its multiplier is attached. DETERMINANT: the normalized "// &
                          "determinant overlap of the two states. SPIN_CHANNEL: the normalized determinant "// &
                          "overlap of the spin channel with the smaller overlap. ORBITAL_PAIR: the cosine "// &
                          "between the pair of orbitals, one in each state, that is closest to orthogonal "// &
                          "(the smallest principal cosine between the occupied subspaces), with the pair "// &
                          "fixed during an inner optimization and re-selected at each multiplier update. "// &
                          "Convergence is always judged on the state overlap.", &
                          usage="MULTIPLIER_CONSTRAINT ORBITAL_PAIR", default_c_val="ORBITAL_PAIR")
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      CALL keyword_create(keyword, __LOCATION__, name="INTRASTATE_COMPONENT", &
                          description="Turn on the intrastate nonorthogonality peanlty.", &
                          usage="INTRASTATE_COMPONENT .FALSE.", default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="STATE_BY_STATE", &
                          description="Optimize electronic states one by one instead of optimizing all simultaneously.", &
                          usage="STATE_BY_STATE .TRUE.", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PRE_ORTHOGONALIZE_STATES", &
                          description="Pre-orthogonalize electronic states before the energy minimization.", &
                          usage="PRE_ORTHOGONALIZE_STATES .TRUE.", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NUMERICAL_GRADIENT", &
                          description="Debugging option. Use numerical, not analytical, gradients in WF optimization.", &
                          usage="NUMERICAL_GRADIENT .FALSE.", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NUMERICAL_GRADIENT_DX", &
                          description="Set dx when evaluating numerical gradients.", &
                          usage="NUMERICAL_GRADIENT_DX 1.0E-5", default_r_val=1.0E-06_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DIPOLE_REFERENCE_TYPE", &
                          description="Type of the reference point to calculate electric "// &
                          "dipole moments using the dipole integrals in the length form.", &
                          enum_c_vals=s2a("COM", "COAC", "USER_DEFINED", "ZERO"), &
                          enum_desc=s2a("Center of mass", &
                                        "Center of atomic charges", &
                                        "User-defined point", &
                                        "Origin of the coordinate system"), &
                          enum_i_vals=(/use_mom_ref_com, &
                                        use_mom_ref_coac, &
                                        use_mom_ref_user, &
                                        use_mom_ref_zero/), &
                          default_i_val=use_mom_ref_com)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="DIPOLE_REFERENCE_POINT", &
                          description="User-defined reference point in the dipole moment calculation.", &
                          usage="DIPOLE_REFERENCE_POINT x y z", &
                          repeats=.FALSE., n_var=3, type_of_var=real_t, unit_str='angstrom')
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      NULLIFY (subsection)
      CALL create_optimizer_section(subsection, optimizer_tidft_pcg)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_penalty_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_WF_guess_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_print_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_tidft_section

! **************************************************************************************************
!> \brief The optimizer section is a collection of keywords that are similar
!>        to all optimization methods (e.g. target error, number of iterations)
!> \param section ...
!> \param optimizer_id   allows to adapt the standard section for specific needs
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   RECURSIVE SUBROUTINE create_optimizer_section(section, optimizer_id)

      TYPE(section_type), POINTER                        :: section
      INTEGER, INTENT(IN)                                :: optimizer_id

      CHARACTER(len=*), PARAMETER :: routineN = 'create_optimizer_section', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: optimizer_type
      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      NULLIFY (section)

      ! choose the name of the section
      SELECT CASE (optimizer_id)
      CASE (optimizer_tidft_pcg)

         CALL section_create(section, __LOCATION__, name="TIDFT_OPTIMIZER_PCG", &
                             description="Controls the PCG optimization of ground- and excited-state MOs.", &
                             n_keywords=9, n_subsections=1, repeats=.FALSE.)
         optimizer_type = optimizer_pcg
         !optimizer_type = optimizer_diis
         !optimizer_type = optimizer_trustr

      CASE DEFAULT
         CPABORT("No default values allowed")
      END SELECT

      NULLIFY (keyword)

      ! add common keywords
      CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER", &
                          description="Maximum number of iterations", &
                          usage="MAX_ITER 100", default_i_val=20)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EPS_ERROR", &
                          description="Target value of the gradient norm", &
                          usage="EPS_ERROR 1.E-6", default_r_val=1.0E-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="WHICH_GRADIENT_NORM", &
                          description="Select which norm of the gradient to use when checking for convergence", &
                          usage="WHICH_GRADIENT_NORM RMS", &
                          default_i_val=dbcsr_norm_maxabsnorm, &
                          enum_c_vals=s2a("RMS", "MAX"), &
                          enum_desc=s2a("Use root mean square norm.", &
                                        "Use max norm."), &
                          enum_i_vals=(/dbcsr_norm_frobenius, dbcsr_norm_maxabsnorm/))
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)


      IF (optimizer_type .EQ. optimizer_pcg .OR. &
          optimizer_type .EQ. optimizer_lin_eq_pcg .OR. &
          optimizer_type .EQ. optimizer_trustr) THEN

         ! CINDY 5/1/23: add weight alpha for QHQ_S precx
         CALL keyword_create(keyword, __LOCATION__, name="PREC_WEIGHT_OVERLAP", &
                     description="Preconditioner weight of the overlap matrix. Should be between 0 and 1.", &
                     usage="PREC_WEIGHT_OVERLAP 0.6", default_r_val=0.5_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="MAX_ITER_OUTER_LOOP", &
                             description="Maximum number of iterations in the outer loop. "// &
                             "Use the outer loop to update the preconditioner and reset the conjugator. "// &
                             "This can speed up convergence significantly.", &
                             usage="MAX_ITER 10", default_i_val=0)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="PRECONDITIONER", &
                             description="Select a preconditioner for the conjugate gradient optimization", &
                             usage="PRECONDITIONER OVERLAP", &
                             default_i_val=tidft_prec_qhq_s, &
                             enum_c_vals=s2a("NONE", "H_S", "QHQ_S", "QHQ_QSQ", "S"), &
                             enum_desc=s2a("Do not use preconditioner.", &
                                           "The Kohn-Sham Hamiltonian plus the overlap. Cholesky inversion is used.", &
                                           "Virtual-virtual block of the Kohn-Sham Hamiltonian plus the overlap. "// &
                                           "Cholesky inversion is used.", &
                                           "Virtual-virtual block of the Kohn-Sham Hamiltonian plus the "// &
                                           "virtual-virtual block of the overlap. Pseudoinverse via diagonalization.", &
                                           "Cholesky inversion of the overlap matrix."), &
                             enum_i_vals=(/tidft_prec_zero, tidft_prec_h_s, tidft_prec_qhq_s, tidft_prec_qhq_qsq, &
                                           tidft_prec_s/))
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)
      ENDIF

      IF (optimizer_type .EQ. optimizer_pcg) THEN

         CALL keyword_create(keyword, __LOCATION__, name="LIN_SEARCH_EPS_ERROR", &
                             description="Target value of the gradient norm during the linear search", &
                             usage="LIN_SEARCH_EPS_ERROR 1.E-2", default_r_val=1.0E-3_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="LIN_SEARCH_STEP_SIZE_GUESS", &
                             description="The size of the first step in the linear search", &
                             usage="LIN_SEARCH_STEP_SIZE_GUESS 0.1", default_r_val=1.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ENDIF

      IF (optimizer_type .EQ. optimizer_pcg .OR. &
          optimizer_type .EQ. optimizer_trustr) THEN

         CALL keyword_create(keyword, __LOCATION__, name="CONJUGATOR", &
                             description="Various methods to compute step directions in the PCG optimization", &
                             usage="CONJUGATOR FLETCHER", &
                             default_i_val=cg_fletcher, &
                             enum_c_vals=s2a( &
                                           "ZERO", &
                                           !"POLAK_RIBIERE", &
                                           !"FLETCHER_REEVES", &
                                           !"HESTENES_STIEFEL", &
                                           "FLETCHER" &
                                           !"LIU_STOREY", &
                                           !"DAI_YUAN", &
                                           !"HAGER_ZHANG" &
                                           ), &
                             enum_desc=s2a( &
                                           "Steepest descent", &
                                           !"Polak and Ribiere", &
                                           !"Fletcher and Reeves", &
                                           !"Hestenes and Stiefel", &
                                           "Fletcher (Conjugate descent)" &
                                           !"Liu and Storey", &
                                           !"Dai and Yuan", &
                                           !"Hager and Zhang" &
                                           ), &
                             enum_i_vals=(/cg_zero, &
                                           !cg_polak_ribiere, &
                                           !cg_fletcher_reeves, &
                                           !cg_hestenes_stiefel, &
                                           cg_fletcher/))
                                           !cg_liu_storey, &
                                           !cg_dai_yuan, &
                                           !cg_hager_zhang/ &
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ENDIF

      IF (optimizer_type .EQ. optimizer_trustr) THEN

         CALL keyword_create(keyword, __LOCATION__, name="ALGORITHM", &
                             description="Selects an algorithm to solve the fixed-radius subproblem", &
                             usage="ALGORITHM CG", &
                             default_i_val=trustr_cauchy, &
                             enum_c_vals=s2a("CG", "CAUCHY", "DOGLEG"), &
                             enum_desc=s2a("Steihaug's iterative CG algorithm that does not invert model Hessian", &
                                           "Compute simple Cauchy point", &
                                           "Dogleg optimizer"), &
                             enum_i_vals=(/trustr_steihaug, trustr_cauchy, trustr_dogleg/))
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="ETA", &
                             description="Must be between 0.0 and 0.25. Rho value below which the "// &
                             "optimization of the model function is not accepted and the "// &
                             "optimization is restarted from the same point but decreased "// &
                             "trust radius. Rho is the ratio of the actual over predicted "// &
                             "change in the objective function", &
                             usage="ETA 0.1", default_r_val=0.25_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="MODEL_GRAD_NORM_RATIO", &
                             description="Stop the fixed-trust-radius (inner) loop optimization "// &
                             "once the ratio of the current norm of the model gradient over the "// &
                             "initial norm drops below this threshold", &
                             usage="MODEL_GRAD_NORM_RATIO 1.E-2", default_r_val=0.01_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="INITIAL_TRUST_RADIUS", &
                             description="Initial trust radius", &
                             usage="INITIAL_TRUST_RADIUS 0.1", default_r_val=0.1_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

         CALL keyword_create(keyword, __LOCATION__, name="MAX_TRUST_RADIUS", &
                             description="Maximum allowed trust radius", &
                             usage="MAX_TRUST_RADIUS 1.0", default_r_val=2.0_dp)
         CALL section_add_keyword(section, keyword)
         CALL keyword_release(keyword)

      ENDIF

   END SUBROUTINE create_optimizer_section

! **************************************************************************************************
!> \brief The section controls penalty methods
!> \param section ...
!> \par History
!>       2022.05 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE create_penalty_section(section)

      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      NULLIFY (section)

      CALL section_create(section, __LOCATION__, name="PENALTY", &
                          description="Add penalty terms to the energy functional.", &
                          n_keywords=3, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)
      ! CINDY 
      CALL keyword_create(keyword, __LOCATION__, name="LOOSE_EPS_SCF", &
                          description="Initial target value of the gradient norm.", &
                          usage="LOOSE_EPS_SCF", default_r_val=1.0E-2_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="NSTEPS_STRENGTH", &
                          variants=(/"MAX_CP_STEP"/), &
                          description="Number of steps in the gradual increase of the intrastate penalty strength "// &
                          "(and the gradual tightening of the convergence criterion).", &
                          usage="NSTEPS_STRENGTH 3", default_i_val=3)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INTERSTATE_NON_PENALTY_STRENGTH", &
                          description="Strength of the interstate nonorthogonality penalty term.", &
                          usage="INTERSTATE_NON_PENALTY_STRENGTH 0.5", &
                          default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! aeph: quadratic penalty coefficient (mu) of the augmented Lagrangian interstate term
      CALL keyword_create(keyword, __LOCATION__, name="AUGMENTED_PENALTY_STRENGTH", &
                          description="Quadratic penalty coefficient (mu) of the augmented Lagrangian "// &
                          "interstate term. Only used when AUGMENTED_COMPONENT is enabled; it is grown "// &
                          "across outer iterations by PENALTY_STRENGTH_INCREASE_FACTOR.", &
                          usage="AUGMENTED_PENALTY_STRENGTH 0.5", &
                          default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="INTRASTATE_NON_PENALTY_STRENGTH", &
                          description="Strength of the intrastate nonorthogonality penalty term.", &
                          usage="INTRASTATE_NON_PENALTY_STRENGTH 0.5", &
                          default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="PENALTY_STRENGTH_INCREASE_FACTOR", &
                          description="Factor that increases the penalty strength.", &
                          usage="PENALTY_STRENGTH_INCREASE_FACTOR 1.1", default_r_val=1.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_NORM_PHI_OFF_DIAG", &
                          description="Maximum allowed deviation from orthogonality for electronic states. "// &
                          "measured as the max norm of the off-diagonal elements of the interstate overlap matrix.", &
                          usage="MAX_NORM_PHI_OFF_DIAG 0.01", default_r_val=0.005_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
                  ! aeph: input threshold for normalized det(phi) "EPS_DET" in penalty section
      CALL keyword_create(keyword, __LOCATION__, name="EPS_DET", &
                          description="Defers ALM if the normalized interstate overlap determinant "// &
                          "det(phi) falls under the set threshold. Prevents near-singular interstates.", &
                          usage="EPS_DET 1.0E-10", &
                          default_r_val=1.0E-10_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      ! aeph: added tau as a input keyword 
      CALL keyword_create(keyword, __LOCATION__, name="AUGMENTED_RATIO", &
                          description="Grow mu only if the constraint |c| failed to decrease by this "// &
                          "factor relative to the previous outer iteration (ALM safeguard ratio).", &
                          usage="AUGMENTED_RATIO 0.25", &
                          default_r_val=0.25_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIPLIER_SAFEGUARD_RATIO", &
                          description="Used with MULTIPLIER_COMPONENT. The interstate penalty strength "// &
                          "is multiplied by PENALTY_STRENGTH_INCREASE_FACTOR only if the current-row "// &
                          "residual max_J |g_IJ| failed to decrease by this ratio relative to the "// &
                          "previous outer iteration.", &
                          usage="MULTIPLIER_SAFEGUARD_RATIO 0.25", &
                          default_r_val=0.25_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIPLIER_CURVATURE_SHIFT", &
                          description="Used with MULTIPLIER_COMPONENT. Estimated excitation energy "// &
                          "(Hartree) of the current state above the ground state. The dual step "// &
                          "uses 2*(lambda - shift) instead of 2*lambda, which is exact for the "// &
                          "two-level model of the collapse. Zero gives the plain multiplier update.", &
                          usage="MULTIPLIER_CURVATURE_SHIFT 0.0", &
                          default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIPLIER_INNER_TOL_FACTOR", &
                          description="Used with MULTIPLIER_COMPONENT. After each outer iteration the "// &
                          "convergence criterion of the next inner optimization is set to this "// &
                          "fraction of the current-row residual, but never looser than the previous "// &
                          "criterion and never tighter than EPS_ERROR. Zero means every inner "// &
                          "optimization is converged to EPS_ERROR.", &
                          usage="MULTIPLIER_INNER_TOL_FACTOR 0.1", &
                          default_r_val=0.1_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIPLIER_MAX_ITER", &
                          description="Used with MULTIPLIER_COMPONENT. Maximum number of outer "// &
                          "(dual update) iterations per state.", &
                          usage="MULTIPLIER_MAX_ITER 50", &
                          default_i_val=50)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_penalty_section

! **************************************************************************************************
!> \brief The section controls WF guess generation for excited states
!> \param section ...
!> \par History
!>       2023.12 created [Hanh D.M. Pham]
!> \author Hanh D.M. Pham
! **************************************************************************************************
   SUBROUTINE create_WF_guess_section(section)

      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      NULLIFY (section)

      CALL section_create(section, __LOCATION__, name="WF_GUESS", &
                          description="Control the generation of initial guess for all possible excited states.", &
                          n_keywords=3, n_subsections=0, repeats=.FALSE.)

      NULLIFY (keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_EXC_ENERGY", &
                          description="maximum allowed excitation energy for the guess (pre-SCF)", &
                          usage="MAX_EXC_ENERGY 100", default_r_val=1000.00_dp, unit_str='eV')
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MIN_EXC_ENERGY", &
                          description="If ground-state MOs are used to generate initial "// &
                          "excited-state MOs, this keyword specifies the minimum energy of "// &
                          "excitations from occupied ground-state orbitals to unoccupied orbitals. "// &
                          "Set to positive non-zero value to eliminate low-energy excitations.", &
                          usage="MIN_EXC_ENERGY 5.0", default_r_val=0.0_dp, unit_str='eV')
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MAX_EXC_ORDER", &
                          description="If ground-state MOs are used to generate initial "// &
                          "excited-state MOs, this keyword specifies the maximum number of "// &
                          "excitations from occupied ground-state orbitals to unoccupied orbitals.", &
                          usage="MAX_EXC_ORDER 2", default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="EXCITE_ALPHA_ONLY", &
                          description="Special case: Generate singlet guesses that excite alpha electrons only.", &
                          usage="EXCITE_ALPHA_ONLY .TRUE.", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword) 

   END SUBROUTINE create_WF_guess_section

! **************************************************************************************************
!> \brief The section controls TIDFT printing
!> \param section ...
!> \par History
!>       2024.07 created [Rustam Z Khaliullin]
!> \author Rustam Z Khaliullin
! **************************************************************************************************
   SUBROUTINE create_print_section(section)

      TYPE(section_type), POINTER                        :: section

      TYPE(section_type), POINTER                        :: subsection
      TYPE(section_type), POINTER                        :: print_key

      CPASSERT(.NOT. ASSOCIATED(section))
      NULLIFY (section, print_key)

      CALL section_create(section, __LOCATION__, name="PRINT", &
                          description="Printing of information during TIDFT runs.", repeats=.FALSE.)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "INTERSTATE_OVERLAP", &
                                       description="Requests printing of the interstate overlap matrix.", &
                                       print_level=debug_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "LOSS_FUNCTION_COMPONENTS", &
                                       description="Requests printing of the components of the loss function.", &
                                       print_level=debug_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      CALL cp_print_key_section_create(print_key, __LOCATION__, "DEBUG_INFO", &
                                       description="Requests printing of detailed information for debugging.", &
                                       print_level=debug_print_level, add_last=add_last_numeric, filename="__STD_OUT__")
      CALL section_add_subsection(section, print_key)
      CALL section_release(print_key)

      NULLIFY (subsection)
      CALL create_nto_analysis_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      NULLIFY (subsection)
      CALL create_mos_molden_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_print_section

END MODULE input_cp2k_tidft
