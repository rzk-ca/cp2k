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
        tidft_prec_qhq_qsq, tidft_mult_constraint_det, tidft_mult_constraint_cosine, &
        tidft_mult_constraint_selective
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
                          description="Constraint function carried by the multipliers of "// &
                          "MULTIPLIER_COMPONENT. DETERMINANT: the signed normalized determinant "// &
                          "overlap g_IJ of the current state with each previous state. "// &
                          "SMALLEST_COSINE: the product over the spin channels of the signed smallest "// &
                          "principal cosine between the occupied subspaces of the two states (or of "// &
                          "its smooth surrogate, see PENALTY%COSINE_SMOOTH), augmented by the "// &
                          "quadratic term PENALTY%COSINE_PENALTY_STRENGTH*c_IJ**2. Both functions "// &
                          "vanish exactly where the states are orthogonal and give the same energies; "// &
                          "the smallest cosine is not suppressed by the other principal cosines of "// &
                          "the pair, so the dual iteration converges also for pairs of states that "// &
                          "differ by two orbitals in one spin channel. SELECTIVE_COSINE: every pair "// &
                          "starts with the DETERMINANT multiplier, so the first inner optimization "// &
                          "of a state is that of DETERMINANT; after each inner optimization a pair "// &
                          "whose overlap is above PENALTY%MAX_NORM_PHI_OFF_DIAG and whose product "// &
                          "P_IJ of the other principal cosines satisfies 2*lambda*P_IJ**2 < "// &
                          "PENALTY%MULTIPLIER_CURVATURE_GUESS (a pair for which the determinant "// &
                          "multiplier converges slowly) is switched to the SMALLEST_COSINE term with "// &
                          "its multiplier converted as mu_J -> mu_J*P_IJ. A switched pair whose "// &
                          "cosine does not decrease by PENALTY%COSINE_FALLBACK_RATIO is returned to "// &
                          "the determinant multiplier for the rest of the state. The interstate "// &
                          "penalty strength is raised only for the sake of determinant pairs and the "// &
                          "strength of the quadratic term does not follow it. Not available with "// &
                          "PENALTY%MULTIPLIER_ACCELERATION.", &
                          usage="MULTIPLIER_CONSTRAINT SMALLEST_COSINE", &
                          default_i_val=tidft_mult_constraint_det, &
                          enum_c_vals=s2a("DETERMINANT", "SMALLEST_COSINE", "SELECTIVE_COSINE"), &
                          enum_desc=s2a("Signed normalized determinant overlap.", &
                                        "Signed smallest principal cosine with quadratic augmentation.", &
                                        "Determinant overlap; smallest cosine only for stagnant pairs."), &
                          enum_i_vals=(/tidft_mult_constraint_det, tidft_mult_constraint_cosine, &
                                        tidft_mult_constraint_selective/))
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

      CALL keyword_create(keyword, __LOCATION__, name="MULTIPLIER_ACCELERATION", &
                          description="Used with MULTIPLIER_COMPONENT. Quasi-Newton dual update for the "// &
                          "pairs of states whose overlap is above the orthogonality threshold. The first "// &
                          "multiplier step is over-relaxed by 1 + MULTIPLIER_CURVATURE_GUESS/(2*lambda*P^2), "// &
                          "where P is the product of all principal cosines of the pair but the smallest; "// &
                          "later steps use an approximate inverse Jacobian of the overlaps with respect to "// &
                          "the multipliers that is corrected (Broyden rank-one update) with the observed "// &
                          "response to each multiplier change, which also captures the coupling between "// &
                          "pairs. After an accelerated step the penalty strength is raised only if the "// &
                          "overlaps responded with the wrong sign.", &
                          usage="MULTIPLIER_ACCELERATION .TRUE.", default_l_val=.FALSE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIPLIER_ACCEL_SINGLE_PAIR", &
                          description="Used with MULTIPLIER_ACCELERATION. Accelerate only one pair per "// &
                          "outer iteration: the previously accelerated pair while its overlap is above "// &
                          "the threshold, otherwise the pair with the largest overlap. The other pairs "// &
                          "take first-order steps, which are close to exact for them, and the inverse "// &
                          "Jacobian is corrected only with the response of the accelerated pair, because "// &
                          "the responses of the other pairs are dominated by their own steps and by the "// &
                          "coupling. FALSE accelerates all pairs above the threshold at once with a full "// &
                          "rank-one update.", &
                          usage="MULTIPLIER_ACCEL_SINGLE_PAIR .TRUE.", default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIPLIER_CURVATURE_GUESS", &
                          description="Used with MULTIPLIER_ACCELERATION. Guess (Hartree) of the curvature "// &
                          "of the energy along the smallest principal cosine of a pair of states, used by "// &
                          "the over-relaxed first multiplier step. Also used with MULTIPLIER_CONSTRAINT "// &
                          "SELECTIVE_COSINE: a pair above the overlap threshold is switched to the "// &
                          "smallest-cosine term if 2*lambda*P_IJ**2 is below this value.", &
                          usage="MULTIPLIER_CURVATURE_GUESS 1.0", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIPLIER_BOOST_MAX", &
                          description="Used with MULTIPLIER_ACCELERATION. Largest over-relaxation factor "// &
                          "of the first accelerated multiplier step.", &
                          usage="MULTIPLIER_BOOST_MAX 100.0", &
                          default_r_val=100.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIPLIER_BOOST_MIN", &
                          description="Used with MULTIPLIER_ACCELERATION. The first multiplier step of a "// &
                          "pair is over-relaxed only if the model factor 1 + MULTIPLIER_CURVATURE_GUESS/"// &
                          "(2*lambda*P^2) is at least this value, i.e. only for pairs with a small second "// &
                          "principal cosine; for the other pairs the plain first-order step is already "// &
                          "close to exact and the curvature guess does not apply.", &
                          usage="MULTIPLIER_BOOST_MIN 2.0", &
                          default_r_val=2.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIPLIER_BOOST_DAMPING", &
                          description="Used with MULTIPLIER_ACCELERATION. The over-relaxation of the first "// &
                          "multiplier step is damped to 1 + damping*(factor - 1): the curvature guess is "// &
                          "uncertain, an undershoot is completed by the following quasi-Newton step whereas "// &
                          "an overshoot perturbs the state.", &
                          usage="MULTIPLIER_BOOST_DAMPING 0.5", &
                          default_r_val=0.5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIPLIER_SECANT_NOISE", &
                          description="Used with MULTIPLIER_ACCELERATION. Smallest relative change of the "// &
                          "overlaps between two outer iterations (norm over the pairs) that is taken as a "// &
                          "genuine response to the multiplier change and used to correct the inverse "// &
                          "Jacobian; smaller changes are attributed to the noise of the inexact inner "// &
                          "optimization.", &
                          usage="MULTIPLIER_SECANT_NOISE 0.05", &
                          default_r_val=0.05_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIPLIER_SECANT_DEPTH", &
                          description="Used with MULTIPLIER_ACCELERATION. Number of previous outer "// &
                          "iterations whose multiplier changes and overlap responses are fitted at once "// &
                          "(multisecant correction of the model inverse Jacobian). One gives the rank-one "// &
                          "Broyden update; n resolves the coupling among n pairs after n updates.", &
                          usage="MULTIPLIER_SECANT_DEPTH 1", &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="MULTIPLIER_STEP_MAX", &
                          description="Used with MULTIPLIER_ACCELERATION. Largest change (Hartree) of a "// &
                          "multiplier in one accelerated step.", &
                          usage="MULTIPLIER_STEP_MAX 1.0", &
                          default_r_val=1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COSINE_PENALTY_STRENGTH", &
                          description="Used with MULTIPLIER_CONSTRAINT SMALLEST_COSINE and "// &
                          "SELECTIVE_COSINE. Strength rho "// &
                          "(Hartree) of the quadratic augmentation rho*c_IJ**2 of the smallest-cosine "// &
                          "constraint of a pair of states. It supplies a stiffness 2*rho along "// &
                          "the smallest principal cosine that does not depend on the other cosines "// &
                          "of the pair. A negative value (the default) means the initial "// &
                          "INTERSTATE_NON_PENALTY_STRENGTH. With SMALLEST_COSINE it is multiplied "// &
                          "together with the interstate penalty strength whenever the safeguard "// &
                          "raises the latter. With SELECTIVE_COSINE it applies to the switched pairs "// &
                          "only and does not follow the interstate penalty strength (see "// &
                          "COSINE_MAX_STRENGTH_STEPS).", &
                          usage="COSINE_PENALTY_STRENGTH 10.0", &
                          default_r_val=-1.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COSINE_FALLBACK_RATIO", &
                          description="Used with MULTIPLIER_CONSTRAINT SELECTIVE_COSINE. A pair that "// &
                          "carries the smallest-cosine term and whose overlap is still above "// &
                          "MAX_NORM_PHI_OFF_DIAG must decrease its constraint value c_IJ by this "// &
                          "factor between two consecutive inner optimizations. Otherwise the term "// &
                          "is taken to conflict with the constraint of another pair: the pair is "// &
                          "returned to the determinant multiplier it would have had and is not "// &
                          "switched again while the current state is optimized (see also "// &
                          "COSINE_MAX_STRENGTH_STEPS).", &
                          usage="COSINE_FALLBACK_RATIO 0.25", &
                          default_r_val=0.25_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COSINE_MAX_STRENGTH_STEPS", &
                          description="Used with MULTIPLIER_CONSTRAINT SELECTIVE_COSINE. Number of "// &
                          "times a pair that fails the COSINE_FALLBACK_RATIO test may instead have "// &
                          "COSINE_PENALTY_STRENGTH multiplied by PENALTY_STRENGTH_INCREASE_FACTOR "// &
                          "before it is returned to the determinant multiplier. The default of zero "// &
                          "returns the pair at the first failure and never changes the strength.", &
                          usage="COSINE_MAX_STRENGTH_STEPS 1", &
                          default_i_val=0)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, __LOCATION__, name="COSINE_SMOOTH", &
                          description="Used with MULTIPLIER_CONSTRAINT SMALLEST_COSINE and "// &
                          "SELECTIVE_COSINE. TRUE: the "// &
                          "smooth surrogate c = eps*[1 + SUM_k (1-lambda_k)**2/lambda_k]**(-1/2) of the "// &
                          "smallest cosine (lambda_k the squared principal cosines of the pair), which "// &
                          "has no kink where the two smallest cosines cross. FALSE: the plain signed "// &
                          "smallest cosine eps*cos(theta_1).", &
                          usage="COSINE_SMOOTH .TRUE.", default_l_val=.TRUE., &
                          lone_keyword_l_val=.TRUE.)
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
