PROGRAM MAIN
    ! Amir Gasmi 11/16/2016
   ! USE FEAInput
    USE Assembly
    USE GaussQuadrature
    USE RHS_2D
    USE Loads_2D
    USE D_Matrix_2D
    USE K_Matrix_2D
    USE Mass_Matrix_2D
    USE BandedLinearSolverCHOLESKEY
    USE EigenValueSolver
    USE Stress_2D
    USE AccuracyAnalysis
    
    
    IMPLICIT NONE
    
    INCLUDE 'omp_lib.h'
    
    ! DECLARE ALL INPUT VARIABLES
    TYPE(Mesh)           :: Meshdata
    TYPE(Material)       :: Materialdata
    TYPE(BCDOF)          :: BCdata
    TYPE(Loads)          :: Loaddata
    TYPE(SimulationData) :: Simdata
    TYPE(Multiplicity)   :: Multiplicitydata ! This is for finding whether nodes in the loaddata belong to more than one element  
    CHARACTER(:), ALLOCATABLE :: OutputPath,LISTFILENAME,InputPath
    CHARACTER(LEN=100) :: ROWFMT
    ! DECLARE ALL LOOP VARIABLES
    INTEGER(IntKi) :: I,J,K,p,q
    
    


    
    !DECLARE ALL THE ASSEMBLY MATRICES
    INTEGER(IntKi), DIMENSION(:,:), ALLOCATABLE  :: ID,NegID
    REAL(ReKi), DIMENSION(:,:), ALLOCATABLE  :: ReducedGlobalStiffnessMatrix, ReducedGlobalMassMatrix,RR,InvReducedGlobalStiffnessMatrix,InvReducedGlobalMassMatrix 
    REAL(ReKi), DIMENSION(:), ALLOCATABLE    :: ReducedGlobalMassSpectral, ReducedRHS, BCVec
    !DECLARE ALL THE EIGENVALUE ANALYSIS MATRICES
    REAL(ReKi), DIMENSION(:,:), ALLOCATABLE  :: EigenMatrix,Zmatrix,GlobalNodalEigenVectorsX,GlobalNodalEigenVectorsY
    REAL(ReKi), DIMENSION(:), ALLOCATABLE    :: EigenValueR,EigenValueI
    INTEGER(IntKi)                           :: ierr
    INTEGER(IntKi), DIMENSION(:), ALLOCATABLE :: int,index,IFAIL, IWORK
    INTEGER(IntKi)                            :: IL, IU, LWORK, M
    REAL(ReKi)                                :: ABSTOL, VL, VU,DLAMCH
    REAL(ReKi), DIMENSION(:), ALLOCATABLE     :: WORK
    
           
    
    
    
    !DECLARE ALL VARIABLES FOR BANDED CHOLESKY SOLVER
    REAL(ReKi), DIMENSION(:,:), ALLOCATABLE    :: ReducedGlobalStiffnessBanded, U
    REAL(ReKi), DIMENSION(:), ALLOCATABLE    ::  R
    INTEGER(IntKi) :: INFO, KD,n
    !DECLARE ALL THE ELEMENT MATRICES
    REAL(ReKi), DIMENSION(:,:), ALLOCATABLE  :: KElement, MassElement, NodalCoordinates
    REAL(ReKi), DIMENSION(:), ALLOCATABLE    :: MassElementSpectral, SurfForce, SurfRHS, BodyForce, BodyRHS, RHS     
    !DECLARE ALL THE LAGRANGE INTERPOLATION AND ITS DERIVATIVES
    REAL(ReKi), DIMENSION(:,:), ALLOCATABLE  :: DL, LMatrix
    REAL(ReKi), DIMENSION(:), ALLOCATABLE    :: DL_Minus_one, DL_Plus_one
    REAL(ReKi), DIMENSION(:), ALLOCATABLE    :: Weights, Roots, QuadraturePoints
  
    !DECLARE MATERIAL PROPERTIES
    REAL(ReKi), DIMENSION(:,:), ALLOCATABLE  :: Dmatrix
    REAL(ReKi)  :: Rho
    ! DECLARE CPU TIME VARIABLES
    REAL(ReKi) :: T1,T2
    ! DECLARE STATIC OUTPUT DATA
    REAL(ReKi), DIMENSION(:,:), ALLOCATABLE    ::GlobalNodalDisplacement, NodalDisplacement, StressDistribution, StressDistributionI, &
                                                 StressAlongEdge
    INTEGER(IntKi)                            :: edge
    
    
    !........ Temporary Variable Declaration ............................
    LOGICAL, PARAMETER :: QUAD = .FALSE.
    CHARACTER(LEN=100) :: ARG
    
    !.......DECLARE ACCURACY ANALYSIS VARIABLES........................
    
    REAL(ReKi) :: PP, Xc, Yc , R1 , R2, Rerror,DispRErr, STrssRErr, TAssembly,TSolver 
    REAL(ReKi), DIMENSION(:,:), ALLOCATABLE :: TrueDisp,TStress
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !OutputPath = 'D:\\Year_2016\\fortran_codes\\FEA\\Common\\Data\\'
    !ListFileName = 'ListFileName.txt'
    !
    !IF(.NOT.(QUAD)) THEN
    !InputPath = 'D:\\Year_2016\\fortran_codes\\FEA\\PlateWithHole\\Spectral\\data\\'
    !ELSE
    !InputPath = 'D:\\Year_2016\\fortran_codes\\FEA\\PlateWithHole\\Quad\\data\\'
    !END IF
    
    
    
   IF (IARGC() /=3 ) THEN
       PRINT*,'THE CORRECT CALL IS DRIVER ListFileName Inputpath Outputpath'
       STOP
   ELSE
       CALL GETARG(1,ARG)
       listFileName = TRIM(ARG) 
       CALL GETARG(2,ARG)
       InputPath = TRIM(ARG)
       CALL GETARG(3,ARG)
       OutputPath = TRIM(ARG)
   END IF
   
       print*, 'ListFileName', ListFileName
       print*,'OutputPath',OutputPath
       print*,'InputPath',InputPath
        
    
    
    !---------------------------------------------------------------------------------------
    !-------------- Beginning Preparing Input Data ----------------------------------------
    !--------------------------------------------------------------------------------------
    
        
    !--------------------------------------- Populating Input file names -----------------------
    
    
    
    CALL PopulateInputFileNames(InputPath,Listfilename,Meshdata,Simdata,Materialdata,Loaddata,BCdata)
    
    
    
    ! ---------------------------------- Populating Mesh data ----------------------------
    
    
    CALL PopulateMeshData(InputPath, Meshdata)
    
    
    
    ! ---------------------------------- Populating Material data ----------------------------
 
    
    CALL PopulateMaterialData(InputPath,Materialdata)

    ! ---------------------------------- Populating Load data ----------------------------

    CALL PopulateLoadData(InputPath, Loaddata)
    !PRINT*,Loaddata%SurfLoads%Node
    !PRINT*,Loaddata%SurfLoads%DOF
    !PRINT*,Loaddata%SurfLoads%Element
    !PRINT*,Loaddata%SurfLoads%Edge
    !PRINT*,Loaddata%SurfLoads%Value
    
    
    ! ---------------------------------- Populating BC data ----------------------------
    
    CALL PopulateBCData(InputPath,BCdata)
    
    ! ---------------------------------- Populating Simulation data ----------------------------
    
    CALL PopulateSimulationData(InputPath, Simdata)
    
    
    
   !------------------------------------------------------------------------------------------------------------------------
   !----------------------    ASSEMBLY OF REDUCED MATRICES & RHS VECTORS --------------------------------------------------
   !-----------------------------------------------------------------------------------------------------------------------
   

    
   
   
    
    IF((Simdata%SimType%StaticAnalysis == .TRUE.) .OR. (Simdata%SimType%EigenValueAnalysis == .TRUE.) ) THEN  ! IS EIGENVALUE ANALYSIS OR STATIC ANALYSIS REQUESTED? IF NOT DO NOTHING.
        
        PRINT*,'................................... BEGIN ASSEMBLY ..............................................'
           ! ------------------------ Allocate all the assembly data ---------------------------------------------------
        ALLOCATE(ID(Meshdata%TotalNumberofDOFperNode,Meshdata%TotalNumberofNodes))
        ALLOCATE(NegID(Meshdata%TotalNumberofDOFperNode,Meshdata%TotalNumberofNodes))
        ALLOCATE(ReducedGlobalStiffnessMatrix(Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize, &
                                                Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize))
        ALLOCATE(ReducedRHS(Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize))
        ALLOCATE(BCVec(BCdata%BCSize))
        
        !----------------------------Allocate all the Element Matrices and Vectors ------------------------------
        ALLOCATE(KElement(Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodesPerElement, &
                            Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodesPerElement))
        ALLOCATE(SurfForce(4*Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodesPerEdge))
        ALLOCATE(SurfRHS(4*Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodesPerEdge))
        ALLOCATE(BodyForce(Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodesPerElement))
        ALLOCATE(BodyRHS(Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodesPerElement))
        ALLOCATE(RHS(Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodesPerElement))
        
        ALLOCATE(NodalCoordinates(Meshdata%TotalNumberOfNodesPerElement,2))
        !------------------------ Allocate Lagrange Arrays -------------------------------------------
        ALLOCATE(DL(Meshdata%TotalNumberOfNodesPerEdge,Meshdata%TotalNumberOfNodesPerEdge))
        ALLOCATE(DL_Minus_one(Meshdata%TotalNumberOfNodesPerEdge))
        ALLOCATE(DL_Plus_one(Meshdata%TotalNumberOfNodesPerEdge))
        ALLOCATE(Weights(Meshdata%TotalNumberOfNodesPerEdge))
        ALLOCATE(roots(Meshdata%TotalNumberOfNodesPerEdge))
        ALLOCATE(QuadraturePoints(Meshdata%TotalNumberOfNodesPerEdge))
        !--------------------- Allocate Material Matrix D-----------------------------
        ALLOCATE(Dmatrix(3,3))
        

        
        
        ! ------------------------- Construct BC Displacement vector and ID matrix -------------------------
        CALL ConstructID(BCdata%BCSize,Meshdata%TotalNumberofDOFperNode,Meshdata%TotalNumberofNodes,BCdata,ID,NegID)
        CALL Construct_BCDOF(BCdata%BCSize,Meshdata%TotalNumberofDOFperNode,Meshdata%TotalNumberofNodes,NegID,BCdata,BCVec)
        
        
        !................. Calculate the Material D matrix......................
        IF (Materialdata%Isotropic) THEN
            IF (Simdata%SimType%PlaneStress)    THEN
                CALL Isotropic_Plane_Stress(Materialdata%Iso%E,Materialdata%Iso%Nu,DMatrix)

            ELSE 
                CALL Isotropic_Plane_Strain(Materialdata%Iso%E,Materialdata%Iso%Nu,DMatrix)
            END IF
        ELSE
            DMatrix = Materialdata%Aniso%D
        END IF

        Rho = Materialdata%Rho
            !.............. Initializing Global matrices to 0.............
        ReducedGlobalStiffnessMatrix = 0
        ReducedRHS = 0
                
        !.............. Applying Loads...........................
        CALL FindMultiplicity(Meshdata,Loaddata,Multiplicitydata)
       
        
        !LOOP OVER ELEMENTS AND CONSTRUCT THE MASS MATRIX AND STIFFNESS MATRIX AND RHS VECTOR
        
        IF ( Simdata%SimType%EigenValueAnalysis == .TRUE. ) THEN ! MASS AT THE SAME TIME WITH STIFFNESS
            
            IF((Simdata%SimType%SpectralElements) .AND. (Simdata%SimType%Quadrature == 1) .AND. (.NOT. (Simdata%SimType%ReducedQuadrature))) THEN ! USE EFFICIENT MATRIX AND VECTOR BUILDING WHEN QUADRATURE POINTS = INTERPOLATION POINTS
                !------------------ Allocate Mass diagonal matrices in case when Spectral elements and quadrature are same
                ALLOCATE(MassElementSpectral(Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodesPerElement))
                ALLOCATE(ReducedGlobalMassSpectral(Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize))
                
                !.............. Initializing Global matrices to 0.............
                ReducedGlobalMassSpectral = 0
                
                        !............... Compute the Gauss-Legendre-Lobbatto Quadrature points and Weights .........
                CALL GaussLegendreLobatto(Meshdata%TotalNumberOfNodesPerEdge-1,QuadraturePoints,Weights)
                !....... Spectral Elements
                Roots = QuadraturePoints
                
                !.................... Compute all Lagrange Derivatives at the different integration points...........
                
                CALL DerivativeLagrangePArray(Meshdata%TotalNumberOfNodesPerEdge,Meshdata%TotalNumberOfNodesPerEdge,Roots,QuadraturePoints,DL)
                CALL DerivativeAllLagrangeP(Meshdata%TotalNumberOfNodesPerEdge,Roots,-1.0d0,DL_Minus_One)
                CALL DerivativeAllLagrangeP(Meshdata%TotalNumberOfNodesPerEdge,Roots,+1.0d0,DL_Plus_One)

                
                DO I =1,Meshdata%TotalNumberofElements ! BEGIN ASSEMBLY
                    
                    
                    !.............Copy Element Nodal Coordinates from Global coordinate using Connectivity matrix
                    
                    DO K =1,Meshdata%TotalNumberOfNodesperElement
                    NodalCoordinates(K,:) = Meshdata%GlobalCoordinates(Meshdata%ConnectivityMatrix(I,K),:)
                    END DO
                    
                    !...............Construct RHS of the Element......................
                    CALL ApplySurfLoad(I,Meshdata,Loaddata,Multiplicitydata,SurfForce)
                    CALL ApplyBodyLoad(I,Meshdata,Loaddata,Multiplicitydata,BodyForce)
                    
                    CALL Eff_Compute_Surface_Force(Meshdata%TotalNumberOfNodesPerEdge,Weights,DL,DL_Minus_one,DL_Plus_one,NodalCoordinates,SurfForce,SurfRHS)
                    
                    CALL Eff_Compute_Body_Force(Meshdata%TotalNumberOfNodesPerEdge,Weights,DL,NodalCoordinates,BodyForce,BodyRHS)
                    CALL Compute_RHS(Meshdata%TotalNumberOfNodesPerEdge,SurfRHS,BodyRHS,RHS) ! SUM THE CONTRIBUTION OF SURFACE AND BODY FORCES.
                     
                    !............. Construct the KElement matrix
                    CALL Eff_Construct_K_2DV2(Meshdata%TotalNumberOfNodesPerEdge,DL,Weights,NodalCoordinates,Dmatrix,KElement) 
                    !............. Construct the MassElement matrix......
                    CALL Eff_Construct_M_2DV2(Meshdata%TotalNumberOfNodesPerEdge,Weights,DL,NodalCoordinates,Rho,MassElementSpectral)
                    !............ Add the element matrices to the global reduced matrices
                    CALL Assemble_K_RHS(I,Meshdata%TotalNumberofElements,Meshdata%TotalNumberofNodesperElement,BCdata%BCSize, &  
                         Meshdata%TotalNumberofDOFperNode,Meshdata%TotalNumberofNodes,                   &
                         ID,NegID,Meshdata%ConnectivityMatrix,BCVec,KElement,RHS, &
                         ReducedGlobalStiffnessMatrix,ReducedRHS,MassElementSpectral,ReducedGlobalMassSpectral)
                    
                END DO ! END ASSEMBLY
                
            
            ELSE IF((Simdata%SimType%Quadrature == 0) .AND. (.NOT. (Simdata%SimType%ReducedQuadrature))) THEN ! USE THE GENERAL APPROACH TO BUILD ELEMENT MATRICES AND VECTORS (here is uniform Elements)
                !............ Allocation of the necessary matrices...................
                ALLOCATE(LMatrix(Meshdata%TotalNumberOfNodesPerEdge,Meshdata%TotalNumberOfNodesPerEdge))
                ALLOCATE(MassElement(Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodesPerElement, &
                            Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodesPerElement))
                ALLOCATE(ReducedGlobalMassMatrix(Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize, &
                                                    Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize))
                
                                !.............. Initializing Global matrices to 0.............
                ReducedGlobalMassMatrix = 0
                
                        !............... Compute the Gauss-Legendre Quadrature points and Weights .........
                CALL GaussLegendre(Meshdata%TotalNumberOfNodesPerEdge,-1.0d0,1.0d0,QuadraturePoints,Weights)
                !....... Non-Spectral Elements Roots are uniformly spaced on the edge
                
                CALL UniformDistributionGenerator(Meshdata%TotalNumberOfNodesPerEdge,-1.0d0,1.0d0,Roots)
                
                !.................... Compute all Lagrange Derivatives at the different integration points...........
                
                CALL DerivativeLagrangePArray(Meshdata%TotalNumberOfNodesPerEdge,Meshdata%TotalNumberOfNodesPerEdge,Roots,QuadraturePoints,DL)
                CALL DerivativeAllLagrangeP(Meshdata%TotalNumberOfNodesPerEdge,Roots,-1.0d0,DL_Minus_One)
                CALL DerivativeAllLagrangeP(Meshdata%TotalNumberOfNodesPerEdge,Roots,+1.0d0,DL_Plus_One)
                
                !.................... Compute all Lagrange polynonials at the different integration points...........
                
                CALL LagrangePArray(Meshdata%TotalNumberOfNodesPerEdge,Meshdata%TotalNumberOfNodesPerEdge,Roots,QuadraturePoints,LMatrix)
                
                
                DO I =1,Meshdata%TotalNumberOfElements ! BEGIN ASSEMBLY
                    
                    
                    !.............Copy Element Nodal Coordinates from Global coordinate using Connectivity matrix
                    
                    DO K =1,Meshdata%TotalNumberOfNodesPerElement
                    NodalCoordinates(K,:) = Meshdata%GlobalCoordinates(Meshdata%ConnectivityMatrix(I,K),:)
                    END DO
                    
                    !...............Construct RHS of the Element......................
                    CALL ApplySurfLoad(I,Meshdata,Loaddata,Multiplicitydata,SurfForce) ! Applying the portion of the surface loading of element I
                    CALL ApplyBodyLoad(I,Meshdata,Loaddata,Multiplicitydata,BodyForce)  ! Applying the portion of the surface loading of element I
                    CALL Compute_Surface_Force(Meshdata%TotalNumberOfNodesPerEdge,Weights,LMatrix,DL,DL_Minus_one,DL_Plus_one,NodalCoordinates,SurfForce,SurfRHS)
                    CALL Compute_Body_Force(Meshdata%TotalNumberOfNodesPerEdge,Weights,LMatrix,DL,NodalCoordinates,BodyForce,BodyRHS)
                    CALL Compute_RHS(Meshdata%TotalNumberOfNodesPerEdge,SurfRHS,BodyRHS,RHS) ! SUM THE CONTRIBUTION OF SURFACE AND BODY FORCES.
                    
                    !............. Construct the KElement matrix
                    CALL  Construct_K_2D(Meshdata%TotalNumberOfNodesPerEdge,LMatrix,DL,Weights,NodalCoordinates,Dmatrix,KElement) 
                    !............. Construct the MassElement matrix......
                    CALL Construct_M_2D(Meshdata%TotalNumberOfNodesPerEdge,Weights,LMatrix,DL,NodalCoordinates,Rho,MassElement)
                     !............ Add the element matrices to the global reduced matrices
                    CALL Assemble_K_M_RHS(I,Meshdata%TotalNumberofElements,Meshdata%TotalNumberofNodesperElement,BCdata%BCSize, &  
                         Meshdata%TotalNumberofDOFperNode,Meshdata%TotalNumberofNodes,                   &
                         ID,NegID,Meshdata%ConnectivityMatrix,BCVec,KElement,RHS, &
                         ReducedGlobalStiffnessMatrix,ReducedRHS,MassElement,ReducedGlobalMassMatrix)
                END DO ! End Assembly
                
            END IF
            
            
            
            
        ELSE ! DO NOT ASSEMBLE MASS MATRIX
            
            T1 = OMP_GET_WTIME()
            
            IF((Simdata%SimType%SpectralElements) .AND. (Simdata%SimType%Quadrature == 1) .AND. (.NOT. (Simdata%SimType%ReducedQuadrature))) THEN ! USE EFFICIENT MATRIX AND VECTOR BUILDING WHEN QUADRATURE POINTS = INTERPOLATION POINTS
                !------------------ Allocate Mass diagonal matrices in case when Spectral elements and quadrature are same
                !ALLOCATE(MassElementSpectral(Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodesPerElement))
                !ALLOCATE(ReducedGlobalMassSpectral(Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize))
                
                
                        !............... Compute the Gauss-Legendre-Lobbatto Quadrature points and Weights .........
                CALL GaussLegendreLobatto(Meshdata%TotalNumberOfNodesPerEdge-1,QuadraturePoints,Weights)
                !....... Spectral Elements
                Roots = QuadraturePoints
                
                !.................... Compute all Lagrange Derivatives at the different integration points...........
                
                CALL DerivativeLagrangePArray(Meshdata%TotalNumberOfNodesPerEdge,Meshdata%TotalNumberOfNodesPerEdge,Roots,QuadraturePoints,DL)
                CALL DerivativeAllLagrangeP(Meshdata%TotalNumberOfNodesPerEdge,Roots,-1.0d0,DL_Minus_One)
                CALL DerivativeAllLagrangeP(Meshdata%TotalNumberOfNodesPerEdge,Roots,+1.0d0,DL_Plus_One)                

            

                
                DO I =1,Meshdata%TotalNumberofElements ! BEGIN ASSEMBLY
                    
                    
                    !.............Copy Element Nodal Coordinates from Global coordinate using Connectivity matrix
                    
                    DO K =1,Meshdata%TotalNumberOfNodesperElement
                    NodalCoordinates(K,:) = Meshdata%GlobalCoordinates(Meshdata%ConnectivityMatrix(I,K),:)
                    END DO
                    
                    !...............Construct RHS of the Element......................
                    CALL ApplySurfLoad(I,Meshdata,Loaddata,Multiplicitydata,SurfForce)
                    CALL ApplyBodyLoad(I,Meshdata,Loaddata,Multiplicitydata,BodyForce)
                    
                    CALL Eff_Compute_Surface_Force(Meshdata%TotalNumberOfNodesPerEdge,Weights,DL,DL_Minus_one,DL_Plus_one,NodalCoordinates,SurfForce,SurfRHS)
                    
                    CALL Eff_Compute_Body_Force(Meshdata%TotalNumberOfNodesPerEdge,Weights,DL,NodalCoordinates,BodyForce,BodyRHS)
                    CALL Compute_RHS(Meshdata%TotalNumberOfNodesPerEdge,SurfRHS,BodyRHS,RHS) ! SUM THE CONTRIBUTION OF SURFACE AND BODY FORCES.
                     
                    !............. Construct the KElement matrix
                    
                    
                    !CALL Eff_Construct_K_2D(Meshdata%TotalNumberOfNodesPerEdge,DL,Weights,NodalCoordinates,Dmatrix,KElement)
                    CALL Eff_Construct_K_2DV2(Meshdata%TotalNumberOfNodesPerEdge,DL,Weights,NodalCoordinates,Dmatrix,KElement)
                     
                    N= Meshdata%TotalNumberOfNodesPerEdge

                    
                    !............. Construct the MassElement matrix......
                    !CALL Eff_Construct_M_2D(Meshdata%TotalNumberOfNodesPerEdge,Weights,DL,NodalCoordinates,Rho,MassElementSpectral)
                    !............ Add the element matrices to the global reduced matrices
                    CALL Assemble_K_RHS(I,Meshdata%TotalNumberofElements,Meshdata%TotalNumberofNodesperElement,BCdata%BCSize, &  
                         Meshdata%TotalNumberofDOFperNode,Meshdata%TotalNumberofNodes,                   &
                         ID,NegID,Meshdata%ConnectivityMatrix,BCVec,KElement,RHS, &
                         ReducedGlobalStiffnessMatrix,ReducedRHS)

                    
                    
                END DO ! END ASSEMBLY
                
            
            ELSE IF((Simdata%SimType%Quadrature == 0) .AND. (.NOT. (Simdata%SimType%ReducedQuadrature))) THEN ! USE THE GENERAL APPROACH TO BUILD ELEMENT MATRICES AND VECTORS (here is uniform Elements)
                !............ Allocation of the necessary matrices...................
                ALLOCATE(LMatrix(Meshdata%TotalNumberOfNodesPerEdge,Meshdata%TotalNumberOfNodesPerEdge))

                
                                !.............. Initializing Global matrices to 0.............

                
                        !............... Compute the Gauss-Legendre Quadrature points and Weights .........
                CALL GaussLegendre(Meshdata%TotalNumberOfNodesPerEdge,-1.0d0,1.0d0,QuadraturePoints,Weights)
                !....... Non-Spectral Elements Roots are uniformly spaced on the edge
                
                CALL UniformDistributionGenerator(Meshdata%TotalNumberOfNodesPerEdge,-1.0d0,1.0d0,Roots)
                
                !.................... Compute all Lagrange Derivatives at the different integration points...........
                
                CALL DerivativeLagrangePArray(Meshdata%TotalNumberOfNodesPerEdge,Meshdata%TotalNumberOfNodesPerEdge,Roots,QuadraturePoints,DL)
                CALL DerivativeAllLagrangeP(Meshdata%TotalNumberOfNodesPerEdge,Roots,-1.0d0,DL_Minus_One)
                CALL DerivativeAllLagrangeP(Meshdata%TotalNumberOfNodesPerEdge,Roots,+1.0d0,DL_Plus_One)
                
                !.................... Compute all Lagrange polynonials at the different integration points...........
                
                CALL LagrangePArray(Meshdata%TotalNumberOfNodesPerEdge,Meshdata%TotalNumberOfNodesPerEdge,Roots,QuadraturePoints,LMatrix)
                
                
                DO I =1,Meshdata%TotalNumberOfElements ! BEGIN ASSEMBLY
                    
                    
                    !.............Copy Element Nodal Coordinates from Global coordinate using Connectivity matrix
                    
                    DO K =1,Meshdata%TotalNumberOfNodesperElement
                    NodalCoordinates(K,:) = Meshdata%GlobalCoordinates(Meshdata%ConnectivityMatrix(I,K),:)
                    END DO
                    
                    !...............Construct RHS of the Element......................
                    CALL ApplySurfLoad(I,Meshdata,Loaddata,Multiplicitydata,SurfForce) ! Applying the portion of the surface loading of element I
                    CALL ApplyBodyLoad(I,Meshdata,Loaddata,Multiplicitydata,BodyForce)  ! Applying the portion of the surface loading of element I
                    CALL Compute_Surface_Force(Meshdata%TotalNumberOfNodesPerEdge,Weights,LMatrix,DL,DL_Minus_one,DL_Plus_one,NodalCoordinates,SurfForce,SurfRHS)
                    CALL Compute_Body_Force(Meshdata%TotalNumberOfNodesPerEdge,Weights,LMatrix,DL,NodalCoordinates,BodyForce,BodyRHS)
                    CALL Compute_RHS(Meshdata%TotalNumberOfNodesPerEdge,SurfRHS,BodyRHS,RHS) ! SUM THE CONTRIBUTION OF SURFACE AND BODY FORCES.
                    
                    !............. Construct the KElement matrix
                    CALL  Construct_K_2D(Meshdata%TotalNumberOfNodesPerEdge,LMatrix,DL,Weights,NodalCoordinates,Dmatrix,KElement) 
                     !............ Add the element matrices to the global reduced matrices
                    CALL Assemble_K_M_RHS(I,Meshdata%TotalNumberofElements,Meshdata%TotalNumberofNodesperElement,BCdata%BCSize, &  
                         Meshdata%TotalNumberofDOFperNode,Meshdata%TotalNumberofNodes,                   &
                         ID,NegID,Meshdata%ConnectivityMatrix,BCVec,KElement,RHS, &
                         ReducedGlobalStiffnessMatrix,ReducedRHS)
                END DO ! End Assembly
                
            END IF  
            
            T2 = OMP_GET_WTIME()  
            
        END IF
        
        PRINT*,'................................... END ASSEMBLY ..............................................'
        TAssembly = T2-T1
        
        PRINT*,'ASSEMBLY TIME :', TAssembly
        
        DEALLOCATE(NodalCoordinates)
    END IF  ! IS EIGENVALUE ANALYSIS OR STATIC ANALYSIS REQUESTED? IF NOT DO NOTHING.
    
    
    
    
    
    
    !................................................................................................................................!
    !............................... LINEAR SOLVER: STATIC ANALYSIS; FIND DISPLACEMENT U.............................................!
    !................................................................................................................................!
    
   
    IF(Simdata%SimType%StaticAnalysis ) THEN
        PRINT*,'.........................BEGIN STATIC ANALYSIS SOLVE........................'
        CALL FindHalfBandWidth(Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize,ReducedGlobalStiffnessMatrix,KD)
        print*,'the Half bandwidth of ReducedGlobalStiffnessMatrix of size',Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize,'is:',KD
        ALLOCATE(ReducedGlobalStiffnessBanded(KD+1,Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize))
        ALLOCATE(U(Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize,1))
        CALL TransformAtoAB(Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize, KD, ReducedGlobalStiffnessMatrix, ReducedGlobalStiffnessBanded)
        
        ! COPY ReducedRHS to U
        U(:,1) = ReducedRHS
        
        
        T2 = OMP_GET_WTIME()  
        CALL DandedSymmPOSSsolver(Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize, KD, 1, ReducedGlobalStiffnessBanded, U, INFO )
        T2 = OMP_GET_WTIME()  
        TSolver = T2-T1
        print*,'INFO OF THE CHOLESKY SOLVER, NONE ZERO THERE IS AN ERROR OR THE SOLVER COULD NOT FIND A SOLUTION:',INFO
        PRINT*,'THE Solver Time:',TSolver
        n = Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
        ALLOCATE(R(n))
        CALL MatrixVectorMultiply(n,n,U(:,1),ReducedGlobalStiffnessMatrix,R)
        PRINT*,'THE RELATIVE ERROR IS:',NORM2(n,R-ReducedRHS)/NORM2(n,ReducedRHS)
        PRINT*,'.........................END STATIC ANALYSIS SOLVE........................'
        
        !PRINT*,'........................ TEST OF THE INVERSE USING CHOLESKY FACTORIZATION FOR A BANDED MATRIX AND INVERSE OF UPPER TRIANGULAR...........'
        !
        !N = Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
        !CALL InverseBandedSymmetricPositive(N, KD,ReducedGlobalStiffnessMatrix,InvReducedGlobalStiffnessMatrix, INFO )
        !ALLOCATE(RR(N,N))
        !CALL MatrixMultiply(N,N,N,ReducedGlobalStiffnessMatrix,InvReducedGlobalStiffnessMatrix,RR)
        !PRINT*,INFO
        !
        !PRINT*,'RR(1,1)=K*k^-1=',RR(1,1),'RR(10,10)=K*k^-1=',RR(10,10),'RR(100,1)=K*k^-1=',RR(100,1)
        !DEALLOCATE(R)
        !DEALLOCATE(InvReducedGlobalStiffnessMatrix)
        !PRINT*,'........................END TEST OF THE INVERSE USING CHOLESKY FACTORIZATION FOR A BANDED MATRIX AND INVERSE OF UPPER TRIANGULAR...........'
    END IF
  
    
    
    
    
    
    
     !................................................................................................................................!
     !............................... EIGENVALUE ANALYSIS; FIND EIGENVALUES AND EIGENVECTORS U.........................................!
     !................................................................................................................................!
    IF((Simdata%SimType%EigenValueAnalysis == .TRUE.)) THEN
        
        
        PRINT*,'........................... BEGIN EIGENVALUE ANALYSIS..........................'
        
        !........................Allocate all the Eigenvalue analysis matrices........................................................
        N =Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
        M = Simdata%OutputReq%EIGEN%NumberOfEigenValues
        ALLOCATE(EigenMatrix(N,N))
        IF((Simdata%OutputReq%EIGEN%ALL) .AND. (Simdata%OutputReq%EIGEN%EigenVectors)) THEN
        ALLOCATE(Zmatrix(N,N))
        ELSE IF (Simdata%OutputReq%EIGEN%EigenVectors) THEN
        ALLOCATE(Zmatrix(N,M))
        END IF
        ALLOCATE(EigenValueR(N))
        
        
        
        
        IF((Simdata%SimType%SpectralElements) .AND. (Simdata%SimType%Quadrature == 1) .AND. (.NOT. (Simdata%SimType%ReducedQuadrature))) THEN ! CONDITIONS FOR MAKING THE MASS MATRIX DIAGONAL MATRIX
            ALLOCATE(int(N))
            ALLOCATE(index(N))
            ALLOCATE(EigenValueI(N))
            CALL CPU_TIME(T1)
            DO I=1,Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
            EigenMatrix(I,:) = ReducedGlobalStiffnessMatrix(I,:)/ReducedGlobalMassSpectral(I)
            END DO
            
        
            !..........................Finding the EigenValues..............................
            IF(Simdata%OutputReq%EIGEN%ALL) THEN ! ALL EIGENVALUES 
                IF(Simdata%OutputReq%EIGEN%EigenVectors) THEN ! INCLUDE IN THE SEARCH THE EIGENVECTORS
                    N = Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
                                   
                    CALL elmhes(N,N,1,N,EigenMatrix,int) ! Transform the General Matrix into a Hessenberg
                    CALL eltran(N,N,1,N,EigenMatrix,int,Zmatrix) ! Find the Transformation matrix useful for finding the EigenVectors
                    CALL hqr2(N,N,1,N,EigenMatrix,EigenValueR,EigenValueI,Zmatrix,ierr) ! Finds the Eigenvalues and Vectors of a Hessenberg Matrix.
                    CALL AscendingOrder(N,EigenValueR,index)
                    print*,'ierr; the error',ierr
                ELSE
                    N = Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
                                   
                    CALL elmhes(N,N,1,N,EigenMatrix,int) ! Transform the General Matrix into a Hessenberg
                    !CALL eltran(N,N,1,N,EigenMatrix,int,Zmatrix) ! Find the Transformation matrix useful for finding the EigenVectors
                    !CALL hqr2(N,N,1,N,EigenMatrix,EigenValueR,EigenValueI,Zmatrix,ierr) ! Finds the Eigenvalues and Vectors of a Hessenberg Matrix.
                    CALL hqr(N,N,1,N,EigenMatrix,EigenValueR,EigenValueI,ierr)  ! Finds all the Eigenvalues of a Hessenberg Matrix.
                    CALL AscendingOrder(N,EigenValueR,index)
                    print*,'ierr; the error',ierr
                END IF
            
            ELSE ! ONLY A FEW EIGENVALUES ARE REQUESTED 
                
                N = Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize   
                ALLOCATE(ReducedGlobalMassMatrix(N,N)) 
                ReducedGlobalMassMatrix = 0
                DO I =1,N
                   ReducedGlobalMassMatrix(I,I) = ReducedGlobalMassSpectral(I)
                END DO
                
                        
                IF(Simdata%OutputReq%EIGEN%EigenVectors) THEN ! INCLUDE IN THE SEARCH THE EIGENVECTORS
                    N = Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize                   
                    
                    VL = 0
                    VU = 0
                    IL = Simdata%OutputReq%EIGEN%LowerIndex
                    IU = Simdata%OutputReq%EIGEN%LowerIndex+Simdata%OutputReq%EIGEN%NumberOfEigenValues-1
                    ABSTOL = 2 *DLAMCH('S')
                    LWORK = max(1,8*N)
                    ALLOCATE(WORK(MAX(1,8*N)))
                    ALLOCATE(IFAIL(N))
                    ALLOCATE(IWORK(5*N))
                    CALL CPU_TIME(T1)
                    CALL DSYGVX( 1, 'V', 'I', 'U', N, ReducedGlobalStiffnessMatrix, N, ReducedGlobalMassMatrix, N, &
                          VL, VU, IL, IU, ABSTOL, M, EigenValueR, Zmatrix, N, WORK, LWORK, IWORK, IFAIL, INFO )
                    CALL CPU_TIME(T2)
                    PRINT*,'Info on the EIGENVALUES OF THE SYSTEM OF MASS STIFFNESS IS:',INFO
                    PRINT*,'THE CPU TIME IT TOOK THE EIGENVALUE ANALYSIS IS:', T2-T1
                    
                    DEALLOCATE(IWORK)
                    DEALLOCATE(WORK)
                    DEALLOCATE(IFAIL)
                    
                    
                    
                ELSE ! All EigenValues without EigenVectors
                     N = Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
                    VL = 0
                    VU = 0
                    IL = Simdata%OutputReq%EIGEN%LowerIndex
                    IU = Simdata%OutputReq%EIGEN%LowerIndex+Simdata%OutputReq%EIGEN%NumberOfEigenValues-1
                    ABSTOL = 2 *DLAMCH('S')
                    LWORK = max(1,8*N)
                    ALLOCATE(WORK(MAX(1,8*N)))
                    ALLOCATE(IFAIL(N))
                    ALLOCATE(IWORK(5*N))
                    CALL CPU_TIME(T1)
                    CALL DSYGVX( 1, 'N', 'I', 'U', N, ReducedGlobalStiffnessMatrix, N, ReducedGlobalMassMatrix, N, &
                          VL, VU, IL, IU, ABSTOL, M, EigenValueR, Zmatrix, N, WORK, LWORK, IWORK, IFAIL, INFO )
                    CALL CPU_TIME(T2)
                    PRINT*,'Info on the EIGENVALUES OF THE SYSTEM OF MASS STIFFNESS IS:',INFO
                    PRINT*,'THE CPU TIME IT TOOK THE EIGENVALUE ANALYSIS IS:', T2-T1
                    
                    DEALLOCATE(IWORK)
                    DEALLOCATE(WORK)
                    DEALLOCATE(IFAIL)
                    
                END IF                
                
            END IF
            CALL CPU_TIME(T2)
            
            PRINT*, 'THE CPU TIME IT TOOK THE EIGENVALUE ANALYSIS IS:', T2-T1
            DEALLOCATE(int)
            
        
        
        ELSE ! CONDITIONS FOR MAKING THE MASS MATRIX A 2-DIMENSIONAL NON-DIAGONAL MATRIX
            

            !CALL CPU_TIME(T1)
            !! ....... multiplying the inverse of the mass matrix with the stiffness matrix...............
            !N = Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
            !CALL FindHalfBandWidth(n,ReducedGlobalMassMatrix,KD)
            !print*,'KD for the Mass is:',KD
            !CALL InverseBandedSymmetricPositive(N, KD,ReducedGlobalMassMatrix,InvReducedGlobalMassMatrix, INFO )
            !CALL MatrixMultiply(N,N,N,ReducedGlobalStiffnessMatrix,InvReducedGlobalMassMatrix,EigenMatrix)
            !PRINT*,'Info on the inverse of the reduced mass matrix:',INFO
            
            
                        !..........................Finding the EigenValues..............................
            IF(Simdata%OutputReq%EIGEN%ALL) THEN ! ALL EIGENVALUES 
                IF(Simdata%OutputReq%EIGEN%EigenVectors) THEN ! INCLUDE IN THE SEARCH THE EIGENVECTORS
                    N = Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
                                   
                    !CALL elmhes(N,N,1,N,EigenMatrix,int) ! Transform the General Matrix into a Hessenberg
                    !CALL eltran(N,N,1,N,EigenMatrix,int,Zmatrix) ! Find the Transformation matrix useful for finding the EigenVectors
                    !CALL hqr2(N,N,1,N,EigenMatrix,EigenValueR,EigenValueI,Zmatrix,ierr) ! Finds the Eigenvalues and Vectors of a Hessenberg Matrix.
                    !CALL AscendingOrder(N,EigenValueR,index)
                    !print*,'ierr; the error',ierr
                    
                    
                    VL = 0
                    VU = 0
                    IL = 1
                    IU = 10
                    ABSTOL = 2 *DLAMCH('S')
                    LWORK = max(1,8*N)
                    ALLOCATE(WORK(MAX(1,8*N)))
                    ALLOCATE(IFAIL(N))
                    ALLOCATE(IWORK(5*N))
                    CALL CPU_TIME(T1)
                    CALL DSYGVX( 1, 'V', 'A', 'U', N, ReducedGlobalStiffnessMatrix, N, ReducedGlobalMassMatrix, N, &
                          VL, VU, IL, IU, ABSTOL, M, EigenValueR, Zmatrix, N, WORK, LWORK, IWORK, IFAIL, INFO )
                    CALL CPU_TIME(T2)
                    PRINT*,'Info on the EIGENVALUES OF THE SYSTEM OF MASS STIFFNESS IS:',INFO
                    PRINT*,'THE CPU TIME IT TOOK THE EIGENVALUE ANALYSIS IS:', T2-T1
                    
                    DEALLOCATE(IWORK)
                    DEALLOCATE(WORK)
                    DEALLOCATE(IFAIL)
                    
                    
                    
                ELSE ! All EigenValues without EigenVectors
                    N = Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
                    VL = 0
                    VU = 0
                    IL = 1
                    IU = 10
                    ABSTOL = 2 *DLAMCH('S')
                    LWORK = max(1,8*N)
                    ALLOCATE(WORK(MAX(1,8*N)))
                    ALLOCATE(IFAIL(N))
                    ALLOCATE(IWORK(5*N))
                    CALL CPU_TIME(T1)
                    CALL DSYGVX( 1, 'N', 'A', 'U', N, ReducedGlobalStiffnessMatrix, N, ReducedGlobalMassMatrix, N, &
                          VL, VU, IL, IU, ABSTOL, M, EigenValueR, Zmatrix, N, WORK, LWORK, IWORK, IFAIL, INFO )
                    CALL CPU_TIME(T2)
                    PRINT*,'Info on the EIGENVALUES OF THE SYSTEM OF MASS STIFFNESS IS:',INFO
                    PRINT*,'THE CPU TIME IT TOOK THE EIGENVALUE ANALYSIS IS:', T2-T1
                    
                    DEALLOCATE(IWORK)
                    DEALLOCATE(WORK)
                    DEALLOCATE(IFAIL)
                    
                    
                END IF
            
            ELSE ! ONLY A FEW EIGENVALUES WITH THEIR CORRESPONDING EIGENVECTORS ARE REQUESTED 
                IF(Simdata%OutputReq%EIGEN%EigenVectors) THEN ! INCLUDE IN THE SEARCH THE EIGENVECTORS
                    N = Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize                   
                    
                    VL = 0
                    VU = 0
                    IL = Simdata%OutputReq%EIGEN%LowerIndex
                    IU = Simdata%OutputReq%EIGEN%LowerIndex+Simdata%OutputReq%EIGEN%NumberOfEigenValues-1
                    ABSTOL = 2 *DLAMCH('S')
                    LWORK = max(1,8*N)
                    ALLOCATE(WORK(MAX(1,8*N)))
                    ALLOCATE(IFAIL(N))
                    ALLOCATE(IWORK(5*N))
                    CALL CPU_TIME(T1)
                    CALL DSYGVX( 1, 'V', 'I', 'U', N, ReducedGlobalStiffnessMatrix, N, ReducedGlobalMassMatrix, N, &
                          VL, VU, IL, IU, ABSTOL, M, EigenValueR, Zmatrix, N, WORK, LWORK, IWORK, IFAIL, INFO )
                    CALL CPU_TIME(T2)
                    PRINT*,'Info on the EIGENVALUES OF THE SYSTEM OF MASS STIFFNESS IS:',INFO
                    PRINT*,'THE CPU TIME IT TOOK THE EIGENVALUE ANALYSIS IS:', T2-T1
                    
                    DEALLOCATE(IWORK)
                    DEALLOCATE(WORK)
                    DEALLOCATE(IFAIL)
                    
                    
                    
                ELSE ! A few EigenValues without EigenVectors
                     N = Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
                    VL = 0
                    VU = 0
                    IL = Simdata%OutputReq%EIGEN%LowerIndex
                    IU = Simdata%OutputReq%EIGEN%LowerIndex+Simdata%OutputReq%EIGEN%NumberOfEigenValues-1
                    ABSTOL = 2 *DLAMCH('S')
                    LWORK = max(1,8*N)
                    ALLOCATE(WORK(MAX(1,8*N)))
                    ALLOCATE(IFAIL(N))
                    ALLOCATE(IWORK(5*N))
                    CALL CPU_TIME(T1)
                    CALL DSYGVX( 1, 'N', 'I', 'U', N, ReducedGlobalStiffnessMatrix, N, ReducedGlobalMassMatrix, N, &
                          VL, VU, IL, IU, ABSTOL, M, EigenValueR, Zmatrix, N, WORK, LWORK, IWORK, IFAIL, INFO )
                    CALL CPU_TIME(T2)
                    PRINT*,'Info on the EIGENVALUES OF THE SYSTEM OF MASS STIFFNESS IS:',INFO
                    PRINT*,'THE CPU TIME IT TOOK THE EIGENVALUE ANALYSIS IS:', T2-T1
                    
                    DEALLOCATE(IWORK)
                    DEALLOCATE(WORK)
                    DEALLOCATE(IFAIL)
                    
                END IF                
                
                
                
                
                
            END IF
            

            
            
            
            
            
        END IF
        PRINT*,'........................... END EIGENVALUE ANALYSIS..........................'
    END IF
    

    
    
    
    
    
    
    
    
    
    
    !.................................. Post Analyis: Outputting Requested Data............................
    IF(Simdata%SimType%StaticAnalysis) THEN
      PRINT*,'...............Outputting the Displacement Solution..........'       
      !................... Output the displacement solution........................
      ALLOCATE(GlobalNodalDisplacement(Meshdata%TotalNumberOfNodes,2))
     DO p=1,Meshdata%TotalNumberofDOFperNode
      DO K=1,Meshdata%TotalNumberofNodes
      IF (ID(p,K) .NE. 0) THEN
        GlobalNodalDisplacement(K,p) = U(ID(p,K),1)
      ELSE
        GlobalNodalDisplacement(K,p) = BCVec(NegID(p,K))
      END IF
      END DO
     END DO
  
    OPEN(UNIT=18,file=OutputPath//'GlobalNodalDisplacement.txt',form='formatted',access='sequential',status='unknown')
    DO K=1,Meshdata%TotalNumberOfNodes
    WRITE(18,'(2(1X,ES23.16))') GlobalNodalDisplacement(K,:)
    END DO 
    CLOSE(18)
    
    
    
    !...............................................................................................................!
    !....................... Accuracy Analysis .....................................................................!
    !...............................................................................................................!
    
    PP = 1
    R1 = 3
    R2 = 23
    Xc = 1
    Yc = 1
    ALLOCATE(TrueDisp(Meshdata%TotalNumberOfNodes,2))
    ALLOCATE(TStress(Meshdata%TotalNumberOfNodes,3))
    CALL TrueDisplacement(Meshdata%TotalNumberOfNodes,PP,Materialdata%Iso%E,Materialdata%Iso%Nu,Xc,Yc,R1,R2,Meshdata%GlobalCoordinates,TrueDisp)
    CALL TrueStress(Meshdata%TotalNumberOfNodes,PP,Materialdata%Iso%E,Materialdata%Iso%Nu,Xc,Yc,R1,R2,Meshdata%GlobalCoordinates,TStress)
    CALL RelativeError(Meshdata%TotalNumberOfNodes,TrueDisp,GlobalNodalDisplacement,DispRErr)
    
    print*,'the relative error between the True displacement and the FEA solution is:',DispRErr
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !DEALLOCATE(GlobalNodalDisplacement)
    
     !OPEN(UNIT=18,file=path//'ReducedGlobalStiffnessMatrix.txt',form='formatted',access='sequential',status='unknown')
     !OPEN(UNIT=19,file=OutputPath//'ReducedForceExt.txt',form='formatted',access='sequential',status='unknown')
     !WRITE(rowfmt,'(A,I6,A)') '(',Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize,'(1X,ES23.16))'
     !DO I=1,Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
     ! WRITE(18,FMT=ROWFMT) ReducedGlobalStiffnessMatrix(I,:)
     ! WRITE(19,'(1(1X,ES23.16))') ReducedRHS(I)
     !END DO
     !ClOSE(18)
     !ClOSE(19)
     
     PRINT*,'...............END Outputting the Displacement Solution..........'
     
     IF(Simdata%OutputReq%STRESS%OUTPUT) THEN
         PRINT*,'............... OUTPUTTING STRESSES ...............................'
         
        IF(Simdata%OutputReq%STRESS%ALL) THEN
            ALLOCATE(NodalDisplacement(Meshdata%TotalNumberOfNodesPerElement,2))
            ALLOCATE(NodalCoordinates(Meshdata%TotalNumberOfNodesPerElement,2))
            ALLOCATE(StressDistribution(Meshdata%TotalNumberOfNodes,3))
            ALLOCATE(StressDistributionI(Meshdata%TotalNumberOfNodesPerElement,3))
            
            IF(Simdata%SimType%Quadrature == 0) THEN ! This is to output stress distribution at nodal points and not at the quadrature points
                CALL DerivativeLagrangePArray(Meshdata%TotalNumberOfNodesPerEdge,Meshdata%TotalNumberOfNodesPerEdge,Roots,Roots,DL)

                !.................... Compute all Lagrange polynonials at the different integration points...........
                
                CALL LagrangePArray(Meshdata%TotalNumberOfNodesPerEdge,Meshdata%TotalNumberOfNodesPerEdge,Roots,Roots,LMatrix)
            END IF
            
            
            
            DO I =1, Meshdata%TotalNumberOfElements
                ! COPY ELEMENT NODAL DISPLACEMENTS FROM GLOBAL DISPLACEMENT
                
                NodalDisplacement = GlobalNodalDisplacement(Meshdata%ConnectivityMatrix(I,:),:)
                ! COPY ELEMENT NODAL COORDINATES FROM GLOBAL COORDINATES
                NodalCoordinates = Meshdata%GlobalCoordinates(Meshdata%ConnectivityMatrix(I,:),:)
               
                
                IF((Simdata%SimType%SpectralElements) .AND. (Simdata%SimType%Quadrature == 1) .AND. (.NOT. (Simdata%SimType%ReducedQuadrature))) THEN
                    CALL Eff_Compute_Stress_2D_Element(Meshdata%TotalNumberOfNodesPerEdge,DL,NodalCoordinates,NodalDisplacement,Dmatrix,StressDistributionI)
                ELSE
                    CALL  Compute_Stress_2D_Element(Meshdata%TotalNumberOfNodesPerEdge,LMatrix,DL,NodalCoordinates,NodalDisplacement,Dmatrix,StressDistributionI)
                END IF
                StressDistribution(Meshdata%ConnectivityMatrix(I,:),:) = StressDistributionI(:,:)
                
            END DO
            OPEN(UNIT=18,file=OutputPath//'StressDistribution.txt',form='formatted',access='sequential',status='unknown')
            DO I =1,Meshdata%TotalNumberOfNodes
                WRITE(18,'(3ES20.13)') StressDistribution(I,:)
            END DO
            CLOSE(18)
            
            !...................................................................................................!
            !......................... Accuracy Analysis .......................................................!
            !...................................................................................................!
            CALL RelativeSTRESSError(Meshdata%TotalNumberOfNodes,TStress,StressDistribution,StrssRErr)
    
            print*,'the relative error between the True Stress and the FEA solution is:',STrssRErr
            
            OPEN(UNIT=18,file=OutputPath//'AccuracyAnalysis.txt',form='formatted',access='APPEND',status='unknown')
            
            WRITE(18,'(2I10,4ES20.13)') Meshdata%TotalNumberOfNodesPerEdge,Meshdata%TotalNumberOfNodes,TAssembly,TSolver,DispRErr,StrssRErr
            
            CLOSE(18)
            
            
            DEALLOCATE(NodalDisplacement)
            DEALLOCATE(NodalCoordinates)
            DEALLOCATE(StressDistribution)
            DEALLOCATE(StressDistributionI)
            
        END IF
        IF(Simdata%OutputReq%STRESS%IsElementEdge) THEN
            ALLOCATE(NodalDisplacement(Meshdata%TotalNumberOfNodesPerElement,2))
            ALLOCATE(NodalCoordinates(Meshdata%TotalNumberOfNodesPerElement,2))
            ALLOCATE(StressAlongEdge(Meshdata%TotalNumberOfNodesPerEdge,3))
            DO I =1, Simdata%OutputReq%STRESS%ElementEdgeSize
                
                ! COPY ELEMENT NODAL DISPLACEMENTS FROM GLOBAL DISPLACEMENT
                NodalDisplacement = GlobalNodalDisplacement(Meshdata%ConnectivityMatrix(Simdata%OutputReq%STRESS%Element(I),:),:)
                ! COPY ELEMENT NODAL COORDINATES FROM GLOBAL COORDINATES
                NodalCoordinates = Meshdata%GlobalCoordinates(Meshdata%ConnectivityMatrix(Simdata%OutputReq%STRESS%Element(I),:),:)
                    IF(Simdata%OutputReq%STRESS%Edge(I)==1) THEN
                        edge =1
                    ELSE IF(Simdata%OutputReq%STRESS%Edge(I)==2) THEN
                        edge = 3
                    ELSE IF(Simdata%OutputReq%STRESS%Edge(I)==3) THEN
                        edge = 4
                    ELSE IF(Simdata%OutputReq%STRESS%Edge(I)==4) THEN
                        edge = 2
                    END IF
                    !         3 (4)
                    !        ---------                   
                    ! 4 (2)  |       | 2 (3)
                    !        |       |
                    !        ---------
                    !         1 (1)
                                    
                
                IF((Simdata%SimType%SpectralElements) .AND. (Simdata%SimType%Quadrature == 1) .AND. (.NOT. (Simdata%SimType%ReducedQuadrature))) THEN
                    

                    CALL Eff_Compute_Stress_2D_Along_Edge(Meshdata%TotalNumberOfNodesPerEdge,DL,DL_Minus_one,DL_Plus_one,NodalCoordinates, &
                                            NodalDisplacement,Dmatrix,edge,StressAlongEdge)
                ELSE
                    CALL  Compute_Stress_2D_Along_Edge(Meshdata%TotalNumberOfNodesPerEdge,LMatrix,DL,DL_Minus_one,DL_Plus_one,NodalCoordinates, &
                                            NodalDisplacement,Dmatrix,edge,StressAlongEdge)
                END IF
                
                WRITE(ARG,'(A,I4,A,I4,A)') 'StressElement',Simdata%OutputReq%STRESS%Element(I),'Edge',Simdata%OutputReq%STRESS%Edge(I),'.txt'
                OPEN(UNIT=18,file=OutputPath//TRIM(ARG),form='formatted',access='sequential',status='unknown')
                DO K =1,Meshdata%TotalNumberOfNodesPerEdge
                    WRITE(18,'(3ES20.13)') StressAlongEdge(K,:)
                END DO
                CLOSE(18)
                
            END DO
            
            
            DEALLOCATE(NodalDisplacement)
            DEALLOCATE(NodalCoordinates)
            DEALLOCATE(StressAlongEdge)
        END IF
             
             
             
             
         
         
     
     
     
     
     
     
        PRINT*,'............. END OUTPUTTING STRESSES ...........................'
     END IF
     
     
     
     
     
     
     
     
     
     
     
    END IF
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
     
     
     IF((Simdata%SimType%EigenValueAnalysis == .TRUE.)) THEN
     
     
     !...................EigenValues + Vectors.........................................
     PRINT*,'................OUTPUTTING EIVENVALUE DATA..............................'
     !IF((Simdata%SimType%SpectralElements) .AND. (Simdata%SimType%Quadrature == 1) .AND. (.NOT. (Simdata%SimType%ReducedQuadrature))) THEN
     !OPEN(UNIT=18,file=OutputPath//'EigenMatrix.txt',form='formatted',access='sequential',status='unknown')
     !WRITE(rowfmt,'(A,I6,A)') '(',Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize,'(1X,ES23.16))'
     !   DO I=1,Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
     !   WRITE(18,FMT=ROWFMT) ReducedGlobalStiffnessMatrix(I,:)/ReducedGlobalMassSpectral(I)
     !   END DO    
     !CLOSE(18)
     !END IF
     
     IF((Simdata%SimType%SpectralElements) .AND. (Simdata%SimType%Quadrature == 1) .AND. (.NOT. (Simdata%SimType%ReducedQuadrature))) THEN ! CONDITIONS FOR MAKING THE MASS MATRIX DIAGONAL MATRIX
     
     IF((Simdata%OutputReq%EIGEN%EigenVectors) .AND. (Simdata%OutputReq%EIGEN%ALL)) THEN
         
      N=   Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
     ALLOCATE(GlobalNodalEigenVectorsX(Meshdata%TotalNumberOfNodes,N))
     ALLOCATE(GlobalNodalEigenVectorsY(Meshdata%TotalNumberOfNodes,N))
     WRITE(rowfmt,'(A,I6,A)') '(',Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize,'(1X,ES23.16))'
     OPEN(UNIT=18,file=OutputPath//'GlobalNodalEigenVectorsX.txt',form='formatted',access='sequential',status='unknown')
     OPEN(UNIT=19,file=OutputPath//'GlobalNodalEigenVectorsY.txt',form='formatted',access='sequential',status='unknown')
      DO I=1,Meshdata%TotalNumberOfNodes
      DO K=1,Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
      IF (ID(1,I) .NE. 0) THEN
          GlobalNodalEigenVectorsX(I,K) = Zmatrix(ID(1,I),index(K))
      ELSE
        GlobalNodalEigenVectorsX(I,K) = BCVec(NegID(1,I))
      END IF
      IF (ID(2,I) .NE. 0) THEN
          GlobalNodalEigenVectorsY(I,K) = Zmatrix(ID(2,I),index(K))
      ELSE
        GlobalNodalEigenVectorsY(I,K) = BCVec(NegID(2,I))
      END IF
      END DO
      WRITE(18,FMT=ROWFMT) GlobalNodalEigenVectorsX(I,:)
      WRITE(19,FMT=ROWFMT) GlobalNodalEigenVectorsY(I,:)
     END DO
     DEALLOCATE(GlobalNodalEigenVectorsX)
     DEALLOCATE(GlobalNodalEigenVectorsY)
     CLOSE(18)
     CLOSE(19)
     ELSE IF (Simdata%OutputReq%EIGEN%EigenVectors) THEN ! A FEW OF EIGENVECTORS ARE REQUESTED

     M=   Simdata%OutputReq%EIGEN%NumberOfEigenValues
     ALLOCATE(GlobalNodalEigenVectorsX(Meshdata%TotalNumberOfNodes,M))
     ALLOCATE(GlobalNodalEigenVectorsY(Meshdata%TotalNumberOfNodes,M))
     WRITE(rowfmt,'(A,I6,A)') '(',M,'(1X,ES23.16))'
     OPEN(UNIT=18,file=OutputPath//'GlobalNodalEigenVectorsX.txt',form='formatted',access='sequential',status='unknown')
     OPEN(UNIT=19,file=OutputPath//'GlobalNodalEigenVectorsY.txt',form='formatted',access='sequential',status='unknown')
      DO I=1,Meshdata%TotalNumberOfNodes
      DO K=1,M
      IF (ID(1,I) .NE. 0) THEN
          GlobalNodalEigenVectorsX(I,K) = Zmatrix(ID(1,I),K)
      ELSE
        GlobalNodalEigenVectorsX(I,K) = BCVec(NegID(1,I))
      END IF
      IF (ID(2,I) .NE. 0) THEN
          GlobalNodalEigenVectorsY(I,K) = Zmatrix(ID(2,I),K)
      ELSE
        GlobalNodalEigenVectorsY(I,K) = BCVec(NegID(2,I))
      END IF
      END DO
      WRITE(18,FMT=ROWFMT) GlobalNodalEigenVectorsX(I,:)
      WRITE(19,FMT=ROWFMT) GlobalNodalEigenVectorsY(I,:)
     END DO
     DEALLOCATE(GlobalNodalEigenVectorsX)
     DEALLOCATE(GlobalNodalEigenVectorsY)
     CLOSE(18)
     CLOSE(19)     
     END IF
     
     
     IF( (Simdata%OutputReq%EIGEN%ALL)) THEN
     OPEN(UNIT=18,file=OutputPath//'EigenValues.txt',form='formatted',access='sequential',status='unknown')
     DO I=1,Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
         WRITE(18,*) EigenValueR(index(I))
     END DO
     CLOSE(18)
     ELSE
     OPEN(UNIT=18,file=OutputPath//'EigenValues.txt',form='formatted',access='sequential',status='unknown')
     DO I=1,Simdata%OutputReq%EIGEN%NumberOfEigenValues
         WRITE(18,*) EigenValueR(I)
     END DO
     CLOSE(18)    
     END IF
     
     
     ELSE ! GENERAL MASS MATRIX AND GENERALIZED EIGENVALUE ANALYSIS.
     IF((Simdata%OutputReq%EIGEN%EigenVectors) .AND. (Simdata%OutputReq%EIGEN%ALL)) THEN
         
      N=   Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
     ALLOCATE(GlobalNodalEigenVectorsX(Meshdata%TotalNumberOfNodes,N))
     ALLOCATE(GlobalNodalEigenVectorsY(Meshdata%TotalNumberOfNodes,N))
     WRITE(rowfmt,'(A,I6,A)') '(',N,'(1X,ES23.16))'
     OPEN(UNIT=18,file=OutputPath//'GlobalNodalEigenVectorsX.txt',form='formatted',access='sequential',status='unknown')
     OPEN(UNIT=19,file=OutputPath//'GlobalNodalEigenVectorsY.txt',form='formatted',access='sequential',status='unknown')
      DO I=1,Meshdata%TotalNumberOfNodes
      DO K=1,N
      IF (ID(1,I) .NE. 0) THEN
          GlobalNodalEigenVectorsX(I,K) = Zmatrix(ID(1,I),K)
      ELSE
        GlobalNodalEigenVectorsX(I,K) = BCVec(NegID(1,I))
      END IF
      IF (ID(2,I) .NE. 0) THEN
          GlobalNodalEigenVectorsY(I,K) = Zmatrix(ID(2,I),K)
      ELSE
        GlobalNodalEigenVectorsY(I,K) = BCVec(NegID(2,I))
      END IF
      END DO
      WRITE(18,FMT=ROWFMT) GlobalNodalEigenVectorsX(I,:)
      WRITE(19,FMT=ROWFMT) GlobalNodalEigenVectorsY(I,:)
     END DO
     DEALLOCATE(GlobalNodalEigenVectorsX)
     DEALLOCATE(GlobalNodalEigenVectorsY)
     CLOSE(18)
     CLOSE(19)
     ELSE IF (Simdata%OutputReq%EIGEN%EigenVectors) THEN ! A FEW OF EIGENVECTORS ARE REQUESTED
     M=   Simdata%OutputReq%EIGEN%NumberOfEigenValues
     ALLOCATE(GlobalNodalEigenVectorsX(Meshdata%TotalNumberOfNodes,M))
     ALLOCATE(GlobalNodalEigenVectorsY(Meshdata%TotalNumberOfNodes,M))
     WRITE(rowfmt,'(A,I6,A)') '(',M,'(1X,ES23.16))'
     OPEN(UNIT=18,file=OutputPath//'GlobalNodalEigenVectorsX.txt',form='formatted',access='sequential',status='unknown')
     OPEN(UNIT=19,file=OutputPath//'GlobalNodalEigenVectorsY.txt',form='formatted',access='sequential',status='unknown')
      DO I=1,Meshdata%TotalNumberOfNodes
      DO K=1,M
      IF (ID(1,I) .NE. 0) THEN
          GlobalNodalEigenVectorsX(I,K) = Zmatrix(ID(1,I),K)
      ELSE
        GlobalNodalEigenVectorsX(I,K) = BCVec(NegID(1,I))
      END IF
      IF (ID(2,I) .NE. 0) THEN
          GlobalNodalEigenVectorsY(I,K) = Zmatrix(ID(2,I),K)
      ELSE
        GlobalNodalEigenVectorsY(I,K) = BCVec(NegID(2,I))
      END IF
      END DO
      WRITE(18,FMT=ROWFMT) GlobalNodalEigenVectorsX(I,:)
      WRITE(19,FMT=ROWFMT) GlobalNodalEigenVectorsY(I,:)
     END DO
     DEALLOCATE(GlobalNodalEigenVectorsX)
     DEALLOCATE(GlobalNodalEigenVectorsY)
     CLOSE(18)
     CLOSE(19)
     END IF
     
     
     IF( (Simdata%OutputReq%EIGEN%ALL)) THEN
         N =Meshdata%TotalNumberofDOFperNode*Meshdata%TotalNumberofNodes-BCdata%BCSize
     OPEN(UNIT=18,file=OutputPath//'EigenValues.txt',form='formatted',access='sequential',status='unknown')
     DO I=1,N
         WRITE(18,*) EigenValueR(I)
     END DO
     CLOSE(18)
     ELSE
     M=   Simdata%OutputReq%EIGEN%NumberOfEigenValues    
     OPEN(UNIT=18,file=OutputPath//'EigenValues.txt',form='formatted',access='sequential',status='unknown')
     DO I=1,M
         WRITE(18,*) EigenValueR(I)
     END DO
     CLOSE(18)    
     END IF         
         
         
         
     END IF
      
     PRINT*,'................END OUTPUTTING EIVENVALUE DATA..............................'
     
     END IF
     
     
     
     
     
    !OPEN(UNIT=19,file=OutputPath//'DLMatrix.txt',form='formatted',access='sequential',status='unknown')
    ! WRITE(rowfmt,'(A,I4,A)') '(',Meshdata%TotalNumberofNodesPerEdge,'(1X,ES23.16))'
    ! DO I=1,Meshdata%TotalNumberofNodesPerEdge
    !  WRITE(19,rowfmt) DL(I,:)
    ! END DO
    ! CLOSE(19)
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    !-------------------------------------------- Clean up ------------------------------------
    IF(ALLOCATED(GlobalNodalDisplacement))  DEALLOCATE(GlobalNodalDisplacement)
    IF(ALLOCATED(ReducedGlobalStiffnessMatrix)) DEALLOCATE(ReducedGlobalStiffnessMatrix)
    IF(ALLOCATED(ReducedGlobalMassMatrix)) DEALLOCATE(ReducedGlobalMassMatrix)
    IF(ALLOCATED(ReducedGlobalMassSpectral)) DEALLOCATE(ReducedGlobalMassSpectral)
    IF(ALLOCATED(ZMatrix)) DEALLOCATE(ZMatrix)    
    
    
    
    
    CALL DestroyMeshData(Meshdata)
    CALL DestroyMaterialData(Materialdata)
    CALL DestroyLoadData(Loaddata)
    CALL DestroyBCData(BCdata)
    CALL DestroySimData(Simdata)
    
    
    STOP
END PROGRAM MAIN
    
    