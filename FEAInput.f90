!=======================================================================================================    
!=======================================================================================================      
!======================================================================================================= 
    MODULE IO
    ! Amir Gasmi 11/14/2016
    USE PrecisionKinds
    IMPLICIT NONE
    
    CONTAINS
    !-----------------------------------------------
    SUBROUTINE FindKeyword(LINE,KEYWORD)
    
    CHARACTER(:), INTENT(IN),ALLOCATABLE :: LINE
    CHARACTER(:),INTENT(OUT),ALLOCATABLE :: KEYWORD
    INTEGER(IntKi) :: POS1,POS2
    INTEGER(IntKi) :: I
    POS1 = 0
    POS2 = 0
    KEYWORD =''
    
    DO I =1,LEN(LINE)-1
        IF( POS1 == 0 ) THEN
            IF (LINE(I:I+1) .EQ. '##') THEN
                POS1 = I
            END IF
        ELSE
            IF(POS2 == 0) THEN
                IF (LINE(I:I+1) .EQ. '##') THEN
                    POS2 = I
                END IF               
            END IF
        END IF
        
    END DO
    IF((POS1>0) .AND. (POS2>POS1))  KEYWORD = LINE(POS1+2:POS2-1)
    
    
    RETURN
    END SUBROUTINE FindKeyword
    !---------------------------------------------
    SUBROUTINE GetLine(unit,line,IOST)
    INTEGER(IntKi),INTENT(IN) :: unit
    CHARACTER(:),INTENT(OUT),ALLOCATABLE :: line
    INTEGER(IntKi),    INTENT(OUT) :: IOST
    CHARACTER(256) :: BUFFER
    INTEGER(IntKi) :: SIZE
    LINE =''
    DO
        READ(UNIT,'(A)',ADVANCE='NO',SIZE=SIZE,IOSTAT=IOST) BUFFER
        IF (IOST>0) THEN
            EXIT
        END IF
        
        LINE = LINE // BUFFER(:SIZE)
        IF(IOST <0) THEN 
            EXIT
        END IF
        
        
    END DO
    RETURN
    END SUBROUTINE GetLine
        !---------------------------------------------
    SUBROUTINE ReadNextRecord(UNIT,VAR)
    INTEGER(IntKi), INTENT(IN) :: UNIT
    CHARACTER(:), ALLOCATABLE :: VAR
    
    CHARACTER(LEN=400) :: TMP
    
    READ(UNIT,*) TMP
    VAR = TRIM(TMP)
    
    RETURN
    END SUBROUTINE ReadNextRecord
    END MODULE IO
!=======================================================================================================    
!=======================================================================================================      
!======================================================================================================= 
    MODULE FEAInput 
    ! Amir Gasmi 11/14/2016
    USE IO
    USE ISO_FORTRAN_ENV,ONLY : IOSTAT_EOR,IOSTAT_END
    IMPLICIT NONE
!-----------------------------------------------------------------------  
! Mesh Data
    TYPE, PUBLIC :: Mesh
        CHARACTER(:),ALLOCATABLE :: MeshFileName
        INTEGER(IntKi) :: TotalNumberOfElements
        INTEGER(IntKi) :: TotalNumberOfNodes
        INTEGER(IntKi) :: TotalNumberOfNodesPerElement
        INTEGER(IntKi) :: TotalNumberOfDOFPerNode
        INTEGER(IntKi) :: TotalNumberOfNodesPerEdge
          
        CHARACTER(:),ALLOCATABLE :: ConnectivityMatrixFileName
        INTEGER(IntKi), ALLOCATABLE, DIMENSION(:,:) :: ConnectivityMatrix
        
        CHARACTER(:),ALLOCATABLE :: CoordinateMatrixFileName
        REAL(ReKi), ALLOCATABLE, DIMENSION(:,:) :: GlobalCoordinates
        
        CHARACTER(:),ALLOCATABLE :: EdgeNodesFileName
        INTEGER(IntKi), ALLOCATABLE, DIMENSION(:,:) :: EdgeNodes
        
    END TYPE Mesh
    
!----------------------------------------------------------------------- 
! Material Data
    TYPE, PUBLIC ::Isotropic
        REAL(ReKI) :: E
        REAL(ReKI) :: nu
    END TYPE Isotropic
!-----------------------------------------------------------------------     
    TYPE, PUBLIC ::Anisotropic
        REAL(ReKI), DIMENSION(3,3) :: D
    END TYPE Anisotropic
!-----------------------------------------------------------------------     
    
    TYPE, PUBLIC :: MaterialSet
        INTEGER(IntKi) :: NodesSize
        INTEGER(IntKi), DIMENSION(:),ALLOCATABLE :: Nodes
        LOGICAL  :: Isotropic ! is the material isotropic nor
        TYPE(Isotropic) :: Iso
        TYPE(Anisotropic) :: Aniso
        REAL(ReKi) :: Rho  ! Mass density
    END TYPE MaterialSet
!----------------------------------------------------------------------- 
    
    TYPE, PUBLIC :: Material
        CHARACTER(:),ALLOCATABLE :: MaterialFileName
        LOGICAL  :: Isotropic ! is the material isotropic
        LOGICAL  :: Homogenous ! is the material uniform or distributed
        TYPE(Isotropic) :: Iso
        TYPE(Anisotropic) :: Aniso
        REAL(ReKi) :: Rho   ! Mass density
        INTEGER(IntKi)  ::MaterialSetsSize 
        TYPE(MaterialSet), DIMENSION(:), ALLOCATABLE :: MaterialSets
    END TYPE Material
!----------------------------------------------------------------------- 
! Loading Data

    
    TYPE, PUBLIC :: SurfaceLoading
        INTEGER(IntKi)              :: SurfSize
        INTEGER(IntKi), ALLOCATABLE ::  Node(:)
        INTEGER(IntKi), ALLOCATABLE ::  DOF(:)
        INTEGER(IntKi), ALLOCATABLE  ::  Edge(:) 
        INTEGER(IntKi), ALLOCATABLE  ::  Element(:) 
        REAL(ReKi), ALLOCATABLE     ::  Value(:)        
    END TYPE SurfaceLoading
!----------------------------------------------------------------------- 
    TYPE, PUBLIC :: BodyLoading
        INTEGER(IntKi)              :: BodySize
        INTEGER(IntKi), ALLOCATABLE ::  Node(:)
        INTEGER(IntKi), ALLOCATABLE ::  DOF(:)
        REAL(ReKi), ALLOCATABLE     ::  Element(:) 
        REAL(ReKi), ALLOCATABLE     ::  Value(:)        
    END TYPE BodyLoading
!----------------------------------------------------------------------- 
    TYPE, PUBLIC :: Loads
        CHARACTER(:),ALLOCATABLE :: LoadsFileName
        TYPE(SurfaceLoading) :: SurfLoads
        TYPE(BodyLoading)    :: BodyLoads
    END TYPE Loads
    
    
 !----------------------------------------------------------------------- 
 ! BC Data   
    TYPE, PUBLIC :: BCDOF
        CHARACTER(:),ALLOCATABLE    :: BCFileName
        INTEGER(IntKi)              :: BCSize 
        INTEGER(IntKi), ALLOCATABLE ::  Node(:)
        INTEGER(IntKi), ALLOCATABLE ::  DOF(:)
        REAL(ReKi), ALLOCATABLE     ::  Value(:)       
    END TYPE BCDOF
    !-----------------------------------------------------------------------     
    ! Simulation Data
    
    TYPE, PUBLIC :: SimulationType
        LOGICAL :: PlaneStress ! is the simulation Plane stress when true, flase simulation is Plane strain
        LOGICAL :: StaticAnalysis ! is static analysis included
        LOGICAL :: EigenValueAnalysis ! is Eigenvalue analysis included
        LOGICAL :: SpectralElements ! is it spectral or regular elements ( roots are Legendre-Lobatto or uniformly spaced)
        LOGICAL :: ReducedQuadrature ! is it full integration or reduced ( quadrature points N or N-1)
        INTEGER(Intki) :: Quadrature ! 0 = Legendre-Gauss 1 = Legendre-Lobatto
        
    END TYPE SimulationType
   !-----------------------------------------------------------------------   
    TYPE, PUBLIC :: STRESS  ! THIS IS STRESSES OR STRAINS
        LOGICAL OUTPUT ! any stresses at all you want to output?
        LOGICAL ALL   ! All stresses?
        LOGICAL IsNodeSet ! Whether you want NodeSets to be outputted
        LOGICAL IsElementEdge ! whether you want certain edges to be outputed
        INTEGER(IntKi)  :: NodeSetSize
        INTEGER(IntKi)  :: ElementEdgeSize
        INTEGER(IntKi), DIMENSION(:), ALLOCATABLE :: NodeSet
        INTEGER(IntKi), DIMENSION(:), ALLOCATABLE :: Element
        INTEGER(IntKi), DIMENSION(:), ALLOCATABLE :: Edge  
    END TYPE STRESS
    
   TYPE, PUBLIC :: EIGEN
       LOGICAL :: ALL ! ALL EIGENVALUES?
       LOGICAL :: EigenVectors ! DO YOU WANT ALSO EIGENVECTORS?
       INTEGER(IntKi) :: NumberOfEigenValues
       INTEGER(IntKi) :: LowerIndex
   END TYPE EIGEN
   
    
    TYPE, PUBLIC :: OutputRequest
        TYPE(EIGEN)  :: EIGEN
        TYPE(STRESS) :: STRAIN
        TYPE(STRESS) :: STRESS
    END TYPE OutputRequest
  !-----------------------------------------------------------------------  
    TYPE, PUBLIC :: SimulationData
        CHARACTER(:),ALLOCATABLE :: SimFileName
        TYPE(SimulationType)     :: SimType
        Type(OutputRequest)      :: OutputReq
    END TYPE SimulationData
      
    
       
    CONTAINS
     !------------------------------------------------------------------------------
    SUBROUTINE PopulateInputFileNames(datapath,Listfilename,Meshdata,Simdata,Materialdata,Loaddata,BCdata)
    CHARACTER(:),ALLOCATABLE, INTENT(IN) :: datapath
    CHARACTER(:),ALLOCATABLE, INTENT(IN) :: Listfilename
    TYPE(Mesh), INTENT(INOUT)            :: Meshdata
    TYPE(SimulationData), INTENT(INOUT)  :: Simdata
    TYPE(Material), INTENT(INOUT)        :: Materialdata
    TYPE(Loads), INTENT(INOUT)           :: Loaddata
    TYPE(BCDOF), INTENT(INOUT)           :: BCdata
    
    
    CHARACTER(:),ALLOCATABLE :: FullPathName
    CHARACTER(:),ALLOCATABLE :: LINE, KEYWORD
    CHARACTER(LEN=400)       :: TMP
    INTEGER(IntKi) :: IOSTAT, I
    INTEGER(IntKi) :: UNIT = 18
    
    ! ------------------------------------ begin READING List of File names --------------------------------------
    FullPathName = trim(datapath) // trim(Listfilename)
    OPEN(UNIT=UNIT,FILE=FullPathName,FORM='FORMATTED',ACCESS='SEQUENTIAL',ACTION='READ',STATUS='OLD',IOSTAT=IOSTAT)
    LINE =''
    KEYWORD =''
    DO
        CALL GetLine(UNIT,LINE,IOSTAT)
        IF(IOSTAT == IOSTAT_END) EXIT
        CALL FindKeyword(LINE,KEYWORD)
        SELECT CASE(KEYWORD)
           
        CASE ('MeshFileName')
            CALL ReadNextRecord(UNIT,Meshdata%MeshFileName)
        CASE('SimFileName')
            CALL ReadNextRecord(UNIT,Simdata%SimFileName)
        CASE('MaterialFileName')
            CALL ReadNextRecord(UNIT, Materialdata%MaterialFileName)
        CASE('LoadsFileName')
            CALL ReadNextRecord(UNIT,Loaddata%LoadsFileName)
        CASE('BCFileName')
            CALL ReadNextRecord(UNIT,BCdata%BCFileName)
            
            
        END SELECT
    END DO
    
    
    
    
    RETURN
    END SUBROUTINE PopulateInputFileNames
 !------------------------------------------------------------------------------   
    SUBROUTINE PopulateMeshData(datapath, Meshdata)
    CHARACTER(:),ALLOCATABLE, INTENT(IN) :: datapath
    TYPE(Mesh), INTENT(INOUT) :: Meshdata
    CHARACTER(:),ALLOCATABLE :: FullPathName
    CHARACTER(:),ALLOCATABLE :: LINE, KEYWORD
    INTEGER(IntKi) :: IOSTAT, I
    INTEGER(IntKi) :: UNIT = 18
   
    ! ------------------------------------ begin Mesh input data --------------------------------------
    FullPathName = trim(datapath) // trim(Meshdata%MeshFileName)
    OPEN(UNIT=UNIT,FILE=FullPathName,FORM='FORMATTED',ACCESS='SEQUENTIAL',ACTION='READ',STATUS='OLD')
    
    DO
        CALL GetLine(UNIT,LINE,IOSTAT)
        IF (IOSTAT ==IOSTAT_END) EXIT
        CALL FindKeyword(LINE,KEYWORD)
        SELECT CASE(KEYWORD)
        CASE ('TotalNumberOfElements')
            READ(UNIT,*) Meshdata%TotalNumberOfElements
        CASE ('TotalNumberOfNodes')
            READ(UNIT,*) Meshdata%TotalNumberOfNodes
        CASE ('TotalNumberOfNodesPerElement')
            READ(UNIT,*) Meshdata%TotalNumberOfNodesPerElement
        CASE ('TotalNumberOfDOFPerNode')
            READ(UNIT,*) Meshdata%TotalNumberOfDOFPerNode
        CASE ('TotalNumberOfNodesPerEdge')
            READ(UNIT,*) Meshdata%TotalNumberOfNodesPerEdge
        CASE ('EdgeNodesFileName')
            CALL ReadNextRecord(UNIT, Meshdata%EdgeNodesFileName)
        CASE ('ConnectivityMatrixFileName')
            CALL ReadNextRecord(UNIT, Meshdata%ConnectivityMatrixFileName)
        CASE ('CoordinateMatrixFileName')
            CALL ReadNextRecord(UNIT, Meshdata%CoordinateMatrixFileName)
        END SELECT      
    END DO
    CLOSE(UNIT)
    
    
     !-------------------------- Reading the Connectivity Data
    FullPathName = trim(datapath) // trim(Meshdata%ConnectivityMatrixFileName)
    OPEN(UNIT=UNIT,FILE=FullPathName,FORM='FORMATTED',ACCESS='SEQUENTIAL',ACTION='READ',STATUS='OLD')
    ALLOCATE(Meshdata%ConnectivityMatrix(Meshdata%TotalNumberOfElements,Meshdata%TotalNumberOfNodesPerElement))
    DO I =1,Meshdata%TotalNumberofElements
    READ(UNIT,*) Meshdata%ConnectivityMatrix(I,:)
    END DO
    ClOSE(UNIT)

    !----------------------------- Reading the Global Coordinates associated to the Connectivity
    FullPathName = trim(datapath) // trim(Meshdata%CoordinateMatrixFileName)
    OPEN(UNIT=UNIT,FILE=FullPathName,FORM='FORMATTED',ACCESS='SEQUENTIAL',ACTION='READ',STATUS='OLD')
    ALLOCATE(Meshdata%GlobalCoordinates(Meshdata%TotalNumberOfNodes,Meshdata%TotalNumberOfDOFPerNode))
    DO I=1,Meshdata%TotalNumberOfNodes
    READ(UNIT,*) Meshdata%GlobalCoordinates(I,:)
    END DO 
    CLOSE(UNIT)
  
  !----------------------------- Reading the Edge nodes
    FullPathName = trim(datapath) // trim(Meshdata%EdgeNodesFileName)
    OPEN(UNIT=UNIT,FILE=FullPathName,FORM='FORMATTED',ACCESS='SEQUENTIAL',ACTION='READ',STATUS='OLD')
    ALLOCATE(Meshdata%EdgeNodes(Meshdata%TotalNumberOfElements*Meshdata%TotalNumberOfNodesPerEdge,4))
    DO I=1,Meshdata%TotalNumberOfElements*Meshdata%TotalNumberOfNodesPerEdge
    READ(UNIT,*) Meshdata%EdgeNodes(I,:)
    END DO 
    CLOSE(UNIT)
  
  

    
    RETURN
    END SUBROUTINE PopulateMeshData
 !------------------------------------------------------------------------------
    SUBROUTINE PopulateSimulationData(datapath, Simdata)
    CHARACTER(:),ALLOCATABLE, INTENT(IN) :: datapath
    TYPE(SimulationData), INTENT(INOUT) :: Simdata
    CHARACTER(:),ALLOCATABLE :: FullPathName
    CHARACTER(:),ALLOCATABLE :: LINE, KEYWORD
    INTEGER(IntKi) :: IOSTAT, I
    INTEGER(IntKi) :: UNIT = 18
    
    ! ------------------------------------ begin Reading Simulation data --------------------------------------
    FullPathName = trim(datapath) // trim(Simdata%SimFileName)
    OPEN(UNIT=UNIT,FILE=FullPathName,FORM='FORMATTED',ACCESS='SEQUENTIAL',ACTION='READ',STATUS='OLD')    
    
    !-------- Read Simtype data --------
    DO
        CALL GetLine(UNIT,LINE,IOSTAT)
        IF ( IOSTAT == IOSTAT_END) EXIT
        CALL FindKeyword(LINE,KEYWORD)
        SELECT CASE (KEYWORD)
        CASE('PlaneStress')
            READ(UNIT,*) Simdata%SimType%PlaneStress
        CASE('StaticAnalysis')
            READ(UNIT,*) Simdata%SimType%StaticAnalysis
        CASE('EigenValueAnalysis')
            READ(UNIT,*) Simdata%SimType%EigenValueAnalysis
        CASE('SpectralElements')
            READ(UNIT,*) Simdata%SimType%SpectralElements
        CASE('ReducedQuadrature')
            READ(UNIT,*) Simdata%SimType%ReducedQuadrature
        CASE('Quadrature')
            READ(UNIT,*) Simdata%SimType%Quadrature
        CASE('EIGEN')
            DO
                CALL GetLine(UNIT,LINE,IOSTAT)
                IF ( IOSTAT == IOSTAT_END) EXIT
                CALL FindKeyword(LINE,KEYWORD)                
                SELECT CASE (KEYWORD)
                CASE('ALL')
                    READ(UNIT,*) Simdata%OutputReq%EIGEN%ALL
                CASE('EigenVectors')
                    READ(UNIT,*) Simdata%OutputReq%EIGEN%EigenVectors
                CASE('NumberOfEigenValues')
                    READ(UNIT,*) Simdata%OutputReq%EIGEN%NumberOfEigenValues
                CASE('LowerIndex')
                    READ(UNIT,*) Simdata%OutputReq%EIGEN%LowerIndex
                CASE DEFAULT
                    BACKSPACE(UNIT)
                    EXIT            
                END SELECT
            END DO
        CASE('STRAIN')
            DO
                CALL GetLine(UNIT,LINE,IOSTAT)
                IF ( IOSTAT == IOSTAT_END) EXIT
                CALL FindKeyword(LINE,KEYWORD)                
                SELECT CASE (KEYWORD)
                CASE('OUTPUT')
                    READ(UNIT,*) Simdata%OutputReq%STRAIN%OUTPUT
                CASE('ALL')
                    READ(UNIT,*) Simdata%OutputReq%STRAIN%ALL
                CASE('IsNodeSet')
                    READ(UNIT,*) Simdata%OutputReq%STRAIN%IsNodeSet
                CASE('IsElementEdge')
                    READ(UNIT,*) Simdata%OutputReq%STRAIN%IsElementEdge
                CASE('NodeSetSize')
                    READ(UNIT,*) Simdata%OutputReq%STRAIN%NodeSetSize
                CASE('ElementEdgeSize')
                    READ(UNIT,*) Simdata%OutputReq%STRAIN%ElementEdgeSize
                CASE('NodeSet')
                    ALLOCATE(Simdata%OutputReq%STRAIN%NodeSet(Simdata%OutputReq%STRAIN%NodeSetSize))
                    DO I=1,Simdata%OutputReq%STRAIN%NodeSetSize
                        READ(UNIT,*) Simdata%OutputReq%STRAIN%NodeSet(I)
                    END DO                    
                CASE('Element')
                    ALLOCATE(Simdata%OutputReq%STRAIN%Element(Simdata%OutputReq%STRAIN%ElementEdgeSize))
                    DO I=1,Simdata%OutputReq%STRAIN%ElementEdgeSize
                        READ(UNIT,*) Simdata%OutputReq%STRAIN%Element(I)
                    END DO
                CASE('Edge')
                    ALLOCATE(Simdata%OutputReq%STRAIN%Edge(Simdata%OutputReq%STRAIN%ElementEdgeSize))
                    DO I=1,Simdata%OutputReq%STRAIN%ElementEdgeSize
                        READ(UNIT,*) Simdata%OutputReq%STRAIN%Edge(I)
                    END DO
                CASE DEFAULT
                    BACKSPACE(UNIT)
                    EXIT            
                END SELECT
                
            END DO
        CASE('STRESS')
            DO
                CALL GetLine(UNIT,LINE,IOSTAT)
                IF ( IOSTAT == IOSTAT_END) EXIT
                CALL FindKeyword(LINE,KEYWORD)
                SELECT CASE (KEYWORD)
                CASE('OUTPUT')
                    READ(UNIT,*) Simdata%OutputReq%STRESS%OUTPUT
                CASE('ALL')
                    READ(UNIT,*) Simdata%OutputReq%STRESS%ALL
                CASE('IsNodeSet')
                    READ(UNIT,*) Simdata%OutputReq%STRESS%IsNodeSet
                CASE('IsElementEdge')
                    READ(UNIT,*) Simdata%OutputReq%STRESS%IsElementEdge
                CASE('NodeSetSize')
                    READ(UNIT,*) Simdata%OutputReq%STRESS%NodeSetSize
                CASE('ElementEdgeSize')
                    READ(UNIT,*) Simdata%OutputReq%STRESS%ElementEdgeSize
                CASE('NodeSet')
                    ALLOCATE(Simdata%OutputReq%STRESS%NodeSet(Simdata%OutputReq%STRESS%NodeSetSize))
                    DO I=1,Simdata%OutputReq%STRESS%NodeSetSize
                        READ(UNIT,*) Simdata%OutputReq%STRESS%NodeSet(I)
                    END DO                    
                CASE('Element')
                    ALLOCATE(Simdata%OutputReq%STRESS%Element(Simdata%OutputReq%STRESS%ElementEdgeSize))
                    DO I=1,Simdata%OutputReq%STRESS%ElementEdgeSize
                        READ(UNIT,*) Simdata%OutputReq%STRESS%Element(I)
                    END DO
                CASE('Edge')
                    ALLOCATE(Simdata%OutputReq%STRESS%Edge(Simdata%OutputReq%STRESS%ElementEdgeSize))
                    DO I=1,Simdata%OutputReq%STRESS%ElementEdgeSize
                        READ(UNIT,*) Simdata%OutputReq%STRESS%Edge(I)
                    END DO                    
                CASE DEFAULT
                    BACKSPACE(UNIT)
                    EXIT                    
                END SELECT
                
            END DO
            
        END SELECT
        
        
    END DO
    
    
    RETURN
    END SUBROUTINE PopulateSimulationData
 
 !------------------------------------------------------------------------------   
    SUBROUTINE PopulateMaterialData(datapath,Materialdata)
    CHARACTER(:),ALLOCATABLE, INTENT(IN) :: datapath
    TYPE(Material), INTENT(INOUT) :: Materialdata
    CHARACTER(:),ALLOCATABLE :: FullPathName
    CHARACTER(:),ALLOCATABLE :: LINE, KEYWORD
    INTEGER(IntKi) :: IOSTAT, I,K
    INTEGER(IntKi) :: UNIT = 18
   
    ! ------------------------------------ begin Reading Material data --------------------------------------
    FullPathName = trim(datapath) // trim(Materialdata%MaterialFileName)
    OPEN(UNIT=UNIT,FILE=FullPathName,FORM='FORMATTED',ACCESS='SEQUENTIAL',ACTION='READ',STATUS='OLD')
    DO
        CALL GetLine(UNIT,LINE,IOSTAT)
        IF (IOSTAT == IOSTAT_END) RETURN
        CALL FindKeyword(LINE,KEYWORD)
        SELECT CASE(KEYWORD)
        CASE('Isotropic')
            READ(UNIT,*) Materialdata%Isotropic
        CASE('Homogenous')
            READ(UNIT,*) Materialdata%Homogenous
        CASE('Iso')
            DO
                CALL GetLine(UNIT,LINE,IOSTAT)
                IF(IOSTAT == IOSTAT_END) RETURN
                CALL FindKeyword(LINE,KEYWORD)
                SELECT CASE(KEYWORD)
                CASE('E')
                    READ(UNIT,*) Materialdata%Iso%E
                CASE('nu')
                    READ(UNIT,*) Materialdata%Iso%nu
                CASE DEFAULT
                    BACKSPACE(UNIT)
                    EXIT
                END SELECT
            END DO
        CASE('Aniso')
            DO
                CALL GetLine(UNIT,LINE,IOSTAT)
                IF(IOSTAT == IOSTAT_END) RETURN
                CALL FindKeyword(LINE,KEYWORD)
                SELECT CASE(KEYWORD)
                CASE('D')
                    DO I=1,3
                        READ(UNIT,*) Materialdata%Aniso%D(I,:)
                    END DO
                CASE DEFAULT
                    BACKSPACE(UNIT)
                    EXIT
                END SELECT
            END DO
        CASE('Rho')
            READ(UNIT,*) Materialdata%Rho
        CASE('MaterialSetsSize')
            READ(UNIT,*) Materialdata%MaterialSetsSize
        CASE('MaterialSets')
            ALLOCATE(Materialdata%MaterialSets(Materialdata%MaterialSetsSize))
            DO I=1,Materialdata%MaterialSetsSize
                DO
                    
                    CALL GetLine(UNIT,LINE,IOSTAT)
                    CALL FindKeyword(LINE,KEYWORD)
                    SELECT CASE(KEYWORD)
                    CASE('NodesSize')
                        READ(UNIT,*) Materialdata%MaterialSets(I)%NodesSize
                    CASE('Nodes')
                        ALLOCATE(Materialdata%MaterialSets(I)%Nodes(Materialdata%MaterialSets(I)%NodesSize))
                        DO K =1,Materialdata%MaterialSets(I)%NodesSize
                            READ(UNIT,*) Materialdata%MaterialSets(I)%Nodes(K)
                        END DO
                    CASE('Isotropic')
                        READ(UNIT,*) Materialdata%MaterialSets(I)%Isotropic
                    CASE('Iso')
                        DO
                            CALL GetLine(UNIT,LINE,IOSTAT)
                            IF(IOSTAT == IOSTAT_END) RETURN
                            CALL FindKeyword(LINE,KEYWORD)
                            SELECT CASE(KEYWORD)
                            CASE('E')
                                READ(UNIT,*) Materialdata%MaterialSets(I)%Iso%E
                            CASE('nu')
                                READ(UNIT,*) Materialdata%MaterialSets(I)%Iso%nu
                            CASE DEFAULT
                                BACKSPACE(UNIT)
                                EXIT
                            END SELECT
                        END DO
                    CASE('Aniso')
                        DO
                            CALL GetLine(UNIT,LINE,IOSTAT)
                            IF(IOSTAT == IOSTAT_END) RETURN
                            CALL FindKeyword(LINE,KEYWORD)
                            SELECT CASE(KEYWORD)
                            CASE('D')
                                DO K=1,3
                                    READ(UNIT,*) Materialdata%MaterialSets(I)%Aniso%D(K,:)
                                END DO
                            CASE DEFAULT
                                BACKSPACE(UNIT)
                                EXIT
                            END SELECT
                        END DO
                    CASE('Rho')
                        READ(UNIT,*) Materialdata%MaterialSets(I)%Rho
                    CASE ('MaterialSetEnd')
                        EXIT
                    END SELECT
                END DO
            END DO
        CASE DEFAULT
            EXIT     
        END SELECT
        
        END DO

        CLOSE(UNIT)
  
    
    
    RETURN
    END SUBROUTINE PopulateMaterialData
!------------------------------------------------------------------------------    
    SUBROUTINE PopulateLoadData(datapath,Loaddata)
    CHARACTER(:),ALLOCATABLE, INTENT(IN) :: datapath
    TYPE(Loads), INTENT(INOUT) :: Loaddata
    CHARACTER(:),ALLOCATABLE :: FullPathName
    CHARACTER(:),ALLOCATABLE :: LINE, KEYWORD
    INTEGER(IntKi) :: IOSTAT, I
    INTEGER(IntKi) :: UNIT = 18

    
    !-------------------------- Reading the BC data
    FullPathName = trim(datapath) // trim(Loaddata%LoadsFileName)
    OPEN(UNIT=UNIT,FILE=FullPathName,FORM='FORMATTED',ACCESS='SEQUENTIAL',ACTION='READ',STATUS='OLD')
    DO
    CALL GetLine(UNIT,LINE,IOSTAT)
    IF (IOSTAT == IOSTAT_END) EXIT
    CALL FindKeyword(Line,KEYWORD)
    SELECT CASE(KEYWORD)
    CASE('SurfLoads')
        READ(UNIT,*) Loaddata%SurfLoads%SurfSize
        IF(loaddata%SurfLoads%SurfSize /= 0 )THEN
            ALLOCATE(Loaddata%SurfLoads%Node(Loaddata%SurfLoads%SurfSize))
            ALLOCATE(Loaddata%SurfLoads%DOF(Loaddata%SurfLoads%SurfSize))
            ALLOCATE(Loaddata%SurfLoads%Edge(Loaddata%SurfLoads%SurfSize))
            ALLOCATE(Loaddata%SurfLoads%Element(Loaddata%SurfLoads%SurfSize))
            ALLOCATE(Loaddata%SurfLoads%Value(Loaddata%SurfLoads%SurfSize))
            DO I = 1, Loaddata%SurfLoads%SurfSize
                READ(UNIT,*)Loaddata%SurfLoads%Node(I),Loaddata%SurfLoads%DOF(I),Loaddata%SurfLoads%Edge(I),Loaddata%SurfLoads%Element(I),Loaddata%SurfLoads%Value(I)
            END DO
        END IF
    CASE('BodyLoads')
        READ(UNIT,*) Loaddata%BodyLoads%BodySize
        IF(loaddata%BodyLoads%BodySize /= 0 )THEN
            ALLOCATE(Loaddata%BodyLoads%Node(Loaddata%BodyLoads%BodySize))
            ALLOCATE(Loaddata%BodyLoads%DOF(Loaddata%BodyLoads%BodySize))
            ALLOCATE(Loaddata%BodyLoads%Element(Loaddata%BodyLoads%BodySize))
            ALLOCATE(Loaddata%BodyLoads%Value(Loaddata%BodyLoads%BodySize))
            DO I = 1, Loaddata%BodyLoads%BodySize
                READ(UNIT,*)Loaddata%BodyLoads%Node(I),Loaddata%BodyLoads%DOF(I),Loaddata%BodyLoads%Value(I)
            END DO
        END IF
    END SELECT
    END DO
    CLOSE(UNIT)

    RETURN 
    END SUBROUTINE PopulateLoadData
!------------------------------------------------------------------------------
    SUBROUTINE PopulateBCData(datapath, BCDOFdata)
    CHARACTER(:),ALLOCATABLE, INTENT(IN) :: datapath
    TYPE(BCDOF), INTENT(INOUT) :: BCDOFdata
    CHARACTER(:),ALLOCATABLE :: FullPathName
    CHARACTER(:),ALLOCATABLE :: LINE, KEYWORD
    INTEGER(IntKi) :: IOSTAT, I
    INTEGER(IntKi) :: UNIT = 18

    
    !-------------------------- Reading the BC data
    FullPathName = trim(datapath) // trim(BCDOFdata%BCFileName)
    OPEN(UNIT=UNIT,FILE=FullPathName,FORM='FORMATTED',ACCESS='SEQUENTIAL',ACTION='READ',STATUS='OLD')
    READ(UNIT,*) BCDOFdata%BCSize
    ALLOCATE(BCDOFdata%Node(BCDOFdata%BCSize))
    ALLOCATE(BCDOFdata%DOF(BCDOFdata%BCSize))
    ALLOCATE(BCDOFdata%Value(BCDOFdata%BCSize))
    DO I =1,BCDOFdata%BCSize
    READ(UNIT,*) BCDOFdata%Node(I),BCDOFdata%DOF(I),BCDOFdata%Value(I)
    END DO
    ClOSE(UNIT)
    
    
    
    RETURN
    END SUBROUTINE PopulateBCData
    
    
    
    !------------------------------------------------------------------------------
    ! Destroying all Arrays in derived types
    !
    SUBROUTINE DestroyMeshData(Meshdata)
    TYPE(Mesh), INTENT(INOUT) :: Meshdata
    
    DEALLOCATE(Meshdata%GlobalCoordinates)
    DEALLOCATE(Meshdata%ConnectivityMatrix)
    DEALLOCATE(Meshdata%EdgeNodes)
    
    RETURN
    END SUBROUTINE DestroyMeshData
        !------------------------------------------------------------
    SUBROUTINE DestroySimData(Simdata)
    TYPE(SimulationData), INTENT(INOUT) :: Simdata
    DEALLOCATE(Simdata%OutputReq%STRAIN%NodeSet)
    DEALLOCATE(Simdata%OutputReq%STRAIN%Element)
    DEALLOCATE(Simdata%OutputReq%STRAIN%Edge)
    DEALLOCATE(Simdata%OutputReq%STRESS%NodeSet)
    DEALLOCATE(Simdata%OutputReq%STRESS%Element)
    DEALLOCATE(Simdata%OutputReq%STRESS%Edge)
 
    
    RETURN
    END SUBROUTINE DestroySimData
    !------------------------------------------------------------
    SUBROUTINE DestroyBCData(BCdata)
    TYPE(BCDOF), INTENT(INOUT) :: BCdata
    
    DEALLOCATE(BCdata%Node)
    DEALLOCATE(BCdata%DOF)
    DEALLOCATE(BCdata%Value)
    
    RETURN
    END SUBROUTINE DestroyBCData
 !------------------------------------------------------------
    SUBROUTINE DestroyLoadData(Loaddata)
    TYPE(Loads), INTENT(INOUT) :: Loaddata
    
    IF(ALLOCATED(Loaddata%SurfLoads%Node)) DEALLOCATE(Loaddata%SurfLoads%Node)
    IF(ALLOCATED(Loaddata%SurfLoads%DOF)) DEALLOCATE(Loaddata%SurfLoads%DOF)
    IF (ALLOCATED(Loaddata%SurfLoads%Value)) DEALLOCATE(Loaddata%SurfLoads%Value)
    IF(ALLOCATED(Loaddata%BodyLoads%Node)) DEALLOCATE(Loaddata%BodyLoads%Node)
    IF(ALLOCATED(Loaddata%BodyLoads%DOF)) DEALLOCATE(Loaddata%BodyLoads%DOF)
    IF(ALLOCATED(Loaddata%BodyLoads%Value)) DEALLOCATE(Loaddata%BodyLoads%Value)    
    RETURN
    END SUBROUTINE DestroyLoadData
 !------------------------------------------------------------
    SUBROUTINE DestroyMaterialData(Materialdata)
    TYPE(Material), INTENT(INOUT) :: Materialdata
    DEALLOCATE(Materialdata%MaterialSets)
    RETURN
    END SUBROUTINE DestroyMaterialData
END MODULE FEAInput
    
 