'''
=======================================================================================
Mercury Simulation Module (MSM): Main Program
=======================================================================================

Developed by:
* Dr. Todd E. Steissberg (ERDC-EL)
* Dr. Billy E. Johnson (ERDC-EL, LimnoTech)
* Dr. Zhonglong Zhang (Portland State University)
* Mr. Mark Jensen (HEC)

This module computes the water quality of a single computational cell. The algorithms 
and structure of this program were adapted from the Fortran 95 version of this module, 
developed by:
* Dr. Billy E. Johnson (ERDC-EL)
* Dr. Zhonglong Zhang (Portland State University)
* Mr. Mark Jensen (HEC)

Version 1.0

Initial Version: June 15, 2021
Last Revision Date: June 15, 2021
'''

'''
module modMercury
  use modGlobal
  use modDLL,    only: TempCorrectionStruct, Arrhenius_TempCorrection, AddIndex
  implicit none
  ! 
  ! temperature correction parameters (21)
  type(TempCorrectionStruct), allocatable, dimension(:)     :: kd23       ! Methylation rate from dissolved HgII to MeHg in water column    (1/day)
  real(R8)                                                  :: kd23_tc
  type(TempCorrectionStruct), allocatable, dimension(:)     :: kdoc23     ! Methylation rate from DOC sorbed HgII to MeHg in water column   (1/day)
  real(R8)                                                  :: kdoc23_tc
  type(TempCorrectionStruct), allocatable, dimension(:)     :: kso42      ! Sediment sulfate reduction rate                                 (1/day)
  real(R8)                                                  :: kso42_tc
  type(TempCorrectionStruct), allocatable, dimension(:)     :: kd32_2     ! Sediment demethylation rate                                     (1/day)
  real(R8)                                                  :: kd32_2_tc
  !
  type(TempCorrectionStruct), allocatable, dimension(:,:)   :: vv         ! Volatilization velocity of Hg0 and MeHg in water column         (m/day)
  real(R8)                                                  :: vv_tc(3)
  !
  real(R8),           allocatable, dimension(:)             :: k12        ! Oxidation rate from Hg0 to HgII in water column                         (1/day)
  real(R8)                                                  :: k12_tc
  real(R8),           allocatable, dimension(:)             :: Ea12       ! Arrhenius activation energy for Hg0 oxidation rate                      (kJ/mol)
  real(R8),           allocatable, dimension(:)             :: Tr12       ! Reference temperature for which Hg0 oxidation rate is reported          (C)
  !
  real(R8),           allocatable, dimension(:)             :: kadpom     ! Adsorption rate of POM sorbed HgII in water column                      (L/�g/day) 
  real(R8)                                                  :: kadpom_tc
  real(R8),           allocatable, dimension(:)             :: kadpom2    ! Adsorption rate of POM sorbed HgII in bed sediment                      (L/�g/day) 
  real(R8)                                                  :: kadpom2_tc
  !
  real(R8),           allocatable, dimension(:)             :: kdapom     ! Desorption rate of POM sorbed HgII in water column                      (1/day) 
  real(R8)                                                  :: kdapom_tc
  real(R8),           allocatable, dimension(:)             :: kdapom2    ! Desorption rate of POM sorbed HgII in bed sediment                      (1/day) 
  real(R8)                                                  :: kdapom2_tc
  !
  real(R8),           allocatable, dimension(:)             :: kadap      ! Adsorption rate of algae sorbed HgII in water column                    (L/�g/day) 
  real(R8)                                                  :: kadap_tc
  real(R8),           allocatable, dimension(:)             :: kdaap      ! Desorption rate of algae sorbed HgII in water column                    (1/day) 
  real(R8)                                                  :: kdaap_tc
  !
  real(R8),           allocatable, dimension(:,:)           :: kadp       ! Adsorption rate of solid j sorbed HgII in water column                  (L/�g/day) 
  real(R8),           allocatable, dimension(:)             :: kadp_tc
  real(R8),           allocatable, dimension(:,:)           :: kadp2      ! Adsorption rate of solid j sorbed HgII in bed sediment                  (L/�g/day) 
  real(R8),           allocatable, dimension(:)             :: kadp2_tc
  !
  real(R8),           allocatable, dimension(:,:)           :: kdap       ! Desorption rate of solid j sorbed HgII in water column                  (1/day) 
  real(R8),           allocatable, dimension(:)             :: kdap_tc
  real(R8),           allocatable, dimension(:,:)           :: kdap2      ! Desorption rate of solid j sorbed HgII in bed sediment                  (1/day) 
  real(R8),           allocatable, dimension(:)             :: kdap2_tc
  !
  real(R8),           allocatable, dimension(:)             :: Eaad       ! Arrhenius activation energy for adsorption rate                         (kJ/mol)
  real(R8),           allocatable, dimension(:)             :: Eade       ! Arrhenius activation energy for desorption rate                         (kJ/mol)
  real(R8),           allocatable, dimension(:)             :: Trade      ! Reference temperature for which adsorption and desorption rate is reported (C)
  !
  ! global parameters (10)
  real(R8),           allocatable, dimension(:)             :: vsap       ! Settling velocity of algae in water column                            (m/d)
  real(R8),           allocatable, dimension(:)             :: vsom       ! Settling velocity of POM in water column                              (m/d)
  real(R8),           allocatable, dimension(:)             :: h2         ! Active Sediment layer thickness                                       (m)
  real(R8),           allocatable, dimension(:)             :: Por        ! Porosity or volume water per volume bed sediment                      (unitless)
  ! 
  real(R8),           allocatable, dimension(:)             :: z2         ! The average bioturbed depth in bed sediments                          (cm)
  real(R8),           allocatable, dimension(:)             :: Db         ! Biodiffusion coefficient representing particle diffusivity in the bed (cm2/d)
  real(R8),           allocatable, dimension(:)             :: Beta       ! Water-side benthic boundary layer mass transfer coefficient           (cm/d)
  real(R8),           allocatable, dimension(:)             :: ps         ! Density of total bed sediments
  !
  real(R8),           allocatable, dimension(:)             :: alpha      ! Coefficient to adjust light extinction coefficient                    (unitless)                                             (1/d) 
  real(R8),           allocatable, dimension(:)             :: res        ! Maximum relative error of Newton-Raphson or Bisection method          (unitless) 
  ! Hg (5)
  real(R8),           allocatable, dimension(:,:)           :: vm         ! Sediment-water mass transfer velocity for HgII and MeHg               (m/d)
  real(R8),           allocatable, dimension(:,:)           :: Dm         ! Molecular diffusivity for HgII and MeHg                               (m2/d)
  real(R8),           allocatable, dimension(:,:)           :: MW         ! Molecular weight of Hg0, HgII and MeHg                                (g/mol)     
  real(R8),           allocatable, dimension(:,:)           :: Hgds       ! Solubility of Hg0, HgII and MeHg                                      (ng/l)
  real(R8),           allocatable, dimension(:,:)           :: KH         ! Henry's constant of Hg0 and MeHg                                      (Pa m3/mol)  
  ! linear partitioning (7)
  real(R8),           allocatable, dimension(:,:)           :: Kdoc       ! Dissolved organic carbon partition coefficient of HgII and MeHg in water column    (L/kg)
  real(R8),           allocatable, dimension(:,:)           :: Kdoc2      ! Dissolved organic carbon partition coefficient of HgII and MeHg in bed sediment    (L/kg)
  real(R8),           allocatable, dimension(:,:)           :: Kpom       ! Organic matter partition coefficient of HgII and MeHg in water column              (L/kg)
  real(R8),           allocatable, dimension(:,:)           :: Kpom2      ! Organic matter partition coefficient of HgII and MeHg in bed sediment              (L/kg)
  real(R8),           allocatable, dimension(:,:)           :: Kap        ! Algae partition coefficient of HgII and MeHg in water column                       (L/kg)
  real(R8),           allocatable, dimension(:,:,:)         :: Kp         ! Solid partition coefficient of HgII and MeHg in water column                       (L/kg)
  real(R8),           allocatable, dimension(:,:,:)         :: Kp2        ! Solid partition coefficient of HgII and MeHg in bed sediment                       (L/kg)
  ! Freundlich and Langmuir equilibrium partitioning (20) 
  real(R8),           allocatable, dimension(:,:)           :: Klpom      ! Langmuir adsorption constant of HgII (MeHg) on organic matter in water column    (L/�g)
  real(R8),           allocatable, dimension(:,:)           :: Klpom2     ! Langmuir adsorption constant of HgII (MeHg) on organic matter in bed sediment    (L/�g)
  real(R8),           allocatable, dimension(:,:)           :: qcpom      ! Adsorption capacity of HgII (MeHg) on POM in water column                        (�g/g)
  real(R8),           allocatable, dimension(:,:)           :: qcpom2     ! Adsorption capacity of HgII (MeHg) on POM in bed sediment                        (�g/g)
  real(R8),           allocatable, dimension(:,:)           :: Kfpom      ! Freundlich adsorption constant of HgII (MeHg) on organic matter in water column  (L/kg)
  real(R8),           allocatable, dimension(:,:)           :: Kfpom2     ! Freundlich adsorption constant of HgII (MeHg) on organic matter in bed sediment  (L/kg)
  real(R8),           allocatable, dimension(:,:)           :: bpom       ! Freundlich exponent of HgII (MeHg) on organic matter in water column             (unitless)
  real(R8),           allocatable, dimension(:,:)           :: bpom2      ! Freundlich exponent of HgII (MeHg) on organic matter in bed sediment             (unitless)
  real(R8),           allocatable, dimension(:,:)           :: Klap       ! Langmuir adsorption constant of HgII (MeHg) on algae in water column             (L/�g)
  real(R8),           allocatable, dimension(:,:)           :: qcap       ! Adsoption capacity of HgII (MeHg) on algae                                       (�g/g)
  real(R8),           allocatable, dimension(:,:)           :: Kfap       ! Freundlich adsorption constant of HgII (MeHg) on algae                           (L/kg)
  real(R8),           allocatable, dimension(:,:)           :: bap        ! Freundlich exponent of HgII (MeHg) on algae                                      (unitless)
  real(R8),           allocatable, dimension(:,:,:)         :: Klp        ! Langmuir adsorption constant of HgII (MeHg) on solid in water column             (L/�g)
  real(R8),           allocatable, dimension(:,:,:)         :: Klp2       ! Langmuir adsorption constant of HgII (MeHg) on solid in bed sediment             (L/�g)
  real(R8),           allocatable, dimension(:,:,:)         :: qcp        ! Adsorption capacity of HgII (MeHg) on solid in water column                      (�g/g)
  real(R8),           allocatable, dimension(:,:,:)         :: qcp2       ! Adsorption capacity of HgII (MeHg) on solid in bed sediment                      (�g/g)
  real(R8),           allocatable, dimension(:,:,:)         :: Kfp        ! Freundlich adsorption constant of HgII (MeHg) on solid in water column           (L/kg)
  real(R8),           allocatable, dimension(:,:,:)         :: Kfp2       ! Freundlich adsorption constant of HgII (MeHg) on solid in bed sediment           (L/kg)
  real(R8),           allocatable, dimension(:,:,:)         :: bp         ! Freundlich exponent of HgII (MeHg) on solid in water column                      (unitless)
  real(R8),           allocatable, dimension(:,:,:)         :: bp2        ! Freundlich exponent of HgII (MeHg) on solid in bed sediment                      (unitless)
  ! transform parameters (14)
  real(R8),           allocatable, dimension(:)             :: kd21       ! Reduction rate from dissolved HgII to Hg0 in water column                        (1/day)
  real(R8),           allocatable, dimension(:)             :: kdoc21     ! Reduction rate from DOC partitioned HgII to Hg0 in water column                  (1/day)
  real(R8),           allocatable, dimension(:,:)           :: I0pht      ! Light intensity at which kd21 and kdoc21 (kd31 and kdoc31) is measured           (w/m2)
  real(R8),           allocatable, dimension(:)             :: kd31       ! Photoreduction rate from dissolved MeHg to Hg0 in water column                   (1/day)
  real(R8),           allocatable, dimension(:)             :: kdoc31     ! Photoreduction rate from DOC partitioned MeHg to Hg0 in water column             (1/day)                               
  real(R8),           allocatable, dimension(:)             :: kd32       ! Demethylation rate from dissolved MeHg to HgII in water                          (1/day)
  real(R8),           allocatable, dimension(:)             :: kdoc32     ! Demethylation rate from DOC sorbed MeHg to HgII in water                         (1/day)
  real(R8),           allocatable, dimension(:)             :: rmso4      ! Ratio of sediment methylation rate and sulfate reduction rate                    (L/mg)                                                                                                
  real(R8),           allocatable, dimension(:)             :: Kso4       ! Half-saturation constant for the effect of sulfate on methylation                (mg-O2/L)                                         
  real(R8),           allocatable, dimension(:)             :: Y12        ! Oxidation yield coefficient from Hg0 to HgII in water column                     (g/g)
  real(R8),           allocatable, dimension(:)             :: Y21        ! Reduction yield coefficient from HgII to Hg0 in water column                     (g/g)
  real(R8),           allocatable, dimension(:)             :: Y23        ! Methylation yield coefficient from HgII to MeHg in water column                  (g/g)
  real(R8),           allocatable, dimension(:)             :: Y31        ! Photodegradation yield coefficient from MeHg to Hg0 in water column              (g/g)
  real(R8),           allocatable, dimension(:)             :: Y32        ! Demethylation yield coefficient from MeHg to HgII in water column                (g/g)
  ! integer parameter (3)
  integer,            allocatable, dimension(:,:)           :: vv_option              ! 1 user defined; 2 compute       
  integer,            allocatable, dimension(:,:)           :: vm_option              ! 1 user defined; 2-5 compute
  integer,            allocatable, dimension(:,:)           :: Hgd_solution_option    ! 1 Newton;       2 Bisection
  !
  ! pathways (25)
  ! Hg0
  real(R8)           :: Hg0_Volatilization            ! Hg0 volatilization in water column    (ng/L/d) 
  real(R8)           :: Hg0_Oxidation                 ! Hg0 oxidation in water column         (ng/L/d)
  ! HgII
  real(R8)           :: HgII_Air_Deposition           ! HgII air deposition in water column   (ng/L/d)  
  real(R8)           :: HgII_Photoreduction           ! HgII photoreduction in water column   (ng/L/d) 
  real(R8)           :: HgII_Methylation              ! HgII methylation in water column      (ng/L/d)
  real(R8)           :: HgII2_Methylation             ! HgII methylation in bed sediment      (ng/L/d)
  real(R8)           :: HgII_Settling                 ! HgII settling from water column       (ng/L/d) 
  real(R8)           :: HgII_Deposition               ! HgII settling to bed sediment         (ng/L/d) 
  real(R8)           :: HgII_Resuspension             ! HgII resuspension to water column     (ng/L/d) 
  real(R8)           :: HgII_Erosion                  ! HgII resuspension from bed sediment   (ng/L/d)
  real(R8)           :: HgII_Burial                   ! HgII burial in bed sediment           (ng/L/d)
  real(R8)           :: HgII_Transfer                 ! HgII transfer to water column         (ng/L/d) 
  real(R8)           :: HgII2_Transfer                ! HgII transfer to bed sediment         (ng/L/d) 
  ! MeHg
  real(R8)           :: MeHg_Air_Deposition           ! MeHg air deposition in water column   (ng/L/d)
  real(R8)           :: MeHg_Volatilization           ! MeHg volatilization in water column   (ng/L/d) 
  real(R8)           :: MeHg_Photoreduction           ! MeHg photoreduction in water column   (ng/L/d) 
  real(R8)           :: MeHg_Demethylation            ! MeHg demethylation in water column    (ng/L/d) 
  real(R8)           :: MeHg2_Demethylation           ! MeHg demethylation in bed sediment    (ng/L/d)
  real(R8)           :: MeHg_Settling                 ! MeHg settling from water column       (ng/L/d)
  real(R8)           :: MeHg_Deposition               ! MeHg settling to bed sediment         (ng/L/d) 
  real(R8)           :: MeHg_Resuspension             ! MeHg resuspension to water column     (ng/L/d)
  real(R8)           :: MeHg_Erosion                  ! MeHg resuspension from bed sediment   (ng/L/d)  
  real(R8)           :: MeHg_Burial                   ! MeHg burial in bed sediment           (ng/L/d)
  real(R8)           :: MeHg_Transfer                 ! MeHg transfer to water column         (ng/L/d) 
  real(R8)           :: MeHg2_Transfer                ! MeHg transfer to bed sediment         (ng/L/d)  
  !
  ! pathway index (25)
  integer            :: Hg0_Volatilization_index            
  integer            :: Hg0_Oxidation_index 
  !               
  integer            :: HgII_Air_deposition_index 
  integer            :: HgII_Photoreduction_index                
  integer            :: HgII_Methylation_index 
  integer            :: HgII2_Methylation_index  
  integer            :: HgII_Settling_index   
  integer            :: HgII_Deposition_index 
  integer            :: HgII_Resuspension_index  
  integer            :: HgII_Erosion_index  
  integer            :: HgII_Burial_index 
  integer            :: HgII_Transfer_index 
  integer            :: HgII2_Transfer_index  
  !          
  integer            :: MeHg_Air_Deposition_index                         
  integer            :: MeHg_Volatilization_index            
  integer            :: MeHg_Photoreduction_index      
  integer            :: MeHg_Demethylation_index 
  integer            :: MeHg2_Demethylation_index                           
  integer            :: MeHg_Settling_index                 
  integer            :: MeHg_Deposition_index
  integer            :: MeHg_Resuspension_index 
  integer            :: MeHg_Erosion_index
  integer            :: MeHg_Burial_index           
  integer            :: MeHg_Transfer_index                                    
  integer            :: MeHg2_Transfer_index                                 
  !
  ! local variables 
  real(R8) :: Cd2(2:3), Cdoc2(2:3)
  !
  real(R8) :: gas_constant = 8.314    ! universal gas constant (J/mole/K)
  real(R8) :: Iav(2:3), SO4_reduction
  real(R8) :: TwaterK, TsedK
  integer  :: imax = 100
  integer  :: i, k
  logical  :: IsWaterCell    
  !
  contains
  !===========================================================================================================================
  ! allocate and initialize all the input parameters with default values
  subroutine InitializeHg() 
    !
    ! Linear equilibrium partition for DOC
    if (use_DOCSorbed(2) .or. use_DOCSorbed(3)) then
      if (allocated(Kdoc)) deallocate(Kdoc) 
      allocate(Kdoc(2:3,nRegion))
      Kdoc = 2.0E5
      !
      if (use_BedSediment) then
        if (allocated(Kdoc2)) deallocate(Kdoc2)
        allocate(Kdoc2(2:3,nRegion))
        Kdoc2 = 2.0E5
      end if
    end if
    !
    ! Linear equilibrium partition for algae, POM and solids  
    if ((use_Equilibrium(2,1) .and. use_AlgaeSorbed(2)) .or. (use_Equilibrium(3,1) .and. use_AlgaeSorbed(3))) then    
      if (allocated(Kap)) deallocate(Kap)
      allocate(Kap(2:3,nRegion))
      Kap = 2.0E5
    end if 
    !
    if ((use_Equilibrium(2,1) .and. use_POMSorbed(2)) .or. (use_Equilibrium(3,1) .and. use_POMSorbed(3))) then
      if (allocated(Kpom)) deallocate(Kpom) 
      allocate(Kpom(2:3,nRegion))
      Kpom = 2.0E5  
    end if
    !
    if ((use_Equilibrium(2,1) .and. use_AnySolidSorbed(2)) .or. (use_Equilibrium(3,1) .and. use_AnySolidSorbed(3))) then    
      if (allocated(Kp)) deallocate(Kp)
      allocate(Kp(2:3,nGS,nRegion))
      Kp = 2.0E5
    end if
    !
    if ((use_Equilibrium(2,2) .and. use_POMSorbed(2)) .or. (use_Equilibrium(3,2) .and. use_POMSorbed(3))) then
      if (allocated(Kpom2)) deallocate(Kpom2) 
      allocate(Kpom2(2:3,nRegion))
      Kpom2 = 2.0E5
    end if
    !
    if ((use_Equilibrium(2,2) .and. use_AnySolidSorbed(2)) .or. (use_Equilibrium(3,2) .and. use_AnySolidSorbed(3))) then     
      if (allocated(Kp2)) deallocate(Kp2)
      allocate(Kp2(2:3,nGS,nRegion))
      Kp2 = 2.0E5      
    end if
    !
    ! Langmuir equilibrium partition for algae, POM and solids
    if ((use_Langmuir(2,1) .and. use_AlgaeSorbed(2)) .or. (use_Langmuir(3,1) .and. use_AlgaeSorbed(3))) then
      if (allocated(Klap)) deallocate(Klap)
      allocate(Klap(2:3, nRegion))
      Klap = 1.0E5
      !
      if (allocated(qcap)) deallocate(qcap)
      allocate(qcap(2:3, nRegion))
      qcap = 1000.0
    end if
    !
    if ((use_Langmuir(2,1) .and. use_POMSorbed(2)) .or. (use_Langmuir(3,1) .and. use_POMSorbed(3))) then    
      if (allocated(Klpom)) deallocate(Klpom)
      allocate(Klpom(2:3, nRegion))
      Klpom = 1.0E5
      !
      if (allocated(qcpom)) deallocate(qcpom)
      allocate(qcpom(2:3, nRegion))
      qcpom = 1000.0
    end if
    !
    if ((use_Langmuir(2,1) .and. use_AnySolidSorbed(2)) .or. (use_Langmuir(3,1) .and. use_AnySolidSorbed(3))) then    
      if (allocated(Klp)) deallocate(Klp)
      allocate(Klp(2:3, nGS,nRegion))
      Klp = 1.0E5
      !
      if (allocated(qcp)) deallocate(qcp)
      allocate(qcp(2:3, nGS,nRegion))
      qcp = 1000.0      
    end if
    !
    if ((use_Langmuir(2,2) .and. use_POMSorbed(2)) .or. (use_Langmuir(3,2) .and. use_POMSorbed(3))) then        
      if (allocated(Klpom2)) deallocate(Klpom2)
      allocate(Klpom2(2:3, nRegion))
      Klpom2 = 1.0E5
      !
      if (allocated(qcpom2)) deallocate(qcpom2)
      allocate(qcpom2(2:3, nRegion))
      qcpom2 = 1000.0
    end if
    !
    if ((use_Langmuir(2,2) .and. use_AnySolidSorbed(2)) .or. (use_Langmuir(3,2) .and. use_AnySolidSorbed(3))) then  
      if (allocated(Klp2)) deallocate(Klp2)
      allocate(Klp2(2:3, nGS,nRegion))
      Klp2 = 1.0E5
      !
      if (allocated(qcp2)) deallocate(qcp2)
      allocate(qcp2(2:3, nGS,nRegion))
      qcp2 = 1000.0     
    end if
    !
    ! Freundlich equilibrium partition for algae, POM and solids
    if ((use_Freundlich(2,1) .and. use_AlgaeSorbed(2)) .or. (use_Freundlich(3,1) .and. use_AlgaeSorbed(3))) then
      if (allocated(Kfap)) deallocate(Kfap)
      allocate(Kfap(2:3, nRegion))
      Kfap = 1.0E2
      !
      if (allocated(bap)) deallocate(bap)
      allocate(bap(2:3, nRegion))
      bap = 1.0
    end if
    !
    if ((use_Freundlich(2,1) .and. use_POMSorbed(2)) .or. (use_Freundlich(3,1) .and. use_POMSorbed(3))) then
      if (allocated(Kfpom)) deallocate(Kfpom)
      allocate(Kfpom(2:3, nRegion))
      Kfpom = 1.0E2
      !
      if (allocated(bpom)) deallocate(bpom)
      allocate(bpom(2:3, nRegion))
      bpom = 1.0
    end if
    !
    if ((use_Freundlich(2,1) .and. use_AnySolidSorbed(2)) .or. (use_Freundlich(3,1) .and. use_AnySolidSorbed(3))) then
      if (allocated(Kfp)) deallocate(Kfp)
      allocate(Kfp(2:3, nGS,nRegion))
      Kfp = 1.0E2
      !
      if (allocated(bp)) deallocate(bp)
      allocate(bp(2:3, nGS,nRegion))
      bp = 1.0
    end if
    !
    if ((use_Freundlich(2,2) .and. use_POMSorbed(2)) .or. (use_Freundlich(3,2) .and. use_POMSorbed(3))) then  
      if (allocated(Kfpom2)) deallocate(Kfpom2)
      allocate(Kfpom2(2:3, nRegion))
      Kfpom2 = 1.0E2
      !
      if (allocated(bpom2)) deallocate(bpom2)
      allocate(bpom2(2:3, nRegion))
      bpom2 = 1.0
    end if
    !
    if ((use_Freundlich(2,2) .and. use_AnySolidSorbed(2)) .or. (use_Freundlich(3,2) .and. use_AnySolidSorbed(3))) then 
      if (allocated(Kfp2)) deallocate(Kfp2)
      allocate(Kfp2(2:3, nGS,nRegion))
      Kfp2 = 1.0E2
      !
      if (allocated(bp2)) deallocate(bp2)
      allocate(bp2(2:3, nGS,nRegion))
      bp2 = 1.0
    end if
    ! 
    ! Non-equilibrium partition for algae, POM and solids
    if (use_NonEquilibrium(1)) then
      if (use_AlgaeSorbed(2)) then
        if (allocated(kadap)) deallocate(kadap)
        allocate(kadap(nRegion))
        kadap = 0.0001
        !
        if (allocated(kdaap)) deallocate(kdaap)
        allocate(kdaap(nRegion))
        kdaap = 0.1
        !
        if (allocated(qcap)) deallocate(qcap)
        allocate(qcap(2:3, nRegion))
        qcap = 1000.0
      end if
      !
      if (use_POMSorbed(2)) then
        if (allocated(kadpom)) deallocate(kadpom)
        allocate(kadpom(nRegion))
        kadpom = 0.0001
        !
        if (allocated(kdapom)) deallocate(kdapom)
        allocate(kdapom(nRegion))
        kdapom = 0.1
        !
        if (allocated(qcpom)) deallocate(qcpom)
        allocate(qcpom(2:3, nRegion))
        qcpom = 1000.0
      end if
      !
      if (use_AnySolidSorbed(2)) then
        if (allocated(kadp)) deallocate(kadp, kadp_tc)
        allocate(kadp(nGS,nRegion), kadp_tc(nGS))
        kadp = 0.0001
        !
        if (allocated(kdap)) deallocate(kdap, kdap_tc)
        allocate(kdap(nGS,nRegion), kdap_tc(nGS))
        kdap = 0.1
        !
        if (allocated(qcp)) deallocate(qcp)
        allocate(qcp(2:3, nGS,nRegion))
        qcp = 1000.0
      end if
    end if
    !
    if (use_NonEquilibrium(2)) then
      if (use_POMSorbed(2)) then
        if (allocated(kadpom2)) deallocate(kadpom2)
        allocate(kadpom2(nRegion))
        kadpom2 = 0.0001
        !
        if (allocated(kdapom2)) deallocate(kdapom2)
        allocate(kdapom2(nRegion))
        kdapom2 = 0.1
        !
        if (allocated(qcpom2)) deallocate(qcpom2)
        allocate(qcpom2(2:3, nRegion))
        qcpom2 = 1000.0
      end if
      !
      if (use_AnySolidSorbed(2)) then
        if (allocated(kadp2)) deallocate(kadp2, kadp2_tc)
        allocate(kadp2(nGS,nRegion), kadp2_tc(nGS))
        kadp2 = 0.0001
        !
        if (allocated(kdap2)) deallocate(kdap2, kdap2_tc)
        allocate(kdap2(nGS,nRegion), kdap2_tc(nGS))
        kdap2 = 0.1
        !
        if (allocated(qcp2)) deallocate(qcp2)
        allocate(qcp2(2:3, nGS,nRegion))
        qcp2 = 1000.0
      end if
    end if
    !
    if (use_NonEquilibrium(1) .or. use_NonEquilibrium(2)) then
      if (use_AlgaeSorbed(2) .or. use_POMSorbed(2) .or. use_AnySolidSorbed(2)) then
        if (allocated(Eaad)) deallocate(Eaad)
        allocate(Eaad(nRegion))
        Eaad = 50.0
        !
        if (allocated(Eade)) deallocate(Eade)
        allocate(Eade(nRegion))
        Eade = 50.0
        !
        if (allocated(TRade)) deallocate(TRade)
        allocate(TRade(nRegion))
        Trade = 20.0
      end if
    end if
    !
    ! Reaction processes 
    if (allocated(vv)) deallocate(vv)
    allocate(vv(3, nRegion))
    vv(1, nRegion)%rc20 = 0.144;  vv(1, nRegion)%theta = 1.024
    vv(3, nRegion)%rc20 = 1.9E-5; vv(3, nRegion)%theta = 1.024
    !
    if (allocated(KH)) deallocate(KH)
    allocate(KH(3, nRegion))
    KH(1, nRegion) = 0.09   
    KH(3, nRegion) = 4.5E-6 
    ! 
    if (allocated(k12)) deallocate(k12)
    allocate(k12(nRegion))
    k12 = 1.0E-3
    !
    if (allocated(Y12)) deallocate(Y12)
    allocate(Y12(nRegion))
    Y12 = 1.0
    !
    if (allocated(Ea12)) deallocate(Ea12)
    allocate(Ea12(nRegion))
    Ea12 = 50.0
    !
    if (allocated(Tr12)) deallocate(Tr12)
    allocate(Tr12(nRegion))
    Tr12 = 20.0
    !
    if (allocated(kd21)) deallocate(kd21)
    allocate(kd21(nRegion))
    kd21 = 0.05
    !
    if (allocated(Y21)) deallocate(Y21)
    allocate(Y21(nRegion))
    Y21 = 1.0
    !
    if (allocated(I0pht)) deallocate(I0pht)
    allocate(I0pht(2:3,nRegion))
    I0pht = 100.0
    !
    if (use_DOCSorbed(2)) then
      if (allocated(kdoc21)) deallocate(kdoc21)
      allocate(kdoc21(nRegion))
      kdoc21 = 0.0
    end if 
    !
    if (allocated(kd23)) deallocate(kd23)
    allocate(kd23(nRegion))
    kd23%rc20 = 1.0E-3; kd23%theta = 1.013
    !
    if (allocated(Y23)) deallocate(Y23)
    allocate(Y23(nRegion)) 
    Y23 = 1.07
    !
    if (use_DOCSorbed(2)) then
      if (allocated(kdoc23)) deallocate(kdoc23)
      allocate(kdoc23(nRegion))
      kdoc23%rc20 = 0.0; kdoc23%theta = 1.013
    end if
    !
    if (use_BedSediment) then  
      if (allocated(kso42)) deallocate(kso42)
      allocate(kso42(nRegion))
      kso42(nRegion)%rc20 = 0.0; kso42(nRegion)%theta = 1.024
      !
      if (allocated(rmso4)) deallocate(rmso4)
      allocate(rmso4(nRegion))   
      rmso4 = 1.0
      !
      if (allocated(kso4)) deallocate(Kso4)
      allocate(Kso4(nRegion))
      Kso4 = 0.5
    end if
    !
    if (allocated(kd31)) deallocate(kd31)
    allocate(kd31(nRegion))
    kd31 = 0.01
    !
    if (allocated(Y31)) deallocate(Y31)
    allocate(Y31(nRegion))
    Y31 = 0.93
    !
    if (use_DOCSorbed(3)) then
      if (allocated(kdoc31)) deallocate(kdoc31)
      allocate(kdoc31(nRegion))
      kdoc31 = 0.0 
    end if
    !
    if (allocated(kd32)) deallocate(kd32)
    allocate(kd32(nRegion))
    kd32 = 0.05
    !
    if (allocated(Y32)) deallocate(Y32)
    allocate(Y32(nRegion))
    Y32 = 0.93
    !
    if (use_DOCSorbed(3)) then
      if (allocated(kdoc32)) deallocate(kdoc32)
      allocate(kdoc32(nRegion))
      kdoc32 = 0.0
    end if
    !
    if (use_BedSediment) then
      if (allocated(kd32_2)) deallocate(kd32_2)
      allocate(kd32_2(nRegion))
      kd32_2%rc20 = 0.2; kd32_2%theta = 1.013  
    end if 
    !
    ! Global parameters
    if (use_AlgaeSorbed(2) .or. use_AlgaeSorbed(3)) then
      if (allocated(vsap)) deallocate(vsap)
      allocate(vsap(nRegion))
      vsap = 0.15
    end if
    !
    if (use_POMSorbed(2) .or. use_POMSorbed(3)) then
      if (allocated(vsom)) deallocate(vsom)
      allocate(vsom(nRegion))
      vsom = 0.1 
    end if
    !
    if (use_BedSediment) then
      if (allocated(h2)) deallocate(h2)
      allocate(h2(nRegion))
      h2 = 0.1
      !   
      if (allocated(Por)) deallocate(Por)
      allocate(Por(nRegion))
      Por = 0.9
      ! 
      if (allocated(z2)) deallocate(z2)
      allocate(z2(nRegion))
      z2 = 5.0
      !
      if (allocated(Db)) deallocate(Db)
      allocate(Db(nRegion))
      Db = 0.03
      !
      if (allocated(beta)) deallocate(beta)
      allocate(beta(nRegion))
      beta = 33.3
      !
      if (allocated(ps)) deallocate(ps)
      allocate(ps(nRegion))
      ps = 2.7
    end if
    ! 
    if (allocated(alpha)) deallocate(alpha)
    allocate(alpha(nRegion))
    alpha = 1.3 
    ! 
    if (((use_Langmuir(2,1) .or. use_Freundlich(2,1) .or. use_Langmuir(2,2) .or. use_Freundlich(2,2)) .and. (use_AlgaeSorbed(2) .or. use_POMSorbed(2) .or. use_AnySolidSorbed(2)))  &     
    .or. ((use_Langmuir(3,1) .or. use_Freundlich(3,1) .or. use_Langmuir(3,2) .or. use_Freundlich(3,2)) .and. (use_AlgaeSorbed(3) .or. use_POMSorbed(3) .or. use_AnySolidSorbed(3)))) then
      if (allocated(res)) deallocate(res)
      allocate(res(nRegion))
      res = 0.001
    end if
    !
    if (allocated(MW)) deallocate(MW)
    allocate(MW(3, nRegion))
    MW(1, nRegion) = 200.59
    MW(2, nRegion) = 271.52
    MW(3, nRegion) = 230.66
    !
    if (allocated(Hgds)) deallocate(Hgds)
    allocate(Hgds(3, nRegion))
    Hgds(1, nRegion) = 5.6E4
    Hgds(2, nRegion) = 2.86E10
    Hgds(3, nRegion) = 1.0E8
    !
    if (use_BedSediment) then 
      if (allocated(Dm)) deallocate(Dm)
      allocate(Dm(2:3,nRegion))
      Dm = 0.0001
      !
      if (allocated(vm)) deallocate(vm)
      allocate(vm(2:3,nRegion))
      vm = 0.01 
    end if
    !
    ! integer parameters
    if (allocated(vv_option)) deallocate(vv_option)
    allocate(vv_option(3, nRegion))
    vv_option = 1
    !
    if (use_BedSediment) then
      if (allocated(vm_option)) deallocate(vm_option)
      allocate(vm_option(2:3,nRegion))
      vm_option = 1
    end if
    !
    if (((use_Langmuir(2,1) .or. use_Freundlich(2,1) .or. use_Langmuir(2,2) .or. use_Freundlich(2,2)) .and. (use_AlgaeSorbed(2) .or. use_POMSorbed(2) .or. use_AnySolidSorbed(2)))  & 
   .or. ((use_Langmuir(3,1) .or. use_Freundlich(3,1) .or. use_Langmuir(3,2) .or. use_Freundlich(3,2)) .and. (use_AlgaeSorbed(3) .or. use_POMSorbed(3) .or. use_AnySolidSorbed(3)))) then 
      if (allocated(Hgd_solution_option)) deallocate (Hgd_solution_option)
      allocate(Hgd_solution_option(2:3, nRegion))
      Hgd_solution_option = 1
    end if
    !
  end subroutine
  !===========================================================================================================================
  ! set a real parameter
  logical function SetHgRealParameter(groupName, name, paramValue)
    character(len=*), intent(in) :: groupName 
    character(len=*), intent(in) :: name
    real(R8),         intent(in) :: paramValue(nRegion)
    integer groupIndex   
    !
    SetHgRealParameter = .true. 
    !
    if (groupName      == 'Element mercury') then
      groupIndex = 1
    else if (groupName == 'Inorganic mercury') then
      groupIndex = 2
    else if (groupName == 'Methylmercury') then
      groupIndex = 3
    end if
    !
    if (name(1:3) == 'Kp_' .or. name(1:4) == 'Kp2_') then
      do i = 1, nGS
        if (name == trim(AddIndex('Kp_',i))) then
          Kp(groupIndex,i,:)  = paramValue
          exit
        else if (name == trim(AddIndex('Kp2_',i))) then
          Kp2(groupIndex,i,:) = paramValue
          exit
        end if
      end do
    else if (name(1:4) == 'Klp_' .or. name(1:5) == 'Klp2_' .or. name(1:4) == 'Kfp_' .or.  & 
             name(1:5) == 'Kfp2_' .or. name(1:3) == 'bp_' .or. name(1:4) == 'bp2_') then
      do i = 1, nGS
        if (name == trim(AddIndex('Klp_',i))) then
          Klp(groupIndex,i,:)  = paramValue
          exit
        else if (name == trim(AddIndex('Klp2_',i))) then
          Klp2(groupIndex,i,:) = paramValue
          exit
        else if (name == trim(AddIndex('Kfp_',i))) then
          Kfp(groupIndex,i,:)  = paramValue
          exit
        else if (name == trim(AddIndex('Kfp2_',i))) then
          Kfp2(groupIndex,i,:) = paramValue
          exit
        else if (name == trim(AddIndex('bp_',i))) then
          bp(groupIndex,i,:)   = paramValue
          exit
        else if (name == trim(AddIndex('bp2_',i))) then
          bp2(groupIndex,i,:)  = paramValue
          exit
        end if
      end do
    else if (name(1:5) == 'kadp_' .or. name(1:5) == 'kdap_'  .or. name(1:4) == 'qcp_' .or.  & 
             name(1:5) == 'qcp2_' .or. name(1:6) == 'kadp2_' .or. name(1:6) == 'kdap2_') then
      do i = 1, nGS
        if (name == trim(AddIndex('kadp_',i))) then
          kadp(i,:) = paramValue
          exit
        else if (name == trim(AddIndex('kdap_',i))) then
          kdap(i,:) = paramValue
          exit
        else if (name == trim(AddIndex('qcp_',i))) then
          qcp(groupIndex,i,:) = paramValue
          exit
        else if (name == trim(AddIndex('kadp2_',i))) then
          kadp2(i,:) = paramValue
          exit
        else if (name == trim(AddIndex('kdap2_',i))) then
          kdap2(i,:) = paramValue
          exit
        else if (name == trim(AddIndex('qcp2_',i))) then
          qcp2(groupIndex,i,:) = paramValue
          exit
        end if
      end do
    else
      select case (name) 
        !
        ! Linear equilibrium partition  
        case ('Kdoc')
          Kdoc(groupIndex,:)     = paramValue
        case ('Kdoc2')
          Kdoc2(groupIndex,:)    = paramValue
        case ('Kap')
          Kap(groupIndex,:)      = paramValue
        case ('Kpom')
          Kpom(groupIndex,:)     = paramValue
        case ('Kpom2')
          Kpom2(groupIndex,:)    = paramValue
        !
        ! Langmuir equilibrium partition
        case ('Klap')
          Klap(groupIndex,:)     = paramValue
        case ('qcap')
          qcap(groupIndex,:)     = paramValue
        case ('Klpom')
          Klpom(groupIndex,:)    = paramValue
        case ('qcpom')
          qcpom(groupIndex,:)    = paramValue
        case ('Klpom2')
          Klpom2(groupIndex,:)   = paramValue
        case ('qcpom2')
          qcpom2(groupIndex,:)   = paramValue
        !
        ! Freundlich-equilibrium partition
        case ('Kfap')
          Kfap(groupIndex,:)     = paramValue
        case ('bap')
          bap(groupIndex,:)      = paramValue
        case ('Kfpom')
          Kfpom(groupIndex,:)    = paramValue
        case ('bpom')
          bpom(groupIndex,:)     = paramValue
        case ('Kfpom2')
          Kfpom2(groupIndex,:)   = paramValue
        case ('bpom2')
          bpom2(groupIndex,:)    = paramValue
        !
        ! Non-equilibrium partition
        case('kadap')
          kadap(:)               = paramValue
        case('kdaap')
          kdaap(:)               = paramValue        
        case('kadpom')
          kadpom(:)              = paramValue
        case('kdapom')
          kdapom(:)              = paramValue
        case('kadpom2')
          kadpom2(:)             = paramValue        
        case('kdapom2')
          kdapom2(:)             = paramValue 
        case('Eaad')
          Eaad(:)                = paramValue
        case('Eade')
          Eade(:)                = paramValue        
        case('Trade')
          Trade(:)               = paramValue
        !
        ! Processes 
        case ('k12')
          k12(:)                 = paramValue
        case ('Y12')
          Y12(:)                 = paramValue
        case ('Ea12')
          Ea12(:)                = paramValue
        case ('Tr12')
          Tr12(:)                = paramValue
        case ('I0pht')
          I0pht(groupIndex,:)    = paramValue
        case ('kd21')
          kd21(:)                = paramValue
        case ('kdoc21')
          kdoc21(:)              = paramValue
        case ('Y21')
          Y21(:)                 = paramValue
        case ('kd23_rc20')
          kd23(:)%rc20           = paramValue
        case ('kd23_theta')
          kd23(:)%theta          = paramValue
        case ('kdoc23_rc20')
          kdoc23(:)%rc20         = paramValue
        case ('kdoc23_theta')
          kdoc23(:)%theta        = paramValue
        case ('Y23')
          Y23(:)                 = paramValue
        case ('kso42_rc20')
          kso42(:)%rc20          = paramValue
        case ('kso42_theta')
          kso42(:)%theta         = paramValue
        case ('rmso4')
          rmso4(:)               = paramValue
        case ('Kso4')
          Kso4(:)                = paramValue
        case ('kd31')
          kd31(:)                = paramValue
        case ('kdoc31')
          kdoc31(:)              = paramValue
        case ('Y31')
          Y31(:)                 = paramValue
        case ('kd32')
          kd32(:)                = paramValue
        case ('kdoc32')
          kdoc32(:)              = paramValue
        case ('Y32')
          Y32(:)                 = paramValue
        case ('kd32_2_rc20')
          kd32_2(:)%rc20         = paramValue
        case ('kd32_2_theta')
          kd32_2(:)%theta        = paramValue  
        case ('vv_rc20')
          vv(groupIndex,:)%rc20  = paramValue
        case ('vv_theta')
          vv(groupIndex,:)%theta = paramValue
        case ('KH')
          KH(groupIndex, :)      = paramValue  
        !  
        ! Global parameters   
        case ('vsap')
          vsap                   = paramValue
        case ('vsom')
          vsom                   = paramValue
        case ('h2')
          h2                     = paramValue
        case ('Por')
          Por                    = paramValue
        !  
        case ('beta')
          beta                   = paramValue
        case ('z2')
          z2                     = paramValue
        case ('Db')
          Db                     = paramValue
        case ('ps')
          ps                     = paramValue   
        case ('res')
          res                    = paramValue 
        case ('alpha')
          alpha                  = paramValue  
        !
        case ('MW')
          MW(groupIndex, :)      = paramValue
        case ('Hgds')
          Hgds(groupIndex,:)     = paramValue
        case ('Dm')
          Dm(groupIndex,:)       = paramValue
        case ('vm')
          vm(groupIndex,:)       = paramValue
        case default
          ! did not find the parameter, return false
          SetHgRealParameter = .false.
      end select       
    end if
    !
  end function
  !===========================================================================================================================
  ! set integer parameter
  logical function SetHgIntegerParameter(groupName, name, paramValue)      
    character(len=*), intent(in) :: groupName
    character(len=*), intent(in) :: name
    integer,          intent(in) :: paramValue(nRegion)
    integer groupIndex
    !
    SetHgIntegerParameter = .true.
    if (groupName      == 'Element mercury') then
      groupIndex = 1
    else if (groupName == 'Inorganic mercury') then
      groupIndex = 2
    else if (groupName == 'Methylmercury') then
      groupIndex = 3
    end if
    select case (name)
      case ('vm Option')
        vm_option(groupIndex,:)            = paramValue
      case ('vv Option')
        vv_option(groupIndex,:)            = paramValue
      case ('Hgd Solution')
        Hgd_solution_option(groupIndex,:)  = paramValue
      case default
        ! did not find the parameter, return false
        SetHgIntegerParameter = .false.
    end select
  end function
  !===========================================================================================================================
  ! set pathway index
  logical function SetHgPathwayIndex(pathwayName, index)
    character(len=*), intent(in) :: pathwayName
    integer                      :: index
    character(len=maxChar)       :: name_tempt, name_tempt2
    !
    SetHgPathwayIndex = .false.
    ! (25)
    select case (pathwayName)
      case ('Atm <--> Hg0') 
        Hg0_Volatilization_index  = index
        SetHgPathwayIndex         = .true.
      case ('Hg0 --> HgII') 
        Hg0_Oxidation_index       = index
        SetHgPathwayIndex         = .true.
      !
      case ('Air --> HgII') 
        HgII_Air_Deposition_index = index
        SetHgPathwayIndex         = .true.
      case ('HgII --> Hg0') 
        HgII_Photoreduction_index = index
        SetHgPathwayIndex         = .true.
      case ('HgII --> MeHg') 
        HgII_Methylation_index    = index
        SetHgPathwayIndex         = .true.
      case ('HgII2 --> MeHg2') 
        HgII2_Methylation_index   = index
        SetHgPathwayIndex         = .true.
      case ('HgII --> Bed') 
        HgII_Settling_index       = index
        SetHgPathwayIndex         = .true.
      case ('HgII --> HgII2') 
        HgII_Deposition_index     = index
        SetHgPathwayIndex         = .true.
      case ('Bed --> HgII') 
        HgII_Resuspension_index   = index
        SetHgPathwayIndex         = .true.
      case ('HgII2 --> HgII') 
        HgII_Erosion_index        = index
        SetHgPathwayIndex         = .true.
      case ('HgII2 burial') 
        HgII_Burial_index         = index
        SetHgPathwayIndex         = .true.
      case ('HgII2 <--> HgII') 
        HgII_Transfer_index       = index
        SetHgPathwayIndex         = .true.
      case ('HgII <--> HgII2') 
        HgII2_Transfer_index      = index
        SetHgPathwayIndex         = .true.
      !
      case ('Air --> MeHg') 
        MeHg_Air_Deposition_index = index
        SetHgPathwayIndex         = .true.
      case ('Atm <--> MeHg') 
        MeHg_Volatilization_index = index
        SetHgPathwayIndex         = .true.
      case ('MeHg --> Hg0') 
        MeHg_Photoreduction_index = index
        SetHgPathwayIndex         = .true.
      case ('MeHg --> HgII') 
        MeHg_Demethylation_index  = index
        SetHgPathwayIndex         = .true.
      case ('MeHg2 --> HgII2') 
        MeHg2_Demethylation_index = index
        SetHgPathwayIndex         = .true.
      case ('MeHg --> Bed') 
        MeHg_Settling_index       = index
        SetHgPathwayIndex         = .true.
      case ('MeHg --> MeHg2') 
        MeHg_Deposition_index     = index
        SetHgPathwayIndex         = .true.
      case ('Bed --> MeHg') 
        MeHg_Resuspension_index   = index
        SetHgPathwayIndex         = .true.
      case ('MeHg2 --> MeHg') 
        MeHg_Erosion_index        = index
        SetHgPathwayIndex         = .true.
      case ('MeHg2 burial') 
        MeHg_Burial_index         = index
        SetHgPathwayIndex         = .true.
      case ('MeHg2 <--> MeHg') 
        MeHg_Transfer_index       = index
        SetHgPathwayIndex         = .true.
      case ('MeHg <--> MeHg2') 
        MeHg2_Transfer_index      = index
        SetHgPathwayIndex         = .true.
      case default
        ! did not find the pathway name, return false
        SetHgPathwayIndex = .false.
      end select
  end function
  !===========================================================================================================================
  ! temperature corrections
  subroutine HgTempCorrection()
    real(R8) :: TRK
    real(R8) :: KL, KG    ! volatilization transfer coefficients
    !
    TwaterK = TwaterC + 273.15
    TsedK   = TsedC   + 273.15
    !
    ! non equilibrium - water column   
    if (use_NonEquilibrium(1)) then
      TRK = Trade(r) + 273.15
      if (use_AlgaeSorbed(2)) then
        kadap_tc = kadap(r) * exp(Eaad(r) * 1000.0 * (TwaterK - TRK) / (gas_constant * TwaterK * TRK))
        kdaap_tc = kdaap(r) * exp(Eade(r) * 1000.0 * (TwaterK - TRK) / (gas_constant * TwaterK * TRK))
      end if
      if (use_POMSorbed(2)) then
        kadpom_tc = kadpom(r) * exp(Eaad(r) * 1000.0 * (TwaterK - TRK) / (gas_constant * TwaterK * TRK))
        kdapom_tc = kdapom(r) * exp(Eade(r) * 1000.0 * (TwaterK - TRK) / (gas_constant * TwaterK * TRK))
      end if
      do k = 1, nGS
        if (use_SolidSorbed(2,k)) then
          kadp_tc(k) = kadp(k,r) * exp(Eaad(r) * 1000.0 * (TwaterK - TRK) / (gas_constant * TwaterK * TRK))
          kdap_tc(k) = kdap(k,r) * exp(Eade(r) * 1000.0 * (TwaterK - TRK) / (gas_constant * TwaterK * TRK))
        end if
      end do
    end if
    !
    ! sediment layer
    if (use_NonEquilibrium(2)) then
      TRK = Trade(r) + 273.15
      if (use_POMSorbed(2)) then
        kadpom2_tc = kadpom2(r) * exp(Eaad(r) * 1000.0 * (TsedK - TRK) / (gas_constant * TsedK * TRK))
        kdapom2_tc = kdapom2(r) * exp(Eade(r) * 1000.0 * (TsedK - TRK) / (gas_constant * TsedK * TRK))
      end if
      do k = 1, nGS
        if (use_SolidSorbed(2,k)) then
          kadp2_tc(k) = kadp2(k,r) * exp(Eaad(r) * 1000.0 * (TsedK - TRK) / (gas_constant * TsedK * TRK))
          kdap2_tc(k) = kdap2(k,r) * exp(Eade(r) * 1000.0 * (TsedK - TRK) / (gas_constant * TsedK * TRK))
        end if
      end do
    end if
    !
    ! volatilization velocity
    do i = 1, 3, 2 
      if (vv_option(i,r) == 2) then
        KL = ka * (32.0 / MW(i,r))**0.25
        KG = 168.0 * wind_speed * (18.0 / MW(i,r))**0.25
        if (KG < 100.0) KG = 100.0
        vv(i,r)%rc20 = 1.0 / (1.0 / KL + gas_constant * TwaterK / KH(i,r) / KG)
      end if
      vv_tc(i) = Arrhenius_TempCorrection(vv(i,r), TwaterC)
    end do
    !
    ! transformations
    TRK     = Tr12(r) + 273.15
    k12_tc  = k12(r) * exp(Ea12(r) * 1000.0 * (TwaterK - TRK) / (gas_constant * TwaterK * TRK))
    !
    kd23_tc = Arrhenius_TempCorrection(kd23(r), TwaterC)
    if (use_DOCSorbed(2)) kdoc23_tc = Arrhenius_TempCorrection(kdoc23(r), TwaterC)
    ! 
    if (use_BedSediment) then    
      kso42_tc  = Arrhenius_TempCorrection(kso42(r), TsedC)
      kd32_2_tc = Arrhenius_TempCorrection(kd32_2(r), TsedC)
    end if
  end subroutine
  !===========================================================================================================================
  ! Compute Linear equilibrium partition concentrations for HgII and MeHg
  subroutine EquilibriumPartitionConc()
    real(R8) :: Rd
    !
    ! water column
    if (IsWaterCell) then
      Rd = 1.0
      if (use_DOCSorbed(i))        Rd = Rd + Kdoc(i,r) * DOC / 1.0E6
      if (use_AlgaeSorbed(i))      Rd = Rd + Kap(i,r)  * Apd / 1.0E6
      if (use_POMSorbed(i))        Rd = Rd + Kpom(i,r) * POM / 1.0E6
      do k = 1, nGS
        if (use_SolidSorbed(i,k))  Rd = Rd + Kp(i,k,r) * Solid(k) / 1.0E6
      end do
      !
      Hgd(i) = Hg(i) / Rd
      if (use_DOCSorbed(i)) then 
        Hgdoc(i) = Kdoc(i,r) * DOC / 1.0E6 / Rd * Hg(i)
      else
        Hgdoc(i) = 0.0
      end if
      if (use_AlgaeSorbed(i)) then
        Hgap(i)  = Kap(i,r) * Apd / 1.0E6 / Rd * Hg(i) 
      else
        Hgap(i)  = 0.0
      end if
      if (use_POMSorbed(i)) then
        Hgpom(i) = Kpom(i,r) * POM / 1.0E6 / Rd * Hg(i) 
      else
        Hgpom(i) = 0.0
      end if
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) then
          Hgp(i,k) = Kp(i,k,r) * Solid(k) / 1.0E6 / Rd * Hg(i)
        else
          Hgp(i,k) = 0.0
        end if
      end do
    else 
      !
      ! sediment layer
      Rd  = Por(r)
      if (use_DOCSorbed(i))        Rd = Rd + Kdoc2(i,r) * DOC2 / 1.0E6 * Por(r)
      if (use_POMSorbed(i))        Rd = Rd + Kpom2(i,r) * POM2 / 1.0E6
      do k = 1, nGS
        if (use_SolidSorbed(i,k))  Rd = Rd + Kp2(i,k,r) * Solid2(k) / 1.0E6
      end do
      !
      Cd2(i) = Por(r) / Rd * Hg2(i) 
      if (use_DOCSorbed(i)) then
        Cdoc2(i) = Kdoc2(i,r) * DOC2 / 1.0E6 * Por(r) / Rd * Hg2(i)
      else
        Cdoc2(i) = 0.0
      end if
      if (use_POMSorbed(i)) then
        Hgpom2(i) = Kpom2(i,r) * POM2 / 1.0E6 / Rd * Hg2(i)
      else
        Hgpom2(i) = 0.0
      end if
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) then
          Hgp2(i,k) = Kp2(i,k,r) * Solid2(k) / 1.0E6 / Rd * Hg2(i) 
        else
          Hgp2(i,k) = 0.0
        end if
      end do
    end if
  end subroutine
  !=========================================================================================================================== 
  ! compute Langmuire equilibrium partition concentrations for HgII and MeHg
  subroutine LangmuirPartitionConc()
    if (IsWaterCell) then
      !
      ! compute Hgd using Newton_Raphson and Bisection 
      if (Hgd_solution_option(i,r) == 1) then
        call NewtonRaphson(Hgd(i))
      else if (Hgd_solution_option(i,r) == 2) then
        call Bisection(Hgd(i))  
      end if
      !
      ! compute HgIIdoc, HgIIap, HgIIpom and HgIIp
      if (use_DOCSorbed(i)) then
        Hgdoc(i)    = DOC / 1.0E6 * Kdoc(i,r) * Hgd(i)
      else
        Hgdoc(i)    = 0.0
      end if
      if (use_AlgaeSorbed(i)) then
        Hgap(i)     = Apd * qcap(i,r) * Klap(i,r) * Hgd(i) / (1.0E3 + Klap(i,r) * Hgd(i))
      else
        Hgap(i)     = 0.0
      end if
      if (use_POMSorbed(i)) then
        Hgpom(i)    = POM * qcpom(i,r) * Klpom(i,r) * Hgd(i) / (1.0E3 + Klpom(i,r) * Hgd(i))
      else
        Hgpom(i)    = 0.0
      end if
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) then
          Hgp(i,k)  = Solid(k) * qcp(i,k,r) * Klp(i,k,r) * Hgd(i) / (1.0E3 + Klp(i,k,r) * Hgd(i))   
        else
          Hgp(i,k)  = 0.0
        end if
      end do
    else
      !
      ! sediment layer    
      if (Hgd_solution_option(i,r) == 1) then
        call NewtonRaphson(Cd2(i))
      else if (Hgd_solution_option(i,r) == 2) then
        call Bisection(Cd2(i))  
      end if
      !
      if (use_DOCSorbed(i)) then
        Cdoc2(i)    = DOC2 / 1.0E6 * Kdoc2(i,r) * Cd2(i)
      else
        Cdoc2(i)    = 0.0 
      end if
      if (use_POMSorbed(i)) then
        Hgpom2(i)   = POM2 * qcpom2(i,r) * Klpom2(i,r) * Cd2(i) / Por(r) / (1.0E3 + Klpom2(i,r) * Cd2(i) / Por(r))
      else
        Hgpom2(i)   = 0.0
      end if
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) then
          Hgp2(i,k) = Solid2(k) * qcp2(i,k,r) * Klp2(i,k,r) * Cd2(i) / Por(r) / (1.0E3 + Klp2(i,k,r) * Cd2(i) / Por(r))  
        else
          Hgp2(i,k) = 0.0
        end if
      end do
    end if
  end subroutine
  !=========================================================================================================================== 
  ! compute Freundlich equilibrium partition concentrations for HgII and MeHg
  subroutine FreundlichPartitionConc()
    ! water column
    if (IsWaterCell) then
      !
      ! compute Hgd using Newton_Raphson or Bisection
      if (Hgd_solution_option(i,r) == 1) then
        call NewtonRaphson(Hgd(i))
      else if (Hgd_solution_option(i,r) == 2) then
        call Bisection(Hgd(i))
      end if
      !
      ! compute Hgdoc, Hgap, Hgpom and Hgp
      if (use_DOCSorbed(i)) then
        Hgdoc(i)    = DOC / 1.0E6 * Kdoc(i,r) * Hgd(i)
      else
        Hgdoc(i)    = 0.0
      end if
      if (use_AlgaeSorbed(i)) then
        Hgap(i)     = Apd * 10**(-3 * bap(i,r)) * Kfap(i,r) * Hgd(i)**bap(i,r)
      else
        Hgap(i)     = 0.0
      end if
      if (use_POMSorbed(i)) then
        Hgpom(i)    = POM * 10**(-3 * bpom(i,r)) * Kfpom(i,r) * Hgd(i)**bpom(i,r)
      else
        Hgpom(i)    = 0.0
      end if
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) then
          Hgp(i,k)  = Solid(k) * 10**(-3 * bp(i,k,r)) * Kfp(i,k,r) * Hgd(i)**bp(i,k,r) 
        else
          Hgp(i,k)  = 0.0
        end if
      end do
    else
      !
      ! sediment layer
      if (Hgd_solution_option(i,r) == 1) then
        call NewtonRaphson(Cd2(i))
      else if (Hgd_solution_option(i,r) == 2) then
        call Bisection(Cd2(i))  
      end if
      !
      if (use_DOCSorbed(i)) then
        Cdoc2(i)    = DOC2 / 1.0E6 * Kdoc2(i,r) * Cd2(i)
      else
        Cdoc2(i)    = 0.0 
      end if
      !
      if (use_POMSorbed(i)) then
        Hgpom2(i)   = POM2 * 10**(-3 * bpom2(i,r)) * Kfpom2(i,r) * (Cd2(i) / Por(r))**bpom2(i,r)
      else
        Hgpom2(i)   = 0.0
      end if
      !
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) then
          Hgp2(i,k) = Solid2(k) * 10**(-3 * bp2(i,k,r)) * Kfp2(i,k,r) * (Cd2(i) / Por(r))**bp2(i,k,r) 
        else
          Hgp2(i,k) = 0.0
        end if
      end do
    end if
  end subroutine
  !=========================================================================================================================== 
  ! compute Non-equilibrium partition concentrations for HgII
  ! only compute at t = 0.
  subroutine NonEquilibriumPartitionConc()
    real(R8) :: Rd
    !
    ! water column 
    if (IsWaterCell) then
      if (abs(t) < 1.0E-10) then
        Rd = 1.0
        if (use_DOCSorbed(2))         Rd = Rd + Kdoc(2,r) * DOC / 1.0E6
        if (use_AlgaeSorbed(2))       Rd = Rd + kadap(r) / kdaap(r) * qcap(2,r) * Apd / 1.0E3
        if (use_POMSorbed(2))         Rd = Rd + kadpom(r) / kdapom(r) * qcpom(2,r) * POM / 1.0E3
        do k = 1, nGS
          if (use_SolidSorbed(2,k))   Rd = Rd + kadp(k,r) / kdap(k,r) * qcp(2,k,r) * Solid(k) / 1.0E3
        end do
        !
        Hgd(2) = Hg(2) / Rd
        if (use_DOCSorbed(2)) then 
          Hgdoc(2)   = Kdoc(2,r) * DOC / 1.0E6 / Rd * Hg(2)
        else
          Hgdoc(2)   = 0.0
        end if
        if (use_AlgaeSorbed(2)) then
          Hgap(2)    = kadap(r) / kdaap(r) * qcap(2,r) * Apd / 1.0E3 / Rd * Hg(2)
        else
          Hgap(2)    = 0.0
        end if
        if (use_POMSorbed(2)) then
          Hgpom(2)   = kadpom(r) / kdapom(r) * qcpom(2,r) * POM / 1.0E3 / Rd * Hg(2)
        else
          Hgpom(2)   = 0.0
        end if
        do k = 1, nGS
          if (use_SolidSorbed(2,k)) then
            Hgp(2,k) = kadp(k,r) / kdap(k,r) * qcp(2,k,r) * Solid(k) / 1.0E3 / Rd * Hg(2)
          else
            Hgp(2,k) = 0.0
          end if
        end do
      end if
    else
      !
      ! sediment layer
      if (abs(t) < 1.0E-10) then
          Rd  = Por(r)
          if (use_DOCSorbed(2))         Rd = Rd + Kdoc2(2,r) * DOC2 / 1.0E6 * Por(r)
          if (use_POMSorbed(2))         Rd = Rd + kadpom2(r) / kdapom2(r) * qcpom2(2,r) * POM2 / 1.0E3
          do k = 1, nGS
            if (use_SolidSorbed(2,k))   Rd = Rd + kadp2(k,r) / kdap2(k,r) * qcp2(2,k,r) * Solid2(k) / 1.0E3
          end do
          !
          Hgd2(2) = 1.0 / Rd * Hg2(2)
          if (use_DOCSorbed(2)) then
            Hgdoc2(2)   = Kdoc2(2,r) * DOC2 / 1.0E6 / Rd * Hg2(2)
          else
            Hgdoc2(2)   = 0.0
          end if
          if (use_POMSorbed(2)) then
            Hgpom2(2)   = kadpom2(r) / kdapom2(r) * qcpom2(2,r) * POM2 / 1.0E3 / Rd * Hg2(2)
          else
            Hgpom2(2)   = 0.0
          end if
          do k = 1, nGS
            if (use_SolidSorbed(2,k)) then
              Hgp2(2,k) = kadp2(k,r) / kdap2(k,r) * qcp2(2,k,r) * Solid2(k) / 1.0E3 / Rd * Hg2(2)
            else
              Hgp2(2,k) = 0.0
            end if
          end do
      end if
      !
                              Cd2(2) = Hgd2(2) * Por(r)
      if (use_DOCSorbed(2)) Cdoc2(2) = Hgdoc2(2) * Por(r)
    end if
  end subroutine
  !=========================================================================================================================== 
  ! function
  double precision function f(x)
    real(R8), intent(in) :: x
    ! 
    f = x
    if (IsWaterCell) then
      f = f - Hg(i)
      if (use_DOCSorbed(i))         f = f + DOC / 1.0E6 * Kdoc(i,r) * x 
      if (use_Langmuir(i,1)) then
        if (use_AlgaeSorbed(i))     f = f + Apd * qcap(i,r) * Klap(i,r) * x / (1.0E3 + Klap(i,r) * x)
        if (use_POMSorbed(i))       f = f + POM * qcpom(i,r) * Klpom(i,r) * x / (1.0E3 + Klpom(i,r) * x)
        do k = 1, nGS
          if (use_SolidSorbed(i,k)) f = f + Solid(k) * qcp(i,k,r) * Klp(i,k,r) * x / (1.0E3 + Klp(i,k,r) * x)
        end do
      else if (use_Freundlich(i,1)) then
        if (use_AlgaeSorbed(i))     f = f + Apd * 10**(-3 * bap(i,r)) * Kfap(i,r) * x**bap(i,r)
        if (use_POMSorbed(i))       f = f + POM * 10**(-3 * bpom(i,r)) * Kfpom(i,r) * x**bpom(i,r)
        do k = 1, nGS
          if (use_SolidSorbed(i,k)) f = f + Solid(k) * 10**(-3 * bp(i,k,r)) * Kfp(i,k,r) * x**bp(i,k,r) 
        end do
      end if
    else 
      !
      ! sediment layer
      f = f - Hg2(i)
      if (use_DOCSorbed(i))         f = f + DOC2 / 1.0E6 * Kdoc2(i,r) * x
      if (use_Langmuir(i,2)) then
        if (use_POMSorbed(i))       f = f + POM2 * qcpom2(i,r) * Klpom2(i,r) * x / Por(r) /  &
                                        (1.0E3 + Klpom2(i,r) * x / Por(r))
        do k = 1, nGS
          if (use_SolidSorbed(i,k)) f = f + Solid2(k) * qcp2(i,k,r) * Klp2(i,k,r) * x /  & 
                                        Por(r) / (1.0E3 + Klp2(i,k,r) * x / Por(r))
        end do
      else if (use_Freundlich(i,2)) then
        if (use_POMSorbed(i))       f = f + POM2 * 10**(-3 * bpom2(i,r)) * Kfpom2(i,r) *  & 
                                        (x / Por(r))**bpom2(i,r)
        do k = 1, nGS
          if (use_SolidSorbed(i,k)) f = f + Solid2(k) * 10**(-3 * bp2(i,k,r)) * Kfp2(i,k,r) *  & 
                                        (x / Por(r))**bp2(i,k,r) 
        end do 
      end if
    end if
  end function 
  !=========================================================================================================================== 
  ! function of df(HgIId)
  double precision function df(x)
    real(R8), intent(in) :: x
    ! 
    df = 1
    if (IsWaterCell) then
      if (use_DOCSorbed(i))            df = df + DOC / 1.0E6 * Kdoc(i,r) 
      if (use_Langmuir(i,1)) then
        if (use_AlgaeSorbed(i))        df = df + Apd * qcap(i,r) * Klap(i,r) * 1.0E3 / (1.0E3 + Klap(i,r) * x)**2.0
        if (use_POMSorbed(i))          df = df + POM * qcpom(i,r) * Klpom(i,r) * 1.0E3 / (1.0E3 + Klpom(i,r) * x)**2.0
        do k = 1, nGS
          if (use_SolidSorbed(i,k))    df = df + Solid(k) * qcp(i,k,r) * Klp(i,k,r) * 1.0E3 /  & 
                                            (1.0E3 + Klp(i,k,r) * x)**2.0
        end do
      else if (use_Freundlich(i,1)) then
        if (use_AlgaeSorbed(i))        df = df + Apd * 10**(-3 * bap(i,r)) * Kfap(i,r) * bap(i,r) * x**(bap(i,r) - 1.0)
        if (use_POMSorbed(i))          df = df + POM * 10**(-3 * bpom(i,r)) * Kfpom(i,r) * bpom(i,r) *  &
                                            x**(bpom(i,r) - 1.0)
        do k = 1, nGS
          if (use_SolidSorbed(i,k))    df = df + Solid(k) * 10**(-3 * bp(i,k,r)) * Kfp(i,k,r) *  &
                                            bp(i,k,r) * x**(bp(i,k,r) - 1.0)
        end do
      end if
    else 
      !
      ! sediment layer
      if (use_DOCSorbed(i))            df = df + DOC2 / 1.0E6 * Kdoc2(i,r)
      if (use_Langmuir(i,2)) then
        if (use_POMSorbed(i))          df = df + POM2 * qcpom2(i,r) * Klpom2(i,r) / Por(r) * 1.0E3 /  &
                                            (1.0E3 + Klpom2(i,r) * x / Por(r))**2.0
        do k = 1, nGS
          if (use_SolidSorbed(i,k))    df = df + Solid2(k) * qcp2(i,k,r) * Klp2(i,k,r) / Por(r) * 1.0E3 /  &
                                            (1.0E3 + Klp2(i,k,r) * x / Por(r))**2.0
        end do
      else if (use_Freundlich(i,2)) then
        if (use_POMSorbed(i))          df = df + POM2 * 10**(-3 * bpom2(i,r)) * Kfpom2(i,r) * bpom2(i,r) *  & 
                                            (x / Por(r))**(bpom2(i,r) - 1.0) / Por(r)
        do k = 1, nGS
          if (use_SolidSorbed(i,k))    df = df + Solid2(k) * 10**(-3 * bp2(i,k,r)) * Kfp2(i,k,r) * bp2(i,k,r) *  & 
                                            (x / Por(r))**(bp2(i,k,r) - 1.0) / Por(r)
        end do    
      end if
    end if
  end function 
  !================================================================================================================================
  subroutine HgPartitions() 
    integer :: j
    ! 
    do i = 2, 3
      do j = 1, 2
        if (j == 1) then
          IsWaterCell = .true.
        else
          IsWaterCell = .false.
        end if
        !
        if (use_Equilibrium(i,j)) then
          call EquilibriumPartitionConc()
        else if (use_Langmuir(i,j)) then 
          if (use_AlgaeSorbed(i) .or. use_POMSorbed(i) .or. use_AnySolidSorbed(i)) then
            call LangmuirPartitionConc()
          else
            ! if only DOCSorbed is selected, use linear equilibrium model
            call EquilibriumPartitionConc()
          end if
        else if (use_Freundlich(i,j)) then
          if (use_AlgaeSorbed(i) .or. use_POMSorbed(i) .or. use_AnySolidSorbed(i)) then
            call FreundlichPartitionConc()
          else
            ! if only DOCSorbed is selected, use linear equilibrium model
            call EquilibriumPartitionConc()
          end if
        else if (i == 2 .and. use_NonEquilibrium(j)) then
          call NonEquilibriumPartitionConc() 
        end if
      end do
    end do
  end subroutine
  !================================================================================================================================  
  subroutine HgPathways()
    real(R8) :: kviscosity
    real(R8) :: Kd2_avg, Tsolid2
    real(R8) :: lambdamax
    !
    !----------------------------------
    ! prepare for pathway computation !
    !----------------------------------
    ! get lambdamax
    lambdamax = alpha(r) * lambda
    do i = 2, 3
      Iav(i) = 1.33 * q_solar / I0pht(i,r) * (1.0 - exp(- lambdamax * depth)) / (lambdamax * depth)  &  
               * (1.0 - 0.56 * cloudiness)
    end do
    if (use_Bedsediment) SO4_reduction = kso42_tc * rmso4(r) * SO42 / (SO42 + Kso4(r)) * SO42 
    !
    ! Hg0 pathways
    !---------------------
    ! Hg0 volatilization !
    !---------------------
    Hg0_Volatilization = vv_tc(1) / depth * (Hg(1) - Hg00 / (KH(1,r) / (gas_constant * TwaterK)))
    !----------------
    ! Hg0 oxidation !
    !----------------
    Hg0_Oxidation = k12_tc * Hg(1)
    ! 
    ! HgII pathways
    !----------------------
    ! HgII air deposition !
    !----------------------
    HgII_Air_Deposition = LHgII * surface_area / volume
    !--------------------------------------------
    ! HgII settling and resuspension and burial !
    !--------------------------------------------
    if (use_POMSorbed(2)) then
      HgII_Settling     = Hgpom(2) * vsom(r) / depth
      if (use_BedSediment)  then
        HgII_Deposition = Hgpom(2) * vsom(r) / h2(r)
        HgII_Burial     = Hgpom2(2) * vb / h2(r)
      end if
    else
      HgII_Settling   = 0.0
      if (use_BedSediment) then
        HgII_Deposition = 0.0
        HgII_Burial     = 0.0
      end if
    end if
    !
    if (use_AlgaeSorbed(2)) then
      HgII_Settling                         = HgII_Settling    + Hgap(2) * vsap(r) / depth
      if (use_BedSediment)  HgII_Deposition = HgII_Deposition  + Hgap(2) * vsap(r) / h2(r)
    end if
    !
    if (use_BedSediment) then
      HgII_Resuspension  = 0.0
      HgII_Erosion       = 0.0
    end if
    do k = 1, nGS
      if (use_SolidSorbed(2,k)) then
        HgII_Settling          = HgII_Settling      + Hgp(2,k)  * vsp(k) / depth
        if (use_BedSediment) then
          HgII_Deposition      = HgII_Deposition    + Hgp(2,k)  * vsp(k) / h2(r)
          HgII_Resuspension    = HgII_Resuspension  + Hgp2(2,k) * vrp(k) / depth
          HgII_Erosion         = HgII_Erosion       + Hgp2(2,k) * vrp(k) / h2(r) 
          HgII_Burial          = HgII_Burial        + Hgp2(2,k) * vb / h2(r)
        end if
      end if
    end do
    !-----------------
    ! HgII reduction !
    !-----------------
    if (use_DOCSorbed(2)) then
      HgII_Photoreduction  = (Hgd(2) * kd21(r) + Hgdoc(2) * kdoc21(r)) * Iav(2)
    else
      HgII_Photoreduction  =  Hgd(2) * kd21(r) * Iav(2)
    end if 
    !-------------------
    ! HgII methylation !
    !-------------------
    if (use_DOCSorbed(2)) then
      HgII_Methylation = Hgd(2) * kd23_tc + Hgdoc(2) * kdoc23_tc  
    else
      HgII_Methylation = Hgd(2) * kd23_tc 
    end if
    if (use_BedSediment) HgII2_Methylation = Cd2(2) * SO4_reduction 
    !
    ! MeHg pathways
    !----------------------
    ! MeHg air deposition !
    !----------------------
    MeHg_Air_Deposition = LMeHg * surface_area / volume
    !--------------------------------------------
    ! MeHg settling and resuspension and burial !
    !--------------------------------------------
    if (use_POMSorbed(3)) then
      MeHg_Settling     = Hgpom(3)  * vsom(r) / depth
      if (use_BedSediment)  then
        MeHg_Deposition = Hgpom(3)  * vsom(r) / h2(r)
        MeHg_Burial     = Hgpom2(3) * vb / h2(r)
      end if
    else
      MeHg_Settling     = 0.0
      if (use_BedSediment) then
        MeHg_Deposition = 0.0
        MeHg_Burial     = 0.0
      end if
    end if
    !
    if (use_AlgaeSorbed(3)) then
      MeHg_Settling                         = MeHg_Settling   + Hgap(3) * vsap(r) / depth
      if (use_BedSediment) MeHg_Deposition  = MeHg_Deposition + Hgap(3) * vsap(r) / h2(r)
    end if
    !
    if (use_BedSediment) then
      MeHg_Resuspension  = 0.0
      MeHg_Erosion       = 0.0
    end if
    !
    do k = 1, nGS
      if (use_SolidSorbed(3,k)) then
        MeHg_Settling        = MeHg_Settling     + Hgp(3,k)  * vsp(k) / depth
        if (use_BedSediment) then
          MeHg_Deposition    = MeHg_Deposition   + Hgp(3,k)  * vsp(k) / h2(r)
          MeHg_Resuspension  = MeHg_Resuspension + Hgp2(3,k) * vrp(k) / depth
          MeHg_Erosion       = MeHg_Erosion      + Hgp2(3,k) * vrp(k) / h2(r) 
          MeHg_Burial        = MeHg_Burial       + Hgp2(3,k) * vb / h2(r)
        end if
      end if
    end do
    !----------------------
    ! MeHg volatilization !
    !----------------------
    MeHg_Volatilization = vv_tc(3) / depth * (Hgd(3) - MeHg0 / (KH(3,r) / (gas_constant * TwaterK)))
    !---------------------
    ! MeHg demethylation !
    !---------------------
    if (use_DOCSorbed(3)) then
      MeHg_Demethylation = (Hgd(3) * kd32(r) + Hgdoc(3) * kdoc32(r)) * Iav(3)   
    else
      MeHg_Demethylation = Hgd(3) * kd32(r) * Iav(3) 
    end if
    if (use_BedSediment) MeHg2_Demethylation = Cd2(3) * kd32_2_tc 
    !----------------------
    ! MeHg photoreduction !
    !----------------------
    if (use_DOCSorbed(3)) then
      MeHg_Photoreduction = (Hgd(3) * kd31(r) + Hgdoc(3) * kdoc31(r)) * Iav(3)
    else
      MeHg_Photoreduction = Hgd(3) * kd31(r) * Iav(3)
    end if
    !------------------------------
    ! HgII and MeHg mass transfer !
    !------------------------------
    if (use_BedSediment) then 
      do i = 2, 3
        if (vm_option(i,r) == 2) then
          !
          ! Formula given by Thibodeaux et al. (2001): bioturbation-driven transport, sorb to sediment
          if (use_Equilibrium(i,2)) then
            Kd2_avg = 0.0
            Tsolid2 = 0.0
            if (use_POMSorbed(i)) then
              Kd2_avg = Kd2_avg + Kpom2(i,r) * POM2
              Tsolid2 = Tsolid2 + POM2
            end if
            do k = 1, nGS
              if (use_SolidSorbed(i,k)) then
                Kd2_avg = Kd2_avg + Kp2(i,k,r) * Solid2(k)
                Tsolid2 = Tsolid2 + Solid2(k)
              end if
            end do
            if (Tsolid2 > 0.0) Kd2_avg  = Kd2_avg / Tsolid2
            if (Kd2_avg > 0.0) vm(i,r)  = 0.01 / (1.0 / beta(r) + z2(r) / (Db(r) * ((1.0 - Por(r)) * ps(r)) * Kd2_avg))
          end if
        else if (vm_option(i,r) == 3) then
          !
          ! Boyer et al.(1994)
          vm(i,r) = Por(r)**3.0 * Dm(i,r) / 0.01
        else if (vm_option(i,r) == 4) then
          !
          ! Di Toro (1981)
          vm(i,r) = 0.19 * Por(r) / MW(i,r)**(2.0 / 3.0)
        else if (vm_option(i,r) == 5) then
          !
          ! Schink and Guinasso (1977)
          ! ustar     = sqrt(9.81 * hydraulic_radius * slope)
          ! Twaterc is used, TsedC mayb be more accurate???
          kviscosity = 1.79E-6 / (1.0 + 0.03368 * TwaterC + 0.000221 * TwaterC**2.0)
          ! vm(i,r)   = ustar * 86400.0 * (Dm(i,r) / 86400.0 / kviscosity)**(2.0 / 3.0) / 24.0
          ! change ustar to shear_velocity
          vm(i,r) = shear_velocity * 86400.0 * (Dm(i,r) / 86400.0 / kviscosity)**(2.0 / 3.0) / 24.0
        end if 
      end do
      !
      HgII_Transfer    = vm(2,r) * (Cd2(2) / Por(r) - Hgd(2)) / depth
      HgII2_Transfer   = vm(2,r) * (Hgd(2) - Cd2(2) / Por(r)) / h2(r)
      if (use_DOCSorbed(2)) then
        HgII_Transfer  = HgII_transfer  + vm(2,r) * (Cdoc2(2) / Por(r) - Hgdoc(2)) / depth
        HgII2_Transfer = HgII2_transfer - vm(2,r) * (Cdoc2(2) / Por(r) - Hgdoc(2)) / h2(r)
      end if
      MeHg_Transfer    = vm(3,r) * (Cd2(3) / Por(r) - Hgd(3)) / depth
      MeHg2_Transfer   = vm(3,r) * (Hgd(3) - Cd2(3) / Por(r)) / h2(r)
      if (use_DOCSorbed(3)) then
        MeHg_Transfer  = MeHg_Transfer  + vm(3,r) * (Cdoc2(3) / Por(r) - Hgdoc(3)) / depth
        MeHg2_Transfer = MeHg2_Transfer - vm(3,r) * (Cdoc2(3) / Por(r) - Hgdoc(3)) / h2(r)
      end if
      !
    end if
  end subroutine
  !=========================================================================================================================== 
  ! compute kinetic rate changes for Hg0
  subroutine Hg0Kinetics()
    !
    dHgdt(1)  = 0.0
    !
    dHgdt(1) = -Hg0_Volatilization - Hg0_Oxidation + HgII_Photoreduction * Y21(r) + MeHg_Photoreduction * Y31(r)    
    !
  end subroutine   
  !=========================================================================================================================== 
  ! compute kinetic rate changes for HgII
  subroutine HgIIKinetics()
    real(R8) :: kad_tmp, Adsorption_Desorption
    !
    dHgdt(2)    = 0.0
    dHgddt(2)   = 0.0
    dHgdocdt(2) = 0.0
    dHgpomdt(2) = 0.0
    dHgapdt(2)  = 0.0
    dHgpdt(2,:) = 0.0
    if (use_BedSediment) then
      dHg2dt(2)    = 0.0  
      dHgd2dt(2)   = 0.0
      dHgdoc2dt(2) = 0.0
      dHgpom2dt(2) = 0.0
      dHgp2dt(2,:) = 0.0
    end if
    !
    ! equilibrium kinetic rate equations
    if (.not. use_NonEquilibrium(1)) then
      dHgdt(2)  = HgII_Air_Deposition - HgII_Settling + Hg0_Oxidation * Y12(r)  &  
                  - HgII_Photoreduction - HgII_Methylation + MeHg_Demethylation * Y32(r) 
      if (use_BedSediment)  dHgdt(2)  = dHgdt(2)  + HgII_Resuspension + HgII_Transfer
    end if
    if (.not. use_NonEquilibrium(2)) then
      if (use_BedSediment) dHg2dt(2) = dHg2dt(2) + HgII_Deposition - HgII_Erosion + HgII2_Transfer  & 
                                       - HgII_Burial- HgII2_Methylation + MeHg2_Demethylation * Y32(r)
    end if
    !
    ! non-equilibrium kinetic rate equations
    if (use_NonEquilibrium(1) .or. use_NonEquilibrium(2)) then
      if (use_NonEquilibrium(1)) then  
        if (Hg(2)>0.0) then
          dHgddt(2)   = dHgddt(2)   + HgII_Air_Deposition * Hgd(2)   / Hg(2)
          if (use_DOCSorbed(2))       dHgdocdt(2) = dHgdocdt(2) + HgII_Air_Deposition * Hgdoc(2) / Hg(2)
          if (use_POMSorbed(2))       dHgpomdt(2) = dHgpomdt(2) + HgII_Air_Deposition * Hgpom(2) / Hg(2)
          if (use_AlgaeSorbed(2))     dHgapdt(2)  = dHgapdt(2)  + HgII_Air_Deposition * Hgap(2)  / Hg(2)
          do k = 1, nGS
            if (use_SolidSorbed(2,k)) dHgpdt(2,k) = dHgpdt(2,k) + HgII_Air_Deposition * Hgp(2,k) / Hg(2)
          end do
        else
          dHgddt(2) = dHgddt(2) + HgII_Air_Deposition
        end if
      end if
      !------------------------------------ 
      ! settling, resuspension and burial !
      !------------------------------------
      if (use_AlgaeSorbed(2)) then     
        if (use_NonEquilibrium(1))                         dHgapdt(2)   = dHgapdt(2)   - Hgap(2) * vsap(r) / depth
        if (use_NonEquilibrium(2) .and. use_POMSorbed(2))  dHgpom2dt(2) = dHgpom2dt(2) + Hgap(2) * vsap(r) / h2(r) 
      end if
      !
      if (use_POMSorbed(2)) then
        if (use_NonEquilibrium(1))                         dHgpomdt(2)  = dHgpomdt(2)  - Hgpom(2) * vsom(r) / depth
        if (use_NonEquilibrium(2))                         dHgpom2dt(2) = dHgpom2dt(2) + (Hgpom(2) * vsom(r) - Hgpom2(2) * vb) / h2(r) 
      end if
      !
      do k = 1, nGS
        if (use_SolidSorbed(2,k)) then 
          if (use_NonEquilibrium(1)) then 
            dHgpdt(2,k)  = dHgpdt(2,k)  - Hgp(2,k)  * vsp(k) / depth
            if (use_BedSediment)     dHgpdt(2,k)  = dHgpdt(2,k)  + Hgp2(2,k) * vrp(k) / depth
          end if
          if (use_NonEquilibrium(2)) dHgp2dt(2,k) = dHgp2dt(2,k) + (Hgp(2,k) * vsp(k) - Hgp2(2,k) * (vrp(k) + vb)) / h2(r)
        end if 
      end do
      !---------------------
      ! HgII mass transfer !
      !---------------------
      if (use_BedSediment) then
        if (use_NonEquilibrium(1)) then  
          dHgddt(2)    = dHgddt(2)    + (Cd2(2) / Por(r) - Hgd(2)) * vm(2,r) / depth
          if (use_DOCSorbed(2))  dHgdocdt(2)  = dHgdocdt(2)  + (Cdoc2(2) / Por(r) - Hgdoc(2)) * vm(2,r) / depth
        end if
        if (use_NonEquilibrium(2)) then   
          dHgd2dt(2)   = dHgd2dt(2)   - (Cd2(2) / Por(r) - Hgd(2)) * vm(2,r) / h2(r)
          if (use_DOCSorbed(2))  dHgdoc2dt(2) = dHgdoc2dt(2) - (Cdoc2(2) / Por(r) - Hgdoc(2)) * vm(2,r) / h2(r)
        end if
      end if
      !----------------
      ! Hg0 oxidation !
      !----------------
      if (use_NonEquilibrium(1))   dHgddt(2)  = dHgddt(2) + k12_tc * Y12(r) * Hg(1)
      !-----------------------
      ! HgII photoreduction !
      !-----------------------
      if (use_NonEquilibrium(1)) then
        dHgddt(2)    = dHgddt(2)   - Hgd(2) * kd21(r) * Iav(2)
        if (use_DOCSorbed(2))    dHgdocdt(2)  = dHgdocdt(2) - Hgdoc(2) * kdoc21(r) * Iav(2)
      end if 
      !-------------------
      ! HgII methylation !
      !-------------------
      if (use_NonEquilibrium(1)) then
        dHgddt(2)    = dHgddt(2)   - Hgd(2) * kd23_tc
        if (use_DOCSorbed(2))    dHgdocdt(2)  = dHgdocdt(2) - Hgdoc(2) * kdoc23_tc 
      end if
      if (use_NonEquilibrium(2)) dHgd2dt(2)   = dHgd2dt(2)  - Cd2(2) * SO4_reduction  
      !---------------------
      ! MeHg demethylation !
      !---------------------
      if (use_NonEquilibrium(1)) then
        dHgddt(2)    = dHgddt(2) + Hgd(3) * kd32(r) * Iav(3) * Y32(r)
        if (use_DOCSorbed(3))    dHgddt(2)    = dHgddt(2) + Hgdoc(3) * kdoc32(r) * Iav(3) * Y32(r) 
      end if
      if (use_NonEquilibrium(2)) dHgd2dt(2)   = dHgd2dt(2) + Cd2(3) * kd32_2_tc * Y32(r)
      !-------------------------------------------------------------------------
      ! non_equilibrium partition                                              !
      !-------------------------------------------------------------------------
      if (use_NonEquilibrium(1)) then
        if (use_AlgaeSorbed(2)) then
          kad_tmp = kadap_tc * 1.0E-3 * (Apd * qcap(2,r) - Hgap(2)) 
          if (Apd * qcap(2,r) > Hgap(2) .and. kad_tmp * dt < 1.0) then
            Adsorption_Desorption = - kdaap_tc * Hgap(2) + kad_tmp * Hgd(2) 
          else if (Apd * qcap(2,r) < Hgap(2)) then
            Adsorption_Desorption = - kdaap_tc * Hgap(2)
          else if (kad_tmp * dt >= 1.0) then
            Adsorption_Desorption = - kdaap_tc * Hgap(2) + Hgd(2) / dt
          end if
          dHgapdt(2) = dHgapdt(2) + Adsorption_Desorption
          dHgddt(2)  = dHgddt(2)  - Adsorption_Desorption 
        end if
        ! HgIIpom must be less than POM * qcpom(r) 
        if (use_POMSorbed(2)) then
          kad_tmp = kadpom_tc * 1.0E-3 * (POM * qcpom(2,r) - Hgpom(2)) 
          if (POM * qcpom(2,r) > Hgpom(2) .and. kad_tmp * dt < 1.0) then
            Adsorption_Desorption = - kdapom_tc * Hgpom(2) + kad_tmp * Hgd(2) 
          else if (POM * qcpom(2,r) <= Hgpom(2)) then
            ! no available adsorption
            Adsorption_Desorption = - kdapom_tc * Hgpom(2)
          else if (kad_tmp * dt >= 1.0) then
            ! maximum adsorptable HgII is HgIId
            Adsorption_Desorption = - kdapom_tc * Hgpom(2) + Hgd(2) / dt
          end if
          dHgpomdt(2)    = dHgpomdt(2) + Adsorption_Desorption
          dHgddt(2)      = dHgddt(2)   - Adsorption_Desorption
        end if
        !
        do k = 1, nGS
          if (use_SolidSorbed(2,k)) then
            kad_tmp = kadp_tc(k) * 1.0E-3 * (Solid(k) * qcp(2,k,r) - Hgp(2,k))
            if (Solid(k) * qcp(2,k,r) > Hgp(2,k) .and. kad_tmp * dt < 1.0) then
              Adsorption_Desorption = - kdap_tc(k) * Hgp(2,k) + kad_tmp * Hgd(2)
            else if (Solid(k) * qcp(2,k,r) <= Hgp(2,k)) then
              Adsorption_Desorption = - kdap_tc(k) * Hgp(2,k)
            else if (kad_tmp * dt >= 1.0) then
              Adsorption_Desorption = - kdap_tc(k) * Hgp(2,k) + Hgd(2) / dt
            end if
            dHgpdt(2,k)  = dHgpdt(2,k) + Adsorption_Desorption
            dHgddt(2)    = dHgddt(2)   - Adsorption_Desorption
          end if
        end do
      end if
      !
      if (use_NonEquilibrium(2)) then 
        if (use_POMSorbed(2)) then   
          kad_tmp = kadpom2_tc / Por(r) * 1.0E-3 * (POM2 * qcpom2(2,r) - Hgpom2(2))
          if (POM2 * qcpom2(2,r) > Hgpom2(2) .and. kad_tmp * dt < 1.0) then 
            Adsorption_Desorption = - kdapom2_tc * Hgpom2(2) + kad_tmp * Cd2(2) 
          else if (POM2 * qcpom2(2,r) <= Hgpom2(2)) then
            Adsorption_Desorption = - kdapom2_tc * Hgpom2(2)
          else if (kad_tmp * dt >= 1.0) then
            Adsorption_Desorption = - kdapom2_tc * Hgpom2(2) + Cd2(2) / dt
          end if
          dHgpom2dt(2) = dHgpom2dt(2) + Adsorption_Desorption
          dHgd2dt(2)   = dHgd2dt(2)   - Adsorption_Desorption
        end if
        !
        do k = 1, nGS
          if (use_SolidSorbed(2,k)) then
            kad_tmp = kadp2_tc(k) * 1.0E-3 / Por(r) * (Solid2(k) * qcp2(2,k,r) - Hgp2(2,k))
            if (Solid2(k) * qcp2(2,k,r) > Hgp2(2,k) .and. kad_tmp * dt < 1.0) then
              Adsorption_Desorption = - kdap2_tc(k) * Hgp2(2,k) + kad_tmp * Cd2(2) 
            else if (Solid2(k) * qcp2(2,k,r) <= Hgp2(2,k)) then
              Adsorption_Desorption = - kdap2_tc(k) * Hgp2(2,k)
            else if (kad_tmp * dt >= 1.0) then
              Adsorption_Desorption = - kdap2_tc(k) * Hgp2(2,k) + Cd2(2) / dt
            end if
            dHgp2dt(2,k)  = dHgp2dt(2,k) + Adsorption_Desorption
            dHgd2dt(2)    = dHgd2dt(2)   - Adsorption_Desorption
          end if
        end do 
        !  
      end if 
    end if
    !
    if (use_BedSediment) then 
      if (use_NonEquilibrium(2)) then    
        Hgd2(2)   = max(Hgd2(2)   + dHgd2dt(2)   * dt / Por(r), 0.0)
        if (use_DOCSorbed(2))       Hgdoc2(2) = max(Hgdoc2(2) + dHgdoc2dt(2) * dt / Por(r), 0.0)
        if (use_POMSorbed(2))       Hgpom2(2) = max(Hgpom2(2) + dHgpom2dt(2) * dt, 0.0)
        do k = 1, nGS 
          if (use_SolidSorbed(2,k)) Hgp2(2,k) = max(Hgp2(2,k) + dHgp2dt(2,k) * dt, 0.0)
        end do
      else
        Hg2(2) = max(Hg2(2) + dHg2dt(2) * dt, 0.0)
      end if
    end if
  end subroutine
  !=========================================================================================================================== 
  ! compute kinetic rate changes for MeHg
  subroutine MeHgKinetics()
    dHgdt(3) = 0.0
    if (use_BedSediment) dHg2dt(3) = 0.0
    !
    dHgdt(3) = MeHg_Air_Deposition - MeHg_Settling - MeHg_Volatilization - MeHg_Demethylation - MeHg_Photoreduction + HgII_Methylation * Y23(r)     
    if (use_BedSediment) then
      dHgdt(3)  = dHgdt(3) + MeHg_resuspension + MeHg_transfer
      dHg2dt(3) = MeHg_Deposition - MeHg_Erosion - MeHg2_Demethylation + MeHg2_Transfer - MeHg_Burial + HgII2_Methylation * Y23(r)
      Hg2(3)    = max(Hg2(3) + dHg2dt(3) * dt, 0.0)
    end if 
  end subroutine    
  !=========================================================================================================================== 
  ! call subroutines   
  subroutine ComputeHgKinetics()
    ! 
    call HgPartitions()
    !
    call HgTempCorrection()
    call HgPathways()
    call Hg0Kinetics()
    call HgIIKinetics() 
    call MeHgKinetics()
  end subroutine
  !===========================================================================================================================
  ! output pathway
  subroutine HgPathwayOutput(na, a)
    integer  :: na
    real(R8) :: a(na)
    !
    if (Hg0_Volatilization_index > 0)           a(Hg0_Volatilization_index)          = Hg0_Volatilization
    if (Hg0_Oxidation_index > 0)                a(Hg0_Oxidation_index)               = Hg0_Oxidation
    if (HgII_Photoreduction_index > 0)          a(HgII_Photoreduction_index)         = HgII_Photoreduction
    if (HgII_Methylation_index > 0)             a(HgII_Methylation_index)            = HgII_Methylation
    if (HgII_Air_Deposition_index > 0)          a(HgII_Air_Deposition_index)         = HgII_Air_Deposition
    if (HgII_Settling_index > 0)                a(HgII_Settling_index)               = HgII_Settling
    if (MeHg_Air_Deposition_index > 0)          a(MeHg_Air_Deposition_index)         = MeHg_Air_Deposition
    if (MeHg_Volatilization_index > 0)          a(MeHg_Volatilization_index)         = MeHg_Volatilization
    if (MeHg_Photoreduction_index > 0)          a(MeHg_Photoreduction_index)         = MeHg_Photoreduction
    if (MeHg_Demethylation_index > 0)           a(MeHg_Demethylation_index)          = MeHg_Demethylation
    if (MeHg_Settling_index > 0)                a(MeHg_Settling_index)               = MeHg_Settling
    !
    if (use_BedSediment) then
      if (HgII2_Methylation_index > 0)          a(HgII2_Methylation_index)           = HgII2_Methylation
      if (HgII_Deposition_index > 0)            a(HgII_Deposition_index)             = HgII_Deposition
      if (HgII_Resuspension_index > 0)          a(HgII_Resuspension_index)           = HgII_Resuspension
      if (HgII_Erosion_index > 0)               a(HgII_Erosion_index)                = HgII_Erosion
      if (HgII_Burial_index > 0)                a(HgII_Burial_index)                 = HgII_Burial
      if (HgII_Transfer_index > 0)              a(HgII_Transfer_index)               = HgII_Transfer
      if (HgII2_Transfer_index > 0)             a(HgII2_Transfer_index)              = HgII2_Transfer
      if (MeHg2_Demethylation_index > 0)        a(MeHg2_Demethylation_index)         = MeHg2_Demethylation
      if (MeHg_Deposition_index > 0)            a(MeHg_Deposition_index)             = MeHg_Deposition
      if (MeHg_Resuspension_index > 0)          a(MeHg_Resuspension_index)           = MeHg_Resuspension
      if (MeHg_Erosion_index > 0)               a(MeHg_Erosion_index)                = MeHg_Erosion
      if (MeHg_Burial_index > 0)                a(MeHg_Burial_index)                 = MeHg_Burial
      if (MeHg_Transfer_index > 0)              a(MeHg_Transfer_index)               = MeHg_Transfer
      if (MeHg2_Transfer_index > 0)             a(MeHg2_Transfer_index)              = MeHg2_Transfer 
    end if  
  end subroutine
  !===========================================================================================================================
  ! Compute derived variables
  subroutine ComputeHgDerivedVariables
    real(R8):: solid_total, solid2_total
    !
    call HgPartitions()
    !
    ! water column 
    if (use_NonEquilibrium(1)) then 
      if (t > 1.0E-10) then
          
        Hg(2) = Hgd(2)
        if (use_DOCSorbed(2))       Hg(2) = Hg(2) + Hgdoc(2)
        if (use_AlgaeSorbed(2))     Hg(2) = Hg(2) + Hgap(2)
        if (use_POMSorbed(2))       Hg(2) = Hg(2) + Hgpom(2)
        do k = 1, nGS
          if (use_SolidSorbed(2,k)) Hg(2) = Hg(2) + Hgp(2,k)
        end do
      end if
    end if
    !
    ! sediment layer
    if (use_BedSediment) then
      if (use_NonEquilibrium(2)) then
        if (t > 1.0E-10) then    
          Hg2(2) = Hgd2(2) * Por(r)
          if (use_DOCSorbed(2)) Hg2(2) = Hg2(2) + Hgdoc2(2) * Por(r)
          if (use_POMSorbed(2)) Hg2(2) = Hg2(2) + Hgpom2(2)
          do k = 1, nGS
            if (use_SolidSorbed(2,k)) Hg2(2) = Hg2(2) + Hgp2(2,k)
          end do
        end if
      else
        Hgd2(2)   = Cd2(2)   / Por(r)
        Hgdoc2(2) = Cdoc2(2) / Por(r)      
      end if
      Hgd2(3)   = Cd2(3)   / Por(r)
      Hgdoc2(3) = Cdoc2(3) / Por(r)
    end if
    !
    do i = 2, 3
      solid_total = 0.0
      Hgpt(i)     = 0.0 
      Hgpts(i)    = 0.0 
      if (use_AlgaeSorbed(i)) then
        solid_total  = solid_total + Apd
        Hgpt(i)      = Hgpt(i) + Hgap(i)
      end if
      if (use_POMSorbed(i)) then         
        solid_total  = solid_total + POM
        Hgpt(i)      = Hgpt(i) + Hgpom(i)
      end if
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) then   
          solid_total  = solid_total + solid(k)
          Hgpt(i)      = Hgpt(i) + Hgp(i,k)
        end if
      end do 
      if (solid_total > 0.0) then
        Hgpts(i) = Hgpt(i) / solid_total * 1.0E3
      else
        Hgpts(i) = 0.0
      end if 
      !
      solid2_total = 0.0
      Hgpt2(i)     = 0.0 
      Hgpts2(i)    = 0.0
      if (use_BedSediment) then
        if (use_POMSorbed(i)) then       
          solid2_total = solid2_total + POM2 
          Hgpt2(i)     = Hgpt2(i) + Hgpom2(i) 
        end if
        do k = 1, nGS
          if (use_SolidSorbed(i,k)) then 
            solid2_total = solid2_total + solid2(k)
            Hgpt2(i)     = Hgpt2(i) + Hgp2(i,k)
          end if
        end do
      end if
      if (solid2_total > 0.0) then 
        Hgpts2(i) = Hgpt2(i) / solid2_total * 1.0E3
      else
        Hgpts2(i) = 0.0
      end if
    end do
  end subroutine
  !=========================================================================================================================== 
  ! Newton-Raphson solution
  subroutine NewtonRaphson(Cd_cmp)
    real(R8) :: Cd_cmp, Cd_new
    real(R8) :: ea
    integer  :: iteration
    ! 
    ! Theory of this method
    ! xn+1 = xn - f(xn)/f(xn)'
    ! 
    ! for two exceptional conditions
    Cd_cmp = 0.0
    if (f(Cd_cmp) == 0.0) return
    !
    if (IsWaterCell) then
      Cd_cmp = Hg(i)
    else
      Cd_cmp = Hg2(i)
    end if
    if (f(Cd_cmp) == 0.0) return
    !
    ! compute under normal conditions
    ! set the initial estimation of Cd
    if ((IsWaterCell .and. use_Langmuir(i,1)) .or. (.not. IsWaterCell .and. use_Langmuir(i,2))) then
      Cd_cmp = 0.0
    end if
    !
    if ((IsWaterCell .and. use_Freundlich(i,1)) .or. (.not. IsWaterCell .and. use_Freundlich(i,2))) then
      ! The smaller Cd_cmp, the better; but non-zero!
      Cd_cmp = Cd_cmp / 10.0
      do while (f(Cd_cmp) > 0.0)
        Cd_cmp = Cd_cmp / 10.0
      end do
    end if
    !
    iteration = 0
    do 
      iteration = iteration + 1
      Cd_new    = Cd_cmp - f(Cd_cmp) / df(Cd_cmp)
      if (Cd_new .ne. 0.0) ea = abs(Cd_new - Cd_cmp) / Cd_new
      Cd_cmp    = Cd_new
      !
      if (iteration >= imax .or. Cd_cmp < 0.0) then
        write(*,*) 'Newton-Raphson method does not converge when computing dissolved HgII concentration. & 
                    Try Bisection method or reduce time step'
        stop
      end if
      !
      if (ea < res(r)) exit
    end do
  end subroutine
  !
  !===========================================================================================================================
  ! Bisection solution
  subroutine Bisection(Cd_cmp)
    real(R8) :: Cd_cmp
    real(R8) :: xlow, xup, xr
    real(R8) :: flow, fup, fr
    integer  :: iteration
    ! 
    xlow = 0.0
    flow = f(xlow)
    if (flow == 0.0) then
      Cd_cmp = xlow
      return
    end if
    !
    if (IsWaterCell) then
      xup = Hg(i)
    else
      xup = Hg2(i) 
    end if 
    fup = f(xup)
    if (fup == 0.0) then
      Cd_cmp = xup
      return
    end if
    !
    if (flow * fup > 0.0) then
      write(*,*) 'Bad dissolved HgII concentration guess & Try a smaller time step'
      stop
    end if
    !
    do iteration = 1, imax
      xr = (xlow + xup) / 2.0
      fr = f(xr)
      if (flow * fr < 0.0) then
        xup = xr
      else if (flow * fr > 0.0) then
        xlow = xr
      else
        Cd_cmp = xr
        exit
      end if
      if (abs(xup - xlow) / xr < res(r)) then
        Cd_cmp = xr
        exit
      end if
    end do
    !
    if (iteration > imax) then
      write(*,*) 'Bisection method does not converge when computing dissolved HgII concentration.'
      stop
    end if 
  end subroutine
end module

'''