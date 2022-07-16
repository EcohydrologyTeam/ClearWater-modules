'''
=======================================================================================
Contaminant Simulation Module (CSM): CSM Functions
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

def ContaminantTempCorrection(TwaterC: float, TsedC: float, nC: int, nSpecies: list):

  # Volatilization transfer coefficients
  TRK = 0.0
  Coef_tc = 0.0
  Coef0_tc = 0.0
  KL = 0.0
  KG = 0.0

  # Compute temperature in Kelvins
  TwaterK = TwaterC + 273.15
  TsedK   = TsedC + 273.15

'''

  do i = 1, nC
    do j = 1, nSpecies(i)
      id = parameter_index(i,j)

      # Non-equilibrium
      if (use_NonEquilibrium(i,1)) then
        TRK = Trade(id,r) + 273.15
        if (use_AlgaeSorbed(i) .or. use_POMSorbed(i) .or. use_AnySolidSorbed(i)) then
          Coef_tc  = exp(Eaad(id,r) * 1000.0 * (TwaterK - TRK) / (gas_constant * TwaterK * TRK))
          Coef0_tc = exp(Eade(id,r) * 1000.0 * (TwaterK - TRK) / (gas_constant * TwaterK * TRK))
        end if
        if (use_AlgaeSorbed(i)) then
          kadap_tc(id) = kadap(id,r) * Coef_tc
          kdeap_tc(id) = kdeap(id,r) * Coef0_tc
        end if
        if (use_POMSorbed(i)) then
          kadpom_tc(id) = kadpom(id,r) * Coef_tc
          kdepom_tc(id) = kdepom(id,r) * Coef0_tc
        end if
        do k = 1, nGS
          if (use_SolidSorbed(i,k)) then
            kadp_tc(id,k) = kadp(id,k,r) * Coef_tc
            kdep_tc(id,k) = kdep(id,k,r) * Coef0_tc
          end if
        end do
      end if
      #
      if (use_NonEquilibrium(i,2)) then
        TRK = Trade(id,r) + 273.15
        if (use_POMSorbed(i) .or. use_AnySolidSorbed(i)) then
          Coef_tc  = exp(Eaad(id,r) * 1000.0 * (TsedK - TRK) / (gas_constant * TsedK * TRK))
          Coef0_tc = exp(Eade(id,r) * 1000.0 * (TsedK - TRK) / (gas_constant * TsedK * TRK))
        end if
        if (use_POMSorbed(i)) then
          kadpom2_tc(id) = kadpom2(id,r) * Coef_tc
          kdepom2_tc(id) = kdepom2(id,r) * Coef0_tc
        end if
        do k = 1, nGS
          if (use_SolidSorbed(i,k)) then
            kadp2_tc(id,k) = kadp2(id,k,r) * Coef_tc
            kdep2_tc(id,k) = kdep2(id,k,r) * Coef0_tc
          end if
        end do
      end if
      # 
      # decay rate - MAF method
      if (use_nOrderDecay(i)) then
        k1d_tc(id) = Arrhenius_TempCorrection(k1d(id,r),    TwaterC)
        if (use_DOCSorbed(i))        k1doc_tc(id)  = Arrhenius_TempCorrection(k1doc(id,r),  TwaterC)
        if (use_AlgaeSorbed(i))      k1ap_tc(id)   = Arrhenius_TempCorrection(k1ap(id,r),   TwaterC)
        if (use_POMSorbed(i))        k1pom_tc(id)  = Arrhenius_TempCorrection(k1pom(id,r),  TwaterC)
        if (use_AnySolidSorbed(i))   k1p_tc(id)    = Arrhenius_TempCorrection(k1p(id,r),    TwaterC)
        if (use_BedSediment) then 
          k1d2_tc(id) = Arrhenius_TempCorrection(k1d2(id,r),   TsedC)
          if (use_DOCSorbed(i))      k1doc2_tc(id) = Arrhenius_TempCorrection(k1doc2(id,r), TsedC)
          if (use_POMSorbed(i))      k1pom2_tc(id) = Arrhenius_TempCorrection(k1pom2(id,r), TsedC)
          if (use_AnySolidSorbed(i)) k1p2_tc(id)   = Arrhenius_TempCorrection(k1p2(id,r),   TsedC)
        end if
      end if
      #
      # hydrolysis rate - AF method
      if (use_Hydrolysis(i)) then
        TRK = Trhyd(id,r) + 273.15
        #
        Coef_tc       = exp(Eahb(id,r) * 1000.0 * (TwaterK - TRK) / (gas_constant * TwaterK * TRK))
        khbd_tc(id)   = khbd(id,r)   * Coef_tc
        if (use_DOCSorbed(i)) khbdoc_tc(id) = khbdoc(id,r) * Coef_tc
        #
        Coef_tc       = exp(Eahn(id,r) * 1000.0 * (TwaterK - TRK) / (gas_constant * TwaterK * TRK))
        khnd_tc(id)   = khnd(id,r)   * Coef_tc
        if (use_DOCSorbed(i)) khndoc_tc(id) = khndoc(id,r) * Coef_tc
        #
        Coef_tc       = exp(Eaha(id,r) * 1000.0 * (TwaterK - TRK) / (gas_constant * TwaterK * TRK))
        khad_tc(id)   = khad(id,r)   * Coef_tc
        if (use_DOCSorbed(i)) khadoc_tc(id) = khadoc(id,r) * Coef_tc
        #
        if (use_BedSediment) then
          Coef_tc        = exp(Eahb(id,r) * 1000.0 * (TsedK - TRK) / (gas_constant * TsedK * TRK))
          khbd2_tc(id)   = khbd2(id,r)   * Coef_tc
          if (use_DOCSorbed(i)) khbdoc2_tc(id) = khbdoc2(id,r) * Coef_tc
          #
          Coef_tc        = exp(Eahn(id,r) * 1000.0 * (TsedK - TRK) / (gas_constant * TsedK * TRK))
          khnd2_tc(id)   = khnd2(id,r)   * Coef_tc
          if (use_DOCSorbed(i)) khndoc2_tc(id) = khndoc2(id,r) * Coef_tc
          #
          Coef_tc        = exp(Eaha(id,r) * 1000.0 * (TsedK - TRK) / (gas_constant * TsedK * TRK))
          khad2_tc(id)   = khad2(id,r)   * Coef_tc
          if (use_DOCSorbed(i)) khadoc2_tc(id) = khadoc2(id,r) * Coef_tc
        end if
      end if
      #
      # photolysis rate - MAF method
      if (use_Photolysis(i)) then
        kphtd_tc(id)   = Arrhenius_TempCorrection(kphtd(id,r),   TwaterC)
        if (use_DOCSorbed(i)) kphtdoc_tc(id) = Arrhenius_TempCorrection(kphtdoc(id,r), TwaterC) 
      end if
      #
      # volatilization velocity - MAF method
      if (use_Volatilization(i) .and. j == 1) then
        if (vv_option(i,r) == 2) then
          KL = ka * (32.0 / MW(id,r))**0.25
          KG = 168.0 * wind_speed * (18.0 / MW(id,r))**0.25
          if (KG < 100.0) KG = 100.0
          vv(i,r)%rc20 = 1.0 / (1.0 / KL + gas_constant * TwaterK / KH(i,r) / KG)
        end if
        #
        vv_tc(i) = Arrhenius_TempCorrection(vv(i,r), TwaterC)
      end if
      #
      # 2nd reaction
      if (use_2ndReaction(i)) then
        TRK = Trer(id,r) + 273.15
        #
        Coef_tc       = exp(Eaer(id,r) * 1000.0 * (TwaterK - TRK) / (gas_constant * TwaterK * TRK))
        kerd_tc(id)   = kerd(id,r)   * Coef_tc
        if (use_DOCSorbed(i))      kerdoc_tc(id) = kerdoc(id,r) * Coef_tc
        if (use_AlgaeSorbed(i))    kerap_tc(id)  = kerap(id,r)  * Coef_tc
        if (use_POMSorbed(i))      kerpom_tc(id) = kerpom(id,r) * Coef_tc
        if (use_AnySolidSorbed(i)) kerp_tc(id)   = kerp(id,r)   * Coef_tc 
        if (use_BedSediment) then
          Coef_tc        = exp(Eaer(id,r) * 1000.0 * (TsedK - TRK) / (gas_constant * TsedK * TRK))
          kerd2_tc(id)   = kerd2(id,r)   * Coef_tc
          if (use_DOCSorbed(i))      kerdoc2_tc(id) = kerdoc2(id,r) * Coef_tc
          if (use_POMSorbed(i))      kerpom2_tc(id) = kerpom2(id,r) * Coef_tc
          if (use_AnySolidSorbed(i)) kerp2_tc(id)   = kerp2(id,r)   * Coef_tc 
        end if
      end if
    end do
  end do
end subroutine


#=========================================================================================================================== 
# Compute equilibrium partition coefficient if necessary
subroutine ContaminantPartitionCoef()
  #
  if (abs(t) < 1.0E-10) then
    do i = 1, nC
      do j = 1, nSpecies(i)
        id = parameter_index(i,j)
        #
        # compute partition coefficients from Kow 
        if (kd_option(id,r) == 2) then
          if (use_DOCSorbed(i)) then
            Kdoc(id,r)  = adoc(id,r) * Kow(id,r)
            if (use_BedSediment) Kdoc2(id,r) = adoc2(id,r) * Kow(id,r)
          end if
          #
          if (use_Equilibrium(i,1)) then
            if (use_AlgaeSorbed(i)) Kap(id,r)  = aap(id,r) * Kow(id,r)
            if (use_POMSorbed(i))   Kpom(id,r) = apom(id,r) * Kow(id,r)
            do k = 1, nGS
              if (use_SolidSorbed(i,k)) Kp(id,k,r) = ap(id,k,r) * Kow(id,r)
            end do
          end if
          #
          if (use_Equilibrium(i,2)) then
            if (use_POMSorbed(i))       Kpom2(id,r) = apom2(id,r) * Kow(id,r)
            do k = 1, nGS
              if (use_SolidSorbed(i,k)) Kp2(id,k,r) = ap2(id,k,r) * Kow(id,r)
            end do
          end if
        end if
        #
      end do
    end do
  end if
end subroutine

#=========================================================================================================================== 
# Compute ionization 
subroutine ContaminantIonization()
  real(R8) :: Coef_tc
  real(R8) :: D
  #
  # compute CHH, COH
  if (use_Ionization_One .or. use_Hydrolysis_one) then
    CHH = 10.0**(-pH)
    COH = 10.0**(pH - 14.0)
  end if    
  #
  do i = 1, nC
    if (nSpecies(i) == 5) then
      #
      # temperature correction for ionization constants - AF method
      Coef_tc = 1000.0 * (TwaterC - Trion(i,r)) / (gas_constant * (TwaterC + 273.15) * (Trion(i,r) + 273.15))
      # change unit of dH from kJ/mol to J/mol
      Ka1_tc(i) = Ka1(i,r) * exp(dHa1(i,r) * Coef_tc)
      Ka2_tc(i) = Ka2(i,r) * exp(dHa2(i,r) * Coef_tc)
      Kb1_tc(i) = Kb1(i,r) * exp(dHb1(i,r) * Coef_tc)
      Kb2_tc(i) = Kb2(i,r) * exp(dHb2(i,r) * Coef_tc)
      if (use_BedSediment) then
        Coef_tc = 1000.0 * (TsedC - Trion(i,r)) / (gas_constant * (TsedC + 273.15) * (Trion(i,r) + 273.15))
        #
        Ka12_tc(i) = Ka1(i,r) * exp(dHa1(i,r) * Coef_tc)
        Ka22_tc(i) = Ka2(i,r) * exp(dHa2(i,r) * Coef_tc)
        Kb12_tc(i) = Kb1(i,r) * exp(dHb1(i,r) * Coef_tc)
        Kb22_tc(i) = Kb2(i,r) * exp(dHb2(i,r) * Coef_tc)
      end if
      #
      # compute ionic fraction for each species
      D = 1.0 + Kb1_tc(i) / COH + Kb1_tc(i) * Kb2_tc(i) / COH / COH   &  
              + Ka1_tc(i) / CHH + Ka1_tc(i) * Ka2_tc(i) / CHH / CHH
      fion(i,1) = 1.0 / D
      fion(i,2) = Kb1_tc(i) / COH / D
      fion(i,3) = Kb1_tc(i) * Kb2_tc(i) / COH / COH / D
      fion(i,4) = Ka1_tc(i) / CHH / D
      fion(i,5) = Ka1_tc(i) * Ka2_tc(i) / CHH / CHH / D
      if (use_BedSediment) then
        D = 1.0 + Kb12_tc(i) / COH + Kb12_tc(i) * Kb22_tc(i) / COH / COH   &
                + Ka12_tc(i) / CHH + Ka12_tc(i) * Ka22_tc(i) / CHH / CHH
        fion2(i,1) = 1.0 / D
        fion2(i,2) = Kb12_tc(i) / COH / D
        fion2(i,3) = Kb12_tc(i) * Kb22_tc(i) / COH / COH / D
        fion2(i,4) = Ka12_tc(i) / CHH / D
        fion2(i,5) = Ka12_tc(i) * Ka22_tc(i) / CHH / CHH / D
      end if
    else if (nSpecies(i) == 1) then
      fion(i,1)  = 1.0
      if (use_BedSediment) fion2(i,1) = 1.0
    end if
  end do
end subroutine

#=========================================================================================================================== 
#------------------------------------------------------------------------------------------------------------------
#$$$$$$$$        Cd2 and Cdoc2 represent dissolved concentration in pore water.
#$$$$$$$$        Cd2_Species and Cdoc2_Species represent dissolved concentration on bed sediment total volume.
#------------------------------------------------------------------------------------------------------------------
# compute Cd_Species, Cdoc_Species, Cap_Species, Cpom_Species, Cp_Species form linear equilibrium partition
subroutine EquilibriumPartitionConc()
  real(R8) :: Rd
  #
  if (IsWaterCell) then
    do j = 1, nSpecies(i)
      id = parameter_index(i,j)
      #
      Rd = 1.0
      if (use_DOCSorbed(i))        Rd = Rd + Kdoc(id,r) * DOC / 1.0E6
      if (use_AlgaeSorbed(i))      Rd = Rd + Kap(id,r)  * Apd / 1.0E6
      if (use_POMSorbed(i))        Rd = Rd + Kpom(id,r) * POM / 1.0E6
      do k = 1, nGS
        if (use_SolidSorbed(i,k))  Rd = Rd + Kp(id,k,r) * Solid(k) / 1.0E6
      end do
      #
      Cd_Species(id) = C(i) * fion(i,j) / Rd
      if (use_AlgaeSorbed(i)) then
        Cap_Species(id)  = Kap(id,r) * Apd / 1.0E6 / Rd * C(i) * fion(i,j)
      else
        Cap_Species(id)  = 0.0
      end if
      if (use_DOCSorbed(i)) then 
        Cdoc_Species(id) = Kdoc(id,r) * DOC / 1.0E6 / Rd * C(i) * fion(i,j)
      else
        Cdoc_Species(id) = 0.0
      end if
      if (use_POMSorbed(i)) then
        Cpom_Species(id) = Kpom(id,r) * POM / 1.0E6 / Rd * C(i) * fion(i,j)
      else
        Cpom_Species(id) = 0.0
      end if
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) then
          Cp_Species(id,k) = Kp(id,k,r) * Solid(k) / 1.0E6 / Rd * C(i) * fion(i,j)
        else
          Cp_Species(id,k) = 0.0
        end if
      end do
    end do
  else
    #
    # sediment layer
    do j = 1, nSpecies(i)
      id = parameter_index(i,j)
      #
      Rd  = Por(r)
      if (use_DOCSorbed(i))        Rd = Rd + Kdoc2(id,r) * DOC2 / 1.0E6 * Por(r)
      if (use_POMSorbed(i))        Rd = Rd + Kpom2(id,r) * POM2 / 1.0E6
      do k = 1, nGS
        if (use_SolidSorbed(i,k))  Rd = Rd + Kp2(id,k,r) * Solid2(k) / 1.0E6
      end do
      #
      Cd2_Species(id) = Por(r) / Rd * C2(i) * fion2(i,j)
      if (use_DOCSorbed(i)) then
        Cdoc2_Species(id) = Kdoc2(id,r) * DOC2 / 1.0E6 * Por(r) / Rd * C2(i) * fion2(i,j)
      else
        Cdoc2_Species(id) = 0.0
      end if
      if (use_POMSorbed(i)) then
        Cpom2_Species(id) = Kpom2(id,r) * POM2 / 1.0E6 / Rd * C2(i) * fion2(i,j)
      else
        Cpom2_Species(id) = 0.0
      end if
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) then
          Cp2_Species(id,k) = Kp2(id,k,r) * Solid2(k) / 1.0E6 / Rd * C2(i) * fion2(i,j)
        else
          Cp2_Species(id,k) = 0.0
        end if
      end do
    end do
  end if
end subroutine
#
#=========================================================================================================================== 
# compute Cd_Species, Cdoc_Species, Cap_Species, Cpom_Species, Cp_Species from Langmuire equilibrium partition
subroutine LangmuirPartitionConc()
  #
  if (IsWaterCell) then
    do j = 1, nSpecies(i)
      id = parameter_index(i,j)
      #
      # compute Cd_Species using Newton-Raphson or Bisection
      if (Cd_solution_option(i,r) == 1) then
        call NewtonRaphson(Cd_Species(id))
      else if (Cd_solution_option(i,r) == 2) then
        call Bisection(Cd_Species(id))
      end if
      #
      # compute Cdoc_Species, Cap_Species, Cpom_Species and Cp_Species from Cd_Species
      if (use_DOCSorbed(i)) then
        Cdoc_Species(id)  = DOC / 1.0E6 * Kdoc(id,r) * Cd_Species(id)
      else
        Cdoc_Species(id)  = 0.0
      end if
      if (use_AlgaeSorbed(i)) then
        Cap_Species(id) = Apd / 1.0E6 * qcap(id,r) * Klap(id,r) * Cd_Species(id) / (1.0 + Klap(id,r) * Cd_Species(id) / 1.0E3)
      else
        Cap_Species(id) = 0.0
      end if
      if (use_POMSorbed(i)) then
        Cpom_Species(id)  = POM / 1.0E6 * qcpom(id,r)  * Klpom(id,r)  * Cd_Species(id) / (1.0 + Klpom(id,r)  * Cd_Species(id)  / 1.0E3)
      else
        Cpom_Species(id)  = 0.0
      end if
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) then
          Cp_Species(id,k)  = Solid(k)  / 1.0E6 * qcp(id,k,r)  * Klp(id,k,r)  * Cd_Species(id) / (1.0 + Klp(id,k,r)  * Cd_Species(id)  / 1.0E3) 
        else
          Cp_Species(id,k)  = 0.0
        end if
      end do
    end do
  else
    #
    # sediment layer
    do j = 1, nSpecies(i)
      id = parameter_index(i,j)
      #
      if (Cd_solution_option(i,r) == 1) then
        call NewtonRaphson(Cd2_Species(id))
      else if (Cd_solution_option(i,r) == 2) then
        call Bisection(Cd2_Species(id))
      end if
      #
      if (use_DOCSorbed(i)) then
        Cdoc2_Species(id) = DOC2 / 1.0E6 * Kdoc2(id,r) * Cd2_Species(id)
      else
        Cdoc2_Species(id) = 0.0 
      end if
      if (use_POMSorbed(i)) then
        Cpom2_Species(id) = POM2 / 1.0E6 * qcpom2(id,r) * Klpom2(id,r) * Cd2_Species(id) / Por(r) / (1.0 + Klpom2(id,r) * Cd2_Species(id) / 1.0E3  / Por(r))
      else
        Cpom2_Species(id) = 0.0
      end if
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) then
          Cp2_Species(id,k) = Solid2(k) / 1.0E6 * qcp2(id,k,r) * Klp2(id,k,r) * Cd2_Species(id) / Por(r) / (1.0 + Klp2(id,k,r) * Cd2_Species(id) / 1.0E3 / Por(r))  
        else
          Cp2_Species(id,k) = 0.0
        end if
      end do
    end do
  end if
  #
end subroutine

#=========================================================================================================================== 
# compute Cd_Species, Cdoc_Species, Cap_Species, Cpom_Species, Cp_Species from Freundlich equilibrium partition
subroutine FreundlichPartitionConc()
  #
  if (IsWaterCell) then
    do j = 1, nSpecies(i)
      id = parameter_index(i,j)
      #
      # compute Cd_Species using Newton-Raphson or Bisection
      if (Cd_solution_option(i,r) == 1) then
        call NewtonRaphson(Cd_Species(id))
      else if (Cd_solution_option(i,r) == 2) then
        call Bisection(Cd_Species(id))
      end if
      #
      # compute Cdoc_Species, Cap_Species, Cpom_Species and Cp_Species from Cd_Species
      if (use_DOCSorbed(i)) then
        Cdoc_Species(id)  = DOC / 1.0E6 * Kdoc(id,r)  * Cd_Species(id)
      else
        Cdoc_Species(id)  = 0.0
      end if
      if (use_AlgaeSorbed(i)) then
        Cap_Species(id) = Apd / 1.0E3 * Kfap(id,r) * Cd_Species(id)**bap(id,r)
      else
        Cap_Species(id) = 0.0
      end if
      if (use_POMSorbed(i)) then
        Cpom_Species(id)  = POM  / 1.0E3 * Kfpom(id,r) * Cd_Species(id)**bpom(id,r)
      else
        Cpom_Species(id)  = 0.0
      end if
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) then
          Cp_Species(id,k)  = Solid(k)  / 1.0E3 * Kfp(id,k,r) * Cd_Species(id)**bp(id,k,r)
        else
          Cp_Species(id,k)  = 0.0
        end if
      end do
    end do
  else
    #
    # sediment layer
    do j = 1, nSpecies(i)
      id = parameter_index(i,j)
      #
      if (Cd_solution_option(i,r) == 1) then
        call NewtonRaphson(Cd2_Species(id))
      else if (Cd_solution_option(i,r) == 2) then
        call Bisection(Cd2_Species(id))
      end if
      #
      # compute Cdoc_Species, Cap_Species, Cpom_Species and Cp_Species from Cd_Species
      if (use_DOCSorbed(i)) then
        Cdoc2_Species(id) = DOC2 / 1.0E6 * Kdoc2(id,r) * Cd2_Species(id)
      else
        Cdoc2_Species(id) = 0.0 
      end if
      if (use_POMSorbed(i)) then
        Cpom2_Species(id) = POM2 / 1.0E3 * Kfpom2(id,r) * (Cd2_Species(id) / Por(r))**bpom2(id,r)
      else
        Cpom2_Species(id) = 0.0
      end if
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) then
          Cp2_Species(id,k) = Solid2(k) / 1.0E3 * Kfp2(id,k,r) * (Cd2_Species(id) / Por(r))**bp2(id,k,r) 
        else
          Cp2_Species(id,k) = 0.0
        end if
      end do
    end do
  end if
end subroutine

#=========================================================================================================================== 
# compute Cd_Species, Cdoc_Species, Cap_Species, Cpom_Species, Cp_Species from non-equilibrium partition
# only conpute non-equilibrium conc at t = 0.
subroutine NonEquilibriumPartitionConc()
  real(R8) :: Rd
  #
  # compute concentration of different phases
  # ionization is always turned off if non-equilibrium is modeled
  id = parameter_index(i,1)
  #
  # set initial values of Cd, Cdoc, Cap, Cpom, Cp if t = 0.0
  # water column
  if (IsWaterCell) then
    if (abs(t) < 1.0E-10) then
      Rd = 1.0
      if (use_DOCSorbed(i))        Rd = Rd + Kdoc(id,r) * DOC / 1.0E6
      if (use_AlgaeSorbed(i))      Rd = Rd + kadap(id,r)  / kdeap(id,r)  * qcap(id,r)  * 1.0E3 * Apd / 1.0E6
      if (use_POMSorbed(i))        Rd = Rd + kadpom(id,r) / kdepom(id,r) * qcpom(id,r) * 1.0E3 * POM / 1.0E6
      do k = 1, nGS
        if (use_SolidSorbed(i,k))  Rd = Rd + kadp(id,k,r) / kdep(id,k,r) * qcp(id,k,r) * 1.0E3 * Solid(k) / 1.0E6
      end do
      #
      Cd(i) = C(i) / Rd
      if (use_DOCSorbed(i)) then 
        Cdoc(i) = Kdoc(id,r) * DOC / 1.0E6 / Rd * C(i)
      else
        Cdoc(i) = 0.0
      end if
      if (use_AlgaeSorbed(i)) then
        Cap(i)  = kadap(id,r) / kdeap(id,r) * qcap(id,r) * 1.0E3 * Apd / 1.0E6 / Rd * C(i)
      else
        Cap(i)  = 0.0
      end if
      if (use_POMSorbed(i)) then
        Cpom(i) = kadpom(id,r) / kdepom(id,r) * qcpom(id,r) * 1.0E3 * POM / 1.0E6 / Rd * C(i)
      else
        Cpom(i) = 0.0
      end if
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) then
          Cp(i,k) = kadp(id,k,r) / kdep(id,k,r) * qcp(id,k,r) * 1.0E3 * Solid(k) / 1.0E6 / Rd * C(i)
        else
          Cp(i,k) = 0.0
        end if
      end do
    end if
    #
    Cd_Species(id) = Cd(i)
    if (use_DOCSorbed(i))         Cdoc_Species(id) = Cdoc(i)
    if (use_AlgaeSorbed(i))       Cap_Species(id)  = Cap(i)
    if (use_POMSorbed(i))         Cpom_Species(id) = Cpom(i)
    do k = 1, nGS
      if(use_SolidSorbed(i,k))    Cp_Species(id,k) = Cp(i,k)
    end do
  else
    #
    # sediment layer
    if (abs(t) < 1.0E-10) then
      Rd  = Por(r)
      if (use_DOCSorbed(i))        Rd = Rd + Kdoc2(id,r) * DOC2 / 1.0E6 * Por(r)
      if (use_POMSorbed(i))        Rd = Rd + kadpom2(id,r) / kdepom2(id,r) * qcpom2(id,r) * 1.0E3 * POM2 / 1.0E6
      do k = 1, nGS
        if (use_SolidSorbed(i,k))  Rd = Rd + kadp2(id,k,r) / kdep2(id,k,r) * qcp2(id,k,r) * 1.0E3 * Solid2(k) / 1.0E6
      end do
      #
      Cd2(i) = 1.0 / Rd * C2(i)
      if (use_DOCSorbed(i)) then
        Cdoc2(i) = Kdoc2(id,r) * DOC2 / 1.0E6 / Rd * C2(i)
      else
        Cdoc2(i) = 0.0
      end if
      if (use_POMSorbed(i)) then
        Cpom2(i) = kadpom2(id,r) / kdepom2(id,r) * qcpom2(id,r) * 1.0E3 * POM2 / 1.0E6 / Rd * C2(i)
      else
        Cpom2(i) = 0.0
      end if
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) then
          Cp2(i,k) = kadp2(id,k,r) / kdep2(id,k,r) * qcp2(id,k,r) * 1.0E3 * Solid2(k) / 1.0E6 / Rd * C2(i)
        else
          Cp2(i,k) = 0.0
        end if
      end do
    end if
    #
    Cd2_Species(id) = Cd2(i) * Por(r)
    if (use_DOCSorbed(i))       Cdoc2_Species(id) = Cdoc2(i) * Por(r)
    if (use_POMSorbed(i))       Cpom2_Species(id) = Cpom2(i)
    do k = 1, nGS
      if (use_SolidSorbed(i,k)) Cp2_Species(id,k) = Cp2(i,k)
    end do
  end if
end subroutine

#=========================================================================================================================== 
# the function of f(Cd)
double precision function f(x)
  real(R8), intent(in) :: x
  #
  f = x
  if (IsWaterCell) then
    f = f - C(i) * fion(i,j)
    if (use_DOCSorbed(i))         f = f + DOC      / 1.0E6 * Kdoc(id,r) * x 
    #
    if (use_Langmuir(i,1)) then
      if (use_AlgaeSorbed(i))     f = f + Apd      / 1.0E3 * qcap(id,r)  * (Klap(id,r)  * x / 1.0E3) / (1.0 + Klap(id,r)  * x / 1.0E3)
      if (use_POMSorbed(i))       f = f + POM      / 1.0E3 * qcpom(id,r) * (Klpom(id,r) * x / 1.0E3) / (1.0 + Klpom(id,r) * x / 1.0E3)
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) f = f + Solid(k) / 1.0E3 * qcp(id,k,r) * (Klp(id,k,r) * x / 1.0E3) / (1.0 + Klp(id,k,r) * x / 1.0E3)
      end do
    else if (use_Freundlich(i,1)) then
      if (use_AlgaeSorbed(i))     f = f + Apd      / 1.0E3 * Kfap(id,r)  * x**bap(id,r)
      if (use_POMSorbed(i))       f = f + POM      / 1.0E3 * Kfpom(id,r) * x**bpom(id,r)
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) f = f + Solid(k) / 1.0E3 * Kfp(id,k,r) * x**bp(id,k,r) 
      end do
    end if
  else
    #
    # sediment layer
    f = f - C2(i) * fion2(i,j)
    if (use_DOCSorbed(i))         f = f + DOC2      / 1.0E6 * Kdoc2(id,r) * x
    if (use_Langmuir(i,2)) then
      if (use_POMSorbed(i))       f = f + POM2      / 1.0E3 * qcpom2(id,r) * (Klpom2(id,r) * x / Por(r) / 1.0E3) / (1.0 + Klpom2(id,r) * x / Por(r) / 1.0E3)
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) f = f + Solid2(k) / 1.0E3 * qcp2(id,k,r) * (Klp2(id,k,r) * x / Por(r) / 1.0E3) / (1.0 + Klp2(id,k,r) * x / Por(r) / 1.0E3)
      end do
    else if (use_Freundlich(i,2)) then
      if (use_POMSorbed(i))       f = f + POM2      / 1.0E3 * Kfpom2(id,r) * (x / Por(r))**bpom2(id,r)
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) f = f + Solid2(k) / 1.0E3 * Kfp2(id,k,r) * (x / Por(r))**bp2(id,k,r) 
      end do
    end if
  end if
end function

#=========================================================================================================================== 
# the function of df(Cd)
double precision function df(x)
  real(R8), intent(in) :: x
  #
  df = 1
  if (IsWaterCell) then
    if (use_DOCSorbed(i))         df = df + DOC      / 1.0E6 * Kdoc(id,r)
    if (use_Langmuir(i,1)) then
      if (use_AlgaeSorbed(i))     df = df + Apd      / 1.0E3 * qcap(id,r)  * (Klap(id,r)  / 1.0E3) / (1.0 + Klap(id,r)  * x / 1.0E3)**2.0
      if (use_POMSorbed(i))       df = df + POM      / 1.0E3 * qcpom(id,r) * (Klpom(id,r) / 1.0E3) / (1.0 + Klpom(id,r) * x / 1.0E3)**2.0
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) df = df + Solid(k) / 1.0E3 * qcp(id,k,r) * (Klp(id,k,r) / 1.0E3) / (1.0 + Klp(id,k,r) * x / 1.0E3)**2.0
      end do
    else if (use_Freundlich(i,1)) then
      if (use_AlgaeSorbed(i))     df = df + Apd      / 1.0E3 * Kfap(id,r)  * bap(id,r)  * x**(bap(id,r)  - 1.0)
      if (use_POMSorbed(i))       df = df + POM      / 1.0E3 * Kfpom(id,r) * bpom(id,r) * x**(bpom(id,r) - 1.0)
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) df = df + Solid(k) / 1.0E3 * Kfp(id,k,r) * bp(id,k,r) * x**(bp(id,k,r) - 1.0)
      end do
    end if
  else
    #
    # sediment layer
    if (use_DOCSorbed(i))         df = df + DOC2      / 1.0E6 * Kdoc2(id,r)
    if (use_Langmuir(i,2)) then
      if (use_POMSorbed(i))       df = df + POM2      / 1.0E3 * qcpom2(id,r) * (Klpom2(id,r) / Por(r) / 1.0E3) / (1.0 + Klpom2(id,r) * x / Por(r) / 1.0E3)**2.0
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) df = df + Solid2(k) / 1.0E3 * qcp2(id,k,r) * (Klp2(id,k,r) / Por(r) / 1.0E3) / (1.0 + Klp2(id,k,r) * x / Por(r) / 1.0E3)**2.0
      end do
    else if (use_Freundlich(i,2)) then
      if (use_POMSorbed(i))       df = df + POM2      / 1.0E3 * Kfpom2(id,r) * bpom2(id,r) * (x / Por(r))**(bpom2(id,r) - 1.0) / Por(r)
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) df = df + Solid2(k) / 1.0E3 * Kfp2(id,k,r) * bp2(id,k,r) * (x / Por(r))**(bp2(id,k,r) - 1.0) / Por(r)
      end do
    end if
  end if    
end function

#=========================================================================================================================== 
# compute Cd_Species, Cdoc_Species, Cap_Species, Cpom_Species, Cp_Species
#         Cd2_Species, Cdoc2_Species, Cpom2_Species, Cp2_Species
subroutine ContaminantPartitions()
  integer :: l
  logical :: ParticleSorbed
  #
  call ContaminantPartitionCoef()
  call ContaminantIonization()
  #
  do i = 1, nC
    do l = 1, 2
      ParticleSorbed = .false.
      if (l == 1) then
        IsWaterCell = .true.
        if (use_AlgaeSorbed(i) .or. use_POMSorbed(i) .or. use_AnySolidSorbed(i)) ParticleSorbed = .true.
      else
        IsWaterCell = .false.
        if (use_POMSorbed(i) .or. use_AnySolidSorbed(i))                         ParticleSorbed = .true.
      end if
      #
      if (use_Equilibrium(i,l)) then
        call EquilibriumPartitionConc()
      else if (use_Langmuir(i,l)) then
        if (ParticleSorbed) then
          call LangmuirPartitionConc()
        else
          # if only DOCSorbed is selected, use linear Equilibrium 
          call EquilibriumPartitionConc()
        end if
      else if (use_Freundlich(i,l)) then
        if (ParticleSorbed) then
          call FreundlichPartitionConc()
        else
          # if only DOCSorbed is selected, use linear Equilibrium model
          call EquilibriumPartitionConc()
        end if
      else if (use_NonEquilibrium(i,l)) then
        call NonEquilibriumPartitionConc()
      end if
    end do
  end do
  #
end subroutine

#=========================================================================================================================== 
# Compute pathways 
subroutine ContaminantPathways()
  real(R8) :: lambdamax
  real(R8) :: kviscosity
  real(R8) :: C_Acid_Hydrolysis, C_Neutral_Hydrolysis, C_Base_Hydrolysis
  real(R8) :: C2_Acid_Hydrolysis, C2_Neutral_Hydrolysis, C2_Base_Hydrolysis
  real(R8) :: Kd2_avg, Tsolid2
  #
  #----------------------------------
  # prepare for pathway computation #
  #----------------------------------
  # get lambdamax
  if (use_Photolysis_one) then 
    lambdamax = alpha(r) * lambda
    Icpht     = 1.33 * q_solar * (1.0 - exp(- lambdamax * depth)) / (lambdamax * depth) * (1.0 - 0.56 * cloudiness)
  end if
  #
  do i = 1, nC
    #
    #-----------------
    # air deposition #
    #-----------------
    C_Air_Deposition(i) = L0(i) * surface_area / volume
    #
    do j = 1, nSpecies(i)
      id = parameter_index(i,j)
      #
      #---------------------------------------
      # settling and resuspension and burial #
      #---------------------------------------
      if (use_POMSorbed(i)) then
        C_C2_Settling(id) = Cpom_Species(id)  * vsom(r) / depth
        if (use_BedSediment)  then
          C2_Settling(id) = Cpom_Species(id)  * vsom(r) / h2(r)
          C2_Burial(id)   = Cpom2_Species(id) * vb / h2(r)
        end if
      else
        C_C2_Settling(id) = 0.0
        if (use_BedSediment) then
          C2_Settling(id) = 0.0
          C2_Burial(id)   = 0.0
        end if
      end if
      #
      if (use_AlgaeSorbed(i)) then
        C_C2_Settling(id) = C_C2_Settling(id) + Cap_Species(id) * vsap(r) / depth
        if (use_BedSediment)  C2_Settling(id) = C2_Settling(id) + Cap_Species(id) * vsap(r) / h2(r)
      end if
      #
      # solids
      if (use_BedSediment) then
        C_C2_Resuspension(id) = 0.0
        C2_C_Resuspension(id) = 0.0
      end if
      #
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) then
          C_C2_Settling(id) = C_C2_Settling(id) + Cp_Species(id,k) * vsp(k) / depth
          if (use_BedSediment) then
            C2_Settling(id)       = C2_Settling(id)       + Cp_Species(id,k)  * vsp(k) / h2(r)
            C_C2_Resuspension(id) = C_C2_Resuspension(id) + Cp2_Species(id,k) * vrp(k) / depth
            C2_C_Resuspension(id) = C2_C_Resuspension(id) + Cp2_Species(id,k) * vrp(k) / h2(r) 
            C2_Burial(id)         = C2_Burial(id)         + Cp2_Species(id,k) * vb / h2(r)
          end if
        end if
      end do
      #
      #----------------
      # mass transfer #
      #----------------
      # get vm
      if (use_BedSediment) then
        #
        if (vm_option(id,r) == 2) then
#           formula given by Thibodeaux et al. (2001): bioturbation-driven transport, sorb to sediment
          if (use_Equilibrium(i,2)) then
            Kd2_avg = 0.0
            Tsolid2 = 0.0
            if (use_POMSorbed(i)) then
              Kd2_avg = Kd2_avg + Kpom2(id,r) * POM2
              Tsolid2 = Tsolid2 + POM2
            end if
            do k = 1, nGS
              if (use_SolidSorbed(i,k)) then
                Kd2_avg = Kd2_avg + Kp2(id,k,r) * Solid2(k)
                Tsolid2 = Tsolid2 + Solid2(k)
              end if
            end do
            if (Tsolid2 > 0.0) Kd2_avg  = Kd2_avg / Tsolid2
            if (Kd2_avg > 0.0) vm(id,r) = 0.01 / (1.0 / beta(r) + z2(r) / (Db(r) * ((1.0 - Por(r)) * ps(r)) * Kd2_avg))
          end if
        else if (vm_option(id,r) == 3) then
          # formula given by Boyer et al.(1994)
          vm(id,r)  = Por(r)**3.0 * Dm(id,r) / 0.01
          # z': characteristic length over which gradient exists at sediment-water interface, 0.01m. 
        else if (vm_option(id,r) == 4) then
          # formula given by Di torio et al.(1981)
          vm(id,r)  = 0.19 * Por(r) / MW(id,r)**(2.0 / 3.0)
        else if (vm_option(id,r) == 5) then
          # formula given by Schink and Guinasso (1977)
          # TwaterC is used here, TsedC may be more accurate. 
          kviscosity = 1.79E-6 / (1.0 + 0.03368 * TwaterC + 0.000221 * TwaterC**2.0)
          vm(id,r)  = shear_velocity * (Dm(id,r) / 86400.0 / kviscosity)**(2.0 / 3.0) / 24.0 * 86400.0
        end if
        #
        C_C2_Transfer(id) = vm(id,r) * (Cd2_Species(id) / Por(r) - Cd_Species(id)) / depth
        C2_C_Transfer(id) = vm(id,r) * (Cd_Species(id) - Cd2_Species(id) / Por(r)) / h2(r)
        if (use_DOCSorbed(i)) then
          C_C2_Transfer(id) = C_C2_Transfer(id) + vm(id,r) * (Cdoc2_Species(id) / Por(r) - Cdoc_Species(id)) / depth
          C2_C_Transfer(id) = C2_C_Transfer(id) - vm(id,r) * (Cdoc2_Species(id) / Por(r) - Cdoc_Species(id)) / h2(r)
        end if
      end if
      #
      #----------------
      # decay         #
      #----------------
      if (use_nOrderDecay(i)) then
        C_Decay(id)  = Cd_Species(id)**nOrder(id,r)  * k1d_tc(id)
        if (use_BedSediment)      C2_Decay(id) = Cd2_Species(id)**nOrder(id,r) * k1d2_tc(id)
        if (use_DOCSorbed(i)) then
          C_Decay(id)  = C_Decay(id)  + Cdoc_Species(id)**nOrder(id,r)  * k1doc_tc(id) 
          if (use_BedSediment)    C2_Decay(id) = C2_Decay(id) + Cdoc2_Species(id)**nOrder(id,r) * k1doc2_tc(id)
        end if
        if (use_AlgaeSorbed(i))   C_Decay(id)  = C_Decay(id)  + Cap_Species(id)**nOrder(id,r)   * k1ap_tc(id)
        if (use_POMSorbed(i)) then
          C_Decay(id)  = C_Decay(id)  + Cpom_Species(id)**nOrder(id,r)  * k1pom_tc(id)
          if (use_BedSediment)    C2_Decay(id) = C2_Decay(id) + Cpom2_Species(id)**nOrder(id,r) * k1pom2_tc(id)
        end if
        #
        do k = 1, nGS
          if (use_SolidSorbed(i,k)) then  
            C_Decay(id)  = C_Decay(id)  + Cp_Species(id,k)**nOrder(id,r)  * k1p_tc(id)
            if (use_BedSediment)  C2_Decay(id) = C2_Decay(id) + Cp2_Species(id,k)**nOrder(id,r) * k1p2_tc(id)
          end if
        end do
      else
        C_Decay(id)  = 0.0
        if (use_BedSediment)      C2_Decay(id) = 0.0
      end if
      #
      #-------------
      # hydrolysis #
      #-------------
      if (use_Hydrolysis(i)) then
        C_Acid_Hydrolysis    = Cd_Species(id) * khad_tc(id) * CHH
        C_Neutral_Hydrolysis = Cd_Species(id) * khnd_tc(id)
        C_Base_Hydrolysis    = Cd_Species(id) * khbd_tc(id) * COH
        if (use_DOCSorbed(i)) then 
          C_Acid_Hydrolysis    = C_Acid_Hydrolysis    + Cdoc_Species(id) * khadoc_tc(id) * CHH
          C_Neutral_Hydrolysis = C_Neutral_Hydrolysis + Cdoc_Species(id) * khndoc_tc(id)
          C_Base_Hydrolysis    = C_Base_Hydrolysis    + Cdoc_Species(id) * khbdoc_tc(id) * COH
        end if
        #
        if (use_BedSediment) then
          C2_Acid_Hydrolysis    = Cd2_Species(id) * khad2_tc(id) * CHH
          C2_Neutral_Hydrolysis = Cd2_Species(id) * khnd2_tc(id)
          C2_Base_Hydrolysis    = Cd2_Species(id) * khbd2_tc(id) * COH
          if (use_DOCSorbed(i)) then
            C2_Acid_Hydrolysis    = C2_Acid_Hydrolysis    + Cdoc2_Species(id) * khadoc2_tc(id) * CHH
            C2_Neutral_Hydrolysis = C2_Neutral_Hydrolysis + Cdoc2_Species(id) * khndoc2_tc(id)
            C2_Base_Hydrolysis    = C2_Base_Hydrolysis    + Cdoc2_Species(id) * khbdoc2_tc(id) * COH
          end if
        end if
        C_Hydrolysis(id)  = C_Acid_Hydrolysis  + C_Neutral_Hydrolysis  + C_Base_Hydrolysis
        if (use_BedSediment) C2_Hydrolysis(id) = C2_Acid_Hydrolysis + C2_Neutral_Hydrolysis + C2_Base_Hydrolysis
      else
        C_Hydrolysis(id)  = 0.0
        if (use_BedSediment) C2_Hydrolysis(id) = 0.0
      end if
      #
      #-------------
      # photolysis #
      #-------------
      if (use_Photolysis(i)) then
        if (use_DOCSorbed(i)) then
          C_Photolysis(id)  = Icpht / I0pht(id,r) * (Cd_Species(id) * kphtd_tc(id) + Cdoc_Species(id) * kphtdoc_tc(id))
        else
          C_Photolysis(id)  = Icpht / I0pht(id,r) * Cd_Species(id)  * kphtd_tc(id)
        end if      
      else
        C_Photolysis(id)    = 0.0
      end if
      #
      #-----------------
      # volatilization #
      #-----------------
      if (use_Volatilization(i) .and. j == 1) then
        C_Volatilization(id) = vv_tc(i) / depth * (Cd_Species(id) - C0(i) / (KH(i,r) / (gas_constant * TwaterK)))
      else
        C_Volatilization(id) = 0.0
      end if
      #
      #-----------------
      # 2nd reaction #
      #-----------------
      if (use_2ndReaction(i)) then 
        C_2ndReaction(id)  = Cd_Species(id)  * Cenv(i)  * kerd_tc(id)
        if (use_BedSediment)    C2_2ndReaction(id) = Cd2_Species(id) * Cenv2(i) * kerd2_tc(id)
        if (use_DOCSorbed(i)) then
          C_2ndReaction(id) = C_2ndReaction(id)  + Cdoc_Species(id)  * Cenv(i)  * kerdoc_tc(id)
          if (use_BedSediment)  C2_2ndReaction(id) = C2_2ndReaction(id) + Cdoc2_Species(id) * Cenv2(i) * kerdoc2_tc(id)
        end if
        if (use_AlgaeSorbed(i)) C_2ndReaction(id)  = C_2ndReaction(id)  + Cap_Species(id)   * Cenv(i)  * kerap_tc(id)
        if (use_POMSorbed(i)) then
          C_2ndReaction(id)  = C_2ndReaction(id)  + Cpom_Species(id)  * Cenv(i)  * kerpom_tc(id)
          if (use_BedSediment) C2_2ndReaction(id) = C2_2ndReaction(id) + Cpom2_Species(id) * Cenv2(i) * kerpom2_tc(id)
        end if
        #
        do k = 1, nGS
          if (use_SolidSorbed(i,k)) then
            C_2ndReaction(id)  = C_2ndReaction(id)  + Cp_Species(id,k)  * Cenv(i)  * kerp_tc(id)
            if (use_BedSediment)  C2_2ndReaction(id) = C2_2ndReaction(id) + Cp2_Species(id,k) * Cenv2(i) * kerp2_tc(id)
          end if
        end do
      else
        C_2ndReaction(id)  = 0.0
        if (use_BedSediment) C2_2ndReaction(id) = 0.0
      end if
      #
      #------------------------------------------
      # transformations
      do k = 1, nTransformProduct
        #
        #-----------------------------------
        # decay transformation #
        #-----------------------------------
        if (nTransform_Decay(i,k) > 0 .and. nTransform_Decay(i,k) <= nC) then
          C_Transform_Decay(id,k)  = C_Decay(id)  * y1(id,k,r)
          if (use_BedSediment) C2_Transform_Decay(id,k) = C2_Decay(id) * y1(id,k,r)
        end if
        #
        #----------------------------
        # hydrolysis transformation #
        #----------------------------
        if (nTransform_Hydrolysis(i,k) > 0 .and. nTransform_Hydrolysis(i,k) <= nC) then
          C_Transform_Hydrolysis(id,k)  = C_Acid_Hydrolysis  * yha(id,k,r) + C_Neutral_Hydrolysis  * yhn(id,k,r) + C_Base_Hydrolysis  * yhb(id,k,r)
          if (use_BedSediment) C2_Transform_Hydrolysis(id,k) = C2_Acid_Hydrolysis * yha(id,k,r) + C2_Neutral_Hydrolysis * yhn(id,k,r) + C2_Base_Hydrolysis * yhb(id,k,r)
        end if
        #
        #----------------------------
        # photolysis transformation #
        #----------------------------
        if (nTransform_Photolysis(i,k) > 0 .and. nTransform_Photolysis(i,k) <= nC)  C_Transform_Photolysis(id,k) = C_Photolysis(id) * ypht(id,k,r)
        #
        #--------------------------------
        # 2nd reaction transformation #
        #--------------------------------
        if (nTransform_Reaction(i,k) > 0 .and. nTransform_Reaction(i,k) <= nC) then
          C_Transform_Reaction(id,k)  = C_2ndReaction(id)  * yer(id,k,r)
          if (use_BedSediment) C2_Transform_Reaction(id,k) = C2_2ndReaction(id) * yer(id,k,r)
        end if
      end do
      #
    end do
    # end of pathway for each species
  end do
  # end of pathway for each contaminant
end subroutine

#=========================================================================================================================== 
# Compute kinetic rate changes 
subroutine ContaminantKinetics()
  real(R8) :: kad_tmp, Adsorption_Desorption
  integer  :: iTransformTo
  #
  # initialize change in contaminants 
  dCdt    = 0.0
  dCddt   = 0.0
  dCdocdt = 0.0
  dCpomdt = 0.0
  dCapdt  = 0.0
  dCpdt   = 0.0
  #
  if (use_BedSediment) then
    dC2dt    = 0.0
    dCd2dt   = 0.0
    dCdoc2dt = 0.0
    dCpom2dt = 0.0
    dCp2dt   = 0.0
  end if
  #
  do i = 1, nC
    #
    if (.not. use_NonEquilibrium(i,1)) dCdt(i) = dCdt(i) + C_Air_Deposition(i)
    #
    do j = 1, nSpecies(i)
      id = parameter_index(i,j)
      #
      # compute kinetic rate equations 
      #
      #-------------------------------------------------------------------------------------------------------------
      # equilibrium: kinetic rate equations of dCdt, dC2dt 
      #-------------------------------------------------------------------------------------------------------------
      if (.not. use_NonEquilibrium(i,1)) then
        dCdt(i) = dCdt(i) - C_C2_Settling(id) - C_Decay(id) - C_Hydrolysis(id) - C_Photolysis(id) - C_Volatilization(id) - C_2ndReaction(id)
        if (use_BedSediment) dCdt(i)  = dCdt(i)  + C_C2_Resuspension(id) + C_C2_Transfer(id) 
      end if 
      if (.not. use_NonEquilibrium(i,2)) then
        if (use_BedSediment) dC2dt(i) = dC2dt(i) + C2_Settling(id) - C2_C_Resuspension(id) + C2_C_Transfer(id) - C2_Burial(id) - C2_Decay(id) - C2_Hydrolysis(id) - C2_2ndReaction(id)
      end if
      #
      #-------------------------------------------------------------------------------------------------------------
      # non-equilibrium: kinetic rate equations of dCddt, dCd2dt, |dCapdt, |dCpomdt, dCpom2dt, |dCpdt, dCp2dt
      #-------------------------------------------------------------------------------------------------------------
      if (use_NonEquilibrium(i,1) .or. use_NonEquilibrium(i,2)) then
        #
        #-----------------
        # air deposition # allocated based on Equilibrium partition. 
        #-----------------
        if (use_NonEquilibrium(i,1)) then
          if (C(i) == 0.0) then
            # all air deposition is given to the dissolved.
            dCddt(i) = dCddt(i) + C_Air_Deposition(i)
          else
            # distribute air deposition according to the current fractions.
            dCddt(i) = dCddt(i) + C_Air_Deposition(i) * Cd_Species(id)   / C(i)
            if (use_DOCSorbed(i))       dCdocdt(i) = dCdocdt(i) + C_Air_Deposition(i) * Cdoc_Species(id) / C(i)
            if (use_AlgaeSorbed(i))     dCapdt(i)  = dCapdt(i)  + C_Air_Deposition(i) * Cap_Species(id)  / C(i)
            if (use_POMSorbed(i))       dCpomdt(i) = dCpomdt(i) + C_Air_Deposition(i) * Cpom_Species(id) / C(i)
            do k = 1, nGS
              if (use_SolidSorbed(i,k)) dCpdt(i,k) = dCpdt(i,k) + C_Air_Deposition(i) * Cp_Species(id,k) / C(i)
            end do
          end if
        end if
        #
        #------------------------------------ 
        # settling, resuspension and burial #
        #------------------------------------
        if (use_AlgaeSorbed(i)) then
          if (use_NonEquilibrium(i,1))                         dCapdt(i)   = dCapdt(i)   - Cap_Species(id) * vsap(r) / depth
          if (use_NonEquilibrium(i,2) .and. use_POMSorbed(i))  dCpom2dt(i) = dCpom2dt(i) + Cap_Species(id) * vsap(r) / h2(r)
        end if
        if (use_POMSorbed(i)) then
          if (use_NonEquilibrium(i,1)) dCpomdt(i)  = dCpomdt(i)  - Cpom_Species(id) * vsom(r) / depth
          if (use_NonEquilibrium(i,2)) dCpom2dt(i) = dCpom2dt(i) + (Cpom_Species(id) * vsom(r) - Cpom2_Species(id) * vb) / h2(r)
        end if
        #
        do k = 1, nGS
          if (use_SolidSorbed(i,k)) then
            if (use_NonEquilibrium(i,1)) then
              dCpdt(i,k) = dCpdt(i,k) - Cp_Species(id,k) * vsp(k) / depth
              if (use_BedSediment)       dCpdt(i,k)  = dCpdt(i,k) + Cp2_Species(id,k) * vrp(k) / depth
            end if
            if (use_NonEquilibrium(i,2)) dCp2dt(i,k) = dCp2dt(i,k) + (Cp_Species(id,k) * vsp(k) - Cp2_Species(id,k) * (vrp(k) + vb)) / h2(r)
          end if 
        end do
        #
        #----------------
        # mass transfer #
        #----------------
        if (use_BedSediment) then
          if (use_NonEquilibrium(i,1)) then
            dCddt(i)   = dCddt(i)   + (Cd2_Species(id)   / Por(r) - Cd_Species(id))   * vm(id,r) / depth
            if (use_DOCSorbed(i)) dCdocdt(i) = dCdocdt(i) + (Cdoc2_Species(id) / Por(r) - Cdoc_Species(id)) * vm(id,r) / depth
          end if
          if (use_NonEquilibrium(i,2)) then
            dCd2dt(i)   = dCd2dt(i)   - (Cd2_Species(id)   / Por(r) - Cd_Species(id))   * vm(id,r) / h2(r)
            if (use_DOCSorbed(i)) dCdoc2dt(i) = dCdoc2dt(i) - (Cdoc2_Species(id) / Por(r) - Cdoc_Species(id)) * vm(id,r) / h2(r)
          end if
        end if
        #
        #--------------------
        # decay             #
        #--------------------
        if (use_nOrderDecay(i)) then
          if (use_NonEquilibrium(i,1)) then
            dCddt(i) = dCddt(i) - Cd_Species(id)**nOrder(id,r) * k1d_tc(id)
            if (use_DOCSorbed(i))       dCdocdt(i) = dCdocdt(i) - Cdoc_Species(id)**nOrder(id,r) * k1doc_tc(id)
            if (use_AlgaeSorbed(i))     dCapdt(i)  = dCapdt(i)  - Cap_Species(id)**nOrder(id,r)  * k1ap_tc(id)
            if (use_POMSorbed(i))       dCpomdt(i) = dCpomdt(i) - Cpom_Species(id)**nOrder(id,r) * k1pom_tc(id)
            do k = 1, nGS
              if (use_SolidSorbed(i,k)) dCpdt(i,k) = dCpdt(i,k) - Cp_Species(id,k)**nOrder(id,r) * k1p_tc(id)
            end do
          end if
          #
          if (use_NonEquilibrium(i,2)) then
            dCd2dt(i) = dCd2dt(i) - Cd2_Species(id)**nOrder(id,r) * k1d2_tc(id)
            if (use_DOCSorbed(i))       dCdoc2dt(i) = dCdoc2dt(i) - Cdoc2_Species(id)**nOrder(id,r) * k1doc2_tc(id)
            if (use_POMSorbed(i))       dCpom2dt(i) = dCpom2dt(i) - Cpom2_Species(id)**nOrder(id,r) * k1pom2_tc(id)
            do k = 1, nGS
              if (use_SolidSorbed(i,k)) dCp2dt(i,k) = dCp2dt(i,k) - Cp2_Species(id,k)**nOrder(id,r) * k1p2_tc(id)
            end do
          end if 
        end if
        #
        #-------------
        # hydrolysis #
        #-------------
        if (use_Hydrolysis(i)) then
          if (use_NonEquilibrium(i,1)) then
            dCddt(i)    = dCddt(i) - Cd_Species(id) * (khbd_tc(id) * COH + khnd_tc(id) + khad_tc(id) * CHH)
            if (use_DOCSorbed(i)) dCdocdt(i)  = dCdocdt(i) - Cdoc_Species(id) * (khbdoc_tc(id) * COH + khndoc_tc(id) + khadoc_tc(id)  * CHH)
          end if
          if (use_NonEquilibrium(i,2)) then
            dCd2dt(i)   = dCd2dt(i) - Cd2_Species(id) * (khbd2_tc(id) * COH + khnd2_tc(id) + khad2_tc(id) * CHH)
            if (use_DOCSorbed(i)) dCdoc2dt(i) = dCdoc2dt(i) - Cdoc2_Species(id) * (khbdoc2_tc(id) * COH + khndoc2_tc(id) + khadoc2_tc(id) * CHH)
          end if
        end if
        #
        #-------------
        # photolysis #
        #-------------
        if (use_Photolysis(i) .and. use_NonEquilibrium(i,1)) then
          dCddt(i)   = dCddt(i) - Cd_Species(id) * kphtd_tc(id) * Icpht / I0pht(id,r)
          if (use_DOCSorbed(i)) dCdocdt(i) = dCdocdt(i) - Cdoc_Species(id) * kphtdoc_tc(id) * Icpht / I0pht(id,r)
        end if
        #
        #-----------------
        # volatilization #
        #-----------------
        if (use_Volatilization(i) .and. use_NonEquilibrium(i,1))  dCddt(i) = dCddt(i) - C_Volatilization(id)
        #
        #-----------------
        # generalized 2nd reaction #
        #-----------------
        if (use_2ndReaction(i)) then
          if (use_NonEquilibrium(i,1)) then
            dCddt(i)   = dCddt(i) - Cd_Species(id) * Cenv(i) * kerd_tc(id)
            if (use_DOCSorbed(i))       dCdocdt(i) = dCdocdt(i) - Cdoc_Species(id) * Cenv(i) * kerdoc_tc(id)
            if (use_AlgaeSorbed(i))     dCapdt(i)  = dCapdt(i)  - Cap_Species(id)  * Cenv(i) * kerap_tc(id)
            if (use_POMSorbed(i))       dCpomdt(i) = dCpomdt(i) - Cpom_Species(id) * Cenv(i) * kerpom_tc(id)
            do k = 1, nGS
              if (use_SolidSorbed(i,k)) dCpdt(i,k) = dCpdt(i,k) - Cp_Species(id,k) * Cenv(i) * kerp_tc(id)
            end do
          end if
          if (use_NonEquilibrium(i,2)) then
            dCd2dt(i)   = dCd2dt(i) - Cd2_Species(id) * Cenv2(i) * kerd2_tc(id)
            if (use_DOCSorbed(i))       dCdoc2dt(i) = dCdoc2dt(i) - Cdoc2_Species(id) * Cenv2(i) * kerdoc2_tc(id)
            if (use_POMSorbed(i))       dCpom2dt(i) = dCpom2dt(i) - Cpom2_Species(id) * Cenv2(i) * kerpom2_tc(id)
            do k = 1, nGS
              if (use_SolidSorbed(i,k)) dCp2dt(i,k) = dCp2dt(i,k) - Cp2_Species(id,k) * Cenv2(i) * kerp2_tc(id)
            end do
          end if 
        end if
        #
        #--------------------------------------------------------------------
        # adsorption and desorption for non-equilibrium partition           #
        #--------------------------------------------------------------------
        # Need to be improved??? Junna.
        if (use_NonEquilibrium(i,1)) then
          if (use_AlgaeSorbed(i)) then
            kad_tmp = kadap_tc(id) * (Apd / 1.0E3 * qcap(id,r) - Cap_Species(id))
            if (kad_tmp > 0.0 .and. kad_tmp * dt < 1.0) then
              Adsorption_DeSorption = - kdeap_tc(id) * Cap_Species(id) + kad_tmp * Cd_Species(id)
            else if (kad_tmp <= 0.0) then
              Adsorption_DeSorption = - kdeap_tc(id) * Cap_Species(id)
            else if (kad_tmp * dt >= 1.0) then
              Adsorption_DeSorption = - kdeap_tc(id) * Cap_Species(id) + Cd_Species(id) / dt
            end if
            dCapdt(i) = dCapdt(i) + Adsorption_DeSorption
            dCddt(i)  = dCddt(i)  - Adsorption_DeSorption
          end if
          #
          if (use_POMSorbed(i)) then
            kad_tmp = kadpom_tc(id) * (POM / 1.0E3 * qcpom(id,r) - Cpom_Species(id))
            if (kad_tmp > 0.0 .and. kad_tmp * dt < 1.0) then
              Adsorption_DeSorption = - kdepom_tc(id) * Cpom_Species(id) + kad_tmp * Cd_Species(id) 
            else if (kad_tmp <= 0.0) then
              # no available adsorption
              Adsorption_DeSorption = - kdepom_tc(id) * Cpom_Species(id)
            else if (kad_tmp * dt >= 1.0) then
              # maximum adsorptable contaminants is Cd_Species(id)
              Adsorption_DeSorption = - kdepom_tc(id) * Cpom_Species(id) + Cd_Species(id) / dt 
            end if
            dCpomdt(i)    = dCpomdt(i) + Adsorption_DeSorption
            dCddt(i)      = dCddt(i)   - Adsorption_DeSorption
          end if
          #
          do k = 1, nGS
            if (use_SolidSorbed(i,k)) then
              kad_tmp = kadp_tc(id,k) * (Solid(k) / 1.0E3 * qcp(id,k,r) - Cp_Species(id,k))
              if (kad_tmp > 0.0 .and. kad_tmp * dt < 1.0) then
                Adsorption_DeSorption = - kdep_tc(id,k) * Cp_Species(id,k) + kad_tmp * Cd_Species(id)
              else if (kad_tmp <= 0.0) then
                Adsorption_DeSorption = - kdep_tc(id,k) * Cp_Species(id,k)
              else if (kad_tmp * dt >= 1.0) then
                Adsorption_DeSorption = - kdep_tc(id,k) * Cp_Species(id,k) + Cd_Species(id) / dt
              end if
              dCpdt(i,k) = dCpdt(i,k) + Adsorption_DeSorption
              dCddt(i)   = dCddt(i)   - Adsorption_DeSorption
            end if
          end do
        end if
        #
        if (use_NonEquilibrium(i,2)) then
          if (use_POMSorbed(i)) then
            kad_tmp = kadpom2_tc(id) / Por(r) * (POM2 / 1.0E3 * qcpom2(id,r) - Cpom2_Species(id))
            if (kad_tmp > 0.0 .and. kad_tmp * dt < 1.0) then 
              Adsorption_DeSorption = - kdepom2_tc(id) * Cpom2_Species(id) + kad_tmp * Cd2_Species(id)
            else if (kad_tmp <= 0.0) then
              Adsorption_DeSorption = - kdepom2_tc(id) * Cpom2_Species(id)
            else if (kad_tmp * dt >= 1.0) then
              Adsorption_DeSorption = - kdepom2_tc(id) * Cpom2_Species(id) + Cd2_Species(id) / dt
            end if
            dCpom2dt(i) = dCpom2dt(i) + Adsorption_DeSorption
            dCd2dt(i)   = dCd2dt(i)   - Adsorption_DeSorption 
          end if
          #
          do k = 1, nGS
            if (use_SolidSorbed(i,k)) then
              kad_tmp = kadp2_tc(id,k) / Por(r) * (Solid2(k) / 1.0E3 * qcp2(id,k,r) - Cp2_Species(id,k))
              if (kad_tmp > 0.0 .and. kad_tmp * dt < 1.0) then
                Adsorption_DeSorption = - kdep2_tc(id,k) * Cp2_Species(id,k) + kad_tmp * Cd2_Species(id)
              else if (kad_tmp <= 0.0) then
                Adsorption_DeSorption = - kdep2_tc(id,k) * Cp2_Species(id,k)
              else if (kad_tmp * dt >= 1.0) then
                Adsorption_DeSorption = - kdep2_tc(id,k) * Cp2_Species(id,k) + Cd2_Species(id) / dt
              end if
              dCp2dt(i,k) = dCp2dt(i,k) + Adsorption_DeSorption
              dCd2dt(i)   = dCd2dt(i)   - Adsorption_DeSorption
            end if
          end do
        end if
        #
      end if
      # 
      #----------------------------------------------
      # add transformations to kinetic rate equation #
      #----------------------------------------------
      # assign all transform products to dissolved part for non-equilibrium
      do k = 1, nTransformProduct
        #
        # decay transform
        iTransformTo = nTransform_Decay(i,k)
        if (iTransformTo > 0 .and. iTransformTo <= nC) then
          if (use_NonEquilibrium(iTransformTo,1)) then
            dCddt(iTransformTo) = dCddt(iTransformTo) + C_Transform_Decay(id,k)
          else
            dCdt(iTransformTo)  = dCdt(iTransformTo)  + C_Transform_Decay(id,k)
          end if
          if (use_BedSediment) then
            if (use_NonEquilibrium(iTransformTo,2)) then
              dCd2dt(iTransformTo) = dCd2dt(iTransformTo) + C2_Transform_Decay(id,k)
            else
              dC2dt(iTransformTo)  = dC2dt(iTransformTo)  + C2_Transform_Decay(id,k)
            end if
          end if
        end if
        #
        # hydrolysis transform
        iTransformTo = nTransform_Hydrolysis(i,k)
        if (iTransformTo > 0 .and. iTransformTo <= nC) then
          if (use_NonEquilibrium(iTransformTo,1)) then
            dCddt(iTransformTo) = dCddt(iTransformTo) + C_Transform_Hydrolysis(id,k)
          else
            dCdt(iTransformTo)  = dCdt(iTransformTo)  + C_Transform_Hydrolysis(id,k)
          end if
          if (use_BedSediment) then
            if (use_NonEquilibrium(iTransformTo,2)) then
              dCd2dt(iTransformTo) = dCd2dt(iTransformTo) + C2_Transform_Hydrolysis(id,k)
            else
              dC2dt(iTransformTo)  = dC2dt(iTransformTo)  + C2_Transform_Hydrolysis(id,k)
            end if
          end if
        end if
        #
        # photolysis transform
        iTransformTo = nTransform_Photolysis(i,k)
        if (iTransformTo > 0 .and. iTransformTo <= nC) then
          if (use_NonEquilibrium(iTransformTo,1)) then
            dCddt(iTransformTo) = dCddt(iTransformTo) + C_Transform_Photolysis(id,k)
          else
            dCdt(iTransformTo)  = dCdt(iTransformTo)  + C_Transform_Photolysis(id,k)
          end if
        end if
        #
        # reaction transform
        iTransformTo = nTransform_Reaction(i,k)
        if (iTransformTo > 0 .and. iTransformTo <= nC) then
          if (use_NonEquilibrium(iTransformTo,1)) then
            dCddt(iTransformTo) = dCddt(iTransformTo) + C_Transform_Reaction(id,k)
          else
            dCdt(iTransformTo)  = dCdt(iTransformTo)  + C_Transform_Reaction(id,k)
          end if
          if (use_BedSediment) then
            if (use_NonEquilibrium(iTransformTo,2)) then
              dCd2dt(iTransformTo) = dCd2dt(iTransformTo) + C2_Transform_Reaction(id,k)
            else
              dC2dt(iTransformTo)  = dC2dt(iTransformTo)  + C2_Transform_Reaction(id,k)
            end if
          end if
        end if
      end do
      #
    end do
    # end for each species
  end do
  # end for each contaminant
  #
  # compute C2, Cd2, Cpom2, Cp2 by adding dC2dt, dCd2dt, dCpom2dt, dCp2dt
  if (use_BedSediment) then
    do i = 1, nC
      if (use_NonEquilibrium(i,2)) then
        Cd2(i) = max(Cd2(i)   + dCd2dt(i)   * dt / Por(r), 0.0)
        if (use_DOCSorbed(i))       Cdoc2(i) = max(Cdoc2(i) + dCdoc2dt(i) * dt / Por(r), 0.0)
        if (use_POMSorbed(i))       Cpom2(i) = max(Cpom2(i) + dCpom2dt(i) * dt, 0.0)
        do k = 1, nGS 
          if (use_SolidSorbed(i,k)) Cp2(i,k) = max(Cp2(i,k) + dCp2dt(i,k) * dt, 0.0)
        end do  
      else
        C2(i) = max(C2(i) + dC2dt(i) * dt, 0.0)
      end if
    end do  
  end if
  #
end subroutine 

#=========================================================================================================================== 
# Call subroutines 
# if ComputeContaminantKinetics() is called right after ComputeContaminantDerivedVariables() for the same WQ cell, 
# no need to call ContaminantPartitions() here.
# ContaminantPartitions() only needs to be called once.
subroutine ComputeContaminantKinetics()
  call ContaminantPartitions()
  call ContaminantTempCorrection()
  call ContaminantPathways()
  call ContaminantKinetics()
end subroutine 

#===========================================================================================================================
# Output pathway
subroutine ContaminantPathwayOutput(na, a)
  integer  :: na
  real(R8) :: a(na)
  #
  do i = 1, nC
    id = parameter_index(i,1)
    if (C_Air_Deposition_index(i) > 0)             a(C_Air_Deposition_index(i))          = C_Air_Deposition(i)
    if (C_Decay_index(i) > 0)                      a(C_Decay_index(i))                   = C_Decay(id)
    if (C_Hydrolysis_index(i) > 0)                 a(C_Hydrolysis_index(i))              = C_Hydrolysis(id)
    if (C_Photolysis_index(i) > 0)                 a(C_Photolysis_index(i))              = C_Photolysis(id)
    if (C_Volatilization_index(i) > 0)             a(C_Volatilization_index(i))          = C_Volatilization(id)
    if (C_2ndReaction_index(i) > 0)                a(C_2ndReaction_index(i))             = C_2ndReaction(id)
    if (C_C2_Settling_index(i) > 0)                a(C_C2_Settling_index(i))             = C_C2_Settling(id)
    do j = 1, nTransformProduct
      if (C_Transform_Decay_index(i,j) > 0)        a(C_Transform_Decay_index(i,j))       = C_Transform_Decay(id,j)
      if (C_Transform_Hydrolysis_index(i,j) > 0)   a(C_Transform_Hydrolysis_index(i,j))  = C_Transform_Hydrolysis(id,j)
      if (C_Transform_Photolysis_index(i,j) > 0)   a(C_Transform_Photolysis_index(i,j))  = C_Transform_Photolysis(id,j)
      if (C_Transform_Reaction_index(i,j) > 0)     a(C_Transform_Reaction_index(i,j))    = C_Transform_Reaction(id,j)
    end do
    #
    if (use_BedSediment) then 
      if (C_C2_Resuspension_index(i) > 0)            a(C_C2_Resuspension_index(i))         = C_C2_Resuspension(id)
      if (C_C2_Transfer_index(i) > 0)                a(C_C2_Transfer_index(i))             = C_C2_Transfer(id)
      #
      if (C2_Decay_index(i) > 0)                     a(C2_Decay_index(i))                  = C2_Decay(id)
      if (C2_Hydrolysis_index(i) > 0)                a(C2_Hydrolysis_index(i))             = C2_Hydrolysis(id)
      if (C2_2ndReaction_index(i) > 0)               a(C2_2ndReaction_index(i))            = C2_2ndReaction(id)
      if (C2_Settling_index(i) > 0)                  a(C2_Settling_index(i))               = C2_Settling(id)
      if (C2_Burial_index(i) > 0)                    a(C2_Burial_index(i))                 = C2_Burial(id)
      if (C2_C_Resuspension_index(i) > 0)            a(C2_C_Resuspension_index(i))         = C2_C_Resuspension(id)
      if (C2_C_Transfer_index(i) > 0)                a(C2_C_Transfer_index(i))             = C2_C_Transfer(id)
      do j = 1, nTransformProduct
        if (C2_Transform_Decay_index(i,j) > 0)       a(C2_Transform_Decay_index(i,j))      = C2_Transform_Decay(id,j)
        if (C2_Transform_Hydrolysis_index(i,j) > 0)  a(C2_Transform_Hydrolysis_index(i,j)) = C2_Transform_Hydrolysis(id,j)
        if (C2_Transform_Reaction_index(i,j) > 0)    a(C2_Transform_Reaction_index(i,j))   = C2_Transform_Reaction(id,j)
      end do
    end if
    #
    if (nSpecies(i) > 1) then
      do j = 2, nSpecies(i)
        id = parameter_index(i,j)
        if (C_Decay_index(i) > 0)                      a(C_Decay_index(i))                   = a(C_Decay_index(i))              + C_Decay(id)
        if (C_Hydrolysis_index(i) > 0)                 a(C_Hydrolysis_index(i))              = a(C_Hydrolysis_index(i))         + C_Hydrolysis(id)
        if (C_Photolysis_index(i) > 0)                 a(C_Photolysis_index(i))              = a(C_Photolysis_index(i))         + C_Photolysis(id)
        if (C_2ndReaction_index(i) > 0)                a(C_2ndReaction_index(i))             = a(C_2ndReaction_index(i))        + C_2ndReaction(id)
        if (C_C2_Settling_index(i) > 0)                a(C_C2_Settling_index(i))             = a(C_C2_Settling_index(i))        + C_C2_Settling(id)
        do k = 1, nTransformProduct
          if (C_Transform_Decay_index(i,k) > 0)        a(C_Transform_Decay_index(i,k))       = a(C_Transform_Decay_index(i,k))       + C_Transform_Decay(id,k)
          if (C_Transform_Hydrolysis_index(i,k) > 0)   a(C_Transform_Hydrolysis_index(i,k))  = a(C_Transform_Hydrolysis_index(i,k))  + C_Transform_Hydrolysis(id,k)
          if (C_Transform_Photolysis_index(i,k) > 0)   a(C_Transform_Photolysis_index(i,k))  = a(C_Transform_Photolysis_index(i,k))  + C_Transform_Photolysis(id,k)
          if (C_Transform_Reaction_index(i,k) > 0)     a(C_Transform_Reaction_index(i,k))    = a(C_Transform_Reaction_index(i,k))    + C_Transform_Reaction(id,k)
        end do
        #
        if (use_BedSediment) then
          if (C_C2_Resuspension_index(i) > 0)            a(C_C2_Resuspension_index(i))         = a(C_C2_Resuspension_index(i))    + C_C2_Resuspension(id)
          if (C_C2_Transfer_index(i) > 0)                a(C_C2_Transfer_index(i))             = a(C_C2_Transfer_index(i))        + C_C2_Transfer(id)
          #
          if (C2_Decay_index(i) > 0)                     a(C2_Decay_index(i))                  = a(C2_Decay_index(i))             + C2_Decay(id)
          if (C2_Hydrolysis_index(i) > 0)                a(C2_Hydrolysis_index(i))             = a(C2_Hydrolysis_index(i))        + C2_Hydrolysis(id)
          if (C2_2ndReaction_index(i) > 0)               a(C2_2ndReaction_index(i))            = a(C2_2ndReaction_index(i))       + C2_2ndReaction(id)
          if (C2_Settling_index(i) > 0)                  a(C2_Settling_index(i))               = a(C2_Settling_index(i))          + C2_Settling(id)
          if (C2_Burial_index(i) > 0)                    a(C2_Burial_index(i))                 = a(C2_Burial_index(i))            + C2_Burial(id)
          if (C2_C_Resuspension_index(i) > 0)            a(C2_C_Resuspension_index(i))         = a(C2_C_Resuspension_index(i))    + C2_C_Resuspension(id)
          if (C2_C_Transfer_index(i) > 0)                a(C2_C_Transfer_index(i))             = a(C2_C_Transfer_index(i))        + C2_C_Transfer(id)
          do k = 1, nTransformProduct
            if (C2_Transform_Decay_index(i,k) > 0)       a(C2_Transform_Decay_index(i,k))      = a(C2_Transform_Decay_index(i,k))      + C2_Transform_Decay(id,k)
            if (C2_Transform_Hydrolysis_index(i,k) > 0)  a(C2_Transform_Hydrolysis_index(i,k)) = a(C2_Transform_Hydrolysis_index(i,k)) + C2_Transform_Hydrolysis(id,k)
            if (C2_Transform_Reaction_index(i,k) > 0)    a(C2_Transform_Reaction_index(i,k))   = a(C2_Transform_Reaction_index(i,k))   + C2_Transform_Reaction(id,k)
          end do
        end if
      end do 
    end if
  end do
end subroutine

#===========================================================================================================================
# Compute derived variables
subroutine ComputeContaminantDerivedVariables
  real(R8) :: Tsolid
  #
  call ContaminantPartitions()
  #
  do i = 1, nC
    #
    # Water column
    if (use_NonEquilibrium(i,1)) then
      #
      # compute C(i)
      if (t > 1.0E-10) then
        C(i) = Cd(i)
        if (use_DOCSorbed(i))         C(i) = C(i) + Cdoc(i)
        if (use_AlgaeSorbed(i))       C(i) = C(i) + Cap(i)  
        if (use_POMSorbed(i))         C(i) = C(i) + Cpom(i)
        do k = 1, nGS
          if (use_SolidSorbed(i,k))   C(i) = C(i) + Cp(i,k)
        end do
      end if
    else
      #
      # compute Cd, Cdoc, Cap, Cpom, Cp
      id        = parameter_index(i,1)
      Cd(i)     = Cd_Species(id)
      Cdoc(i)   = Cdoc_Species(id)
      Cap(i)    = Cap_Species(id)
      Cpom(i)   = Cpom_Species(id)
      do k = 1, nGS 
        Cp(i,k) = Cp_Species(id,k)
      end do
      do j = 2, nSpecies(i)
        id = parameter_index(i,j)
        Cd(i) = Cd(i) + Cd_Species(id)
        if (use_DOCSorbed(i))       Cdoc(i) = Cdoc(i) + Cdoc_Species(id)
        if (use_AlgaeSorbed(i))     Cap(i)  = Cap(i)  + Cap_Species(id)
        if (use_POMSorbed(i))       Cpom(i) = Cpom(i) + Cpom_Species(id)
        do k = 1, nGS
          if (use_SolidSorbed(i,k)) Cp(i,k) = Cp(i,k) + Cp_Species(id,k)
        end do
      end do
      #
      # compute Cion
      if (use_Ionization(i)) then
        Cion(i) = C(i)  * (1.0 - fion(i,1))
      else
        Cion(i) = 0.0
      end if
    end if
    #
    # compute sorbed concentration on total solids: Cpt, Cpts
    Cpt(i) = 0.0
    Tsolid = 0.0
    if (use_AlgaeSorbed(i)) then
      Cpt(i) = Cpt(i) + Cap(i)
      Tsolid = Tsolid + Apd
    end if
    if (use_POMSorbed(i)) then
      Cpt(i) = Cpt(i) + Cpom(i)
      Tsolid = Tsolid + POM
    end if
    do k = 1, nGS
      if (use_SolidSorbed(i,k)) then
        Cpt(i) = Cpt(i) + Cp(i,k)
        Tsolid = Tsolid + Solid(k)
      end if
    end do
    if (Tsolid > 0.0) then
      Cpts(i) = Cpt(i) / (Tsolid / 1.0E6)
    else
      Cpts(i) = 0.0
    end if
    #
    # Sediment layer
    if (use_BedSediment) then
      if (use_NonEquilibrium(i,2)) then
        #
        # compute C2(i)
        if (t > 1.0E-10) then
          C2(i) = Cd2(i) * Por(r)
          if (use_DOCSorbed(i))       C2(i) = C2(i) + Cdoc2(i) * Por(r)
          if (use_POMSorbed(i))       C2(i) = C2(i) + Cpom2(i)
          do k = 1, nGS
            if (use_SolidSorbed(i,k)) C2(i) = C2(i) + Cp2(i,k)
          end do
        end if
      else
        #
        # compute Cd2, Cdoc2, Cpom2, Cp2
        id = parameter_index(i,1)
        Cd2(i)     = Cd2_Species(id)   / Por(r)
        Cdoc2(i)   = Cdoc2_Species(id) / Por(r) 
        Cpom2(i)   = Cpom2_Species(id)
        do k = 1, nGS
          Cp2(i,k) = Cp2_Species(id,k)
        end do
        do j = 2, nSpecies(i)
          id = parameter_index(i,j)
          Cd2(i) = Cd2(i) + Cd2_Species(id) / Por(r)
          if (use_DOCSorbed(i))       Cdoc2(i) = Cdoc2(i) + Cdoc2_Species(id) / Por(r)
          if (use_POMSorbed(i))       Cpom2(i) = Cpom2(i) + Cpom2_Species(id)
          do k = 1, nGS
            if (use_SolidSorbed(i,k)) Cp2(i,k) = Cp2(i,k) + Cp2_Species(id,k)
          end do
        end do
        #
        # compute Cion2
        if (use_Ionization(i)) then
          Cion2(i) = C2(i) * (1.0 - fion2(i,1))
        else
          Cion2(i) = 0.0
        end if
      end if
      #
      # compute sorbed concentration on total solids: Cpt2, Cpts2
      Cpt2(i) = 0.0
      Tsolid  = 0.0
      if (use_POMSorbed(i)) then
        Cpt2(i) = Cpt2(i) + Cpom2(i)
        Tsolid  = Tsolid  + POM2
      end if
      do k = 1, nGS
        if (use_SolidSorbed(i,k)) then
          Cpt2(i) = Cpt2(i) + Cp2(i,k)
          Tsolid  = Tsolid  + Solid2(k)
        end if
      end do
      if (Tsolid > 0.0) then
        Cpts2(i) = Cpt2(i) / (Tsolid / 1.0E6)
      else
        Cpts2(i) = 0.0
      end if
    end if
  end do
end subroutine

#=========================================================================================================================== 
# Newton-Raphson solution
subroutine NewtonRaphson(Cd_cmp)
  real(R8) :: Cd_cmp, Cd_new
  real(R8) :: ea
  integer  :: iteration
  # 
  # Theory of this method
  # xn+1 = xn - f(xn)/f(xn)'
  # 
  # deal with two exceptional conditions
  Cd_cmp = 0.0
  if (f(Cd_cmp) == 0.0) return
  #
  if (IsWaterCell) then
    Cd_cmp = C(i)  * fion(i,j)
  else
    Cd_cmp = C2(i) * fion2(i,j)
  end if
  if (f(Cd_cmp) == 0.0) return
  #
  # compute under normal conditions
  # set the initial estimation of Cd
  if (IsWaterCell .and. use_Langmuir(i,1) .or. .not. IsWaterCell .and. use_Langmuir(i,2)) then
    Cd_cmp = 0.0
  end if
  #
  if (IsWaterCell .and. use_Freundlich(i,1) .or. .not. IsWaterCell .and. use_Freundlich(i,2)) then
    # The smaller Cd_cmp, the better; but non-zero#
    Cd_cmp = Cd_cmp / 10.0
    do while (f(Cd_cmp) > 0.0)
      Cd_cmp = Cd_cmp / 10.0
    end do
  end if
  #
  iteration = 0
  do 
    iteration = iteration + 1
    Cd_new    = Cd_cmp - f(Cd_cmp) / df(Cd_cmp)
    if (Cd_new .ne. 0.0) ea = abs(Cd_new - Cd_cmp) / Cd_new
    Cd_cmp    = Cd_new
    #
    if (iteration >= imax .or. Cd_cmp < 0.0) then
      # It will be determined by RAS GUI??? 
      write(*,*) 'Newton-Raphson method does not converge when computing dissolved contaminant concentration. Try Bisection method or reduce time step'
      stop
    end if
    #
    if (ea < res(r)) exit
  end do
end subroutine

#===========================================================================================================================
# Bisection solution
subroutine Bisection(Cd_cmp)
  real(R8) :: Cd_cmp
  real(R8) :: xlow, xup, xr
  real(R8) :: flow, fup, fr
  integer  :: iteration
  #
  xlow = 0.0
  flow = f(xlow)
  if (flow == 0.0) then
    Cd_cmp = xlow
    return
  end if

  if (IsWaterCell) then
    xup = C(i) * fion(i,j)
  else
    xup = C2(i) * fion2(i,j)
  end if 
  fup  = f(xup)
  if (fup == 0.0) then
    Cd_cmp = xup
    return
  end if

  if (flow * fup > 0.0) then
    write(*,*) 'Bad dissolved contaminant concentration guess & Try a smaller time step'
    stop
  end if
  
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
  
  if (iteration > imax) then
    write(*,*) 'Bisection method does not converge when computing dissolved contaminant concentration.'
    stop
  end if 
end subroutine

'''