!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2018  DFTB+ developers group                                                      !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> update the SCC hamiltonian
module dftbp_hamiltonian
  use dftbp_accuracy, only : dp, lc
  use dftbp_assert
  use dftbp_commontypes, only : TOrbitals
  use dftbp_periodic, only : TNeighbourList
  use dftbp_potentials, only : TPotentials
  use dftbp_shift, only : add_shift, total_shift
  use dftbp_spin, only : getSpinShift
  use dftbp_spinorbit, only : getDualSpinOrbitShift
  use dftbp_dftbplusu, only : TDftbU
  use dftbp_message, only : error
  use dftbp_thirdorder, only : TThirdOrder
  use dftbp_solvation, only : TSolvation
  use dftbp_environment, only : TEnvironment
  use dftbp_scc, only : TScc
  use dftbp_dispersions, only : TDispersionIface

  implicit none

  private
  public :: resetExternalPotentials, getSccHamiltonian, mergeExternalPotentials
  public :: setUpExternalElectricField, resetInternalPotentials, addChargePotentials
  public :: addBlockChargePotentials, TRefExtPot

  !> Container for external potentials
  type :: TRefExtPot
    real(dp), allocatable :: atomPot(:,:)
    real(dp), allocatable :: shellPot(:,:,:)
    real(dp), allocatable :: potGrad(:,:)
  end type TRefExtPot

contains


  !> Sets the external potential components to zero
  subroutine resetExternalPotentials(refExtPot, potential)

    !> Reference external potential (usually set via API)
    type(TRefExtPot), intent(in) :: refExtPot

    !> Potential contributions
    type(TPotentials), intent(inout) :: potential

    if (allocated(refExtPot%atomPot)) then
      potential%extAtom(:,:) = refExtPot%atomPot
    else
      potential%extAtom(:,:) = 0.0_dp
    end if
    if (allocated(refExtPot%shellPot)) then
      potential%extShell(:,:,:) = refExtPot%shellPot
    else
      potential%extShell(:,:,:) = 0.0_dp
    end if
    potential%extBlock(:,:,:,:) = 0.0_dp
    if (allocated(refExtPot%potGrad)) then
      potential%extGrad(:,:) = refExtPot%potGrad
    else
      potential%extGrad(:,:) = 0.0_dp
    end if

  end subroutine resetExternalPotentials


  !> Merges atomic and shell resolved external potentials into blocked one
  subroutine mergeExternalPotentials(orb, species, potential)

    !> Atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species for atoms
    integer, intent(in) :: species(:)

    !> Potential energy contributions
    type(TPotentials), intent(inout) :: potential

    call total_shift(potential%extShell, potential%extAtom, orb, species)
    call total_shift(potential%extBlock, potential%extShell, orb, species)

  end subroutine mergeExternalPotentials


  !> Sets up electric external field
  subroutine setUpExternalElectricField(tEfield, tTimeDepEField, tPeriodic, EFieldStrength,&
      & EFieldVector, EFieldOmega, EFieldPhase, neighbourList, nNeighbourSK, iCellVec,&
      & img2CentCell, cellVec, deltaT, iGeoStep, coord0Fold, coord, extAtomPot, extPotGrad, EField,&
      & absEField)

    !> Whether electric field should be considered at all
    logical, intent(in) :: tEfield

    !> Is there an electric field that varies with geometry step during MD?
    logical, intent(in) :: tTimeDepEField

    !> Is this a periodic geometry
    logical, intent(in) :: tPeriodic

    !> What is the field strength
    real(dp), intent(in) :: EFieldStrength

    !> What is the field direction
    real(dp), intent(in) :: EFieldVector(:)

    !> Is there an angular frequency for the applied field
    real(dp), intent(in) :: EFieldOmega

    !> What is the phase of the field
    integer, intent(in) :: EFieldPhase

    !> Atomic neighbours
    type(TNeighbourList), intent(in) :: neighbourList

    !> Number of neighbours for each atom
    integer, intent(in) :: nNeighbourSK(:)

    !> Index for unit cells
    integer, intent(in) :: iCellVec(:)

    !> Image atom to central cell atom number
    integer, intent(in) :: img2CentCell(:)

    !> Vectors (in units of the lattice constants) to cells of the lattice
    real(dp), intent(in) :: cellVec(:,:)

    !> Time step in MD
    real(dp), intent(in) :: deltaT

    !> Number of the geometry step
    integer, intent(in) :: iGeoStep

    !> Atomic coordinates in central cell
    real(dp), allocatable, intent(in) :: coord0Fold(:,:)

    !> all coordinates
    real(dp), intent(in) :: coord(:,:)

    !> Potentials on atomic sites
    real(dp), intent(inout) :: extAtomPot(:)

    !> Gradient of potential on atomic sites with respect of nucleus positions. Shape: (3, nAtom)
    real(dp), intent(inout) :: extPotGrad(:,:)

    !> Resulting electric field
    real(dp), intent(out) :: EField(:)

    !> Magnitude of the field
    real(dp), intent(out) :: absEField

    integer :: nAtom
    integer :: iAt1, iAt2, iNeigh
    character(lc) :: tmpStr

    if (.not. tEField) then
      EField(:) = 0.0_dp
      absEField = 0.0_dp
      return
    end if

    nAtom = size(nNeighbourSK)

    Efield(:) = EFieldStrength * EfieldVector
    if (tTimeDepEField) then
      Efield(:) = Efield * sin(EfieldOmega * deltaT * real(iGeoStep + EfieldPhase, dp))
    end if
    absEfield = sqrt(sum(Efield**2))
    if (tPeriodic) then
      do iAt1 = 1, nAtom
        do iNeigh = 1, nNeighbourSK(iAt1)
          iAt2 = neighbourList%iNeighbour(iNeigh, iAt1)
          ! overlap between atom in central cell and non-central cell
          if (iCellVec(iAt2) /= 0) then
            ! component of electric field projects onto vector between cells
            if (abs(dot_product(cellVec(:, iCellVec(iAt2)), EfieldVector)) > epsilon(1.0_dp)) then
              write(tmpStr, "(A, I0, A, I0, A)") 'Interaction between atoms ', iAt1, ' and ',&
                  & img2CentCell(iAt2), ' crosses the saw-tooth discontinuity in the electric&
                  & field.'
              call error(tmpStr)
            end if
          end if
        end do
      end do
      do iAt1 = 1, nAtom
        extAtomPot(iAt1) = extAtomPot(iAt1) + dot_product(coord0Fold(:, iAt1), Efield)
      end do
    else
      do iAt1 = 1, nAtom
        extAtomPot(iAt1) = extAtomPot(iAt1) + dot_product(coord(:, iAt1), Efield)
      end do
    end if
    extPotGrad(:,:) = extPotGrad + spread(EField, 2, nAtom)

  end subroutine setUpExternalElectricField


  !> Reset internal potential related quantities
  subroutine resetInternalPotentials(tDualSpinOrbit, xi, orb, species, potential)

    !> Is dual spin orbit being used (block potentials)
    logical, intent(in) :: tDualSpinOrbit

    !> Spin orbit constants if required
    real(dp), allocatable, intent(in) :: xi(:,:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> chemical species
    integer, intent(in) :: species(:)

    !> potentials in the system
    type(TPotentials), intent(inout) :: potential

    @:ASSERT(.not. tDualSpinOrbit .or. allocated(xi))

    potential%intAtom(:,:) = 0.0_dp
    potential%intShell(:,:,:) = 0.0_dp
    potential%intBlock(:,:,:,:) = 0.0_dp
    potential%orbitalBlock(:,:,:,:) = 0.0_dp
    potential%iOrbitalBlock(:,:,:,:) = 0.0_dp
    if (tDualSpinOrbit) then
      call getDualSpinOrbitShift(potential%iOrbitalBlock, xi, orb, species)
    end if

  end subroutine resetInternalPotentials


  !> Add potentials coming from point charges.
  subroutine addChargePotentials(env, sccCalc, updateScc, qInput, q0, chargePerShell, orb, species,&
      & neighbourList, img2CentCell, spinW, solvation, thirdOrd, potential, dispersion)

    !> Environment settings
    type(TEnvironment), intent(inout) :: env

    !> SCC module internal variables
    type(TScc), intent(inout) :: sccCalc

    !> Whether the charges in the scc calculator should be updated before obtaining the potential
    logical, intent(in) :: updateScc
    
    !> Input atomic populations
    real(dp), intent(in) :: qInput(:,:,:)

    !> reference atomic occupations
    real(dp), intent(in) :: q0(:,:,:)

    !> charges per atomic shell
    real(dp), intent(in) :: chargePerShell(:,:,:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> species of all atoms
    integer, target, intent(in) :: species(:)

    !> neighbours to atoms
    type(TNeighbourList), intent(in) :: neighbourList

    !> map from image atom to real atoms
    integer, intent(in) :: img2CentCell(:)

    !> spin constants
    real(dp), intent(in), allocatable :: spinW(:,:,:)

    !> Solvation mode
    class(TSolvation), allocatable, intent(inout) :: solvation

    !> third order SCC interactions
    type(TThirdOrder), allocatable, intent(inout) :: thirdOrd

    !> Potentials acting
    type(TPotentials), intent(inout) :: potential

    !> Dispersion interactions object
    class(TDispersionIface), allocatable, intent(in) :: dispersion

    ! local variables
    real(dp), allocatable :: atomPot(:,:)
    real(dp), allocatable :: shellPot(:,:,:)
    integer, pointer :: pSpecies0(:)
    integer :: nAtom, nSpin

    nAtom = size(qInput, dim=2)
    nSpin = size(qInput, dim=3)
    pSpecies0 => species(1:nAtom)

    allocate(atomPot(nAtom, nSpin))
    allocate(shellPot(orb%mShell, nAtom, nSpin))

    if (updateScc) then
      call sccCalc%updateCharges(env, qInput, q0, orb, species)
    end if
    call sccCalc%updateShifts(env, orb, species, neighbourList%iNeighbour, img2CentCell)
    call sccCalc%getShiftPerAtom(atomPot(:,1))
    call sccCalc%getShiftPerL(shellPot(:,:,1))

    if (allocated(dispersion)) then
      call dispersion%addPotential(atomPot(:,1))
    end if

    potential%intAtom(:,1) = potential%intAtom(:,1) + atomPot(:,1)
    potential%intShell(:,:,1) = potential%intShell(:,:,1) + shellPot(:,:,1)

    if (allocated(thirdOrd)) then
      call thirdOrd%updateCharges(pSpecies0, neighbourList, qInput, q0, img2CentCell, orb)
      call thirdOrd%getShifts(atomPot(:,1), shellPot(:,:,1))
      potential%intAtom(:,1) = potential%intAtom(:,1) + atomPot(:,1)
      potential%intShell(:,:,1) = potential%intShell(:,:,1) + shellPot(:,:,1)
    end if

    if (allocated(solvation)) then
      call solvation%updateCharges(env, pSpecies0, neighbourList, qInput, q0, img2CentCell, orb)
      call solvation%getShifts(atomPot(:,1), shellPot(:,:,1))
      potential%intAtom(:,1) = potential%intAtom(:,1) + atomPot(:,1)
      potential%intShell(:,:,1) = potential%intShell(:,:,1) + shellPot(:,:,1)
    end if

    if (nSpin /= 1 .and. allocated(spinW)) then
      call getSpinShift(shellPot, chargePerShell, species, orb, spinW)
      potential%intShell = potential%intShell + shellPot
    end if

    call total_shift(potential%intShell, potential%intAtom, orb, species)
    call total_shift(potential%intBlock, potential%intShell, orb, species)

  end subroutine addChargePotentials


  !> Add potentials comming from on-site block of the dual density matrix.
  subroutine addBlockChargePotentials(qBlockIn, qiBlockIn, dftbU, tImHam, species, orb, potential)

    !> block input charges
    real(dp), allocatable, intent(in) :: qBlockIn(:,:,:,:)

    !> imaginary part
    real(dp), allocatable, intent(in) :: qiBlockIn(:,:,:,:)

    !> is this a +U calculation
    type(TDftbU), intent(in), allocatable :: dftbU

    !> does the hamiltonian have an imaginary part in real space?
    logical, intent(in) :: tImHam

    !> chemical species of all atoms
    integer, intent(in) :: species(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

    !> potentials acting in system
    type(TPotentials), intent(inout) :: potential

    if (allocated(dftbU)) then
      if (tImHam) then
        call dftbU%getDftbUShift(potential%orbitalBlock, potential%iorbitalBlock, qBlockIn,&
            & qiBlockIn, species,orb)
      else
        call dftbU%getDftbUShift(potential%orbitalBlock, qBlockIn, species, orb)
      end if
      potential%intBlock = potential%intBlock + potential%orbitalBlock
    end if

  end subroutine addBlockChargePotentials


  !> Returns the Hamiltonian for the given scc iteration
  subroutine getSccHamiltonian(H0, over, nNeighbourSK, neighbourList, species, orb, iSparseStart,&
      & img2CentCell, potential, isREKS, ham, iHam)

    !> non-SCC hamiltonian (sparse)
    real(dp), intent(in) :: H0(:)

    !> overlap (sparse)
    real(dp), intent(in) :: over(:)

    !> Number of atomic neighbours
    integer, intent(in) :: nNeighbourSK(:)

    !> list of atomic neighbours
    type(TNeighbourList), intent(in) :: neighbourList

    !> species of atoms
    integer, intent(in) :: species(:)

    !> atomic orbital information
    type(TOrbitals), intent(in) :: orb

    !> Index for atomic blocks in sparse data
    integer, intent(in) :: iSparseStart(:,:)

    !> image atoms to central cell atoms
    integer, intent(in) :: img2CentCell(:)

    !> potential acting on sustem
    type(TPotentials), intent(in) :: potential

    !> Is this DFTB/SSR formalism
    logical, intent(in) :: isREKS

    !> resulting hamitonian (sparse)
    real(dp), intent(inout) :: ham(:,:)

    !> imaginary part of hamiltonian (if required, signalled by being allocated)
    real(dp), allocatable, intent(inout) :: iHam(:,:)

    integer :: nAtom

    nAtom = size(orb%nOrbAtom)

    if (.not. isREKS) then
      ham(:,:) = 0.0_dp
      ham(:,1) = h0
    end if

    call add_shift(ham, over, nNeighbourSK, neighbourList%iNeighbour, species, orb, iSparseStart,&
        & nAtom, img2CentCell, potential%intBlock)

    if (allocated(iHam)) then
      iHam(:,:) = 0.0_dp
      call add_shift(iHam, over, nNeighbourSK, neighbourList%iNeighbour, species, orb,&
          & iSparseStart, nAtom, img2CentCell, potential%iorbitalBlock)
    end if

  end subroutine getSccHamiltonian


end module dftbp_hamiltonian
