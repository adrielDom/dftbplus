!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'
#:include 'error.fypp'

module dftbp_loscorrection
  use dftbp_assert
  use dftbp_accuracy
  use dftbp_commontypes
  use dftbp_message
  use dftbp_blasroutines
  use dftbp_scc, only : TScc
  use dftbp_transcharges, only : transq
  implicit none
  public

  !> Data type for localised orbital scaling correction
  type :: TLOSCorrection

    !> energy penalty function
    character :: penalty

    !> parameters of penalty function
    real(dp), allocatable :: PenaltyParam(:)

    !> atomic spin constants
    real(dp), allocatable :: spinW(:)

    !> SCF-LOSC Calculation?
    logical :: tSCF

  end type TLOSCorrection

contains

  !> This subroutine calculates the LO scaling correction term of total energy
  subroutine getElosc(e_losc, sccCalc, nNeighbor, iNeighbor, img2CentCell, ind,&
      & S, cc, filling, uu, species, LOSC)

    !> LOSC energy (output)
    real(dp), intent(out) :: e_losc(:)

    !> Self-consistent charge module settings
    type(TScc), intent(in) :: sccCalc

    integer,  intent(in)  :: nNeighbor(:)
    integer,  intent(in)  :: iNeighbor(0:,:)
    integer,  intent(in)  :: img2CentCell(:)
    integer,  intent(in)  :: ind(:)
    real(dp), intent(in)  :: S(:,:), cc(:,:,:)
    real(dp), intent(in)  :: filling(:,:)
    real(dp), intent(in)  :: uu(:,:,:)

    !> chemical species
    integer, intent(in) :: species(:)

    !> Container for LOSC calculation data
    type(TLOSCorrection), allocatable, intent(in) :: LOSC

    real(dp), allocatable :: gamma(:,:)
    real(dp), allocatable :: lambdaQQ(:,:), lambdaSqrQQ(:,:)
    real,     allocatable :: occNr(:,:)

    integer :: nAtom, nOrb, nSpin
    integer :: iAt1, iAt2, iAt2f, iNeigh, iSp1
    real(dp) :: rTmp, rTmpSpin
    real(dp) :: term1, term2, term3, term4

    nAtom = size(e_losc)
    nOrb = size(filling, dim=1)
    nSpin = size(filling, dim=2)

    allocate(gamma(nAtom, nAtom))
    allocate(lambdaQQ(nAtom, nAtom))
    allocate(lambdaSqrQQ(nAtom, nAtom))

    allocate(occNr(nOrb, nSpin))
    occNr(:,:) = real(filling(:,:), 4)

    e_losc = 0.0_dp

    term1 = 0.0_dp
    term2 = 0.0_dp
    term3 = 0.0_dp
    term4 = 0.0_dp

    call sccCalc%getAtomicGammaMatrix(gamma, iNeighbor, img2CentCell)
    call getlambdaQQ(lambdaQQ,uu,ind,S,cc,occNr,nNeighbor,iNeighbor,img2CentCell)
    call getlambdaSqrQQ(lambdaSqrQQ,uu,ind,S,cc,occNr,nNeighbor,iNeighbor,img2CentCell)


    do iAt1 = 1, nAtom
      do iNeigh = 0, nNeighbor(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)

        rTmp = gamma(iAt1,iAt2f) * lambdaQQ(iAt1, iAt2f)&
            & - gamma(iAt1,iAt2f) * lambdaSqrQQ(iAt1, iAt2f)

        ! term1 and term2 just for debugging
        term1 = term1 + 0.5_dp * gamma(iAt1, iAt2f) * lambdaQQ(iAt1, iAt2f)
        term2 = term2 + 0.5_dp * gamma(iAt1, iAt2f) * lambdaSqrQQ(iAt1, iAt2f)

        e_losc(iAt1) = e_losc(iAt1) + 0.5_dp * rTmp

        if (iAt1 /= iAt2f) then
          e_losc(iAt2f) = e_losc(iAt2f) + 0.5_dp * rTmp
        end if

      end do
      if (nSpin /= 1) then
        iSp1 = species(iAt1)
        rTmpSpin = LOSC%spinW(iSp1) * lambdaQQ(iAt1, iAt1)&
            & - LOSC%spinW(iSp1) * lambdaSqrQQ(iAt1, iAt1)
        e_losc(iAt1) = e_losc(iAt1) + 0.5_dp * rTmpSpin

        term3 = term3 + 0.5_dp * LOSC%spinW(iSp1) * lambdaQQ(iAt1, iAt1)
        term4 = term4 + 0.5_dp * LOSC%spinW(iSp1) * lambdaSqrQQ(iAt1, iAt1)
      end if
    end do

    print *, "term1", term1
    print *, "term2", term2
    print *, "term3", term3
    print *, "term4", term4

  end subroutine getElosc


  subroutine getOrbitallets(uu, ener, nNeighbor, iNeighbor, img2CentCell, ind,&
      & S, cc, filling, rat, LOSC)
    real(dp), intent(out) :: uu(:,:,:)
    real(dp), intent(in) :: ener(:,:)
    integer, intent(in)   :: nNeighbor(:)
    integer,  intent(in)  :: iNeighbor(0:,:)
    integer,  intent(in)  :: img2CentCell(:)
    integer,  intent(in)  :: ind(:)
    real(dp), intent(in)  :: S(:,:), cc(:,:,:)
    real(dp), intent(in)  :: filling(:,:)
    real(dp), intent(in)  :: rat(:,:)

    !> Container for LOSC calculation data
    type(TLOSCorrection), allocatable, intent(in) :: LOSC

    logical :: tFconverged
    integer :: iSpin, pp, ss, qq, rr, ii, iFloop, nFloop
    integer :: nOrb, nSpin
    real (dp) :: D0, D1, D2, C0, C1
    real (dp) :: Aq, Bq, Ar, Br
    real (dp) :: wqs, wrs
    real (dp) :: fqr, f_der1, f_der2, fqr_tmp, F, F_new
    real (dp) :: phi(20), phi_new, phi_tmp, delta_phi
    real (dp), allocatable :: norm(:,:), U_tmp(:,:)

    @:ASSERT(allocated(LOSC))

    nOrb = size(uu, dim=2)
    nSpin = size(uu, dim=3)

    allocate(norm(nOrb, nSpin))
    allocate(U_tmp(nOrb, nSpin))

    tFconverged = .false.
    nFloop = 100
    uu = 0.0_dp
    F = 0.0_dp

    !initializing the LO orbitals as the KS orbitals
    do pp = 1, nOrb
      uu(pp,pp,:) = 1.0_dp
    end do

    iFloop = 1
    loopF: do while (iFloop <= nFloop)
      print *, "============================"
      print *, "i =", iFloop
      print *, "============================"

      F_new = 0.0_dp
      do iSpin = 1, nSpin
        do qq = 1, nOrb - 1
          do rr = qq + 1, nOrb

            do ii = 1, 20
              phi(ii) = rand(0) * 6.28
            end do

            print *, "spin =", iSpin
            print *, "q =", qq
            print *, "r =", rr
            print *, "----------------"

            call getLODipole(Aq, qq, iSpin, rat, ind, S, cc, uu)
            call getLODipoleSqr(Bq, qq, iSpin, rat, ind, S, cc, uu)
            call getLODipole(Ar, rr, iSpin, rat, ind, S, cc, uu)
            call getLODipoleSqr(Br, rr, iSpin, rat, ind, S, cc, uu)

            C0 = Aq - Bq**2 + Ar - Br**2
            C1 = 0.5_dp*( Bq**2 + Br**2 )

            D0 = 0.0_dp
            D1 = 0.0_dp
            D2 = 0.0_dp
            do ss = 1, nOrb

              wqs = penalty(abs(ener(qq, iSpin) - ener(ss, iSpin)), LOSC)
              wrs = penalty(abs(ener(rr, iSpin) - ener(ss, iSpin)), LOSC)

              D0 = D0 + 0.5 * (wqs + wrs) * (uu(qq, ss, iSpin)**2 + uu(rr, ss, iSpin)**2)
              D1 = D1 + (wqs - wrs) * uu(qq, ss, iSpin) * uu(rr, ss, iSpin)
              D2 = D2 + 0.5 * (wqs - wrs) * (uu(qq, ss, iSpin)**2 - uu(rr, ss, iSpin)**2)
            end do

            fqr = C0 + D0 + C1 * dsin(phi(1))**2 +&
                & D1 * dsin(phi(1)) + D2 * dcos(phi(1))
            phi_new = phi(1)
            do ii = 1, 20

              f_der1 = 2 * C1 * dsin(phi(ii)) * dcos(phi(ii)) +&
                  & D1 * dcos(phi(ii)) - D2 * dsin(phi(ii))
              f_der2 = 2 * C1 - 4 * C1 * dsin(phi(ii))**2 -&
                  & D1 * dsin(phi(ii)) - D2 * dcos(phi(ii))

              delta_phi = -1.0_dp * f_der1/f_der2
              phi_tmp = phi(ii) + delta_phi

              fqr_tmp = C0 + D0 + C1 * dsin(phi_tmp)**2 + D1 * dsin(phi_tmp)&
                  & + D2 * dcos(phi_tmp)
              if (fqr_tmp < fqr) then
                fqr = fqr_tmp
                phi_new = phi_tmp
              end if
            end do

            F_new = F_new + fqr
            U_tmp(:, iSpin) = uu(qq, :, iSpin) * dcos(0.5 * phi_new) +&
                & uu(rr, :, iSpin) * dsin(0.5 * phi_new)

            uu(rr, :, iSpin) = -1.0_dp * uu(qq, :, iSpin) * dsin(0.5 * phi_new) +&
                & uu(rr, :, iSpin) * dcos(0.5 * phi_new)
            uu(qq, :, iSpin) = U_tmp(:, iSpin)

            norm(:,:) = sqrt(sum(uu**2, dim=2))

            print *, "phi:", phi
            print *, "phi_new:", phi_new
            print *, "fqr:", fqr
            print *, "Uq:", uu(qq,:,iSpin)
            print *, "norm:", norm(qq,iSpin)
            print *, "Ur:", uu(rr,:,iSpin)
            print *, "norm:", norm(rr,iSpin)

          end do
        end do
      end do

      phi = phi_new
      F_new = F_new/(nOrb - 1)
      print *, "F:", F_new

      tFconverged = (abs(F_new - F) < 0.00001_dp)

      if (tFconverged) then
        exit loopF
      end if

      iFloop = iFloop + 1
      F = F_new

    end do loopF

  contains

    function penalty(x,losc) result(w)
      real(dp), intent(in) :: x
      type(TLOSCorrection), intent(in) :: losc
      real(dp) :: w

      real(dp) :: ener_frac, tmp
      real(dp) :: R0, eps0, gamma, eta

      ! parameters of the penalty function
      R0 = losc%PenaltyParam(1)
      eps0 = losc%PenaltyParam(2)
      gamma = losc%PenaltyParam(3)
      eta = losc%PenaltyParam(4)

      ener_frac = x / eps0
      tmp = 1 - exp( -1.0_dp * ener_frac**eta )

      if (x < eps0) then
        w = tmp * R0**2
      else
        w = tmp * R0**2 * ener_frac**gamma
      end if

    end function penalty

  end subroutine getOrbitallets


  subroutine getLODipole(B, qq, sigma, rat, ind, S, cc, uu)
    real(dp), intent(out)          :: B
    integer,  intent(in)           :: qq
    integer,  intent(in)           :: sigma
    real(dp), intent(in)           :: rat(:,:)
    integer,  intent(in)           :: ind(:)
    real(dp), intent(in)           :: S(:,:), cc(:,:,:)
    real(dp), intent(in)           :: uu(:,:,:)

    integer :: nOrb, ss, tt
    real(dp) :: d_st

    nOrb = size(cc, dim=1)

    B = 0.0_dp

    do ss = 1, nOrb
      do tt = 1, nOrb
        call getTransDipole(d_st, ss, tt, sigma, rat, ind, S, cc)
        B = B + uu(qq, ss, sigma) * uu(qq, tt, sigma) * d_st
      end do
    end do

  end subroutine getLODipole


  subroutine getLODipoleSqr(A, qq, sigma, rat, ind, S, cc, uu)
    real(dp), intent(out)          :: A
    integer,  intent(in)           :: qq
    integer,  intent(in)           :: sigma
    real(dp), intent(in)           :: rat(:,:)
    integer,  intent(in)           :: ind(:)
    real(dp), intent(in)           :: S(:,:), cc(:,:,:)
    real(dp), intent(in)           :: uu(:,:,:)

    integer :: nOrb, ss, tt
    real(dp) :: d2_st

    nOrb = size(cc, dim=1)

    A = 0.0_dp

    do ss = 1, nOrb
      do tt = 1, nOrb
        call getTransDipoleSqr(d2_st, ss, tt, sigma, rat, ind, S, cc)
        A = A + uu(qq, ss, sigma) * uu(qq, tt, sigma) * d2_st
      end do
    end do

  end subroutine getLODipoleSqr

  subroutine getTransDipole(d_st, ss, tt, sigma, X_A, ind, S, cc)
    real(dp), intent(out)          :: d_st
    integer,  intent(in)           :: ss, tt
    integer,  intent(in)           :: sigma
    real(dp), intent(in)           :: X_A(:,:)
    integer,  intent(in)           :: ind(:)
    real(dp), intent(in)           :: S(:,:), cc(:,:,:)

    integer :: nOrb, nAtom, iAt
    logical :: updwn
    real(dp), allocatable :: stimc(:,:,:)
    real(dp), allocatable :: qst(:), R_A(:)

    nOrb = size(cc, dim=1)
    nAtom = size(X_A, dim=2)
    allocate(stimc(norb, norb, 2))
    allocate(qst(nAtom))
    allocate(R_A(nAtom))

    call symm(stimc(:,:,sigma), "L", S, cc(:,:,sigma))
    updwn = (sigma == 1)

    qst(:) = transq(ss, tt, ind, updwn, stimc, cc)
    R_A(:) = sqrt(sum(X_A**2, dim=1))

    d_st = 0.0_dp
    do iAt = 1, nAtom
      d_st = d_st + R_A(iAt)*qst(iAt)
    end do
!    d_st = matmul(R_A,qst)

  end subroutine getTransDipole


  subroutine getTransDipoleSqr(d2_st, ss, tt, sigma, X_A, ind, S, cc)
    real(dp), intent(out)          :: d2_st
    integer,  intent(in)           :: ss, tt
    integer,  intent(in)           :: sigma
    real(dp), intent(in)           :: X_A(:,:)
    integer,  intent(in)           :: ind(:)
    real(dp), intent(in)           :: S(:,:), cc(:,:,:)

    integer :: nOrb, nAtom, iAt
    logical :: updwn
    real(dp), allocatable :: stimc(:,:,:)
    real(dp), allocatable :: qst(:), R_ASqr(:)

    nOrb = size(cc, dim=1)
    nAtom = size(X_A, dim=2)
    allocate(stimc(norb, norb, 2))
    allocate(qst(nAtom))
    allocate(R_ASqr(nAtom))

    call symm(stimc(:,:,sigma), "L", S, cc(:,:,sigma))
    updwn = (sigma == 1)

    qst(:) = transq(ss, tt, ind, updwn, stimc, cc)
    R_ASqr(:) = sum(X_A**2, dim=1)
!    d2_st = matmul(R_ASqr,qst)
    d2_st = 0.0_dp
    do iAt = 1, nAtom
      d2_st = d2_st + R_ASqr(iAt)*qst(iAt)
    end do

  end subroutine getTransDipoleSqr

  !> Get the localized occupation matrix \lambda_{pq,sigma} = \sum_s n_{s,sigma} U_{ps,sigma}U_{qs,sigma}
  subroutine getlambda(lambda, pp, sigma, qq, uu, occNr)
    real(dp), intent(out)          :: lambda
    integer,  intent(in)           :: pp, qq
    integer,  intent(in)           :: sigma
    real(dp), intent(in)           :: uu(:,:,:)
    real, intent(in)               :: occNr(:,:)

    integer :: ss, nOrb

    nOrb = size(uu, dim=2)
    lambda = 0.0_dp

    do ss = 1, nOrb
      lambda = lambda + uu(pp,ss,sigma) * uu(qq,ss,sigma) * occNr(ss,sigma)
    end do

  end subroutine getlambda

  !> Returns the LO charges Q_{p,\sigma}(A) = \sum_{st} U_{ps,\sigma}U_{pt,\sigma} q_A^{st,\sigma}
  subroutine getLOcharge(Q, pp, sigma, uu, ind, S, cc)
    real(dp), intent(out)          :: Q(:)
    integer,  intent(in)           :: pp
    integer,  intent(in)           :: sigma
    real(dp), intent(in)           :: uu(:,:,:)
    integer,  intent(in)           :: ind(:)
    real(dp), intent(in)           :: S(:,:), cc(:,:,:)

    integer :: ss, tt, nOrb, nAtom, iAt
    logical :: updwn
    real(dp), allocatable :: stimc(:,:,:)
    real(dp), allocatable :: qst(:)

    nOrb = size(uu, dim=2)
    nAtom = size(Q)
    allocate(stimc(norb, norb, 2))
    allocate(qst(nAtom))

    Q = 0.0_dp

    call symm(stimc(:,:,sigma), "L", S, cc(:,:,sigma))
    updwn = (sigma == 1)

    do ss = 1, nOrb
      do tt = 1, nOrb
        qst(:) = transq(ss, tt, ind, updwn, stimc, cc)
        Q(:) = Q(:) + uu(pp,ss,sigma) * uu(pp,tt,sigma) * qst(:)

        print *, "s,t", ss, tt
        print *, "q_st", qst(:)
      end do
    end do

  end subroutine getLOcharge


  subroutine getlambdaQQ(lambdaQQ, uu, ind, S, cc,&
      & occNr, nNeighbor, iNeighbor, img2CentCell)
    real(dp), intent(out)          :: lambdaQQ(:,:)
    real(dp), intent(in)           :: uu(:,:,:)
    integer,  intent(in)           :: ind(:)
    real(dp), intent(in)           :: S(:,:), cc(:,:,:)
    real,     intent(in)           :: occNr(:,:)
    integer,  intent(in)           :: nNeighbor(:)
    integer,  intent(in)           :: iNeighbor(0:,:)
    integer,  intent(in)           :: img2CentCell(:)

    integer :: pp, nOrb, nAtom, nSpin
    integer :: iAt1, iAt2, iAt2f, iNeigh, iSpin
    real(dp), allocatable :: Q(:)
    real(dp):: lambda

    lambdaQQ = 0.0_dp
    nAtom = size(lambdaQQ, dim=1)
    nOrb  = size(occNr, dim=1)
    nSpin = size(occNr, dim=2)

    allocate(Q(nAtom))

    do iAt1 = 1, nAtom
      do iNeigh = 0, nNeighbor(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        do iSpin = 1, nSpin
          do pp = 1, nOrb

            call getLOcharge(Q,pp,iSpin,uu,ind,S,cc)
            call getlambda(lambda,pp,pp,iSpin,uu,occNr)

            ! For debugging only
            if (iAt1 == 1) then
              print *, "-------------"
              print *, "orb", pp
              print *, "lambda", lambda
              print *, "Q", Q
            end if

            lambdaQQ(iAt1,iAt2f) = lambdaQQ(iAt1,iAt2f) +&
                & lambda * Q(iAt1) * Q(iAt2f)

            if (iAt1 /= iAt2f) then
              lambdaQQ(iAt2f,iAt1) = lambdaQQ(iAt2f,iAt1) +&
                  & lambda * Q(iAt1) * Q(iAt2f)
            end if

          end do
        end do

      end do
    end do

  end subroutine getlambdaQQ


  subroutine getlambdaSqrQQ(lambdaSqrQQ, uu, ind, S, cc,&
      & occNr, nNeighbor, iNeighbor, img2CentCell)
    real(dp), intent(out)          :: lambdaSqrQQ(:,:)
    real(dp), intent(in)           :: uu(:,:,:)
    integer,  intent(in)           :: ind(:)
    real(dp), intent(in)           :: S(:,:), cc(:,:,:)
    real, intent(in)               :: occNr(:,:)
    integer,  intent(in)           :: nNeighbor(:)
    integer,  intent(in)           :: iNeighbor(0:,:)
    integer,  intent(in)           :: img2CentCell(:)

    integer :: pp, qq, nOrb, nAtom, nSpin
    integer :: iAt1, iAt2, iAt2f, iNeigh, iSpin
    real(dp), allocatable :: Q_p(:), Q_q(:)
    real(dp):: lambda

    lambdaSqrQQ = 0.0_dp
    nAtom = size(lambdaSqrQQ, dim=1)
    nOrb  = size(occNr, dim=1)
    nSpin = size(occNr, dim=2)

    allocate(Q_p(nAtom))
    allocate(Q_q(nAtom))

    do iAt1 = 1, nAtom
      do iNeigh = 0, nNeighbor(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        do iSpin = 1, nSpin
          do pp = 1, nOrb
            do qq = 1, nOrb

              call getLOcharge(Q_p,pp,iSpin,uu,ind,S,cc)
              call getLOcharge(Q_q,qq,iSpin,uu,ind,S,cc)

              call getlambda(lambda,pp,qq,iSpin,uu,occNr)

              lambdaSqrQQ(iAt1,iAt2f) = lambdaSqrQQ(iAt1,iAt2f) +&
                  & lambda**2 * Q_p(iAt1) * Q_q(iAt2f)

              if (iAt1 /= iAt2f) then
                lambdaSqrQQ(iAt2f,iAt1) = lambdaSqrQQ(iAt2f,iAt1) +&
                  & lambda**2 * Q_p(iAt1) * Q_q(iAt2f)
              end if

            end do
          end do
        end do

      end do
    end do

  end subroutine getlambdaSqrQQ


end module dftbp_loscorrection
