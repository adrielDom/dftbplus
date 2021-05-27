!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

module dftbp_loscorrection
# include "assert.h"
# include "allocate.h"
  use Accuracy
  use CommonTypes
  use Message
  use BlasRoutines
  use SCC, only : getAtomicGammaMatrix
  use LinearResponseCommon, only : transq
  implicit none
  public

  !> Data type for localised orbital scaling correction
  type :: TLOSCorrection

    !> energy penalty function
    character :: penalty

    !> parameters of penalty function
    real(dp), allocatable :: PenaltyParam(:)

    !> SCF-LOSC Calculation?
    logical :: tSCF

  end type TLOSCorrection

contains

  !> Get the LOSC term of total energy
  subroutine getElosc(e_losc, nNeighbor, iNeighbor, img2CentCell, ind,&
      & S, cc, filling, uu)
    real(dp), intent(out) :: e_losc(:)
    integer,  intent(in)  :: nNeighbor(:)
    integer,  intent(in)  :: iNeighbor(0:,:)
    integer,  intent(in)  :: img2CentCell(:)
    integer,  intent(in)  :: ind(:)
    real(dp), intent(in)  :: S(:,:), cc(:,:,:)
    real(dp), intent(in)  :: filling(:,:)
    real(dp), intent(in)  :: uu(:,:)

    real(dp), allocatable :: gamma(:,:)
    real(dp), allocatable :: lambdaQQ(:,:), lambdaSqrQQ(:,:)
    real,     allocatable :: occNr(:)

    integer :: nAtom, nOrb
    integer :: iAt1, iAt2, iAt2f, iNeigh
    real(dp) :: rTmp, term1, term2

    nAtom = size(e_losc)

    ALLOCATE_(gamma, (nAtom, nAtom))
    ALLOCATE_(lambdaQQ, (nAtom, nAtom))
    ALLOCATE_(lambdaSqrQQ, (nAtom, nAtom))

    ALLOCATE_(occNr, (size(filling,dim=1)))
    occNr = real(filling(:,1), 4)

    e_losc = 0.0_dp

    term1 = 0.0_dp
    term2 = 0.0_dp

    call getAtomicGammaMatrix(gamma, iNeighbor, img2CentCell)
    call getlambdaQQ(lambdaQQ,uu,ind,S,cc,occNr,nNeighbor,iNeighbor,img2CentCell)
    call getlambdaSqrQQ(lambdaSqrQQ,uu,ind,S,cc,occNr,nNeighbor,iNeighbor,img2CentCell)


    do iAt1 = 1, nAtom
      do iNeigh = 0, nNeighbor(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)

        rTmp = gamma(iAt1,iAt2) * lambdaQQ(iAt1,iAt2)&
            & - gamma(iAt1,iAt2) * lambdaSqrQQ(iAt1,iAt2)

! term1 and term2 just for debugging purposes
        term1 = term1 + 0.5_dp * gamma(iAt1,iAt2) * lambdaQQ(iAt1,iAt2)
        term2 = term2 + 0.5_dp * gamma(iAt1,iAt2) * lambdaSqrQQ(iAt1,iAt2)

        e_losc(iAt1) = e_losc(iAt1) + 0.5_dp * rTmp

        if (iAt1 /= iAt2f) then
          e_losc(iAt2f) = e_losc(iAt2f) + 0.5_dp * rTmp
        end if

      end do
    end do

    print *, "term1", term1
    print *, "term2", term2

  end subroutine getElosc


  subroutine getOrbitallets(uu, ener, nNeighbor, iNeighbor, img2CentCell, ind,&
      & S, cc, filling, rat)
    real(dp), intent(out) :: uu(:,:)
    real(dp), intent(in) :: ener(:)
    integer, intent(in)   :: nNeighbor(:)
    integer,  intent(in)  :: iNeighbor(0:,:)
    integer,  intent(in)  :: img2CentCell(:)
    integer,  intent(in)  :: ind(:)
    real(dp), intent(in)  :: S(:,:), cc(:,:,:)
    real(dp), intent(in)  :: filling(:,:)
    real(dp), intent(in)  :: rat(:,:)

    logical :: tFconverged
    integer :: pp, ss, qq, rr, ii, iFloop, nFloop
    integer :: nOrb
    real (dp) :: D0, D1, D2, C0, C1
    real (dp) :: Aq, Bq, Ar, Br
    real (dp) :: wqs, wrs
    real (dp) :: fqr, f_der1, f_der2, fqr_tmp, F, F_new
    real (dp) :: phi(20), phi_new, phi_tmp, delta_phi
    real (dp), allocatable :: norm(:), U_tmp(:)

    nOrb = size(uu, dim=2)
    ALLOCATE_(norm, (nOrb))
    ALLOCATE_(U_tmp, (nOrb))
    tFconverged = .false.
    nFloop = 100
    uu = 0.0_dp
    F = 0.0_dp

    !initializing the LO orbitals as the KS orbitals
    do pp = 1, nOrb
      do ss = pp, nOrb
        if (pp == ss) then
          uu(pp,ss) = 1.0_dp
        end if
      end do
    end do

    iFloop = 1
    loopF: do while (iFloop <= nFloop)
      print *, "============================"
      print *, "i =", iFloop
      print *, "============================"
      F_new = 0.0_dp
      do qq = 1, nOrb - 1
        do rr = qq + 1, nOrb

          do ii = 1, 20
            phi(ii) = rand(0)*6.28
          end do

          print *, "q =", qq
          print *, "r =", rr
          print *, "----------------"

          call getLODipole(Aq, qq, rat, ind, S, cc, uu)
          call getLODipoleSqr(Bq, qq, rat, ind, S, cc, uu)
          call getLODipole(Ar, rr, rat, ind, S, cc, uu)
          call getLODipoleSqr(Br, rr, rat, ind, S, cc, uu)

          C0 = Aq - Bq**2 + Ar - Br**2
          C1 = 0.5_dp*( Bq**2 + Br**2 )

          D0 = 0.0_dp
          D1 = 0.0_dp
          D2 = 0.0_dp
          do ss = 1, nOrb

            wqs = penalty(abs(ener(qq)-ener(ss)))
            wrs = penalty(abs(ener(rr)-ener(ss)))

            D0 = D0 + 0.5 * (wqs + wrs) * (uu(qq,ss)**2 + uu(rr,ss)**2)
            D1 = D1 + (wqs - wrs) * uu(qq,ss) * uu(rr,ss)
            D2 = D2 + 0.5 * (wqs - wrs) * (uu(qq,ss)**2 - uu(rr,ss)**2)
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
          U_tmp(:) = uu(qq,:)*dcos(0.5*phi_new) + uu(rr,:)*dsin(0.5*phi_new)
          uu(rr,:) = -1.0_dp*uu(qq,:)*dsin(0.5*phi_new) + uu(rr,:)*dcos(0.5*phi_new)

          uu(qq,:) = U_tmp(:)

          norm(:) = sqrt(sum(uu**2, dim=2))

          print *, "phi:", phi
          print *, "phi_new:", phi_new
          print *, "fqr:", fqr
          print *, "Uq:", uu(qq,:)
          print *, "norm:", norm(qq)
          print *, "Ur:", uu(rr,:)
          print *, "norm:", norm(rr)

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

    function penalty(x) result(w)
      real(dp), intent(in) :: x
      real(dp) :: w

      real(dp) :: ener_frac, tmp
      real(dp) :: R0, eps0, gamma, eta

      ! parameters of the penalty function
      R0 = 5.1_dp !Bohr (2.7 AA)
      eps0 = 0.92_dp !Hartree (2.5 eV)
      gamma = 2.0_dp
      eta = 3.0_dp

      ener_frac = x / eps0
      tmp = 1 - exp( -1.0_dp * ener_frac**eta )

      if (x < eps0) then
        w = tmp * R0**2
      else
        w = tmp * R0**2 * ener_frac**gamma
      end if

    end function penalty

  end subroutine getOrbitallets


  subroutine getLODipole(B, qq, rat, ind, S, cc, uu)
    real(dp), intent(out)          :: B
    integer,  intent(in)           :: qq
    real(dp), intent(in)           :: rat(:,:)
    integer,  intent(in)           :: ind(:)
    real(dp), intent(in)           :: S(:,:), cc(:,:,:)
    real(dp), intent(in)           :: uu(:,:)

    integer :: nOrb, ss, tt
    real(dp) :: d_st

    nOrb = size(cc, dim=1)

    B = 0.0_dp

    do ss = 1, nOrb
      do tt = 1, nOrb

        call getTransDipole(d_st, ss, tt, rat, ind, S, cc)
        B = B + uu(qq,ss)*uu(qq,tt)*d_st
      end do
    end do

  end subroutine getLODipole


  subroutine getLODipoleSqr(A, qq, rat, ind, S, cc, uu)
    real(dp), intent(out)          :: A
    integer,  intent(in)           :: qq
    real(dp), intent(in)           :: rat(:,:)
    integer,  intent(in)           :: ind(:)
    real(dp), intent(in)           :: S(:,:), cc(:,:,:)
    real(dp), intent(in)           :: uu(:,:)

    integer :: nOrb, ss, tt
    real(dp) :: d2_st

    nOrb = size(cc, dim=1)

    A = 0.0_dp

    do ss = 1, nOrb
      do tt = 1, nOrb

        call getTransDipoleSqr(d2_st, ss, tt, rat, ind, S, cc)
        A = A + uu(qq,ss)*uu(qq,tt)*d2_st
      end do
    end do

  end subroutine getLODipoleSqr

  subroutine getTransDipole(d_st, ss, tt, X_A, ind, S, cc)
    real(dp), intent(out)          :: d_st
    integer,  intent(in)           :: ss, tt
    real(dp), intent(in)           :: X_A(:,:)
    integer,  intent(in)           :: ind(:)
    real(dp), intent(in)           :: S(:,:), cc(:,:,:)

    integer :: nOrb, nAtom, iAt
    real(dp), allocatable :: stimc(:,:,:)
    real(dp), allocatable :: qst(:), R_A(:)

    nOrb = size(cc, dim=1)
    nAtom = size(X_A, dim=2)
    ALLOCATE_(stimc, (norb, norb, 2))
    ALLOCATE_(qst, (nAtom))
    ALLOCATE_(R_A, (nAtom))

    call symm(stimc(:,:,1), "L", S, cc(:,:,1), "U")
    call transq(ss, tt, ind, .true., stimc, cc, qst)
    R_A(:) = sqrt(sum(X_A**2, dim=1))

    d_st = 0.0_dp
    do iAt = 1, nAtom
      d_st = d_st + R_A(iAt)*qst(iAt)
    end do
!    d_st = matmul(R_A,qst)

  end subroutine getTransDipole


  subroutine getTransDipoleSqr(d2_st, ss, tt, X_A, ind, S, cc)
    real(dp), intent(out)          :: d2_st
    integer,  intent(in)           :: ss, tt
    real(dp), intent(in)           :: X_A(:,:)
    integer,  intent(in)           :: ind(:)
    real(dp), intent(in)           :: S(:,:), cc(:,:,:)

    integer :: nOrb, nAtom, iAt
    real(dp), allocatable :: stimc(:,:,:)
    real(dp), allocatable :: qst(:), R_ASqr(:)

    nOrb = size(cc, dim=1)
    nAtom = size(X_A, dim=2)
    ALLOCATE_(stimc, (norb, norb, 2))
    ALLOCATE_(qst, (nAtom))
    ALLOCATE_(R_ASqr, (nAtom))

    call symm(stimc(:,:,1), "L", S, cc(:,:,1), "U")
    call transq(ss, tt, ind, .true., stimc, cc, qst)
    R_ASqr(:) = sum(X_A**2, dim=1)
!    d2_st = matmul(R_ASqr,qst)
    d2_st = 0.0_dp
    do iAt = 1, nAtom
      d2_st = d2_st + R_ASqr(iAt)*qst(iAt)
    end do

  end subroutine getTransDipoleSqr

  !> Get the localized occupation matrix \lambda_{pq} = \sum_s n_s U_{ps}U_{qs}
  subroutine getlambda(lambda_pq, pp, qq, uu, occNr)
    real(dp), intent(out)          :: lambda_pq(:)
    integer,  intent(in)           :: pp, qq
    real(dp), intent(in)           :: uu(:,:)
    real, intent(in)               :: occNr(:)

    integer :: ss, nOrb
    real :: occNrSpin(2)

    nOrb = size(uu, dim=2)
    lambda_pq = 0.0_dp

    do ss = 1, nOrb
      ! CHECK THIS! for spin polarization it should exist a second term involving spin constants W
      if (occNr(ss) > 1) then
        occNrSpin(1) = 1.0
        occNrSpin(2) = occNr(ss) - 1.0
      else
        occNrSpin(1) = occNr(ss)
        occNrSpin(2) = 0.0
      end if

      lambda_pq(:) = lambda_pq(:) + &
          & uu(pp,ss)*uu(qq,ss)*occNrSpin(:)
    end do

  end subroutine getlambda

  !> Returns the transition LO charges Q_p(A) = \sum_{st} U_{ps}U_{pt} q_A^{st}
  subroutine transLocQ(Q,pp,uu,ind,S,cc)
    real(dp), intent(out)          :: Q(:)
    integer,  intent(in)           :: pp
    real(dp), intent(in)           :: uu(:,:)
    integer,  intent(in)           :: ind(:)
    real(dp), intent(in)           :: S(:,:), cc(:,:,:)

    integer :: ss, tt, nOrb, nAtom, iAt
    real(dp), allocatable :: stimc(:,:,:)
    real(dp), allocatable :: qst(:)

    nOrb = size(uu, dim=2)
    nAtom = size(Q)
    ALLOCATE_(stimc, (norb, norb, 2))
    ALLOCATE_(qst, (nAtom))

    Q = 0.0_dp

    call symm(stimc(:,:,1), "L", S, cc(:,:,1), "U")

    do iAt = 1, nAtom
      print *, "**************"
      print *, "atom:", iAt

      do ss = 1, nOrb
        do tt = 1, nOrb
          call transq(ss, tt, ind, .true., stimc, cc, qst)
          Q(iAt) = Q(iAt) + uu(pp,ss) * uu(pp,tt) * qst(iAt)

          print *, "s,t", ss, tt
          print *, "q_st", qst(iAt)
        end do
      end do
    end do

  end subroutine transLocQ


  subroutine getlambdaQQ(lambdaQQ, uu, ind, S, cc,&
      & occNr, nNeighbor, iNeighbor, img2CentCell)
    real(dp), intent(out)          :: lambdaQQ(:,:)
    real(dp), intent(in)           :: uu(:,:)
    integer,  intent(in)           :: ind(:)
    real(dp), intent(in)           :: S(:,:), cc(:,:,:)
    real,     intent(in)           :: occNr(:)
    integer,  intent(in)           :: nNeighbor(:)
    integer,  intent(in)           :: iNeighbor(0:,:)
    integer,  intent(in)           :: img2CentCell(:)

    integer :: pp, nOrb, nAtom
    integer :: iAt1, iAt2, iAt2f, iNeigh, iSpin
    real(dp), allocatable :: Q(:)
    real(dp):: lambda(2)

    lambdaQQ = 0.0_dp
    nAtom = size(lambdaQQ, dim=1)
    nOrb  = size(occNr)

    ALLOCATE_(Q, (nAtom))

    do iAt1 = 1, nAtom
      do iNeigh = 0, nNeighbor(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        do iSpin = 1, 2
          do pp = 1, nOrb

            call transLocQ(Q,pp,uu,ind,S,cc)
            call getlambda(lambda,pp,pp,uu,occNr)

            if (iAt1 == 1) then
              print *, "-------------"
              print *, "orb", pp
              print *, "lambda", lambda
              print *, "Q", Q
            end if

            lambdaQQ(iAt1,iAt2) = lambdaQQ(iAt1,iAt2) +&
                & lambda(iSpin) * Q(iAt1) * Q(iAt2)
          end do
        end do

      end do
    end do

  end subroutine getlambdaQQ


  subroutine getlambdaSqrQQ(lambdaSqrQQ, uu, ind, S, cc,&
      & occNr, nNeighbor, iNeighbor, img2CentCell)
    real(dp), intent(out)          :: lambdaSqrQQ(:,:)
    real(dp), intent(in)           :: uu(:,:)
    integer,  intent(in)           :: ind(:)
    real(dp), intent(in)           :: S(:,:), cc(:,:,:)
    real, intent(in)               :: occNr(:)
    integer,  intent(in)           :: nNeighbor(:)
    integer,  intent(in)           :: iNeighbor(0:,:)
    integer,  intent(in)           :: img2CentCell(:)

    integer :: pp, qq, nOrb, nAtom
    integer :: iAt1, iAt2, iAt2f, iNeigh, iSpin
    real(dp), allocatable :: Q_p(:), Q_q(:)
    real(dp):: lambda(2)

    lambdaSqrQQ = 0.0_dp
    nAtom = size(lambdaSqrQQ, dim=1)
    nOrb  = size(occNr)

    ALLOCATE_(Q_p, (nAtom))
    ALLOCATE_(Q_q, (nAtom))

    do iAt1 = 1, nAtom
      do iNeigh = 0, nNeighbor(iAt1)
        iAt2 = iNeighbor(iNeigh, iAt1)
        iAt2f = img2CentCell(iAt2)
        do iSpin = 1, 2
          do pp = 1, nOrb
            do qq = 1, nOrb

              call transLocQ(Q_p,pp,uu,ind,S,cc)
              call transLocQ(Q_q,qq,uu,ind,S,cc)

              call getlambda(lambda,pp,qq,uu,occNr)

              lambdaSqrQQ(iAt1,iAt2) = lambdaSqrQQ(iAt1,iAt2) +&
                  & lambda(iSpin)**2 * Q_p(iAt1) * Q_q(iAt2)

            end do
          end do
        end do

      end do
    end do

  end subroutine getlambdaSqrQQ


end module dftbp_loscorrection
