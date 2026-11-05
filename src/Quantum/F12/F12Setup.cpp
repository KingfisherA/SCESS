#include <Quantum/Q_PCH.h>
#include "F12.h"
//
#include "Dahlia/Tools/Threading.h"
//
#include "Quantum/HF/HF.h"
#include "Quantum/Molecule.h"
#include "Quantum/Integrals/IntegralTransforms.h"
#include "Quantum/Integrals/Integrals.h"
#include "Quantum/Integrals/IntegralsF12.h"
#include "Tensor/TensorEinsum.h"
#include "Tensor/TensorOperations.h"
#include "Tensor/TensorSlice.h"

//
//
namespace SCESS::F12 {
    //
    //
    using namespace Eigen;
    using Integrals::F12::AOSpace;
    using Integrals::F12::Operator2e;
    //
    //
    namespace {
        //chemists' (p q | r s) -> physicists' <p r | q s>
        TEITensor ChemistsToPhysicists(const TEITensor &chemists) {
            enum { p, q, r, s };
            return Evaluate(Einsum<Indices<p, q, r, s>, Indices<p, r, q, s> >(chemists));
        }
    }
    //
    //
    F12Calculation PrepareF12Calculation(const std::string &obsName,
                                         const std::string &cabsName,
                                         const Molecule &mol,
                                         const DH::Thread *thread,
                                         std::stringstream *console,
                                         bool buildQMCIntegrals) {
        if (mol.GetNoElectrons() % 2 != 0)
            throw std::runtime_error("F12 methods require a closed-shell molecule (an even number of electrons)");

        auto ShouldStop = [thread]() { return thread && thread->ShouldStop(); };
        auto atoms = mol.GetAtomsLibInt();

        F12Calculation calc;
        calc.basisSet = F12BasisSet(obsName, cabsName, atoms);
        const F12BasisSet &basis = calc.basisSet;
        const Index nOBS = (Index) basis.GetNoHFFuncs();
        const Index nCABS = (Index) basis.GetNoCABSFuncs();
        const Index N = nOBS + nCABS;

        if (console)
            *console << "--- F12 basis sets ---\n"
                    << "  OBS  = " << obsName << " (" << nOBS << " functions)\n"
                    << "  CABS = " << cabsName << " (" << nCABS << " functions)\n"
                    << "Evaluating overlap, kinetic, nuclear & OBS eri integrals...\n";

        //-- AO integrals for the reference --
        const MatrixXd S = Integrals::F12::CalculateOverlapMatrix(basis);
        const MatrixXd T = Integrals::F12::CalculateKineticMatrix(basis);
        const MatrixXd V = Integrals::F12::CalculateNuclearMatrix(basis, atoms);

        Integrals::MolInts obsInts;
        obsInts.S = S.topLeftCorner(nOBS, nOBS);
        obsInts.T = T.topLeftCorner(nOBS, nOBS);
        obsInts.V = V.topLeftCorner(nOBS, nOBS);
        obsInts.eri = Integrals::F12::OBSTwoElectronRepulsionTensor(basis);
        if (ShouldStop())
            return calc;

        if (console)
            *console << "Obtaining the closed-shell reference in the OBS...\n";
        //Convergence on the change in the density matrix. Machine precision is not reachable for larger basis
        //sets, and 1e-10 leaves errors in the orbitals far below those of the F12 approximations
        constexpr double SCF_CONVERGENCE = 1e-10;
        const uint32_t electronPairs = mol.GetNoElectrons() / 2;
        calc.hfWfn = HF::UHF(obsInts,
                             electronPairs,
                             electronPairs,
                             mol.CalculateVnn(),
                             SCF_CONVERGENCE,
                             100,
                             thread,
                             console);
        if (ShouldStop())
            return calc;

        calc.idxs = GenerateF12Idxs(basis, calc.hfWfn);
        const Index nOcc = calc.idxs.nOcc;
        F12MOIntegrals &ints = calc.ints;

        //Occupied orbitals over the OBS + CABS functions (they have no CABS components)
        MatrixXd C_occ = MatrixXd::Zero(N, nOcc);
        C_occ.topRows(nOBS) = calc.hfWfn.C_Alpha.leftCols(nOcc);

        //-- Fock matrix over the OBS + CABS functions --
        if (console)
            *console << "Building the Fock matrix over the OBS + CABS...\n";
        //(i nu | j sig) for the occupied i, j. The exchange matrix follows from it and once the CABS is known it
        //is completed into the <ij|pq> integrals
        TEITensor halfERI = Integrals::F12::HalfTransformTwoElectron(basis, Operator2e::Coulomb,
                                                                    C_occ, AOSpace::All, C_occ, AOSpace::All);
        MatrixXd K_AO = MatrixXd::Zero(N, N);
        for (Index nu = 0; nu < N; nu++)
            for (Index sig = 0; sig < N; sig++)
                for (Index i = 0; i < nOcc; i++)
                    K_AO(nu, sig) += halfERI(i, nu, i, sig);

        const MatrixXd D = C_occ.topRows(nOBS) * C_occ.topRows(nOBS).transpose();
        const MatrixXd J_AO = Integrals::F12::CoulombMatrix(basis, D);
        const MatrixXd F_AO = T + V + 2.0 * J_AO - K_AO;
        if (ShouldStop())
            return calc;

        //-- Canonical CABS --
        if (console)
            *console << "Constructing and canonicalising the CABS...\n";
        const MatrixXd Q = Integrals::Transforms::BuildCABSCoefficients(S, F_AO, calc.hfWfn.C_Alpha, nOBS, nCABS);

        ints.h = Q.transpose() * (T + V) * Q;
        ints.F = Q.transpose() * F_AO * Q;
        ints.K = Q.transpose() * K_AO * Q;

        //-- Integrals with an occupied pair, <ij|pq> --
        if (console)
            *console << "Evaluating the eri, f12, f12^2, f12/r12 & |del f12|^2 integrals with an occupied pair...\n";
        ints.eri = ChemistsToPhysicists(Integrals::F12::FinishTransform(halfERI, Q, Q));
        halfERI = TEITensor();

        auto PairIntegrals = [&](Operator2e op) {
            return ChemistsToPhysicists(Integrals::F12::FinishTransform(
                Integrals::F12::HalfTransformTwoElectron(basis, op, C_occ, AOSpace::All, C_occ, AOSpace::All), Q, Q));
        };
        ints.f12 = PairIntegrals(Operator2e::F12);
        if (ShouldStop())
            return calc;
        ints.f12Sq = PairIntegrals(Operator2e::F12Squared);
        if (ShouldStop())
            return calc;
        ints.f12DelSq = PairIntegrals(Operator2e::DelF12Squared);
        if (ShouldStop())
            return calc;
        ints.f12Repulsion = PairIntegrals(Operator2e::F12Coulomb);
        if (ShouldStop())
            return calc;

        //-- OBS-only ERIs --
        ints.eriOBS = Integrals::Transforms::TransformTEI(obsInts.eri, calc.hfWfn.C_Alpha, true);
        obsInts.eri = TEITensor();

        //-- <mx|pq> for CCSD-QMC[MP2-F12] --
        if (buildQMCIntegrals) {
            if (console)
                *console << "Evaluating the <mx|pq> integrals...\n";
            //(m p | x q) with p, q over the OBS
            const MatrixXd Q_CABS = Q.rightCols(nCABS);
            const MatrixXd &C_OBS = calc.hfWfn.C_Alpha;
            ints.eri_oxpq = ChemistsToPhysicists(Integrals::F12::FinishTransform(
                Integrals::F12::HalfTransformTwoElectron(basis, Operator2e::Coulomb,
                                                         C_occ, AOSpace::OBS, Q_CABS, AOSpace::OBS),
                C_OBS, C_OBS));
        }

        //Corrects the basis set error of the reference energy by singly exciting into the CABS virtuals
        calc.E_CABSSingles = EvaluateCABSSinglesCorrection(MatrixToTensor(ints.F), calc.idxs);
        if (console)
            *console << std::setprecision(12)
                    << "E(HF) = " << calc.hfWfn.E_HF << " hartrees\n"
                    << "E(CABS singles correction) = " << calc.E_CABSSingles << " hartrees\n\n";

        return calc;
    }
    //
    F12MOIntegrals ExtractF12MOIntegrals(const Integrals::MolInts &moInts,
                                         const Integrals::F12::F12Integrals &moF12Ints,
                                         const MatrixXd &F,
                                         F12Idxs f12Idxs) {
        enum { I, J, P, Q, R, S, M, X };
        const Index nOcc = f12Idxs.nOcc, nOBS = f12Idxs.nOBS, nCABS = f12Idxs.nCABS, N = nOBS + nCABS;
        TensorSliceIndices idxs = {
            {I, {0, nOcc}}, {J, {0, nOcc}}, {M, {0, nOcc}},
            {P, {0, N}}, {Q, {0, N}},
            {R, {0, nOBS}}, {S, {0, nOBS}},
            {X, {nOBS, nCABS}}
        };

        F12MOIntegrals ints;
        ints.h = moInts.T + moInts.V;
        ints.F = F;
        ints.K = MatrixXd::Zero(N, N);
        for (Index p = 0; p < N; p++)
            for (Index q = 0; q < N; q++)
                for (Index i = 0; i < nOcc; i++)
                    ints.K(p, q) += moInts.eri(i, i, p, q);

        ints.eriOBS = Slice<Indices<R, S, R, S> >(idxs, moInts.eri);
        ints.eri = Slice<Indices<I, J, P, Q> >(idxs, moInts.eri);
        ints.f12 = Slice<Indices<I, J, P, Q> >(idxs, moF12Ints.f12);
        ints.f12Sq = Slice<Indices<I, J, P, Q> >(idxs, moF12Ints.f12Sq);
        ints.f12DelSq = Slice<Indices<I, J, P, Q> >(idxs, moF12Ints.f12DelSq);
        ints.f12Repulsion = Slice<Indices<I, J, P, Q> >(idxs, moF12Ints.f12Repulsion);
        ints.eri_oxpq = Slice<Indices<M, X, R, S> >(idxs, moInts.eri);
        return ints;
    }
    //
    //
}
