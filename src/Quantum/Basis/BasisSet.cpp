#include <Quantum/Q_PCH.h>
#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/replace.hpp>
//
#include "Resources/AssetFolder.h"
#include <libint2/basis.h>
#include <libint2/config.h>
//
#include "Quantum/Molecule.h"
#include "Quantum/Basis/BasisSet.h"
#include "Quantum/HF/HF.h"
#include "Utilities/Utilities.h"

//
//
namespace SCESS
{
    //
    //
    namespace
    {
        //Display order of the basis set families (folders of Assets/Basis); any other family is listed after these
        const std::vector<std::string> FAMILY_ORDER = {"Pople", "Dunning", "Karlsruhe", "F12", BasisLibrary::CABS_FAMILY};
        //Family of basis sets placed directly in Assets/Basis
        constexpr const char *OTHER_FAMILY = "Other";

        size_t FamilyRank(const std::string &family)
        {
            const auto it = std::ranges::find(FAMILY_ORDER, family);
            return (size_t) std::distance(FAMILY_ORDER.begin(), it);
        }

        //Lower case with the 6-31G* style polarisation shorthand written as the 6-31G(d) style used by the file names
        std::string CanonicalName(const std::string &name)
        {
            std::string canonical = boost::algorithm::to_lower_copy(name);
            boost::algorithm::replace_all(canonical, "**", "(d,p)");
            boost::algorithm::replace_all(canonical, "*", "(d)");
            return canonical;
        }

        std::vector<BasisSetInfo> ScanBasisSets()
        {
            namespace fs = std::filesystem;
            std::vector<BasisSetInfo> basisSets;
            const fs::path basisPath = Res::AssetsFolder::Get().GetResourcePath("Basis").string();

            std::error_code ec;
            for (fs::recursive_directory_iterator itr(basisPath, ec), end; !ec && itr != end; itr.increment(ec))
            {
                if (!itr->is_regular_file() || itr->path().extension() != ".json")
                    continue;

                const fs::path parent = itr->path().parent_path();
                basisSets.push_back({itr->path().stem().string(),
                                     fs::equivalent(parent, basisPath) ? OTHER_FAMILY : parent.filename().string(),
                                     itr->path()});
            }

            std::ranges::sort(basisSets, [](const BasisSetInfo &a, const BasisSetInfo &b)
            {
                const size_t rankA = FamilyRank(a.family), rankB = FamilyRank(b.family);
                if (rankA != rankB)
                    return rankA < rankB;
                if (a.family != b.family)
                    return a.family < b.family;
                return CanonicalName(a.name) < CanonicalName(b.name);
            });

            return basisSets;
        }
    }
    //
    const std::vector<BasisSetInfo> &BasisLibrary::List()
    {
        //The assets folder must have been initialised before the first call
        static const std::vector<BasisSetInfo> basisSets = ScanBasisSets();
        return basisSets;
    }
    //
    std::vector<std::string> BasisLibrary::ListFamilies()
    {
        std::vector<std::string> families;
        for (const BasisSetInfo &info : List())
            if (std::ranges::find(families, info.family) == families.end())
                families.push_back(info.family);
        return families;
    }
    //
    const BasisSetInfo *BasisLibrary::Find(const std::string &name)
    {
        const std::string canonical = CanonicalName(name);
        for (const BasisSetInfo &info : List())
            if (CanonicalName(info.name) == canonical)
                return &info;
        return nullptr;
    }
    //
    const BasisSetInfo &BasisLibrary::Get(const std::string &name)
    {
        const BasisSetInfo *info = Find(name);
        if (!info)
            throw std::runtime_error("Cannot locate the basis set " + name + " in Assets/Basis");
        return *info;
    }
    //
    std::string BasisLibrary::DefaultCABS(const std::string &obsName)
    {
        const BasisSetInfo *obs = Find(obsName);
        if (!obs)
            return "";

        //cc-pVXZ-F12 -> cc-pVXZ-F12-OptRI; cc-pVXZ and aug-cc-pVXZ -> aug-cc-pVXZ-OptRI
        std::string candidate;
        const std::string lower = boost::algorithm::to_lower_copy(obs->name);
        if (lower.starts_with("cc-pv") && lower.ends_with("-f12"))
            candidate = obs->name + "-OptRI";
        else if (lower.starts_with("aug-cc-pv") && lower.size() == std::string("aug-cc-pvxz").size())
            candidate = obs->name + "-OptRI";
        else if (lower.starts_with("cc-pv") && lower.size() == std::string("cc-pvxz").size())
            candidate = "aug-" + obs->name + "-OptRI";

        const BasisSetInfo *cabs = candidate.empty() ? nullptr : Find(candidate);
        return cabs ? cabs->name : "";
    }
    //
    nlohmann::json BasisLibrary::LoadJSON(const std::string &name)
    {
        const BasisSetInfo &info = Get(name);
        std::ifstream jsonFile(info.path);
        if (!jsonFile)
            throw std::runtime_error("Cannot open the basis set file " + info.path.string());
        return nlohmann::json::parse(jsonFile);
    }
    //
    std::vector<libint2::Shell> BasisLibrary::LoadShells(const std::string &name,
                                                         const std::vector<libint2::Atom> &atoms)
    {
        const nlohmann::json basisSet = LoadJSON(name);
        const nlohmann::json &elements = basisSet.at("elements");

        std::vector<libint2::Shell> shells;
        for (const libint2::Atom &atom : atoms)
        {
            const std::string Z = std::to_string(atom.atomic_number);
            if (!elements.contains(Z) || !elements.at(Z).contains("electron_shells"))
                throw std::runtime_error("The basis set " + name + " has no functions for the element Z = " + Z);

            for (const nlohmann::json &shell : elements.at(Z).at("electron_shells"))
            {
                const nlohmann::json &angMomenta = shell.at("angular_momentum");
                const nlohmann::json &exponents = shell.at("exponents");
                const nlohmann::json &coefficients = shell.at("coefficients");
                const bool pure = shell.value("function_type", "gto") != "gto_cartesian";

                //As in BasisSet::LoadAndAddBasisSet, every coefficient array is a separate contracted shell
                if (angMomenta.size() > 1 && angMomenta.size() != coefficients.size())
                    throw std::runtime_error("Basis set shell has mismatched angular momenta and coefficient sets");

                for (size_t j = 0; j < coefficients.size(); j++)
                {
                    const int l = angMomenta[angMomenta.size() > 1 ? j : 0].get<int>();
                    if (l > LIBINT_MAX_AM)
                        throw std::runtime_error("The basis set " + name + " has functions of angular momentum " +
                                                 std::to_string(l) + " but libint2 was built with a maximum of " +
                                                 std::to_string(LIBINT_MAX_AM));

                    libint2::svector<double> alpha, coeffs;
                    for (size_t k = 0; k < coefficients[j].size(); k++)
                    {
                        const double coefficient = std::stod(coefficients[j][k].get<std::string>());
                        if (coefficient == 0.0)
                            continue;
                        alpha.push_back(std::stod(exponents[k].get<std::string>()));
                        coeffs.push_back(coefficient);
                    }

                    shells.emplace_back(alpha, libint2::svector<libint2::Shell::Contraction>{{l, pure, coeffs}},
                                        std::array<double, 3>{atom.x, atom.y, atom.z});
                }
            }
        }

        return shells;
    }
    //
    //
    BasisSet::BasisSet(const std::string &basisSetName, const Molecule &mol)
    {
        LoadAndAddBasisSet(basisSetName, mol);
    }
    //
    std::vector<std::string> BasisSet::ListBasisSets()
    {
        std::vector<std::string> names;
        for (const BasisSetInfo &info : BasisLibrary::List())
            names.push_back(info.name);
        return names;
    }
    //
    void BasisSet::LoadAndAddBasisSet(const std::string &basisSetName, const Molecule &mol)
    {
        const nlohmann::json basisSet = BasisLibrary::LoadJSON(basisSetName);

        std::vector<Eigen::Vector3i> angularMomenta;
        auto GenerateAngularMomenta = [&angularMomenta](uint32_t angMom)
        {
            angularMomenta.clear();
            angularMomenta.reserve((angMom + 1) * (angMom + 2) / 2);

            for (uint32_t x = 0; x <= angMom; x++)
                for (uint32_t y = 0; y <= angMom; y++)
                    for (uint32_t z = 0; z <= angMom; z++)
                    {
                        if (x + y + z == angMom)
                            angularMomenta.emplace_back(x, y, z);
                    }
        };

        for (const Atom &a : mol.GetAtoms())
        {
            const std::string Z = std::to_string(a.GetCharge());
            if (!basisSet.at("elements").contains(Z) || !basisSet.at("elements").at(Z).contains("electron_shells"))
                throw std::runtime_error("The basis set " + basisSetName + " has no functions for the element Z = " + Z);
            const nlohmann::json &shells = basisSet.at("elements").at(Z).at("electron_shells");
            for (size_t i = 0; i < shells.size(); i++)
            {
                const nlohmann::json &shell = shells[i];
                const nlohmann::json &angMomenta = shell["angular_momentum"];
                const nlohmann::json &exponents = shell["exponents"];
                const nlohmann::json &coefficients = shell["coefficients"];

                // Every coefficient array is a separate contracted function. For Pople SP-type shells
                // (several angular momenta) the j-th array belongs to the j-th angular momentum; for
                // general contractions (one angular momentum) every array shares that angular momentum.
                if (angMomenta.size() > 1 && angMomenta.size() != coefficients.size())
                    throw std::runtime_error("Basis set shell has mismatched angular momenta and coefficient sets");

                for (size_t j = 0; j < coefficients.size(); j++)
                {
                    int angularMomentum = angMomenta[angMomenta.size() > 1 ? j : 0].get<int>();
                    GenerateAngularMomenta((uint32_t)angularMomentum);

                    const nlohmann::json &coeffCollection = coefficients[j];
                    for (const Eigen::Vector3i &angMom : angularMomenta)
                    {
                        CGF cgf(a.GetPos());
                        for (size_t k = 0; k < coeffCollection.size(); k++)
                        {
                            double coefficient = std::stod(coeffCollection[k].get<std::string>());
                            if (coefficient == 0.0)
                                continue;
                            double exponent = std::stod(exponents[k].get<std::string>());
                            cgf.Add(angMom, PGF(exponent, coefficient));
                        }
                        mBasisFuncs.push_back(cgf);
                    }
                }
            }
        }
    }
    //
    void BasisSet::AddCGF(const CGF &cgf)
    {
        mBasisFuncs.push_back(cgf);
    }
    //
    size_t BasisSet::GetNoPrimitives() const
    {
        size_t numPrimitives = 0;

        for (const CGF &cgf : mBasisFuncs)
            numPrimitives += cgf.Size();

        return numPrimitives;
    }
    //
    void F12BasisSet::AddBasisSets(const std::string &hfBasisSetName, const std::string &cabsBasisSetName,
                                   const std::vector<libint2::Atom> &molecule)
    {
        mShells.clear();
        mHFShellsStartIdx = 0;
        mShells = BasisLibrary::LoadShells(hfBasisSetName, molecule);

        mCABSShellsSetStartIdx = mShells.size();

        const std::vector<libint2::Shell> cabsShells = BasisLibrary::LoadShells(cabsBasisSetName, molecule);
        mShells.insert(mShells.end(), cabsShells.begin(), cabsShells.end());

        mNumOfFuncs = 0;
        mMaxL = 0;
        mNumPrimitives = 0;
        for (auto &&shell : mShells)
        {
            mNumOfFuncs += shell.size();
            mNumPrimitives = std::max(shell.nprim(), mNumPrimitives);
            for (auto &&c : shell.contr)
                mMaxL = std::max(c.l, mMaxL);
        }
    }
    //
    std::vector<size_t> F12BasisSet::ShellToBasisIdx() const
    {
        std::vector<size_t> result(mShells.size(), 0);

        size_t counter = 0;
        for (size_t i = 0; i < mShells.size(); i++)
        {
            result[i] = counter;
            counter += mShells[i].size();
        }

        return result;
    }
    //
    //
}
