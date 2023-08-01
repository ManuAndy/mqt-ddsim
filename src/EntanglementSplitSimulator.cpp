#include "EntanglementSplitSimulator.hpp"

#include <cassert>
#include <cmath>

template<class Config>
bool EntanglementSplitSimulator<Config>::combine(std::unique_ptr<dd::Package<Config>>& dd, const qc::Qubit x, const qc::Qubit y) {
    Slice lower = uniontable.at(std::max(x, y));
    Slice upper = uniontable.at(std::min(x, y));

    if (lower != upper) {
        Slice combined = Slice::combine(dd, upper, lower);

        // update all affected qubits in the uniontable
        for (auto i = combined.start; i <= combined.end; i++) {
            uniontable.at(i) = combined;
        }

        return true;
    } else {
        return false;
    }
}

template<class Config>
std::map<std::string, std::size_t> EntanglementSplitSimulator<Config>::simulateSlicing(const std::size_t shots, const size_t nqubits) {
    auto sliceDD = std::make_unique<dd::Package<Config>>(nqubits);

    for (qc::Qubit i = 0; i < nqubits; ++i) {
        // initialize map and vectors
        auto slice = Slice(sliceDD, i, i);
        uniontable.push_back(slice);
    }
    // iterate over operations and apply them to corresponding qubits
    for (const auto& op: *CircuitSimulator<Config>::qc) {
        apply(sliceDD, op);
        sliceDD->garbageCollect();
    }
    std::cout << "length " << uniontable.at(0).nqubits << "\n";
    // finally combine all the dds
    for (qc::Qubit i = 0; i < nqubits; i ++) {
        combine(sliceDD, 0, i);
    }
    auto edge = uniontable.at(0).edge;
    std::cout << "length " << uniontable.at(0).nqubits << "\n";
    auto result = Simulator<Config>::measureAllNonCollapsing(shots, edge);
    return result;
}

template<class Config>
void EntanglementSplitSimulator<Config>::apply(std::unique_ptr<dd::Package<Config>>& sliceDD, const std::unique_ptr<qc::Operation>& op) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    if (reinterpret_cast<qc::StandardOperation*>(op.get()) != nullptr) { // TODO change control and target if wrong direction
        qc::Qubit firstTarget = op->getTargets().front();

        // combine targets
        for (const auto& target: op->getTargets()) {
            combine(sliceDD, target, firstTarget); // TODO omit the call where target == firstTarget
        }
        // combine controls
        for (const auto& control: op->getControls()) {
            combine(sliceDD, control.qubit, firstTarget);
        }

        // retrieve combined slice and its edge
        // TODO this could be done more elegantly
        Slice combined = uniontable.at(firstTarget);

        // apply gate
        combined.apply(sliceDD, op);

        for (auto i = combined.start; i <= combined.end; i++) {
            uniontable.at(i) = combined;
        }
    } else {
        throw std::invalid_argument("Only StandardOperations are supported for now.");
    }
}

template<class Config>
std::map<std::string, std::size_t> EntanglementSplitSimulator<Config>::simulate(const std::size_t shots) {
    auto nqubits                = CircuitSimulator<Config>::getNumberOfQubits();
    auto result = simulateSlicing(shots, nqubits);
    return result;
}

template class EntanglementSplitSimulator<dd::DDPackageConfig>;
