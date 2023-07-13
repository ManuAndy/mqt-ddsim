#include "EntanglementSplitSimulator.hpp"

#include <cassert>
#include <cmath>
#include <taskflow/taskflow.hpp>


template<class Config>
qc::VectorDD EntanglementSplitSimulator<Config>::simulateSlicing(size_t nqubits, size_t controls) {
    std::vector<Slice> sliceVector;
    //map qubit -> pointer of slice
    std::vector<Slice> mp;
    std::unique_ptr<dd::Package<Config>> sliceDD;
    sliceDD = std::make_unique<dd::Package<Config>>(CircuitSimulator<Config>::getNumberOfQubits());
    for (qc::Qubit i = 0; i < nqubits; ++i) {
        // initialize map and vectors
        mp.push_back(Slice(sliceDD, i, i, controls));
    }
    // iterate over operations and apply them to corresponding qubits
    // map and vectors can be used inside apply function as well, so deciding
    // combining slices there shouldn't be a problem
    for (const auto& op: *CircuitSimulator<Config>::qc) {
        qc::Qubit firstQubit = 0;
        if (op->isUnitary()) {
            std::set<qc::Qubit> usedQubits = op->getUsedQubits();
            bool atLeastOneTarget = false;
            for (const auto& qubit: usedQubits) {
                if (!atLeastOneTarget) {
                    firstQubit = qubit;
                    atLeastOneTarget = true;
                }
                else {
                    mp[firstQubit] = mp[firstQubit].combine(sliceDD, mp, mp[qubit]);
                }
            }
            if (atLeastOneTarget) {
                for (const auto& qubit: usedQubits) {
                    mp[qubit] = mp[firstQubit];
                }
            }
            [[maybe_unused]] auto x = mp[firstQubit].apply(sliceDD, op);
        }
    }
    //sliceDD->garbageCollect();
    // finally combine all the dds
    auto result = mp[0];
    size_t size = mp.size();
    for (size_t i = 1; i < size; ++i) {
        if (&mp[i] != mp.data()) {
            result = mp[0].combine(sliceDD, mp, mp[i]);
        }
        sliceDD->incRef(result.edge);
    }
    return result.edge;
}

template<class Config>
bool EntanglementSplitSimulator<Config>::Slice::apply(std::unique_ptr<dd::Package<Config>>& sliceDD, const std::unique_ptr<qc::Operation>& op) {
    bool isSplitOp = false;
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    if (reinterpret_cast<qc::StandardOperation*>(op.get()) != nullptr) { // TODO change control and target if wrong direction
        qc::Targets  opTargets{};
        qc::Controls opControls{};

        // check targets
        bool targetInSplit      = false;
        size_t targetIdx = 0;
        bool targetInOtherSplit = false;
        for (const auto& target: op->getTargets()) {
            if (start <= target && target <= end) {
                opTargets.push_back(target);
                targetInSplit = true;
                targetIdx = target;
            } else {
                // merge slices into 1
                if (targetInSplit) {

                }
                targetInOtherSplit = true;
            }
        }

        if (targetInSplit && targetInOtherSplit && !op->getControls().empty()) {
            throw std::invalid_argument("Multiple Targets that are in different slices are not supported at the moment");
        }

        // check controls
        for (const auto& control: op->getControls()) {
            if (start <= control.qubit && control.qubit <= end) {
                opControls.emplace(qc::Control{control.qubit, control.type});
            } else { // other controls are set to the corresponding value
                if (targetInSplit) {
                    isSplitOp              = true;
                    const bool nextControl = getNextControl();
                    if ((control.type == qc::Control::Type::Pos && !nextControl) || // break if control is not activated
                        (control.type == qc::Control::Type::Neg && nextControl)) {
                        return true;
                    }
                }
            }
        }

        if (targetInOtherSplit && !opControls.empty()) { // control slice for split
            if (opControls.size() > 1) {
                throw std::invalid_argument("Multiple controls in control slice of operation are not supported at the moment");
            }

            isSplitOp          = true;
            const bool control = getNextControl();
            for (const auto& c: opControls) {
                sliceDD->decRef(edge); // TODO incref and decref could be integrated in delete edge
                edge = sliceDD->deleteEdge(edge, static_cast<dd::Qubit>(c.qubit), control ? (c.type == qc::Control::Type::Pos ? 0 : 1) : (c.type == qc::Control::Type::Pos ? 1 : 0));
                sliceDD->incRef(edge);
            }
        } else if (targetInSplit) { // target slice for split or operation in split
            const auto&           param = op->getParameter();
            qc::StandardOperation newOp(nqubits, opControls, opTargets, op->getType(), param, start);
            sliceDD->decRef(edge);
            edge = sliceDD->multiply(dd::getDD(&newOp, sliceDD), edge, static_cast<dd::Qubit>(start));
            sliceDD->incRef(edge);
        }
    } else {
        throw std::invalid_argument("Only StandardOperations are supported for now.");
    }
    return isSplitOp;
}

template<class Config>
std::map<std::string, std::size_t> EntanglementSplitSimulator<Config>::simulate(std::size_t shots) {
    auto nqubits    = CircuitSimulator<Config>::getNumberOfQubits();
    auto                                 result       = simulateSlicing(nqubits, 0);
    return Simulator<Config>::measureAllNonCollapsing(shots, result);
}

template class EntanglementSplitSimulator<dd::DDPackageConfig>;
