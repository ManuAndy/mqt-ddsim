#ifndef ENTANGLEMENTS_SPLIT_SIMULATOR_HPP
#define ENTANGLEMENTS_SPLIT_SIMULATOR_HPP

#include "CircuitOptimizer.hpp"
#include "CircuitSimulator.hpp"
#include "QuantumComputation.hpp"
#include "dd/Export.hpp"
#include "dd/Operations.hpp"
#include "dd/Package.hpp"

#include <complex>
#include <memory>

template<class Config = dd::DDPackageConfig>
class EntanglementSplitSimulator: public CircuitSimulator<Config> {
public:

    EntanglementSplitSimulator(std::unique_ptr<qc::QuantumComputation>&& qc_,
                         const ApproximationInfo&                  approxInfo_):
            CircuitSimulator<Config>(std::move(qc_), approxInfo_) {
        // remove final measurements
        qc::CircuitOptimizer::removeFinalMeasurements(*(CircuitSimulator<Config>::qc));
    }

    explicit EntanglementSplitSimulator(std::unique_ptr<qc::QuantumComputation>&& qc_):
            EntanglementSplitSimulator(std::move(qc_), {}) {}

    EntanglementSplitSimulator(std::unique_ptr<qc::QuantumComputation>&& qc_,
                         const ApproximationInfo&                  approxInfo_,
                         const std::uint64_t                       seed_):
            CircuitSimulator<Config>(std::move(qc_), approxInfo_, seed_) {
        // remove final measurements
        qc::CircuitOptimizer::removeFinalMeasurements(*(CircuitSimulator<Config>::qc));
    }

    /** [copied from Simulator.hpp]
     * Run the simulation in the (derived) class.
     * @param shots number of shots to take from the final quantum state
     * @return a map from the strings representing basis states to the number of times they have been measured
     */
    std::map<std::string, std::size_t> simulate(const std::size_t shots) override;

private:

    class Slice {

    public:
        qc::Qubit    start;
        qc::Qubit    end;
        qc::Qubit    nqubits;
        qc::VectorDD edge{};

        explicit Slice(std::unique_ptr<dd::Package<Config>>& dd, const qc::Qubit start_, const qc::Qubit end_):
                start(start_), end(end_), nqubits(end - start + 1), edge(dd->makeZeroState(static_cast<dd::QubitCount>(nqubits), start)) {
            dd->incRef(edge);
        }

        explicit Slice(std::unique_ptr<dd::Package<Config>>& dd, qc::VectorDD edge_, const qc::Qubit start_, const qc::Qubit end_):
                start(start_), end(end_), nqubits(end - start + 1), edge(edge_) {
            dd->incRef(edge);
        }

        static Slice combine(std::unique_ptr<dd::Package<Config>>& dd, Slice& upper, Slice& lower) {
            if (upper.end + 1 != lower.start) {
                throw std::invalid_argument("Only adjacent slices can be combined for now.");
            }
            auto edge = dd->kronecker(lower.edge, upper.edge, false);

            return Slice(dd, edge, upper.start, lower.end);
        }

        void apply(std::unique_ptr<dd::Package<Config>>& sliceDD, const std::unique_ptr<qc::Operation>& op) {
            const auto&       param = op->getParameter();
            const auto&  opControls = op->getControls();
            const auto&   opTargets = op->getTargets();
            qc::StandardOperation newOp(nqubits, opControls, opTargets, op->getType(), param, start);
            sliceDD->decRef(edge);
            edge = sliceDD->multiply(dd::getDD(&newOp, sliceDD), edge, static_cast<dd::Qubit>(start));
            sliceDD->incRef(edge);
        }

        friend bool operator < (Slice const& a, Slice const& b) {
            return (a.start < b.start) || (a.end < b.end);
        }

        friend bool operator == (Slice const& a, Slice const& b) {
            return (a.start == b.start) || (a.end == b.end);
        }

        friend bool operator != (Slice const& a, Slice const& b) {
            return ! (a == b);
        }
    };

    // The uniontable maps every qubit to its slice where it is contained in (second entry in the tuple); the first entry
    // denotes its relative position wihtin the slice
    std::vector<Slice> uniontable{};

    /**
     * @brief Combines the slices of two different qubits x and y and updates all related datastructures
     *
     * @param dd
     * @param x
     * @param y
     * @return true if two slices where combined
     * @return false if qubits were already in the same slice
     */
    bool combine(std::unique_ptr<dd::Package<Config>>& dd, const qc::Qubit x, const qc::Qubit y);

    qc::VectorDD simulateSlicing(const size_t nqubits);

    void apply(std::unique_ptr<dd::Package<Config>>& sliceDD, const std::unique_ptr<qc::Operation>& op);
};

#endif //ENTANGLEMENTS_SPLIT_SIMULATOR_HPP
