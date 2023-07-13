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

    std::map<std::string, std::size_t> simulate(std::size_t shots) override;


    template<class ReturnType = dd::ComplexValue>
    [[nodiscard]] std::vector<ReturnType> getVectorFromHybridSimulation() const {
        if (CircuitSimulator<Config>::getNumberOfQubits() >= 60) {
            // On 64bit system the vector can hold up to (2^60)-1 elements, if memory permits
            throw std::range_error("getVector only supports less than 60 qubits.");
        }
        if constexpr (std::is_same_v<ReturnType, decltype(finalAmplitudes)>) {
            return finalAmplitudes;
        }
        std::vector<ReturnType> amplitudes;
        std::transform(finalAmplitudes.begin(), finalAmplitudes.end(), std::back_inserter(amplitudes), [](std::complex<dd::fp> x) -> ReturnType { return {x.real(), x.imag()}; });
        return amplitudes;
    }

private:
    std::vector<std::complex<dd::fp>> finalAmplitudes{};


    qc::VectorDD simulateSlicing(size_t nqubits, std::size_t controls);

    class Slice {
    protected:
        qc::Qubit nextControlIdx = 0;

        std::size_t getNextControl() {
            std::size_t idx = 1UL << nextControlIdx;
            nextControlIdx++;
            return controls & idx;
        }

    public:
        qc::Qubit    start;
        qc::Qubit    end;
        std::size_t  controls;
        qc::Qubit    nqubits;
        qc::VectorDD edge{};

        explicit Slice(std::unique_ptr<dd::Package<Config>>& dd, const qc::Qubit start_, const qc::Qubit end_, const std::size_t controls_):
                start(start_), end(end_), controls(controls_), nqubits(end - start + 1), edge(dd->makeZeroState(static_cast<dd::QubitCount>(nqubits), start_)) {
            dd->incRef(edge);
        }

        explicit Slice(std::unique_ptr<dd::Package<Config>>& dd, qc::VectorDD edge_, const qc::Qubit start_, const qc::Qubit end_, const std::size_t controls_):
                start(start_), end(end_), controls(controls_), nqubits(end - start + 1), edge(edge_) {
            dd->incRef(edge);
        }

        Slice combine(std::unique_ptr<dd::Package<Config>>& dd, std::vector<Slice> mp, Slice lower) {
            if ((lower.start <= end && lower.start >= start) || (lower.end <= end && lower.end >= start)) {
                return *this;
            }
                dd->kronecker(edge, lower.edge);
                // update data structures
                for (qc::Qubit i = lower.start; i <= lower.end; ++i) {
                    mp[i] = *this;
                }
                size_t startMin;
                size_t endMax;
                if (lower.start < start) {
                    startMin = lower.start;
                    endMax = end;
                } else {
                    startMin = start;
                    endMax = lower.end;
                }
                // TODO: improve this part. For now we combine all the qubits between two slices to avoid problems.
                for (qc::Qubit i = startMin; i < endMax; ++i) {
                    if (&mp[i] != this && &mp[i] != this) {
                        dd->kronecker(edge, mp[i].edge);
                        for (qc::Qubit j = mp[i].start; j < mp[i].end; ++j) {
                            // assign old qubits to the new slice in order not to use kronecker more than once for same slice
                            mp[j] = *this;
                        }
                        mp[i] = *this;
                    }
                }
                return Slice(dd, startMin, endMax, controls);
        }

        // returns true if this operation was a split operation
        bool apply(std::unique_ptr<dd::Package<Config>>& sliceDD, const std::unique_ptr<qc::Operation>& op);
    };
};

#endif //ENTANGLEMENTS_SPLIT_SIMULATOR_HPP
