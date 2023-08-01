#include "EntanglementSplitSimulator.hpp"
#include "CircuitSimulator.hpp"


#include <gtest/gtest.h>
#include <memory>

enum {
NumberOfTimeMeasurements = 1
};

template <
        class result_t   = std::chrono::milliseconds,
        class clock_t    = std::chrono::steady_clock,
        class duration_t = std::chrono::milliseconds
        >
auto since(std::chrono::time_point<clock_t, duration_t> const& start)
{
    return std::chrono::duration_cast<result_t>(clock_t::now() - start);
}

using namespace qc::literals;
/*TEST(EntanglementSplitSimulatorDDSIM, TrivialTest) {
    auto quantumComputation = [] {
        auto qc = std::make_unique<qc::QuantumComputation>(2);
        qc->h(0);         // 1| h q[0];
        qc->x(1, {0_pc}); // 2| cx q[0], q[1];
        // qc->i(1); // some dummy operation
        return qc;
    };

    EntanglementSplitSimulator ddsim(quantumComputation());

    auto resultAmp = ddsim.simulate(128);
    for (const auto& entry: resultAmp) {
        std::cout << "resultAmp[" << entry.first << "] = " << entry.second << "\n";
    }

    ASSERT_EQ(resultAmp.size(), 2);
    auto it = resultAmp.find("00");
    ASSERT_TRUE(it != resultAmp.end());
    EXPECT_NEAR(static_cast<double>(it->second), 64, 32);
    it = resultAmp.find("11");
    ASSERT_TRUE(it != resultAmp.end());
    EXPECT_NEAR(static_cast<double>(it->second), 64, 32);
}

using namespace qc::literals;
TEST(EntanglementSplitSimulatorDDSIM, grover5) {
    auto qc = std::make_unique<qc::QuantumComputation>("circuits/grover_5.qasm");

    EntanglementSplitSimulator ddsim(std::move(qc));

    std::vector<int64_t> times;
    std::map<std::string, std::size_t> resultAmp;
    for (size_t i = 0; i < NumberOfTimeMeasurements; ++i) {
        qc = std::make_unique<qc::QuantumComputation>("circuits/grover_5.qasm");

        EntanglementSplitSimulator ddsimNew(std::move(qc));
        auto start     = std::chrono::steady_clock::now();
        resultAmp = ddsimNew.simulate(512);
        times.push_back(since(start).count());
    }
    int64_t sum = 0;
    for (int64_t time : times) {
        sum += time;
    }
    std::cout << "Elapsed(ms)=" << sum / NumberOfTimeMeasurements << "\n";
    for (const auto& entry: resultAmp) {
        std::cout << "resultAmp[" << entry.first << "] = " << entry.second << "\n";
    }

    EXPECT_NEAR(resultAmp.size(), 16, 2);
    auto it = resultAmp.find("1000000");
    EXPECT_NEAR(static_cast<double>(it->second), 16, 40);
    it = resultAmp.find("1000001");
    EXPECT_NEAR(static_cast<double>(it->second), 8, 40);
    it = resultAmp.find("1000100");
    EXPECT_NEAR(static_cast<double>(it->second), 24, 40);
    it = resultAmp.find("1000101");
    EXPECT_NEAR(static_cast<double>(it->second), 8, 40);
    it = resultAmp.find("1000100");
    EXPECT_NEAR(static_cast<double>(it->second), 16, 40);
    it = resultAmp.find("1001101");
    EXPECT_NEAR(static_cast<double>(it->second), 8, 40);
    it = resultAmp.find("1001100");
    EXPECT_NEAR(static_cast<double>(it->second), 6, 40);
    it = resultAmp.find("1001101");
    EXPECT_NEAR(static_cast<double>(it->second), 12, 40);
    it = resultAmp.find("1100000");
    EXPECT_NEAR(static_cast<double>(it->second), 280, 100);
    it = resultAmp.find("1100001");
    EXPECT_NEAR(static_cast<double>(it->second), 12, 40);
    it = resultAmp.find("1100100");
    EXPECT_NEAR(static_cast<double>(it->second), 4, 40);
    it = resultAmp.find("1100101");
    EXPECT_NEAR(static_cast<double>(it->second), 20, 40);
    it = resultAmp.find("1101000");
    EXPECT_NEAR(static_cast<double>(it->second), 32, 40);
    it = resultAmp.find("1101001");
    EXPECT_NEAR(static_cast<double>(it->second), 24, 40);
    it = resultAmp.find("1101100");
    EXPECT_NEAR(static_cast<double>(it->second), 8, 40);
    it = resultAmp.find("1101101");
    EXPECT_NEAR(static_cast<double>(it->second), 6, 40);
}

using namespace qc::literals;
TEST(EntanglementSplitSimulatorDDSIM, grover5_simple) {
    auto qc = std::make_unique<qc::QuantumComputation>("circuits/grover_5.qasm");

    CircuitSimulator ddsim(std::move(qc));

    std::vector<int64_t> times;
    std::map<std::string, std::size_t> resultAmp;
    for (size_t i = 0; i < NumberOfTimeMeasurements; ++i) {
        qc = std::make_unique<qc::QuantumComputation>("circuits/grover_5.qasm");

        CircuitSimulator ddsimNew(std::move(qc));
        auto start     = std::chrono::steady_clock::now();
        resultAmp = ddsimNew.simulate(512);
        times.push_back(since(start).count());
    }
    int64_t sum = 0;
    for (int64_t time : times) {
        sum += time;
    }
    std::cout << "Elapsed(ms)=" << sum / NumberOfTimeMeasurements << "\n";
    for (const auto& entry: resultAmp) {
        std::cout << "resultAmp[" << entry.first << "] = " << entry.second << "\n";
    }

    EXPECT_NEAR(resultAmp.size(), 16, 2);
    auto it = resultAmp.find("1000000");
    EXPECT_NEAR(static_cast<double>(it->second), 16, 40);
    it = resultAmp.find("1000001");
    EXPECT_NEAR(static_cast<double>(it->second), 8, 40);
    it = resultAmp.find("1000100");
    EXPECT_NEAR(static_cast<double>(it->second), 24, 40);
    it = resultAmp.find("1000101");
    EXPECT_NEAR(static_cast<double>(it->second), 8, 40);
    it = resultAmp.find("1000100");
    EXPECT_NEAR(static_cast<double>(it->second), 16, 40);
    it = resultAmp.find("1001101");
    EXPECT_NEAR(static_cast<double>(it->second), 8, 40);
    it = resultAmp.find("1001100");
    EXPECT_NEAR(static_cast<double>(it->second), 6, 40);
    it = resultAmp.find("1001101");
    EXPECT_NEAR(static_cast<double>(it->second), 12, 40);
    it = resultAmp.find("1100000");
    EXPECT_NEAR(static_cast<double>(it->second), 280, 100);
    it = resultAmp.find("1100001");
    EXPECT_NEAR(static_cast<double>(it->second), 12, 40);
    it = resultAmp.find("1100100");
    EXPECT_NEAR(static_cast<double>(it->second), 4, 40);
    it = resultAmp.find("1100101");
    EXPECT_NEAR(static_cast<double>(it->second), 20, 40);
    it = resultAmp.find("1101000");
    EXPECT_NEAR(static_cast<double>(it->second), 32, 40);
    it = resultAmp.find("1101001");
    EXPECT_NEAR(static_cast<double>(it->second), 24, 40);
    it = resultAmp.find("1101100");
    EXPECT_NEAR(static_cast<double>(it->second), 8, 40);
    it = resultAmp.find("1101101");
    EXPECT_NEAR(static_cast<double>(it->second), 6, 40);
}
*/
using namespace qc::literals;
TEST(EntanglementSplitSimulatorDDSIM, grover15) {
    auto qc = std::make_unique<qc::QuantumComputation>("circuits/grover_15.qasm");

    EntanglementSplitSimulator ddsim(std::move(qc));

    std::vector<int64_t> times;
    std::map<std::string, std::size_t> resultAmp;
    for (size_t i = 0; i < NumberOfTimeMeasurements; ++i) {
        qc = std::make_unique<qc::QuantumComputation>("circuits/grover_15.qasm");

        EntanglementSplitSimulator ddsimNew(std::move(qc));
        auto start     = std::chrono::steady_clock::now();
        resultAmp = ddsimNew.simulate(512);
        times.push_back(since(start).count());
    }
    int64_t sum = 0;
    for (int64_t time : times) {
        sum += time;
    }
    std::cout << "Elapsed(ms)=" << sum / NumberOfTimeMeasurements << "\n";
}

TEST(EntanglementSplitSimulatorDDSIM, grover15_simple) {
    auto qc = std::make_unique<qc::QuantumComputation>("circuits/grover_15.qasm");

    CircuitSimulator ddsim(std::move(qc));

    std::vector<int64_t> times;
    std::map<std::string, std::size_t> resultAmp;
    for (size_t i = 0; i < NumberOfTimeMeasurements; ++i) {
        qc = std::make_unique<qc::QuantumComputation>("circuits/grover_15.qasm");
        CircuitSimulator ddsimNew(std::move(qc));
        auto start     = std::chrono::steady_clock::now();
        resultAmp = ddsimNew.simulate(512);
        times.push_back(since(start).count());
    }
    int64_t sum = 0;
    for (int64_t time : times) {
        sum += time;
    }
    std::cout << "Elapsed(ms)=" << sum / NumberOfTimeMeasurements << "\n";
}


/*using namespace qc::literals;
TEST(EntanglementSplitSimulatorDDSIM, grover7) {
    auto qc = std::make_unique<qc::QuantumComputation>("circuits/grover_7.qasm");

    EntanglementSplitSimulator ddsim(std::move(qc));
    auto start     = std::chrono::steady_clock::now();
    auto resultAmp = ddsim.simulate(512);
    std::cout << "Elapsed(ms)=" << since(start).count() << std::endl;
    for (const auto& entry: resultAmp) {
        std::cout << "resultAmp[" << entry.first << "] = " << entry.second << "\n";
    }
} */

/*using namespace qc::literals;
TEST(EntanglementSplitSimulatorDDSIM, grover7_simple) {
    auto qc = std::make_unique<qc::QuantumComputation>("circuits/grover_7.qasm");

    CircuitSimulator ddsim(std::move(qc));
    auto start     = std::chrono::steady_clock::now();
    auto resultAmp = ddsim.simulate(512);
    std::cout << "Elapsed(ms)=" << since(start).count() << std::endl;
    for (const auto& entry: resultAmp) {
        std::cout << "resultAmp[" << entry.first << "] = " << entry.second << "\n";
    }
} */
