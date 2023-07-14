#include "EntanglementSplitSimulator.hpp"

#include <gtest/gtest.h>
#include <memory>

using namespace qc::literals;
TEST(EntanglementSplitSimulatorDDSIM, TrivialTest) {
    auto quantumComputation = [] {
        auto qc = std::make_unique<qc::QuantumComputation>(2);
        qc->x(0);         // 1| h q[0];
        qc->h(1, {0_pc}); // 2| cx q[0], q[1];
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
