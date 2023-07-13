#include "EntanglementSplitSimulator.hpp"

#include <gtest/gtest.h>
#include <memory>

using namespace qc::literals;
TEST(EntanglementSplitSimulatorDDSIM, TrivialTest) {
    auto quantumComputation = [] {
        auto qc = std::make_unique<qc::QuantumComputation>(3);
        qc->h(0);
        qc->x(1, {0});
        //qc->i(1); // some dummy operations
       // qc->i(1);
        return qc;
    };

    EntanglementSplitSimulator ddsim(quantumComputation());

    auto resultAmp = ddsim.simulate(8192);
    for (const auto& entry: resultAmp) {
        std::cout << "resultAmp[" << entry.first << "] = " << entry.second << "\n";
    }

    ASSERT_EQ(resultAmp.size(), 4);
    auto it = resultAmp.find("0000");
    ASSERT_TRUE(it != resultAmp.end());
    EXPECT_NEAR(static_cast<double>(it->second), 2048, 128);
    it = resultAmp.find("0010");
    ASSERT_TRUE(it != resultAmp.end());
    EXPECT_NEAR(static_cast<double>(it->second), 2048, 128);
    it = resultAmp.find("0100");
    ASSERT_TRUE(it != resultAmp.end());
    EXPECT_NEAR(static_cast<double>(it->second), 2048, 128);
    it = resultAmp.find("1110");
    ASSERT_TRUE(it != resultAmp.end());
    EXPECT_NEAR(static_cast<double>(it->second), 2048, 128);
}
