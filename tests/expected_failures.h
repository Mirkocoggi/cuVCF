#pragma once
#include "harness.h"
#include <set>
#include <algorithm>
namespace test {
struct KnownFailure {
    std::string backend, case_id, assertion, reason;
};
inline void apply_expected_failures(std::vector<Case>& cases) {
#ifdef CUVCF_TEST_GPU
    const std::string backend="gpu";
#else
    const std::string backend="cpu";
#endif
    // Populate only after a native reproduction, with a matching README ledger entry.
    // {"cpu", "write/example", "exact/assertion/key", "issue reference and cause"},
    const std::vector<KnownFailure> known{};
    std::set<std::string> seen;
    for(const auto& failure:known) {
        if(failure.backend!=backend) continue;
        if(failure.assertion.empty() || failure.reason.empty() || !seen.insert(failure.case_id).second)
            throw std::runtime_error("invalid expected-failure entry: "+failure.case_id);
        auto found=std::find_if(cases.begin(),cases.end(),[&](const Case& c){return c.id==failure.case_id;});
        if(found==cases.end()) throw std::runtime_error("stale expected-failure entry: "+failure.case_id);
        found->expected_failure=failure.assertion;
    }
}
}
