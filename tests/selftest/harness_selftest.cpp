#include "support/harness.h"
int main(int argc, char** argv) {
    using namespace test;
    std::vector<Case> cases;
    cases.push_back({"harness/outcomes", [] {
        std::ostringstream sink;
        expect_equal(int(execute({"pass", [] {}}, sink)), int(Result::pass), "pass");
        auto assertion = [] { throw Assertion("value", "wrong value"); };
        expect_equal(int(execute({"fail", assertion}, sink)), int(Result::fail), "failure");
        expect_equal(int(execute({"known", assertion, false, "value"}, sink)), int(Result::xfail), "xfail");
        expect_equal(int(execute({"fixed", [] {}, false, "value"}, sink)), int(Result::xpass), "xpass");
        expect_equal(int(execute({"unrelated", assertion, false, "other"}, sink)), int(Result::fail), "unrelated");
        expect_equal(int(execute({"setup", [] { throw std::runtime_error("setup"); }, false, "value"}, sink)), int(Result::fail), "setup");
    }});
    cases.push_back({"harness/exit_status", [] {
        char name[]="nested-run"; char* args[]={name}; std::ostringstream sink;
        expect_equal(run({{"probe", [] { throw Assertion("probe", "intentional"); }}},1,args,sink),1,"failure_exit");
        expect_equal(run({{"probe", [] {},false,"probe"}},1,args,sink),1,"xpass_exit");
    }});
    cases.push_back({"harness/cleanup", [] {
        std::filesystem::path path;
        { TempDir tmp; path = tmp.file("input"); tmp.write("input", "fixture"); require(std::filesystem::exists(path), "created", "fixture missing"); }
        require(!std::filesystem::exists(path), "cleanup", "fixture left behind");
    }});
    return run(cases, argc, argv);
}
