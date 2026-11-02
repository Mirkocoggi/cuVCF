#pragma once
#include <filesystem>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>

namespace test {
struct Assertion : std::runtime_error {
    std::string key;
    Assertion(std::string k, const std::string& detail)
        : std::runtime_error(k + ": " + detail), key(std::move(k)) {}
};
inline void require(bool ok, const std::string& key, const std::string& detail) {
    if (!ok) throw Assertion(key, detail);
}
template<class A, class B> void expect_equal(const A& actual, const B& expected, const std::string& key) {
    if (!(actual == expected)) {
        std::ostringstream out; out << "expected " << expected << ", got " << actual;
        throw Assertion(key, out.str());
    }
}
class TempDir {
    std::filesystem::path path_;
public:
    TempDir() {
        auto pattern = (std::filesystem::temp_directory_path() / "cuvcf-tests-XXXXXX").string();
        if (!mkdtemp(pattern.data())) throw std::runtime_error("cannot create test directory");
        path_ = pattern;
    }
    ~TempDir() { std::error_code ec; std::filesystem::remove_all(path_, ec); }
    std::filesystem::path file(const std::string& name) const { return path_ / name; }
    std::string write(const std::string& name, const std::string& data) const {
        auto path = file(name); std::ofstream out(path, std::ios::binary);
        out.exceptions(std::ios::badbit | std::ios::failbit); out << data; out.close();
        return path.string();
    }
    TempDir(const TempDir&) = delete;
    TempDir& operator=(const TempDir&) = delete;
};
struct Case {
    std::string id;
    std::function<void()> run;
    bool extended = false;
    // Only this exact assertion is eligible; exceptions/setup errors never count as xfail.
    std::string expected_failure;
    Case(std::string name, std::function<void()> fn, bool extra=false, std::string failure={})
        : id(std::move(name)), run(std::move(fn)), extended(extra), expected_failure(std::move(failure)) {}
};
enum class Result { pass, fail, xfail, xpass };
inline Result execute(const Case& c, std::ostream& log) {
    try {
        c.run();
        if (!c.expected_failure.empty()) { log << "XPASS " << c.id << '\n'; return Result::xpass; }
        log << "PASS " << c.id << '\n'; return Result::pass;
    } catch (const Assertion& e) {
        if (!c.expected_failure.empty() && e.key == c.expected_failure) {
            log << "XFAIL " << c.id << " " << e.what() << '\n'; return Result::xfail;
        }
        log << "FAIL " << c.id << " " << e.what() << '\n'; return Result::fail;
    } catch (const std::exception& e) {
        log << "ERROR " << c.id << " " << e.what() << '\n'; return Result::fail;
    }
}
inline int run(const std::vector<Case>& cases, int argc, char** argv, std::ostream& output=std::cout) {
    const auto started=std::chrono::steady_clock::now();
    bool extended = false, list = false; std::string filter;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--extended") extended = true;
        else if (arg == "--list") list = true;
        else if (arg == "--filter" && i + 1 < argc) filter = argv[++i];
        else { std::cerr << "Usage: " << argv[0] << " [--list] [--extended] [--filter substring]\n"; return 2; }
    }
    int passed = 0, failed = 0, known = 0, selected = 0;
    for (const auto& c : cases) {
        if ((c.extended && !extended) || c.id.find(filter) == std::string::npos) continue;
        ++selected;
        if (list) { output << c.id << '\n'; continue; }
        switch (execute(c, output)) {
        case Result::pass: ++passed; break;
        case Result::xfail: ++known; break;
        default: ++failed; break;
        }
    }
    if (!selected) { std::cerr << "No cases selected\n"; return 2; }
    if (!list) output << passed << " passed, " << known << " expected failures, " << failed << " failed\n";
    if (!list) output << "Elapsed: " << std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count() << " s\n";
    return failed ? 1 : 0;
}
}
