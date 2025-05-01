#include "quickMath.h"
using namespace quickMath;

int main() {
    _initQuickMath();

    const int N = 1000000;
    std::vector<double> angles(N);
    std::vector<double> log_inputs(N);
    std::vector<double> unlog_inputs(N);
    std::vector<double> pow_bases(N), pow_exponents(N);
    std::vector<double> asin_inputs(N), acos_inputs(N);
    std::vector<double> sin_inputs(N), cos_inputs(N), tan_inputs(N);
    std::vector<double> x_unit(N), y_unit(N); // for atan2

    for (int i = 0; i < N; ++i) {
        double t = static_cast<double>(i) / static_cast<double>(N - 1);
        log_inputs[i] = 0.01 + t * (1000.0 - 0.01);
        unlog_inputs[i] = -3.0 + t * 7.0;
        pow_bases[i] = 0.1 + t * 999.9;
        pow_exponents[i] = -2.0 + 4.0 * t;
        angles[i] = -4.0 * _internal::_PI + t * 8.0 * _internal::_PI;
        sin_inputs[i] = cos_inputs[i] = tan_inputs[i] = angles[i];
        x_unit[i] = std::cos(angles[i]);
        y_unit[i] = std::sin(angles[i]);
        asin_inputs[i] = y_unit[i]; // in [-1,1]
        acos_inputs[i] = x_unit[i];
    }

    std::vector<double> qlog(N), slog(N);
    std::vector<double> qunlog(N), sunlog(N);
    std::vector<double> qpow(N), spow(N);
    std::vector<double> qsin(N), ssin(N);
    std::vector<double> qcos(N), scos(N);
    std::vector<double> qtan(N), stan(N);
    std::vector<double> qasin(N), sasin(N);
    std::vector<double> qacos(N), sacos(N);
    std::vector<double> qatan2(N), satan2(N);

    auto time = [](auto&& func) {
        auto start = std::chrono::high_resolution_clock::now();
        func();
        return std::chrono::duration<double, std::milli>(
            std::chrono::high_resolution_clock::now() - start).count();
    };

    std::cout << std::fixed << std::setprecision(16);
    std::cout << "\n--- Mixed Function Timing (Realistic Usage) ---\n";

    std::cout << "quickMath total time: " << time([&]() {
        for (int i = 0; i < N; ++i) {
            qlog[i]    = quickLog(log_inputs[i]);
            qunlog[i]  = quickUnlog(unlog_inputs[i]);
            qpow[i]    = quickPow(pow_bases[i], pow_exponents[i]);
            qsin[i]    = quickSin(sin_inputs[i]);
            qcos[i]    = quickCos(cos_inputs[i]);
            qtan[i]    = quickTan(tan_inputs[i]);
            qasin[i]   = quickAsin(asin_inputs[i]);
            qacos[i]   = quickAcos(acos_inputs[i]);
            qatan2[i]  = quickAtan2(y_unit[i], x_unit[i]);
        }
    }) << " ms\n";

    std::cout << "cmath total time: " << time([&]() {
        for (int i = 0; i < N; ++i) {
            slog[i]    = std::log10(log_inputs[i]);
            sunlog[i]  = std::pow(10.0, unlog_inputs[i]);
            spow[i]    = std::pow(pow_bases[i], pow_exponents[i]);
            ssin[i]    = std::sin(sin_inputs[i]);
            scos[i]    = std::cos(cos_inputs[i]);
            stan[i]    = std::tan(tan_inputs[i]);
            sasin[i]   = std::asin(asin_inputs[i]);
            sacos[i]   = std::acos(acos_inputs[i]);
            satan2[i]  = std::atan2(y_unit[i], x_unit[i]);
        }
    }) << " ms\n";

    // ---------------------- Error Analysis Per Function ---------------------
    auto rel_error = [](double a, double b) -> double {
        return (std::abs(b) > 1e-15) ? std::abs(a - b) / std::abs(b) : 0.0;
    };

    auto abs_error = [](double a, double b) -> double {
        return std::abs(a - b);
    };

    auto analyze = [&](const std::vector<double>& q, const std::vector<double>& s,
                       std::string name, bool useRel = true) {
        double sum = 0.0;
        int count = 0;
        for (int i = 0; i < N; ++i) {
            if (std::isfinite(q[i]) && std::isfinite(s[i])) {
                double err = useRel ? rel_error(q[i], s[i]) : abs_error(q[i], s[i]);
                sum += err;
                ++count;
            }
        }
        std::cout << (useRel ? "Avg Rel Error" : "Avg Abs Error")
                  << " (" << name << "): " << (sum / count) << "\n";
    };

    std::cout << "\n--- Mixed Function Errors (Per Function) ---\n";
    analyze(qlog,    slog,   "quickLog vs std::log10", false);
    analyze(qunlog,  sunlog, "quickUnlog vs std::pow10");
    analyze(qpow,    spow,   "quickPow vs std::pow");
    analyze(qsin,    ssin,   "quickSin vs std::sin");
    analyze(qcos,    scos,   "quickCos vs std::cos");
    analyze(qtan,    stan,   "quickTan vs std::tan");
    analyze(qasin,   sasin,  "quickAsin vs std::asin");
    analyze(qacos,   sacos,  "quickAcos vs std::acos");
    analyze(qatan2,  satan2, "quickAtan2 vs std::atan2");

    return 0;
}