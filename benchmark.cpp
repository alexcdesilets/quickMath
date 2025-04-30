#include "quickMath.h"

using namespace quickMath;

int main() {
    _initQuickMath();

    const int N = 1000000;
    std::vector<double> angles(N);
    std::vector<double> log_inputs(N);
    std::vector<double> unlog_inputs(N);
    std::vector<double> pow_bases(N);
    std::vector<double> pow_exponents(N);
    std::vector<double> tan_safe_angles(N);
    std::vector<double> unit_circle_x(N), unit_circle_y(N); // for atan2
    std::vector<double> acos_inputs(N), asin_inputs(N);

    std::cout << "Preparing " << N << " test inputs..." << std::endl;
    for (int i = 0; i < N; ++i) {
        double t = static_cast<double>(i) / static_cast<double>(N - 1);
        log_inputs[i] = 0.01 + t * (1000.0 - 0.01);
        unlog_inputs[i] = -3.0 + t * 7.0;
        pow_bases[i] = 0.1 + t * 999.9;
        pow_exponents[i] = -2.0 + 4.0 * t;
        angles[i] = -4.0 * quickMath::_internal::_PI + t * (8.0 * quickMath::_internal::_PI);
        tan_safe_angles[i] = -1.4 * quickMath::_internal::_PI + t * (2.8 * quickMath::_internal::_PI);
        unit_circle_x[i] = std::cos(angles[i]);
        unit_circle_y[i] = std::sin(angles[i]);
        asin_inputs[i] = std::sin(t * quickMath::_internal::_PI - quickMath::_internal::_PI / 2.0);  // in [-1, 1]
        acos_inputs[i] = std::cos(t * quickMath::_internal::_PI);               // in [-1, 1]
    }
    std::cout << "Inputs prepared." << std::endl;

    std::vector<double> qlog_results(N), stdlog_results(N);
    std::vector<double> qunlog_results(N), stdpow10_results(N);
    std::vector<double> quickpow_results(N), stdpow_results(N);
    std::vector<double> qsin_results(N), stdsin_results(N);
    std::vector<double> qcos_results(N), stdcos_results(N);
    std::vector<double> qtan_results(N), stdtan_results(N);
    std::vector<double> qasin_results(N), stdasin_results(N);
    std::vector<double> qacos_results(N), stdacos_results(N);
    std::vector<double> qatan2_results(N), stdatan2_results(N);

    std::cout << std::fixed << std::setprecision(16);
    std::cout << "\n--- Benchmarking ---" << std::endl;

    auto time = [](auto&& func) {
        auto start = std::chrono::high_resolution_clock::now();
        func();
        return std::chrono::duration<double, std::milli>(std::chrono::high_resolution_clock::now() - start).count();
    };

    std::cout << "quickLog time:         " << time([&]() {
        for (int i = 0; i < N; ++i) qlog_results[i] = quickLog(log_inputs[i]);
    }) << " ms\n";

    std::cout << "std::log10 time:       " << time([&]() {
        for (int i = 0; i < N; ++i) stdlog_results[i] = std::log10(log_inputs[i]);
    }) << " ms\n";

    std::cout << "quickUnlog time:       " << time([&]() {
        for (int i = 0; i < N; ++i) qunlog_results[i] = quickUnlog(unlog_inputs[i]);
    }) << " ms\n";

    std::cout << "std::pow(10,x) time:   " << time([&]() {
        for (int i = 0; i < N; ++i) stdpow10_results[i] = std::pow(10.0, unlog_inputs[i]);
    }) << " ms\n";

    std::cout << "quickPow time:         " << time([&]() {
        for (int i = 0; i < N; ++i) quickpow_results[i] = quickPow(pow_bases[i], pow_exponents[i]);
    }) << " ms\n";

    std::cout << "std::pow time:         " << time([&]() {
        for (int i = 0; i < N; ++i) stdpow_results[i] = std::pow(pow_bases[i], pow_exponents[i]);
    }) << " ms\n";

    std::cout << "quickSin time:         " << time([&]() {
        for (int i = 0; i < N; ++i) qsin_results[i] = quickSin(angles[i]);
    }) << " ms\n";

    std::cout << "std::sin time:         " << time([&]() {
        for (int i = 0; i < N; ++i) stdsin_results[i] = std::sin(angles[i]);
    }) << " ms\n";

    std::cout << "quickCos time:         " << time([&]() {
        for (int i = 0; i < N; ++i) qcos_results[i] = quickCos(angles[i]);
    }) << " ms\n";

    std::cout << "std::cos time:         " << time([&]() {
        for (int i = 0; i < N; ++i) stdcos_results[i] = std::cos(angles[i]);
    }) << " ms\n";

    std::cout << "quickTan time:         " << time([&]() {
        for (int i = 0; i < N; ++i) qtan_results[i] = quickTan(tan_safe_angles[i]);
    }) << " ms\n";

    std::cout << "std::tan time:         " << time([&]() {
        for (int i = 0; i < N; ++i) stdtan_results[i] = std::tan(tan_safe_angles[i]);
    }) << " ms\n";

    std::cout << "quickAsin time:        " << time([&]() {
        for (int i = 0; i < N; ++i) qasin_results[i] = quickAsin(asin_inputs[i]);
    }) << " ms\n";

    std::cout << "std::asin time:        " << time([&]() {
        for (int i = 0; i < N; ++i) stdasin_results[i] = std::asin(asin_inputs[i]);
    }) << " ms\n";

    std::cout << "quickAcos time:        " << time([&]() {
        for (int i = 0; i < N; ++i) qacos_results[i] = quickAcos(acos_inputs[i]);
    }) << " ms\n";

    std::cout << "std::acos time:        " << time([&]() {
        for (int i = 0; i < N; ++i) stdacos_results[i] = std::acos(acos_inputs[i]);
    }) << " ms\n";

    std::cout << "quickAtan2 time:       " << time([&]() {
        for (int i = 0; i < N; ++i) qatan2_results[i] = quickAtan2(unit_circle_y[i], unit_circle_x[i]);
    }) << " ms\n";

    std::cout << "std::atan2 time:       " << time([&]() {
        for (int i = 0; i < N; ++i) stdatan2_results[i] = std::atan2(unit_circle_y[i], unit_circle_x[i]);
    }) << " ms\n";

    // --- Error Analysis ---
    double log_error_sum_abs = 0.0;
    double unlog_error_sum_rel = 0.0, pow_error_sum_rel = 0.0;
    double sin_error_sum_rel = 0.0, cos_error_sum_rel = 0.0, tan_error_sum_rel = 0.0;
    double asin_error_sum_rel = 0.0, acos_error_sum_rel = 0.0, atan2_error_sum_rel = 0.0;
    int unlog_error_count = 0, pow_error_count = 0;
    int sin_error_count = 0, cos_error_count = 0, tan_error_count = 0;
    int asin_error_count = 0, acos_error_count = 0, atan2_error_count = 0;

    for (int i = 0; i < N; ++i) {
        if (std::isfinite(qlog_results[i]) && std::isfinite(stdlog_results[i]))
            log_error_sum_abs += std::abs(qlog_results[i] - stdlog_results[i]);

        if (std::isfinite(qunlog_results[i]) && std::isfinite(stdpow10_results[i])) {
            double s = std::abs(stdpow10_results[i]);
            if (s > 1e-15) {
                unlog_error_sum_rel += std::abs(qunlog_results[i] - stdpow10_results[i]) / s;
                ++unlog_error_count;
            }
        }

        if (std::isfinite(quickpow_results[i]) && std::isfinite(stdpow_results[i])) {
            double s = std::abs(stdpow_results[i]);
            if (s > 1e-15) {
                pow_error_sum_rel += std::abs(quickpow_results[i] - stdpow_results[i]) / s;
                ++pow_error_count;
            }
        }

        if (std::isfinite(qsin_results[i]) && std::isfinite(stdsin_results[i])) {
            double s = std::abs(stdsin_results[i]);
            if (s > 1e-15) {
                sin_error_sum_rel += std::abs(qsin_results[i] - stdsin_results[i]) / s;
                ++sin_error_count;
            }
        }

        if (std::isfinite(qcos_results[i]) && std::isfinite(stdcos_results[i])) {
            double s = std::abs(stdcos_results[i]);
            if (s > 1e-15) {
                cos_error_sum_rel += std::abs(qcos_results[i] - stdcos_results[i]) / s;
                ++cos_error_count;
            }
        }

        if (std::isfinite(qtan_results[i]) && std::isfinite(stdtan_results[i])) {
            double s = std::abs(stdtan_results[i]);
            if (s > 1e-6) {
                tan_error_sum_rel += std::abs(qtan_results[i] - stdtan_results[i]) / s;
                ++tan_error_count;
            }
        }

        if (std::isfinite(qasin_results[i]) && std::isfinite(stdasin_results[i])) {
            double s = std::abs(stdasin_results[i]);
            if (s > 1e-15) {
                asin_error_sum_rel += std::abs(qasin_results[i] - stdasin_results[i]) / s;
                ++asin_error_count;
            }
        }

        if (std::isfinite(qacos_results[i]) && std::isfinite(stdacos_results[i])) {
            double s = std::abs(stdacos_results[i]);
            if (s > 1e-15) {
                acos_error_sum_rel += std::abs(qacos_results[i] - stdacos_results[i]) / s;
                ++acos_error_count;
            }
        }

        if (std::isfinite(qatan2_results[i]) && std::isfinite(stdatan2_results[i])) {
            double s = std::abs(stdatan2_results[i]);
            if (s > 1e-15) {
                atan2_error_sum_rel += std::abs(qatan2_results[i] - stdatan2_results[i]) / s;
                ++atan2_error_count;
            }
        }
    }

    std::cout << "\n--- Average Errors ---" << std::endl;
    std::cout << "Avg Abs Error (quickLog vs std::log10): " << (log_error_sum_abs / N) << std::endl;
    if (unlog_error_count)
        std::cout << "Avg Rel Error (quickUnlog vs std::pow10): " << (unlog_error_sum_rel / unlog_error_count) << std::endl;
    if (pow_error_count)
        std::cout << "Avg Rel Error (quickPow vs std::pow): " << (pow_error_sum_rel / pow_error_count) << std::endl;
    if (sin_error_count)
        std::cout << "Avg Rel Error (quickSin vs std::sin): " << (sin_error_sum_rel / sin_error_count) << std::endl;
    if (cos_error_count)
        std::cout << "Avg Rel Error (quickCos vs std::cos): " << (cos_error_sum_rel / cos_error_count) << std::endl;
    if (tan_error_count)
        std::cout << "Avg Rel Error (quickTan vs std::tan): " << (tan_error_sum_rel / tan_error_count) << std::endl;
    if (asin_error_count)
        std::cout << "Avg Rel Error (quickAsin vs std::asin): " << (asin_error_sum_rel / asin_error_count) << std::endl;
    if (acos_error_count)
        std::cout << "Avg Rel Error (quickAcos vs std::acos): " << (acos_error_sum_rel / acos_error_count) << std::endl;
    if (atan2_error_count)
        std::cout << "Avg Rel Error (quickAtan2 vs std::atan2): " << (atan2_error_sum_rel / atan2_error_count) << std::endl;

    return 0;
}
