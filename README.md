
# quickMath

**Author:** Alexander Desilets  
**Email:** AlexCDesilets@gmail.com  

## License

This project is licensed under the GNU General Public License v3.0 (GPLv3). You are free to use, modify, and redistribute the code under the terms of that license.

However, if you wish to use this library in a **proprietary or closed-source commercial project**, a **separate commercial license is required**.

### 🔒 Commercial Licensing

To use `quickMath` in a closed-source or commercial setting without releasing your own code under GPLv3, please contact:

**Alexander Desilets**  
📧 Email: AlexCDesilets@gmail.com

You’ll receive a commercial license that grants you the right to use `quickMath` without the obligations of GPLv3.

## Overview

`quickMath` is a C++ header-only library that provides **fast approximations** for common math functions using **lookup tables**. It is designed for **performance-critical** applications such as real-time graphics or simulations where **maximum speed** is required and **ultra-high precision** is not the priority.

The design is inspired by slide-rule thinking: **fast but less accurate** than native `std::` functions. Many functions achieve **~10⁻¹⁰ to 10⁻¹¹ relative error**.

---

## Key Features

- Fast approximations of common math functions: `log`, `pow`, `exp`, `sqrt`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan2`
- Uses **precomputed caches** and **linear/bilinear interpolation**
- Cache presets: `smallCache`, `defaultCache`, `bigCache` with tradeoffs in memory and accuracy
- Suitable for embedded or performance-bound systems

---

## Build Instruction

You must match "g++ main.cpp -o main.exe -std=c++23 -O3 -march=native -ffast-math" in your compiler of choice to receive maximum benefit

---

## Initialization

You **must call** `_initQuickMath()` **once** before using any of the `quick...` functions:

```cpp
quickMath::_initQuickMath(quickMath::initVals::defaultCache);
```

This builds internal caches. It's **not thread-safe**—initialize from the main thread.

---

## Accuracy Notes

| Function      | Typical Accuracy (default cache) | Notes |
|---------------|----------------------------------|-------|
| `quickLog`    | ~10⁻¹¹                            | Very good |
| `quickExp`    | ~10⁻¹¹                            | Good |
| `quickPow`    | ~10⁻¹⁰                            | Good |
| `quickSin`    | ~10⁻⁸                             | Good |
| `quickCos`    | ~10⁻¹⁰                            | Very good |
| `quickTan`    | Poor near asymptotes             | Use `quickSin / quickCos` |
| `quickAsin`   | ~10⁻⁶                             | Moderate |
| `quickAcos`   | ~10⁻³                             | Poor near edges |
| `quickAtan2`  | ~10⁻⁴                             | Interpolated from grid |

---

## Notes

- The cache size significantly affects accuracy and memory usage.
- Caches are cleared and rebuilt if `_initQuickMath()` is called again with a different size.
- Use `#define QUICKMATH_VERBOSE` before including the header to see cache generation logs.

---

## Warnings

- `quickTan`, `quickAcos`, and `quickAtan2` are **less accurate**, especially near function discontinuities.
- All approximations assume **finite input**. Invalid input (NaN/Inf) returns NaN.
- The library is **not thread-safe** during initialization.

