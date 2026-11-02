# Verification scripts

Independent checks backing the findings in [CODE_REVIEW.md](../../CODE_REVIEW.md). The Python scripts port SCESS algorithms to NumPy and compare them against PySCF. The C++ files are minimal programs that compile against the real SCESS and Dahlia headers.

Requirements: Python 3 with `numpy` and `pyscf` (tested with PySCF 2.14), and a C++20 compiler.

Run everything from this directory, because `mp2f12_port.py` and `gamma_scan.py` import `f12_ints.py` from it.

| Script | Checks | Review section |
|---|---|---|
| `basis_check.py` | Effect of the basis-loader contraction merge on H₂O/STO-3G RHF | §3.1 |
| `cc_check.py` | CCSD (`Wmbej` factor) and CCSD(T) (permutations, E[4] bra) vs PySCF | §5 |
| `ci_check.py` | FCI Slater–Condon rules and phase algorithm vs PySCF (LiH) | §6 |
| `notation_check.py` | Chemists'/physicists' mismatch in CCSD-QMC[MP2-F12] (H₂) | §7.3 |
| `cabs_check2.py` | CABS canonicalisation and CABS singles (H₂O; edit `geom` for H₂) | §8.1 |
| `f12_ints.py` | Gaussian-geminal integrals from PySCF range-separated ERIs; run it directly for self-tests | §8 |
| `mp2f12_port.py` | MP2-F12 port, Fock trace bug, independent B, optimised energies | §8 |
| `gamma_scan.py` | Geminal-exponent scan against the MRCC value in `Tests.cpp` | §8.2 |
| `einsum_trace.cpp` | `Einsum` with a repeated index sums over both positions instead of taking the trace | §8.2 |
| `selfjoin.cpp` | `DH::Thread::MarkAsFinished()` self-join aborts the process | §2.2 |

## Usage

```sh
python3 basis_check.py
python3 cc_check.py
python3 ci_check.py
python3 notation_check.py
python3 cabs_check2.py
python3 f12_ints.py                 # self-tests, ~1 s
python3 mp2f12_port.py H2           # ~5 s
python3 mp2f12_port.py H2O          # ~3.5 min
python3 gamma_scan.py               # ~30 s

c++ -std=c++20 -I../../src -I../../Libraries/Eigen einsum_trace.cpp -o einsum_trace && ./einsum_trace
c++ -std=c++20 -I../../../Dahlia/src -DDH_API= selfjoin.cpp -o selfjoin && ./selfjoin
```

Both bugs are now fixed. `Einsum` rejects repeated indices within one operand (traces) at compile time, so `einsum_trace.cpp` is expected to fail to compile with the `static_assert` "Repeated indices (traces) are not supported by Einsum". `selfjoin` exits cleanly.

The deterministic results are also checked by the regression tests, `SCESS test` (run from `bin/<platform>/<config>`).
