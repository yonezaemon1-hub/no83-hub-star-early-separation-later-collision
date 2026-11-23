# Research No. 83 — public preprint v1.0

**Early-Time Edge Identifiability and a Later Observation Collision in Hub-Star Kuramoto Networks**

Author: Ryutaro Yonezu. Version: v1.0 (24 November 2026). All rights reserved.

The manuscript establishes graph-uniform early-time separation and a computer-assisted five-vertex collision at a later observation time. See the PDF for exact assumptions, proof, bounds, and limitations. Independent AI mathematical audit: PASS for the stated results and code inspection; not human peer reviewed. Novelty: INCONCLUSIVE. References: PARTIAL.

The author reported a successful Windows compile and run; its output was transcribed in the private working record. The original native Windows log and its original-file hash were not provided, so original-file provenance remains INCONCLUSIVE. The independent AI auditor inspected the MPFR source and used separate calculations; it did not execute the author's MPFR binary.

## Reproduce the certificate

Requires a C++ compiler, MPFR and GMP development libraries, and Python 3 (for the optional reduction check). From a directory containing the two source files:

```sh
g++ -O2 -std=c++17 no83_stage29_cert.cpp -lmpfr -lgmp -o no83_cert
./no83_cert
python3 no83_reduced_check.py
```

The C++ code computes interval enclosures for discrete midpoint iterates. The manuscript supplies the separate analytic bound needed to infer signs for the exact ODE solutions. See `SHA256SUMS.txt` for public file hashes.
