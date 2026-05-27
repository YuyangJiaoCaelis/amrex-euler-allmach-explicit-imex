# Research Baseline - 2026-05-27

This repository is the clean AMReX Euler research baseline for work after Report 1.

It contains the current structured implementation of:

- Explicit O2 HLLC
- Explicit O2 Low-Mach Corrected HLLC-P
- IMEX T1/S2 BDLTV20
- Report 1 reproduction scripts for Riemann, Gresho vortex, periodic advection blob, and same-gamma shock-density-bubble tests

Use this repository as the starting point for future research tasks, including Report 2 development, hardware timing, MPI/CUDA preparation, solver refinement, and additional validation.

The historical parent folder contains archived design notes, review logs, old prototypes, and report-writing material. Those files are useful for reference, but new code development should start from this baseline.
