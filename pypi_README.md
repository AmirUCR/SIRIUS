<img width="175" alt="SIRIUS Logo" src="https://github.com/user-attachments/assets/c9e8c503-9cdb-41fe-b060-0f5e1aa78760">

# Introduction
SIRIUS (_<ins>S</ins>ystematische <ins>I</ins>dentifikation <ins>R</ins>edundanter, <ins>I</ins>dentisch <ins>U</ins>ebersetzter <ins>S</ins>equenzen_) is a synthetic biology tool leveraging Google OR-Tools integer programming to design genetic sequences with the shortest and fewest possible homologous fragments between pairs within minutes.

- Design _n_ gene sequences all translating to a given protein _P_
- Effectively synthesize sequences with maximal, optimal divergence
- Written in pure C++

**Overview of the SIRIUS workflow.** **Step (1)** The input to SIRIUS is a protein sequence of interest _P_ and the desired number _n_ of synonymous DNA sequences to be designed; **Step (2)** SIRIUS solves an integer linear program with the objective function detailed in its publication, and millions of variables and constraints; **Step (3)** SIRIUS produces of _n_ synonymous DNA sequences that encode _P_ with the fewest and shortest common subsequences.

# Documentation
You may find the documentation for SIRIUS at its [GitHub Wiki](https://github.com/ucrbioinfo/SIRIUS/wiki).

# Support
If you run into any issues or have suggestions for SIRIUS, please report them on our GitHub Issues tracker. It's the fastest way to get support and helps us improve SIRIUS for everyone. You may also email the authors at their provided e-mail addresses on the publication directly.

# About
SIRIUS has been developed and is maintained by <ins>Amir</ins>sadra Mohseni, and Stefano Lonardi at the University of California, Riverside.
