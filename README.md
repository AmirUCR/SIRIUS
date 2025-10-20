<img width="175" alt="SIRIUS Logo" src="https://github.com/user-attachments/assets/c9e8c503-9cdb-41fe-b060-0f5e1aa78760">

# Introduction
[Work in progress – Major release update coming late October] 

SIRIUS (_<ins>S</ins>ystematische <ins>I</ins>dentifikation <ins>R</ins>edundanter, <ins>I</ins>dentisch <ins>U</ins>ebersetzter <ins>S</ins>equenzen_) is a synthetic biology tool leveraging Google OR-Tools mixed-integer programming to design genetic sequences with the shortest and fewest possible homologous fragments between pairs within minutes.

- Design _n_ gene sequences all translating to a given protein _P_
- Synthesize sequences with maximal, optimal divergence
- Written in pure C++
  
<img width="1310" height="362" alt="image" src="https://github.com/user-attachments/assets/f06cae9f-d7f9-44a0-85d5-5097706e0590" />

__Overview of the SIRIUS workflow.__ __Step 1__: The input to SIRIUS is a protein sequence of interest (_P_) and
the desired number (_n_) of synonymous DNA sequences to be designed. __Step 2__: SIRIUS formulates and
solves an integer linear program (ILP) containing millions of variables and constraints. Illustration shows
the objective function described in the Methods section, together with the codon choices for each amino
acid in the example peptide _P_ from Step 1. __Step 3__: The resulting feasible or optimal solution consists
of n output synonymous DNA sequences with the fewest and shortest homologous fragments between all
sequence pairs. Highlights indicate homologous fragments between any pair.

# Documentation
You may find the documentation for SIRIUS at its [GitHub Wiki](https://github.com/AmirUCR/SIRIUS/wiki).

# Support
If you run into any issues or have suggestions for SIRIUS, please report them on our GitHub Issues tracker. It's the fastest way to get support and helps us improve SIRIUS for everyone. You may also email the authors at their provided e-mail addresses on the publication directly.

# About
SIRIUS has been developed and is maintained by <ins>Amir</ins>sadra Mohseni, and Stefano Lonardi at the University of California, Riverside.
