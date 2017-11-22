# MusicProfileHMMs
Python routines for multiple sequence alignment, generating profile HMMs, focused on music sequences.

These are my standard set of routines for data-driven analysis of unidimensional music sequences. They contain functions for pairwise alignment (based on BioPython), multiple sequence alignment (using MAFFT or progressive alignment with T-COFFEE), substitution matrix creation from alignments (based on BioPython), consensus generation from an MSA (majority voting, or data fusion), profile HMM modelling of an MSA using Krogh's original model (1994) and HMMER.

