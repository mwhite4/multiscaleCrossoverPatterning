# multiscaleCrossoverPatterning

Custom code (Matlab functions/scripts) associated with publication:

White, M. A., Weiner, B., Chu, L., Lim, G., and Kleckner, N. (2024).  Crossover Interference Mediates Multiscale Patterning Along Meiotic Chromosomes. bioRxiv. https://doi.org/10.1101/2024.01.28.577645.

August 2025:
The above manuscript has been accepted, in principle, for publication in Nature Communications

White, M. A., Weiner, B., Lim, G., Prentiss, M. Chu, L., and Kleckner, N. (2025).  Crossover Interference Mediates Multiscale Patterning Along Meiotic Chromosomes. Nature Communications, in press.


Abstract

Meiotic crossover interference is a one-dimensional spatial patterning process that produces evenly-spaced crossovers.  Quantitative analysis of diagnostic molecules along budding yeast chromosomes reveals that this process sets up two interdigitated patterns, of shorter and longer periodicity, by "two-tiered" patterning.  Both tiers comprise clustered assemblies of three types of molecules ("triads") representing the three major components of meiotic chromosomes (crossover recombination, axes, and the synaptonemal complex).  One tier of triads occurs at sites of majority (“canonical”) crossovers.  Second tier triads are more widely spaced but also exhibit interference, dependent on the same functions as canonical crossover interference.  Diverse lines of evidence suggest that second tier triads arise at sites of previously mysterious "minority" crossovers.  Finally, conserved protein remodeler Pch2/TRIP13 modulates the abundance of triad components, specifically in longer periodicity triads, dynamically in real time.  Potential roles of triad structure, mechanisms of two-tiered patterning, and the nature of minority crossovers are discussed.

## Software Dependencies
MATLAB: tested on MATLAB_R2023b running on MAC Sonoma 14.4.1

Mathematica

## Installation Guide
No specific installation is required.  Matlab functions (.m files) should be added to the MATLAB file path.
Mathematica scripts (.nb files) may need to be unzipped.

## Demo
- Data associated with the above paper has been deposited in a public repository.  It also serves the purposes of a demo dataset and can be accessed via https://doi.org/10.7910/DVN/5LEWYF.
  
- File MP_outputAnalysisOf100ControlSimulations.pdf shows an example output for simulting and analyzing data using MP_MathematicaScriptForSimulatingAndAnalyzingControlDataset.nb.  Each row is a separate simulated intensity profile.  The first column is the simulated intensity profile.  The middle column is the corresponding fourier transform in k-space.  The third column is the corresponding fourier transform in position space.

## Instructions for use
Each function contains instructions for use.
