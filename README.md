# p5cq
Automatically exported from code.google.com/p/p5cq

Super old processing project, dumped from google code.

It's readme used to be like this:

Processing implementation for the constant Q transform

The constant Q transform method to transform a discrete Fourier transform (DFT) into a constant Q transform, where Q is the ratio of center frequency to bandwidth. This allow to have a discrete signal represented in the frequency domain by a fixed amount of bins per octave. For example setting it to 12 bins per octave each bin represents a semitone in the equally tempered musical scale, so the frequency domain representation is much more "musical" than the standard DFT. 
This implementation is based on the work by Judith C. Brown & Miller S. Puckette http://www.wellesley.edu/Physics/brown/pubs/effalgV92P2698-P2701.pdf
and taking as starting point the matlab implementation by Benjamin Blankertz http://wwwmath1.uni-muenster.de/logik/org/staff/blankertz/constQ/constQ.html


It works with the current release of processing (2.2.1)

it needs ESS audio library to work.
download it from here.
http://www.tree-axis.com/Ess/
