______________________________
# Hydroelastic Stability Analysis of an Inverted T Hydrofoil

Matlab/Octave code for the calculation of Eignemodes / Eigenvectors of an inverted "T" hydrofoil at high speed.
______________________________
## description

This code uses a lumped-mass Finite Element Method tightly coupled with a Potential flow solution for the fluid in order to calculate the resonance frequency for linear perturbations to the foil while moving at high speed through the water.

It is quasi quasi 3d
- structure: 2d bending coupled with 1d rotation
- fluid: 2d aerfoil slices (no 3d effects)

You can validate the foil structural dynamics by removing the fluid and testing it in various configuration (cantilever beam, free beam, etc.)

______________________________
## Forcing frequency/phase-shift diagram - getfr.m

I think it is just a pretty standard Bode diagram.
Its plotting the frequency response of the FSI system as a whole.
If you force the system with some external force at a frequency (the forcing frequency)...
then the system will respond with a particular amplitude that is some ratio of the force you are applying (the amplitude ratio)... and it will also respond with a phase that is different from the forcing frequency (the phase shift)

When the system resonates in particular modes... then the phase shift will be in-line with the external force (they will start to resonate together)... and the amplitude ratio will spike... that is why you see those spikes in both the amplitude ratio and phase diagrams... those are the points where the foil will resonate with any external force...

The type and nature of the external force is not important... because you can assume that in a practical application there will always be some external force being applied to the system (small disturbances or pertrubations in the air/water it doesn't matter)... what is important is that the system will pick up on those and resonate with them at particular frequencies

If it were very important... then you could generate a kind of frequency/amplitude diagram of those external forces and then you could get an idea of the likelihood and severity of the FSI system (foil) being influenced by those frequencies and the magnitude of them...

So for example, if the foil resonates in mode 1 at say 1 hz... and if the amplitude ratio is say 3... then you could make an engineering assessment of the environment this foil will operate within... is it likely to encounter external forces at a frequency of 1Hz and, if so, how much energy do those forces have.  So the foil may be cutting through the ocean at 30 knots (15 m/s) through ocean swell of wavelength around 15m.... in this case there is a rather large external force being encoutered at around 1Hz... we could estimate the magnitude of that external force on the foil (or parts of the foil) assuming it is static.  So maybe the external wave forces on the static foil are say 100kN... then assuming that will resonate with the foil at mode 1... then due to the amplitude ratio at that resonance point (amplitude ratio of 3) ... then you could use that to apply your engineering factors of safety... you have to at least triple the static forces to account for possible resoonance at that frequency





