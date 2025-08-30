# NCS
Library for numerically solving BCS gap equation for noncentrosymmetrical 1D superconductor and finding critical temperature, chemical potential and nodal points.

# 
This was orginally part of my master thesis writen under supervision of doc. hab. Grzegorz Haran. 

# Results

[![Plot Preview](plots/tc_1.1_1_2-eps-converted-to.pdf)](plots/tc_1.1_1_2-eps-converted-to.pdf)

# TO DO
1. Still a major refactoring is in order.
2. Upgrade setting up of GL quadrature.
3. Some Cuda kernel optimizations.
4. Add documentation.
5. Present derivation of the equation on english.
6. Derrive more general equation to support more then 2 bands or plasmonic component.
7. Include corrections to FD distribution.
8. Add regression tests
9. Bring back other d(k) vector directions
10. Revisit results and universality