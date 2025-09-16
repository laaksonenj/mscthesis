# Master's thesis

This repository contains the source files for my master's thesis titled "L2 Convergence of the p- Finite Element Method for the 2D Poisson Problem with a Dirac Delta Load".

The thesis considers the partial differential equation
```math
-\Delta u = \delta_{x_0},
```
where $\delta_{x_0}$ is the [Dirac delta function](https://en.wikipedia.org/wiki/Dirac_delta_function), over an arbitrary two-dimensional polygonal convex domain.
The equation is known to have a solution, and the research question is whether the p-version of the [finite element method](https://en.wikipedia.org/wiki/Finite_element_method)
can be used to approximate the solution. The research question is answered affirmatively in two parts.
The first part is theoretical in nature, and it provides a theoretical bound for the approximation error.
The second part is computational where the error is computed exactly and compared to the theoretical bound.
The results show that the theoretical error bound matches closely with the observed bound in the general case.

The computational results were obtained with a custom-made finite element solver whose source code can be found under the `src` directory.
The solver is implemented in C++, and it depends on a few open-source third-party libraries.
The Eigen library is used for all linear algebra computations.
The GMP library is used to provide the underlying multiprecision data types to minimize the effects of floating-point errors.
OpenMP is used to parallelize a couple of embarrassingly parallel loops.
