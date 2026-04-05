# aramappings 0.2.0

Updated the package due to important changes in impoted packages

- Removed functions and code related to solver glpkAPI, since it has been 
archived.
- Use new function CVXR::psolve, with solver Clarabel instead of ECOS
- The clarabel package has been updated as well, but this did not require any
updates


# aramappings 0.1.3

Removed package dependencies involving data sets. The data sets used in tests
and examples are now included in the package.


# aramappings 0.1.2

Fixed installation instructions in README and vignette


# aramappings 0.1.1

Fixed CRAN submission issues:

- Reduced size of images <br>
- Increased test tolerance in test-ara_unconstrained_l1.R


# aramappings 0.1.0

This is the first release of **aramappings**

- Adaptable Radial Axes (ARA) mapping functions:
    -   `ara_unconstrained_L2()`
    -   `ara_unconstrained_L1()`
    -   `ara_unconstrained_Linf()`
    -   `ara_exact_L2()`
    -   `ara_exact_L1()`
    -   `ara_exact_Linf()`
    -   `ara_ordered_L2()`
    -   `ara_ordered_L1()`
    -   `ara_ordered_Linf()` <br>
<br>
- Plotting functions:
    -   `draw_ara_plot_2d_standardized()` <br>
<br>
- Vignettes:
    - [Introduction to **aramappings**](https://manuelrubio.github.io/aramappings/articles/intro_to_aramappings.html)
