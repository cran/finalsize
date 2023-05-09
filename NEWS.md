# finalsize 0.2

This is the second release of _finalsize_, and includes:

1. The package's C++ code for the solver algorithms has been moved from source files under `src/` into headers under `inst/include/`
2. The package includes a package header under `inst/include` called `finalsize.h`, which allows other Rcpp packages to link to _finalsize_ and reuse the solver algorithms provided here.
3. Three helper functions have been added to:
    - Calculate the effective basic reproductive number $R_{eff}$ in a population with heterogeneous social contacts and susceptibility to infection
    - Convert from the basic reproductive number provided by the user $R_0$ to the transmission rate (denoted by $\lambda$), given the population's social contacts and susceptibility structure
    - Convert from a user-provided transmission rate to the basic reproductive number given the population's social contacts and susceptibility structure
4. Two new vignettes have been added which are intended to serve as contextual background information for users:
    - A guide to constructing susceptibility matrices
    - A comparison with a simple susceptible-infectious-recovered epidemic model
5. Package infrastructure has been modified:
    - The package now specifies C++17 as standard
    - The `Readme.Rmd` uses auto-rendering of the package and repository name to `Readme.md` via an updated `render-readme` Github Actions workflow
    - The Cpplint workflow now includes linting and checking using Cppcheck for header files
6. Updated `NEWS.md` file to track changes to the package.

# finalsize 0.1.0.9000

* Added a `NEWS.md` file to track changes to the package.

# finalsize 0.1

Initial release of _finalsize_, an R package to calculate the final size of an epidemic in a population with demographic variation in social contacts and in susceptibility to infection.

This release includes:
1. A choice of equation solver functions.
2. 100% code coverage,
3. A basic usage vignette, and two advanced vignettes,
4. Example data from the POLYMOD dataset obtained using the _socialmixr_ R package,
5. Workflows to render the vignettes and README as a website.