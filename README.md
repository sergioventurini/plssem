# plssem
###### Current release: 0.3.0
###### Stata version required: at least 15.1
Stata package for Structural Equation Modeling with Partial Least Squares (PLS-SEM).

# Installation note    

To install `plssem` directly from GitHub you need to use the `github` Stata module. You can install the latest version of the `github` command by executing the following code in your Stata session:

    net install github, from("https://haghish.github.io/github/")

Then, you can install `plssem` simply using the following code in Stata:

    github install sergioventurini/plssem

Alternatively, you can install the package manually downloading the files from this GitHub repository and placing it in your Stata `PERSONAL` directory (if you don't know what this is, run the `adopath` command).

The `examples.do` file contains many examples taken from the literature as well as some simulated data examples. All the example data sets are placed in the separate `data` directory in this repository.

# Authors
Sergio Venturini, Department of Decision Sciences, Università Bocconi, Milan, Italy

E-mail: sergio.venturini@unibocconi.it

Mehmet Mehmetoglu, Department of Psychology, Norwegian University of Science and Technology, Trondheim, Norway

E-mail: mehmetm@svt.ntnu.no

# Bugs
In case you find any bug, please send us an e-mail or open an issue on GitHub.

# Citation    
You can cite the `plssem` package as:

Venturini, S., Mehmetoglu, M. (2019). plssem: A Stata Package for Structural Equation Modeling with Partial Least Squares. Journal of Statistical Software, 88(8)1--35

Paper webpage: https://www.jstatsoft.org/article/view/v088i08

GitHub repository: https://github.com/sergioventurini/plssem

# Copyright
This software is distributed under the GPL-3 license (see LICENSE file).
