# cutoff

This is an R package that implements the method used in Trang et al. (2015). It fits a finite mixture model (Schlattman 2009)
to a bimodal distribution using the Expectation-Maximization algorithm (Do and Batzoglou 2008). Confidence interval of the
mixture parameter is found using the method of Oakes (1999). The fitted finite mixture model is then used to calculate a cutoff
value that separates the data in two groups, given a type-1 error to belong to one of the two modes.

### References
* Do, C. B., and S. Batzoglou. 2008. What is the expectation maximization algorithm? Nat Biotechnol 26:897–899.
* Oakes, D. 1999. Direct calculation of the information matrix via the EM algorithm. J R Statist Soc B 61:479–482.
* Schlattmann, P. 2009. Medical Applications of Finite Mixture Models. Springer Verlag.
* Trang, N. V., M. Choisy, N. T. Nakagomi, N. T. U. Chinh, Y. H. Doan, T. Yamashiro, J. E. Bryant, et al. 2015. Determination of
cut-off cycle threshold values in routine RT–PCR assays to assist differential diagnosis of norovirus in children hospitalized
for acute gastroenteritis. Epidemiol Infect. In press.

# Installation
To install cutoff2, run the following code

``` devtools::install_github("DRJP/cutoff/cutoff2@penalties", build_manual=TRUE, build_vignettes=TRUE) ```
