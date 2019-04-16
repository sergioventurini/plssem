*!plssem_classes version 0.3.0
*!Written 15Jan2019
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 15.1

mata:

struct plssem_struct {
	real scalar n								// number of observations used in the analysis
	real scalar niter						// number of iterations to reach convergence
	real scalar converged				// indicator that the algorithm did converge
	real rowvector diff					// row vector of the max outer weights differences
	real matrix loadings				// matrix of estimated outer loadings
	real matrix xloadings				// matrix of estimated cross loadings
	real matrix path						// matrix of estimated path coefficients
	real matrix path_v					// matrix of estimated path coefficient variances
	real matrix inner_v					// diagonal matrix of inner relation variances
	real matrix scores					// matrix of estimated latent scores
	real matrix total_effects		// matrix of estimated total effects
	real matrix inner_weights		// matrix of estimated inner weights
	real matrix outer_weights		// matrix of estimated outer weights
	real matrix evo							// matrix of the outer weights evolution
	real rowvector r2						// row vector of the structural model's R2's
	real rowvector r2_a					// row vector of the structural model's adjusted
															// R2's
	real matrix M								// measurement model's adjacency matrix
	real matrix S								// structural model's adjacency matrix (if any)
	real matrix e								// measurement model's residuals
	real matrix f								// structural model's residuals, if any
	real matrix R_c							// matrix of consistent latent variable
															// correlations (PLSc)
	real matrix loadings_c			// matrix of consistent outer loadings (PLSc)
	real matrix xloadings_c			// matrix of consistent cross loadings (PLSc)
	real matrix path_c					// matrix of consistent path coefficients (PLSc)
	real matrix scores_c				// matrix of estimated latent scores (PLSc)
	real matrix total_effects_c	// matrix of estimated total effects (PLSc)
	real rowvector r2_c					// row vector of the structural model's R2's
															// (PLSc)
	real rowvector r2_a_c				// row vector of the structural model's adjusted
															// R2's (PLSc)
	real matrix e_c							// measurement model's residuals (PLSc)
	real matrix f_c							// structural model's residuals, if any (PLSc)
}

struct plssem_struct_boot {
	real scalar reps						// number of bootstrap replications
	real matrix loadings_reps		// matrix of outer loadings bootstrap replications
	real matrix loadings_bs			// matrix of outer loadings bootstrap estimates
	real matrix loadings_v			// matrix of outer loadings bootstrap variances
	real matrix xloadings_reps	// matrix of cross loadings bootstrap replications
	real matrix xloadings_bs		// matrix of cross loadings bootstrap estimates
	real matrix xloadings_v			// matrix of cross loadings bootstrap variances
	real matrix path_reps				// matrix of path coefficients bootstrap
															// replications
	real matrix path_bs					// matrix of path coefficients bootstrap estimates
	real matrix path_v					// matrix of path coefficients bootstrap variances
}

struct plssem_struct_rebus {
	real scalar niter						// number of iterations to reach convergence
	real colvector rebus_class	// column vector of the final class memberships
	real colvector touse_vec		// column vector for the observations used
	real scalar rN0							// indicator that any of the classes is empty
	real scalar rN_lte_5				// indicator that any of the classes has 5
															// observations or less
	real scalar nclass					// number of classes
	real scalar maxiter					// maximum number of REBUS-PLS iterations
	real scalar stop						// REBUS-PLS stopping criterion
}

struct plssem_struct_matrix {
	real matrix mat							// generic real matrix
}

struct plssem_struct_fimix {
	real scalar niter						// number of iterations to reach convergence
	real colvector fimix_class	// column vector of the final class memberships
	real colvector touse_vec		// column vector of the observations used
	real scalar rN0							// indicator that any of the classes is empty
	real scalar rN_lte_5				// indicator that any of the classes has 5
															// observations or less
	real scalar nclass					// number of classes
	real scalar maxiter					// maximum number of FIMIX iterations
	real scalar stop						// FIMIX stopping criterion
	real scalar restart					// number of restarts of the EM algorithm
	real rowvector rho					// row vector of mixture proportions
	struct plssem_struct_matrix ///
		colvector path						// column vector structure with path coefficients
	real colvector classfreq		// column vector of class frequencies
	real colvector ll						// column vector of incomplete loglikelihood values
	real colvector ll_c					// column vector of complete loglikelihood values
	real colvector ll_restart		// column vector of incomplete loglikelihood
															//   values across the different restarts
	real matrix ic							// real matrix of information criteria
}

struct plssem_struct_gas {
	real scalar niter						// number of iterations to reach convergence
	real colvector gas_class		// column vector of the final class memberships
	real colvector touse_vec		// column vector of the observations used
	real scalar nclass					// number of classes
	real scalar popsize					// number of individuals per generation
	real scalar ngen						// number of generations to run
}

end
