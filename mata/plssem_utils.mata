*!plssem_utils version 0.3.0
*!Written 15Jan2019
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

version 15.1

mata:

real matrix scale(real matrix X, |real scalar biased, real rowvector scale, ///
	real rowvector center)
{
	/* Description:
		 ------------
		 Function that standardizes the columns of a real matrix X
	*/
	
	/* Arguments:
		 ----------
		 - X			--> real matrix whose columns must be scaled
		 - biased	--> (optional) real scalar; if not zero the biased sample
									variances are computed
		 - scale	--> (optional) real rowvector containing the column scales
									(after centering)
		 - center	--> (optional) real rowvector containing the column centers
	*/
	
	/* Returned value:
		 ---------------
		 - Xs 		--> real matrix with columns that have been centered and/or scaled
	*/

	real matrix cen, sc, sc_inv, Xc, Xcs
	real scalar n, factor
	
	n = rows(X)
	if (center == J(1, 0, .)) {
		cen = J(n, 1, mean(X))
	}
	else {
		cen = J(n, 1, center)
	}
	if (scale == J(1, 0, .)) {
		sc = sqrt(diag(variance(X)))
	}
	else {
		sc = diag(scale)
	}
	if (biased == .) {
		biased = 0
	}
	if (biased) {
		factor = sqrt((n - 1)/n)
		sc = sc :* factor
	}
	Xc = X - cen
	sc_inv = invsym(sc)
	Xcs = Xc * sc_inv

	return(Xcs)
}

real matrix unscale(real matrix Xs, real rowvector scale, ///
	real rowvector center)
{
	/* Description:
		 ------------
		 Function that unstandardizes the columns of a real matrix Xs
	*/
	
	/* Arguments:
		 ----------
		 - Xs			--> real matrix whose columns are standardized
		 - scale	--> real rowvector containing the column scales
									(after centering)
		 - center	--> real rowvector containing the column centers
	*/
	
	/* Returned value:
		 ---------------
		 - X 			--> real matrix with columns that have been unstandardized
	*/

	real matrix cen, sc, Xunscaled, X
	real scalar n
	
	n = rows(Xs)
	cen = J(n, 1, center)
	sc = diag(scale)
	Xunscaled = Xs * sc
	X = Xunscaled + cen

	return(X)
}

real rowvector sd(real matrix X, |real scalar biased)
{
	/* Description:
		 ------------
		 Function that computes the standard deviation of columns of X
	*/
	
	/* Arguments:
		 ----------
		 - X			--> real matrix of data
		 - biased	--> (optional) real scalar; if not zero the biased sample variance
									is computed
	*/
	
	/* Returned value:
		 ---------------
		 - sc --> real rowvector containing the X column standard deviations
	*/
	
	real rowvector sc
	real scalar n, bias, factor
	
	n = rows(X)

	sc = sqrt(diagonal(variance(X)))'
	
	if (biased == . | biased == 0) {
		bias = 0
	}
	if (bias) {
		factor = sqrt((n - 1)/n)
		sc = sc:*factor
	}

	return(sc)
}

real colvector which(real matrix X, |string scalar minmax)
{
	/* Description:
		 ------------
		 Function that determines the index of the (first) minimum or maximum of a
		 numeric matrix
	*/
	
	/* Arguments:
		 ----------
		 - X			--> real matrix of data
		 - minmax	--> (optional) string scalar for minimum ("min") or maximum
									("max"); deafult is "min"
	*/
	
	/* Returned value:
		 ---------------
		 - index --> real vector containing the indexes at which min or max occur
	*/
	
	real scalar i, j, N, w
	real colvector index
	
	pragma unset i		// this is needed to avoid useless warning messages
	pragma unset w		// this is needed to avoid useless warning messages

	N = rows(X)
	
	if (args() == 1) {
		minmax = "min"
	}

	index = J(0, 1, .)
	for (j = 1; j <= N; j++) {
		if (minmax == "min") {
			minindex(X[j, .], 1, i, w)
		}
		else if (minmax == "max") {
			maxindex(X[j, .], 1, i, w)
		}
		else {
			_error(3000, "the 'minmax' argument value is unrecognized")
		}
		index = (index \ i[1])
	}
	
	return(index)
}

real matrix rdirichlet(real scalar n, real rowvector alpha)
{
	/* Description:
		 ------------
		 Function that generates random deviates from a Dirichlet distribution
	*/
	
	/* Arguments:
		 ----------
		 - n 			--> real scalar providing the number of deviates to draw
		 - alpha	--> real rowvector containing the Dirichlet parameters
	*/
	
	/* Returned value:
		 ---------------
		 - out		--> real colvector containing the generated Dirichlet deviates
	*/
	
	real matrix out
	
	out = rgamma(n, 1, alpha, 1)
	
	return(out :/ rowsum(out))
}

real colvector euclidean_dist(real matrix X, real scalar idx)
{
	/* Description:
		 ------------
		 Function that computes the dissimilarities among a set of observations
	*/
	
	/* Arguments:
		 ----------
		 - X					--> real matrix containing the data
		 - idx 				--> real scalar providing the row index of X for which to
											calculate the dissimilarities
	*/
	
	/* Returned value:
		 ---------------
		 - D					--> vector of dissimilarities
	*/
	
	real colvector D
	real scalar n
	
	n = rows(X)
	D = sqrt(rowsum((X - J(n, 1, 1)*X[idx, .]):^2))
	
	return(D)
}

real scalar prod(real vector v)
{
	/* Description:
		 ------------
		 Function that computes the product of a vector's elements
	*/
	
	/* Arguments:
		 ----------
		 - v					--> vector of real elements
	*/
	
	/* Returned value:
		 ---------------
		 - res 				--> product of the v vector's elements
	*/

	real scalar i, nr, res
	
	nr = length(v)
	
	res = 1
	for (i = 1; i <= nr; i++) {
		res = res*v[i]
	}
	
	return(res)
}

void print_struct_matrix(struct plssem_struct_matrix colvector A)
{
	/* Description:
		 ------------
		 Function that prints the elements of a struct of plssem_struct_matrix
		 objects
	*/
	
	/* Arguments:
		 ----------
		 - A					--> struct of plssem_struct_matrix objects in a colvector
											layout
	*/
	
	/* Returned value:
		 ---------------
		 - void 			--> no object returned
	*/

	real scalar i, j, nr, nc
	
	nr = rows(A)
	nc = cols(A)
	
	for (i = 1; i <= nr; i++) {
		for (j = 1; j <= nc; j++) {
			A[i, j].mat
		}
	}
}

void cleanup()
{
	/* Description:
		 ------------
		 Function that cleans up Mata from temporary objects 
	*/
	
	/* Arguments:
		 ----------
		 - void				--> no argument passed
	*/
	
	/* Returned value:
		 ---------------
		 - void 			--> no object returned
	*/

	string colvector names
	real scalar i
	
	names = direxternal("__*")

	for (i = 1; i <= rows(names); i++) {
		rmexternal(names[i])
	}
}

end
