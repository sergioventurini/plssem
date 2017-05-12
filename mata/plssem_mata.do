* set trace on

	/* TO DO:
		 1. 
	*/

	/* ISSUES:
		 1. 
	*/

*!plssem_mata version 0.1.0
*!Written 09May2017
*!Written by Sergio Venturini and Mehmet Mehmetoglu
*!The following code is distributed under GNU General Public License version 3 (GPL-3)

capture mata: mata drop scale()
version 10
mata:
real matrix scale(real matrix X, |real rowvector center, real rowvector scale)
{
	/* Description:
		 ------------
		 Function that standardizes the columns of a real matrix X.
	*/
	
	/* Arguments:
		 ----------
		 - X			--> real matrix whose columns must be scaled.
		 - center	--> (optional) real rowvector containing the column centers.
		 - scale	--> (optional) real rowvector containing the column scales
									(after centering).
	*/
	
	/* Returned value:
		 ---------------
		 - Xs --> real matrix with columns that have been centered and/or scaled.
	*/

	real matrix cen, sc, Xc, Xcs
	real scalar n
	
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
	Xc = X - cen
	sc_inv = luinv(sc)
	Xcs = Xc * sc_inv

	return(Xcs)
}
end

capture mata: mata drop sd()
version 10
mata:
real rowvector sd(real matrix X)
{
	/* Description:
		 ------------
		 Function that computes the standard deviation of columns of X.
	*/
	
	/* Arguments:
		 ----------
		 - X --> real matrix of data.
	*/
	
	/* Returned value:
		 ---------------
		 - sd --> real rowvector containing the X column standard deviations.
	*/
	
	return(sqrt(diagonal(variance(X)))')
}
end

capture mata: mata drop which()
version 10
mata:
real vector which(real matrix X, |string scalar minmax)
{
	/* Description:
		 ------------
		 Function that determines the index of the (first) minimum or maximum of a
		 numeric matrix.
	*/
	
	/* Arguments:
		 ----------
		 - X			--> real matrix of data.
		 - minmax	--> (optional) string scalar for minimum ("min") or maximum
									("max"); deafult is "min"
	*/
	
	/* Returned value:
		 ---------------
		 - index --> real vector containing the indexes at which min or max occur.
	*/
	
	real scalar i, j, N, w
	real vector index
	
	pragma unset i // this is needed to avoid useless warning messages
	pragma unset w // this is needed to avoid useless warning messages

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
		index = (index \ i)
	}
	
	return(index)
}
end

capture mata: mata drop rebus_cm()
version 10
mata:
real colvector rebus_cm(real matrix e, real matrix f, real matrix loads, ///
	real rowvector r2)
{
	/* Description:
		 ------------
		 Function that computes the REBUS-PLS closeness measure.
	*/
	
	/* Arguments:
		 ----------
		 - e 			--> real matrix containing the PLS-SEM measurement model residuals.
		 - f			--> real matrix containing the PLS-SEM structural model residuals.
		 - loads	--> real matrix containing the PLS-SEM communalities.
		 - r2			--> real rowvector containing the PLS-SEM R-squares.
	*/
	
	/* Returned value:
		 ---------------
		 - cm			--> real colvector containing the closeness measures for all
									observations.
	*/
	
	real scalar N
	real colvector cm, left, right
	real matrix com, e2, f2, r2_mat
	
	N = rows(e)
	
	e2 = e:^2
	com = diag(rowsum(loads:^2))
	left = rowsum(e2 * luinv(com), .)
	left = left / colsum(left)

	f2 = f:^2
	r2_mat = diag(r2)
	right = rowsum(f2 * luinv(r2_mat), .)
	right = right / colsum(right)

	cm = (N - 2)*sqrt(left :* right)
	
	return(cm)
}
end

capture mata: mata drop rebus_gqi()
version 10
mata:
real scalar rebus_gqi(real matrix ind, real matrix lat, real matrix e, ///
	real matrix f, real colvector cl)
{
	/* Description:
		 ------------
		 Function that computes the Group Quality Index (GQI) for a REBUS solution.
	*/
	
	/* Arguments:
		 ----------
		 - ind		--> real matrix containing the indicator values stacked by REBUS
									group.
		 - lat		--> real matrix containing the latent variable values stacked by
									REBUS group.
		 - e 			--> real matrix containing the PLS-SEM measurement model residuals.
		 - f			--> real matrix containing the PLS-SEM structural model residuals.
		 - cl			--> real colvector containing the REBUS classification.
	*/
	
	/* Returned value:
		 ---------------
		 - gqi		--> real scalar providing the Group Quality Index (GQI) value.
	*/
	
	real scalar N, Ncl, gqi
	real colvector Nk, left, right
	real matrix e2, f2, e2_k, f2_k, ind_k, lat_k, ind_k_den, lat_k_den
	
	N = rows(ind)
	Ncl = max(cl)
	e2 = e:^2
	f2 = f:^2
	Nk = J(0, 1, .)
	left = J(0, 1, .)
	right = J(0, 1, .)
	for (k = 1; k <= Ncl; k++) {
		Nk = (Nk \ rows(cl[selectindex(cl :== k), 1]))
		e2_k = e2[selectindex(cl :== k), .]
		f2_k = f2[selectindex(cl :== k), .]
		ind_k = ind[selectindex(cl :== k), .]
		lat_k = lat[selectindex(cl :== k), .]
		ind_k_den = diagonal(crossdev(ind_k, mean(ind_k), ind_k, mean(ind_k)))
		lat_k_den = diagonal(crossdev(lat_k, mean(lat_k), lat_k, mean(lat_k)))
		left = (left \ (1 - mean(colsum(e2_k)' :/ ind_k_den)))
		right = (right \ (1 - mean(colsum(f2_k)' :/ lat_k_den)))
	}
	// sqrt((Nk' * (left :* right))/N) // this is the calculation according
																		 // to Sanchez's plspm R package

	gqi = sqrt((Nk' * left)*(Nk' * right))/N	// this is the calculation according
																						// to Trinchera's PhD thesis (2007)

	return(gqi)
}
end

********************************************************************************

mata: mata mosave scale(), dir(PERSONAL) replace
mata: mata mosave sd(), dir(PERSONAL) replace
mata: mata mosave which(), dir(PERSONAL) replace
mata: mata mosave rebus_cm(), dir(PERSONAL) replace
mata: mata mosave rebus_gqi(), dir(PERSONAL) replace
