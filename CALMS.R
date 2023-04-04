# It takes as input:
#   - Y: a response vector
#   - phi: a design matrix
#   - lam: a tuning parameter controlling the degree of sparsity in the solution
#   - etaini: an initial value for the optimization procedure
#   - admmAbsTol: absolute tolerance for the ADMM algorithm
#   - admmRelTol: relative tolerance for the ADMM algorithm
#   - Gamma: a parameter controlling the weight assigned to non-zero coefficients in the penalty term
# The output is a list containing:
#   - one.step.x: the solution obtained after one step of the CMAL algorithm
#   - x: the final solution obtained after the CMAL algorithm
#   - t2-t1: the time taken for the first step of the algorithm (multi_signal_lasso)
#   - t2-t1+t4-t3: the total time taken for the algorithm

# Some Remarks:
#   multi_signal_lasso is called to perform the first step of the CMAL algorithm. Its inputs are Y, phi, lam, etaini, admmAbsTol and admmRelTol. The output is x, which is stored in the variable x.
#   The weights for the penalty term are computed using apply and cbind functions. The resulting vector wei contains the minimum values of the abs(x) and abs(x-1) columns, where x is the solution obtained from the first step of the CMAL algorithm.
#   If all elements of wei are zero, then the solution is already binary and no further processing is needed. The function returns a list with one.step.x and x equal to the same binary solution.
#   If there are non-zero elements in wei, then the indices ind and indi are computed based on the values of wei.
#   Next, the design matrix phi2 is computed by weighting the columns of phi corresponding to non-zero elements in wei.
#   Depending on whether ind contains a single or multiple indices, the response vector Y2 is computed using matrix operations.


# The function constrained_MDSL is called to perform the second step of the CMAL algorithm. Its inputs are:
#   - Y2: the modified response vector after removing the effect of the zero-weighted design matrix columns
#   - phi2: the modified design matrix
#   - lam: the tuning parameter controlling the degree of sparsity in the solution
#   - rho: a parameter for the ADMM algorith
#   - CMALweights: the weights for the penalty term
#   - etaini: an initial value for the optimization procedure
#   - admmAbsTol: absolute tolerance for the ADMM algorithm
#   - admmRelTol: relative tolerance for the ADMM algorithm
# Finally, the resulting solution x2 is combined with the first step solution one.step.x to obtain the final solution. The output is a list containing one.step.x, x, and two time measurements t2-t1 and t2-t1+t4-t3.





