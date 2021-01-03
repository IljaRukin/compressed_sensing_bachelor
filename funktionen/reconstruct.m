function [result,fit_error,comp_time] = reconstruct(pic_samples,indices,dim,trafo,solver,sigma)

addpath('./solver/');
addpath('./solver/FPC_AS_v1.21\src');
addpath('./solver/NESTA_v1.1');
addpath('./solver/spgl1-1.9');
addpath('./solver/FPC_AS_v1.21/prob_gen/classes');

%%%%% normalize image: mean=0 & variance=1
mean = sum(sum(sum(pic_samples)))/numel(pic_samples);
pic_samples = pic_samples - mean;
variance = sum(sum(sum(pic_samples.^2)))/numel(pic_samples);
pic_samples = pic_samples/sqrt(variance);

%select transformations (basis)
if strcmp(trafo,'dft2')
    AA = @(x) fft2(x)/sqrt(prod(dim));
    AT = @(x) ifft2(x)*sqrt(prod(dim));
elseif strcmp(trafo,'dct2')
    AA = @(x) dct2(x);
    AT = @(x) idct2(x);
elseif strcmp(trafo,'cdft')
    AA = @(x) compress_coeff(fft2(x))/sqrt(prod(dim));
    AT = @(x) ifft2(decompress_coeff(x))*sqrt(prod(dim));
elseif strcmp(trafo,'haar')
    divisibility = 0;divider=2;
    while rem(dim,divider)==0
        divisibility = divisibility +1; divider = divider*2;
    end
    %%%
    AA = @(x) wave_2d_nonstandard(x,'haar',divisibility);
    AT = @(x) iwave_2d_nonstandard(x,'haar',divisibility);
elseif strcmp(trafo,'db04')
    divisibility = 0;divider=2;
    while rem(dim,divider)==0
        divisibility = divisibility +1; divider = divider*2;
    end
    %%%
    AA = @(x) wave_2d_nonstandard(x,'db04',divisibility);
    AT = @(x) iwave_2d_nonstandard(x,'db04',divisibility);
else
    fprintf('wrong transformation basis ! \n chose one of the following: spgl1 nesta fpcas');
end

%first guess
expanded = zeros(dim);
expanded(indices) = pic_samples;
x0 = AT(expanded);
x0 = reshape(x0, prod(dim),1);

%solvers with documented parameters
if strcmp(solver,'spgl1')
    
    opA = @(x,mode) func_basis_sparse(x,indices,dim,AA,AT,mode);
    opts = spgSetParms('verbosity',0,'bpTol',1e-6,'decTol',1e-4,'optTol',1e-6);
    sigma = sigma*sqrt(numel(pic_samples));
    
    tic
    x = spg_bpdn(opA,pic_samples,sigma,opts);
    comp_time=toc;
    
%SPGL1
%https://www.cs.ubc.ca/~mpf/spgl1/index.html
%
%Solves the basis pursuit denoise (BPDN) problem
%min_X ||X||_1  subject to  ||A X - B||_2 <= SIGMA,
%
% [x, r, g, info] = spg_bpdn(A, b, sigma, options)
%
% INPUTS
% ======
% A        m x n-transform matrix (basis)
%          can be an explicit matrix or an operator y = A(x,mode)
%          (mode=1 basis / mode=2 adjoint)
% b        m-vector of measured sample points
% sigma    noise estimate for signal.
% x0       first guess vector. [default: zero vector]
% opts
%          fid :           [ positive integer        |     1 ]     file id for log output (1=matlab screen)
%          verbosity :     [ integer: 1, 2, or 3     |     3 ]     0=quiet, 1=some output, 2=more output
%          iterations :    [ positive integer        |  10*m ]     Max. number of iterations
%          nPrevVals :     [ positive integer        |    10 ]     
%          bpTol :         [ positive scalar         | 1e-06 ]     Tolerance for identifying a basis pursuit solution.
%          lsTol :         [ positive scalar         | 1e-06 ]     
%          optTol :        [ positive scalar         | 1e-04 ]     Optimality tolerance (default is 1e-4).
%          decTol :        [ positive scalar         | 1e-04 ]     Larger decTol means more frequent Newton updates.
%          stepMin :       [ positive scalar         | 1e-16 ] 
%          stepMax :       [ positive scalar         | 1e+05 ] 
%          rootMethod :    [ 1=linear, 2=quadratic   |     2 ] 
%          activeSetIt :   [ positive integer        |   Inf ] 
%          subspaceMin :   [ 0=no, 1=yes             |     0 ] 
%          iscomplex :     [ 0=no, 1=yes, NaN=auto   |   NaN ] 
%          maxMatvec :     [ positive integer        |   Inf ] 
%          weights :       [ vector                  |     1 ] 
%          project :       [ projection function     |    @()] 
%          primal_norm :   [ primal norm eval fun    |    @()] 
%          dual_norm :     [ dual norm eval fun      |    @()] 
%
% OUTPUTS
% =======
% x        solution of the problem
% r        residual, r = b - Ax
% g        gradient, g = -A'r
% info     structure with the following information:
%          .tau     final value of tau (see sigma above)
%          .rNorm   two-norm of the optimal residual
%          .rGap    relative duality gap (an optimality measure)
%          .gNorm   Lagrange multiplier of (LASSO)
%          .stat    = 1 found a BPDN solution
%                   = 2 found a BP sol'n; exit based on small gradient
%                   = 3 found a BP sol'n; exit based on small residual
%                   = 4 found a LASSO solution
%                   = 5 error: too many iterations
%                   = 6 error: linesearch failed
%                   = 7 error: found suboptimal BP solution
%                   = 8 error: too many matrix-vector products
%          .time    total solution time (seconds)
%          .nProdA  number of multiplications with A
%          .nProdAt number of multiplications with A'
elseif strcmp(solver,'nesta')
    
    opAA = @(x) func_basis_sparse(x,indices,dim,AA,AT,1);
    opAT = @(x) func_basis_sparse(x,indices,dim,AA,AT,2);
    muf = 1e-4; %%%used for approximating l1 norm (smoothed l1 norm); mu->0 produces l1 norm mu->1 for faster solution
    sigma = sigma*sqrt(numel(pic_samples));
    
    opts = [];
    opts.Verbose = 0; %mute output
    opts.xplug = x0; %guess x
    opts.maxintiter = 5; % number of continuation steps. default is 5
    opts.TOlVar = 1e-04; % tolerance for the stopping criteria
    opts.TypeMin = 'L1'; %minimize l1 norm
    
    tic
    x = NESTA(opAA,opAT,pic_samples,muf,sigma,opts);
    comp_time=toc;
    
%Nesta
%https://statweb.stanford.edu/~candes/nesta/
%
% Solves a L1 minimization problem under a quadratic constraint using the
% Nesterov algorithm, with continuation:
%min_x || U x ||_1 s.t. ||y - Ax||_2 <= delta
%
%
%[x,niter,residuals,outputData,opts]= NESTA(opAA,opAT,pic_samples,mu,delta,opts);
%
% INPUTS
% ======
%
% AA/AT    m x n-transform matrix (basis) and it's adjoint
%          can be an explicit matrix or an operator
% b        m-vector of measured sample points
% muf      value of mu at the last continuation step.
%          used for approximating l1 norm; mu->0 produces l1 norm mu->1 for faster solution
% sigma    noise estimate for signal.
% opts
%          xplug        the first guess for the primal prox-function, and
%                       also the initial point for xk.  By default, xplug = At(b)
%          U/Ut         operator to modify the l1 norm (and it's adjoint ?)
%          normU        norm of operator U
%          MaxIntIter	number of continuation steps. default is 5
%          maxiter      max number of iterations in an inner loop. default is 10,000
%          TolVar       tolerance for the stopping criteria
%          stopTest     which stopping criteria to apply
%                       (1: stop when the relative change in the objective
%                       function is less than TolVar /
%                       2 : stop with the l_infinity norm of difference in the xk variable is less than TolVar)
%          TypeMin      Wich norm to minimize: ('L1': l1 norm / 'tv': smoothed version of the total-variation norm)
%          Verbose      how much output to display
%                       (0 or false: little output / 1 or true: output every iteration /
%                       number p greater than 1: output every pth iteration)
%          fid          id for log output
%                       (1: Matlab screen / file-id of a file opened with fopen: log output redirected to this file)
%          errFcn       Error Function or true coefficient vector norm(xk-x_true), it's value is displayed every iteration.
%          outFcn       Error Function or true coefficient vector norm(xk-x_true), it's value is stored to outputData.
%experimental features:
%          AAtinv - this is an experimental new option.  AAtinv
%                   is the inverse of AA^*.  This allows the use of a 
%                   matrix A which is not a projection, but only
%                   for the noiseless (i.e. delta = 0) case.
%                   If the SVD of A is U*S*V', then AAtinv = U*(S^{-2})*U'.
%          USV - another experimental option.  This supercedes
%                   the AAtinv option, so it is recommended that you
%                   do not define AAtinv.  This allows the use of a matrix
%                   A which is not a projection, and works for the
%                   noisy ( i.e. delta > 0 ) case.
%          USV should contain three fields: 
%          USV.U  is the U from [U,S,V] = svd(A)
%                   likewise, opts.USV.S and opts.USV.V are S and V
%                   from svd(A).  S may be a matrix or a vector.
%
% OUTPUTS
% =======
%
% xk             solution
% niter          number of iterations
% residuals      first column is the residual at every step,
%                second column is the value of f_mu at every step
% outputData	 output from opts.outFcn, if supplied.
% opts           the structure containing the options that were used    
elseif strcmp(solver,'fpcas')
    
    A = A_operator(@(x) func_basis_sparse(x,indices,dim,AA,AT,1), @(x) func_basis_sparse(x,indices,dim,AA,AT,2));
    mu = sigma; %not documented well, but it works this way !!!
    M =[];
    opts=[];
    opts.sub_mxitr = 40; opts.gtol = 1e-4; opts.gtol_scale_x = 1e-8; %most important parameters
    opts.x0 = x0; %guess x
    opts.record = -1; %no output before termination

    tic
    x = FPC_AS(prod(dim),A,pic_samples,mu,M,opts);
    comp_time = toc;
    
%FPC AS
%https://www.caam.rice.edu/~optimization/L1/FPC_AS/
%https://www.caam.rice.edu/~optimization/L1/FPC_AS/FPC_AS_v1.0_Manual.pdf
%
%Solve BPDN using an active set continuation algorithm using shrinkage and sub-optimization.
%min mu*||x||_1 + 0.5*||Ax - b||_M^2,
%
% [x,Out] = FPC_AS(n,A,b,mu,M,opts)
%
% INPUTS
% ======
%
% AA/AT    m x n-transform matrix (basis) and it's adjoint
%          can be an explicit matrix or an operator
% b        m-vector of measured sample points
% M        positive definite matrix or the empty matrix []. (M=[]: sets as an identity matrix)
%          Note: maximum eigenvalue of A' M A should be <=1
%          automaticlly scale matrix: opts.scale_A = 1 as theta*A by factor theta
% mu       l1 regularization parameter (tradeoff between l1 and accuracy)
% opts
%           manual says, that the algorithm is highly tuned and shouldn't be modified.
%           recomended parameters for high accuracy on noisy data are:
%           opts.sub_mxitr = 10; opts.gtol = 1e-3; opts.gtol_scale_x =1e-6;
%   'x0':   initial solution
%           default: - , valid range: [-Inf, Inf]
%   ’init’: methods of initialization, integer
%           default: 2, valid range:{0,1,2}
%   'xs':   exact solution whose purpose is only for comparison
%           default: - , valid range: [-Inf, Inf]
%   ’tol_eig’: tolerance for eigenvalues
%           default: 0.0001, valid range: [0, 1]
%   ’scale_A’: if 1 scale the matrix A so that max of eigs(A*AT) equals 1
%           default: 0, valid range:{0, 1}
%   ’eigs_mxitr’: max number of iterations for eigs(A*AT), integer
%           default: 20, valid range: [1, 100]
%   ’eps’: machine accuarcy
%           default: 1e-16, valid range: [0, 1]
%   ’zero’: minimal magnitude of x
%           default: 1e-08, valid range: [0, 1e+10]
%   ’dynamiczero’: set the thresholding level dynamically or not, integer
%           default: 0, valid range:{0, 1}
%   'minK': an estimate of the number of nonzero components in optimal solution. 
%           default: m/2, valid range: [1, n]
%   ’tauD’: a parameter for shrinkage
%           default: min(1.999,-1.665*m/n + 2.665), valid range: [0, 100000]
%   ’tau_min’: minimal tau
%           default: 0.0001, valid range: [0, 100000]
%   ’taumax’: minimal tau
%           default: 1000, valid range: [0, 100000]
%   'mxitr': ax number of iterations, integer
%           default: 1000, valid range: [1, 100000]
%   'gtol': termination criterion on ``crit2'', the maximum norm of sub-gradient
%           default: 1e-08, valid range: [0, 1]
%   'gtol_scale_x': termination criterion on ``crit2'' scaled by max(norm(x), 1). 
%           default: 1e-12, vaild range: [0, 1]
%   ’f_rel_tol’: Tolerance on the relative changes of function value
%           default: 1e-12, valid range: [0, 1]
%   'f_value_tol': Tolerance on the optimal objective value, stop if f(x) <= f_value_tol
%           default: 0, vaild range: [0, inf]
%   ’lsmxitr’: max number of iterations of of line search subroutine, integer
%           default: 5, valid range: [2, 100]
%   ’gamma’: a parameter for the nonmonotone line search
%           default: 0.85, valid range: [0, 1]
%   ’c’: a parameter for Armijo condition
%           default: 0.001, valid range: [0, 1
%   ’beta’: a parameter for decreasing the step size in the nonmonotone line search
%           default: 0.5, valid range: [0, 1
%   ’eta’: a parameter for decreasing mu
%           default: 0.5, valid range: [0, 1]
%   ’eta_rho’: a parameter for decreasing eta
%           default: 0.5, valid range: [0, 1]
%   ’eta_min’: minimal eta
%           default: 0.001, valid range: [0, 1]
%   ’eta_max’: maximal eta
%           default: 0.8, valid range: [0, 1]
%   ’max_itr_mu’: control the decreasing of eta, iteration number between the changes of mu, integer
%           default: 3, valid range: [1, 100]
%   ’sub_mxitr’: max number of iterations for doing sub-optimization, integer
%           default: 80, valid range: [1, 100000]
%   'lbfgs_m': storage number of L-BFGS. 
%           default: 5, valid range: [1, 100]
%   ’kappa_g_d’: tolerance for checking whether do sub-optimization or not
%           default: 10, valid range: [1, 1000]
%   ’kappa_rho’: a parameter for increasing kappa_g_d
%           default: 10, valid range: [1, 1000]
%   ’tol_start_sub’: Tolerance of starting sub-optimization
%           default: 1e-06, valid range: [0, 1]
%   ’min_itr_shrink’: min number of iterations between two sub-optimization, integer
%           default: 3, valid range: [1, 1000]
%   ’max_itr_shrink’: max number of iterations between two sub-optimization, integer
%           default: 20, valid range: [1, 1000]
%   'record': print level, -1=quiet, 0=some output, 1=more output. 
%           default: 0 valid range: {-1, 0, 1}
%   'PrintOptions': print options, 0=quiet, 1=output 
%           default: 0, valid range: {0, 1}
%
%%%%%not documented in manual ???
%   'sub_opt_meth': choice of sub-optimization methods.
%           default: 'lbfgs', valid values: {'lbfgs','lbfgsb','pcg'};
%   'zero': a lower bound of minimal magnitude of nonzero components of optimal solution
%           default: 1e-08, valid range: [0, 1e+10]
%   'dynamic_zero': on/off switch for setting `zero' dynamically
%           default: 0, valid range: {0, 1}
%
% OUTPUTS
% =======
%
% x             solution of the problem
% Out           structure eith the following information
%    .cpu       -the total amount of CPU time 
%    .exit      -exit status, 1='normal', 10='exceed max iterations'
%    .mesg      -a string describe the detailed exit status
%    .itr       -the number of iterations
%    .f         -exit function value, mu*||x||_1 + 0.5*||Ax - b||_M^2
%    .nrm1x     -exit l1 norm ||x||_1
%    .rNorm     -exit l2 discrepancy term ||Ax - b||_M
%    .g         -exit gradient of 0.5*||Ax - b||_M^2
%
%    .zero      -final thresholding value for computing termination criteria 
%    .crit2     -violation of optimality, norm(g(T)-mu,'inf'), where T = union(nz_x, z_xa),
%                nz_x = x>Out.zero, the set of nonzero components whose magnitude of x 
%                       are larger than Out.zero
%                z_xa = ~nz_x & (||g|-mu| > Out.zero); the set of nonzero components 
%                       whose magnitude of |g|-mu are larger than Out.zero
%
%    .nnz_x     -number of elements such that x_i > Out.zero
%    .nCont     -number of continuation steps taken
%    .nSubOpt   -number of sub-optimization problems solved
%    .nProdA    -total number of A*x    
%    .nProdAt   -total number of A'*x
%    .nfe       -number of A*x called by shrinkage
%    .nge       -number of A'*x called by shrinkage
%    .nfe_sub   -number of A*x called by sub-optimization
%    .nge_sub   -number of A'*x called by sub-optimization
%
%    .opts      -options used
%
%    --> The following fields are available if an optimal solution opts.xs 
%    is provided in the input. opts.xs is compared to x after x is truncated
%    in the way that x_i=0 if x_i <Out.zero. The comparison results are given in
%
%    .sgn       -number of nonzero components of x that have different sign compared to those of opts.xs,
%    .miss      -number of (missed) components that are zero in x but nonzero in opts.xs
%    .over      -number of (overshot) components that are nonzero in x but zero in opts.xs
%
%    They are computed as 
%         jnt  = union(nz_xs, nz_x); nz_xs = abs(xs) > opts.eps;
%         sgn  = nnz(sign(x(jnt))~=sign(xs(jnt))); 
%         miss = nnz(nz_xs&(~nz_x))
%         over = nnz(nz_x&(~nz_xs))
%
else
    fprintf('wrong solver name: \n chose \"spgl1\" or \"nesta\" or \"fpcas\" as solver')
end

%rearrange output
x = reshape(x,dim);
%transform backwards (reconstruct whole image from coefficients)
pic_reconstructed = real(AT(x));

%calculate fit error
fit_error = sum(sum(sum((pic_samples-pic_reconstructed(indices)).^2)))/sum(sum(sum(pic_samples)));

%reverse normalization
result = pic_reconstructed*sqrt(variance) + mean;
result = real(result);

end