%
% Perform Armijo Line Search
% the input are
%   fun  = the function to minimize with the derivative of the function
%   xk   = starting point
%   d    = search direction
%
%  optional extra arguments
%
%    tau   = 0.5    reduction parameter for Armijo algorithm
%    c1    = 10-3   parameter for Armijo test
%    f0    = f(0)   for non monotone line search
%    lmin  = 1E-8   minimum step length
%    lmax  = 1E+8   maximum step length
%
function lambda = ArmijoLineSearch( fun, xk, f0, Df0, d, varargin )

  opts = { 'c1', 'tau', 'lmin', 'lmax', 'minDF0', 'f0' } ;
  defs = { 1E-3,   0.5,  1E-20,   1E+8,    -1E+6,  f0 } ;

  OPTS = parseArgs( opts, defs, varargin{:} ) ;

  lambdaMin = OPTS.lmin   ;
  lambdaMax = OPTS.lmax   ;
  tau       = OPTS.tau    ;
  c1        = OPTS.c1     ;
  f0        = OPTS.f0     ;
  minDF0    = OPTS.minDF0 ;

  if c1 <= 0 || c1 >=1
    fprintf( 1, 'ArmijoLineSearch: Bad parameter c1 = %g\n', c1 ) ;
    lambda = NaN ;
    return ;
  end

  if any(isnan(xk)) && any(isinf(xk))
    disp(xk) ;
    fprintf( 1, 'ArmijoLineSearch: Bad initial point!\n' ) ;
    lambda = NaN ;
    return ;
  end

  if isnan(Df0) || Df0 >= 0
    fprintf( 1, 'ArmijoLineSearch: Bad search direction Df0 = %g\n', Df0 ) ;
    lambda = NaN ;
    return ;
  end

  c1Df0   = c1 * Df0 ;
  dumpMin = 2 ;
  dumpMax = 100 ;

  % if derivative is too strong use Simple Armijo
  if Df0 < minDF0
    lambda = ArmijoSimpleSearch( fun, xk, f0, Df0, d, varargin{:} ) ;
    return ;
  end

  lambda  = 1 ;
  fLambda = feval( fun, xk + lambda * d ) ;

  if fLambda <= f0 + lambda*c1Df0
    % forward search
    while fLambda <= f0 + lambda*c1Df0
      if lambda > lambdaMax
        fprintf( 1, 'ArmijoLineSearch forward search failed! fLambda=%g lambda=%g\n', fLambda, lambda ) ;
        lambda = rev*lambda ;
        return ;
      end
      lambda  = lambda*2 ;
      fLambda = feval( fun, xk + lambda * d ) ;
    end ;
    lambda = lambda/2 ;
  else
    % backward search
    lambdap = lambda ;
    while fLambda > f0 + lambda*c1Df0
      if lambda < lambdaMin
        fprintf( 1, 'ArmijoLineSearch backward search failed! fLambda=%g lambda=%g\n', fLambda, lambda ) ;
        lambda = rev*lambda ;
        return ;
      end
      if lambda > lambdap
        lambdaTmp = cubic( f0, Df0, fLambda, lambda, fLambdap, lambdap ) ;
      else
        lambdaTmp = quadratic( f0, Df0, fLambda, lambda ) ;
      end ;
      lambdap  = lambda ;
      fLambdap = fLambda ;
      lambda   = max( min( lambdaTmp, lambda/dumpMin ), lambda/dumpMax ) ;
      fLambda  = feval( fun, xk + lambda * d ) ;
    end ;
  end ;
end
%                         _           _   _      
%    __ _ _   _  __ _  __| |_ __ __ _| |_(_) ___ 
%   / _` | | | |/ _` |/ _` | '__/ _` | __| |/ __|
%  | (_| | |_| | (_| | (_| | | | (_| | |_| | (__ 
%   \__, |\__,_|\__,_|\__,_|_|  \__,_|\__|_|\___|
%      |_|                                       
%      
function alpha = quadratic( f0, Df0, fp, p )
  alpha = Df0 * p^2 / ( 2*(f0+Df0*p-fp) ) ;
end
%              _     _      
%    ___ _   _| |__ (_) ___ 
%   / __| | | | '_ \| |/ __|
%  | (__| |_| | |_) | | (__ 
%   \___|\__,_|_.__/|_|\___|
% 
function alpha = cubic( f0, Df0, fl, l, fp, p )
  % f0 + Df0 * x + a * x^2 + b * x^3
  bf1 = fl - (f0+l*Df0) ;
  bf2 = fp - (f0+p*Df0) ;
  bf3 = (p*l)^2*(p-l) ;
  a   = (p^3 * bf1 - l^3 * bf2) / bf3 ;
  b   = (l^2 * bf2 - p^2 * bf1) / bf3 ;

  if abs(b) < 1E-10
    % points disposed near a parabola
    alpha = -Df0 / (2*a) ;
  else
    delta = a^2 - 3*b*Df0 ;
    if delta < 0
      % complex root, using quadratic again
      alpha = quadratic( f0, Df0, fl, l ) ;
    else
      % regular case, using cubic
      alpha = (sqrt(delta)-a)/(3*b) ;
    end;
  end;
end
