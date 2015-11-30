%
% Perform Wolfe Line Search
%
% the input are
%   fun  = the function to minimize with the derivative of the function
%   xk   = starting point
%   d    = search direction
%
%  optional extra arguments
%
%    tau   = 0.5    reduction parameter for Armijo algorithm
%    c1    = 10E-3  parameter for Wolfe test
%    c2    = 0.5    parameter for Wolfe test
%    f0    = f(0)   for non monotone line search
%    lmin  = 1E-8   minimum step length
%    lmax  = 1E+8   maximum step length
%

function l = WolfeLineSearch( fun, xk, f0, Df0, d, varargin )

  opts = { 'c1', 'c2', 'lepsi', 'lmin', 'lmax', 'minDF0', 'f0' };
  defs = { 1E-3,  0.49,  1E-10,  1E-20,   1E+8, -1E+6,    f0 } ;

  OPTS = parseArgs( opts, defs, varargin{:} ) ;

  lMin    = OPTS.lmin ;
  lMax    = OPTS.lmax  ;
  lepsi   = OPTS.lepsi ;
  c1      = OPTS.c1 ;
  c2      = OPTS.c2 ;
  minDF0  = OPTS.minDF0 ;
  f0      = OPTS.f0 ;
  dumpMin = 2    ;
  dumpMax = 100  ;

  l = 1 ;

  if c1 <= 0 || c1 >=c2 || c2 >= 1
    l = NaN ;
    fprintf( 1, 'WolfeLineSearch: Bad parameters c1 = %g c2=%g\n', c1, c2) ;
    return ;
  end

  if isnan(Df0) || Df0 >= 0
    l = NaN ;
    fprintf( 1, 'WolfeLineSearch: Bad search direction Df0 = %g\n', Df0 ) ;
    return ;
  end

  % if derivative is too strong use Simple Armijo
  if Df0 < minDF0
    l = ArmijoSimpleSearch( fun, xk, f0, Df0, d, varargin{:} ) ;
    return ;
  end
  
  c1Df0 = c1 * Df0 ;
  c2Df0 = c2 * Df0 ;
  while l > lMin
    fl = feval( fun, xk + l * d ) ;
    if fl <= f0 + l*c1Df0 ;
      % goto ZOOM
      break ;
    else
      if l == 1
        lTmp = quadratic( f0, Df0, fl, l ) ;
      else
        lTmp = cubic( f0, Df0, fl, l, fp, p ) ;
      end ;
      p  = l ;
      fp = fl ;
      l  = max( min( lTmp, l/dumpMin ), l/dumpMax ) ;
    end ;
  end ;

  if l < lMin
    fprintf( 'WolfeLineSearch: failed lambda too short %g\n', l ) ;
    l = NaN ; % Failed search
    return ;
  end ;

  %   ________   ___  __  __ 
  %  |__  / _ \ / _ \|  \/  |
  %    / / | | | | | | |\/| |
  %   / /| |_| | |_| | |  | |
  %  /____\___/ \___/|_|  |_|
  %
  [f,Df] = feval( fun, xk + l * d ) ;
  Dfl = Df*d ;
  if Dfl >= c2Df0 ;
    return ; % found Wolfe point
  end ;

  if l == 1
    % forward search
    while l <= lMax 
      p  = l ;
      fp = fl ;
      l  = 2*l ;
      fl = feval( fun, xk + l * d ) ;
      if fl > f0 + l*c1Df0
        p  = l ;
        fp = fl ;
        % goto REFINE
        break ;
      end ;
      [f,Df] = feval( fun, xk + l * d ) ;
      Dfl = Df * d ;
      if Dfl >= c2Df0 ;
        return ; % found Wolfe point
      end
    end
  end

  if l > lMax
    fprintf( 'WolfeLineSearch (zoom): failed lambda too short %g\n', l ) ;
    l = NaN ; % failed search
    return ;
  end

  %   ____  _____ _____ ___ _   _ _____ 
  %  |  _ \| ____|  ___|_ _| \ | | ____|
  %  | |_) |  _| | |_   | ||  \| |  _|  
  %  |  _ <| |___|  _|  | || |\  | |___ 
  %  |_| \_\_____|_|   |___|_| \_|_____|
  %
  lLo  = l ;
  fLo  = fl ;
  DfLo = Dfl ;

  lHi  = p ;
  fHi  = fp ;

  Delta = lHi - lLo ;
  
  while abs(Delta) > lepsi

    deltaLambda = quadratic( fLo, DfLo, fHi, Delta ) ;
    if deltaLambda > 0
      deltaLambda = min( Delta/dumpMin, max( Delta/dumpMax, deltaLambda )) ;
    else
      deltaLambda = max( Delta/dumpMin, min( Delta/dumpMax, deltaLambda )) ;
    end
    l  = lLo + deltaLambda ;
    fl = feval( fun, xk + l * d  ) ;

    if fl > f0 + l*c1Df0 || fl > fLo
      lHi = l ;
      fHi = fl ;
    else
      [f,Df] = feval( fun, xk + l * d ) ;
      Dfl = Df * d ;
      if Dfl >= c2Df0;
        return ; % found Wolfe point
      end ;

      lLo  = l ;
      fLo  = fl ;
      DfLo = Dfl ;
    end
    Delta = lHi - lLo ;
  end
  fprintf( 1, 'WolfeLineSearch (refine): failed refine lHi=%g lLo=%g DfLo=%g Df0=%g\n', lHi, lLo, DfLo, Df0 ) ;
  l = NaN ;
  return ; % failed search
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

