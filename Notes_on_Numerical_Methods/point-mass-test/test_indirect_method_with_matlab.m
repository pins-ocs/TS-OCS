
function [z,elapsed,ok] = test_indirect_method_with_matlab(auxdata_in)
  global NF_eval NJF_eval auxdata ;
  auxdata = auxdata_in ;
 
  NF_eval  = 0 ;
  NJF_eval = 0 ;

  % compute guess solution
  z0 = indirect_method_guess_solution( auxdata ) ;

  opt  = optimoptions('lsqnonlin','Display','iter', ...
                      'Jacobian','on', 'DerivativeCheck', 'off', ...
                      'TolX', 1e-8, 'TolFun', 1e-8, ...
                      'MaxFunEvals',1000000,'MaxIter',5000) ;
  tic
  auxdata.epsilon = 1e-1 ;
  [z1,resnorm,residual,exitflag,output] = lsqnonlin(@f_model,z0,[],[],opt) ;
  auxdata.epsilon = 1e-3 ;
  [z2,resnorm,residual,exitflag,output] = lsqnonlin(@f_model,z1,[],[],opt) ;
  auxdata.epsilon = auxdata_in.epsilon ;
  [z,resnorm,residual,exitflag,output] = lsqnonlin(@f_model,z2,[],[],opt) ;
  elapsed = toc ;
  
  if exitflag > 0 && exitflag <= 3
    ok = output.iterations ;
  else
    ok = -1 ;
    stop ;
  end

end

function varargout = f_model( x )
  global auxdata ;
  varargout{1} = indirect_method_F( x, auxdata ) ;
  if nargout > 1
    varargout{2} = indirect_method_JF( x, auxdata ) ;
  end
end
