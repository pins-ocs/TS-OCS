function [z,elapsed,ok] = test_direct_method_with_matlab(auxdata_in)

  global NF_eval NJF_eval auxdata ;

  auxdata = auxdata_in ;
 
  NF_eval  = 0 ;
  NJF_eval = 0 ;

  [lb,ub,cl,cu] = direct_method_bound(auxdata) ; 

  % compute guess solution
  z0 = direct_method_guess_solution( auxdata ) ;

  opt = optimoptions('fmincon','Display','iter', ...
                     'GradObj','on', 'GradConstr','on',...
                     'TolX', 1e-8, 'TolFun', 1e-8, ...
                     'MaxFunEvals',1000000,'MaxIter',5000) ;
  tic
  [z,fval,exitflag,output] = fmincon(@f_model,z0,[],[],[],[],lb,ub,@f_constr,opt) ;
  elapsed = toc ;
  
  if exitflag == 1
    ok = output.iterations ;
  else
    ok = -1 ;
  end

end

%% 
% map the indices with the corresponding index in the spase matrix
function varargout = f_model(z)
  global auxdata ;
  varargout{1} = direct_method_objective(z,auxdata);
  if nargout > 1
    varargout{2} = direct_method_gradient(z,auxdata);
  end
end

% map the indices with the corresponding index in the spase matrix
function varargout = f_constr(z)
  global auxdata ;
  tmp = direct_method_constraints(z,auxdata);
  idx = 2*auxdata.N+3 ;
  varargout{1} = -tmp(idx+1:end) ;
  varargout{2} = tmp(1:idx) ;
  if nargout > 2
    tmp = direct_method_constraints_jacobian(z,auxdata);
    varargout{3} = -tmp(idx+1:end,:).' ;
    varargout{4} = tmp(1:idx,:).' ;
  end
end
