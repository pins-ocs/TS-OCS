function [z,elapsed,ok] = test_indirect_method_with_strscne(auxdata_in)

  global auxdata ;
  addpath('../scripts') ;
  
  auxdata = auxdata_in ;
 
  % compute guess solution
  z0 = indirect_method_guess_solution( auxdata ) ;
  
  % compute guess solution
  lb = -1000*ones(auxdata.nvars,1) ;
  ub = 1000*ones(auxdata.nvars,1) ;
  
  options = TRESNEI();

  options.max_itns  = 100  ;
  options.max_feval = 20000 ;
  options.delta     = 1     ;
  options.jacobian  = 'on'  ;
  options.tol_F     = 1e-8  ;
  options.tol_opt   = 1e-8  ;
  options.output    = 0     ;

  tic
  auxdata.epsilon = 1e-1  ;
  [z1,ierr,output] = TRESNEI(z0,[auxdata.nvars,0],'indirect_method_f_model',lb,ub,options) ; 
  if ierr == 0
    ok = output.itns ;
    auxdata.epsilon = 1e-3  ;
    [z2,ierr,output] = TRESNEI(z1,[auxdata.nvars,0],'indirect_method_f_model',lb,ub,options) ;  
    if ierr == 0
      ok = ok + output.itns ;
      auxdata.epsilon = auxdata_in.epsilon  ;
      [z,ierr,output] = TRESNEI(z2,[auxdata.nvars,0],'indirect_method_f_model',lb,ub,options) ;  
      ok = ok + output.itns ;
    end
  end
  elapsed = toc ;

  if ierr == 0
    ok = output.itns ;
  else
    ok = -1 ;
  end

end
