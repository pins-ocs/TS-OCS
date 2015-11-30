function [z,elapsed,ok] = test_indirect_method_with_strscne(auxdata_in)

  global auxdata ;
  addpath('../scripts') ;
  
  auxdata = auxdata_in ;

  % compute guess solution
  z0 = indirect_method_guess_solution( auxdata ) ;
  
  % compute guess solution
  lb = -1000*ones(auxdata.nvars,1) ;
  ub = 1000*ones(auxdata.nvars,1) ;
  
  tol   = [1e-8 1e-8] ;
  parms = [10000,10000,-2,0] ;

  tic
  auxdata.epsilon = 1e-2 ;
  [z1,ierr,output,history,grad]=STRSCNE(z0,@fun,tol,lb,ub,parms,@Jfun) ;
  ok = output(1) ;
  if ierr == 0
    auxdata.epsilon = 1e-4 ;
    [z2,ierr,output,history,grad]=STRSCNE(z1,@fun,tol,lb,ub,parms,@Jfun) ;
    ok = ok + output(1) ;
    auxdata.epsilon = auxdata_in.epsilon ;
    if ierr == 0
      [z,ierr,output,history,grad]=STRSCNE(z2,@fun,tol,lb,ub,parms,@Jfun) ;
      ok = ok + output(1) ;
    end
  end
  elapsed = toc ;

  if ierr ~= 0
    ok = -1 ;
  end

end

function RES = fun( x )
  global auxdata ;
  RES = indirect_method_F( x, auxdata ) ;
end

function RES = Jfun( x )
  global auxdata ;
  RES = indirect_method_JF( x, auxdata ) ;
end
