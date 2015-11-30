%   Copyright (C) 2011 Enrico Bertolazzi
% 
%   This program is free software; you can redistribute it and/or
%   modify it under the terms of the GNU General Public License
%   as published by the Free Software Foundation; either version 2
%   of the License, or (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details. 
% 
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
%   USA.
%
function OPTIONS = parseArgs( opts, defs, varargin )
  % OPTIONS = parseArgs(OPTS,DEFS,'PARAM1',VALUE1,...) creates a copy of OLDOPTS
  % with the named parameters altered with the specified values.
  % initialize to default values

  if length(opts) ~= length(defs)
    error('parseArgs','number of options differents of number of defaulty values') ;
  end
  
  for k=1:length(defs)
    OPTIONS.(opts{k}) = defs{k} ;
  end ;

  for k=1:2:nargin-2
    opt = varargin{k} ;
    if ~ischar(opt)
      error('parseArgs:NoPropName','Expected argument %d to be a string property name.', k);
    else
      if any(strcmp(opt,opts))
        value = varargin{k+1} ;
        OPTIONS.(opt) = value ;
      end
    end
  end
end
