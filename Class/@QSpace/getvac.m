function A=getvac(A,varargin)
% function A=getvac(A [,k])
%
%    reduce QSpace A to rank-2 Id tensor which only has
%    vacuum state / scalar representation with each index.
%    if k is specified, itag{k} will be used.
%
% Options
%
%    -1d   reduce to rank-1 vector (NB! only permitted for
%          for scalar representation with all q-lables zero!)
%
% Wb,Oct18,14

% outsourced to QSpace/getId.m

  getopt('init',varargin);
     d1flag=getopt('-1d');
  args=getopt('get_remaining');

  A=getId(A,[],args{:});

  if d1flag
     A=contract(getIdentityQS(A,1,A,1),'23',A,'12*');
  end

end

