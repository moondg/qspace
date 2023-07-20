function A=getIdentity(varargin)
% function A=getIdentity([args])
%
%    Wrapper routine to getIdentityQS
%    that returns QSpace object rather than plain structure.
%    
%    This wrapper supports a tweak in the trailing itag 
%    for the new fused state space: if represent, it can 
%    be specified as cell string { 'tag', 'x' } with some
%    extensions string 'x'; in this case trailing marks
%    such as ['*] with 'tag' are moved to the end of 'x'
%    instead, prior to concatening them.
%
% Wb,May28,16

% adapted from QSpace//contract.m

% adapted usage that permits to move trailing marks in t{1} ['*]
% to end of t{2}; an even number of marks cancel each other // Wb,Jul03,23
  if iscell(varargin{end}), t=varargin{end};
     if ~isempty(t) && ischar(t{1})
        n=0; t=regexprep(t,'([\*]+)$(?@n=n+length($1);)','');
        m=0; t=regexprep(t,'(['']+)$(?@m=m+length($1);)',''); t=[t{:}];
        if mod(m,2), t=[t '''']; end
        if mod(n,2), t=[t '*' ]; end
     end
  end

  A=class(getIdentityQS(varargin{:}),'QSpace');

end

