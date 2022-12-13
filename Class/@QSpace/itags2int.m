function kk=itags2int(A)
% function k=itags2int(A)
%
%    Extract integer indices carried with itags,
%    assuming the itag format '[non-number]##[non-number]'
%    where outgoing indices (i.e., which carry a trailing
%    conjugate flag `*') are returned as negative integers.
%    If no valid number is encountered, nan is returned
%    for that index.
%
% Wb,Nov19,21

% adapted from itag2iter()
% see also Hamilton1d/privte/get_kidx.m which returns the same data
% for a single QSpace A. // Wb,Nov19,21

  kk=cell(size(A)); n=numel(A);

  for i=1:n
     if ~isfield(A(i).info,'itags'), kk{i}=[]; continue; end
     it=A(i).info.itags; r=numel(it); e=0;
     x=regexprep(it,'^[^\d]*(\d+)[^\d]*$','$1');
     for j=1:r
        x{j}=str2num(x{j});
        if numel(x{j})==1, if it{j}(end)=='*', x{j}=-x{j}; end
        else x{j}=nan; end
     end
     kk{i}=[x{:}];
  end
  if n==1, kk=kk{1}; end

end

