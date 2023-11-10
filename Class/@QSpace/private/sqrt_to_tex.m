function str=sqrt_to_tex(str,xflag)
% function str=sqrt_to_tex(str [,xflag])
% Wb,Oct19,23

% tags: srd_to_tex, surd_to_tex

  srd=char(8730);
  str=regexprep(str,'sqrt\(([^()]*)\)',[srd '$1']);

  srds='\surd';
  if nargin>1 && xflag, srds=['\' srds]; end
  srds=['?$' srds '$?'];

  for i=fliplr(find(str==srd))
      str=[str(1:i-1) srds str(i+1:end)];
  end

end

