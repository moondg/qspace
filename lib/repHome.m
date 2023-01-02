function fname=repHome(fname)
% Function: fname=repHome(fname)
% Wb,2006 ; Wb,Nov05,10

% this is only used to shorten long paths to personalized tags
% feel free to adapt this to your own needs, otheriwse ignore
% Wb,Aug04,17

% [fname,i]=system_getline(['pwdprompt ''' fname '''']);

  U=getenv('USER');

% e.g. see also matlab command 'cto lma'
  qq = { 
     ['/data/.*'    U '/Matlab(|/Data)'], '$LMA'
     ['/scratch/.*' U '/Matlab/Data'   ], '$SCR'
     ['/data/.*'    U '/bin'           ], '$LBIN'
     ['/scratch/.*' U '/bin'           ], '$SBIN'
  };

  qx={'RC_STORE', '$RCS'
      'RC_SYNC',  '$RSY' 
      'HOME',     '~' 
  };
  for i=1:size(qx,1)
     q=getenv(qx{i});
     if ~isempty(q), qq(end+1,:)={q,qx{i,2}}; end
  end

  for i=1:size(qq,1)
     fname=regexprep(fname,qq{i,:});
  end

end

