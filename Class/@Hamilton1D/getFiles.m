function [ff,Iout]=getFiles(HAM)
% function [ff,Iout]=getFiles(HAM)
%
%   Get data file listing (not including *_info.mat)
%
%   If second output argument is requested, this uses the file
%   time stamps as estimate (first guess) to find kc.
%   Note also that in case the DMRG is stored in work space,
%   no files (and thus also no file stamps) are available.
%
% Adapted from @Hamilton1D/private/setupStorage()
% Wb,Sep11,15

% WRN! rather unreliable/fragile in test sessions
% where files may have been resaved! // Wb,May31,19

  L=numel(HAM.mpo);

  if ~isempty(HAM.store), s=HAM.store;
     if nargout>1, ff={}; Iout=struct;
     else ff=dir(''); end
     return
  end

  if isempty(HAM.mat), wberr('invalid storage specification'); end
  mat=HAM.mat;

  ff=dir(sprintf('%s_*.mat',mat));
  if isempty(ff)
     wblog('WRN','no files found');
     if nargout>1, ff={}; Iout=struct; end
     return
  end

  i=find([ff.isdir]);
  if ~isempty(i), disp(ff(i(1)))
     wberr('file listing includes directories !? (L=%g)',L); 
  end

  for i=1:numel(ff)
     k=0; regexp(ff(i).name,'.*_([0-9]+)(?@k=str2num($1);)\.mat');
     ff(i).k=k;
  end

  i=find([ff.k]>0);
  [kk,is]=sort([ff(i).k]); ff=ff(i(is));

  if ~isequal(kk,1:L), q=[min(kk), max(kk), numel(kk), L];
     wberr('unexpected DMRG file listing (%g .. %g; %g/%g)',q); 
  end

  if nargout<2, return; end

  tt=[ff.datenum];
  kc=find(tt==max(tt));

  k1=kc(1); k2=kc(end);

  if     k1<= 3, sdir=+1;
  elseif k2>L-3, sdir=-1; else
     q=[ mean(tt(1:k1-1)), mean(tt(k2+1:end)) ];
     sdir=-sign(diff(q));
  end

  if k1<k2
     if k1<=3, kc=k2;
     elseif k2>L-3, kc=k1; else
        q=[ mean(tt(1:k1-1)), mean(tt(k2+1:end)) ];
        if diff(q)>0
             kc=k1;
        else kc=k2; end
     end
  end

  Iout=struct('folder',ff(1).folder,...
    'bytes',sum([ff.bytes]),'tt',[ff.datenum],'kc',kc,'sdir',sdir);
  ff={ff.name};

end 

