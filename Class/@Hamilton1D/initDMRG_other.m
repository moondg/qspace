function Iout=initDMRG_other(HAM,dmrg0,varargin)
% function Iout=initDMRG_other(HAM,dmrg0, [,opts])
%
%    initialize A based on other DMRG run.
%
% Options
%
%   --RCStore  also rsync RCStore from earlier run
%
% Wb,Sep02,20

  getopt('INIT',varargin);
     if getopt('-v'), vflag=2;
     elseif getopt('-q'), vflag=0; else vflag=1; end
     RCflag=getopt('--RCStore');
  getopt('check_error');

  if ~ischar(dmrg0) || length(dmrg0)<6, dmrg0, 
     error('Wb:ERR','invalid usage (unexpected DMRG reference)'); end

  [ff,I2]=get_files(HAM.mat); nf=numel(ff); L=nf-1;

  if ~isempty(regexp(dmrg0,'[a-zA-Z\/]'))
     [f0,I0]=get_files(dmrg0); L0=numel(f0)-1;

     if L~=L0, error('Wb:ERR',...
       'invalid usage (system length mismatch L=%g/%g',L0,L); end
     check_same_system(HAM,I0.HAM);

     if isequal(ff(1).folder,f0(1).folder) && isequal(ff(1).name,f0(1).name)
        wblog('WRN','got same DMRG data !? (return)'); 
        return
     end

     if vflag
        fprintf(1,['\n>> initializing with existing DMRG data\n' ... 
        '        %s\n   -> %s\n\n'],repHome(dmrg0),repHome(HAM.mat));
     end

     if RCflag
        q=[f0(1).folder '/../RCStore']; 
        if exist(q,'dir')
           fprintf(1,'>> rsyncing RCSore %s -> %s\n\n',...
              repHome(q),repHome(ff(1).folder));
           system(['rsync -var ' q ' ' ff(1).folder '/../ | wc']);
        else
           fprintf(1,['>> rsyncing RCSore does not exist - ' ... 
           'ignore\n\n'], repHome(q));
        end
     end

     for k=1:nf
        mat0=[f0(k).folder '/' f0(k).name];
        mat2=[ff(k).folder '/' ff(k).name];

        fprintf(1,'%5g/%g %s \r',k,L,mat2);

        X0=load(mat0);
        X2=load(mat2);

        if k<=L
           X2.info.itags=X0.info.itags;
           X0.info=X2.info;
        else
           X2.NPsi=X0.NPsi;
           X0=X2;
        end

        save(mat2,'-struct','X0');
     end
     fprintf(1,'\n');
  else
     error('Wb:ERR','[to be cont''d'); 
  end

  kc=getCurrentSite(HAM);
  if kc>1
     wblog('WRN','got current site kc=%g => 1',kc); 
     setCurrentSite(HAM,1);
  end

end

% -------------------------------------------------------------------- %

function check_same_system(HAM,HM2)

   q=rmfield([HAM.info.IS, HM2.info.IS],'istr');
   if ~isequal(q(1),q(2))
       wbdie('got system mismatch (different HAM.IS)'); end

   if ~isequal(HAM.oez,HM2.oez)
       wbdie('got system mismatch (different HAM.oez)'); end

    q={ struct(HAM.ops), struct(HM2.ops) };
    for i=1:numel(q), for j=1:numel(q{i})
        q{i}(j).op=q{i}(j).op.Q;
    end, end

   if ~isequal(q{1},q{2})
       wbdie('got system mismatch (different HAM.ops)'); end

   q=[length(HAM.mpo), length(HM2.mpo)];
   if diff(q) 
       wbdie('system length mismatch (L=%g/%g)',q);
   end
end

% -------------------------------------------------------------------- %

function [ff,I]=get_files(mat)

  i=[mat '_*.mat']; ff=dir(i); ff(find([ff.isdir]))=[];
  p=regexprep(mat,'.*/','');
  i=find(cellfun(@numel, regexp({ff.name},[p '_' '(info|\d+)\.mat'])));
  ff=ff(i); 

  if isempty(ff)
   % D=regexprep(mat,'/[^\/]*$','');
   % if exist(D,'dir')
   %      I=dir(D); I=I(1);
   %      I.name=p;
   % else error('Wb:ERR','invalid directory ''%s''',p); end
   % return
     error('Wb:ERR','invalid usage (missing %s files)',mat);
  end

  if isempty(regexp(ff(end).name,'_info'))
     error('Wb:ERR','invalid HAM.mat'); end

  fI=ff(end); fI=[fI.folder '/' fI.name];
  I=load(fI); I.info_mat=fI;

  L=[numel(ff)-1, numel(I.HAM.mpo)];
  if diff(L), error('Wb:ERR',sprintf('invalid HAM (L=%g/%g)',L)); end

end

% -------------------------------------------------------------------- %

