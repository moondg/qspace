function [kc,Ic]=getCurrentSite(HAM,varargin)
% function [kc,Ic]=getCurrentSite(HAM [,opts])
% Wb,Apr12,14 ; Wb,Mar06,21

  getopt('INIT',varargin);
     if     getopt('-V'), vflag=3;
     elseif getopt('-v'), vflag=2;
     elseif getopt('-q'), vflag=0; else vflag=1; end
  getopt('check_error');

  L=numel(HAM.mpo);

  [ff,Ix]=getFiles(HAM);

  sdir=Ix.sdir; kc=-1;

  [~,kk]=sort(abs((1:L)-(Ix.kc+0.25*sdir))); nx=0; l=0;
  for k=kk, l=l+1;
     q=load_dmrg_data(HAM,k,'info'); itags{k}=q.itags;
     if ~isempty(q.itags)
        i=find(cellfun(@numel, regexp(q.itags(1:3),'\*$')),1);
        if numel(i)>1, wberr('got invalid MPS setting'); 
        elseif isempty(i)
           kc=k; break
        elseif k==1 && i==1, kc=k; sdir=-1; break;
        elseif k==L && i==2, kc=k; sdir=+1; break;
        end
     else nx=nx+1; end
  end

  msg='';
  if nx, if kc>0, kc=-kc; else kc=-Ix.kc; end
     msg=sprintf('%g/%g sites unitialized',nx,L);
     if l<L, msg=[msg sprintf(' (%g checked)',l)]; end
     if nargout<2 && vflag
        wblog('WRN','MPS not yet initialized (%g/%g)',nx,L);
     end
  elseif kc<0, wberr('failed to identify current site'); end

  if nargout>1
     Ic=add2struct('-',kc,sdir,'active=0',msg);

     tt=sort(Ix.tt);
     Ic.time=datestr(tt(end));
     Ic.dt=24*3600*(tt(end)-tt(1))/L;

     Ic.active=0;
     dt=max(1/(24*3600),mean(diff(tt))); dt(2)=(now-tt(end))/dt;
     if dt(2)<L, Ic.active=1; end
     Ic.folder=Ix.folder; Ic.bytes=Ix.bytes;
  end

  if nargout || ~vflag, return; end

  I=load_dmrg_data(HAM,1,'info');
  if isfield(I,'isw')
     if isempty(I.isw), isw=1; else isw=I.isw+1; end
  else isw=0; end

  if     sdir>0, sdir='(>>)';
  elseif sdir<0, sdir='(<<)'; else sdir=''; end

  sx='';
  if vflag && kc==floor(kc) ...
     Ic=load_dmrg_data(HAM,kc,'info');
     q=findstrc(Ic.itags,'\*$','-n');
     if isequal(q,[0 0 0 1])
        if vflag>2
           q=load_dmrg_data(HAM,kc,'AK');
           d=getDimQS(q.AK); fprintf(1,[
             '\n   Global symmetry QPsi = [%s] (%d multiplets, ' ... 
             '%d states)\n'],mat2str(uniquerows(q.Q{4})),d(:,end));
        else
           sx='got global Psi index';
        end
     elseif ~isequal(q,[0 0 0]), disp(q)
        wblog('WRN','invalid current site !?'); 
     end
  end

  if kc<0
     if ~isempty(sx), sx=[', ' sx]; end
     fprintf(1,'   Current site not yet specified (kc=%g%s)\n',kc,sx);

  elseif round(kc)==kc
     if ~isempty(sx), sx=[' (' sx ')']; end
     fprintf(1,'   Current sweep.site is kc=%g.%g %s%s\n',isw,kc,sdir,sx);

  elseif round(2*kc)/2==kc,
     if isempty(sdir), sdir='(unknown sweepdir)'; end
     ss={ '\e[35m', sprintf('kc=%g:%g %s',isw,kc,sdir),'\e[0m','' };
     if ~isempty(sx), ss{end}=[', ' sx]; end
     printfc(['   ' ss{1} 'Current sweep.bond is ' ss{2:end} '\n']);
  else
     ss={ '\e[35mWRN ', sprintf('(kc=%g',kc),'',' !?)\e[0m\n' };
     if ~isempty(sx), ss{end-1}=[', ' sx]; end
     printfc(['   ' ss{1} 'failed to determine current bond ' ss{2:end} ]);
  end

  clear kc

end

