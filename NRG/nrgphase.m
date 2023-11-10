function [ph,I]=nrgphase(NRG,varargin)
% function [ph,I]=nrgphase(NRG [,opts])
%
%    determine phase shifts from last NRG iteration.
%    by specifying either the NRG-tag or the Hamiltonian
%    in the K(ept) basis of given NRG-iteration.
%
%    NB! it is assumed that first quantum number of the
%    HK QSpace refers to charge count (use option 'iqn', otherwise).
%
% Options
%
%   'iqn',..  which quantum number(s) to consider particle-numbers (1)
%             NB! this decides on particle/hole character of phase shifts
%   'eps',..  energy threshold relative to dE0 for degeneracy
%   '-v'      verbose
%   '-q'      quiet
%   '-p'      plot results
%
% Examples
%
%    [ph,I]=nrgphase([k,]'NRG/NRG'); % without k, pick even/odd with lowest E0
%    [ph,I]=nrgphase(HK);
%
% Wb,Mar08,09

  if nargin<1
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  if ~ischar(NRG)
     if ~isa(NRG,'QSpace') && ~isfield(NRG,'Q')
        wbdie('invalid usage (invalid NRG space)');
     end
     HK=NRG;
  else
     if nargin>1 && isnumeric(varargin{1})
        k=varargin{1}; varargin=varargin(2:end);
     else
        f=[NRG '_info.mat'];
        if ~exist(f,'file')
           wbdie('invalid NRG space (file %s does not exist)',f);
        end

        load(f,'E0'); n=length(E0);
        if E0(end)<E0(end-1)
        k=n-1; else k=n-2; end
     end

     f=sprintf('%s_%02d.mat',NRG,k);
     if ~exist(f,'file')
        wbdie('invalid NRG space (file %s does not exist)',f);
     end

     s=load(f,'HK','HT');
     if isempty(s.HK.Q), HK=s.HT; else HK=s.HK; end

  end

  getopt('init',varargin);
     iqn  =getopt('iqn',1);
     eps  =getopt('eps',0.02);
     pflag=getopt('-p');
     vflag=getopt('-v');
     qflag=getopt('-q');
  getopt('check_error');

  [EQD,qq]=getEQdata(QSpace(HK),'iq',iqn);
  if isempty(EQD), wbdie('empty QSpace !?'); end

  if EQD(1)~=0
     wblog('WRN','lowest energy is unequal zero !?? (%g)',EQD(1));
     EQD(:,1)=EQD(:,1)-EQD(1); n=numel(HK.data);
     for i=1:n, HK.data{i}=HK.data{i}-EQD(1); end
  end

  D0=EQD(:,end);
  ee=EQD(:,1);
  i0=find(ee<ee(1)+eps); II=[]; Er=[];

  while 1
     [ph,e0,EQD,Er]=nrgphase_1(EQD,Er,i0,eps,qflag);
     if isempty(ph), break; end

     q=add2struct('-',ph,e0);
     if ~isempty(II), II(end+1)=q; else II=q; end
  end

  ph=cat(1,II.ph);
  e0=cat(1,II.e0); eps(2)=min(e0)/10;

  D1=EQD(:,end); EQD(:,end)=D0;

  m=EQD(i0,end);
  if numel(i0)==1 && m>1
    if mod(length(ph),m), wbdie(...
      'expecting multiple of ground state degeneracy (%g/%g=?)',numel(ph),m);
    end
    e=norm(diff(reshape(ph,m,[]),1));
    if e<1E-4
       wblog(' * ','removing ground state degeneracy (g=%g)',m); 
       ph_=ph; ph=ph(1:m:end);
       e0_=e0; e0=e0(1:m:end);
    else
       wblog('ERR','failed to deal with ground state degeneracy (g=%g)',m); 
    end
  end

  if ~qflag
     wblog(' * ','got %g phase shifts',numel(ph));
     fprintf(1,'\n   ##  phase/pi  energy e0\n  %s\n',repmat('-',1,30));
     fprintf(1,'% 4d. %8.5g %8.4g\n',[(1:numel(ph))',ph,e0]');
     fprintf(1,'\n');
  end

  if numel(ph)>4,
     wblog('WRN','got range of phases (match_fermisea may be slow)');
     keyboard
  end

  EQD(:,end)=1;
  [EQD,m,Im]=match_fermisea(EQD,ph,e0,eps(2),qflag);

  dn=D0-EQD(:,end);  EQD=[ EQD(:,1:end-1), D0, D1, EQD(:,end) ];
  I=add2struct('-',ph,e0,Im,'ph_?','e0_?',EQD,eps);
  if ~pflag || isempty(ph), return; end

ah=smaxis(2,2,'tag',mfilename); header('%M'); addt2fig Wb

  plotQSpectra(QSpace(HK),'ah',ah(1,:),'-x')

setax(ah(2,1)); xl=[0 min(120,length(ee))];

  plot(ee,'o-','Disp','matched'); xlim(xl); sms(4); hold on; grid on
  ylabel('energy')

  i=find(EQD(:,end));
  if ~isempty(i)
     h=plot(i,EQD(i,1),'ro','Disp','non-matched'); sms(h,4);
     legdisp('Location','SouthEast','erase');
  end

  for i=1:numel(ph)
    postext(0.04,0.99-0.06*i,'\phi_%g=%6.5g,  E_%g=%6.5g',i,ph(i),i,e0(i),...
    {'BackgroundC','w','Margin',0.01,'FontSize',10});
  end

setax(ah(1,iqn))

  i=find(dn);
  if ~isempty(i)
     h=plot(EQD(i,2),EQD(i,1),'o','Color',[1 .5 0],'Disp','matched');
     sms(h,4);

     h=plot(Im.qx+EQD(i0(1),2),Im.ex,'bx','MarkerS',6,...
      'Disp','fermisea reference data');
     legdisp('Orientation','Horizontal','Location','North','dx',[-0.03 0.07]);

     Er=EQD(find(EQD(:,end),1)); ymark(Er);
     wblog(' * ','matched all multiplets up to E=%.4g',Er);
  end

setax(ah(2,2))

  h=plot(qq,'.-'); ytight(1.1); xlim(xl); hold on; grid on
  i=1:numel(h); i(iqn)=[]; 
  set(h(i),'Marker','o'); sms(h(i),4); % 'LineSt','none',
  mv2back(blurl(h(i))); set(h(i),'LineW',1);
  label('state index','Q-labels');

end

% -------------------------------------------------------------------- %
function [ph,e0,EQD,Er]=nrgphase_1(EQD,Er,i0,eps,qflag)

   i0=i0(find(EQD(i0,end),1)); ph=[]; e0=[];
   if isempty(i0), return; end

 % NB! charge excitation may occur across channels!
 % strategy:
 % * take particle-excitation (dq=+1)
 % * match with all hole-excitations (dq=-1) with energy<e0_approx
 %   to find a match for e0_exact
   jp=find(EQD(:,2)==(EQD(i0,2)+1) & EQD(:,3),1);

   if isempty(Er)
      jn=find(EQD(:,2)==(EQD(i0,2)-1) & EQD(:,3),1);
      Er=sum(EQD([jp,jn]));
   end

   jn=find(EQD(:,2)==(EQD(i0,2)-1) & EQD(:,3) & EQD(:,1)<(Er+eps));
   j0=find(EQD(:,2)== EQD(i0,2)    & EQD(:,3)); E0=EQD(j0,1);

   if isempty(jn), return; end

   for i=1:numel(jn)
      e0=sum(EQD([jp,jn(i)]));
      [dE,is]=sort(abs(E0-e0));
      if i>1 && dE(1)<q(3) || i==1
         q=[i, is(1), dE(1)];
      end
   end
   jn=jn(q(1)); j0=j0(q(2)); dE=q(3); e0=sum(EQD([jp,jn]));

   if ~qflag, wblog('==>',...
      '%8.5g + %8.5g = %8.5g (@%.3g)',EQD([jp,jn]),e0,abs(EQD(j0)-e0)); 
   end

   if dE>=eps
      return
   end

   ph=EQD(jp)/e0;

   i=[jn,j0,jp]; EQD(i,end)=EQD(i,end)-1;
end

% -------------------------------------------------------------------- %

function [EQD,m,I]=match_fermisea(EQD,ph,e0,eps,qflag)

       if numel(ph)<=4, Er=5*min(e0);
   elseif numel(ph)<=6, Er=4*min(e0);
   elseif numel(ph)<=8, Er=3*min(e0); else Er=2*min(e0); end

   [ex,qx]=fermisea(ph,e0,'Er',Er);

   nx=0; m=0; e2=Inf; err2=0; D0=EQD(:,end);
   for i=1:numel(ex)
      ix=find(EQD(:,2)==qx(i) & EQD(:,3)); if isempty(ix), continue; end
      [dE,is]=sort(abs(EQD(ix,1)-ex(i)));
      if dE(1)<eps
         ix=ix(is(1)); EQD(ix,3)= EQD(ix,3)-1; m=m+1;
         err2=err2+(EQD(ix,1)-ex(i))^2;
      else 
         e2=min(e2,ex(i));
         if ex(i)<2*e0, nx=nx+1; end
      end
   end

   if nargout>2, I=add2struct('-',ex,qx,nx,e2,err2); end

   if ~qflag
   wblog(' * ','found %g (%d/%d) matches (@ %.4g up to %.4g)',...
     m,length(ex),numel(find(D0-EQD(:,end))), ...
     sqrt(err2)/(m*min(e0)),e2/min(e0));
   end

end

% -------------------------------------------------------------------- %

