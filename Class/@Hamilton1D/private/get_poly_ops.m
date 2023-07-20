function [Sab,Ix]=get_poly_ops(S,J,varargin)
% function Y=get_poly_ops(S,n [,sstr])
% function Y=get_poly_ops(S,J [,sstr,opts])
%
%    get interaction bilinears for powers of S.S interactions
%    acting across two sites 1 and 2.
%
%    Usage #1 (with the second argument a single integer):
%    this returns paired-up operators
%    having Y(2k-3).Y(2k-2)' = (S.S')^k for k=2:n.
%
%    Usage #2 (with a vector J):
%    this returns the single pari of operators
%    having Y(1).Y(2)' = \sum_{k=1:n} J(k) * (S.S')^k 
%
% Options
%
%    sstr   string for S (for logging, and for HAM.ops returned)
%    --tr   subtract trace, i.e., make traceless (usage #2 only)
%
% Wb,Apr02,20

% adpated from former get_aklt_ops.m

% -------------------------------------------------------------- %
% NB! (S.S)^kterms (as well as Heisenberg interactions between   %
% different sites with different spin S require the contraction  %
% Sa.Sb' with Sa!=Sb (!)  // Wb,Aug27,15                         %
% -------------------------------------------------------------- %

  if nargin<2
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  getopt('init',varargin);
     notr=getopt('--notr');
  sstr=getopt('get_last','S');

  if numel(sstr)==1
       s1=sstr; s2='';
  else s1='S'; s2=[' ' sstr]; end

  if notr, o={'--notr'}; else o={}; end

  S2=contract(S,3,S,'3*','1423'); % (a a',b'b) -> (a b,a'b')

  if numel(J)==1 && J>=1 && J==round(J), n=J;
     if n>9, wbdie('check value for n=%g !?',n); end
     SK=S2; s1=sprintf('(%s,%s)',s1,s1);
     for k=2:n
        SK=contract(SK,'34',S2,'12');
        s=sprintf('%s^%g%s',s1,k,s2);
        [Sa,Sb,Ix]=splitH4_SdotS(SK,'-v','istr',s,o{:});

        Sab(1,k-1)=init_ops(Sa,[s '_a'],'~hconj','-q');
        Sab(2,k-1)=init_ops(Sb,[s '_b'],'~hconj','-q');
     end
     Sab=Sab(:);
  else
     SK=S2; X2=J(1)*SK; n=numel(J);
     for k=2:n
        SK=contract(SK,'34',S2,'12'); X2=X2+J(k)*SK;
     end

     s=sprintf(',%.3g',J); s=sprintf('J%s=[%s]%s',s1,s(2:end),s2);
     [Sa,Sb,Ix]=splitH4_SdotS(X2,'-v','istr',s,o{:});

     Sab(1)=init_ops(Sa,[s '_a'],'~hconj','-q');
     Sab(2)=init_ops(Sb,[s '_b'],'~hconj','-q');
  end

end

% -------------------------------------------------------------------- %

