function qfmt=get_q_fmt(qtype,r,m)
% introduce compact notation for SU(N), Sp(2n), etc
% Wb,Apr02,15

% see also QSpace::getqfmt() // Wb,Dec14,15

   q=strread(qtype,'%s','delimiter',',;')';

   if isempty(q)
      qfmt=repmat(' %2g',1,m);
      if r>1, qfmt=repmat([ ' ;' qfmt],1,r); qfmt=qfmt(4:end);
      else qfmt=qfmt(2:end); end
      return
   end

   for i=1:numel(q)
      j=regexp(q{1,i},'\d');
      if isempty(j), q{2,i}='%2g';
      else
         w=q{1,i}(1:j-1); n=str2num(q{1,i}(j:end));
         switch w
            case 'SU', q{2,i}=repmat('%X',1,n-1);
            case 'Sp', q{2,i}=repmat('%X',1,n/2);
            case 'SO', q{2,i}=repmat('%X',1,floor(n/2));
            case {'Z','P'}, q{2,i}='%2g';
            otherwise wbdie('unexpected symmetry %s',q{1,i});
         end
      end
   end
   q(3,:)={' '}; q{3,end}='';

   qfmt=[q{2:3,:}];

   if nargin>1
      qfmt=repmat({qfmt},2,r);
      qfmt(2,:)={' ; '}; qfmt{end}='';
      qfmt=[qfmt{:}];
   end

end

