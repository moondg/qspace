function get_data_from_info(HAM)
% function get_data_from_info(HAM)
%
%    In case matlab session was not saved in tst_Hamilton1D,
%    (re)extract data from info fields
%
% Wb,Sep12,15

   L=numel(HAM.mpo);

   q=struct('E0',[],'e0',[],'se',[],'rr',[],'r2',[],'NK',[],'Dk',[]);

   I=load_dmrg_data(HAM,2,'info');
   ns=ceil(size(I.dsw,1)/2);

   IL=repmat(q,[L-1,2,ns]);
   Il=repmat(q,[L-1,2]);

   for k=1:L-1
      fprintf(1,'  k=%2g/%g ... \r',k,L);

      I=load_dmrg_data(HAM,k,'info');

      ns=size(I.dsw,1);

      for i=1:ns
         if mod(i,2), j=1; else j=2; end

         if k>1 && k<L
            is=floor((i+1)/2);
            if mod(i,2), j=1; else j=2; end
         else
            is=i;
            if k==1, j=2; else j=1; end
         end

         q.NK=I.dsw(i,1);
         q.se=I.dsw(i,2);
         q.E0=I.dsw(i,3);
         q.r2=I.dsw(i,5);

         if iscell(I.ek)
              x=I.ek{is,j};
         else x=I.ek; end

         if numel(x)>1, x=sum(x(:,end));
         elseif isempty(x), x=nan; end

         q.e0=x;

         if i==1
            Dk=getDimQS(I.Psi);
            q.rr=I.rr;
         end

         q.Dk=Dk;
         if q.Dk(1)~=q.NK
            q.Dk(1)=q.NK; q.Dk(end)=nan;
         end

         IL(k,j,is)=q;
         if i+2>ns, Il(k,j)=q; end
      end
   end
   fprintf(1,'\r%60s\r','');

   assignin('caller','IL',IL);
   assignin('caller','Il',Il);

   NK=getdatafield(IL,'NK','-q'); s=size(NK);
   n1=max(reshape(NK,s(1)*s(2),[]),[],1);

   assignin('caller','NKEEP',n1);
   assignin('caller','Nkeep',max(n1));
   assignin('caller','isw',numel(n1));
   assignin('caller','nsw',numel(n1));

   rtol=getopt(HAM.info,'sweep','rtol',1E-15);
   assignin('caller','rtol',rtol);

end

