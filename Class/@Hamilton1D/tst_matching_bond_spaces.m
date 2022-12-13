function check_matching_bond_spaces(HAM)
% function check_matching_bond_spaces(HAM)
%
%    can it happen that through truncation certain symmetry
%    sector is present in A(k) but not A(k+1)?
%
% Wb,Apr24,14

  L=numel(HAM.mpo);
  Xl=load_dmrg_data(HAM,1);

  for k=2:L
     Xk=load_dmrg_data(HAM,k);
     [i1,i2,I]=matchIndex(Xl.AK.Q{2},Xk.AK.Q{1});
     if ~isempty(I.ix1) || ~isempty(I.ix2), w=[0 0];
        if ~isempty(I.ix1), x1=getsub(Xl.AK,I.ix1), w(1)=norm(x1)^2; end
        if ~isempty(I.ix2), x2=getsub(Xk.AK,I.ix2), w(2)=norm(x2)^2; end
        wblog(' * ','k=%g => {%s}(%.8g), {%s}(%.8g)\N',k-0.5,...
          vec2str(I.ix1),w(1),vec2str(I.ix2),w(2));
     end
     Xl=Xk;
  end

end

