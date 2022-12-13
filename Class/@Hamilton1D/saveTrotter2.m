function saveTrotter2(HAM,HMX,tag)
% function saveTrotter2(HAM,HMX,tag)
% Wb,Aug15,19

% outsourced from tst_tdDMRG.m // Wb,Aug15,19

  if isempty(tag), error('Wb:ERR','\n   ERR invalid usage'); end

  m=[HAM.mat tag];
  wblog('I/O','saving TR2 to %s',repHome(m));

  for k=0:numel(HAM.mpo)
	 if k
		  S=HMX.user.trotter.data(k);
	 else S=rmfield(HMX.user.trotter,'data');
	 end
	 save2(sprintf('%s_%03g.mat',m,k),'-q','-f','-struct','S');
  end

end

