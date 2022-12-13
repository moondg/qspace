function [HAM,t]=acceptTrotter2Step(HAM)
% function [HAM,t]=acceptTrotter2Step(HAM)
% Wb,Jun04,19

  L=length(HAM.mpo);

  for k=L:-1:1
	 Ir=load_trotter_data(HAM,k);

	 Ir.At=itagrep(Ir.Af,'F','K');

	 ut=HAM.user.trotter.mpo(k);
	 if k<L
		  Q={ut,Xm}; if k>1, Ir.Xm=Xm; end
	 else Q=ut; end
	 Xm=contract(Ir.Af,'!1*',{Ir.At,Q}); % RR',m(po)

	 HAM=save_trotter_data(HAM,Ir,k);
  end

  q=HAM.user.trotter.info; 
  t = q.ttot + q.dt;
  HAM.user.trotter.info.ttot = t;

end

