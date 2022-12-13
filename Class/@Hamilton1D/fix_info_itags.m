function fix_info_itags(HAM,kk)
% function fix_info_itags(HAM,kk)
% fix in case X.info.itags was missed to update
% Wb,Mar03,17

  for k=kk
     fprintf(1,'\r   k = %2g  \r',k);
     q=load_dmrg_data(HAM,k);

   % ----------------------------- %
     q.info.itags=q.AK.info.itags;
     save_dmrg_data(HAM,q,k);
  end
  fprintf(1,'\r%60s\r','');

end

