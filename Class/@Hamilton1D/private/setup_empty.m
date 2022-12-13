function [HAM]=setup_empty(varargin)
% function [HAM]=setup_empty(L)
% e.g. called at the very beginning of each private/setup_*.m routines
% outsourced from Hamilton1D.m // Wb,May31,19
% EMPTY ; SETUP_EMPTY

  getopt('init',varargin);
     L=getopt('L',0);
     M=getopt('M',[]);
  getopt('check_error');

  HAM.info=struct('istr','','host',mlinfo,'IS',struct,'mpo',[]);

  HAM.oez=struct;
  HAM.ops=struct;

    q=struct(...
      'M',M, ...
      'iop',[],...
      'stype',1,...
      'user',[]...
    );

  HAM.mpo=repmat(q,1,L);

  HAM.store='';
  HAM.mat='';
  HAM.user=[];

end

