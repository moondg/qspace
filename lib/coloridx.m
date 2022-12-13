function coloridx(varargin)
% function coloridx(varargin)
% Wb,Oct25,20

  getopt('init',varargin);
     if getopt('--reset')
        set(gca,'ColorOrderIndex',1);
        hold on
     end
  getopt('check_error');

end

