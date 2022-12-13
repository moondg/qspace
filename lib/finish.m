% function finish()
% ------------------------------------------------------------------- %
% NB! finish.m is called by matlab command quit e.g. see 'help quit'
% NB! do not make this a function (e.g. see 'save' below!)
% ------------------------------------------------------------------- %
% NB! to not keep this in $ML root directory, as the latter will
% not be included in MCC path (hence this would not find finish.m)
% Wb,Nov16,15
% ------------------------------------------------------------------- %

% if exist('batchmode') if batchmode, return; end; end

% eval(sprintf('! pstime.pl %g',getpid));
% system(sprintf('pstime.pl -Vml %g',getpid)); % Wb,Mar28,14

  fprintf(1,'\n>> Exiting MatLab (finish.m) ...\n>> Bye.\n');

return

  button = questdlg(['Howdy from your finish.m script' 10 'Ready to quit?'], ...
                    'Exit Dialog','Yes','No','Save','Yes');
  switch button
    case 'Yes',
    case 'No',
      quit cancel;
    case 'Maybe'
      printfc('\nwhy? '); why; printfc('\n');
      quit cancel;
    case 'Save'
      save 
  end

