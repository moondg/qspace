function s=num2text(n)

  if n<0 || n~=round(n)
     wberr('%s requires uint number', mfilename); end

  switch n
     case  0, s='zero';
     case  1, s='one';
     case  2, s='two';
     case  3, s='three';
     case  4, s='four';
     case  5, s='five';
     case  6, s='six';
     case  7, s='seven';
     case  8, s='eight';
     case  9, s='nine';
     case 10, s='ten';
     case 11, s='eleven';
     case 12, s='twelve';
     otherwise, s=sprintf('%d',n);
  end

return

  switch n
    case 13, s='thirteen';
    case 14, s='fourteen';
    case 15, s='fifteen';
    case 16, s='sixteen';
    case 17, s='seventeen';
    case 18, s='eighteen';
    case 19, s='nineteen';
  end

