function dd = field2arr(cc, fname, sub)
% Function: field2arr
%    extract specific NUMERIC field of array
%    preserving the dimensions of the 
%
% Usage: dd = field2arr(stru_arr, fieldname [, sub])
%
%    stru_arr    array of specific structure
%    fieldname   field name
%    sub         data sub specification, e.g. ':,1' selects the first column
%
% Wb,Nov15,05

  if nargin<2, eval(['help ' mfilename]); return; end

  if nargin==2
     eval(sprintf('dd=cat(1,cc.%s);',fname));

     if prod(size(dd))==prod(size(cc))
        dd=reshape(dd,size(cc));

        return
     else
        sz=size(cc); dd=nan(sz);
        cmd=sprintf('if ~isempty(cc(i).%s), dd(i)=cc(i).%s; end',fname,fname);
        for i=1:prod(sz), eval(cmd); end
     end
  else
     sz=size(cc); dd=zeros(sz);
     cmd=sprintf('dd(i)=cc(i).%s(%s);',fname, sub)

     for i=1:prod(sz)
        try
           eval(cmd);
        catch l
           cmd
           wblog(1,'ERR','Invalid sub index specification `%s''', sub);
           rethrow(l);
        end
     end
  end

return

