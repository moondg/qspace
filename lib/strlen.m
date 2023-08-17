function [n,dN]=strlen(s)
% function [n,dn]=strlen(s)
%
%    Return length n of string s.
%
%    Here n represents the actual number of characters,
%    being aware of multibyte utf-characters, etc.
%
%    The offset due to utf-like characters is returned as dn.
%    E.g., dn must add to string format width specifications
%    such as '%##s' in fprintf commands, in order to the
%    get the correct width shown.
%
% Wb,Aug15,23

% http://www.zedwood.com/article/cpp-utf8-strlen-function
% see also clab/c.c // Wb,Aug15,23


   n=length(s); dN=0;

   cref=32768; % 2^15 32768 / last bit set
   if all(s<cref), return; end % no special character expected

   qm=256*uint32([
       hex2dec('0xE0') hex2dec('0xC0')
       hex2dec('0xF0') hex2dec('0xE0')
       hex2dec('0xF8') hex2dec('0xF0')
   ]);

   i=1; m=size(qm,1);

 % NB! matlab has 2-byte representation for each character
   while i<=n, c=uint32(s(i)); dn=1;
      if c>255
         for j=1:m
            if bitand(c,qm(j,1)) == qm(j,2)
             % largest set of bits set to '1110 ...' indicates
             % number of characters required
               dn=1+j; dN=dN+j; break
            end
         end
      end
      i=i+1+dn;
   end

   n=n-dN;

end

