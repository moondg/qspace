%  Usage: [kk,I] = wbhash(dd,N [,opts]);
%
%     Compute keys kk using given hash for each row of the
%     input array using given hash specified by the hash
%     function (uses offset and multiplier) together with
%     the final dimension N of the hash key (hash mod N).
%     The hash function is given by
%
%         unsigned long h=offset;
%         for (i=0; i<len; ++i) {
%            h = mult*h + s[i];
%         }
%
%     Here s represents a row in dd of length len.
%     Default values for the hash function are
%
%         offset = 5381 (prime)
%         mult = 33 = 2^5 + 1;
%
%     which correspond to the hash function 'djb2'
%     as suggested by Daniel B. Bernstein Having mult=33,
%     the hash value is evaluated through h = (h<<5)+h+s[i].
%
%  Options
%
%   For algorithm based on multiplier
%     'mult',..   multiplier (33)
%     'offset',.. offset (5381)
%
%   Alternatively the following flags use default string hashes
%     '-stl'      uses hash<std::string> string_hash(char(dd)).
%     '-wbb'      adapted versions of algorithms in boost/hash.cpp
%     '-ump'      test mode for unordered_map
%
%  Wb,Jul01,14
