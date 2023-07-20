% function contractQS(A,..B,..)
%
%     performs pairwise contraction of QSpace tensors
%     [see usage #1 and #2 below), while also properly taking care of
%     the underlying Clebsch Gordan coefficient spaces (if present).
%
%     Each QSpace can be used as is, or as its 'conjugate' where the
%     conjugate of a QSpace A, i.e. conj(A) is defined as the QSpace
%
%     a) with all `arrows' reversed
%     b) keeping the SAME qlabels [for reversing individual
%        arrows, see getIdentitQS(..,'-0')]
%     c) and complex conjugation of all data{*}
%        if applicable, i.e. a given QSpace is complex
%
%  Note that because of (a), the specification of conjugation flags
%  (`conj-flags') is also important for QSpaces with all-real matrix
%  elements.
%
%  Usage #1: S=contractQS(A, ica, B, icb [, perm, OPTS ]);
%
%     Plain contraction of a single pair of QSpaces, A and B,
%     with respect to given explicitly specified  sets of
%     contraction indices ica and icb, which can be specified
%
%       - numerically (e.g. [1 2]),
%       - or as strings (e.g. '1,2', or '1 2')
%       - or as compact strings (e.g. '1,2', or '12')
%
%     The last 'compact format' is only possible / unique,
%     of course, if the tensors A and B do not have more than
%     9 legs (which basically never occurs), such that the
%     contraction indices reamin in the single digits
%     (this can be further relaxed, though, by extending the
%     digital range to letter, i.e. using a-z after 0-9).
%
%     The recommended way to specify conj-flag with usage #1
%     is together with the contraction indices in string notation!
%     For example,
%
%       contractQS(A,[1 3],B,[1 3],'conjA')  is equivalent to
%       contractQS(A,'1,3;*',B,[1 3])        is equivalent to
%       contractQS(A,'13*',B,'13')
%
%     (Deprecated) options specific to usage #1:
%
%       'conjA'  use complex conjugate of A for contraction
%       'conjB'  use complex conjugate of B for contraction
%
%  Usage #2: S=contractQS({A,{B,{..,C}}},... [, perm, OPTS ]);
%
%     Generalized 'cell-contraction' of tensors: when encountering
%     a cell, the content of that cell is contracted first, before
%     using its results. This allows the specification of an entire
%     patter of pairwise contractions based on a nested cell structure
%     where the lowest-level contractions are performed first.
%     If an optional permutation [perm] is specified as an explicit
%     index array (non-string), it is applied to the final result only.
%
%     Cell contractions are furthermore based on QSpace 'itags'
%     i.e. string labels for indices with up to 8chars, and which
%     are specified in QSpace.info.itags. This offers automated
%     contraction ('auto-contraction') of pairs of tensors solely
%     based on matching itags! Uusage #2 therefore does not (also)
%     allow explicit specification of contraction indices as in
%     usage #1.  Therefore itags (plus conj-flags) must be unique.
%
%     Itags must always also contain individual conjugate flags
%     (this represent the bare minimum that must be specified with
%     each QSpacein v3): the conjugate flags on individual indices
%     (legs) of a tensor determine  whether that index (leg) is
%     in- or out-going, with the convention that
%
%         all out-going indices have a trailing * in their itags!
%
%     For example, an A-tensor with L(eft), R(ight), s(=local)
%     indices may have itags A.info.itags={'L','R*',s'}
%     assuming (L,R,s) index order.
%
%     In usage #2, for every operator additional optional strings
%     can be specified, appearing right after the affected tensor
%     e.g. QSpace A:
%
%       A,'!ij'  do not contract indices specified by ij (in compact
%                format) despite they share common matching itags.
%       A,'*'    apply overall (complex) conjugation on given input
%                tensor A (see early comments above)
%       A,'!ij*' both of the above in a single instruction
%                the conjugate flag '*' always trailing.
%
%       A,'ij'   
%       A,'ij*'  explicitly specify indices to contract
%                this is required only in the presence of degenerate itags
%                i.e., the case when identical itags that appear on multiple
%                legs including the same direction (conjugate flag).
%
%     Itags may be set or adapted on the fly (this is performed
%     prior to the auto-contraction together wit the above directions)
%
%       X,'--op:<tag>[:<opl=op>]'
%          specify itags for given (e.g. local) operator.
%
%          The last option considers X an operator, and hence assumes
%          operator itags '<tag>;<tag>*[;opl*]' for QSpace X;
%          the default operator label is `op', but may be changed
%          by specifying a trailing ':<opl>' as indicated above.
%          As a safeguard, this issues a warning, if existing itags are
%          overwritten. This is relevant e.g. for local operators that
%          are applied to a very specific site with associated itag.
%
%          [11/24/2018] the specified <tag> may now also represent
%          a regular expression (regexp), recognized by non-alphabetic
%          special characters, while ignoring conj-flags (without
%          special characters, the specified <tag> is taken as is!).
%          This usage then searches for a matching itag in the paired
%          up QSpace (cell) in the contraction.
%
%          Ex. Consider A some QSpace with a single local index
%          that starts with `s', e.g., like 's010' for site 10;
%          then contract(A,Xloc,'-op:^s') will autocontract
%          the local operator Xloc to the correct local index in A
%          (here the regex `^s' indicates `starts with s');
%          An alternative operator itag may still be specified
%          by adding a trailing ':opl' as indicated earlier.
%
%       A,'--itag:s/pat/rep/[gi]
%
%          replace/modify existing itags on the fly for a particular
%          recursive level of the cell contraction based on regular
%          expressions (regex) using ECMAScript grammar (cf. C++/regex).
%          The trailing flags enable case insensitive replacement [i]
%          and global replacement [g] of all possible matches.
%          The syntax is much analogous to perl regex.
%
%     An adaptation of usage #2 can also be used
%     for plain sequential contractions
%
%         S=contractQS(A [,flagsA],B [,flagsB],C,... [, perm, OPTS ])
%
%     which is equivalent to
%     S=contractQS({A [,flagsA],{B [,flagsB], {C,...}}}, [,perm,OPTS ]).
%     That is, by grouping A*(B*(C*...)), sequential contractions start
%     from the right end onwards to the beginning of the set.
%     Non-contracted indices are collected in the order they appear.
%
%  The remaining trailing OPTS are
%
%     perm  permutation to be applied to the final object;
%           NB! [06/02/2019] this permutation can be shorter
%           than the rank of the resulting QSpace; in this case
%           it only affects the leading range of indices.
%
%     '-v'  debug mode that shows level of cell contraction
%           together with actual contractions performed.
%           Internally, degenerate itags are frequently flagged
%           in order to make them unique and thus to differentiate them;
%           when printed, the flagged bits are formatted as <itag>⏐#
%           using the utf character `⏐' to indicate that the subsequent
%           number # is not part of the bare itag string.
%
%  Mixed usage of #2 and #1 is not possible.
%  Usage #2 is the typically recommended way because autocontraction
%  makes it far easier to perform entire contractions networks
%  without having to manually track and specify index locations.
%
%  AW (C) May 2010-2023

% -------------------------------------------------------------------- %
% CHANGE LOG:
% -------------------------------------------------------------------- %
% [07/14/2023] cell-contraction now also permits
%     explicit specification of index to contract for particular
%     QSpace on compact notation such as '12' equivalent to [1 2]
%     rather than just excluding indices that could be contracted
%     based on matching itags.
%
% [07/17/2023] introduced '--itag:s/pat/rep/[gi]
% [07/17/2023] replaced '-op: by '--op:'
%     yet with '-op:' still permitted for backward compatibility

