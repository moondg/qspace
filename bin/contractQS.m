% Usage: contractQS(A,..B,..)
%
%     performs pairwise contraction of tensors
%     e.g. A and B in QSpace format while also taking care of the
%     underlying Clebsch Gordan coefficient spaces (if present).
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
%  Usage #2: S=contractQS({A,{B,C}},... [, perm, OPTS ]);
%
%     Generalized 'cell-contraction' of tensors: when encountering
%     a cell, the content of that cell is contracted first, before
%     using its results. This allows the specification of an entire
%     patter of pairwise contractions based on nested structures
%     where the lowest-level contractions are performed first.
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
%       A,'!ij*' both of the above (with * always trailing).
%
%       X,'-op:<tag>[:<opl=op>]' specify itags for given (e.g. local) operator.
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
%     An adaptation of usage #2 also be used for a plain sequential
%     contractions where
%         S=contractQS(A,B,C,... [, perm, OPTS ]); is equivalent to
%         S=contractQS({A,{B, {C,...}}}, [,perm,OPTS ]),
%     i.e. the sequential contractions are started from the end
%     onwards to the beginning of the set.
%
%  The remaining trailing OPTS are
%
%     perm  permutation to be applied to the final object;
%           NB! [06/02/2019] this permutation can be shorter
%           than the rank of the resulting QSpace; in this case
%           it only affects the leading range of indices.
%
%     '-v'  verbose mode that shows level of cell contraction
%           together with actual contractions performed.
%
%  Mixed usage of #2 and #1 is not possible.
%
%  Note that mex files do not allow to return class objects,
%  hence S is returned as a QSpace structure. To get a class
%  object, use QSpace(contractQS(...)), or equivalently,
%  if contract() is properly defined as a wrapper routine
%  within matlab's Class/@QSpace (see MPS Pack), having a QSpace
%  input A, this may be shortended to contract(A,...).
%
%  AW (C) May10 ; Aug12 ; Dec14 ; Sep16
