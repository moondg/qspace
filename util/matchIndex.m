%  Usage: [Ia,Ib,I] = matchIndex(QA, QB [,OPTS])
%
%     QA or QB are matrizes whos rows are matched; consequently, 
%     their number of columns (row dimension) must be equal.
%
%  Options
%
%      '-s'    sort indizes with respect to A
%      '-q'    quiet mode, i.e. do not issue warning when NaN's
%              are encountered.
%      'eps',..specify epsilon within which to accept match (0)
%      '-f'    use float instead of double (gets rid of numerical
%              noise in case of double precision input data)
%
%  Output
%      Ia      index into QA (QA.Q)
%      Ib      index into QB (QB.Q)
%
%      I       an info structure with elements
%      .n      number of matches [=length(Ia)]
%      .dim    number of input records [ size(QA,1), size(QB,1) ]
%      .ix1    entries in QA not matched in QB
%      .ix2    entries in QB not matched in QA
%
%  AW (C) 2006 ; Feb10,11
