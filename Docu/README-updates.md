
The following contains a history log of more significant changes.
<br>Newest entries are shown first.

## Code update to tensor permute / HPTT

`[12/06/2024]` 
The permutation of indices of arbitrary n-dimensional arrays
was significantly updated with improved cache performance.
To be specific, this concerns the permute of the dense
arrays in X.data{} for some QSpace X, may it be implicit
within mex-routines as in contractions,
or explicit by using `permuteQS`.

Alternatively, QSpace now also permits to link the 
**`high-performance tensor transpose (HPTT)`** library.
When that library is compiled locally and the environmental
variable HPTT_ROOT is set to point to the root path
of that library, then the Makefile will link `-lhptt`
and use it for all tensor permutes, instead (this assumes
that paths $(HPTT_ROOT)/include/ [hptt.h] and $(HPTT_ROOT)/bin
[libhptt.so or similar] exist).
Thanks to Dai-Wei Qu and Wei Li (ITP Chinese Academy of Sciences, Beijing)
and Changkai Zhang (LMU Munich) for pointing out this library.


## Git-repo QSpace-v4-osx12-Monterey retired

`[10/17/2024]` 
The git repository `QSpace-v4-osx12-Monterey`
compiled for (Intel-based) macOS  Monterey has been retired,
and will no longer be maintained. It is replaced
by the new git repository **`QSpace-v4-osx15-Sequoia`**
which is compiled on the macOS ARM architecture
M3 within the Rosetta environment.

