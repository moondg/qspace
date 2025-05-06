The following contains a history log of more significant changes.
with newest entries shown first.

## (re)reading ENV variables in QSpace mex-routines 

`[04/15/2025]`
Previously QSpace mex-routines checked / re-read QSpace specific
environmental variables with every single call to the mex-routine.
This has been reduced to re-checking at most once per second
which thus is only relevant if a particular mex-routine is called
rapidly many times within a second. The change was motivated
by an observed slow-down in this regard for Linux environments.
Acknowledgment Seung-Sup Lee / Geng-Dong Zhou (Seoul University).


## Added symmetries  `Aspin(:)` and `SU2spin(:)`

`[03/03/2025]`
The spin symmetry now can also be specified individually
across channels (rather than just the total spin), if the
inter-channel interaction acts trivially on the spin sector.
For this purpose, `getLocalSpace` now permits adding a
trailing `(:)` string to the `Aspin` or `SU2spin` symmetry specification.
Note that the `(:)` flag was already implicit for the case of
a non-abelian channel symmetry like SU(N).
Thanks to Geng-Dong Zhou (Peking University, Beijing) and Seung-Sup Lee 
(Seoul National University, South Korea) for pointing this out.


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

