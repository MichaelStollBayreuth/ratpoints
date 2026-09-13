# The bugs the single-thread review found, and their fixes

Branch `review-bugfixes` off `v2.3` at ca5cc45, 2026-09-13.  TODO item 18.
The bugs come from the review of 2026-09-12 (`review-2026-09-12/REVIEW-REPORT.md`
in the repository root, untracked, section 1; the verdicts of the two skeptics
per finding are in `verdicts-final.md` there).  This file records what was
done, what was decided along the way, and how it was checked.  It stays on
the branch; the merge strips it.

## The five bugs

**B1. Degree 1 read past the Sturm arrays** (`sturm.c`).  The division loop
`for(k = 2; k <= degree; k++)` never runs for degree 1, so `k` is left at 2
and the sign-count loop reads `sturm_degs[2]` and `sturm[2][...]` out of
arrays of length 2; SIGSEGV every time.  Reached by every `y^2 = ax + b`
and by every `y^2 = c2 x^2 + c1 x` (reversal case 1 lowers the degree) or
`c2 = 0`.  Fix: `if(k > degree) k = degree;` after the loop.  For degree
2 and up the loop always leaves through the `d2 == 0` break with
`k <= degree` (`sturm_degs` strictly decreases and `sturm_degs[degree]` is
forced to 0), so the clamp is dead code there; for degree 1 it gives the
chain `f, f'`, one real root, and `ivlocal[1 + (degree>>1)]` has the right
size because `{f > 0}` is a half-line.  Present in 2.2.3.

**B2. Points at infinity dropped when the positivity region misses the
search domain** (`find_points.c`, the Sturm step).  `_ratpoints_compute_sturm`
returns 0 both for "f negative everywhere" and for "the intersection with
the domain is empty", and `find_points_work` returned 0 for both before it
had looked at the points at infinity.  The first case is harmless (even
degree, negative leading coefficient, no point at infinity), the second
loses the point of every odd-degree curve whose real branch starts beyond
the height bound (`y^2 = x^3 + k` with `|k|` above about `H^3`), two
points of an even-degree curve with a square leading coefficient under
`-k`/`-l`/`-u` (with reversal allowed the case cannot arise), and after a
reversal the lost point is affine: `./ratpoints '0 1 0 0 -1000000000' 45`
printed nothing where the curve has `(0 : 0 : 1)`.  Fix: split the test on
the return value; for 0, return at once only when there is no point at
infinity (`!((degree & 1) || lcfsq)`), otherwise set `sturm_empty` and
return after the point-at-infinity block.  That placement matters: the
finder's version (always fall through) ran the whole set-up on the 90
test1 curves that stop here, +0.11% of the suite.  The verbose message no
longer claims the polynomial is always negative.  Present in 2.2.3.

**B3. `num_bits[16]` read uninitialised** (`get_2adic_info`).  The early
return for `db == 0` (no residue class of the denominator mod 16 admits a
numerator) skipped the block that fills `num_bits`; the caller then read
all sixteen entries.  By default the answer was still right (`db == 0`
forces the Jacobi test on, which keeps the checked denominator loop, and
that loop sifts nothing when `den_bits == 0`), but the chosen parameters
and the `-v` output were nondeterministic, and with `-j -F 0` (or the API
equivalents, both public) the plain loop tests nothing but `num_bits` and
printed nondeterministic points outside the height bound after a run 5000
to 10000 times too long.  Two fixes: `num_bits` is zeroed before the early
return, and `find_points_work` returns right after `get_2adic_info` when
`den_bits == 0` -- there is no affine point then, and none at infinity
either, since an odd degree always leaves the class `v_2(b) >= 4`
admissible, so the degree is even and the leading coefficient is a
non-square mod 16.  (The return is guarded by `!point_at_infty` anyway.)
Present in 2.2.3 (the first half).

**B4. `run_shape` put U on its floor for `num_none`** (`run_shape`).
`which_bits == num_none` only says that no *odd* denominator has an
admissible numerator; the even ones are still sieved, with `num_odd`
forced in `sift()`.  `nums = 0.0` made `u_words` 1 while `n_denom` stayed
right, so `phase_2_offset` gave 0 (no second stage), the table term of
`prime_key` became of order 1e6 (phase 1 took the smallest primes), the
stage-3 rule saw S near 0 (no third stage), and `adapt_primes` returned at
once on `sp2 <= sp1` (no correction).  41 of the 900 test1 curves that
sieve; 25 to 40 per cent on each of them, -0.34% of test1 and -0.66% of
testhigh in instructions.  Fix: delete the special case; `nums *= 0.5` is
exact for these curves (every sieved denominator is even).  New in 2.3
(`run_shape` is new), so nothing to backport.

**B5. Input fields used as outputs** (`find_points_work`, `sturm.c`).
The review named `max_forbidden` (the number of forbidden divisors found
was written into it; only a negative value is reset to the default, so
for a caller that fills `args` once and loops over curves the test latched
off after the first curve without forbidden divisors: 1.9x slower at
16383, 2.3x at 200000) and `num_inter` (the intersected interval count was
written back, so the next curve's domain was intersected with the previous
curve's positivity region: 26% of the test1 curves lose points, silently).
Looking for the whole class turned up more: the *chosen* `sp1` and `sp2`
were written back (so the second curve of such a loop ran with the first
curve's primes, and with adaptation off, since an explicit value disables
it), `num_primes` was normalised in place (which turns the default into a
hard limit), `b_low`, `b_high`, `array_size`, `sturm` were normalised in
place, and the derived "do not reverse" decision was stored in the caller's
`RATPOINTS_NO_REVERSE` bit.  The documentation's loop example suggests
exactly the reuse that trips over all of this; `rptest.c` and `main.c`
happen to refill everything per call, which is why nothing had noticed.

Fix, in three parts.  (1) `find_points_work` is now a wrapper: it saves the
input fields and the caller's `num_inter` intervals, calls the search
(`find_points_work_1`), and puts them back; what the search used is
reported in three new output fields `sp1_used`, `sp2_used`, `sp3_used`
(`sp3` becomes a working field; it had only been documented as an output
since 2026-09-10, unreleased).  Left as the search made them, on purpose:
`cof` and `degree` when the polynomial was reversed or a leading zero
dropped (`RATPOINTS_REVERSED` reports it; the pair still describes the
polynomial worked with), and the flag bits that report on the run.
(2) The derived no-reversal condition sets a private bit
`RATPOINTS_NO_REVERSE_AUTO`, cleared on entry by the input mask, and the
two tests check both bits.  (3) The write into `max_forbidden` is deleted
(nothing read the value).  `main.c` is unaffected: it prints the search
intervals from the caller's copy *before* the call.  The documentation
now says which fields are inputs, that they come back unchanged, and what
the two exceptions are.  For 2.2.4 the same wrapper (without the new
fields, which would break the ABI) and the same documentation note.

## Hygiene from the finders' cut list

Four left shifts of negative values, `i << RBA_SHIFT` and
`i << (RBA_SHIFT+1)` in `sift.c` (twice each) and `nr << (del-der)`,
`nl << (der-del)` in `sturm.c`, are multiplications now (same instruction;
these were the only UBSan reports on the four suites).  `sieve_spec
ssp[args->sp2]` and `long bp_list[args->sp3_max]` (four sites) had length 0
under `-n 0 -N 0` / `-p 1`, which is undefined; they get at least 1.  The
divisor enumeration `t = *div0 * p; if(t <= b_high)` compares
`*div0 <= b_high / p` instead, so it cannot wrap for `b_high` above about
4e15.  `RATPOINTS_DEFAULT_NUM_PRIMES` and `RATPOINTS_ARRAY_SIZE` got the
`#ifndef` guard their neighbours have.  Not done: `sieves0` is still not
`const` (it is assigned to non-const pointers; no runtime effect), and the
"2^60" comment the finders mention could not be found.

## Tests

`test3.sh` / `testbase3` / `make test3`: fifteen invocations of
`./ratpoints`, one or more per bug, with `-q` so the output is the point
list only.  Every point list was checked against an independent
brute-force search in PARI/GP over all coprime `(a, b)` with `|a|, b <= H`
(`scratchpad/brute.gp`: `F(a,b) = sum c_k a^k b^(e-k)`, `e` the degree
rounded up to even, `issquare`, plus the points at infinity by hand);
identical sets in all ten cases with points.  The `-k` and `-l -u` cases
print the two points at infinity as expected, the `num_bits` cases print
nothing (and the 300000 one now runs in a millisecond), the `run_shape`
curve's points are the ones in `testbase`.

`rptest -O` / `make test1once`: `rptest` fills the input fields once before
the loop over the 1008 curves (only `cof` and `degree` per curve) and, after
every call, prints a line for each input field whose value differs from
what was set; `make test1once` compares the output with `testbase`.  Zero
lines, output identical.

The existing suites: `test1`, `test1many`, `testdegrees`, `test2`,
`testhigh`, `testhighmany` byte-identical to their references, on the
default build and, for `test1`, `test1once`, `test3`, `testdegrees`, on the
128-bit, 64-bit, SSE, `RATPOINTS_CHUNK=1` and `USE_LONG_IN_PHASE_2` builds
(scratch copy).  valgrind on the debug build: 0 errors on eight
reproducers including the former uninitialised read.

## Backport to 2.2.4

B1, B2, B3 (both halves), B5 (wrapper without new fields, the private
reversal bit, the dead write), the four shifts, the divisor overflow; not
B4 (no `run_shape`), not the VLAs (2.2.3 has `ssp[args->sp2]` too -- check),
not the header guards (`RATPOINTS_ARRAY_SIZE` exists there, check).  Plus
`test3` (the same script and reference; 2.2.3's output for these must be
checked to be the same point lists), the documentation (`ratpoints-doc-2.2.tex`
change log and the field semantics), the version strings.
