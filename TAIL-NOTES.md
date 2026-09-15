# The tail leg: sieving a range's last bit arrays without padding them to a whole chunk

Branch `tail-leg` off `v2.3` at 9ca70c8, 2026-09-15.  The review's P2
(`review-2026-09-12/REVIEW-REPORT.md`; verdict F9 in `verdicts-final.md`,
both skeptics confirmed and corrected it), the first item of the sieve
group; TODO item 25.  Measured on the i7-1355U, pinned to P-core 0, in core
cycles, paired and interleaved with the base, medians of the per-round
ratios, outputs compared with the references before anything was timed
(`pair.sh` in `review-2026-09-12/`).

## What was wrong

`sift()` (find_points.c) cut the numerator interval of a denominator into
blocks of `array_size` bit arrays and rounded the last block up to a
multiple of `RATPOINTS_CHUNK` (16): `_ratpoints_sift0` then sieved the
padding with every prime of the first phase, the scan of the second phase
walked it, and a loop after the first phase zeroed it.  The verdicts
counted the padding at 14.2 per cent of all bit arrays swept in `make
test1` (height 16383; 11.2 in test1many, 13.3 in testdegrees), 80 per cent
at height 1000, 90 at height 200, 1.4 at 200000: a block at 16383 is 41
bit arrays long on average, so the last chunk is mostly padding, and at
height 1000 the whole interval of a denominator is 3 bit arrays in a chunk
of 16.

## What was done

(filled in as the steps land)

## Measurements

(filled in from the chain)
