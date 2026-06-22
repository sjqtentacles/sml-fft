# sml-fft

Fast Fourier transforms in pure Standard ML — `fft`/`ifft` over complex
arrays, a real-input `rfft`, and FFT-based linear `convolve` — built on
[`sml-complex`](https://github.com/sjqtentacles/sml-complex). Power-of-two
lengths use an iterative **radix-2 Cooley-Tukey** transform; every other
length uses **Bluestein's** chirp-z algorithm, so all sizes run in
`O(n log n)`. No FFI, no external dependencies, and **deterministic**,
byte-identically under both [MLton](http://mlton.org/) and
[Poly/ML](https://www.polyml.org/).

## Status

- 73 assertions, green on MLton and Poly/ML.
- Basis-library only; deterministic across compilers.
- Vendors `sml-complex` (Layout B), so the repo builds standalone.

## Install

With [`smlpkg`](https://github.com/diku-dk/smlpkg):

```
smlpkg add github.com/sjqtentacles/sml-fft
smlpkg sync
```

Include the MLB from your own (it pulls in the vendored `sml-complex`):

```
local
  $(SML_LIB)/basis/basis.mlb
  lib/github.com/sjqtentacles/sml-fft/... (via smlpkg)
in
  ...
end
```

This brings `structure Fft` (and the vendored `Complex`) into scope.

## Quick start

```sml
val signal = Array.fromList [1.0, 2.0, 3.0, 4.0, 5.0, 0.0, ~1.0, 2.0]

val spectrum = Fft.rfft signal            (* Complex.t array, length 8 *)
val mags     = Array.map Complex.abs spectrum

(* round-trips up to floating-point rounding *)
val back = Fft.ifft (Fft.fft (Array.map (fn x => Complex.complex (x, 0.0)) signal))

(* linear convolution via the FFT: a moving sum of width 3 *)
val conv = Fft.convolve (Array.fromList [1.0, 2.0, 3.0, 4.0, 5.0],
                         Array.fromList [1.0, 1.0, 1.0])
(* [1, 3, 6, 9, 12, 9, 5] *)
```

## API (`signature FFT`)

```sml
val fft      : Complex.t array -> Complex.t array   (* forward DFT          *)
val ifft     : Complex.t array -> Complex.t array   (* inverse, 1/n-scaled  *)
val rfft     : real array -> Complex.t array        (* DFT of a real signal *)
val convolve : real array * real array -> real array (* linear convolution  *)
```

The forward transform is unnormalized:

```
fft  x : X[k] = sum_{j=0}^{n-1} x[j] * exp(-2*pi*i*j*k/n)
ifft X : x[j] = (1/n) sum_{k=0}^{n-1} X[k] * exp(+2*pi*i*j*k/n)
```

so `ifft (fft x) = x` up to rounding. `convolve (a, b)` returns the linear
convolution of length `length a + length b - 1` (empty if either input is
empty). An empty array transforms to an empty array.

### More transforms

```sml
val irfft : Complex.t array -> real array              (* inverse real FFT     *)
val fft2  : Complex.t array array -> Complex.t array array  (* 2D forward DFT   *)
val ifft2 : Complex.t array array -> Complex.t array array  (* 2D inverse DFT   *)
val dct   : real array -> real array                   (* DCT-II (unnormalized)*)
val idct  : real array -> real array                   (* DCT-III, its inverse *)
```

- `irfft` is the real part of `ifft`, so `irfft (rfft x) = x` up to rounding.
- `fft2`/`ifft2` treat a `Complex.t array array` as a row-major matrix and
  transform along rows then columns (separable). `ifft2` is
  `1/(rows*cols)`-normalized, so `ifft2 (fft2 m) = m`. Rectangular and empty
  matrices are handled.
- `dct` is the unnormalized type-II transform
  `X[k] = sum_j x[j] cos(pi*(2j+1)*k/(2n))`; `idct` is the type-III transform
  scaled by `1/n`, the exact inverse, so `idct (dct x) = x`. For example
  `dct [1,1,1,1] = [4,0,0,0]`.

### Conventions

- All functions are pure: each copies its input and returns a fresh array;
  inputs are never mutated.
- Power-of-two lengths take the in-place radix-2 path; other lengths take the
  Bluestein path (a single power-of-two convolution of a chirped signal).
- `ifft` is the conjugate trick `conj (fft (conj X)) / n`, so the forward
  kernel is the only numeric code.
- Floating-point results are compared with an explicit epsilon in the tests
  (`approx`/`approxC`, `eps = 1e-9`) — never string- or structural-equality,
  since `Real.toString` differs between compilers and FFTs only recover values
  up to rounding.

## Build & test

```
make test        # MLton
make test-poly   # Poly/ML
make all-tests   # both
make example     # build + run examples/demo.sml
make clean
```

Both compilers run the same strict-TDD suite (`test/test.sml`), whose oracle is
a naive `O(n^2)` DFT and a naive linear convolution. Highlights:

- **DFT equivalence:** `fft` matches the naive DFT for sizes `1..16`, including
  non-power-of-two lengths (`3, 5, 6, 7, 12`) that exercise Bluestein.
- **Inversion:** `ifft (fft x) ~= x` (and `fft (ifft X) ~= X`).
- **Linearity:** `fft(a*x + y) = a*fft(x) + fft(y)`.
- **Known transforms:** a constant signal maps to a delta at bin 0; a pure
  cosine maps to peaks of magnitude `n/2` at bins `±k`.
- **Parseval:** `sum |x[j]|^2 = (1/n) sum |X[k]|^2`.
- **Convolution theorem:** `convolve (a, b)` matches the direct convolution.

### Poly/ML note

CI builds Poly/ML 5.9.1 from source rather than using the Ubuntu package
(Poly/ML 5.7.1), whose X86 code generator crashes (`asGenReg raised while
compiling`) on heavy real-arithmetic code. See `.github/workflows/ci.yml`.

## License

MIT — see [LICENSE](LICENSE).
