(* demo.sml

   A tiny tour of the `Fft` library: transform a small real signal, print the
   per-bin magnitudes, and confirm the round-trip `ifft (fft x)` recovers the
   input. Build and run with `make example`. *)

structure C = Complex
structure F = Fft

val pi = Math.pi

fun line s = print (s ^ "\n")

(* Fixed-width-ish real formatting, identical across compilers. *)
fun fmt x =
  let
    val s = Real.fmt (StringCvt.FIX (SOME 4)) x
  in
    s
  end

(* A length-8 signal: DC offset + a cosine at bin 1 + a cosine at bin 3. *)
val n = 8
val signal =
  Array.tabulate
    (n, fn j =>
      let val t = real j
      in 1.0
         + 2.0 * Math.cos (2.0 * pi * 1.0 * t / real n)
         + 1.0 * Math.cos (2.0 * pi * 3.0 * t / real n)
      end)

val () = line "Input signal (real):"
val () =
  Array.appi
    (fn (j, v) => line ("  x[" ^ Int.toString j ^ "] = " ^ fmt v))
    signal

val spectrum = F.rfft signal

val () = line ""
val () = line "Spectrum magnitudes |X[k]|:"
val () =
  Array.appi
    (fn (k, z) => line ("  |X[" ^ Int.toString k ^ "]| = " ^ fmt (C.abs z)))
    spectrum

(* Round-trip back to the time domain. *)
val recovered = F.ifft spectrum
val maxErr =
  Array.foldli
    (fn (j, z, acc) =>
       Real.max (acc, Real.abs (C.re z - Array.sub (signal, j))))
    0.0
    recovered

val () = line ""
val () = line ("Round-trip max error (ifft (fft x) vs x): " ^ fmt maxErr)

(* A small linear convolution via the FFT: a moving sum of width 3. *)
val a = Array.fromList [1.0, 2.0, 3.0, 4.0, 5.0]
val box = Array.fromList [1.0, 1.0, 1.0]
val conv = F.convolve (a, box)

val () = line ""
val () = line "convolve([1,2,3,4,5], [1,1,1]) (moving sum, width 3):"
val () =
  line
    ("  ["
     ^ String.concatWith ", "
         (List.map fmt (Array.foldr (op ::) [] conv))
     ^ "]")
