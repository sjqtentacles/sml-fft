(* test.sml

   Strict-TDD suite for `Fft`. The transform is built on floating-point
   `real`/`Complex.t`, so every comparison goes through an explicit epsilon
   tolerance (`approx`/`approxC` and their array lifts) rather than string or
   structural equality -- `Real.toString` differs between MLton and Poly/ML,
   and FFT round-trips only recover the input up to rounding. The same `eps`
   is shared by every check so failures are legible and identical across both
   compilers.

   The oracle for correctness is a naive O(n^2) DFT (`naiveDft`) and a naive
   O(n*m) linear convolution (`naiveConv`); the fast transform is required to
   agree with them within `eps`. *)

structure Tests =
struct

  structure C = Complex
  structure F = Fft

  val eps = 1E~9
  val pi  = Math.pi

  (* ---- tolerance helpers ---- *)
  fun approx (a, b) = Real.abs (a - b) <= eps
  fun approxC (z, w) =
    approx (C.re z, C.re w) andalso approx (C.im z, C.im w)

  (* Elementwise complex-array comparison within eps (lengths must match). *)
  fun approxArrC (xs, ys) =
    Array.length xs = Array.length ys andalso
    let
      fun loop i =
        i >= Array.length xs orelse
        (approxC (Array.sub (xs, i), Array.sub (ys, i)) andalso loop (i + 1))
    in
      loop 0
    end

  (* Elementwise real-array comparison within eps (lengths must match). *)
  fun approxArrR (xs, ys) =
    Array.length xs = Array.length ys andalso
    let
      fun loop i =
        i >= Array.length xs orelse
        (approx (Array.sub (xs, i), Array.sub (ys, i)) andalso loop (i + 1))
    in
      loop 0
    end

  (* ---- builders ---- *)
  fun carr xs = Array.fromList (List.map (fn (a, b) => C.complex (a, b)) xs)
  fun rcarr xs = Array.fromList (List.map (fn a => C.complex (a, 0.0)) xs)
  fun rarr xs = Array.fromList (xs : real list)

  (* ---- oracles ---- *)

  (* Naive O(n^2) forward DFT: X[k] = sum_j x[j] * exp(-2*pi*i*k*j/n). *)
  fun naiveDft (x : C.t array) =
    let
      val n = Array.length x
      fun bin k =
        let
          fun loop (j, acc) =
            if j >= n then acc
            else
              let
                val ang = ~2.0 * pi * real (k * j) / real n
                val w = C.complex (Math.cos ang, Math.sin ang)
              in
                loop (j + 1, C.add (acc, C.mul (Array.sub (x, j), w)))
              end
        in
          loop (0, C.complex (0.0, 0.0))
        end
    in
      Array.tabulate (n, bin)
    end

  (* Naive linear convolution of two real sequences (length la+lb-1). *)
  fun naiveConv (a : real array, b : real array) =
    let
      val la = Array.length a
      val lb = Array.length b
      val outLen = la + lb - 1
      fun out k =
        let
          fun loop (j, acc) =
            if j >= la then acc
            else
              let val kj = k - j
              in
                if kj >= 0 andalso kj < lb
                then loop (j + 1, acc + Array.sub (a, j) * Array.sub (b, kj))
                else loop (j + 1, acc)
              end
        in
          loop (0, 0.0)
        end
    in
      Array.tabulate (outLen, out)
    end

  (* Naive separable 2D DFT oracle: naive DFT along rows, then columns. *)
  fun naiveDft2 (m : C.t array array) =
    let
      val rows = Array.length m
    in
      if rows = 0 then Array.fromList []
      else
        let
          val cols = Array.length (Array.sub (m, 0))
        in
          if cols = 0 then Array.tabulate (rows, fn _ => Array.fromList [])
          else
            let
              val rowT = Array.tabulate (rows, fn r => naiveDft (Array.sub (m, r)))
              val result =
                Array.tabulate (rows, fn _ =>
                  Array.array (cols, C.complex (0.0, 0.0)))
              fun doCol c =
                if c >= cols then ()
                else
                  let
                    val col =
                      naiveDft (Array.tabulate (rows, fn r =>
                        Array.sub (Array.sub (rowT, r), c)))
                    fun writeBack r =
                      if r >= rows then ()
                      else
                        (Array.update (Array.sub (result, r), c,
                                       Array.sub (col, r));
                         writeBack (r + 1))
                  in
                    writeBack 0;
                    doCol (c + 1)
                  end
            in
              doCol 0;
              result
            end
        end
    end

  (* Elementwise comparison of two matrices (array of rows) within eps. *)
  fun approxArrArrC (xs, ys) =
    Array.length xs = Array.length ys andalso
    let
      fun loop i =
        i >= Array.length xs orelse
        (approxArrC (Array.sub (xs, i), Array.sub (ys, i)) andalso loop (i + 1))
    in
      loop 0
    end

  (* Naive DCT-II oracle (direct sum): X[k] = sum_j x[j] cos(pi*(2j+1)*k/(2n)). *)
  fun naiveDct (x : real array) =
    let
      val n = Array.length x
      fun bin k =
        let
          fun loop (j, acc) =
            if j >= n then acc
            else
              let
                val ang = pi * real (2 * j + 1) * real k / real (2 * n)
              in
                loop (j + 1, acc + Array.sub (x, j) * Math.cos ang)
              end
        in
          loop (0, 0.0)
        end
    in
      Array.tabulate (n, bin)
    end

  (* Energy helpers for Parseval. *)
  fun energyR (x : real array) =
    Array.foldl (fn (v, acc) => acc + v * v) 0.0 x
  fun energyC (x : C.t array) =
    Array.foldl (fn (z, acc) => acc + C.re z * C.re z + C.im z * C.im z) 0.0 x

  (* A fixed, deterministic real test signal of length n. *)
  fun signal n =
    Array.tabulate
      (n, fn j =>
        let val t = real j
        in 0.5 + Math.sin (0.7 * t) + 0.3 * Math.cos (1.9 * t + 0.4) end)

  fun runAll () =
    let
      (* ---- fft == naive DFT, several sizes incl. non-power-of-two ---- *)
      val () = Harness.section "fft equals naive DFT (within eps)"
      val () =
        List.app
          (fn n =>
            let
              val x = signal n
              val xc = Array.tabulate (n, fn j => C.complex (Array.sub (x, j), 0.0))
            in
              Harness.check
                ("n = " ^ Int.toString n)
                (approxArrC (F.fft xc, naiveDft xc))
            end)
          [1, 2, 3, 4, 5, 6, 7, 8, 12, 16]

      (* ---- ifft (fft x) ~= x ---- *)
      val () = Harness.section "ifft (fft x) = x (within eps)"
      val () =
        List.app
          (fn n =>
            let
              val xc = Array.tabulate (n, fn j =>
                         C.complex (Array.sub (signal n, j),
                                    Math.cos (0.3 * real j)))
            in
              Harness.check
                ("round-trip n = " ^ Int.toString n)
                (approxArrC (F.ifft (F.fft xc), xc))
            end)
          [1, 2, 3, 4, 6, 8, 9, 16]

      (* ---- fft (ifft X) ~= X as well ---- *)
      val () =
        let
          val xc = carr [(1.0, 0.0), (2.0, ~1.0), (0.0, 3.0), (~1.0, 0.5),
                         (2.0, 2.0), (0.5, ~0.5)]
        in
          Harness.check "fft (ifft X) = X (n = 6)"
            (approxArrC (F.fft (F.ifft xc), xc))
        end

      (* ---- linearity: fft(a*x + y) = a*fft(x) + fft(y) ---- *)
      val () = Harness.section "Linearity"
      val () =
        let
          val n = 8
          val a = 2.5
          val x = Array.tabulate (n, fn j => C.complex (Math.sin (0.6 * real j),
                                                        Math.cos (0.2 * real j)))
          val y = Array.tabulate (n, fn j => C.complex (real j - 3.0,
                                                        Math.sin (1.1 * real j)))
          val combined =
            Array.tabulate (n, fn j =>
              C.add (C.scale (a, Array.sub (x, j)), Array.sub (y, j)))
          val lhs = F.fft combined
          val fx = F.fft x
          val fy = F.fft y
          val rhs =
            Array.tabulate (n, fn k =>
              C.add (C.scale (a, Array.sub (fx, k)), Array.sub (fy, k)))
        in
          Harness.check "fft(a*x + y) = a*fft(x) + fft(y)" (approxArrC (lhs, rhs))
        end

      (* ---- known transform: constant signal -> delta at bin 0 ---- *)
      val () = Harness.section "Known transforms"
      val () =
        let
          val n = 8
          val ones = Array.tabulate (n, fn _ => C.complex (1.0, 0.0))
          val X = F.fft ones
          val bin0 = approxC (Array.sub (X, 0), C.complex (real n, 0.0))
          fun restZero () =
            let
              fun loop k =
                k >= n orelse
                (approxC (Array.sub (X, k), C.complex (0.0, 0.0)) andalso loop (k + 1))
            in
              loop 1
            end
        in
          Harness.check "constant -> delta at bin 0" (bin0 andalso restZero ())
        end

      val () =
        let
          (* A pure cosine cos(2*pi*k0*j/n) has peaks of magnitude n/2 at
             bins k0 and n-k0, and ~0 elsewhere. *)
          val n = 16
          val k0 = 3
          val cosine =
            Array.tabulate
              (n, fn j => C.complex (Math.cos (2.0 * pi * real k0 * real j / real n), 0.0))
          val X = F.fft cosine
          val peakLo = approx (C.abs (Array.sub (X, k0)), real n / 2.0)
          val peakHi = approx (C.abs (Array.sub (X, n - k0)), real n / 2.0)
          fun othersZero () =
            let
              fun loop k =
                k >= n orelse
                (if k = k0 orelse k = n - k0 then loop (k + 1)
                 else approx (C.abs (Array.sub (X, k)), 0.0) andalso loop (k + 1))
            in
              loop 0
            end
        in
          Harness.check "single cosine -> peaks at +/-k0"
            (peakLo andalso peakHi andalso othersZero ())
        end

      (* ---- Parseval's theorem: sum|x|^2 = (1/n) sum|X|^2 ---- *)
      val () = Harness.section "Parseval's theorem"
      val () =
        List.app
          (fn n =>
            let
              val x = signal n
              val xc = Array.tabulate (n, fn j => C.complex (Array.sub (x, j), 0.0))
              val X = F.fft xc
            in
              Harness.check
                ("n = " ^ Int.toString n)
                (approx (energyC xc, energyC X / real n))
            end)
          [4, 6, 8, 12]

      val () =
        Harness.check "Parseval (real energy form, n = 8)"
          (let
             val x = signal 8
             val xc = Array.tabulate (8, fn j => C.complex (Array.sub (x, j), 0.0))
           in
             approx (energyR x, energyC (F.fft xc) / 8.0)
           end)

      (* ---- rfft matches fft of the real-embedded signal ---- *)
      val () = Harness.section "rfft (real input)"
      val () =
        let
          val xs = [1.0, 2.0, 3.0, 4.0, 5.0, 0.0, ~1.0]   (* length 7 *)
          val x = rarr xs
          val xc = rcarr xs
        in
          Harness.check "rfft x = fft (embed x)" (approxArrC (F.rfft x, F.fft xc))
        end

      (* ---- convolution theorem: convolve = naive linear convolution ---- *)
      val () = Harness.section "Convolution theorem"
      val () =
        let
          val a = rarr [1.0, 2.0, 3.0]
          val b = rarr [0.0, 1.0, 0.5, ~2.0]
        in
          Harness.check "convolve(a,b) = naive (3 x 4)"
            (approxArrR (F.convolve (a, b), naiveConv (a, b)))
        end
      val () =
        let
          val a = rarr [2.0, ~1.0, 0.0, 4.0, 3.0]
          val b = rarr [1.0, 1.0, 1.0]   (* boxcar / moving sum *)
        in
          Harness.check "convolve(a,b) = naive (5 x 3)"
            (approxArrR (F.convolve (a, b), naiveConv (a, b)))
        end
      val () =
        let
          (* Equal, power-of-two-ish lengths landing on a non-trivial pad. *)
          val a = signal 6
          val b = signal 6
        in
          Harness.check "convolve(a,b) = naive (6 x 6)"
            (approxArrR (F.convolve (a, b), naiveConv (a, b)))
        end
      val () =
        let
          val a = rarr [1.0, ~2.0, 3.0, 4.0, ~1.0, 2.0, 0.5, ~0.5]
          val b = rarr [0.5, 0.25]
        in
          Harness.check "convolve(a,b) = naive (8 x 2)"
            (approxArrR (F.convolve (a, b), naiveConv (a, b)))
        end

      (* ---- irfft recovers a real signal ---- *)
      val () = Harness.section "irfft (inverse real FFT)"
      val () =
        List.app
          (fn n =>
            let
              val x = signal n
            in
              Harness.check
                ("irfft (rfft x) = x, n = " ^ Int.toString n)
                (approxArrR (F.irfft (F.rfft x), x))
            end)
          [1, 2, 3, 4, 5, 6, 8, 9, 16]

      (* ---- 2D transforms: fft2 == naive separable DFT, and round-trip ---- *)
      val () = Harness.section "fft2 / ifft2 (separable 2D)"
      val () =
        List.app
          (fn (rows, cols) =>
            let
              val m =
                Array.tabulate (rows, fn r =>
                  Array.tabulate (cols, fn c =>
                    C.complex (0.5 + Math.sin (0.7 * real r + 0.3 * real c),
                               Math.cos (0.2 * real r - 0.5 * real c))))
            in
              Harness.check
                ("fft2 = naive 2D DFT, " ^ Int.toString rows ^ "x"
                 ^ Int.toString cols)
                (approxArrArrC (F.fft2 m, naiveDft2 m));
              Harness.check
                ("ifft2 (fft2 m) = m, " ^ Int.toString rows ^ "x"
                 ^ Int.toString cols)
                (approxArrArrC (F.ifft2 (F.fft2 m), m))
            end)
          [(1, 1), (2, 3), (3, 2), (4, 4), (2, 5), (6, 3)]
      val () =
        Harness.check "fft2 of empty matrix is empty"
          (Array.length (F.fft2 (Array.fromList [])) = 0)

      (* ---- DCT-II / DCT-III ---- *)
      val () = Harness.section "dct / idct"
      val () =
        List.app
          (fn n =>
            let
              val x = signal n
            in
              Harness.check
                ("dct = naive DCT-II, n = " ^ Int.toString n)
                (approxArrR (F.dct x, naiveDct x));
              Harness.check
                ("idct (dct x) = x, n = " ^ Int.toString n)
                (approxArrR (F.idct (F.dct x), x))
            end)
          [1, 2, 3, 4, 5, 6, 8, 12]
      val () =
        Harness.check "dct [1,1,1,1] = [4,0,0,0]"
          (approxArrR (F.dct (rarr [1.0, 1.0, 1.0, 1.0]),
                       rarr [4.0, 0.0, 0.0, 0.0]))

      (* ---- edge cases ---- *)
      val () = Harness.section "Edge cases"
      val () =
        Harness.check "fft of empty array is empty"
          (Array.length (F.fft (Array.fromList [])) = 0)
      val () =
        Harness.check "fft of singleton is identity"
          (approxArrC (F.fft (carr [(7.0, ~3.0)]), carr [(7.0, ~3.0)]))
    in
      ()
    end

  fun run () = (Harness.reset (); runAll (); Harness.run ())
end
