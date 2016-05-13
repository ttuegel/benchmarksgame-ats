staload "./emmintrin.sats"

implement m128d_load (a : double, b : double) : m128d = let
  var a = a
  var b = b
  in _mm_loadh_pd(_mm_load_sd(addr@a), addr@b)
  end

implement m128d_store (m : m128d) : @(double, double) = let
  var a = 0.0
  var b = 0.0
  val () = _mm_storel_pd (addr@a, m)
  val () = _mm_storeh_pd (addr@b, m)
  in @(a, b)
  end