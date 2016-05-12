#include "share/atspre_define.hats"
#include "share/atspre_staload.hats"

staload "libc/SATS/stdlib.sats"
staload "libc/SATS/math.sats"

typedef vec3 = @{ x = double, y = double, z = double }

fun vec3_add(a : vec3, b : vec3) : vec3 =
  @{ x = a.x + b.x
   , y = a.y + b.y
   , z = a.z + b.z
   }
overload + with vec3_add

fun vec3_sub(a : vec3, b : vec3) : vec3 =
  @{ x = a.x - b.x
   , y = a.y - b.y
   , z = a.z - b.z
   }
overload - with vec3_sub

fun vec3_scale(s : double, a : vec3) : vec3 =
  @{ x = s * a.x
   , y = s * a.y
   , z = s * a.z
   }
overload * with vec3_scale

typedef system = @{ pos = arrszref(vec3)
                  , vel = arrszref(vec3)
                  , mass = arrszref(double)
                  , n = int
                  }

val dp = 365.24 : double
val dt = 0.01 : double

val sol =
  @{ pos = @{ x = 0.0, y = 0.0, z = 0.0 } : vec3
   , vel = @{ x = 0.0, y = 0.0, z = 0.0 } : vec3
   , mass = 4 * M_PI ** 2
   }

val jupiter =
  @{ pos = @{ x = 4.84143144246472090e+00
            , y = ~1.16032004402742839e+00
            , z= ~1.03622044471123109e-01
            } : vec3
   , vel = @{ x = 1.66007664274403694e-03 * dp
            , y = 7.69901118419740425e-03 * dp
            , z = ~6.90460016972063023e-05 * dp
            } : vec3
   , mass = 9.54791938424326609e-04 * sol.mass
   }

val saturn =
  @{ pos = @{ x = 8.34336671824457987e+00
            , y = 4.12479856412430479e+00
            , z = ~4.03523417114321381e-01
            } : vec3
   , vel = @{ x = ~2.76742510726862411e-03 * dp
            , y = 4.99852801234917238e-03 * dp
            , z = 2.30417297573763929e-05 * dp
            } : vec3
   , mass = 2.85885980666130812e-04 * sol.mass
   }

val uranus =
  @{ pos = @{ x = 1.28943695621391310e+01
            , y = ~1.51111514016986312e+01
            , z = ~2.23307578892655734e-01
            } : vec3
   , vel = @{ x = 2.96460137564761618e-03 * dp
            , y = 2.37847173959480950e-03 * dp
            , z = ~2.96589568540237556e-05 * dp
            } : vec3
   , mass = 4.36624404335156298e-05 * sol.mass
   }

val neptune =
  @{ pos = @{ x = 1.53796971148509165e+01
            , y = ~2.59193146099879641e+01
            , z = 1.79258772950371181e-01
            } : vec3
   , vel = @{ x = 2.68067772490389322e-03 * dp
            , y = 1.62824170038242295e-03 * dp
            , z = ~9.51592254519715870e-05 * dp
            } : vec3
   , mass = 5.15138902046611451e-05 * sol.mass
   }

fun offset (@{ vel = v, mass = m, n = n, ... } : system) : vec3 = let
  fun loop (i : int, off : vec3) :<cloref1> vec3 =
    if i < n
      then loop (succ(i), off + m[i] * v[i])
      else off
  in (~1.0 / sol.mass) * loop(0, @{ x = 0.0, y = 0.0, z = 0.0 })
  end

fun init () : system = let
  val ms = (arrszref)$arrpsz(sol.mass, jupiter.mass, saturn.mass, uranus.mass, neptune.mass)
  val vs = (arrszref)$arrpsz(sol.vel, jupiter.vel, saturn.vel, uranus.vel, neptune.vel)
  val xs = (arrszref)$arrpsz(sol.pos, jupiter.pos, saturn.pos, uranus.pos, neptune.pos)
  val sys = @{ pos = xs, vel = vs, mass = ms, n = 5 }
  val () = vs[0] := offset(sys)
  in @{ pos = xs, vel = vs, mass = ms, n = 5 }
  end

fun squared (@{ x = x, y = y, z = z } : vec3) = x * x + y * y + z * z

fun kinetic(@{ mass = m, vel = v, n = n, ... } : system) : double = let
  fun loop(i : size_t, T : double) :<cloref1> double =
    if i < n
      then let
        val dT = 0.5 * m[i] * squared(v[i])
        in loop(succ(i), T + dT)
      end
      else T
  in loop(i2sz(0), 0.0)
end

fun potential
  (@{ mass = m, vel = v, pos = x, n = n, ... } : system) : double = let
  fun loop(V : double, i : size_t, j : size_t) :<cloref1> double =
    if j < n
      then let
        val r = x[i] - x[j]
        val dV = (~1.0) * m[i] * m[j] / sqrt_double(squared(r))
        in loop(V + dV, i, succ(j))
        end
      else if i < n
        then loop(V, succ(i), succ(succ(i)))
        else V
  in loop(0.0, i2sz(0), i2sz(1))
  end

fun energy(sys : system) : double =
  kinetic(sys) + potential(sys)

fun kinematic(@{ vel = v, pos = x, n = n, ... } : system) : void = let
  fun loop(i : size_t) :<cloref1> void =
    if i < n
      then ( x[i] := x[i] + dt * v[i]
           ; loop(succ i) )
      else ()
  in loop(i2sz(0))
  end

fun dynamic(@{ mass = m, pos = x, vel = v, n = n, ... } : system) : void = let
  fun loop(i : size_t, j : size_t) :<cloref1> void =
    if j < n
      then let
        val rij = x[i] - x[j]
        val rr = squared(rij)
        val r = sqrt_double(rr)
        val mag = dt / (r * rr)
        val () = v[i] := v[i] - mag * m[j] * rij
        val () = v[j] := v[j] + mag * m[i] * rij
        in loop(i, succ(j))
        end
      else if i < n
        then loop(succ(i), succ(succ(i)))
        else ()
  in loop(i2sz(0), i2sz(1))
  end

fun advance(sys : system) : void = ( dynamic(sys) ; kinematic(sys) )

fun run(i : int, sys : system) : void =
  if i > 0
    then ( advance(sys) ; run(i - 1, sys) )
    else ()

implement main (argc, argv) = let
  val sys = init()
  in if (1 < argc)
    then (
      $extfcall(void, "printf", "%.9f\n", energy(sys));
      run(atoi(argv[1]), sys);
      $extfcall(void, "printf", "%.9f\n", energy(sys));
      0
    )
    else 1
end