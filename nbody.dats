#include "share/atspre_define.hats"
#include "share/atspre_staload.hats"

staload "libc/SATS/math.sats"

typedef pos = @{ x = double, y = double, z = double }

typedef vel = @{ x = double, y = double, z = double }

typedef body = @{ pos = pos, vel = vel, mass = double }

fun compute_offset (system : arrszref(body)) : (double, double, double) = let
  val n = system.size()
  fun loop (off : (double, double, double), i : size_t)
    :<cloref1> (double, double, double) =
      if i < n
        then let
          val v = system[i].vel
          val m = system[i].mass
          val off1 = ( off.0 + v.x * m
                     , off.1 + v.y * m
                     , off.2 + v.z * m
                     )
          in loop (off1, i + i2sz(1))
        end
        else off
  in loop ((0.0, 0.0, 0.0), i2sz(0))
end

(* Offset Sol's momentum so that the system has no net momentum. *)
fun apply_offset (sol : body, off : (double, double, double)) : body =
  @{ pos = sol.pos
   , vel = @{ x = sol.vel.x - off.0 / sol.mass
            , y = sol.vel.y - off.1 / sol.mass
            , z = sol.vel.z - off.2 / sol.mass
            }
   , mass = sol.mass
   }

val pi = 3.141592653589793 : double

val sol =
  @{ pos = @{ x = 0.0, y = 0.0, z = 0.0 } : pos
   , vel = @{ x = 0.0, y = 0.0, z = 0.0 } : vel
   , mass = 4 * pi ** 2
   } : body

val dp = 365.24 : double

val jupiter =
  @{ pos = @{ x = 4.84143144246472090e+00
            , y = ~1.16032004402742839e+00
            , z= ~1.03622044471123109e-01
            } : pos
   , vel = @{ x = 1.66007664274403694e-03 * dp
            , y = 7.69901118419740425e-03 * dp
            , z = ~6.90460016972063023e-05 * dp
            } : vel
   , mass = 9.54791938424326609e-04 * sol.mass
   } : body

val saturn =
  @{ pos = @{ x = 8.34336671824457987e+00
            , y = 4.12479856412430479e+00
            , z = ~4.03523417114321381e-01
            } : pos
   , vel = @{ x = ~2.76742510726862411e-03 * dp
            , y = 4.99852801234917238e-03 * dp
            , z = 2.30417297573763929e-05 * dp
            } : vel
   , mass = 2.85885980666130812e-04 * sol.mass
   } : body

val uranus =
  @{ pos = @{ x = 1.28943695621391310e+01
            , y = ~1.51111514016986312e+01
            , z = ~2.23307578892655734e-01
            } : pos
   , vel = @{ x = 2.96460137564761618e-03 * dp
            , y = 2.37847173959480950e-03 * dp
            , z = ~2.96589568540237556e-05 * dp
            } : vel
   , mass = 4.36624404335156298e-05 * sol.mass
   } : body

val neptune =
  @{ pos = @{ x = 1.53796971148509165e+01
            , y = ~2.59193146099879641e+01
            , z = 1.79258772950371181e-01
            } : pos
   , vel = @{ x = 2.68067772490389322e-03 * dp
            , y = 1.62824170038242295e-03 * dp
            , z = ~9.51592254519715870e-05 * dp
            } : vel
   , mass = 5.15138902046611451e-05 * sol.mass
   } : body

fun init() : arrszref(body) = let
  val system = arrszref_make_elt<body>(i2sz(5), sol)
  val () = system[1] := jupiter
  val () = system[2] := saturn
  val () = system[3] := uranus
  val () = system[4] := neptune
  val off = compute_offset(system)
  val () = system[0] := apply_offset(sol, off)
  in system
end

fun squared(vec : @{ x = double, y = double, z = double }) = let
  val x = vec.x
  val y = vec.y
  val z = vec.z
  in x * x + y * y + z * z
end

fun kinetic(system : arrszref(body)) : double = let
  val n = system.size()
  fun loop(T : double, i : size_t) : double =
    if i < n
      then let
        val dT = 0.5 * system[i].mass * squared(system[i].vel)
        in loop(T + dT, i + i2sz(1))
      end
      else T
  in loop(0.0, i2sz(0))
end

fun potential(system : arrszref(body)) : double = let
  val n = system.size()
  fun loop(V : double, i : size_t, j : size_t) : double =
    if j < n
      then let
        val ri = system[i].pos
        val rj = system[j].pos
        val rij = @{ x = ri.x - rj.x
                   , y = ri.y - rj.y
                   , z = ri.z - rj.z
                   }
        val dV = (~1.0) * system[i].mass * system[j].mass / sqrt_double(squared(rij))
        in loop(V + dV, i, j + i2sz(1))
      end
      else if i < n
        then loop(V, i + i2sz(1), i + i2sz(2))
        else V
  in loop(0.0, i2sz(0), i2sz(1))
end

fun energy(system : arrszref(body)) : double =
  kinetic(system) + potential(system)

val dt = 0.01 : double

fun kinematic(system : arrszref(body)) : void = let
  val n = system.size()
  fun loop(i : size_t) :<cloref1> void =
    if i < n
      then let
        val r = system[i].pos
        val v = system[i].vel
        in (
          system[i] := @{ pos = @{ x = r.x + v.x * dt
                                 , y = r.y + v.y * dt
                                 , z = r.z + v.z * dt
                                 }
                        , vel = v
                        , mass = system[i].mass
                        };
          loop(i + i2sz(1))
        )
      end
      else ()
  in loop(i2sz(0))
end

fun dynamic(system : arrszref(body)) : void = let
  val n = system.size()
  fun loop(i : size_t, j : size_t) :<cloref1> void =
    if j < n
      then let
        val ri = system[i].pos
        val rj = system[j].pos
        val rij = @{ x = ri.x - rj.x
                   , y = ri.y - rj.y
                   , z = ri.z - rj.z
                   }
        val rr = squared(rij)
        val r = sqrt_double(rr)
        val mag = dt / (rr * r)
        val vi = system[i].vel
        val vj = system[j].vel
        val mi = system[i].mass
        val mj = system[j].mass
        val () = system[i] := @{ pos = ri
                               , vel = @{ x = vi.x - rij.x * mj * mag
                                        , y = vi.y - rij.y * mj * mag
                                        , z = vi.z - rij.z * mj * mag
                                        }
                               , mass = mi
                               }
        val () = system[j] := @{ pos = rj
                               , vel = @{ x = vj.x + rij.x * mi * mag
                                        , y = vj.y + rij.y * mi * mag
                                        , z = vj.z + rij.z * mi * mag
                                        }
                               , mass = mj
                               }
        in loop(i, j + i2sz(1))
      end
      else if i < n
        then loop(i + i2sz(1), i + i2sz(2))
        else ()
  in loop(i2sz(0), i2sz(1))
end

fun advance(system : arrszref(body)) : void =
  ( dynamic(system) ; kinematic(system) )

fun run(system : arrszref(body), i : int) : void =
  if i > 0
    then ( advance(system) ; run(system, i - 1) )
    else ()

implement main0 () = let
  val system = init()
  in (
    $extfcall(void, "printf", "%.9f\n", energy(system));
    run(system, 50000000);
    $extfcall(void, "printf", "%.9f\n", energy(system));
  )
end