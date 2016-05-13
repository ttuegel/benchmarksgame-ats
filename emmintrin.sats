%{#
#include "emmintrin.h"
%}

typedef m128d = $extype"__m128d"

fun _mm_add_pd (m128d, m128d) : m128d = "mac#"
overload + with _mm_add_pd

fun _mm_sub_pd (m128d, m128d) : m128d = "mac#"
overload - with _mm_sub_pd

fun _mm_sqrt_pd (m128d) : m128d = "mac#"

fun _mm_load_sd (ptr) : m128d = "mac#"
fun _mm_loadh_pd (m128d, ptr) : m128d = "mac#"
fun m128d_load (double, double) : m128d

fun _mm_storel_pd (ptr, m128d) : void = "mac#"
fun _mm_storeh_pd (ptr, m128d) : void = "mac#"
fun m128d_store (m128d) : @(double, double)