#pragma once
#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <math.h>
#if defined(__aarch64__) || defined(_M_ARM64)
    #include <arm_neon.h>
#elif defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    #include <xmmintrin.h>
#endif

inline double drand()
{
	return (double)rand()/(double)RAND_MAX;
}

inline void vrand(double *v)
{
	v[0] = 2.0*drand()-1.0;
	v[1] = 2.0*drand()-1.0;
	v[2] = 2.0*drand()-1.0;
}

inline void vset(double x, double y, double z, double *v)
{
	v[0] = x;
	v[1] = y;
	v[2] = z;
}

inline void vcopy(double *v, double *vr)
{
	vr[0] = v[0];
	vr[1] = v[1];
	vr[2] = v[2];
}

inline double vmag(double *v)
{
	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

inline double vmagsq(double *v)
{
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

inline double vmagsoft(double *v, double r_soft)
{
	return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + r_soft);
}

inline double vmagsqsoft(double *v, double r_soft)
{
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + r_soft;
}

inline void vscale(double *v, double s, double *vr)
{
	vr[0] = v[0]*s;
	vr[1] = v[1]*s;
	vr[2] = v[2]*s;
}

inline void vnorm(double *v)
{
	vscale(v,1.0/vmag(v),v);
}

inline void vneg(double *v)
{
	v[0] = -v[0];
	v[1] = -v[1];
	v[2] = -v[2];
}

inline void vadd(double *v1, double *v2, double *vr)
{
	vr[0] = v1[0] + v2[0];
	vr[1] = v1[1] + v2[1];
	vr[2] = v1[2] + v2[2];
}

inline void vscaleadd(double *v, double s, double *vr)
{
	vr[0] = v[0]*s+vr[0];
	vr[1] = v[1]*s+vr[1];
	vr[2] = v[2]*s+vr[2];
}

inline void vsub(double *v1, double *v2, double *vr)
{
	vr[0] = v1[0] - v2[0];
	vr[1] = v1[1] - v2[1];
	vr[2] = v1[2] - v2[2];
}

inline double vdot(double *v1, double *v2)
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

inline void vcross(double *v1, double *v2, double *vr)
{
	double v[3];
	v[0] = v1[1]*v2[2] - v1[2]*v2[1];
	v[1] = v1[2]*v2[0] - v1[0]*v2[2];
	v[2] = v1[0]*v2[1] - v1[1]*v2[0];
	vcopy(v,vr);
}

inline void vmean(double *v1, double *v2, double *vr)
{
	vadd(v1,v2,vr);
	vscale(vr,1.0/2.0,vr);
}

inline bool vequal(double *v1, double *v2)
{
	return ((v1[0]==v2[0]) && (v1[1]==v2[1]) &&  (v1[2]==v2[2]));
}

inline void vmin(double *v1, double *v2, double *vr)
{
	vr[0] = (v1[0]<v2[0])?v1[0]:v2[0];
	vr[1] = (v1[1]<v2[1])?v1[1]:v2[1];
	vr[2] = (v1[2]<v2[2])?v1[2]:v2[2];
}

inline void vmax(double *v1, double *v2, double *vr)
{
	vr[0] = (v1[0]>v2[0])?v1[0]:v2[0];
	vr[1] = (v1[1]>v2[1])?v1[1]:v2[1];
	vr[2] = (v1[2]>v2[2])?v1[2]:v2[2];
}

inline double vmaxel(double *v)
{
	double m;
	m = (v[0]>v[1])?v[0]:v[1];
	return (m>v[2])?m:v[2];
}

inline double vmaxsub(double *v1, double *v2)
{
	double sub[3];
	vsub(v1,v2,sub);
	return vmaxel(sub);
}

inline void vmaxboundcube(double *p_min_in, double *p_max_in, double *p_min_out, double *p_max_out)
{
	double halfwidth;
	double mean[3],delta[3];
	halfwidth = vmaxsub(p_max_in,p_min_in)/2.0;
	vmean(p_max_in,p_min_in,mean);
	vset(halfwidth,halfwidth,halfwidth,delta);
	vcopy(mean,p_min_out);
	vsub(p_min_out,delta,p_min_out);
	vcopy(mean,p_max_out);
	vadd(p_max_out,delta,p_max_out);
}

inline void vrot_zy(double *v, double angle, double *vr)
{
	double ca = cos(angle);
	double sa = sin(angle);
	double vt[3];

	vt[2] = v[2]*ca - v[1]*sa;
	vt[1] = v[1]*ca + v[2]*sa;
	vt[0] = v[0];

	vcopy(vt,vr);
}

// Float versions for BHTree

inline void vsetf(float x, float y, float z, float *v)
{
	v[0] = x; v[1] = y; v[2] = z;
}

inline void vcopyf(float *v, float *vr)
{
	vr[0] = v[0]; vr[1] = v[1]; vr[2] = v[2];
}

inline void vaddf(float *v1, float *v2, float *vr)
{
	vr[0] = v1[0]+v2[0]; vr[1] = v1[1]+v2[1]; vr[2] = v1[2]+v2[2];
}

inline void vsubf(float *v1, float *v2, float *vr)
{
	vr[0] = v1[0]-v2[0]; vr[1] = v1[1]-v2[1]; vr[2] = v1[2]-v2[2];
}

inline void vscalef(float *v, float s, float *vr)
{
	vr[0] = v[0]*s; vr[1] = v[1]*s; vr[2] = v[2]*s;
}

inline float vmagsqf(float *v)
{
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

inline float vmagsqsoftf(float *v, float r_soft)
{
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + r_soft;
}

inline void vmeanf(float *v1, float *v2, float *vr)
{
	vr[0] = (v1[0]+v2[0])*0.5f;
	vr[1] = (v1[1]+v2[1])*0.5f;
	vr[2] = (v1[2]+v2[2])*0.5f;
}

inline bool vequalf(float *v1, float *v2)
{
	return (v1[0]==v2[0]) && (v1[1]==v2[1]) && (v1[2]==v2[2]);
}

inline void vscaleaddf(float *v, float s, double *vr)
{
	vr[0] += v[0] * s;
	vr[1] += v[1] * s;
	vr[2] += v[2] * s;
}

inline void vminf(float *v1, float *v2, float *vr)
{
	vr[0] = (v1[0]<v2[0])?v1[0]:v2[0];
	vr[1] = (v1[1]<v2[1])?v1[1]:v2[1];
	vr[2] = (v1[2]<v2[2])?v1[2]:v2[2];
}

inline void vmaxf(float *v1, float *v2, float *vr)
{
	vr[0] = (v1[0]>v2[0])?v1[0]:v2[0];
	vr[1] = (v1[1]>v2[1])?v1[1]:v2[1];
	vr[2] = (v1[2]>v2[2])?v1[2]:v2[2];
}

inline float fast_rsqrtf(float x)
{
#if defined(__aarch64__) || defined(_M_ARM64)
    // -------------------------------------------------------------------------
    // ARM64 NEON Path (Jetson Orin Nano / Apple Silicon / ARM Servers)
    // -------------------------------------------------------------------------
    // 1. Get initial 8-bit reciprocal square root estimate
    float32x4_t val = vdupq_n_f32(x);
    float32x4_t est = vrsqrteq_f32(val);

    // 2. Exact Newton-Raphson step via native NEON instruction (vrsqrtsq_f32)
    //    vrsqrtsq_f32(a, b) computes: (3 - a * b) / 2
    float32x4_t step = vrsqrtsq_f32(val, vmulq_f32(est, est));
    float32x4_t refined = vmulq_f32(est, step);

    return vgetq_lane_f32(refined, 0);

#elif defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    // -------------------------------------------------------------------------
    // x86 / x64 SSE Path (MSVC / GCC / Clang on Intel & AMD)
    // -------------------------------------------------------------------------
    __m128 val = _mm_set_ss(x);
    __m128 est = _mm_rsqrt_ss(val);
    __m128 half = _mm_set_ss(0.5f);
    __m128 three = _mm_set_ss(3.0f);

    // One Newton-Raphson iteration: est = est * (3 - x * est^2) * 0.5
    __m128 est2 = _mm_mul_ss(est, est);
    __m128 xe2 = _mm_mul_ss(val, est2);
    __m128 diff = _mm_sub_ss(three, xe2);
    __m128 refined = _mm_mul_ss(_mm_mul_ss(est, diff), half);

    return _mm_cvtss_f32(refined);

#else
    // -------------------------------------------------------------------------
    // Portable Scalar Fallback (Any other architecture)
    // -------------------------------------------------------------------------
    return 1.0f / std::sqrt(x);
#endif
}

inline float fast_r3_inv(float dsq)
{
	float rsqrt = fast_rsqrtf(dsq);
	return rsqrt * rsqrt * rsqrt;
}
