#pragma once

#include <iostream>
#include <string>
#include <SDL3/SDL.h>
#include "glad.h"
#include <vector>
#include <thread>
#include <random>
#include <algorithm>
#include "VectorMath.h"
#include "BHTree.h"
#include "ThreadPool.h"

struct Camera
{
	double pos[3];
	double lookAt[3];
	double phi,theta;
};

class Simulation
{
	int N_Bodies;
	int N_Systems;
	std::vector<int> N_System_Bodies;

	double G;
	double FDE;
	double dt;
	double t;
	double r_soft;
	double BH_Opening_Theta;
	double DisplayScale;
	double accel_sq_color_thresh;
	// Exposure for the HDR tone-map post pass (out = 1 - exp(-exposure * L)).
	// Higher lifts faint/lone particles; cores saturate gracefully toward white.
	double ToneMapExposure;

	bool Gravity_P2P;
	bool Gravity_Oct;
	bool Record_Video;
	bool Data_Log;
	bool Info_Display;

	bool CamOrbit;
	double CamOrbitTheta;

	// Camera_lookAt_System: when >= 0 the camera's look-at point is retargeted
	// every frame onto that system's central body, so a moving galaxy stays
	// centred in frame. Only the look-at point moves; the camera's spherical
	// offset (phi, theta, r) is preserved, so the viewing angle and zoom are
	// unchanged and remain under user control. -1 disables the feature, leaving
	// the fixed Camera_lookAt point in effect.
	int CamFollowSystem;
	double EndTime;

	int DisplayWidth, DisplayHeight;

	std::vector<double> halo_vc;
	std::vector<double> halo_rc_sq;
	// Halo truncation radius. <= 0 means untruncated (enclosed mass grows without
	// bound as M_halo(r) ~ Vc^2 * r). When set, the enclosed mass is frozen at
	// M_halo(Rh) beyond Rh and the force falls off as 1/r^2, which is what
	// Salo & Laurikainen (2000) sect 2.2 do (Rh = 400 arcsec for the primary,
	// equal to its disc truncation Rd).
	std::vector<double> halo_rh;
	std::vector<double> halo_M_rh;   // precomputed Vc^2 * Rh^3/(Rh^2+Rc^2)
	std::vector<int> halo_central;
	std::vector<int> body_system;
	std::vector<double> halo_center;

	// Salo & Laurikainen (2000) treat each rigid halo as an INERTIAL body whose
	// centre is integrated under gravity ("the disc back-action is taken into
	// account in the halo motion"), rather than re-pinned to the live particle
	// barycentre each step. These hold that centre's velocity, current/previous
	// acceleration (for velocity-Verlet), and its inertial mass.
	std::vector<double> halo_vel;
	std::vector<double> halo_acc;
	std::vector<double> halo_acc_prev;
	std::vector<double> halo_mass;   // inertial mass of the halo (0 = no halo)

	std::vector<double> mass;

	std::vector<double> pos_data;
	std::vector<double> vel_data;
	std::vector<double> acc_data;
	std::vector<double> acc_prev_data;

	std::vector<double*> pos;
	std::vector<double*> vel;
	std::vector<double*> acc;
	std::vector<double*> acc_prev;

	std::vector<bool> has_gravity;

	// Dissipative "sticky particle" gas (Salo & Laurikainen 2000 sect 2.1).
	// Gas particles gravitate exactly like stars (they live in the Octree via the
	// arrays above); is_gas[i] additionally marks them for inelastic gas-gas
	// collisions handled by ProcessGasCollisions(). char, not bool, so the
	// collision loop can index it without the vector<bool> proxy.
	std::vector<char> is_gas;

	std::vector<double> pos_sq;
	std::vector<double> vel_sq;
	std::vector<double> acc_sq;

	std::vector<float> pos_f;
	std::vector<int> sortedIdx;
	std::vector<int> sortTemp;
	std::vector<uint32_t> mortonCodes;
	int numActiveBodies;

	double totalKE, totalPE, totalE;
	double totalP, totalL, comSpeed;
	double virialRatio;
	std::vector<double> body_pot;

	bool Remove_Halo_Monopole;       // cancel net force of the rigid analytic halos

	// Gas collision parameters (see is_gas above). On each gas-gas impact the
	// component of the relative velocity along the line of centres becomes
	// -gas_alpha times its prior value (gas_alpha = 0 is the paper's fully
	// inelastic value); two gas particles collide when within gas_radius sum.
	// Gas collision parameters (see is_gas above). Collisions are handled by a
	// kinetic Monte-Carlo (DSMC-style) step: within each collision cell, pairs are
	// selected stochastically at the physical rate n*sigma*v_rel, and on each
	// accepted impact the component of the relative velocity along a STOCHASTICALLY
	// sampled line of centres becomes -gas_alpha times its prior value (gas_alpha=0
	// is the paper's fully inelastic value). Sampling the impact geometry (rather
	// than reading it off instantaneous positions) makes the cooling isotropic, so
	// in-plane sigma_r cools as well as sigma_z -- unlike a geometric-overlap test,
	// which in a thin disc is biased toward vertical impacts. sigma = pi*(2*radius)^2.
	double gas_alpha;                // restitution coefficient (default 0)
	double gas_radius;               // per-gas-particle collision radius, code units
	double gas_cell_size;            // DSMC collision-cell edge, code units (default 6)
	double gas_softening;            // gravitational softening LENGTH for gas sinks,
	                                 // code units; 0 = use r_soft (star softening).
	                                 // Larger than r_soft suppresses gas self-gravity
	                                 // fragmentation (balls) while stars keep sharp
	                                 // (small-softening) arms. See SoftSq().
	bool gasEnabled;                 // true iff any body is tagged gas
	std::mt19937 gasRng;             // RNG for the Monte-Carlo collision sampling

	// Warmup / initialization phase.
	// When InitializationTime > 0 the simulation starts at t = -InitializationTime
	// with the systems dynamically ISOLATED: gravity acts only within a system,
	// cross-system halo terms are off, and each system's bulk velocity is
	// withheld. This lets each galaxy relax out of its initial particle-noise
	// transients before the interaction begins. At t >= 0 the stored bulk
	// velocities are applied and full N-body coupling resumes.
	// Zero (the default) disables all of this: the run starts at t = 0 exactly
	// as before.
	double InitializationTime;
	bool warmupActive;               // true while t < 0 and warmup is enabled
	bool bulkVelocityApplied;        // one-shot guard for the t=0 transition
	std::vector<double> system_bulk_vel;   // 3 per system, withheld during warmup

	Camera Cam;

	// Particle rendering (modern GL, instanced)
	GLuint particleVAO, particleShapeVBO, particlePosVBO, particleColorVBO;
	GLuint particleShader;
	float *posBuf;
	float *clrBuf;

	// Octree wireframe rendering
	GLuint octreeVAO, octreeVBO;
	GLuint octreeShader;
	std::vector<float> octreeVerts;

	// FPS overlay
	GLuint fpsVAO, fpsVBO;
	GLuint fpsShader;

	// Off-screen FBO for recording
	GLuint recordFBO, recordColorTex, recordDepthRBO;

	// HDR accumulation buffer (RGBA16F): particles are additively blended here so
	// core densities can exceed 1.0 without clipping, then a tone-map fullscreen
	// pass compresses it into the LDR final target (screen or record FBO).
	GLuint hdrFBO, hdrColorTex, hdrDepthRBO;
	GLuint tonemapShader, tonemapVAO;

	int winWidth, winHeight;

	BHTree Octree;

	FILE *DataLog;

	// TEMPORARY m=2 bar/arm diagnostic: periodically log the m=2 Fourier amplitude
	// A2(R) and phase of each disc's surface density, so bar (strong, radially
	// coherent A2 with near-constant phase) and spiral (A2 with winding phase) can
	// be measured. Enabled by the BarDiagnostic script command; writes
	// bar_diagnostic.csv. Slated for removal once the generalized diagnostic system
	// lands.
	FILE *barLog = nullptr;
	int barLogEvery = 0;        // steps between rows; 0 = disabled
	long barStepCount = 0;

	ThreadPool *pool;

	// Scale factor s such that the halo acceleration on a body is
	//   a_vec = s * (halo_centre - position)
	// Continuous at r = Rh by construction: the truncated branch evaluates to
	// Vc^2/(Rh^2+Rc^2) there, matching the cored-isothermal branch exactly.
	inline double HaloScale(int sys, double rsq) const {
		double rh = halo_rh[sys];
		if (rh <= 0.0 || rsq <= rh*rh)
			return halo_vc[sys]*halo_vc[sys] / (rsq + halo_rc_sq[sys]);
		return halo_M_rh[sys] / (rsq * sqrt(rsq));
	}

	// External halo potential energy per unit mass, Phi(r), consistent with
	// HaloScale (a_vec = HaloScale * (centre - pos) = -grad Phi). Integrating
	// dPhi/dr = HaloScale * r:
	//   inside  (r <= Rh): cored isothermal, matched at Rh to the outer branch
	//   outside (r  > Rh): Keplerian -M_rh/r, with Phi(infinity) = 0
	// For an untruncated halo (Rh <= 0) the potential has no zero at infinity
	// (it diverges logarithmically), so we reference Phi(0) = 0 instead.
	inline double HaloPotential(int sys, double rsq) const {
		double rc_sq = halo_rc_sq[sys];
		double vc2 = halo_vc[sys]*halo_vc[sys];
		double rh = halo_rh[sys];
		if (rh > 0.0) {
			double rh_sq = rh*rh;
			double M_rh = halo_M_rh[sys];
			if (rsq > rh_sq)
				return -M_rh / sqrt(rsq);
			return 0.5*vc2*log((rsq + rc_sq)/(rh_sq + rc_sq)) - M_rh/rh;
		}
		return 0.5*vc2*log((rsq + rc_sq)/rc_sq);
	}

	// Softening SQUARED (epsilon^2) for the gravitational force ON body i, added to
	// |dx|^2 by the force kernels. Softening is SINK-based: gas particles feel a
	// smoother (larger gas_softening) potential so their self-gravity cannot fragment
	// into sub-softening clumps, while stars keep the small r_soft for sharp arms.
	// The two possible squared values are precomputed by UpdateSofteningSq() whenever
	// r_soft/gas_softening change, so this hot-path lookup is just a branch + read.
	double star_soft_sq = 0.0;       // r_soft^2 (collisionless particles)
	double gas_soft_sq = 0.0;        // gas_softening^2, or r_soft^2 if gas_softening<=0
	inline double SoftSq(int i) const {
		return is_gas[i] ? gas_soft_sq : star_soft_sq;
	}
	// Recompute the cached squared softenings from r_soft / gas_softening.
	void UpdateSofteningSq() {
		star_soft_sq = r_soft * r_soft;
		gas_soft_sq  = (gas_softening > 0.0) ? gas_softening * gas_softening
		                                     : star_soft_sq;
	}

	void Allocate();
	void LoadScript(const std::string &path);
	void InitGL();
	void CreateRecordFBO(int width, int height);
	void CreateHDRFBO(int width, int height);
	GLuint CompileShader(const char *vertSrc, const char *fragSrc);
	void CalcAccelRangeP2P(int iStart, int iEnd);
	void CalcAccelRangeOct(int iStart, int iEnd);
	void ZeroAccelerationRange(int iStart, int iEnd);
	void ComputeHaloCenters();
	void IntegrateHaloCenters();
	void CalcAccelerations();
	void CalcLeapFrogPositionsRange(int iStart, int iEnd);
	void CalcLeapFrogPositions();
	void CalcLeapFrogVelocitiesAndOutputsRange(int iStart, int iEnd);
	void CalcLeapFrogVelocitiesAndOutputs();
	void CalcOutputsRange(int iStart, int iEnd);
	void CalcOutputs();
	void CalcSystemQuantities();
	void ProcessGasCollisions();
	void LogBarDiagnostic();
	// Zero one system's net momentum. Called by the procedural generators only,
	// where random particle phases leave a residual of order v_c/sqrt(N) that is
	// pure sampling noise. Systems built from explicit `Body` state vectors do not
	// call it: there the net momentum is physical (real ephemeris data).
	// Must be called BEFORE the system's bulk velocity is applied.
	void ZeroNetMomentum(int system);
	void ApplyBulkVelocity(int system);
	void ApplyBulkVelocities();
	void UpdateCameraFollow();
	void BuildOctreeForSystem(int sys, int &outFirst, int &outCount);
	void CalcAccelIsolated();
	void BuildOctreeVerts(int nodeIdx);
public:
	bool DrawOctree = false;
	bool multiThreading = true;
	int numThreads = 4;

	static void ParseDisplaySize(const std::string &scriptPath, int &width, int &height);

	Simulation(const std::string &scriptPath);
	~Simulation();

	// h_r is the exponential disc scale length, a required physical input. It is
	// independent of the truncation radius R because R/h_r is not universal:
	// Salo & Laurikainen (2000) truncate M51a at 4 h_r and M51b at 7.3 h_r.
	// Caller must ensure 0 < h_r < R.
	void LoadGalaxyDiscState(int system, double *sysPos, double *sysVel, double *discNormal, double M_central, double M_disc, double R, double Ri, double h_r, double Q, double haloVc, double haloRc, double haloRh, double sigmaZratio = 0.7, double gasMass = 0.0, double gasFraction = 0.0);
	void LoadSphericalUniverseState(int system, double *sysPos, double *sysVel, double M, double R, double H, double haloVc, double haloRc, double haloRh, double gasMass = 0.0, double gasFraction = 0.0);
	void BuildOctree();
	void Step();
	void CamMove(double d_phi, double d_theta, double d_r);
	void CamShift(double dx, double dy, double dz);
	void ReSizeGL(int width, int height);
	void DrawGL();
	void DrawInfo(double fps);
	void ReadFramePixels(uint8_t *rgbOut);
	void BlitToScreen();
	void SaveState();
	bool ReadState();

	bool GetRecordVideo() const { return Record_Video; }

	// Whether the frame about to be presented should be written to the video.
	// During the warmup phase (t < 0) the systems are isolated and held at rest,
	// so those frames are setup rather than simulation output.
	//
	// The comparison carries a half-step tolerance because t is accumulated as
	// t += dt and so never lands exactly on zero (from -2.0 at dt = 0.0005 the
	// closest value is -1.65e-13). The tolerance includes that frame as t = 0
	// while still excluding the step before it.
	bool ShouldRecordFrame() const {
		if (InitializationTime <= 0.0) return true;
		return t >= -0.5 * dt;
	}
	double GetTime() const { return t; }
	double GetEndTime() const { return EndTime; }
	int GetDisplayWidth() const { return DisplayWidth; }
	int GetDisplayHeight() const { return DisplayHeight; }
	bool GetCamOrbit() const { return CamOrbit; }
	double GetCamOrbitTheta() const { return CamOrbitTheta; }
	bool GetInfoDisplay() const { return Info_Display; }
};
