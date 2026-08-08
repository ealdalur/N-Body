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
	std::vector<int> halo_central;
	std::vector<int> body_system;
	std::vector<double> halo_center;

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

	std::vector<double> pos_sq;
	std::vector<double> vel_sq;
	std::vector<double> acc_sq;

	std::vector<float> pos_f;
	std::vector<int> sortedIdx;
	std::vector<int> sortTemp;
	std::vector<uint32_t> mortonCodes;
	int numActiveBodies;

	double totalKE, totalPE, totalE;
	std::vector<double> body_pot;

	bool Remove_Halo_Monopole;       // cancel net force of the rigid analytic halos

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

	int winWidth, winHeight;

	BHTree Octree;

	FILE *DataLog;

	ThreadPool *pool;

	void Allocate();
	void LoadScript(const std::string &path);
	void InitGL();
	void CreateRecordFBO(int width, int height);
	GLuint CompileShader(const char *vertSrc, const char *fragSrc);
	void CalcAccelRangeP2P(int iStart, int iEnd);
	void CalcAccelRangeOct(int iStart, int iEnd);
	void ZeroAccelerationRange(int iStart, int iEnd);
	void PinCentralBodies();
	void ComputeHaloCenters();
	void CalcAccelerations();
	void CalcLeapFrogPositionsRange(int iStart, int iEnd);
	void CalcLeapFrogPositions();
	void CalcLeapFrogVelocitiesAndOutputsRange(int iStart, int iEnd);
	void CalcLeapFrogVelocitiesAndOutputs();
	void CalcOutputsRange(int iStart, int iEnd);
	void CalcOutputs();
	void CalcEnergy();
	void RemoveHaloMonopole();
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
	void LoadGalaxyDiscState(int system, double *sysPos, double *sysVel, double *discNormal, double M, double Mfrac, double R, double Ri, double h_r, double Q, double haloVc, double haloRc);
	void LoadSphericalUniverseState(int system, double *sysPos, double *sysVel, double M, double R, double H, double haloVc, double haloRc);
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
