#pragma once

#include <iostream>
#include <SDL3/SDL.h>
#include "glad.h"
#include <vector>
#include <thread>
#include <random>
#include <algorithm>
#include "VectorMath.h"
#include "BHTree.h"
#include "ThreadPool.h"

const int N_BODIES = 800000;
const int N_SYSTEMS = 2;
const int N_SYSTEM_BODIES[N_SYSTEMS] = {600000, 200000};

// const int N_BODIES = 100000;
// const int N_SYSTEMS = 1;
// const int N_SYSTEM_BODIES[N_SYSTEMS] = {N_BODIES};

const int N_STATES = 6;

const bool GRAVITY_P2P = false;
const bool GRAVITY_OCT = true;

const bool SOLVER_RK4		= false;
const bool SOLVER_LEAP_FROG = true;

const bool DATA_LOG		= false;
const bool RECORD_VIDEO	= true;

struct Camera
{
	double pos[3];
	double phi,theta;
};

class Simulation
{
	double G;
	double FDE;
	double dt;
	double r_soft;
	double BH_Opening_Theta;

	double halo_vc[N_SYSTEMS];
	double halo_rc_sq[N_SYSTEMS];
	int halo_central[N_SYSTEMS];
	int body_system[N_BODIES];

	double mass[N_BODIES];

	double states[N_BODIES*N_STATES];
	double states_e[N_BODIES*N_STATES];
	double states_d[4][N_BODIES*N_STATES];

	double *pos[N_BODIES];
	double *pos_t[N_BODIES];
	double *vel[N_BODIES];
	double *acc[N_BODIES];
	double *acc_t[N_BODIES];
	double *acc_prev[N_BODIES];

	bool has_gravity[N_BODIES];

	double pos_sq[N_BODIES];
	double vel_sq[N_BODIES];
	double acc_sq[N_BODIES];

	float pos_f[N_BODIES * 3];
	int sortedIdx[N_BODIES];
	uint32_t mortonCodes[N_BODIES];
	int numActiveBodies;

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

	void InitGL();
	void CreateRecordFBO(int width, int height);
	GLuint CompileShader(const char *vertSrc, const char *fragSrc);
	void CalcAccelRangeP2P(int iStart, int iEnd);
	void CalcAccelRangeOct(int iStart, int iEnd);
	void PrepareDerivativeDataRange(double *s, double *s_d, int iStart, int iEnd);
	void CalcDerivatives(double *s, double *s_d);
	void CalcRK4StateEstimateRange(double *s_est, double *s_curr, double *s_d, double scalar, double dt, int iStart, int iEnd);
	void CalcRK4StateEstimate(double *s_est, double *s_curr, double *s_d, double scalar, double dt);
	void CalcLeapFrogPositionsRange(int iStart, int iEnd);
	void CalcLeapFrogPositions();
	void CalcLeapFrogVelocitiesRange(int iStart, int iEnd);
	void CalcLeapFrogVelocities();
	void CalcLeapFrogVelocitiesAndOutputsRange(int iStart, int iEnd);
	void CalcLeapFrogVelocitiesAndOutputs();
	void CalcOutputsRange(int iStart, int iEnd);
	void CalcOutputs();
	void BuildOctreeVerts(int nodeIdx);
public:
	bool DrawOctree = false;
	bool multiThreading = true;
	int numThreads = 4;

	Simulation();
	~Simulation();

	void LoadGalaxyDiscState(int system, double *sysPos, double *sysVel, double M, double Mfrac, double R, double Ri, double Vtol, double haloVc, double haloRc);
	void LoadSphericalUniverseState(double M, double R, double V);
	void LoadDefaultState();
	void BuildOctree();
	void Step();
	void CamMove(double d_phi, double d_theta, double d_r);
	void ReSizeGL(int width, int height);
	void DrawGL();
	void DrawFPS(double fps);
	void ReadFramePixels(uint8_t *rgbOut);
	void SaveState();
	bool ReadState();
};
