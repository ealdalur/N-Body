#include <iostream>
#include <SDL3/SDL.h>
#include <SDL3/SDL_opengl.h>
#include <vector>
#include <thread>
#include <random>
#include "VectorMath.h"
#include "BHTree.h"

const int N_BODIES = 40000;
const int N_SYSTEMS = 2;
const int N_SYSTEM_BODIES[N_SYSTEMS] = {30000, 10000};

// const int N_BODIES = 100000;
// const int N_SYSTEMS = 1;
// const int N_SYSTEM_BODIES[N_SYSTEMS] = {N_BODIES};

const int N_STATES = 6;

const bool GRAVITY_P2P = false;
const bool GRAVITY_OCT = true;

const bool SOLVER_RK4		= false;
const bool SOLVER_LEAP_FROG = true;

const bool DATA_LOG		= false;
const bool RECORD_IMG	= false;

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
	
	Camera Cam;
	GLuint dList;

	BHTree *Octree;

	FILE *DataLog;

	std::vector<std::thread> threads;
	
	void InitGL();
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
	void CalcOutputsRange(int iStart, int iEnd);
	void CalcOutputs();
public:
	bool DrawOctree = false;
	bool multiThreading = true;
	int numThreads = 4;
	
	Simulation();
	~Simulation();
	
	void LoadGalaxyDiscState(int system, double *sysPos, double *sysVel, double M, double Mfrac, double R, double Ri, double Vtol);
	void LoadSphericalUniverseState(double M, double R, double V);
	void LoadDefaultState();
	void BuildOctree();
	void Step();
	void CamMove(double d_phi, double d_theta, double d_r);
	void SetPerspective(double fovY, double aspect, double zNear, double zFar);
	void ReSizeGL(int width, int height);
	void DrawGL();
	void DrawOctant(BHTree *Oct);
	void SaveState();
	bool ReadState();
};