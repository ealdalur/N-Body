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

const int N_STATES = 6;

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

	bool Gravity_P2P;
	bool Gravity_Oct;
	bool Record_Video;
	bool Data_Log;

	bool CamOrbit;
	double CamOrbitTheta;
	double EndTime;

	int DisplayWidth, DisplayHeight;

	std::vector<double> halo_vc;
	std::vector<double> halo_rc_sq;
	std::vector<int> halo_central;
	std::vector<int> body_system;
	std::vector<double> halo_center;

	std::vector<double> mass;

	std::vector<double> states;
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
	void CalcDerivatives();
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

	static void ParseDisplaySize(const std::string &scriptPath, int &width, int &height);

	Simulation(const std::string &scriptPath);
	~Simulation();

	void LoadGalaxyDiscState(int system, double *sysPos, double *sysVel, double *discNormal, double M, double Mfrac, double R, double Ri, double Vtol, double haloVc, double haloRc);
	void LoadSphericalUniverseState(int system, double *sysPos, double *sysVel, double M, double R, double H, double haloVc, double haloRc);
	void BuildOctree();
	void Step();
	void CamMove(double d_phi, double d_theta, double d_r);
	void CamShift(double dx, double dy, double dz);
	void ReSizeGL(int width, int height);
	void DrawGL();
	void DrawFPS(double fps);
	void ReadFramePixels(uint8_t *rgbOut);
	void BlitToScreen();
	void SaveState();
	bool ReadState();

	bool GetRecordVideo() const { return Record_Video; }
	double GetTime() const { return t; }
	double GetEndTime() const { return EndTime; }
	int GetDisplayWidth() const { return DisplayWidth; }
	int GetDisplayHeight() const { return DisplayHeight; }
	bool GetCamOrbit() const { return CamOrbit; }
	double GetCamOrbitTheta() const { return CamOrbitTheta; }
};
