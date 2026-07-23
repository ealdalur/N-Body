#include "Simulation.h"
#include <cstring>
#include <fstream>
#include <sstream>

Simulation::Simulation(const std::string &scriptPath)
{
	multiThreading = true;
	int hardwareConcurrency = std::thread::hardware_concurrency();
	std::cout << "Hardware thread contexts: " << hardwareConcurrency << std::endl;
	numThreads = hardwareConcurrency - 2;
	if (numThreads < 4) numThreads = 4;
	std::cout << "Number of compute threads: " << numThreads << std::endl;

	pool = new ThreadPool(numThreads);

	N_Bodies = 0;
	N_Systems = 0;
	G = 1.0;
	FDE = 0.0;
	dt = 0.0005;
	t = 0.0;
	r_soft = 0.1;
	BH_Opening_Theta = 0.5;
	DisplayScale = 1.0;
	Gravity_P2P = false;
	Gravity_Oct = true;
	Solver_RK4 = false;
	Solver_LeapFrog = true;
	Record_Video = false;
	Data_Log = false;

	LoadScript(scriptPath);
	Allocate();

	posBuf = new float[N_Bodies * 3];
	clrBuf = new float[N_Bodies * 4];

	InitGL();

	CalcOutputs();

	if (Gravity_Oct)
		BuildOctree();

	CalcDerivatives(states.data(), states_d0.data());
	CalcOutputs();

	if (Data_Log)
	{
		DataLog = fopen("Simulation.dat", "wb");
		fwrite(&N_Bodies, sizeof(N_Bodies), 1, DataLog);
	}
}

Simulation::~Simulation()
{
	glDeleteVertexArrays(1, &particleVAO);
	glDeleteBuffers(1, &particleShapeVBO);
	glDeleteBuffers(1, &particlePosVBO);
	glDeleteBuffers(1, &particleColorVBO);
	glDeleteProgram(particleShader);
	glDeleteVertexArrays(1, &octreeVAO);
	glDeleteBuffers(1, &octreeVBO);
	glDeleteProgram(octreeShader);
	glDeleteVertexArrays(1, &fpsVAO);
	glDeleteBuffers(1, &fpsVBO);
	glDeleteProgram(fpsShader);
	if (recordFBO) {
		glDeleteFramebuffers(1, &recordFBO);
		glDeleteTextures(1, &recordColorTex);
		glDeleteRenderbuffers(1, &recordDepthRBO);
	}
	delete pool;
	delete[] posBuf;
	delete[] clrBuf;

	if (Data_Log) fclose(DataLog);
}

void Simulation::Allocate()
{
	halo_vc.resize(N_Systems, 0.0);
	halo_rc_sq.resize(N_Systems, 0.0);
	halo_central.resize(N_Systems, 0);
	halo_center.resize(N_Systems * 3, 0.0);
	body_system.resize(N_Bodies, 0);

	mass.resize(N_Bodies, 0.0);

	states.resize(N_Bodies * N_STATES, 0.0);
	states_e.resize(N_Bodies * N_STATES, 0.0);
	states_d0.resize(N_Bodies * N_STATES, 0.0);
	states_d1.resize(N_Bodies * N_STATES, 0.0);
	states_d2.resize(N_Bodies * N_STATES, 0.0);
	states_d3.resize(N_Bodies * N_STATES, 0.0);

	pos.resize(N_Bodies, nullptr);
	pos_t.resize(N_Bodies, nullptr);
	vel.resize(N_Bodies, nullptr);
	acc.resize(N_Bodies, nullptr);
	acc_t.resize(N_Bodies, nullptr);
	acc_prev.resize(N_Bodies, nullptr);

	for (int i = 0; i < N_Bodies; i++) {
		int ki = i * N_STATES;
		pos[i] = states.data() + ki;
		vel[i] = states.data() + ki + 3;
		acc[i] = states_d0.data() + ki + 3;
		acc_prev[i] = states_d1.data() + ki + 3;
	}

	has_gravity.resize(N_Bodies, true);

	pos_sq.resize(N_Bodies, 0.0);
	vel_sq.resize(N_Bodies, 0.0);
	acc_sq.resize(N_Bodies, 0.0);

	pos_f.resize(N_Bodies * 3, 0.0f);
	sortedIdx.resize(N_Bodies, 0);
	sortTemp.resize(N_Bodies, 0);
	mortonCodes.resize(N_Bodies, 0);
}

void Simulation::LoadScript(const std::string &path)
{
	std::string searchPaths[] = { path, "../" + path, "../../" + path };
	std::ifstream file;
	std::string resolvedPath;

	for (const auto &p : searchPaths) {
		file.open(p);
		if (file.is_open()) { resolvedPath = p; break; }
	}

	if (!file.is_open()) {
		std::cerr << "Error: Could not open script file: " << path << std::endl;
		std::cerr << "Searched:" << std::endl;
		for (const auto &p : searchPaths)
			std::cerr << "  " << p << std::endl;
		exit(1);
	}

	std::cout << "Loading script: " << resolvedPath << std::endl;

	std::string line;
	while (std::getline(file, line)) {
		// Strip comments
		size_t commentPos = line.find('#');
		if (commentPos != std::string::npos)
			line = line.substr(0, commentPos);

		std::istringstream iss(line);
		std::string key;
		if (!(iss >> key)) continue;

		if (key == "G") {
			iss >> G;
		} else if (key == "FDE") {
			iss >> FDE;
		} else if (key == "dt") {
			iss >> dt;
		} else if (key == "r_soft") {
			iss >> r_soft;
		} else if (key == "BH_Opening_Theta") {
			iss >> BH_Opening_Theta;
		} else if (key == "DisplayScale") {
			iss >> DisplayScale;
		} else if (key == "Solver") {
			std::string val;
			iss >> val;
			Solver_RK4 = (val == "RK4");
			Solver_LeapFrog = (val == "LeapFrog");
		} else if (key == "Gravity") {
			std::string val;
			iss >> val;
			Gravity_P2P = (val == "P2P");
			Gravity_Oct = (val == "Octree");
		} else if (key == "DataLog") {
			int val; iss >> val;
			Data_Log = (val != 0);
		} else if (key == "RecordVideo") {
			int val; iss >> val;
			Record_Video = (val != 0);
		} else if (key == "N_SystemBodies") {
			N_System_Bodies.clear();
			int n;
			while (iss >> n) {
				N_System_Bodies.push_back(n);
			}
			N_Systems = (int)N_System_Bodies.size();
			N_Bodies = 0;
			for (int i = 0; i < N_Systems; i++)
				N_Bodies += N_System_Bodies[i];
		} else if (key == "GalaxyDisc") {
			int system;
			double px, py, pz, vx, vy, vz;
			double nx, ny, nz;
			double M, Mfrac, R, Ri, Vtol, haloVc, haloRc;
			iss >> system >> px >> py >> pz >> vx >> vy >> vz >> nx >> ny >> nz >> M >> Mfrac >> R >> Ri >> Vtol >> haloVc >> haloRc;

			if (states.empty()) Allocate();

			double sysPos[3] = {px, py, pz};
			double sysVel[3] = {vx, vy, vz};
			double discNormal[3] = {nx, ny, nz};
			LoadGalaxyDiscState(system, sysPos, sysVel, discNormal, M, Mfrac, R, Ri, Vtol, haloVc, haloRc);
		} else if (key == "SphericalUniverse") {
			int system;
			double px, py, pz, vx, vy, vz;
			double M, R, H, haloVc, haloRc;
			iss >> system >> px >> py >> pz >> vx >> vy >> vz >> M >> R >> H >> haloVc >> haloRc;

			if (states.empty()) Allocate();

			double sysPos[3] = {px, py, pz};
			double sysVel[3] = {vx, vy, vz};
			LoadSphericalUniverseState(system, sysPos, sysVel, M, R, H, haloVc, haloRc);
		} else if (key == "Body") {
			int system;
			double px, py, pz, vx, vy, vz, m;
			iss >> system >> px >> py >> pz >> vx >> vy >> vz >> m;

			if (states.empty()) Allocate();

			// Find the next available slot in this system
			int sysIdx = 0;
			for (int i = 0; i < system; i++) sysIdx += N_System_Bodies[i];

			// Find first uninitialized body in this system (by checking mass == 0)
			int bi = -1;
			for (int i = 0; i < N_System_Bodies[system]; i++) {
				if (mass[sysIdx + i] == 0.0) { bi = sysIdx + i; break; }
			}
			if (bi >= 0) {
				pos[bi][0] = px; pos[bi][1] = py; pos[bi][2] = pz;
				vel[bi][0] = vx; vel[bi][1] = vy; vel[bi][2] = vz;
				mass[bi] = m;
				has_gravity[bi] = true;
				body_system[bi] = system;
			}
		} else if (key == "Camera") {
			double px, py, pz;
			iss >> px >> py >> pz;
			vset(px, py, pz, Cam.pos);
			double r = vmag(Cam.pos);
			double phi = acos(py / r);
			double theta = atan2(pz, px);
			if (sin(phi) < 1e-6)
				theta = M_PI/2;
			Cam.phi = phi - M_PI/2;
			Cam.theta = -(theta - M_PI/2);
		} else {
			std::cerr << "Warning: Unknown script key: " << key << std::endl;
		}
	}

	std::cout << "N_Bodies: " << N_Bodies << std::endl;
	std::cout << "N_Systems: " << N_Systems << std::endl;
	for (int i = 0; i < N_Systems; i++)
		std::cout << "  System " << i << ": " << N_System_Bodies[i] << " bodies" << std::endl;
}

void Simulation::LoadGalaxyDiscState(int system, double *sysPos, double *sysVel, double *discNormal, double M, double Mfrac, double R, double Ri, double Vtol, double haloVc, double haloRc){

	double m,r,theta,vm,m_orbit;
	double p[3],v[3];

	// Build orthonormal frame from disc normal
	double n[3]; vcopy(discNormal, n); vnorm(n);
	double u[3], w[3];
	// Pick a vector not parallel to n to seed the frame
	double seed[3] = {1.0, 0.0, 0.0};
	if (fabs(n[0]) > 0.9) vset(0.0, 1.0, 0.0, seed);
	vcross(n, seed, u); vnorm(u);
	vcross(n, u, w); vnorm(w);

	int sysIdx = 0;
	for (int i=0; i<system; i++) sysIdx += N_System_Bodies[i];

	halo_vc[system] = haloVc;
	halo_rc_sq[system] = haloRc * haloRc;
	halo_central[system] = sysIdx;

	for (int i=0; i<N_System_Bodies[system]; i++)
		body_system[sysIdx+i] = system;

	m = Mfrac*M/(N_System_Bodies[system]-1);
	mass[sysIdx] = M;
	has_gravity[sysIdx] = true;
	vcopy(sysPos, pos[sysIdx]);
	vcopy(sysVel, vel[sysIdx]);

	for (int i=1; i<N_System_Bodies[system]; i++)
	{
		mass[sysIdx+i] = m;
		has_gravity[sysIdx+i] = true;

		r = (R-Ri)*drand() + Ri;
		theta = 2*M_PI*drand();

		// In-plane position
		double ct = cos(theta), st = sin(theta);
		double p_plane[3];
		p_plane[0] = r*ct*u[0] + r*st*w[0];
		p_plane[1] = r*ct*u[1] + r*st*w[1];
		p_plane[2] = r*ct*u[2] + r*st*w[2];

		// Height scatter scales with r (~5% of r, matching old spherical coord behavior)
		double h = r*(M_PI/64.0)*(2*drand()-1);
		p[0] = p_plane[0] + h*n[0];
		p[1] = p_plane[1] + h*n[1];
		p[2] = p_plane[2] + h*n[2];

		m_orbit = M + m*(N_System_Bodies[system]-1)*(r*r/(R*R));
		double vc_sq = G*m_orbit/r + haloVc*haloVc*r*r/(r*r + haloRc*haloRc);
		vm = sqrt(vc_sq);
		vm = (-1.0+Vtol*(2*drand()-1))*vm;

		// Orbital velocity: perpendicular to radial (in-plane), within disc plane
		vcopy(p_plane,v); vnorm(v);
		vcross(v,n,v); vnorm(v);
		vscale(v,1.0*vm,v);
		vadd(p,pos[sysIdx],pos[sysIdx+i]);
		vadd(v,vel[sysIdx],vel[sysIdx+i]);
	}
}

void Simulation::LoadSphericalUniverseState(int system, double *sysPos, double *sysVel, double M, double R, double H, double haloVc, double haloRc) {

	double m,r,theta,phi;
	double p[3],v[3];

	int sysIdx = 0;
	for (int i=0; i<system; i++) sysIdx += N_System_Bodies[i];

	halo_vc[system] = haloVc;
	halo_rc_sq[system] = haloRc * haloRc;
	halo_central[system] = sysIdx;

	m = M/N_System_Bodies[system];

	std::mt19937 generator;
	std::normal_distribution<double> normal(0.0, 1.0);

	for (int i=0; i<N_System_Bodies[system]; i++)
	{
		mass[sysIdx+i] = m;
		has_gravity[sysIdx+i] = true;
		body_system[sysIdx+i] = system;

		r = R*pow(drand(),(1.0/3.0));
		theta = 2*M_PI*drand();
		phi = acos(2*drand()-1);
		vset(r*sin(phi)*cos(theta),r*cos(phi),r*sin(phi)*sin(theta),p);

		double pp[3];
		pp[0] = normal(generator);
		pp[1] = normal(generator);
		pp[2] = normal(generator);
		vscale(pp,10.0,pp);
		vadd(p,pp,p);

		// Hubble flow: each particle gets velocity v = H * r (radial outward)
		vscale(p, H, v);

		vadd(p, sysPos, pos[sysIdx+i]);
		vadd(v, sysVel, vel[sysIdx+i]);
	}
}

void Simulation::CalcAccelRangeP2P(int iStart, int iEnd) {

	double a[3];
	double r_halo[3];

	for (int i=iStart; i<=iEnd; i++)
	{
		vscaleadd(pos_t[i],FDE,acc_t[i]);

		for (int j=0; j<N_Bodies; j++)
		{
			if (j != i)
			{
				vsub(pos_t[j],pos_t[i],a);
				double dsq = vmagsqsoft(a, r_soft);
				double r3_inv = 1.0 / (dsq * sqrt(dsq));
				vscaleadd(a, G * mass[j] * r3_inv, acc_t[i]);
			}
		}

		int sys = body_system[i];
		double *hc = &halo_center[sys * 3];
		r_halo[0] = hc[0] - pos_t[i][0];
		r_halo[1] = hc[1] - pos_t[i][1];
		r_halo[2] = hc[2] - pos_t[i][2];
		double rsq = vmagsq(r_halo);
		if (rsq > 1e-10) {
			double halo_scale = halo_vc[sys] * halo_vc[sys] / (rsq + halo_rc_sq[sys]);
			vscaleadd(r_halo, halo_scale, acc_t[i]);
		}

		for (int s = 0; s < N_Systems; s++) {
			if (s == sys || halo_vc[s] == 0.0) continue;
			double *hc2 = &halo_center[s * 3];
			r_halo[0] = hc2[0] - pos_t[i][0];
			r_halo[1] = hc2[1] - pos_t[i][1];
			r_halo[2] = hc2[2] - pos_t[i][2];
			rsq = vmagsq(r_halo);
			double halo_scale = halo_vc[s] * halo_vc[s] / (rsq + halo_rc_sq[s]);
			vscaleadd(r_halo, halo_scale, acc_t[i]);
		}
	}
}

void Simulation::CalcAccelRangeOct(int iStart, int iEnd) {

	double a[3];
	double r_halo[3];
	float pf[3];

	for (int i=iStart; i<=iEnd; i++)
	{
		int bi = sortedIdx[i];
		vscaleadd(pos_t[bi],FDE,acc_t[bi]);

		pf[0] = (float)pos_t[bi][0];
		pf[1] = (float)pos_t[bi][1];
		pf[2] = (float)pos_t[bi][2];
		Octree.CalcAcceleration(pf, bi, (float)G, (float)r_soft, (float)(BH_Opening_Theta * BH_Opening_Theta), a);
		vadd(acc_t[bi],a,acc_t[bi]);

		int sys = body_system[bi];
		double *hc = &halo_center[sys * 3];
		r_halo[0] = hc[0] - pos_t[bi][0];
		r_halo[1] = hc[1] - pos_t[bi][1];
		r_halo[2] = hc[2] - pos_t[bi][2];
		double rsq = vmagsq(r_halo);
		if (rsq > 1e-10) {
			double halo_scale = halo_vc[sys] * halo_vc[sys] / (rsq + halo_rc_sq[sys]);
			vscaleadd(r_halo, halo_scale, acc_t[bi]);
		}

		for (int s = 0; s < N_Systems; s++) {
			if (s == sys || halo_vc[s] == 0.0) continue;
			double *hc2 = &halo_center[s * 3];
			r_halo[0] = hc2[0] - pos_t[bi][0];
			r_halo[1] = hc2[1] - pos_t[bi][1];
			r_halo[2] = hc2[2] - pos_t[bi][2];
			rsq = vmagsq(r_halo);
			double halo_scale = halo_vc[s] * halo_vc[s] / (rsq + halo_rc_sq[s]);
			vscaleadd(r_halo, halo_scale, acc_t[bi]);
		}
	}
}

void Simulation::PrepareDerivativeDataRange(double *s, double *s_d, int iStart, int iEnd) {

	int ki;

	for (int i=iStart; i<=iEnd; i++)
	{
		ki = i*N_STATES;

		pos_t[i] = s+ki;
		acc_t[i] = s_d+ki+3;

		vcopy(s+ki+3,s_d+ki);
		vset(0.0,0.0,0.0,acc_t[i]);
	}
}

void Simulation::PinCentralBodies()
{
	int sysIdx = 0;
	for (int sys = 0; sys < N_Systems; sys++) {
		if (halo_vc[sys] > 0.0) {
			int ci = halo_central[sys];
			double cx = 0.0, cy = 0.0, cz = 0.0;
			double cvx = 0.0, cvy = 0.0, cvz = 0.0;
			double total_m = 0.0;
			for (int i = 0; i < N_System_Bodies[sys]; i++) {
				double mi = mass[sysIdx + i];
				cx += mi * pos[sysIdx + i][0];
				cy += mi * pos[sysIdx + i][1];
				cz += mi * pos[sysIdx + i][2];
				cvx += mi * vel[sysIdx + i][0];
				cvy += mi * vel[sysIdx + i][1];
				cvz += mi * vel[sysIdx + i][2];
				total_m += mi;
			}
			double inv_m = 1.0 / total_m;
			pos[ci][0] = cx * inv_m;
			pos[ci][1] = cy * inv_m;
			pos[ci][2] = cz * inv_m;
			vel[ci][0] = cvx * inv_m;
			vel[ci][1] = cvy * inv_m;
			vel[ci][2] = cvz * inv_m;
		}
		sysIdx += N_System_Bodies[sys];
	}
}

void Simulation::ComputeHaloCenters()
{
	int sysIdx = 0;
	for (int sys = 0; sys < N_Systems; sys++) {
		double cx = 0.0, cy = 0.0, cz = 0.0, total_m = 0.0;
		for (int i = 0; i < N_System_Bodies[sys]; i++) {
			double mi = mass[sysIdx + i];
			cx += mi * pos_t[sysIdx + i][0];
			cy += mi * pos_t[sysIdx + i][1];
			cz += mi * pos_t[sysIdx + i][2];
			total_m += mi;
		}
		double inv_m = 1.0 / total_m;
		halo_center[sys * 3 + 0] = cx * inv_m;
		halo_center[sys * 3 + 1] = cy * inv_m;
		halo_center[sys * 3 + 2] = cz * inv_m;
		sysIdx += N_System_Bodies[sys];
	}
}

void Simulation::CalcDerivatives(double *s, double *s_d)
{
	if (multiThreading) {
		int chunk_size = N_Bodies / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_Bodies-1) : (iStart + chunk_size - 1);
			pool->submit([this, s, s_d, iStart, iEnd]() { PrepareDerivativeDataRange(s, s_d, iStart, iEnd); });
		}
		pool->waitAll();

		ComputeHaloCenters();

		if (Gravity_P2P) {
			int accel_chunk = N_Bodies / numThreads;
			for (int i = 0; i < numThreads; ++i) {
				int iStart = i * accel_chunk;
				int iEnd = (i == numThreads - 1) ? (N_Bodies-1) : (iStart + accel_chunk - 1);
				pool->submit([this, iStart, iEnd]() { CalcAccelRangeP2P(iStart, iEnd); });
			}
			pool->waitAll();
		}
		if (Gravity_Oct) {
			int accel_chunk = numActiveBodies / numThreads;
			for (int i = 0; i < numThreads; ++i) {
				int iStart = i * accel_chunk;
				int iEnd = (i == numThreads - 1) ? (numActiveBodies-1) : (iStart + accel_chunk - 1);
				pool->submit([this, iStart, iEnd]() { CalcAccelRangeOct(iStart, iEnd); });
			}
			pool->waitAll();
		}
	} else {
		PrepareDerivativeDataRange(s, s_d, 0, N_Bodies-1);
		ComputeHaloCenters();
		if (Gravity_P2P) CalcAccelRangeP2P(0, N_Bodies-1);
		if (Gravity_Oct) CalcAccelRangeOct(0, numActiveBodies-1);
	}
}

static uint32_t expandBits(uint32_t v)
{
	v = (v | (v << 16)) & 0x030000FF;
	v = (v | (v <<  8)) & 0x0300F00F;
	v = (v | (v <<  4)) & 0x030C30C3;
	v = (v | (v <<  2)) & 0x09249249;
	return v;
}

static uint32_t morton3D(float x, float y, float z, float *p_min, float *p_max)
{
	float inv = 1.0f / (p_max[0] - p_min[0]);
	uint32_t ix = (uint32_t)((x - p_min[0]) * inv * 1023.0f);
	uint32_t iy = (uint32_t)((y - p_min[1]) * inv * 1023.0f);
	uint32_t iz = (uint32_t)((z - p_min[2]) * inv * 1023.0f);
	if (ix > 1023) ix = 1023;
	if (iy > 1023) iy = 1023;
	if (iz > 1023) iz = 1023;
	return (expandBits(ix) << 2) | (expandBits(iy) << 1) | expandBits(iz);
}

void Simulation::BuildOctree()
{
	int numActive = 0;

	for (int i = 0; i < N_Bodies; i++)
	{
		if (has_gravity[i])
		{
			pos_f[i*3]   = (float)pos[i][0];
			pos_f[i*3+1] = (float)pos[i][1];
			pos_f[i*3+2] = (float)pos[i][2];
			sortedIdx[numActive++] = i;
		}
	}

	if (numActive == 0) { numActiveBodies = 0; return; }

	float p_min[3], p_max[3];
	int first = sortedIdx[0];
	p_min[0] = pos_f[first*3]; p_min[1] = pos_f[first*3+1]; p_min[2] = pos_f[first*3+2];
	p_max[0] = p_min[0]; p_max[1] = p_min[1]; p_max[2] = p_min[2];

	for (int i = 1; i < numActive; i++)
	{
		int bi = sortedIdx[i];
		for (int j = 0; j < 3; j++) {
			if (pos_f[bi*3+j] < p_min[j]) p_min[j] = pos_f[bi*3+j];
			if (pos_f[bi*3+j] > p_max[j]) p_max[j] = pos_f[bi*3+j];
		}
	}

	// Make bounding box cubic
	float maxDim = p_max[0] - p_min[0];
	if (p_max[1] - p_min[1] > maxDim) maxDim = p_max[1] - p_min[1];
	if (p_max[2] - p_min[2] > maxDim) maxDim = p_max[2] - p_min[2];
	maxDim *= 1.001f;
	float center[3] = {(p_min[0]+p_max[0])*0.5f, (p_min[1]+p_max[1])*0.5f, (p_min[2]+p_max[2])*0.5f};
	p_min[0] = center[0] - maxDim*0.5f; p_max[0] = center[0] + maxDim*0.5f;
	p_min[1] = center[1] - maxDim*0.5f; p_max[1] = center[1] + maxDim*0.5f;
	p_min[2] = center[2] - maxDim*0.5f; p_max[2] = center[2] + maxDim*0.5f;

	for (int i = 0; i < numActive; i++) {
		int bi = sortedIdx[i];
		mortonCodes[bi] = morton3D(pos_f[bi*3], pos_f[bi*3+1], pos_f[bi*3+2], p_min, p_max);
	}

	// Radix sort by Morton code (3 passes over 10-bit chunks)
	for (int shift = 0; shift < 30; shift += 10) {
		int counts[1024] = {};
		for (int i = 0; i < numActive; i++)
			counts[(mortonCodes[sortedIdx[i]] >> shift) & 0x3FF]++;
		int total = 0;
		for (int i = 0; i < 1024; i++) {
			int c = counts[i]; counts[i] = total; total += c;
		}
		for (int i = 0; i < numActive; i++) {
			int bi = sortedIdx[i];
			sortTemp[counts[(mortonCodes[bi] >> shift) & 0x3FF]++] = bi;
		}
		memcpy(sortedIdx.data(), sortTemp.data(), numActive * sizeof(int));
	}

	Octree.Reset(p_min, p_max);

	for (int i = 0; i < numActive; i++)
	{
		int bi = sortedIdx[i];
		Octree.InsertBody(pos_f.data() + bi*3, (float)mass[bi], bi);
	}

	Octree.CalcMasses();
	numActiveBodies = numActive;
}

void Simulation::CalcOutputsRange(int iStart, int iEnd) {

	for (int i=iStart; i<=iEnd; i++)
	{
		pos_sq[i] = vdot(pos[i],pos[i]);
		vel_sq[i] = vdot(vel[i],vel[i]);
		acc_sq[i] = vdot(acc[i],acc[i]);
	}
}

void Simulation::CalcOutputs() {

	if (multiThreading) {
		int chunk_size = N_Bodies / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_Bodies-1) : (iStart + chunk_size - 1);
			pool->submit([this, iStart, iEnd]() { CalcOutputsRange(iStart, iEnd); });
		}
		pool->waitAll();
	} else {
		CalcOutputsRange(0, N_Bodies-1);
	}
}

void Simulation::CalcRK4StateEstimateRange(double *s_est, double *s_curr, double *s_d, double scalar, double dt, int iStart, int iEnd) {
	for (int i=iStart; i<=iEnd; i++)
	{
		s_est[i] = s_curr[i] + scalar * s_d[i] * dt;
	}
}

void Simulation::CalcRK4StateEstimate(double *s_est, double *s_curr, double *s_d, double scalar, double dt) {

	if (multiThreading) {
		int chunk_size = N_Bodies*N_STATES / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_Bodies*N_STATES-1) : (iStart + chunk_size - 1);
			pool->submit([this, s_est, s_curr, s_d, scalar, dt, iStart, iEnd]() {
				CalcRK4StateEstimateRange(s_est, s_curr, s_d, scalar, dt, iStart, iEnd);
			});
		}
		pool->waitAll();
	} else {
		CalcRK4StateEstimateRange(s_est, s_curr, s_d, scalar, dt, 0, N_Bodies*N_STATES-1);
	}
}

void Simulation::CalcLeapFrogPositionsRange(int iStart, int iEnd) {
	for (int i=iStart; i<=iEnd; i++)
	{
		vcopy(acc[i],acc_prev[i]);

		vscaleadd(vel[i],dt,pos[i]);
		vscaleadd(acc[i],0.5*dt*dt,pos[i]);
	}
}

void Simulation::CalcLeapFrogPositions() {

	if (multiThreading) {
		int chunk_size = N_Bodies / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_Bodies-1) : (iStart + chunk_size - 1);
			pool->submit([this, iStart, iEnd]() { CalcLeapFrogPositionsRange(iStart, iEnd); });
		}
		pool->waitAll();
	} else {
		CalcLeapFrogPositionsRange(0, N_Bodies-1);
	}
}

void Simulation::CalcLeapFrogVelocitiesAndOutputsRange(int iStart, int iEnd) {
	double a[3];

	for (int i=iStart; i<=iEnd; i++)
	{
		vadd(acc[i],acc_prev[i],a);
		vscaleadd(a,0.5*dt,vel[i]);

		pos_sq[i] = vdot(pos[i],pos[i]);
		vel_sq[i] = vdot(vel[i],vel[i]);
		acc_sq[i] = vdot(acc[i],acc[i]);
	}
}

void Simulation::CalcLeapFrogVelocitiesAndOutputs() {

	if (multiThreading) {
		int chunk_size = N_Bodies / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_Bodies-1) : (iStart + chunk_size - 1);
			pool->submit([this, iStart, iEnd]() { CalcLeapFrogVelocitiesAndOutputsRange(iStart, iEnd); });
		}
		pool->waitAll();
	} else {
		CalcLeapFrogVelocitiesAndOutputsRange(0, N_Bodies-1);
	}
}

void Simulation::CalcLeapFrogVelocitiesRange(int iStart, int iEnd) {
	double a[3];

	for (int i=iStart; i<=iEnd; i++)
	{
		vadd(acc[i],acc_prev[i],a);
		vscaleadd(a,0.5*dt,vel[i]);
	}
}

void Simulation::CalcLeapFrogVelocities() {

	if (multiThreading) {
		int chunk_size = N_Bodies / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_Bodies-1) : (iStart + chunk_size - 1);
			pool->submit([this, iStart, iEnd]() { CalcLeapFrogVelocitiesRange(iStart, iEnd); });
		}
		pool->waitAll();
	} else {
		CalcLeapFrogVelocitiesRange(0, N_Bodies-1);
	}
}

void Simulation::Step()
{
	if (Solver_RK4)
	{
		const double c[4] = {0.0, 1.0/2.0, 1.0/2.0, 1.0};
		const double b[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};

		double *sd[4] = {states_d0.data(), states_d1.data(), states_d2.data(), states_d3.data()};

		CalcDerivatives(states.data(), sd[0]);
		for (int k=1; k<4; k++)
		{
			CalcRK4StateEstimate(states_e.data(), states.data(), sd[k-1], c[k], dt);
			CalcDerivatives(states_e.data(), sd[k]);
		}

		for (int k=0; k<4; k++)
		{
			CalcRK4StateEstimate(states.data(), states.data(), sd[k], b[k], dt);
		}
	}

	if (Solver_LeapFrog)
	{
		CalcLeapFrogPositions();
		PinCentralBodies();
		CalcDerivatives(states.data(), states_d0.data());
		CalcLeapFrogVelocitiesAndOutputs();
	}
	else
	{
		CalcOutputs();
	}

	if (Gravity_Oct)
	{
		BuildOctree();
	}

	if (Data_Log)
	{
		fwrite(pos_sq.data(), N_Bodies*sizeof(double), 1, DataLog);
		fwrite(vel_sq.data(), N_Bodies*sizeof(double), 1, DataLog);
		fwrite(acc_sq.data(), N_Bodies*sizeof(double), 1, DataLog);
	}

	t += dt;
}

void Simulation::CamMove(double d_phi, double d_theta, double d_r)
{

	double r = vmag(Cam.pos);
	double phi = acos(Cam.pos[1]/r);
	double theta = atan2(Cam.pos[2],Cam.pos[0]);

	// Near poles, atan2(~0,~0) is degenerate — recover theta from Cam.theta
	if (sin(phi) < 1e-6)
		theta = -(Cam.theta) + M_PI/2;

	phi += d_phi;
	if (phi < M_PI/180) phi = M_PI/180;
	if (phi > M_PI*(1.0-1.0/180)) phi = M_PI*(1.0-1.0/180);
	theta += d_theta;
	r += d_r;

	vset(
			r*sin(phi)*cos(theta),
			r*cos(phi),
			r*sin(phi)*sin(theta),
			Cam.pos
		);

	Cam.phi = phi - M_PI/2;
	Cam.theta = -(theta - M_PI/2);

}

GLuint Simulation::CompileShader(const char *vertSrc, const char *fragSrc)
{
	GLuint vert = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(vert, 1, &vertSrc, nullptr);
	glCompileShader(vert);

	GLuint frag = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(frag, 1, &fragSrc, nullptr);
	glCompileShader(frag);

	GLuint prog = glCreateProgram();
	glAttachShader(prog, vert);
	glAttachShader(prog, frag);
	glLinkProgram(prog);

	glDeleteShader(vert);
	glDeleteShader(frag);
	return prog;
}

void Simulation::InitGL()
{
	recordFBO = 0;
	recordColorTex = 0;
	recordDepthRBO = 0;

	glDisable(GL_DEPTH_TEST);
	glDepthMask(GL_FALSE);
	glDisable(GL_CULL_FACE);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);
	glEnable(GL_BLEND);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	// --- Particle shader (instanced billboards, expansion on GPU) ---
	const char *particleVert = R"(
		#version 330 core
		layout(location = 0) in vec2 aOffset;
		layout(location = 1) in vec3 aPos;
		layout(location = 2) in vec4 aColor;
		out vec4 fColor;
		uniform mat4 uVP;
		uniform vec3 uRight;
		uniform vec3 uUp;
		void main() {
			vec3 worldPos = aPos + aOffset.x * uRight + aOffset.y * uUp;
			gl_Position = uVP * vec4(worldPos, 1.0);
			fColor = aColor;
		}
	)";

	const char *particleFrag = R"(
		#version 330 core
		in vec4 fColor;
		out vec4 FragColor;
		void main() {
			FragColor = fColor;
		}
	)";

	particleShader = CompileShader(particleVert, particleFrag);

	// Star shape: 2 quads (rotated 45 deg), each as 2 triangles = 12 vertices
	const float d = 0.5f;
	const float s45 = 0.7071067f;
	float shapeData[24];
	shapeData[0]=-d; shapeData[1]=-d;
	shapeData[2]=-d; shapeData[3]= d;
	shapeData[4]= d; shapeData[5]= d;
	shapeData[6]=-d; shapeData[7]=-d;
	shapeData[8]= d; shapeData[9]= d;
	shapeData[10]=d; shapeData[11]=-d;
	float c0x=-d*s45-(-d)*s45, c0y=-d*s45+(-d)*s45;
	float c1x=-d*s45-d*s45,    c1y=-d*s45+d*s45;
	float c2x=d*s45-d*s45,     c2y=d*s45+d*s45;
	float c3x=d*s45-(-d)*s45,  c3y=d*s45+(-d)*s45;
	shapeData[12]=c0x; shapeData[13]=c0y;
	shapeData[14]=c1x; shapeData[15]=c1y;
	shapeData[16]=c2x; shapeData[17]=c2y;
	shapeData[18]=c0x; shapeData[19]=c0y;
	shapeData[20]=c2x; shapeData[21]=c2y;
	shapeData[22]=c3x; shapeData[23]=c3y;

	glGenVertexArrays(1, &particleVAO);
	glGenBuffers(1, &particleShapeVBO);
	glGenBuffers(1, &particlePosVBO);
	glGenBuffers(1, &particleColorVBO);

	glBindVertexArray(particleVAO);

	glBindBuffer(GL_ARRAY_BUFFER, particleShapeVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(shapeData), shapeData, GL_STATIC_DRAW);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, nullptr);
	glEnableVertexAttribArray(0);
	glVertexAttribDivisor(0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, particlePosVBO);
	glBufferData(GL_ARRAY_BUFFER, N_Bodies * 3 * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
	glEnableVertexAttribArray(1);
	glVertexAttribDivisor(1, 1);

	glBindBuffer(GL_ARRAY_BUFFER, particleColorVBO);
	glBufferData(GL_ARRAY_BUFFER, N_Bodies * 4 * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
	glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 0, nullptr);
	glEnableVertexAttribArray(2);
	glVertexAttribDivisor(2, 1);

	glBindVertexArray(0);

	// --- Octree wireframe shader ---
	const char *octreeVert = R"(
		#version 330 core
		layout(location = 0) in vec3 aPos;
		uniform mat4 uVP;
		void main() {
			gl_Position = uVP * vec4(aPos, 1.0);
		}
	)";
	const char *octreeFrag = R"(
		#version 330 core
		out vec4 FragColor;
		uniform vec4 uColor;
		void main() {
			FragColor = uColor;
		}
	)";
	octreeShader = CompileShader(octreeVert, octreeFrag);

	glGenVertexArrays(1, &octreeVAO);
	glGenBuffers(1, &octreeVBO);
	glBindVertexArray(octreeVAO);
	glBindBuffer(GL_ARRAY_BUFFER, octreeVBO);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
	glEnableVertexAttribArray(0);
	glBindVertexArray(0);

	// --- FPS overlay shader ---
	const char *fpsVert = R"(
		#version 330 core
		layout(location = 0) in vec2 aPos;
		uniform mat4 uProj;
		void main() {
			gl_Position = uProj * vec4(aPos, 0.0, 1.0);
		}
	)";
	const char *fpsFrag = R"(
		#version 330 core
		out vec4 FragColor;
		uniform vec4 uColor;
		void main() {
			FragColor = uColor;
		}
	)";
	fpsShader = CompileShader(fpsVert, fpsFrag);

	glGenVertexArrays(1, &fpsVAO);
	glGenBuffers(1, &fpsVBO);
	glBindVertexArray(fpsVAO);
	glBindBuffer(GL_ARRAY_BUFFER, fpsVBO);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, nullptr);
	glEnableVertexAttribArray(0);
	glBindVertexArray(0);
}

void Simulation::CreateRecordFBO(int width, int height)
{
	if (recordFBO) {
		glDeleteFramebuffers(1, &recordFBO);
		glDeleteTextures(1, &recordColorTex);
		glDeleteRenderbuffers(1, &recordDepthRBO);
	}

	glGenFramebuffers(1, &recordFBO);
	glBindFramebuffer(GL_FRAMEBUFFER, recordFBO);

	glGenTextures(1, &recordColorTex);
	glBindTexture(GL_TEXTURE_2D, recordColorTex);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, recordColorTex, 0);

	glGenRenderbuffers(1, &recordDepthRBO);
	glBindRenderbuffer(GL_RENDERBUFFER, recordDepthRBO);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, width, height);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, recordDepthRBO);

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void Simulation::ReSizeGL(int width, int height)
{
	if (height == 0) height = 1;
	winWidth = width;
	winHeight = height;
	glViewport(0, 0, width, height);

	if (Record_Video)
		CreateRecordFBO(width, height);
}

void Simulation::BuildOctreeVerts(int nodeIdx)
{
	BHNode &n = Octree.GetNode(nodeIdx);
	float x0 = (float)n.p_min[0], y0 = (float)n.p_min[1], z0 = (float)n.p_min[2];
	float x1 = (float)n.p_max[0], y1 = (float)n.p_max[1], z1 = (float)n.p_max[2];

	float edges[] = {
		x0,y0,z0, x1,y0,z0,  x1,y0,z0, x1,y1,z0,  x1,y1,z0, x0,y1,z0,  x0,y1,z0, x0,y0,z0,
		x0,y0,z1, x1,y0,z1,  x1,y0,z1, x1,y1,z1,  x1,y1,z1, x0,y1,z1,  x0,y1,z1, x0,y0,z1,
		x0,y0,z0, x0,y0,z1,  x1,y0,z0, x1,y0,z1,  x1,y1,z0, x1,y1,z1,  x0,y1,z0, x0,y1,z1
	};
	octreeVerts.insert(octreeVerts.end(), edges, edges + 72);

	for (int i = 0; i < 8; i++) {
		if (n.Octants[i] != -1) {
			BuildOctreeVerts(n.Octants[i]);
		}
	}
}

void Simulation::DrawGL()
{
	if (Record_Video)
		glBindFramebuffer(GL_FRAMEBUFFER, recordFBO);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	double cphi = cos(Cam.phi), sphi = sin(Cam.phi);
	double ctheta = cos(Cam.theta), stheta = sin(Cam.theta);

	float r00=(float)ctheta,       r10=(float)(sphi*stheta),  r20=(float)(cphi*stheta);
	float r01=0.0f,                r11=(float)cphi,           r21=(float)(-sphi);
	float r02=(float)(-stheta),    r12=(float)(sphi*ctheta),  r22=(float)(cphi*ctheta);

	float cx=(float)Cam.pos[0], cy=(float)Cam.pos[1], cz=(float)Cam.pos[2];
	float tx = -(r00*cx + r01*cy + r02*cz);
	float ty = -(r10*cx + r11*cy + r12*cz);
	float tz = -(r20*cx + r21*cy + r22*cz);

	float view[16] = {
		r00, r10, r20, 0,
		r01, r11, r21, 0,
		r02, r12, r22, 0,
		tx,  ty,  tz,  1
	};

	float aspect = (float)winWidth / (float)winHeight;
	float fov = 45.0f * (float)M_PI / 180.0f;
	float camDist = sqrtf(cx*cx + cy*cy + cz*cz);
	float zNear = 0.1f, zFar = camDist * 10.0f;
	float f = 1.0f / tanf(fov / 2.0f);
	float proj[16] = {
		f/aspect, 0, 0, 0,
		0, f, 0, 0,
		0, 0, (zFar+zNear)/(zNear-zFar), -1,
		0, 0, 2*zFar*zNear/(zNear-zFar), 0
	};

	float vp[16];
	for (int c = 0; c < 4; c++) {
		for (int r = 0; r < 4; r++) {
			vp[c*4+r] = 0;
			for (int k = 0; k < 4; k++) {
				vp[c*4+r] += proj[k*4+r] * view[c*4+k];
			}
		}
	}

	float right[3] = { r00, 0.0f, r02 };
	float up[3] = { r10, r11, r12 };

	// --- Update particle data ---
	//const float alpha_base = 100000.0f / N_Bodies;
	//const float alpha_clamped = (alpha_base > 0.2f) ? 0.2f : ((alpha_base < 0.02f) ? 0.02f : alpha_base);
	const float alpha_clamped = 0.2f;
	
	std::vector<int> sysIndices(N_Systems);
	sysIndices[0] = 0;
	for (int j = 1; j < N_Systems; j++) sysIndices[j] = sysIndices[j-1] + N_System_Bodies[j-1];

	float ds = (float)DisplayScale;
	for (int i = 0; i < N_Bodies; i++) {
		posBuf[i*3+0] = (float)pos[i][0] * ds;
		posBuf[i*3+1] = (float)pos[i][1] * ds;
		posBuf[i*3+2] = (float)pos[i][2] * ds;

		bool sysBody = false;
		for (int j = 0; j < N_Systems; j++) {
			if (i == sysIndices[j]) { sysBody = true; break; }
		}

		if (sysBody) {
			clrBuf[i*4+0] = 0.0f; clrBuf[i*4+1] = 1.0f; clrBuf[i*4+2] = 0.0f; clrBuf[i*4+3] = 1.0f;
		} else {
			float r = (float)cbrt(acc_sq[i] / 1000000.0);
			float b = 1.0f - r;
			r = (r < 0.3f) ? 0.3f : r;
			b = (b < 0.3f) ? 0.3f : b;
			clrBuf[i*4+0] = r; clrBuf[i*4+1] = 0.3f; clrBuf[i*4+2] = b; clrBuf[i*4+3] = alpha_clamped;
		}
	}

	glUseProgram(particleShader);
	glUniformMatrix4fv(glGetUniformLocation(particleShader, "uVP"), 1, GL_FALSE, vp);
	glUniform3fv(glGetUniformLocation(particleShader, "uRight"), 1, right);
	glUniform3fv(glGetUniformLocation(particleShader, "uUp"), 1, up);

	glBindVertexArray(particleVAO);
	glBindBuffer(GL_ARRAY_BUFFER, particlePosVBO);
	glBufferData(GL_ARRAY_BUFFER, N_Bodies * 3 * sizeof(float), posBuf, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, particleColorVBO);
	glBufferData(GL_ARRAY_BUFFER, N_Bodies * 4 * sizeof(float), clrBuf, GL_DYNAMIC_DRAW);
	glDrawArraysInstanced(GL_TRIANGLES, 0, 12, N_Bodies);
	glBindVertexArray(0);

	// --- Draw octree wireframe ---
	if (DrawOctree) {
		octreeVerts.clear();
		BuildOctreeVerts(0);

		glUseProgram(octreeShader);
		glUniformMatrix4fv(glGetUniformLocation(octreeShader, "uVP"), 1, GL_FALSE, vp);
		glUniform4f(glGetUniformLocation(octreeShader, "uColor"), 0.0f, 0.0f, 1.0f, 0.5f);

		glBindVertexArray(octreeVAO);
		glBindBuffer(GL_ARRAY_BUFFER, octreeVBO);
		glBufferData(GL_ARRAY_BUFFER, octreeVerts.size() * sizeof(float), octreeVerts.data(), GL_DYNAMIC_DRAW);
		glDrawArrays(GL_LINES, 0, (GLsizei)(octreeVerts.size() / 3));
		glBindVertexArray(0);
	}

	glUseProgram(0);
	glFlush();
}

void Simulation::DrawFPS(double fps)
{
	static const GLubyte font5x7[][7] = {
		{0x00,0x00,0x00,0x00,0x00,0x00,0x00}, // ' ' 32
		{0x04,0x04,0x04,0x04,0x04,0x00,0x04}, // '!' 33
		{0x0A,0x0A,0x00,0x00,0x00,0x00,0x00}, // '"' 34
		{0x0A,0x0A,0x1F,0x0A,0x1F,0x0A,0x0A}, // '#' 35
		{0x04,0x0F,0x14,0x0E,0x05,0x1E,0x04}, // '$' 36
		{0x18,0x19,0x02,0x04,0x08,0x13,0x03}, // '%' 37
		{0x08,0x14,0x14,0x08,0x15,0x12,0x0D}, // '&' 38
		{0x04,0x04,0x00,0x00,0x00,0x00,0x00}, // '\'' 39
		{0x02,0x04,0x08,0x08,0x08,0x04,0x02}, // '(' 40
		{0x08,0x04,0x02,0x02,0x02,0x04,0x08}, // ')' 41
		{0x00,0x04,0x15,0x0E,0x15,0x04,0x00}, // '*' 42
		{0x00,0x04,0x04,0x1F,0x04,0x04,0x00}, // '+' 43
		{0x00,0x00,0x00,0x00,0x00,0x04,0x08}, // ',' 44
		{0x00,0x00,0x00,0x1F,0x00,0x00,0x00}, // '-' 45
		{0x00,0x00,0x00,0x00,0x00,0x00,0x04}, // '.' 46
		{0x00,0x01,0x02,0x04,0x08,0x10,0x00}, // '/' 47
		{0x0E,0x11,0x13,0x15,0x19,0x11,0x0E}, // '0' 48
		{0x04,0x0C,0x04,0x04,0x04,0x04,0x0E}, // '1' 49
		{0x0E,0x11,0x01,0x06,0x08,0x10,0x1F}, // '2' 50
		{0x0E,0x11,0x01,0x06,0x01,0x11,0x0E}, // '3' 51
		{0x02,0x06,0x0A,0x12,0x1F,0x02,0x02}, // '4' 52
		{0x1F,0x10,0x1E,0x01,0x01,0x11,0x0E}, // '5' 53
		{0x06,0x08,0x10,0x1E,0x11,0x11,0x0E}, // '6' 54
		{0x1F,0x01,0x02,0x04,0x08,0x08,0x08}, // '7' 55
		{0x0E,0x11,0x11,0x0E,0x11,0x11,0x0E}, // '8' 56
		{0x0E,0x11,0x11,0x0F,0x01,0x02,0x0C}, // '9' 57
		{0x00,0x00,0x04,0x00,0x04,0x00,0x00}, // ':' 58
		{0x00,0x00,0x04,0x00,0x04,0x04,0x08}, // ';' 59
		{0x02,0x04,0x08,0x10,0x08,0x04,0x02}, // '<' 60
		{0x00,0x00,0x1F,0x00,0x1F,0x00,0x00}, // '=' 61
		{0x08,0x04,0x02,0x01,0x02,0x04,0x08}, // '>' 62
		{0x0E,0x11,0x01,0x02,0x04,0x00,0x04}, // '?' 63
		{0x0E,0x11,0x17,0x15,0x17,0x10,0x0E}, // '@' 64
		{0x0E,0x11,0x11,0x1F,0x11,0x11,0x11}, // 'A' 65
		{0x1E,0x11,0x11,0x1E,0x11,0x11,0x1E}, // 'B' 66
		{0x0E,0x11,0x10,0x10,0x10,0x11,0x0E}, // 'C' 67
		{0x1E,0x11,0x11,0x11,0x11,0x11,0x1E}, // 'D' 68
		{0x1F,0x10,0x10,0x1E,0x10,0x10,0x1F}, // 'E' 69
		{0x1F,0x10,0x10,0x1E,0x10,0x10,0x10}, // 'F' 70
		{0x0E,0x11,0x10,0x17,0x11,0x11,0x0F}, // 'G' 71
		{0x11,0x11,0x11,0x1F,0x11,0x11,0x11}, // 'H' 72
		{0x0E,0x04,0x04,0x04,0x04,0x04,0x0E}, // 'I' 73
		{0x07,0x02,0x02,0x02,0x02,0x12,0x0C}, // 'J' 74
		{0x11,0x12,0x14,0x18,0x14,0x12,0x11}, // 'K' 75
		{0x10,0x10,0x10,0x10,0x10,0x10,0x1F}, // 'L' 76
		{0x11,0x1B,0x15,0x15,0x11,0x11,0x11}, // 'M' 77
		{0x11,0x19,0x15,0x13,0x11,0x11,0x11}, // 'N' 78
		{0x0E,0x11,0x11,0x11,0x11,0x11,0x0E}, // 'O' 79
		{0x1E,0x11,0x11,0x1E,0x10,0x10,0x10}, // 'P' 80
		{0x0E,0x11,0x11,0x11,0x15,0x12,0x0D}, // 'Q' 81
		{0x1E,0x11,0x11,0x1E,0x14,0x12,0x11}, // 'R' 82
		{0x0E,0x11,0x10,0x0E,0x01,0x11,0x0E}, // 'S' 83
		{0x1F,0x04,0x04,0x04,0x04,0x04,0x04}, // 'T' 84
		{0x11,0x11,0x11,0x11,0x11,0x11,0x0E}, // 'U' 85
		{0x11,0x11,0x11,0x11,0x11,0x0A,0x04}, // 'V' 86
		{0x11,0x11,0x11,0x15,0x15,0x1B,0x11}, // 'W' 87
		{0x11,0x11,0x0A,0x04,0x0A,0x11,0x11}, // 'X' 88
		{0x11,0x11,0x0A,0x04,0x04,0x04,0x04}, // 'Y' 89
		{0x1F,0x01,0x02,0x04,0x08,0x10,0x1F}, // 'Z' 90
	};

	char lines[2][32];
	snprintf(lines[0], sizeof(lines[0]), "T=%.4F", t);
	snprintf(lines[1], sizeof(lines[1]), "FPS:%d", (int)(fps + 0.5));

	int maxLen = 0;
	for (int l = 0; l < 2; l++) {
		int len = (int)strlen(lines[l]);
		if (len > maxLen) maxLen = len;
	}

	std::vector<float> verts;
	float lineHeight = 10.0f;
	int scale = 1;

	for (int l = 0; l < 2; l++) {
		float x = (float)winWidth - (float)maxLen * 6.0f - 4.0f;
		float y = (float)winHeight - 12.0f - l * lineHeight;

		for (const char *p = lines[l]; *p; p++) {
			char c = *p;
			if (c < 32 || c > 90) { x += 6 * scale; continue; }
			const GLubyte *glyph = font5x7[c - 32];
			for (int row = 0; row < 7; row++) {
				for (int col = 0; col < 5; col++) {
					if (glyph[6 - row] & (1 << (4 - col))) {
						float bx = x + col * scale, by = y + row * scale;
						float bx2 = bx + scale, by2 = by + scale;
						verts.push_back(bx);  verts.push_back(by);
						verts.push_back(bx2); verts.push_back(by);
						verts.push_back(bx2); verts.push_back(by2);
						verts.push_back(bx);  verts.push_back(by);
						verts.push_back(bx2); verts.push_back(by2);
						verts.push_back(bx);  verts.push_back(by2);
					}
				}
			}
			x += 6 * scale;
		}
	}

	if (verts.empty()) return;

	float orthoProj[16] = {
		2.0f/winWidth, 0, 0, 0,
		0, 2.0f/winHeight, 0, 0,
		0, 0, -1, 0,
		-1, -1, 0, 1
	};

	glUseProgram(fpsShader);
	glUniformMatrix4fv(glGetUniformLocation(fpsShader, "uProj"), 1, GL_FALSE, orthoProj);
	glUniform4f(glGetUniformLocation(fpsShader, "uColor"), 1.0f, 1.0f, 1.0f, 0.8f);

	glBindVertexArray(fpsVAO);
	glBindBuffer(GL_ARRAY_BUFFER, fpsVBO);
	glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(float), verts.data(), GL_DYNAMIC_DRAW);
	glDrawArrays(GL_TRIANGLES, 0, (GLsizei)(verts.size() / 2));
	glBindVertexArray(0);
	glUseProgram(0);
}

void Simulation::ReadFramePixels(uint8_t *rgbOut)
{
	glBindFramebuffer(GL_READ_FRAMEBUFFER, recordFBO);
	glReadPixels(0, 0, winWidth, winHeight, GL_RGB, GL_UNSIGNED_BYTE, rgbOut);

	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	glBlitFramebuffer(0, 0, winWidth, winHeight, 0, 0, winWidth, winHeight,
					  GL_COLOR_BUFFER_BIT, GL_NEAREST);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	int rowBytes = winWidth * 3;
	for (int y = 0; y < winHeight / 2; y++) {
		uint8_t *top = rgbOut + y * rowBytes;
		uint8_t *bot = rgbOut + (winHeight - 1 - y) * rowBytes;
		for (int x = 0; x < rowBytes; x++) {
			uint8_t tmp = top[x];
			top[x] = bot[x];
			bot[x] = tmp;
		}
	}
}

void Simulation::SaveState()
{
	FILE *StateFile;

	StateFile = fopen("State.dat", "wb");

	fwrite(&G,			sizeof(double), 1, StateFile);
	fwrite(&FDE,		sizeof(double), 1, StateFile);
	fwrite(&dt,			sizeof(double), 1, StateFile);
	fwrite(&r_soft,		sizeof(double), 1, StateFile);

	fwrite(Cam.pos,		sizeof(double), 3, StateFile);
	fwrite(&Cam.phi,	sizeof(double), 1, StateFile);
	fwrite(&Cam.theta,	sizeof(double), 1, StateFile);

	fwrite(states.data(),		sizeof(double), N_Bodies*N_STATES,	StateFile);
	fwrite(mass.data(),			sizeof(double), N_Bodies,			StateFile);

	fclose(StateFile);
}

bool Simulation::ReadState()
{
	FILE *StateFile;
	size_t result;

	StateFile = fopen("State.dat", "rb");
	if (StateFile == nullptr)
	{
        return false;
    }

	result = fread(&G,			sizeof(double), 1, StateFile);
	result = fread(&FDE,		sizeof(double), 1, StateFile);
	result = fread(&dt,			sizeof(double), 1, StateFile);
	result = fread(&r_soft,		sizeof(double), 1, StateFile);

	result = fread(Cam.pos,		sizeof(double), 3, StateFile);
	result = fread(&Cam.phi,	sizeof(double), 1, StateFile);
	result = fread(&Cam.theta,	sizeof(double), 1, StateFile);

	result = fread(states.data(),		sizeof(double), N_Bodies*N_STATES,	StateFile);
	result = fread(mass.data(),			sizeof(double), N_Bodies,			StateFile);
	(void)result;

	fclose(StateFile);

	return true;
}
