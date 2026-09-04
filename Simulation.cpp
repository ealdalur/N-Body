#include "Simulation.h"
#include <cstring>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <algorithm>

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
	accel_sq_color_thresh = 1000000.0;
	ToneMapExposure = 1.2;
	Gravity_P2P = false;
	Gravity_Oct = true;
	Record_Video = false;
	Data_Log = false;
	Info_Display = true;
	totalKE = 0.0;
	totalPE = 0.0;
	totalE = 0.0;
	totalP = 0.0;
	totalL = 0.0;
	comSpeed = 0.0;
	virialRatio = -1.0;
	Remove_Halo_Monopole = true;
	// Gas defaults: fully inelastic (paper's alpha=0) and a collision radius of
	// 0.155 code units = 0.0005 primary Rd (Salo & Laurikainen's gas radius).
	// gasEnabled flips on only if a system actually declares gas particles.
	gas_alpha = 0.0;
	gas_radius = 0.155;
	gas_cell_size = 6.0;
	gas_softening = 0.0;
	gasEnabled = false;
	gasRng.seed(20250820u);       // fixed seed: reproducible collision sampling
	InitializationTime = 0.0;
	warmupActive = false;
	bulkVelocityApplied = true;   // nothing withheld unless warmup is enabled
	CamOrbit = false;
	CamOrbitTheta = 0.0;
	CamFollowSystem = -1;         // -1 = feature off, fixed Camera_lookAt
	vset(0.0, 0.0, 0.0, Cam.lookAt);
	EndTime = -1.0;
	DisplayWidth = 1280;
	DisplayHeight = 720;

	LoadScript(scriptPath);
	Allocate();
	UpdateSofteningSq();   // cache r_soft^2 / gas_softening^2 for the force kernels

	if (warmupActive) {
		std::cout << "Warmup enabled: starting at t = " << t
		          << ", systems isolated (no cross-system gravity, bulk"
		          << " velocities withheld) until t = 0" << std::endl;
	}

	posBuf = new float[N_Bodies * 3];
	clrBuf = new float[N_Bodies * 4];

	InitGL();

	if (warmupActive) {
		CalcAccelIsolated();
	} else {
		if (Gravity_Oct)
			BuildOctree();
		CalcAccelerations();
	}
	CalcOutputs();

	// Centre on the followed system before the first frame is drawn.
	UpdateCameraFollow();

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
	if (hdrFBO) {
		glDeleteFramebuffers(1, &hdrFBO);
		glDeleteTextures(1, &hdrColorTex);
		glDeleteRenderbuffers(1, &hdrDepthRBO);
	}
	glDeleteProgram(tonemapShader);
	glDeleteVertexArrays(1, &tonemapVAO);
	delete pool;
	delete[] posBuf;
	delete[] clrBuf;

	if (Data_Log) fclose(DataLog);
	if (barLog) fclose(barLog);
}

void Simulation::Allocate()
{
	halo_vc.resize(N_Systems, 0.0);
	halo_rc_sq.resize(N_Systems, 0.0);
	halo_rh.resize(N_Systems, 0.0);
	halo_M_rh.resize(N_Systems, 0.0);
	halo_central.resize(N_Systems, 0);
	halo_center.resize(N_Systems * 3, 0.0);
	halo_vel.resize(N_Systems * 3, 0.0);
	halo_acc.resize(N_Systems * 3, 0.0);
	halo_acc_prev.resize(N_Systems * 3, 0.0);
	halo_mass.resize(N_Systems, 0.0);
	system_bulk_vel.resize(N_Systems * 3, 0.0);
	body_system.resize(N_Bodies, 0);

	mass.resize(N_Bodies, 0.0);

	pos_data.resize(N_Bodies * 3, 0.0);
	vel_data.resize(N_Bodies * 3, 0.0);
	acc_data.resize(N_Bodies * 3, 0.0);
	acc_prev_data.resize(N_Bodies * 3, 0.0);

	pos.resize(N_Bodies, nullptr);
	vel.resize(N_Bodies, nullptr);
	acc.resize(N_Bodies, nullptr);
	acc_prev.resize(N_Bodies, nullptr);

	for (int i = 0; i < N_Bodies; i++) {
		pos[i] = pos_data.data() + i * 3;
		vel[i] = vel_data.data() + i * 3;
		acc[i] = acc_data.data() + i * 3;
		acc_prev[i] = acc_prev_data.data() + i * 3;
	}

	has_gravity.resize(N_Bodies, true);
	is_gas.resize(N_Bodies, 0);   // 0 = collisionless star; set per system at load

	pos_sq.resize(N_Bodies, 0.0);
	vel_sq.resize(N_Bodies, 0.0);
	acc_sq.resize(N_Bodies, 0.0);

	pos_f.resize(N_Bodies * 3, 0.0f);
	sortedIdx.resize(N_Bodies, 0);
	sortTemp.resize(N_Bodies, 0);
	mortonCodes.resize(N_Bodies, 0);
	body_pot.resize(N_Bodies, 0.0);
}

void Simulation::ParseDisplaySize(const std::string &path, int &width, int &height)
{
	width = 1280;
	height = 720;

	std::string searchPaths[] = { path, "../" + path, "../../" + path };
	std::ifstream file;
	for (const auto &p : searchPaths) {
		file.open(p);
		if (file.is_open()) break;
	}
	if (!file.is_open()) return;

	std::string line;
	while (std::getline(file, line)) {
		size_t commentPos = line.find('#');
		if (commentPos != std::string::npos)
			line = line.substr(0, commentPos);
		std::istringstream iss(line);
		std::string key;
		if (!(iss >> key)) continue;
		if (key == "Display") {
			iss >> width >> height;
			break;
		}
	}
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
			// Gravitational softening LENGTH (code units) for the collisionless
			// component; the force kernels use its square (SoftSq). Sets the smallest
			// resolved gravitational scale ~ r_soft.
			iss >> r_soft;
		} else if (key == "Gas_Softening") {
			// Gravitational softening LENGTH (code units) applied to GAS sinks only
			// (SoftSq). Set larger than r_soft to smooth gas self-gravity and stop it
			// fragmenting into balls, while stars keep the sharp r_soft. 0 = use r_soft.
			iss >> gas_softening;
		} else if (key == "BH_Opening_Theta") {
			iss >> BH_Opening_Theta;
		} else if (key == "DisplayScale") {
			iss >> DisplayScale;
		} else if (key == "accel_sq_color_thresh") {
			iss >> accel_sq_color_thresh;
		} else if (key == "ToneMapExposure") {
			iss >> ToneMapExposure;
		} else if (key == "Gravity") {
			std::string val;
			iss >> val;
			Gravity_P2P = (val == "P2P");
			Gravity_Oct = (val == "Octree");
		} else if (key == "DataLog") {
			int val; iss >> val;
			Data_Log = (val != 0);
		} else if (key == "BarDiagnostic") {
			// TEMPORARY: steps between m=2 bar/arm CSV rows; 0 (or absent) disables.
			iss >> barLogEvery;
		} else if (key == "RecordVideo") {
			int val; iss >> val;
			Record_Video = (val != 0);
		} else if (key == "Info_Display") {
			int val; iss >> val;
			Info_Display = (val != 0);
		} else if (key == "RemoveHaloMonopole") {
			int val; iss >> val;
			Remove_Halo_Monopole = (val != 0);
		} else if (key == "Gas_Restitution") {
			// Restitution coefficient for gas-gas collisions (0 = fully inelastic,
			// the paper's value). Only matters if a system declares gas particles.
			iss >> gas_alpha;
		} else if (key == "Gas_Radius") {
			// Per-gas-particle collision radius in code units (default 0.155 =
			// 0.0005 primary Rd). Sets the collision cross-section sigma=pi*(2r)^2,
			// hence the collision rate; calibrate against the sigma_gas diagnostic.
			iss >> gas_radius;
		} else if (key == "Gas_Cell_Size") {
			// Edge length (code units) of the DSMC collision cells. Must be large
			// enough to contain several in-plane gas neighbours (so in-plane
			// collisions occur) yet small compared with the disc/arm scale. The
			// collision RATE is independent of this (rate ~ local density); it only
			// sets the spatial locality of collision partners.
			iss >> gas_cell_size;
		} else if (key == "InitializationTime") {
			iss >> InitializationTime;
			if (InitializationTime > 0.0) {
				// Start before zero; the systems evolve in isolation until t=0.
				t = -InitializationTime;
				warmupActive = true;
				bulkVelocityApplied = false;
			}
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
			double M_central, M_disc, R, Ri, h_r, Q, haloVc, haloRc, haloRh;
			if (!(iss >> system >> px >> py >> pz >> vx >> vy >> vz
			          >> nx >> ny >> nz >> M_central >> M_disc >> R >> Ri >> h_r
			          >> Q >> haloVc >> haloRc >> haloRh)) {
				std::cerr << "Error: malformed GalaxyDisc command." << std::endl;
				std::cerr << "Expected 19 values: system posX posY posZ velX velY velZ"
				          << " normalX normalY normalZ centralMass discMass R Ri h_r Q"
				          << " haloVc haloRc haloRh" << std::endl;
				std::cerr << "  (optional 20th value: sigma_z/sigma_r ratio, default 0.7;" << std::endl;
				std::cerr << "   optional 21st/22nd: gas disc mass and gas fraction (0..1), default 0)" << std::endl;
				std::cerr << "  got: " << line << std::endl;
				exit(1);
			}
			// Optional 20th field: ratio of vertical to radial velocity dispersion,
			// sigma_z/sigma_r, used for the isothermal-sheet vertical structure.
			// Defaults to 0.7 (Salo & Laurikainen 2000) when omitted, so existing
			// 19-column scripts are unaffected. operator>> writes 0 on failure since
			// C++11, so read into a temp and keep the default when no token is present.
			double sigmaZratio = 0.7;
			{ double tmp; if (iss >> tmp) sigmaZratio = tmp; }
			// Optional 21st/22nd fields: gas disc mass and gas FRACTION -- the
			// fraction (0..1) of this system's bodies that are dissipative gas
			// (the rest are collisionless stars). Both default 0 (no gas). They are
			// positional, so to set them the 20th sigma_z field must also be present.
			// A fraction (rather than an absolute count) is resolution-independent:
			// change N_SystemBodies and the gas count scales automatically. The gas
			// particles are the LAST round(gasFraction*N) bodies of this system.
			double gasMass = 0.0;
			double gasFraction = 0.0;
			{ double tmp; if (iss >> tmp) gasMass = tmp; }
			{ double tmp; if (iss >> tmp) gasFraction = tmp; }
			if (h_r <= 0.0) {
				std::cerr << "Error: GalaxyDisc scale length (h_r) must be > 0, got "
				          << h_r << std::endl;
				exit(1);
			}
			if (haloRh < 0.0) {
				std::cerr << "Error: GalaxyDisc haloRh must be >= 0 (0 = untruncated), got "
				          << haloRh << std::endl;
				exit(1);
			}
			if (h_r >= R) {
				std::cerr << "Error: GalaxyDisc scale length (h_r = " << h_r
				          << ") must be less than the outer radius R = " << R
				          << std::endl;
				exit(1);
			}
			if (sigmaZratio < 0.0) {
				std::cerr << "Error: GalaxyDisc sigma_z/sigma_r ratio must be >= 0, got "
				          << sigmaZratio << std::endl;
				exit(1);
			}
			if (M_disc <= 0.0) {
				std::cerr << "Error: GalaxyDisc disc mass must be > 0, got "
				          << M_disc << std::endl;
				exit(1);
			}
			if (M_central < 0.0) {
				std::cerr << "Error: GalaxyDisc central body mass must be >= 0, got "
				          << M_central << std::endl;
				exit(1);
			}
			if (gasFraction < 0.0 || gasFraction > 1.0 || gasMass < 0.0) {
				std::cerr << "Error: GalaxyDisc gas mass must be >= 0 and gas fraction"
				          << " in [0,1], got gasMass=" << gasMass
				          << " gasFraction=" << gasFraction << std::endl;
				exit(1);
			}
			if ((gasFraction > 0.0) != (gasMass > 0.0)) {
				std::cerr << "Error: GalaxyDisc gas mass and gas fraction must both"
				          << " be positive or both zero, got gasMass=" << gasMass
				          << " gasFraction=" << gasFraction << std::endl;
				exit(1);
			}

			if (pos_data.empty()) Allocate();

			double sysPos[3] = {px, py, pz};
			double sysVel[3] = {vx, vy, vz};
			double discNormal[3] = {nx, ny, nz};
			LoadGalaxyDiscState(system, sysPos, sysVel, discNormal, M_central, M_disc, R, Ri, h_r, Q, haloVc, haloRc, haloRh, sigmaZratio, gasMass, gasFraction);
		} else if (key == "SphericalUniverse") {
			int system;
			double px, py, pz, vx, vy, vz;
			double M, R, H, haloVc, haloRc, haloRh;
			if (!(iss >> system >> px >> py >> pz >> vx >> vy >> vz
			          >> M >> R >> H >> haloVc >> haloRc >> haloRh)) {
				std::cerr << "Error: malformed SphericalUniverse command." << std::endl;
				std::cerr << "Expected >=13 values: system posX posY posZ velX velY velZ"
				          << " mass radius H haloVc haloRc haloRh [gasMass] [gasFraction]" << std::endl;
				std::cerr << "  got: " << line << std::endl;
				exit(1);
			}

			// Optional 14th/15th fields: dissipative gas mass and gas FRACTION -- the
			// fraction (0..1) of this system's bodies that are gas (the rest are
			// collisionless). Both default 0 (no gas), so existing 13-field scripts are
			// unaffected. When set, M is the COLLISIONLESS mass budget and gasMass the
			// gas budget (total = M + gasMass); the gas are the LAST round(gasFraction*N)
			// bodies. Same convention as GalaxyDisc.
			double gasMass = 0.0;
			double gasFraction = 0.0;
			{ double tmp; if (iss >> tmp) gasMass = tmp; }
			{ double tmp; if (iss >> tmp) gasFraction = tmp; }
			if (gasFraction < 0.0 || gasFraction > 1.0 || gasMass < 0.0) {
				std::cerr << "Error: SphericalUniverse gas mass must be >= 0 and gas"
				          << " fraction in [0,1], got gasMass=" << gasMass
				          << " gasFraction=" << gasFraction << std::endl;
				exit(1);
			}
			if ((gasFraction > 0.0) != (gasMass > 0.0)) {
				std::cerr << "Error: SphericalUniverse gas mass and gas fraction must"
				          << " both be positive or both zero." << std::endl;
				exit(1);
			}

			if (pos_data.empty()) Allocate();

			double sysPos[3] = {px, py, pz};
			double sysVel[3] = {vx, vy, vz};
			LoadSphericalUniverseState(system, sysPos, sysVel, M, R, H, haloVc, haloRc, haloRh, gasMass, gasFraction);
		} else if (key == "Body") {
			int system;
			double px, py, pz, vx, vy, vz, m;
			iss >> system >> px >> py >> pz >> vx >> vy >> vz >> m;

			if (pos_data.empty()) Allocate();

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
		} else if (key == "Display") {
			iss >> DisplayWidth >> DisplayHeight;
		} else if (key == "End_Time") {
			iss >> EndTime;
		} else if (key == "Camera_Orbit") {
			double theta;
			iss >> theta;
			CamOrbit = true;
			CamOrbitTheta = theta;
		} else if (key == "Camera_lookAt_System") {
			int sys;
			if (!(iss >> sys)) {
				std::cerr << "Error: Camera_lookAt_System requires a system index."
				          << std::endl;
				std::cerr << "  got: " << line << std::endl;
				exit(1);
			}
			if (sys < 0) {
				std::cerr << "Error: Camera_lookAt_System index must be >= 0, got "
				          << sys << std::endl;
				exit(1);
			}
			// N_Systems is not known until N_SystemBodies is parsed, so the upper
			// bound is checked after the script is fully read.
			CamFollowSystem = sys;
		} else if (key == "Camera_lookAt") {
			double lx, ly, lz;
			iss >> lx >> ly >> lz;
			vset(lx, ly, lz, Cam.lookAt);
		} else if (key == "Camera") {
			double px, py, pz;
			iss >> px >> py >> pz;
			vset(px, py, pz, Cam.pos);
		} else {
			std::cerr << "Warning: Unknown script key: " << key << std::endl;
		}
	}

	double rel[3];
	vsub(Cam.pos, Cam.lookAt, rel);
	double r = vmag(rel);
	double phi = acos(rel[1] / r);
	double theta = atan2(rel[2], rel[0]);
	if (sin(phi) < 1e-6)
		theta = M_PI/2;
	Cam.phi = phi - M_PI/2;
	Cam.theta = -(theta - M_PI/2);

	std::cout << "N_Bodies: " << N_Bodies << std::endl;
	std::cout << "N_Systems: " << N_Systems << std::endl;
	for (int i = 0; i < N_Systems; i++)
		std::cout << "  System " << i << ": " << N_System_Bodies[i] << " bodies" << std::endl;

	// Camera_lookAt_System can appear before N_SystemBodies, so the upper bound
	// is only knowable here.
	if (CamFollowSystem >= 0) {
		if (CamFollowSystem >= N_Systems) {
			std::cerr << "Error: Camera_lookAt_System index " << CamFollowSystem
			          << " is out of range; there are " << N_Systems
			          << " system(s) (valid indices 0.." << N_Systems - 1 << ")"
			          << std::endl;
			exit(1);
		}
		std::cout << "Camera following system " << CamFollowSystem
		          << " (look-at retargeted to its central body each frame)"
		          << std::endl;
	}
}

void Simulation::LoadGalaxyDiscState(int system, double *sysPos, double *sysVel, double *discNormal, double M_central, double M_disc, double R, double Ri, double h_r, double Q, double haloVc, double haloRc, double haloRh, double sigmaZratio, double gasMass, double gasFraction){

	double r,theta;
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
	// Halo truncation. <= 0 leaves the halo untruncated. When set, the enclosed
	// mass is frozen at M_halo(Rh) beyond Rh, so precompute that constant.
	halo_rh[system] = haloRh;
	halo_M_rh[system] = (haloRh > 0.0)
	    ? haloVc*haloVc * haloRh*haloRh*haloRh / (haloRh*haloRh + haloRc*haloRc)
	    : 0.0;
	halo_central[system] = sysIdx;

	// The halo centre is now an inertial body integrated under gravity (see
	// IntegrateHaloCenters), not re-derived from the particle barycentre. Seed it
	// at the galaxy's start position, at rest (the bulk velocity is applied at
	// t=0 with the particles). Its inertial mass is the total (truncated) halo
	// mass; for an untruncated halo we use the mass enclosed within the disc
	// radius R as a finite proxy.
	halo_center[system*3+0] = sysPos[0];
	halo_center[system*3+1] = sysPos[1];
	halo_center[system*3+2] = sysPos[2];
	halo_vel[system*3+0] = halo_vel[system*3+1] = halo_vel[system*3+2] = 0.0;
	halo_acc[system*3+0] = halo_acc[system*3+1] = halo_acc[system*3+2] = 0.0;
	halo_acc_prev[system*3+0] = halo_acc_prev[system*3+1] = halo_acc_prev[system*3+2] = 0.0;
	if (haloVc > 0.0)
		halo_mass[system] = (haloRh > 0.0)
		    ? halo_M_rh[system]
		    : haloVc*haloVc * R*R*R / (R*R + haloRc*haloRc);
	else
		halo_mass[system] = 0.0;

	for (int i=0; i<N_System_Bodies[system]; i++)
		body_system[sysIdx+i] = system;

	// Gas ("sticky") particles are the LAST nGas bodies of this system's sub-array,
	// so the central body (index sysIdx, first) becomes gas only if nGas == N. The
	// count is derived from the gas FRACTION so it scales with the particle count.
	// The disc mass M_disc is the STAR budget and gasMass the GAS budget; both trace
	// the same exponential profile, so the disc's gravitating mass is their sum.
	int Nsys = N_System_Bodies[system];
	int nGas = (int)llround(gasFraction * Nsys);
	if (nGas < 0) nGas = 0;
	if (nGas > Nsys) nGas = Nsys;
	int gasStart = sysIdx + Nsys - nGas;                 // first gas index
	int nStarDisc = (Nsys - nGas - 1 > 0) ? (Nsys - nGas - 1) : 0;
	double m_star = (nStarDisc > 0) ? M_disc / nStarDisc : 0.0;
	double m_gas  = (nGas > 0) ? gasMass / nGas : 0.0;
	if (nGas > 0) gasEnabled = true;

	bool centralIsGas = (nGas >= Nsys);
	mass[sysIdx] = centralIsGas ? m_gas : M_central;
	is_gas[sysIdx] = centralIsGas ? 1 : 0;
	has_gravity[sysIdx] = true;
	vcopy(sysPos, pos[sysIdx]);

	// Record this system's bulk velocity. Particles are built WITHOUT it so that
	// ZeroNetMomentum() below sees only the sampling noise; the bulk motion is
	// added afterwards (or, under warmup, withheld until t = 0). The bulk
	// POSITION is always applied -- the systems must be spatially separated.
	system_bulk_vel[system*3+0] = sysVel[0];
	system_bulk_vel[system*3+1] = sysVel[1];
	system_bulk_vel[system*3+2] = sysVel[2];

	vset(0.0, 0.0, 0.0, vel[sysIdx]);

	// Exponential disc: surface density Sigma(r) = Sigma_0 * exp(-r/h_r).
	//
	// The radial NUMBER density is the surface density times the area of a ring,
	//   p(r) dr = Sigma(r) * 2*pi*r dr  ~  r * exp(-r/h_r) dr
	// which is exactly a Gamma(shape=2, scale=h_r) distribution. A Gamma with
	// integer shape 2 is the sum of two independent exponentials, so
	//   r = -h_r * (ln u1 + ln u2) = -h_r * ln(u1*u2)
	// samples it exactly, with no numerical inversion. Draws outside [Ri, R] are
	// rejected (~10% for a disc truncated at 4 h_r).
	//
	// h_r is supplied by the caller and validated at parse time (0 < h_r < R).
	// It is independent of R because R/h_r is not a universal ratio: Salo &
	// Laurikainen (2000) truncate M51a at 4 h_r and M51b at 7.3 h_r. A disc
	// truncated at 4 h_r retains ~91% of its mass.

	// Enclosed mass fraction for the exponential disc:
	// M_enc(r) / M_disc = [1 - (1 + r/h_r)*exp(-r/h_r)] / [1 - (1 + R/h_r)*exp(-R/h_r)]
	double enc_denom = 1.0 - (1.0 + R/h_r)*exp(-R/h_r);

	// Gaussian RNG for velocity dispersion
	std::mt19937 gen(42 + system);
	std::normal_distribution<double> normal(0.0, 1.0);

	// Local disc kinematics at radius rr. Factored into a lambda so the
	// asymmetric-drift correction below can finite-difference d ln(Sigma*sigma_r^2)
	// with respect to R. Returns:
	//   vc_sq   circular speed squared (enclosed disc mass + halo)
	//   Omega   angular speed vc/rr
	//   Sigma   exponential surface density, normalized to the TRUNCATED disc mass
	//   kappa   epicyclic frequency from the full rotation curve
	//   sig_r   Toomre-Q radial dispersion; sig_phi = sig_r*kappa/(2*Omega)
	// Total disc mass driving the rotation curve and surface density: stars plus
	// gas (both sampled from the same exponential profile), so the disc gravity
	// and Toomre-Q dispersions see the full disc, not just the star budget.
	double M_disc_total = M_disc + gasMass;
	double haloRc_sq = haloRc * haloRc;
	auto discProps = [&](double rr, double &vc_sq, double &Omega, double &Sigma,
	                     double &kappa, double &sig_r, double &sig_phi){
		double enc = (1.0 - (1.0 + rr/h_r)*exp(-rr/h_r)) / enc_denom;
		double morb = M_central + M_disc_total*enc;
		vc_sq = G*morb/rr + haloVc*haloVc*rr*rr/(rr*rr + haloRc_sq);
		Omega = sqrt(vc_sq) / rr;
		// Sigma_0 = M_disc_total / (2*pi*h_r^2 * enc_denom): the /enc_denom factor
		// makes the surface density integrate to the total disc mass over the
		// truncated disc, consistent with dMenc_dr (which carries the same factor).
		Sigma = (M_disc_total / (2.0*M_PI*h_r*h_r*enc_denom)) * exp(-rr/h_r);
		double dMenc_dr = M_disc_total * (rr/(h_r*h_r))*exp(-rr/h_r) / enc_denom;
		double dvc_sq_dr = -G*morb/(rr*rr) + G*dMenc_dr/rr
		                   + haloVc*haloVc * 2.0*rr*haloRc_sq
		                     / ((rr*rr + haloRc_sq)*(rr*rr + haloRc_sq));
		double kappa_sq = dvc_sq_dr/rr + 2.0*vc_sq/(rr*rr);
		if (kappa_sq < 0.0) kappa_sq = 4.0*Omega*Omega;
		kappa = sqrt(kappa_sq);
		sig_r = Q * 3.36 * G * Sigma / kappa;
		sig_phi = sig_r * kappa / (2.0*Omega);
	};

	for (int i=1; i<N_System_Bodies[system]; i++)
	{
		// Star vs gas: the last nGas bodies (index >= gasStart) are gas particles,
		// carrying the gas per-particle mass; the rest are collisionless stars.
		is_gas[sysIdx+i] = (sysIdx+i >= gasStart) ? 1 : 0;
		mass[sysIdx+i] = is_gas[sysIdx+i] ? m_gas : m_star;
		has_gravity[sysIdx+i] = true;

		// Sample radius from the exponential disc profile, r ~ Gamma(2, h_r),
		// rejecting draws outside the disc. drand() is guarded against returning
		// exactly 0, which would give log(0).
		do {
			double u1 = drand(), u2 = drand();
			if (u1 < 1e-300) u1 = 1e-300;
			if (u2 < 1e-300) u2 = 1e-300;
			r = -h_r * log(u1 * u2);
		} while (r < Ri || r > R);
		theta = 2*M_PI*drand();

		// In-plane position
		double ct = cos(theta), st = sin(theta);
		double p_plane[3];
		p_plane[0] = r*ct*u[0] + r*st*w[0];
		p_plane[1] = r*ct*u[1] + r*st*w[1];
		p_plane[2] = r*ct*u[2] + r*st*w[2];

		// Local kinematics at this radius (circular speed, dispersions, kappa).
		double vc_sq, Omega, Sigma, kappa, sigma_r, sigma_phi;
		discProps(r, vc_sq, Omega, Sigma, kappa, sigma_r, sigma_phi);

		// Vertical structure: self-gravitating isothermal sheet (Spitzer 1942),
		//   rho(z) = rho0 * sech^2(z/z0),  z0 = sigma_z^2 / (2*pi*G*Sigma),
		// with sigma_z/sigma_r = sigmaZratio (default 0.7, Salo & Laurikainen 2000,
		// sect 2.2; set per-disc via the optional 20th GalaxyDisc field). Sample z
		// by inverting the CDF F(z) = (1 + tanh(z/z0))/2, i.e. z = z0*atanh(2u-1).
		double sigma_z = sigmaZratio * sigma_r;
		double z0 = sigma_z*sigma_z / (2.0*M_PI*G*Sigma);
		double uz = drand();
		if (uz < 1e-6) uz = 1e-6;
		if (uz > 1.0 - 1e-6) uz = 1.0 - 1e-6;
		double zc = (z0 > 0.0)
		    ? z0 * 0.5 * log((1.0 + (2.0*uz - 1.0)) / (1.0 - (2.0*uz - 1.0)))
		    : 0.0;
		p[0] = p_plane[0] + zc*n[0];
		p[1] = p_plane[1] + zc*n[1];
		p[2] = p_plane[2] + zc*n[2];

		// Asymmetric drift: the mean streaming (tangential) velocity is LOWER than
		// the circular speed by the pressure support of the random motions
		// (Binney & Tremaine, asymmetric drift equation):
		//   vc^2 - <v_phi>^2 = sigma_r^2 [ sigma_phi^2/sigma_r^2 - 1
		//                                  - d ln(Sigma sigma_r^2)/d ln R ]
		// The log-derivative is taken numerically from discProps().
		double v_stream;
		{
			double hs = 1.0e-3 * r;
			double vcs, Om, ka, sp, Sg1, sr1, Sg2, sr2;
			discProps(r + hs, vcs, Om, Sg1, ka, sr1, sp);
			discProps(r - hs, vcs, Om, Sg2, ka, sr2, sp);
			double f1 = Sg1*sr1*sr1, f2 = Sg2*sr2*sr2;
			double dln_f = (log(f1) - log(f2)) / (log(r + hs) - log(r - hs));
			double bracket = (sigma_phi*sigma_phi)/(sigma_r*sigma_r) - 1.0 - dln_f;
			double v_stream_sq = vc_sq - sigma_r*sigma_r*bracket;
			v_stream = (v_stream_sq > 0.0) ? sqrt(v_stream_sq) : 0.0;
		}

		// Gaussian random velocities about the drift-corrected streaming speed.
		double v_tan  = -(v_stream + sigma_phi * normal(gen));
		double v_rad  = sigma_r * normal(gen);
		double v_vert = sigma_z * normal(gen);

		// Tangential direction: perpendicular to radial, in disc plane
		double t_hat[3], r_hat[3];
		vcopy(p_plane, r_hat); vnorm(r_hat);
		vcross(r_hat, n, t_hat); vnorm(t_hat);

		v[0] = v_tan * t_hat[0] + v_rad * r_hat[0] + v_vert * n[0];
		v[1] = v_tan * t_hat[1] + v_rad * r_hat[1] + v_vert * n[1];
		v[2] = v_tan * t_hat[2] + v_rad * r_hat[2] + v_vert * n[2];

		vadd(p, pos[sysIdx], pos[sysIdx+i]);
		vcopy(v, vel[sysIdx+i]);
	}

	// Random particle phases leave a net momentum of order v_c/sqrt(N). Remove it
	// now, while the particles still carry only their internal motion.
	ZeroNetMomentum(system);

	// Apply the bulk motion, unless warmup is deferring it to t = 0.
	if (!warmupActive) ApplyBulkVelocity(system);
}

void Simulation::LoadSphericalUniverseState(int system, double *sysPos, double *sysVel, double M, double R, double H, double haloVc, double haloRc, double haloRh, double gasMass, double gasFraction) {

	double r,theta,phi;
	double p[3],v[3];

	int sysIdx = 0;
	for (int i=0; i<system; i++) sysIdx += N_System_Bodies[i];

	halo_vc[system] = haloVc;
	halo_rc_sq[system] = haloRc * haloRc;
	// Halo truncation. <= 0 leaves the halo untruncated. When set, the enclosed
	// mass is frozen at M_halo(Rh) beyond Rh, so precompute that constant.
	halo_rh[system] = haloRh;
	halo_M_rh[system] = (haloRh > 0.0)
	    ? haloVc*haloVc * haloRh*haloRh*haloRh / (haloRh*haloRh + haloRc*haloRc)
	    : 0.0;
	halo_central[system] = sysIdx;

	// Inertial halo centre (see IntegrateHaloCenters / LoadGalaxyDiscState).
	halo_center[system*3+0] = sysPos[0];
	halo_center[system*3+1] = sysPos[1];
	halo_center[system*3+2] = sysPos[2];
	halo_vel[system*3+0] = halo_vel[system*3+1] = halo_vel[system*3+2] = 0.0;
	halo_acc[system*3+0] = halo_acc[system*3+1] = halo_acc[system*3+2] = 0.0;
	halo_acc_prev[system*3+0] = halo_acc_prev[system*3+1] = halo_acc_prev[system*3+2] = 0.0;
	if (haloVc > 0.0)
		halo_mass[system] = (haloRh > 0.0)
		    ? halo_M_rh[system]
		    : haloVc*haloVc * R*R*R / (R*R + haloRc*haloRc);
	else
		halo_mass[system] = 0.0;

	// Split into a collisionless budget (M) and a dissipative gas budget (gasMass).
	// The gas particles are the LAST nGas bodies of this system (same convention as
	// GalaxyDisc); since every particle is sampled independently from the same
	// uniform sphere + Hubble flow, that slice is just a spatially random subset.
	int Nsys = N_System_Bodies[system];
	int nGas = (int)llround(gasFraction * Nsys);
	if (nGas < 0) nGas = 0;
	if (nGas > Nsys) nGas = Nsys;
	int gasStart = sysIdx + Nsys - nGas;                 // first gas index
	int nColl = Nsys - nGas;
	double m_coll = (nColl > 0) ? M / nColl : 0.0;
	double m_gas  = (nGas > 0) ? gasMass / nGas : 0.0;
	if (nGas > 0) gasEnabled = true;

	// Bulk velocity is recorded and applied after the momentum correction; see
	// LoadGalaxyDiscState for the rationale.
	system_bulk_vel[system*3+0] = sysVel[0];
	system_bulk_vel[system*3+1] = sysVel[1];
	system_bulk_vel[system*3+2] = sysVel[2];

	std::mt19937 generator;
	std::normal_distribution<double> normal(0.0, 1.0);

	for (int i=0; i<N_System_Bodies[system]; i++)
	{
		bool isGas = (sysIdx + i) >= gasStart && nGas > 0;
		mass[sysIdx+i] = isGas ? m_gas : m_coll;
		is_gas[sysIdx+i] = isGas ? 1 : 0;
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
		vcopy(v, vel[sysIdx+i]);
	}

	ZeroNetMomentum(system);
	if (!warmupActive) ApplyBulkVelocity(system);
}

void Simulation::CalcAccelRangeP2P(int iStart, int iEnd) {

	double a[3];
	double r_halo[3];

	for (int i=iStart; i<=iEnd; i++)
	{
		vscaleadd(pos[i],FDE,acc[i]);

		double pot = 0.0;
		double softI = SoftSq(i);          // sink-based: constant over the j-loop
		for (int j=0; j<N_Bodies; j++)
		{
			if (j != i)
			{
				vsub(pos[j],pos[i],a);
				double dsq = vmagsqsoft(a, softI);
				double r_inv = 1.0 / sqrt(dsq);
				double r3_inv = r_inv / dsq;
				vscaleadd(a, G * mass[j] * r3_inv, acc[i]);
				pot += -G * mass[j] * r_inv;
			}
		}
		body_pot[i] = pot * mass[i];

		int sys = body_system[i];
		double *hc = &halo_center[sys * 3];
		r_halo[0] = hc[0] - pos[i][0];
		r_halo[1] = hc[1] - pos[i][1];
		r_halo[2] = hc[2] - pos[i][2];
		double rsq = vmagsq(r_halo);
		if (rsq > 1e-10) {
			double halo_scale = HaloScale(sys, rsq);
			vscaleadd(r_halo, halo_scale, acc[i]);
		}

		for (int s = 0; s < N_Systems; s++) {
			if (s == sys || halo_vc[s] == 0.0) continue;
			double *hc2 = &halo_center[s * 3];
			r_halo[0] = hc2[0] - pos[i][0];
			r_halo[1] = hc2[1] - pos[i][1];
			r_halo[2] = hc2[2] - pos[i][2];
			rsq = vmagsq(r_halo);
			double halo_scale = HaloScale(s, rsq);
			vscaleadd(r_halo, halo_scale, acc[i]);
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
		vscaleadd(pos[bi],FDE,acc[bi]);

		pf[0] = (float)pos[bi][0];
		pf[1] = (float)pos[bi][1];
		pf[2] = (float)pos[bi][2];
		double pot;
		Octree.CalcAcceleration(pf, bi, (float)G, (float)SoftSq(bi), (float)(BH_Opening_Theta * BH_Opening_Theta), a, &pot);
		vadd(acc[bi],a,acc[bi]);
		body_pot[bi] = pot * mass[bi];

		int sys = body_system[bi];
		double *hc = &halo_center[sys * 3];
		r_halo[0] = hc[0] - pos[bi][0];
		r_halo[1] = hc[1] - pos[bi][1];
		r_halo[2] = hc[2] - pos[bi][2];
		double rsq = vmagsq(r_halo);
		if (rsq > 1e-10) {
			double halo_scale = HaloScale(sys, rsq);
			vscaleadd(r_halo, halo_scale, acc[bi]);
		}

		for (int s = 0; s < N_Systems; s++) {
			if (s == sys || halo_vc[s] == 0.0) continue;
			double *hc2 = &halo_center[s * 3];
			r_halo[0] = hc2[0] - pos[bi][0];
			r_halo[1] = hc2[1] - pos[bi][1];
			r_halo[2] = hc2[2] - pos[bi][2];
			rsq = vmagsq(r_halo);
			double halo_scale = HaloScale(s, rsq);
			vscaleadd(r_halo, halo_scale, acc[bi]);
		}
	}
}

void Simulation::ZeroAccelerationRange(int iStart, int iEnd) {

	for (int i=iStart; i<=iEnd; i++)
	{
		vset(0.0,0.0,0.0,acc[i]);
	}
}

void Simulation::ComputeHaloCenters()
{
	int sysIdx = 0;
	for (int sys = 0; sys < N_Systems; sys++) {
		double cx = 0.0, cy = 0.0, cz = 0.0, total_m = 0.0;
		for (int i = 0; i < N_System_Bodies[sys]; i++) {
			double mi = mass[sysIdx + i];
			cx += mi * pos[sysIdx + i][0];
			cy += mi * pos[sysIdx + i][1];
			cz += mi * pos[sysIdx + i][2];
			total_m += mi;
		}
		double inv_m = 1.0 / total_m;
		halo_center[sys * 3 + 0] = cx * inv_m;
		halo_center[sys * 3 + 1] = cy * inv_m;
		halo_center[sys * 3 + 2] = cz * inv_m;
		sysIdx += N_System_Bodies[sys];
	}
}

void Simulation::IntegrateHaloCenters()
{
	// Acceleration of each halo centre, treated as an inertial body (Salo &
	// Laurikainen 2000: "the disc back-action is taken into account in the halo
	// motion"). Two momentum-conserving contributions:
	//   1. Disc back-reaction. Every particle the halo pulls pulls back on it
	//      (Newton's 3rd law). The net force the halo exerts on all gravitating
	//      particles is fh_s = sum_i m_i * HaloScale_s * (hc_s - pos_i); the
	//      reaction on the halo is -fh_s. This REPLACES RemoveHaloMonopole: rather
	//      than deleting that net force from the particles, the halo recoils under
	//      it, so the halo is never dragged around by a tidal tail.
	//   2. Halo-halo. Each halo centre falls in every other halo's field exactly
	//      as a particle there would. Beyond both truncation radii the pair is
	//      equal-and-opposite point masses, so total momentum is conserved.
	// Together with the particle forces this reproduces the analytic relative
	// orbit (compute_M51.py) while leaving only the real, weak disc friction.
	int NS = N_Systems;
	std::vector<double> fh(NS * 3, 0.0);

	for (int i = 0; i < N_Bodies; i++) {
		if (!has_gravity[i]) continue;
		double mi = mass[i];
		for (int s = 0; s < NS; s++) {
			if (halo_vc[s] == 0.0) continue;
			double *hc = &halo_center[s * 3];
			double r_halo[3] = { hc[0]-pos[i][0], hc[1]-pos[i][1], hc[2]-pos[i][2] };
			double rsq = vmagsq(r_halo);
			if (rsq <= 1e-10) continue;
			double scale = HaloScale(s, rsq);
			fh[s*3+0] += mi * scale * r_halo[0];
			fh[s*3+1] += mi * scale * r_halo[1];
			fh[s*3+2] += mi * scale * r_halo[2];
		}
	}

	for (int s = 0; s < NS; s++) {
		if (halo_vc[s] == 0.0 || halo_mass[s] <= 0.0) {
			halo_acc[s*3+0] = halo_acc[s*3+1] = halo_acc[s*3+2] = 0.0;
			continue;
		}
		double inv = 1.0 / halo_mass[s];
		double a[3] = { -fh[s*3+0]*inv, -fh[s*3+1]*inv, -fh[s*3+2]*inv };
		double *hcs = &halo_center[s * 3];
		for (int o = 0; o < NS; o++) {
			if (o == s || halo_vc[o] == 0.0) continue;
			double *hco = &halo_center[o * 3];
			double d[3] = { hco[0]-hcs[0], hco[1]-hcs[1], hco[2]-hcs[2] };
			double rsq = vmagsq(d);
			if (rsq <= 1e-10) continue;
			double scale = HaloScale(o, rsq);
			a[0] += scale * d[0];
			a[1] += scale * d[1];
			a[2] += scale * d[2];
		}
		halo_acc[s*3+0] = a[0];
		halo_acc[s*3+1] = a[1];
		halo_acc[s*3+2] = a[2];
	}
}

void Simulation::CalcAccelerations()
{
	if (multiThreading) {
		int chunk_size = N_Bodies / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_Bodies-1) : (iStart + chunk_size - 1);
			pool->submit([this, iStart, iEnd]() { ZeroAccelerationRange(iStart, iEnd); });
		}
		pool->waitAll();

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
		ZeroAccelerationRange(0, N_Bodies-1);
		if (Gravity_P2P) CalcAccelRangeP2P(0, N_Bodies-1);
		if (Gravity_Oct) CalcAccelRangeOct(0, numActiveBodies-1);
	}

	// Halo centres are inertial bodies: derive their acceleration from the disc
	// back-reaction and the other halos. This replaces the old RemoveHaloMonopole
	// band-aid -- the net halo force is now balanced by the halo recoiling rather
	// than being deleted from the particles, so momentum is conserved and the halo
	// is never dragged off the nucleus by tidal debris.
	IntegrateHaloCenters();
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

void Simulation::CalcSystemQuantities() {
	totalKE = 0.0;
	double selfPE = 0.0;   // particle-particle self-gravity (each pair counted twice)
	double extPE  = 0.0;   // external halo + FDE potential (counted once per body)

	// First pass: energy, total mass, mass-weighted position (COM), momentum
	double M = 0.0;
	double com[3] = {0.0, 0.0, 0.0};   // Sigma m_i * r_i
	double P[3]   = {0.0, 0.0, 0.0};   // Sigma m_i * v_i (total momentum)
	for (int i = 0; i < N_Bodies; i++) {
		totalKE += 0.5 * mass[i] * vel_sq[i];
		selfPE  += body_pot[i];

		// External potential energy of this body. Halo: matches the force model --
		// post-warmup a body feels every system's halo, but during warmup systems
		// are isolated and feel only their OWN halo (cross-halo terms omitted).
		// FDE: a = FDE*pos => Phi = -0.5*FDE*|pos|^2 about the origin.
		double phi = 0.0;
		if (warmupActive) {
			int s = body_system[i];
			if (halo_vc[s] != 0.0) {
				double *hc = &halo_center[s * 3];
				double d[3] = { pos[i][0]-hc[0], pos[i][1]-hc[1], pos[i][2]-hc[2] };
				phi += HaloPotential(s, vdot(d, d));
			}
		} else {
			for (int s = 0; s < N_Systems; s++) {
				if (halo_vc[s] == 0.0) continue;
				double *hc = &halo_center[s * 3];
				double d[3] = { pos[i][0]-hc[0], pos[i][1]-hc[1], pos[i][2]-hc[2] };
				phi += HaloPotential(s, vdot(d, d));
			}
		}
		if (FDE != 0.0)
			phi += -0.5 * FDE * vdot(pos[i], pos[i]);
		extPE += mass[i] * phi;

		M += mass[i];
		com[0] += mass[i] * pos[i][0];
		com[1] += mass[i] * pos[i][1];
		com[2] += mass[i] * pos[i][2];
		P[0] += mass[i] * vel[i][0];
		P[1] += mass[i] * vel[i][1];
		P[2] += mass[i] * vel[i][2];
	}
	selfPE *= 0.5;   // each pair counted twice in the body_pot summation
	totalPE = selfPE + extPE;
	totalE = totalKE + totalPE;

	totalP = vmag(P);

	// Center of mass position and velocity
	double Rcom[3] = {com[0]/M, com[1]/M, com[2]/M};
	double Vcom[3] = {P[0]/M, P[1]/M, P[2]/M};
	comSpeed = vmag(Vcom);

	// Second pass: angular momentum about the COM.
	// L = Sigma m_i * (r_i - Rcom) x v_i  (absolute velocities; the Vcom cross
	// term vanishes because Sigma m_i (r_i - Rcom) = 0)
	double L[3] = {0.0, 0.0, 0.0};
	for (int i = 0; i < N_Bodies; i++) {
		double r[3], l[3];
		vsub(pos[i], Rcom, r);
		vcross(r, vel[i], l);
		L[0] += mass[i] * l[0];
		L[1] += mass[i] * l[1];
		L[2] += mass[i] * l[2];
	}
	totalL = vmag(L);

	// Virial ratio 2K/|U| using the FULL potential energy (self-gravity + halo +
	// FDE). 1 = virialized equilibrium, 2 = marginally bound (E=0), >2 unbound.
	// Sentinel -1 when |U| is ~0 (no binding yet) so the display can show "---".
	double absU = fabs(totalPE);
	virialRatio = (absU > 0.0) ? (2.0 * totalKE / absU) : -1.0;
}

// Dissipative "sticky particle" gas collisions (Salo & Laurikainen 2000 sect 2.1),
// implemented as a kinetic Monte-Carlo (DSMC-style, Bird 1994) step -- conceptually
// the same scheme the paper uses. Rather than detecting instantaneous geometric
// overlaps (which in a thin disc are almost all VERTICAL, so they cool sigma_z but
// barely touch sigma_r), we bin the gas into collision cells and, per cell, select
// collision pairs STOCHASTICALLY at the physical rate 0.5*N*(N-1)*sigma*v_rel*dt/V
// (No-Time-Counter method), accepting each candidate with probability
// |v_rel|/v_rel_max. The line of centres for each accepted impact is SAMPLED from
// the hard-sphere distribution instead of read off instantaneous positions, which
// makes the cooling isotropic -- so sigma_r relaxes toward equilibrium along with
// sigma_z. On impact the line-of-centres velocity component is reversed and scaled
// by gas_alpha (0 = fully inelastic, the paper's value) via a momentum-conserving
// impulse. Operator-split after the gravity kick; modifies velocities only. Stars
// never collide. Runs during warmup too, so isolated gas reaches its cool
// equilibrium (sigma_gas ~ 5-10 km/s) before t=0.
void Simulation::ProcessGasCollisions()
{
	if (!gasEnabled || gas_radius <= 0.0 || gas_cell_size <= 0.0) return;

	const double PI = 3.14159265358979323846;
	const double L = gas_cell_size;
	const double inv_cell = 1.0 / L;
	const double sigma = PI * (2.0*gas_radius) * (2.0*gas_radius);   // pi*(2r)^2
	const double floorLen = 2.0 * gas_radius;   // min bounding-box edge per axis
	const long long B = 1LL << 20;
	auto cellKey = [B](long long ix, long long iy, long long iz) -> long long {
		return ((ix + B) & 0x1FFFFF) | (((iy + B) & 0x1FFFFF) << 21)
		     | (((iz + B) & 0x1FFFFF) << 42);
	};

	// Bin gas particles into collision cells.
	std::unordered_map<long long, std::vector<int>> grid;
	for (int i = 0; i < N_Bodies; i++) {
		if (!is_gas[i]) continue;
		double *p = pos[i];
		grid[cellKey((long long)floor(p[0]*inv_cell),
		             (long long)floor(p[1]*inv_cell),
		             (long long)floor(p[2]*inv_cell))].push_back(i);
	}
	if (grid.empty()) return;

	std::uniform_real_distribution<double> U01(0.0, 1.0);

	for (auto &kv : grid) {
		std::vector<int> &cell = kv.second;
		int Nc = (int)cell.size();
		if (Nc < 2) continue;

		// Cell mean velocity, an upper bound on pairwise |v_rel| (2 * max deviation),
		// and the bounding box of the occupied volume. Using the bounding box rather
		// than a fixed cubic L^3 lets the local density n = Nc/V reflect a thin disc:
		// a pancake of points has a small box, so n is not underestimated (which
		// would spuriously suppress the collision rate).
		double vm[3] = {0,0,0};
		double lo[3] = { 1e300, 1e300, 1e300 }, hi[3] = { -1e300, -1e300, -1e300 };
		for (int idx : cell)
			for (int k = 0; k < 3; k++) {
				vm[k] += vel[idx][k];
				if (pos[idx][k] < lo[k]) lo[k] = pos[idx][k];
				if (pos[idx][k] > hi[k]) hi[k] = pos[idx][k];
			}
		double invN = 1.0 / Nc;
		for (int k = 0; k < 3; k++) vm[k] *= invN;
		double maxdev = 0.0;
		for (int idx : cell) {
			double a0 = vel[idx][0]-vm[0], a1 = vel[idx][1]-vm[1], a2 = vel[idx][2]-vm[2];
			double d = sqrt(a0*a0 + a1*a1 + a2*a2);
			if (d > maxdev) maxdev = d;
		}
		double vrel_max = 2.0 * maxdev;
		if (vrel_max <= 0.0) continue;
		double Vc = std::max(hi[0]-lo[0], floorLen)
		          * std::max(hi[1]-lo[1], floorLen)
		          * std::max(hi[2]-lo[2], floorLen);

		// NTC candidate count, with the fractional part carried stochastically.
		double Mreal = 0.5 * (double)Nc * (Nc - 1) * sigma * vrel_max * dt / Vc;
		long Mcand = (long)Mreal;
		if (U01(gasRng) < (Mreal - (double)Mcand)) Mcand++;
		long Mmax = 4L * Nc;                        // safety cap against tiny Vc
		if (Mcand > Mmax) Mcand = Mmax;

		for (long c = 0; c < Mcand; c++) {
			int a = (int)(U01(gasRng) * Nc);
			int b2 = (int)(U01(gasRng) * (Nc - 1));
			if (b2 >= a) b2++;                      // distinct partner
			if (a  >= Nc) a  = Nc - 1;
			if (b2 >= Nc) b2 = Nc - 1;
			int i = cell[a], j = cell[b2];
			double *vi = vel[i], *vj = vel[j];
			double g[3] = { vi[0]-vj[0], vi[1]-vj[1], vi[2]-vj[2] };
			double grel = sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
			if (grel <= 0.0) continue;
			if (U01(gasRng) * vrel_max > grel) continue;   // NTC acceptance ~ |v_rel|

			// Hard-sphere line of centres n̂: angle theta from -ĝ has
			// p(theta) ~ sin(theta)cos(theta), so mu = cos(theta) = sqrt(U); azimuth
			// uniform. Built in a frame about ĝ, giving g·n̂ = -mu*grel < 0.
			double gh[3] = { g[0]/grel, g[1]/grel, g[2]/grel };
			double seed[3] = { 1.0, 0.0, 0.0 };
			if (fabs(gh[0]) > 0.9) { seed[0]=0.0; seed[1]=1.0; seed[2]=0.0; }
			double e1[3] = { seed[1]*gh[2]-seed[2]*gh[1],
			                 seed[2]*gh[0]-seed[0]*gh[2],
			                 seed[0]*gh[1]-seed[1]*gh[0] };
			double e1n = sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]);
			e1[0]/=e1n; e1[1]/=e1n; e1[2]/=e1n;
			double e2[3] = { gh[1]*e1[2]-gh[2]*e1[1],
			                 gh[2]*e1[0]-gh[0]*e1[2],
			                 gh[0]*e1[1]-gh[1]*e1[0] };
			double mu = sqrt(U01(gasRng));
			double st = sqrt(std::max(0.0, 1.0 - mu*mu));
			double eps = 2.0 * PI * U01(gasRng);
			double ce = cos(eps), se = sin(eps);
			double nrm[3];
			for (int k = 0; k < 3; k++)
				nrm[k] = -mu*gh[k] + st*(ce*e1[k] + se*e2[k]);

			double dvn = g[0]*nrm[0] + g[1]*nrm[1] + g[2]*nrm[2];   // = -mu*grel < 0
			double mi = mass[i], mj = mass[j];
			double mred = (mi*mj) / (mi + mj);
			double Jn = -(1.0 + gas_alpha) * dvn * mred;
			double ji = Jn/mi, jj = Jn/mj;
			vi[0] += ji*nrm[0]; vi[1] += ji*nrm[1]; vi[2] += ji*nrm[2];
			vj[0] -= jj*nrm[0]; vj[1] -= jj*nrm[1]; vj[2] -= jj*nrm[2];
		}
	}
}

// TEMPORARY bar/arm diagnostic. For each system, measure the m = 1..6 Fourier
// amplitudes A_m(R) and phases of the disc's mass distribution vs cylindrical
// radius. A strong, radially-coherent m=2 in the INNER disc (near-constant phase)
// is a bar; whichever m dominates the OUTER disc is the arm number (m=2 two-arm,
// m=3 three-arm, m=4 four-arm, ...), and a phase that winds with radius marks a
// spiral. The disc plane is found from the system's own mass-weighted angular
// momentum, so no orientation input is needed. Written to bar_diagnostic.csv. To
// be removed with the diagnostic cleanup.
void Simulation::LogBarDiagnostic() {
	if (barLogEvery <= 0) return;
	if ((barStepCount++ % barLogEvery) != 0) return;

	const int NR = 20;                 // radial bins per disc
	const int MMAX = 6;                // Fourier orders m = 1..MMAX logged
	if (!barLog) {
		barLog = fopen("bar_diagnostic.csv", "w");
		if (!barLog) { barLogEvery = 0; return; }
		fprintf(barLog, "t,system,bin,r_lo,r_hi,n,"
		        "A1,A2,A3,A4,A5,A6,ph1,ph2,ph3,ph4,ph5,ph6\n");
	}

	int base = 0;
	for (int sys = 0; sys < N_Systems; sys++) {
		int Nsys = N_System_Bodies[sys];
		if (Nsys < 2) { base += Nsys; continue; }

		// Mass-weighted centre, mean velocity, and angular momentum -> disc normal.
		double c[3] = {0,0,0}, vc[3] = {0,0,0}, mtot = 0.0;
		for (int i = 0; i < Nsys; i++) {
			int gi = base + i; double mi = mass[gi];
			for (int k = 0; k < 3; k++) { c[k] += mi*pos[gi][k]; vc[k] += mi*vel[gi][k]; }
			mtot += mi;
		}
		if (mtot <= 0.0) { base += Nsys; continue; }
		for (int k = 0; k < 3; k++) { c[k] /= mtot; vc[k] /= mtot; }
		double L[3] = {0,0,0};
		for (int i = 0; i < Nsys; i++) {
			int gi = base + i; double mi = mass[gi];
			double dx = pos[gi][0]-c[0], dy = pos[gi][1]-c[1], dz = pos[gi][2]-c[2];
			double dvx = vel[gi][0]-vc[0], dvy = vel[gi][1]-vc[1], dvz = vel[gi][2]-vc[2];
			L[0] += mi*(dy*dvz - dz*dvy);
			L[1] += mi*(dz*dvx - dx*dvz);
			L[2] += mi*(dx*dvy - dy*dvx);
		}
		double Ln = sqrt(L[0]*L[0]+L[1]*L[1]+L[2]*L[2]);
		double nrm[3] = { 0.0, 1.0, 0.0 };
		if (Ln > 0.0) { nrm[0]=L[0]/Ln; nrm[1]=L[1]/Ln; nrm[2]=L[2]/Ln; }

		// In-plane basis (u,w) for the azimuth phi = atan2(ip.w, ip.u).
		double seed[3] = {1.0, 0.0, 0.0};
		if (fabs(nrm[0]) > 0.9) { seed[0]=0.0; seed[1]=1.0; seed[2]=0.0; }
		double u[3] = { nrm[1]*seed[2]-nrm[2]*seed[1],
		                nrm[2]*seed[0]-nrm[0]*seed[2],
		                nrm[0]*seed[1]-nrm[1]*seed[0] };
		double un = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]); u[0]/=un; u[1]/=un; u[2]/=un;
		double w[3] = { nrm[1]*u[2]-nrm[2]*u[1],
		                nrm[2]*u[0]-nrm[0]*u[2],
		                nrm[0]*u[1]-nrm[1]*u[0] };

		// Cylindrical radius of every particle; bin linearly over [0, Rmax].
		double Rmax = 0.0;
		std::vector<double> Rp(Nsys), phip(Nsys);
		for (int i = 0; i < Nsys; i++) {
			int gi = base + i;
			double dp[3] = { pos[gi][0]-c[0], pos[gi][1]-c[1], pos[gi][2]-c[2] };
			double z = dp[0]*nrm[0]+dp[1]*nrm[1]+dp[2]*nrm[2];
			double ipx = dp[0]-z*nrm[0], ipy = dp[1]-z*nrm[1], ipz = dp[2]-z*nrm[2];
			double R = sqrt(ipx*ipx+ipy*ipy+ipz*ipz);
			Rp[i] = R;
			phip[i] = atan2(ipx*w[0]+ipy*w[1]+ipz*w[2], ipx*u[0]+ipy*u[1]+ipz*u[2]);
			if (R > Rmax) Rmax = R;
		}
		if (Rmax <= 0.0) { base += Nsys; continue; }
		double dR = Rmax / NR;

		// Per bin: mass, count, and the m=1..MMAX Fourier sums of the mass
		// distribution. A_m = |sum m e^{i m phi}| / sum m ; the m that dominates
		// the OUTER disc is the arm number (m=2 two-arm, m=3 three-arm, ...), while
		// a strong radially-coherent m=2 in the INNER disc is the bar.
		std::vector<double> sM(NR,0.0);
		std::vector<long> cnt(NR,0);
		std::vector<double> sC(NR*(MMAX+1),0.0), sS(NR*(MMAX+1),0.0);
		for (int i = 0; i < Nsys; i++) {
			int b = (int)(Rp[i] / dR); if (b >= NR) b = NR-1; if (b < 0) b = 0;
			double mi = mass[base + i];
			sM[b] += mi;
			cnt[b]++;
			// cos(m phi)/sin(m phi) built by Chebyshev recurrence from m=1.
			double cph = cos(phip[i]), sph = sin(phip[i]);
			double cm = cph, sm = sph;               // m = 1
			for (int m = 1; m <= MMAX; m++) {
				sC[b*(MMAX+1)+m] += mi*cm;
				sS[b*(MMAX+1)+m] += mi*sm;
				double cn = cm*cph - sm*sph;         // (m+1) angle addition
				double sn = sm*cph + cm*sph;
				cm = cn; sm = sn;
			}
		}
		for (int b = 0; b < NR; b++) {
			if (sM[b] <= 0.0) continue;
			double A[MMAX+1], ph[MMAX+1];
			for (int m = 1; m <= MMAX; m++) {
				double C = sC[b*(MMAX+1)+m], S = sS[b*(MMAX+1)+m];
				A[m] = sqrt(C*C + S*S) / sM[b];
				ph[m] = atan2(S, C) / m;             // pattern orientation for order m
			}
			fprintf(barLog,
			        "%.6f,%d,%d,%.3f,%.3f,%ld,"
			        "%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,"
			        "%.5f,%.5f,%.5f,%.5f,%.5f,%.5f\n",
			        t, sys, b, b*dR, (b+1)*dR, cnt[b],
			        A[1],A[2],A[3],A[4],A[5],A[6],
			        ph[1],ph[2],ph[3],ph[4],ph[5],ph[6]);
		}
		base += Nsys;
	}
	fflush(barLog);
}

void Simulation::ZeroNetMomentum(int system) {
	// Remove one system's net momentum. The procedural generators sample particle
	// azimuths randomly and never correct the sum, leaving a residual of order
	// v_c/sqrt(N). Nothing damps it, so it carries the system off the origin at
	// constant velocity forever -- for a 100k-particle disc at 220 km/s that is
	// ~1 code unit per time unit.
	//
	// Called from the generators themselves, so it applies exactly to the systems
	// that need it. Systems built from explicit `Body` state vectors never reach
	// here, which is correct: their net momentum is physical (real ephemerides,
	// where the Sun's recoil balances the planets) and zeroing it would corrupt
	// the trajectories.
	//
	// PRECONDITION: the system's particles must not yet carry the bulk velocity.
	// Otherwise this would subtract the intended bulk motion along with the noise.
	//
	// Per system rather than globally, because a global correction cancels only
	// the total: each individual system would keep its own residual drift, and
	// under warmup the isolated systems would wander away from their specified
	// separation before t = 0.
	int sysIdx = 0;
	for (int i = 0; i < system; i++) sysIdx += N_System_Bodies[i];
	int n = N_System_Bodies[system];

	double cvx = 0.0, cvy = 0.0, cvz = 0.0, total_m = 0.0;
	for (int i = 0; i < n; i++) {
		double mi = mass[sysIdx + i];
		cvx += mi * vel[sysIdx + i][0];
		cvy += mi * vel[sysIdx + i][1];
		cvz += mi * vel[sysIdx + i][2];
		total_m += mi;
	}
	if (total_m <= 0.0) return;

	double inv_m = 1.0 / total_m;
	double v0[3] = {cvx*inv_m, cvy*inv_m, cvz*inv_m};
	for (int i = 0; i < n; i++) {
		vel[sysIdx + i][0] -= v0[0];
		vel[sysIdx + i][1] -= v0[1];
		vel[sysIdx + i][2] -= v0[2];
	}
}

void Simulation::ApplyBulkVelocity(int system) {
	// Add the system's stored bulk velocity to every one of its particles.
	// A uniform shift leaves all internal relative motion untouched, so the disc
	// structure is preserved exactly; only the system as a whole starts moving.
	int sysIdx = 0;
	for (int i = 0; i < system; i++) sysIdx += N_System_Bodies[i];
	double *vb = &system_bulk_vel[system * 3];
	for (int i = 0; i < N_System_Bodies[system]; i++) {
		vel[sysIdx + i][0] += vb[0];
		vel[sysIdx + i][1] += vb[1];
		vel[sysIdx + i][2] += vb[2];
	}
	// The inertial halo centre shares the bulk motion, so it starts at t=0 moving
	// with its galaxy (its velocity was withheld during warmup, like the particles').
	halo_vel[system*3+0] += vb[0];
	halo_vel[system*3+1] += vb[1];
	halo_vel[system*3+2] += vb[2];
}

void Simulation::ApplyBulkVelocities() {
	// End of warmup: hand every system the bulk velocity that was withheld at
	// load time.
	for (int sys = 0; sys < N_Systems; sys++)
		ApplyBulkVelocity(sys);
	bulkVelocityApplied = true;

	std::cout << "t=0: warmup complete, applying bulk velocities and enabling"
	          << " cross-system gravity" << std::endl;
}

void Simulation::BuildOctreeForSystem(int sys, int &outFirst, int &outCount) {
	// Build the tree from ONE system's gravitating bodies only.
	//
	// This is the reliable way to isolate the systems. Spatial separation is NOT
	// sufficient with a shared tree: the root node spans every body, and the
	// Barnes-Hut acceptance test (d*d <= theta_sq * dsq in BHTree.cpp) gets
	// EASIER to satisfy as separation grows, so a distant system is more likely
	// to be swallowed as a single accepted multipole, not less. Excluding its
	// bodies from the tree entirely is a structural guarantee that needs no
	// assumption about distance or opening angle.
	int sysStart = 0;
	for (int s = 0; s < sys; s++) sysStart += N_System_Bodies[s];

	int numActive = 0;
	for (int i = 0; i < N_System_Bodies[sys]; i++) {
		int bi = sysStart + i;
		if (!has_gravity[bi]) continue;
		pos_f[bi*3]   = (float)pos[bi][0];
		pos_f[bi*3+1] = (float)pos[bi][1];
		pos_f[bi*3+2] = (float)pos[bi][2];
		sortedIdx[numActive++] = bi;
	}

	outFirst = 0;
	outCount = numActive;
	if (numActive == 0) return;

	float p_min[3], p_max[3];
	int first = sortedIdx[0];
	p_min[0] = pos_f[first*3]; p_min[1] = pos_f[first*3+1]; p_min[2] = pos_f[first*3+2];
	p_max[0] = p_min[0]; p_max[1] = p_min[1]; p_max[2] = p_min[2];
	for (int i = 1; i < numActive; i++) {
		int bi = sortedIdx[i];
		for (int j = 0; j < 3; j++) {
			if (pos_f[bi*3+j] < p_min[j]) p_min[j] = pos_f[bi*3+j];
			if (pos_f[bi*3+j] > p_max[j]) p_max[j] = pos_f[bi*3+j];
		}
	}

	float maxDim = p_max[0] - p_min[0];
	if (p_max[1] - p_min[1] > maxDim) maxDim = p_max[1] - p_min[1];
	if (p_max[2] - p_min[2] > maxDim) maxDim = p_max[2] - p_min[2];
	maxDim *= 1.001f;
	float center[3] = {(p_min[0]+p_max[0])*0.5f, (p_min[1]+p_max[1])*0.5f,
	                   (p_min[2]+p_max[2])*0.5f};
	for (int j = 0; j < 3; j++) {
		p_min[j] = center[j] - maxDim*0.5f;
		p_max[j] = center[j] + maxDim*0.5f;
	}

	// Morton sort for cache-coherent insertion, same as BuildOctree()
	for (int i = 0; i < numActive; i++) {
		int bi = sortedIdx[i];
		mortonCodes[bi] = morton3D(pos_f[bi*3], pos_f[bi*3+1], pos_f[bi*3+2],
		                           p_min, p_max);
	}
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
	for (int i = 0; i < numActive; i++) {
		int bi = sortedIdx[i];
		Octree.InsertBody(pos_f.data() + bi*3, (float)mass[bi], bi);
	}
	Octree.CalcMasses();
}

void Simulation::CalcAccelIsolated() {
	// Warmup force evaluation: each system feels ONLY its own gravity and its
	// OWN halo. No cross-system particle forces, no cross-system halo terms.
	ZeroAccelerationRange(0, N_Bodies-1);
	ComputeHaloCenters();

	int sysStart = 0;
	for (int sys = 0; sys < N_Systems; sys++) {
		int n = N_System_Bodies[sys];

		if (Gravity_Oct) {
			int first, count;
			BuildOctreeForSystem(sys, first, count);
			// sortedIdx[0..count) now holds only this system's bodies, and the
			// tree contains only those bodies, so this is exactly self-gravity.
			// Threaded the same way as the normal path: the tree is read-only
			// here and each body writes only its own acc[], so this is safe.
			if (multiThreading && count >= numThreads && numThreads > 1) {
				int chunk = count / numThreads;
				for (int th = 0; th < numThreads; ++th) {
					int kStart = th * chunk;
					int kEnd = (th == numThreads - 1) ? (count - 1) : (kStart + chunk - 1);
					pool->submit([this, kStart, kEnd]() {
						for (int k = kStart; k <= kEnd; k++) {
							int bi = sortedIdx[k];
							float pf[3] = {(float)pos[bi][0], (float)pos[bi][1],
							               (float)pos[bi][2]};
							double a[3], pot;
							Octree.CalcAcceleration(pf, bi, (float)G, (float)SoftSq(bi),
							        (float)(BH_Opening_Theta*BH_Opening_Theta), a, &pot);
							vadd(acc[bi], a, acc[bi]);
							body_pot[bi] = pot * mass[bi];
							vscaleadd(pos[bi], FDE, acc[bi]);
						}
					});
				}
				pool->waitAll();
			} else {
				for (int k = 0; k < count; k++) {
					int bi = sortedIdx[k];
					float pf[3] = {(float)pos[bi][0], (float)pos[bi][1], (float)pos[bi][2]};
					double a[3], pot;
					Octree.CalcAcceleration(pf, bi, (float)G, (float)SoftSq(bi),
					                        (float)(BH_Opening_Theta*BH_Opening_Theta),
					                        a, &pot);
					vadd(acc[bi], a, acc[bi]);
					body_pot[bi] = pot * mass[bi];
					vscaleadd(pos[bi], FDE, acc[bi]);
				}
			}
		} else {
			// P2P restricted to this system's index range
			for (int i = sysStart; i < sysStart + n; i++) {
				vscaleadd(pos[i], FDE, acc[i]);
				double pot = 0.0, d[3];
				double softI = SoftSq(i);      // sink-based: constant over the j-loop
				for (int j = sysStart; j < sysStart + n; j++) {
					if (j == i) continue;
					vsub(pos[j], pos[i], d);
					double dsq = vmagsqsoft(d, softI);
					double r_inv = 1.0 / sqrt(dsq);
					vscaleadd(d, G * mass[j] * (r_inv/dsq), acc[i]);
					pot += -G * mass[j] * r_inv;
				}
				body_pot[i] = pot * mass[i];
			}
		}

		// Own halo only -- the cross-halo loop is deliberately omitted.
		if (halo_vc[sys] > 0.0) {
			double *hc = &halo_center[sys * 3];
			for (int i = sysStart; i < sysStart + n; i++) {
				double r_halo[3];
				r_halo[0] = hc[0] - pos[i][0];
				r_halo[1] = hc[1] - pos[i][1];
				r_halo[2] = hc[2] - pos[i][2];
				double rsq = vmagsq(r_halo);
				if (rsq > 1e-10) {
					double s = HaloScale(sys, rsq);
					vscaleadd(r_halo, s, acc[i]);
				}
			}
		}

		sysStart += n;
	}

	// Halo monopole removal, applied PER SYSTEM during warmup. The normal
	// global version would transfer spurious momentum between isolated systems.
	if (Remove_Halo_Monopole) {
		sysStart = 0;
		for (int sys = 0; sys < N_Systems; sys++) {
			int n = N_System_Bodies[sys];
			if (halo_vc[sys] > 0.0) {
				double fh[3] = {0.0, 0.0, 0.0};
				double total_m = 0.0;
				double *hc = &halo_center[sys * 3];
				for (int i = sysStart; i < sysStart + n; i++) {
					if (!has_gravity[i]) continue;
					double mi = mass[i];
					total_m += mi;
					double r_halo[3];
					r_halo[0] = hc[0] - pos[i][0];
					r_halo[1] = hc[1] - pos[i][1];
					r_halo[2] = hc[2] - pos[i][2];
					double rsq = vmagsq(r_halo);
					if (rsq > 1e-10) {
						double s = HaloScale(sys, rsq);
						fh[0] += mi * s * r_halo[0];
						fh[1] += mi * s * r_halo[1];
						fh[2] += mi * s * r_halo[2];
					}
				}
				if (total_m > 0.0) {
					double inv_m = 1.0 / total_m;
					double aCorr[3];
					vscale(fh, inv_m, aCorr);
					for (int i = sysStart; i < sysStart + n; i++) {
						if (!has_gravity[i]) continue;
						acc[i][0] -= aCorr[0];
						acc[i][1] -= aCorr[1];
						acc[i][2] -= aCorr[2];
					}
				}
			}
			sysStart += n;
		}
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

	// Drift the inertial halo centres with the same velocity-Verlet step. Held
	// static during warmup (systems isolated and pinned in place); the barycentre
	// tracking in CalcAccelIsolated keeps halo_center on the relaxing disc until
	// t=0, after which it evolves purely under gravity.
	if (!warmupActive) {
		for (int s = 0; s < N_Systems; s++) {
			if (halo_mass[s] <= 0.0) continue;
			for (int k = 0; k < 3; k++) {
				halo_acc_prev[s*3+k] = halo_acc[s*3+k];
				halo_center[s*3+k] += halo_vel[s*3+k]*dt + 0.5*halo_acc[s*3+k]*dt*dt;
			}
		}
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

	// Kick the inertial halo centres (velocity-Verlet second half).
	if (!warmupActive) {
		for (int s = 0; s < N_Systems; s++) {
			if (halo_mass[s] <= 0.0) continue;
			for (int k = 0; k < 3; k++)
				halo_vel[s*3+k] += 0.5*(halo_acc[s*3+k] + halo_acc_prev[s*3+k])*dt;
		}
	}
}

void Simulation::Step()
{
	CalcLeapFrogPositions();

	// End of warmup. Checked BEFORE the force evaluation so the first step at
	// t >= 0 already sees the coupled system and the applied bulk velocities.
	if (warmupActive && t >= 0.0) {
		warmupActive = false;
		if (!bulkVelocityApplied) ApplyBulkVelocities();
	}

	if (warmupActive) {
		// Systems isolated: each builds its own tree and feels only its own
		// gravity and halo. CalcAccelIsolated does its own tree construction
		// per system, so the shared BuildOctree() below is skipped.
		CalcAccelIsolated();
	} else {
		// The octree must be rebuilt from the current positions before forces are
		// evaluated against it, otherwise every step computes forces from the
		// previous step's tree -- a body is attracted toward where mass used to be.
		// The constructor already initialized in this order; Step() did not, so this
		// aligns the two. Correctness fix rather than a drift fix: it cut the net
		// spurious gravity force ~30x when measured, but that error is largely
		// direction-randomizing and cancels in the sum, so its contribution to
		// center-of-mass drift was small (a ~12% change in COM velocity growth).
		if (Gravity_Oct)
		{
			BuildOctree();
		}

		CalcAccelerations();
	}

	CalcLeapFrogVelocitiesAndOutputs();

	// Operator-split gas dissipation: after the gravity velocity kick, resolve
	// inelastic gas-gas collisions (no-op unless a system declared gas particles).
	// Runs during warmup too, so an isolated gas disc relaxes to its cool
	// equilibrium (Salo & Laurikainen: sigma_gas ~ 5-10 km/s) before t = 0.
	ProcessGasCollisions();

	CalcSystemQuantities();

	LogBarDiagnostic();

	if (Data_Log)
	{
		fwrite(pos_sq.data(), N_Bodies*sizeof(double), 1, DataLog);
		fwrite(vel_sq.data(), N_Bodies*sizeof(double), 1, DataLog);
		fwrite(acc_sq.data(), N_Bodies*sizeof(double), 1, DataLog);
	}

	t += dt;

	// Retarget the camera after positions have advanced, so the followed body is
	// centred using its current position rather than last frame's.
	UpdateCameraFollow();
}

void Simulation::UpdateCameraFollow()
{
	// Retarget the look-at point onto the followed system's centre (its inertial
	// halo centre, or its centre of mass if the system has no halo).
	//
	// The camera position is moved by the SAME delta as the look-at point, which
	// keeps the relative offset vector (Cam.pos - Cam.lookAt) exactly unchanged.
	// That offset is what defines the spherical coordinates phi, theta and r, so
	// the viewing angle and zoom are untouched and the W/A/S/D/J/L controls and
	// Camera_Orbit continue to work against the moving target.
	//
	// Recomputing the position from stored angles instead would fight the user's
	// input and accumulate drift through the acos/atan2 round trip.
	if (CamFollowSystem < 0) return;

	int sys = CamFollowSystem;
	double target[3];
	if (halo_mass[sys] > 0.0) {
		// Follow the inertial halo centre, not the central-body particle. The halo
		// centre is a smooth dynamical point -- its shot noise is damped by the
		// large halo mass and by double integration -- and it stays locked on the
		// galaxy core rather than drifting with tidal debris. Following the single
		// (now free, low-mass) central particle instead made the view jitter.
		target[0] = halo_center[sys*3+0];
		target[1] = halo_center[sys*3+1];
		target[2] = halo_center[sys*3+2];
	} else {
		// No halo: fall back to the system's mass-weighted centre of mass.
		int sysIdx = 0;
		for (int s = 0; s < sys; s++) sysIdx += N_System_Bodies[s];
		double cx = 0.0, cy = 0.0, cz = 0.0, m_tot = 0.0;
		for (int i = 0; i < N_System_Bodies[sys]; i++) {
			double mi = mass[sysIdx + i];
			cx += mi * pos[sysIdx + i][0];
			cy += mi * pos[sysIdx + i][1];
			cz += mi * pos[sysIdx + i][2];
			m_tot += mi;
		}
		double inv = (m_tot > 0.0) ? 1.0 / m_tot : 0.0;
		target[0] = cx * inv; target[1] = cy * inv; target[2] = cz * inv;
	}

	double delta[3];
	vsub(target, Cam.lookAt, delta);

	vadd(Cam.lookAt, delta, Cam.lookAt);
	vadd(Cam.pos, delta, Cam.pos);
}

void Simulation::CamMove(double d_phi, double d_theta, double d_r)
{
	double rel[3];
	vsub(Cam.pos, Cam.lookAt, rel);

	double r = vmag(rel);
	double phi = acos(rel[1]/r);
	double theta = atan2(rel[2],rel[0]);

	if (sin(phi) < 1e-6)
		theta = -(Cam.theta) + M_PI/2;

	phi += d_phi;
	if (phi < M_PI/180) phi = M_PI/180;
	if (phi > M_PI*(1.0-1.0/180)) phi = M_PI*(1.0-1.0/180);
	theta += d_theta;
	r += d_r;

	vset(
			Cam.lookAt[0] + r*sin(phi)*cos(theta),
			Cam.lookAt[1] + r*cos(phi),
			Cam.lookAt[2] + r*sin(phi)*sin(theta),
			Cam.pos
		);

	Cam.phi = phi - M_PI/2;
	Cam.theta = -(theta - M_PI/2);
}

void Simulation::CamShift(double dx, double dy, double dz)
{
	Cam.lookAt[0] += dx;
	Cam.lookAt[1] += dy;
	Cam.lookAt[2] += dz;
	Cam.pos[0] += dx;
	Cam.pos[1] += dy;
	Cam.pos[2] += dz;
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
	hdrFBO = 0;
	hdrColorTex = 0;
	hdrDepthRBO = 0;

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
		out vec2 vLocal;
		uniform mat4 uVP;
		uniform vec3 uRight;
		uniform vec3 uUp;
		uniform float uSize;
		void main() {
			// aOffset is a unit quad in [-0.5, 0.5]; uSize scales the billboard to
			// its world size, while vLocal (kept unit) drives a size-independent
			// Gaussian in the fragment shader.
			vec2 off = aOffset * uSize;
			vec3 worldPos = aPos + off.x * uRight + off.y * uUp;
			gl_Position = uVP * vec4(worldPos, 1.0);
			fColor = aColor;
			vLocal = aOffset;
		}
	)";

	const char *particleFrag = R"(
		#version 330 core
		in vec4 fColor;
		in vec2 vLocal;
		out vec4 FragColor;
		void main() {
			// Radial Gaussian falloff turns the flat billboard into a soft round
			// glow, so lone particles read as faint stars and overlaps blend
			// smoothly. vLocal is in [-0.5, 0.5]; the steep constant makes the
			// glow vanish (~0.02) by the quad edge so no straight edges show, which
			// is why a single quad suffices (the old 8-point star is no longer
			// needed). Additive blending (SRC_ALPHA, ONE) accumulates the
			// alpha-weighted colour into the HDR buffer; the tone-map pass
			// compresses it later.
			float r2 = dot(vLocal, vLocal);
			float fall = exp(-16.0 * r2);
			FragColor = vec4(fColor.rgb, fColor.a * fall);
		}
	)";

	particleShader = CompileShader(particleVert, particleFrag);

	// Unit quad (2 triangles = 6 vertices), corners at +/-0.5. The fragment
	// shader's Gaussian falloff shapes it into a round glow, so a single quad is
	// enough; uSize (set at draw time) controls its world size.
	const float d = 0.5f;
	float shapeData[12];
	shapeData[0]=-d; shapeData[1]=-d;
	shapeData[2]=-d; shapeData[3]= d;
	shapeData[4]= d; shapeData[5]= d;
	shapeData[6]=-d; shapeData[7]=-d;
	shapeData[8]= d; shapeData[9]= d;
	shapeData[10]=d; shapeData[11]=-d;

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

	// --- Tone-map post-process shader ---
	// Fullscreen triangle generated from gl_VertexID (no vertex buffer needed).
	// Reads the HDR accumulation texture and compresses it into displayable [0,1]
	// with an exponential curve: bright cores roll off toward white instead of
	// hard-clipping, while faint lone particles are lifted into visibility.
	const char *tonemapVert = R"(
		#version 330 core
		out vec2 vUV;
		void main() {
			vec2 p = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
			vUV = p;
			gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
		}
	)";
	const char *tonemapFrag = R"(
		#version 330 core
		in vec2 vUV;
		out vec4 FragColor;
		uniform sampler2D uHDR;
		uniform float uExposure;
		void main() {
			vec3 hdr = texture(uHDR, vUV).rgb;
			vec3 mapped = 1.0 - exp(-uExposure * hdr);
			FragColor = vec4(mapped, 1.0);
		}
	)";
	tonemapShader = CompileShader(tonemapVert, tonemapFrag);
	glGenVertexArrays(1, &tonemapVAO);
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

void Simulation::CreateHDRFBO(int width, int height)
{
	if (hdrFBO) {
		glDeleteFramebuffers(1, &hdrFBO);
		glDeleteTextures(1, &hdrColorTex);
		glDeleteRenderbuffers(1, &hdrDepthRBO);
	}

	glGenFramebuffers(1, &hdrFBO);
	glBindFramebuffer(GL_FRAMEBUFFER, hdrFBO);

	// RGBA16F so additive accumulation of overlapping particles can exceed 1.0
	// without clipping; the tone-map pass compresses this range for display.
	glGenTextures(1, &hdrColorTex);
	glBindTexture(GL_TEXTURE_2D, hdrColorTex);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, width, height, 0, GL_RGBA, GL_FLOAT, nullptr);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, hdrColorTex, 0);

	glGenRenderbuffers(1, &hdrDepthRBO);
	glBindRenderbuffer(GL_RENDERBUFFER, hdrDepthRBO);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, width, height);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, hdrDepthRBO);

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void Simulation::ReSizeGL(int width, int height)
{
	if (height == 0) height = 1;
	winWidth = width;
	winHeight = height;
	glViewport(0, 0, width, height);

	// The HDR buffer is always needed: particles render into it every frame
	// regardless of whether we are recording.
	CreateHDRFBO(width, height);

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
	// Particles are additively accumulated into the HDR buffer first, then a
	// tone-map pass compresses that into the LDR final target (screen or record
	// FBO). This keeps on-screen and recorded output identical.
	glBindFramebuffer(GL_FRAMEBUFFER, hdrFBO);
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
	// The per-particle alpha is an accumulation weight, not a coverage value:
	// additive blending sums color*alpha into the HDR buffer and ToneMapExposure
	// controls overall brightness. Only the ratio between classes matters here --
	// system-centre markers (5.0) accumulate 5x faster than ordinary stars (1.0).
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
			clrBuf[i*4+0] = 0.0f; clrBuf[i*4+1] = 1.0f; clrBuf[i*4+2] = 0.0f; clrBuf[i*4+3] = 5.0f;
		} else {
			float r = (float)cbrt(acc_sq[i] / accel_sq_color_thresh);
			float b = 1.0f - r;
			r = (r < 0.3f) ? 0.3f : r;
			b = (b < 0.3f) ? 0.3f : b;
			clrBuf[i*4+0] = r; clrBuf[i*4+1] = 0.3f; clrBuf[i*4+2] = b; clrBuf[i*4+3] = 1.0f;
		}
	}

	glUseProgram(particleShader);
	glUniformMatrix4fv(glGetUniformLocation(particleShader, "uVP"), 1, GL_FALSE, vp);
	glUniform3fv(glGetUniformLocation(particleShader, "uRight"), 1, right);
	glUniform3fv(glGetUniformLocation(particleShader, "uUp"), 1, up);
	// World size of a particle billboard. Chosen with the frag-shader Gaussian
	// (exp(-16 r^2)) so the visible glow roughly matches the old 8-point star while
	// staying fully contained inside the single quad (no clipped straight edges).
	glUniform1f(glGetUniformLocation(particleShader, "uSize"), 1.8f);

	glBindVertexArray(particleVAO);
	glBindBuffer(GL_ARRAY_BUFFER, particlePosVBO);
	glBufferData(GL_ARRAY_BUFFER, N_Bodies * 3 * sizeof(float), posBuf, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, particleColorVBO);
	glBufferData(GL_ARRAY_BUFFER, N_Bodies * 4 * sizeof(float), clrBuf, GL_DYNAMIC_DRAW);
	glDrawArraysInstanced(GL_TRIANGLES, 0, 6, N_Bodies);
	glBindVertexArray(0);

	// --- Tone-map the HDR scene into the LDR final target ---
	// Final target is the record FBO when recording, else the default framebuffer.
	// Overlays (octree wireframe, HUD text) are drawn afterwards in LDR on top of
	// the tone-mapped image so they are not themselves compressed/dimmed.
	GLuint finalFBO = Record_Video ? recordFBO : 0;
	glBindFramebuffer(GL_FRAMEBUFFER, finalFBO);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glDisable(GL_BLEND);   // tone-map fully replaces the target
	glUseProgram(tonemapShader);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, hdrColorTex);
	glUniform1i(glGetUniformLocation(tonemapShader, "uHDR"), 0);
	glUniform1f(glGetUniformLocation(tonemapShader, "uExposure"), (float)ToneMapExposure);
	glBindVertexArray(tonemapVAO);
	glDrawArrays(GL_TRIANGLES, 0, 3);
	glBindVertexArray(0);
	glEnable(GL_BLEND);    // restore additive blend for the overlays below

	// --- Draw octree wireframe ---
	if (DrawOctree) {
		if (Gravity_Oct) {
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
	}

	glUseProgram(0);
	glFlush();
}

void Simulation::DrawInfo(double fps)
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

	const int NUM_LINES = 9;
	const char *labels[NUM_LINES] = {"FPS", "T", "KE", "PE", "E", "VIR", "P", "L", "VCOM"};
	double vals[NUM_LINES] = {
		(double)(int)(fps + 0.5), t, totalKE, totalPE, totalE,
		(virialRatio < 0.0 ? 0.0 : virialRatio), totalP, totalL, comSpeed
	};

	// Numeric portion of each line (sign, if any, is part of these strings)
	char nums[NUM_LINES][32];
	snprintf(nums[0], sizeof(nums[0]), "%d",   (int)(fps + 0.5));
	snprintf(nums[1], sizeof(nums[1]), "%.4F", t);
	snprintf(nums[2], sizeof(nums[2]), "%.4E", totalKE);
	snprintf(nums[3], sizeof(nums[3]), "%.4E", totalPE);
	snprintf(nums[4], sizeof(nums[4]), "%.4E", totalE);
	if (virialRatio < 0.0)
		snprintf(nums[5], sizeof(nums[5]), "---");
	else
		snprintf(nums[5], sizeof(nums[5]), "%.4F", virialRatio);
	snprintf(nums[6], sizeof(nums[6]), "%.4E", totalP);
	snprintf(nums[7], sizeof(nums[7]), "%.4E", totalL);
	snprintf(nums[8], sizeof(nums[8]), "%.4E", comSpeed);

	// Longest label sets the right-alignment width for all labels
	int maxLabel = 0;
	for (int l = 0; l < NUM_LINES; l++) {
		int len = (int)strlen(labels[l]);
		if (len > maxLabel) maxLabel = len;
	}

	// Compose: right-aligned label, space, '=', then a space before the
	// number (omitted when the number is negative so the '-' takes that slot)
	char lines[NUM_LINES][32];
	for (int l = 0; l < NUM_LINES; l++) {
		const char *pad = (vals[l] < 0.0) ? "" : " ";
		snprintf(lines[l], sizeof(lines[l]), "%*s =%s%s",
		         maxLabel, labels[l], pad, nums[l]);
	}

	int maxLen = 0;
	for (int l = 0; l < NUM_LINES; l++) {
		int len = (int)strlen(lines[l]);
		if (len > maxLen) maxLen = len;
	}

	std::vector<float> verts;
	float lineHeight = 10.0f;
	int scale = 1;

	for (int l = 0; l < NUM_LINES; l++) {
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
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(0, 0, winWidth, winHeight, GL_RGB, GL_UNSIGNED_BYTE, rgbOut);

	// Blit to default framebuffer for on-screen display
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	glBlitFramebuffer(0, 0, winWidth, winHeight, 0, 0, winWidth, winHeight,
					  GL_COLOR_BUFFER_BIT, GL_NEAREST);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	// Flip vertically (glReadPixels returns bottom-up)
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

void Simulation::BlitToScreen()
{
	glBindFramebuffer(GL_READ_FRAMEBUFFER, recordFBO);
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	glBlitFramebuffer(0, 0, winWidth, winHeight, 0, 0, winWidth, winHeight,
					  GL_COLOR_BUFFER_BIT, GL_NEAREST);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
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

	fwrite(pos_data.data(),		sizeof(double), N_Bodies*3,	StateFile);
	fwrite(vel_data.data(),		sizeof(double), N_Bodies*3,	StateFile);
	fwrite(mass.data(),			sizeof(double), N_Bodies,	StateFile);

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

	result = fread(pos_data.data(),		sizeof(double), N_Bodies*3,	StateFile);
	result = fread(vel_data.data(),		sizeof(double), N_Bodies*3,	StateFile);
	result = fread(mass.data(),			sizeof(double), N_Bodies,	StateFile);
	(void)result;

	fclose(StateFile);

	UpdateSofteningSq();   // r_soft may have changed on restore

	return true;
}
