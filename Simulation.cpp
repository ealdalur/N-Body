#include "Simulation.h"
#include <cstring>

Simulation::Simulation()
{
	multiThreading = true;
	int hardwareConcurrency = std::thread::hardware_concurrency();
	std::cout << "Hardware thread contexts: " << hardwareConcurrency << std::endl;
	numThreads = hardwareConcurrency - 2;
	if (numThreads < 4) numThreads = 4;
	std::cout << "Number of compute threads: " << numThreads << std::endl;

	pool = new ThreadPool(numThreads);

	posBuf = new float[N_BODIES * 3];
	clrBuf = new float[N_BODIES * 4];

	InitGL();

	int ki;
	for (int i=0; i<N_BODIES; i++)
	{
		ki = i*N_STATES;

		pos[i] = states+ki;
		vel[i] = states+ki+3;
		acc[i] = states_d[0]+ki+3;
		acc_prev[i] = states_d[1]+ki+3;

		has_gravity[i] = true;
	}

	if (!ReadState())
	{
		LoadDefaultState();
	}

	CalcOutputs();

	BuildOctree();
	DrawOctree = false;

	CalcDerivatives(states,states_d[0]);
	CalcOutputs();

	if (DATA_LOG)
	{
		DataLog = fopen("Simulation.dat", "wb");
		fwrite(&N_BODIES, sizeof(N_BODIES), 1, DataLog);
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

	if (DATA_LOG) fclose(DataLog);
}

void Simulation::LoadGalaxyDiscState(int system, double *sysPos, double *sysVel, double M, double Mfrac, double R, double Ri, double Vtol, double haloVc, double haloRc){

	double m,r,theta,phi,vm,m_orbit;
	double p[3],v[3];
	double y[3] = {0.0, 1.0, 0.0};

	int sysIdx = 0;
	for (int i=0; i<system; i++) sysIdx += N_SYSTEM_BODIES[i];

	halo_vc[system] = haloVc;
	halo_rc_sq[system] = haloRc * haloRc;
	halo_central[system] = sysIdx;

	for (int i=0; i<N_SYSTEM_BODIES[system]; i++)
		body_system[sysIdx+i] = system;

	m = Mfrac*M/(N_SYSTEM_BODIES[system]-1);
	mass[sysIdx] = M;
	has_gravity[sysIdx] = true;
	vcopy(sysPos, pos[sysIdx]);
	vcopy(sysVel, vel[sysIdx]);

	for (int i=1; i<N_SYSTEM_BODIES[system]; i++)
	{
		mass[sysIdx+i] = m;
		has_gravity[sysIdx+i] = true;
		
		r = (R-Ri)*drand() + Ri;
		
		phi = M_PI/2 + M_PI/64*(2*drand()-1);
		theta = 2*M_PI*drand();
		vset(r*sin(phi)*cos(theta),r*cos(phi),r*sin(phi)*sin(theta),p);

		m_orbit = M + m*(N_SYSTEM_BODIES[system]-1)*(r*r/(R*R));
		double vc_sq = G*m_orbit/r + haloVc*haloVc*r*r/(r*r + haloRc*haloRc);
		vm = sqrt(vc_sq);
		vm = (-1.0+Vtol*(2*drand()-1))*vm;

		vcopy(p,v); vnorm(v);
		vcross(v,y,v); vnorm(v);
		vscale(v,1.0*vm,v);
		vadd(p,pos[sysIdx],pos[sysIdx+i]);
		vadd(v,vel[sysIdx],vel[sysIdx+i]);
	}
}

void Simulation::LoadSphericalUniverseState(double M, double R, double V) {

	double m,r,theta,phi;
	double p[3],v[3];

	m = M/N_BODIES;
	double pp[3];

	std::mt19937 generator;
	std::normal_distribution<double> normal(0.0, 1.0);

	for (int i=0; i<N_BODIES; i++)
	{
		mass[i] = m + 0.0*m*normal(generator);
		has_gravity[i] = true;

		r = R*pow(drand(),(1.0/3.0));
		theta = 2*M_PI*drand();
		phi = acos(2*drand()-1);
		vset(r*sin(phi)*cos(theta),r*cos(phi),r*sin(phi)*sin(theta),p);

		pp[0] = normal(generator);
		pp[1] = normal(generator);
		pp[2] = normal(generator);
		vscale(pp,10.0,pp);
		vadd(p,pp,p);

		/*vcopy(p,v);
		vnorm(v);
		r = vmag(p);
		vscale(v,0.1*r,v);
		vrand(pp);
		vadd(v,pp,v);*/
		vrand(v);
		vscale(v,V,v);

		vcopy(p,pos[i]);
		vcopy(v,vel[i]);
	}
}

void Simulation::LoadDefaultState()
{

	G = 1.0;
	FDE = 0.0*1.0;
	dt = 0.0005;	
	r_soft = 0.1;
	BH_Opening_Theta = 0.5;
	double pos[3];
	double vel[3];

	FDE = 0.0;
	vset(0.0, 0.0, 0.0, pos);
	vset(0.0, 0.0, 0.0, vel);
	LoadGalaxyDiscState(0, pos, vel, 1.0e7, 0.5, 250.0, 25.0, 0.1, 200.0, 50.0);
	vset( 300.0, 0.0, -500.0, pos);
	vset(-150.0, 0.0, 0.0, vel);
	LoadGalaxyDiscState(1, pos, vel, 6.0e6, 0.5, 125.0, 25.0, 0.1, 155.0, 25.0);
	
	// FDE = 1.0;
	// LoadSphericalUniverseState(2.1e7, 210.0, 200.0);

	vset(0.0,500.0,500.0,Cam.pos);
	Cam.phi = atan2(-Cam.pos[1],Cam.pos[2]);
	Cam.theta = 0;
}

void Simulation::CalcAccelRangeP2P(int iStart, int iEnd) {

	double a[3];
	double r_halo[3];

	for (int i=iStart; i<=iEnd; i++)
	{
		vscaleadd(pos_t[i],FDE,acc_t[i]);

		for (int j=0; j<N_BODIES; j++)
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
		int ci = halo_central[sys];
		if (i != ci) {
			vsub(pos_t[ci], pos_t[i], r_halo);
			double rsq = vmagsq(r_halo);
			double halo_scale = halo_vc[sys] * halo_vc[sys] / (rsq + halo_rc_sq[sys]);
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
		int ci = halo_central[sys];
		if (bi != ci) {
			vsub(pos_t[ci], pos_t[bi], r_halo);
			double rsq = vmagsq(r_halo);
			double halo_scale = halo_vc[sys] * halo_vc[sys] / (rsq + halo_rc_sq[sys]);
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

void Simulation::CalcDerivatives(double *s, double *s_d)
{
	if (multiThreading) {
		int chunk_size = N_BODIES / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_BODIES-1) : (iStart + chunk_size - 1);
			pool->submit([this, s, s_d, iStart, iEnd]() { PrepareDerivativeDataRange(s, s_d, iStart, iEnd); });
		}
		pool->waitAll();

		int accel_chunk = numActiveBodies / numThreads;
		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * accel_chunk;
			int iEnd = (i == numThreads - 1) ? (numActiveBodies-1) : (iStart + accel_chunk - 1);
			if (GRAVITY_P2P) pool->submit([this, iStart, iEnd]() { CalcAccelRangeP2P(iStart, iEnd); });
			if (GRAVITY_OCT) pool->submit([this, iStart, iEnd]() { CalcAccelRangeOct(iStart, iEnd); });
		}
		pool->waitAll();
	} else {
		PrepareDerivativeDataRange(s, s_d, 0, N_BODIES-1);
		if (GRAVITY_P2P) CalcAccelRangeP2P(0, N_BODIES-1);
		if (GRAVITY_OCT) CalcAccelRangeOct(0, numActiveBodies-1);
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
	const double MAX_OCT_DST_SQ = 1000000.0;
	int numActive = 0;

	for (int i = 0; i < N_BODIES; i++)
	{
		if ((has_gravity[i]) && (pos_sq[i] < MAX_OCT_DST_SQ))
		{
			pos_f[i*3]   = (float)pos[i][0];
			pos_f[i*3+1] = (float)pos[i][1];
			pos_f[i*3+2] = (float)pos[i][2];
			sortedIdx[numActive++] = i;
		}
	}

	float p_min[3], p_max[3];
	vcopyf(pos_f + sortedIdx[0]*3, p_min);
	vcopyf(pos_f + sortedIdx[0]*3, p_max);
	for (int i = 1; i < numActive; i++)
	{
		vminf(pos_f + sortedIdx[i]*3, p_min, p_min);
		vmaxf(pos_f + sortedIdx[i]*3, p_max, p_max);
	}

	float halfwidth;
	float sub[3], mean[3], delta[3];
	vsubf(p_max, p_min, sub);
	halfwidth = sub[0];
	if (sub[1] > halfwidth) halfwidth = sub[1];
	if (sub[2] > halfwidth) halfwidth = sub[2];
	halfwidth *= 0.5f;
	vmeanf(p_max, p_min, mean);
	vsetf(halfwidth, halfwidth, halfwidth, delta);
	vsubf(mean, delta, p_min);
	vaddf(mean, delta, p_max);

	for (int i = 0; i < numActive; i++)
	{
		int bi = sortedIdx[i];
		mortonCodes[bi] = morton3D(pos_f[bi*3], pos_f[bi*3+1], pos_f[bi*3+2], p_min, p_max);
	}

	// Radix sort by Morton code (3 passes over 10-bit chunks)
	static int sortTemp[N_BODIES];
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
		memcpy(sortedIdx, sortTemp, numActive * sizeof(int));
	}

	Octree.Reset(p_min, p_max);

	for (int i = 0; i < numActive; i++)
	{
		int bi = sortedIdx[i];
		Octree.InsertBody(pos_f + bi*3, (float)mass[bi], bi);
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

void Simulation::CalcOutputs()
{
	if (multiThreading) {
		int chunk_size = N_BODIES / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_BODIES-1) : (iStart + chunk_size - 1);
			pool->submit([this, iStart, iEnd]() { CalcOutputsRange(iStart, iEnd); });
		}
		pool->waitAll();
	} else {
		CalcOutputsRange(0, N_BODIES-1);
	}
}

void Simulation::CalcRK4StateEstimateRange(double *s_est, double *s_curr, double *s_d, double scalar, double dt, int iStart, int iEnd) {
	
	for (int i=iStart; i<=iEnd; i++)
	{
		s_est[i] = s_curr[i] + dt*scalar*s_d[i];
	}
}

void Simulation::CalcRK4StateEstimate(double *s_est, double *s_curr, double *s_d, double scalar, double dt) {

	if (multiThreading) {
		int chunk_size = N_BODIES*N_STATES / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_BODIES*N_STATES-1) : (iStart + chunk_size - 1);
			pool->submit([this, s_est, s_curr, s_d, scalar, dt, iStart, iEnd]() {
				CalcRK4StateEstimateRange(s_est, s_curr, s_d, scalar, dt, iStart, iEnd);
			});
		}
		pool->waitAll();
	} else {
		CalcRK4StateEstimateRange(s_est, s_curr, s_d, scalar, dt, 0, N_BODIES*N_STATES-1);
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
		int chunk_size = N_BODIES / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_BODIES-1) : (iStart + chunk_size - 1);
			pool->submit([this, iStart, iEnd]() { CalcLeapFrogPositionsRange(iStart, iEnd); });
		}
		pool->waitAll();
	} else {
		CalcLeapFrogPositionsRange(0, N_BODIES-1);
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
		int chunk_size = N_BODIES / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_BODIES-1) : (iStart + chunk_size - 1);
			pool->submit([this, iStart, iEnd]() { CalcLeapFrogVelocitiesAndOutputsRange(iStart, iEnd); });
		}
		pool->waitAll();
	} else {
		CalcLeapFrogVelocitiesAndOutputsRange(0, N_BODIES-1);
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
		int chunk_size = N_BODIES / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_BODIES-1) : (iStart + chunk_size - 1);
			pool->submit([this, iStart, iEnd]() { CalcLeapFrogVelocitiesRange(iStart, iEnd); });
		}
		pool->waitAll();
	} else {
		CalcLeapFrogVelocitiesRange(0, N_BODIES-1);
	}
}

void Simulation::Step()
{
	if (SOLVER_RK4)
	{
		const double c[4] = {0.0, 1.0/2.0, 1.0/2.0, 1.0};
		const double b[4] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};

		CalcDerivatives(states,states_d[0]);
		for (int k=1; k<4; k++)
		{
			CalcRK4StateEstimate(states_e,states,states_d[k-1],c[k],dt);
			CalcDerivatives(states_e,states_d[k]);
		}

		for (int k=0; k<4; k++)
		{
			CalcRK4StateEstimate(states,states,states_d[k],b[k],dt);
		}
	}

	if (SOLVER_LEAP_FROG)
	{
		CalcLeapFrogPositions();
		CalcDerivatives(states,states_d[0]);
		CalcLeapFrogVelocitiesAndOutputs();
	}
	else
	{
		CalcOutputs();
	}

	if (GRAVITY_OCT)
	{
		BuildOctree();
	}

	if (DATA_LOG)
	{
		fwrite(pos_sq, N_BODIES*sizeof(double), 1, DataLog);
		fwrite(vel_sq, N_BODIES*sizeof(double), 1, DataLog);
		fwrite(acc_sq, N_BODIES*sizeof(double), 1, DataLog);
	}

	//CamMove(-0.0002/3.0,0.0,-0.5/3.0);
}

void Simulation::CamMove(double d_phi, double d_theta, double d_r)
{

	double r = vmag(Cam.pos);
	double phi = acos(Cam.pos[1]/r);
	double theta = atan2(Cam.pos[2],Cam.pos[0]);
	
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
	// Quad 1 triangle vertices in 2D billboard space
	float shapeData[24];
	// Quad 1, tri 1: (-d,-d), (-d,d), (d,d)
	shapeData[0]=-d; shapeData[1]=-d;
	shapeData[2]=-d; shapeData[3]= d;
	shapeData[4]= d; shapeData[5]= d;
	// Quad 1, tri 2: (-d,-d), (d,d), (d,-d)
	shapeData[6]=-d; shapeData[7]=-d;
	shapeData[8]= d; shapeData[9]= d;
	shapeData[10]=d; shapeData[11]=-d;
	// Quad 2 (rotated 45 degrees): corners rotated
	float c0x=-d*s45-(-d)*s45, c0y=-d*s45+(-d)*s45;  // (-d,-d) rotated 45
	float c1x=-d*s45-d*s45,    c1y=-d*s45+d*s45;      // (-d, d) rotated 45
	float c2x=d*s45-d*s45,     c2y=d*s45+d*s45;       // ( d, d) rotated 45
	float c3x=d*s45-(-d)*s45,  c3y=d*s45+(-d)*s45;    // ( d,-d) rotated 45
	// Quad 2, tri 1
	shapeData[12]=c0x; shapeData[13]=c0y;
	shapeData[14]=c1x; shapeData[15]=c1y;
	shapeData[16]=c2x; shapeData[17]=c2y;
	// Quad 2, tri 2
	shapeData[18]=c0x; shapeData[19]=c0y;
	shapeData[20]=c2x; shapeData[21]=c2y;
	shapeData[22]=c3x; shapeData[23]=c3y;

	glGenVertexArrays(1, &particleVAO);
	glGenBuffers(1, &particleShapeVBO);
	glGenBuffers(1, &particlePosVBO);
	glGenBuffers(1, &particleColorVBO);

	glBindVertexArray(particleVAO);

	// Attribute 0: shape offsets (per-vertex, from shape VBO)
	glBindBuffer(GL_ARRAY_BUFFER, particleShapeVBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(shapeData), shapeData, GL_STATIC_DRAW);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, nullptr);
	glEnableVertexAttribArray(0);
	glVertexAttribDivisor(0, 0);

	// Attribute 1: instance position (per-instance)
	glBindBuffer(GL_ARRAY_BUFFER, particlePosVBO);
	glBufferData(GL_ARRAY_BUFFER, N_BODIES * 3 * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
	glEnableVertexAttribArray(1);
	glVertexAttribDivisor(1, 1);

	// Attribute 2: instance color (per-instance)
	glBindBuffer(GL_ARRAY_BUFFER, particleColorVBO);
	glBufferData(GL_ARRAY_BUFFER, N_BODIES * 4 * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
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

	if (RECORD_VIDEO)
		CreateRecordFBO(width, height);
}

void Simulation::BuildOctreeVerts(int nodeIdx)
{
	BHNode &n = Octree.GetNode(nodeIdx);
	float x0 = (float)n.p_min[0], y0 = (float)n.p_min[1], z0 = (float)n.p_min[2];
	float x1 = (float)n.p_max[0], y1 = (float)n.p_max[1], z1 = (float)n.p_max[2];

	// 12 edges of the cube as line segments (24 vertices)
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
	if (RECORD_VIDEO)
		glBindFramebuffer(GL_FRAMEBUFFER, recordFBO);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// View = Rx(-phi) * Ry(-theta) * T(-cam)
	// Rx(-phi): rotates around X by -phi
	// Ry(-theta): rotates around Y by -theta
	double phi = Cam.phi, theta = Cam.theta;
	double cp = cos(-phi), sp = sin(-phi);
	double ct = cos(-theta), st = sin(-theta);

	// R = Rx * Ry (column-major)
	// Row 0: (ct, 0, -st)
	// Row 1: (sp*st, cp, sp*ct) -- wait, let me just expand properly
	// Rx(-phi) = {{1,0,0},{0,cp,-sp},{0,sp,cp}}  (using cp=cos(-phi), sp=sin(-phi))
	// Ry(-theta) = {{ct,0,st},{0,1,0},{-st,0,ct}}
	// R = Rx * Ry:
	//   row0: ct, 0, st
	//   row1: sp*(-st), cp, sp*ct  -- wait no: Rx*Ry[row1] = {0*ct + cp*0 + (-sp)*(-st), 0*0+cp*1+(-sp)*0, 0*st+cp*0+(-sp)*ct}
	//       = {sp*st, cp, -sp*ct}  (using sp=sin(-phi)=-sin(phi), so sp*st = -sin(phi)*(-sin(theta)) = sin(phi)*sin(theta))
	//   row2: {sp*ct+...}  let me just compute it directly

	// Rx(-phi) * Ry(-theta), element by element:
	// [0][0]=ct      [0][1]=0   [0][2]=st
	// [1][0]=-sp*st  [1][1]=cp  [1][2]=sp*ct   (note: actually sp*(-st) for first but sp=sin(-phi))
	// [2][0]=-cp*st  [2][1]=-sp [2][2]=cp*ct   (wait, let me redo)

	// Standard: Rx(a) = {{1,0,0},{0,cos(a),-sin(a)},{0,sin(a),cos(a)}}
	// Ry(b) = {{cos(b),0,sin(b)},{0,1,0},{-sin(b),0,cos(b)}}
	// With a=-phi, b=-theta:
	// Let's just use cp=cos(phi), sp=sin(phi), ct=cos(theta), st=sin(theta)
	// and compute Rx(-phi)*Ry(-theta) directly:
	// Rx(-phi) = {{1,0,0},{0,cos(phi),sin(phi)},{0,-sin(phi),cos(phi)}}
	// Ry(-theta) = {{cos(theta),0,-sin(theta)},{0,1,0},{sin(theta),0,cos(theta)}}
	// Product R = Rx(-phi)*Ry(-theta):
	// R[0] = {ct, 0, -st}
	// R[1] = {sp*st, cp, sp*ct}  -- wait: row1 = {0*ct+cp*0+sp*st, 0*0+cp*1+sp*0, 0*(-st)+cp*0+sp*ct} = {sp*st, cp, sp*ct}
	//   Actually: Rx(-phi) row1 = {0, cos(phi), sin(phi)}
	//   Ry(-theta) col0={cos(theta), 0, sin(theta)}, col1={0,1,0}, col2={-sin(theta),0,cos(theta)}
	//   R[1][0] = 0*ct + cos(phi)*0 + sin(phi)*st = sin(phi)*sin(theta)
	//   R[1][1] = 0*0 + cos(phi)*1 + sin(phi)*0 = cos(phi)
	//   R[1][2] = 0*(-st) + cos(phi)*0 + sin(phi)*ct = sin(phi)*cos(theta)
	//   Hmm that gives the same sign issue. Let me just re-derive properly.

	// Use the ORIGINAL angles from old code:
	// The old code did: glRotated(-phi_deg, 1,0,0); glRotated(-theta_deg, 0,1,0); glTranslated(-cam)
	// where phi_deg = Cam.phi * 180/PI, theta_deg = Cam.theta * 180/PI
	// So the view rotation is Rx(-Cam.phi) * Ry(-Cam.theta)

	double cphi = cos(Cam.phi), sphi = sin(Cam.phi);
	double ctheta = cos(Cam.theta), stheta = sin(Cam.theta);

	// Rx(-phi) * Ry(-theta), rows of the 3x3 rotation part:
	// Row 0: { cos(theta), 0, -sin(theta) }
	// Row 1: { -sin(phi)*(-sin(theta)), cos(phi), -sin(phi)*cos(theta) }
	//       = { sin(phi)*sin(theta), cos(phi), -sin(phi)*cos(theta) }
	// Row 2: { -cos(phi)*(-sin(theta)), -(-sin(phi)), -cos(phi)*... }
	//   Let me expand Rx(-phi) = {{1,0,0},{0,cos(phi),sin(phi)},{0,-sin(phi),cos(phi)}}
	//   Ry(-theta) = {{cos(theta),0,-sin(theta)},{0,1,0},{sin(theta),0,cos(theta)}}
	//   R = Rx(-phi) * Ry(-theta):
	//   R[0][0..2] = {cos(theta), 0, -sin(theta)}  (first row of Rx is {1,0,0})
	//   R[1][0] = 0*cos(theta) + cos(phi)*0 + sin(phi)*sin(theta) = sin(phi)*sin(theta)
	//   R[1][1] = 0*0 + cos(phi)*1 + sin(phi)*0 = cos(phi)
	//   R[1][2] = 0*(-sin(theta)) + cos(phi)*0 + sin(phi)*cos(theta) = sin(phi)*cos(theta)
	//   R[2][0] = 0*cos(theta) + (-sin(phi))*0 + cos(phi)*sin(theta) = cos(phi)*sin(theta)
	//   R[2][1] = 0*0 + (-sin(phi))*1 + cos(phi)*0 = -sin(phi)
	//   R[2][2] = 0*(-sin(theta)) + (-sin(phi))*0 + cos(phi)*cos(theta) = cos(phi)*cos(theta)

	// Column-major view matrix (rotation part in columns, translation in column 3)
	float r00=(float)ctheta,       r10=(float)(sphi*stheta),  r20=(float)(cphi*stheta);
	float r01=0.0f,                r11=(float)cphi,           r21=(float)(-sphi);
	float r02=(float)(-stheta),    r12=(float)(sphi*ctheta),  r22=(float)(cphi*ctheta);

	// Translation = R * (-camPos)
	float cx=(float)Cam.pos[0], cy=(float)Cam.pos[1], cz=(float)Cam.pos[2];
	float tx = -(r00*cx + r01*cy + r02*cz);
	float ty = -(r10*cx + r11*cy + r12*cz);
	float tz = -(r20*cx + r21*cy + r22*cz);

	// Column-major 4x4
	float view[16] = {
		r00, r10, r20, 0,
		r01, r11, r21, 0,
		r02, r12, r22, 0,
		tx,  ty,  tz,  1
	};

	// Projection matrix (perspective, column-major)
	float aspect = (float)winWidth / (float)winHeight;
	float fov = 45.0f * (float)M_PI / 180.0f;
	float zNear = 0.1f, zFar = 10000.0f;
	float f = 1.0f / tanf(fov / 2.0f);
	float proj[16] = {
		f/aspect, 0, 0, 0,
		0, f, 0, 0,
		0, 0, (zFar+zNear)/(zNear-zFar), -1,
		0, 0, 2*zFar*zNear/(zNear-zFar), 0
	};

	// VP = proj * view (column-major)
	float vp[16];
	for (int c = 0; c < 4; c++) {
		for (int r = 0; r < 4; r++) {
			vp[c*4+r] = 0;
			for (int k = 0; k < 4; k++) {
				vp[c*4+r] += proj[k*4+r] * view[c*4+k];
			}
		}
	}

	// Billboard vectors (row 0 and row 1 of the view rotation = camera right and up in world space)
	float right[3] = { r00, (float)0.0f, r02 };  // {cos(theta), 0, -sin(theta)}
	float up[3] = { r10, r11, r12 };              // {sin(phi)*sin(theta), cos(phi), sin(phi)*cos(theta)}

	// --- Update particle data ---
	const float alpha_base = 10000.0f / N_BODIES;
	const float alpha_clamped = (alpha_base > 0.2f) ? 0.2f : ((alpha_base < 0.02f) ? 0.02f : alpha_base);

	int sysIndices[N_SYSTEMS];
	sysIndices[0] = 0;
	for (int j = 1; j < N_SYSTEMS; j++) sysIndices[j] = sysIndices[j-1] + N_SYSTEM_BODIES[j-1];

	for (int i = 0; i < N_BODIES; i++) {
		posBuf[i*3+0] = (float)pos[i][0];
		posBuf[i*3+1] = (float)pos[i][1];
		posBuf[i*3+2] = (float)pos[i][2];

		bool sysBody = false;
		for (int j = 0; j < N_SYSTEMS; j++) {
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

	// Upload and draw particles (instanced: 12 verts per star shape, N_BODIES instances)
	glUseProgram(particleShader);
	glUniformMatrix4fv(glGetUniformLocation(particleShader, "uVP"), 1, GL_FALSE, vp);
	glUniform3fv(glGetUniformLocation(particleShader, "uRight"), 1, right);
	glUniform3fv(glGetUniformLocation(particleShader, "uUp"), 1, up);

	glBindVertexArray(particleVAO);
	// Buffer orphaning: re-allocate before uploading to avoid implicit sync stalls
	glBindBuffer(GL_ARRAY_BUFFER, particlePosVBO);
	glBufferData(GL_ARRAY_BUFFER, N_BODIES * 3 * sizeof(float), posBuf, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, particleColorVBO);
	glBufferData(GL_ARRAY_BUFFER, N_BODIES * 4 * sizeof(float), clrBuf, GL_DYNAMIC_DRAW);
	glDrawArraysInstanced(GL_TRIANGLES, 0, 12, N_BODIES);
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

	char fpsText[16];
	snprintf(fpsText, sizeof(fpsText), "FPS:%d", (int)(fps + 0.5));

	std::vector<float> verts;
	float x = (float)winWidth - (float)strlen(fpsText) * 6.0f - 4.0f;
	float y = (float)winHeight - 12.0f;
	int scale = 1;

	for (const char *p = fpsText; *p; p++) {
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

	if (verts.empty()) return;

	// Orthographic projection
	float proj[16] = {
		2.0f/winWidth, 0, 0, 0,
		0, 2.0f/winHeight, 0, 0,
		0, 0, -1, 0,
		-1, -1, 0, 1
	};

	glUseProgram(fpsShader);
	glUniformMatrix4fv(glGetUniformLocation(fpsShader, "uProj"), 1, GL_FALSE, proj);
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
	// Read pixels from the FBO
	glBindFramebuffer(GL_READ_FRAMEBUFFER, recordFBO);
	glReadPixels(0, 0, winWidth, winHeight, GL_RGB, GL_UNSIGNED_BYTE, rgbOut);

	// Blit FBO to the default framebuffer for on-screen display
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	glBlitFramebuffer(0, 0, winWidth, winHeight, 0, 0, winWidth, winHeight,
					  GL_COLOR_BUFFER_BIT, GL_NEAREST);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);

	// glReadPixels returns bottom-up; flip vertically
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

	fwrite(states,		sizeof(double), N_BODIES*N_STATES,	StateFile);
	fwrite(mass,		sizeof(double), N_BODIES,			StateFile);
	fwrite(has_gravity, sizeof(bool),	N_BODIES,			StateFile);

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

	result = fread(states,		sizeof(double), N_BODIES*N_STATES,	StateFile);
	result = fread(mass,		sizeof(double), N_BODIES,			StateFile);
	result = fread(has_gravity,	sizeof(bool),	N_BODIES,			StateFile);
	
	fclose(StateFile);

	return true;
}