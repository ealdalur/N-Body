#include "Simulation.h"

Simulation::Simulation()
{
	multiThreading = true;
	int hardwareConcurrency = std::thread::hardware_concurrency();
	std::cout << "Hardware thread contexts: " << hardwareConcurrency << std::endl;
	if (hardwareConcurrency > 8) {
		numThreads = 8;
	} else {
		numThreads = 4;
	}
	std::cout << "Number of compute threads: " << numThreads << std::endl;

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
	
	double p[3];
	vset(0.0,0.0,0.0,p);
	Octree = new BHTree(p,p,nullptr);
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
	glDeleteLists(dList,2);
	delete Octree;

	if (DATA_LOG) fclose(DataLog);
}

void Simulation::LoadGalaxyDiscState(int system, double *sysPos, double *sysVel, double M, double Mfrac, double R, double Ri, double Vtol){
	
	double m,r,theta,phi,vm,m_orbit;
	double p[3],v[3];
	double y[3] = {0.0, 1.0, 0.0};

	int sysIdx = 0;
	for (int i=0; i<system; i++) sysIdx += N_SYSTEM_BODIES[i];

	m = Mfrac*M/(N_SYSTEM_BODIES[system]-1);
	mass[sysIdx] = M;
	has_gravity[sysIdx] = true;
	//vset(0.0,0.0,0.0,pos[sysIdx]);
	//vset(0.0,0.0,0.0,vel[sysIdx]);
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
		vm = sqrt(G*m_orbit/r);
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
	double pos[3];
	double vel[3];

	FDE = 0.0;
	vset(0.0, 0.0, 0.0, pos);
	vset(0.0, 0.0, 0.0, vel);
	LoadGalaxyDiscState(0, pos, vel, 1.0e7, 0.5, 250.0, 25.0, 0.1);
	vset( 300.0, 0.0, -500.0, pos);
	vset(-150.0, 0.0, 0.0, vel);
	LoadGalaxyDiscState(1, pos, vel, 6.0e6, 0.5, 125.0, 25.0, 0.1);
	
	// FDE = 1.0;
	// LoadSphericalUniverseState(2.1e7, 210.0, 200.0);

	vset(0.0,500.0,500.0,Cam.pos);
	Cam.phi = atan2(-Cam.pos[1],Cam.pos[2]);
	Cam.theta = 0;
}

void Simulation::CalcAccelRangeP2P(int iStart, int iEnd) {

	double a[3];
	double r;

	for (int i=iStart; i<=iEnd; i++)
		{
			vscaleadd(pos_t[i],FDE,acc_t[i]);

			for (int j=0; j<N_BODIES; j++)
			{
				if ((j != i) && (has_gravity[j]))
				{
					vsub(pos_t[j],pos_t[i],a);
					r = vmagsoft(a,r_soft);
					vscale(a,G/(r*r*r),a);
					if (has_gravity[j]) vscaleadd(a,mass[j],acc_t[i]);
				}
			}
		}
}

void Simulation::CalcAccelRangeOct(int iStart, int iEnd) {

	double a[3];

	for (int i=iStart; i<=iEnd; i++)
	{
		vscaleadd(pos_t[i],FDE,acc_t[i]);

		Octree->CalcAcceleration(pos_t[i],i,G,r_soft,a);
		vadd(acc_t[i],a,acc_t[i]);
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
		threads.clear();
		int chunk_size = N_BODIES / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_BODIES-1) : (iStart + chunk_size - 1);
			threads.emplace_back(&Simulation::PrepareDerivativeDataRange, this, s, s_d, iStart, iEnd);
		}

		for (auto& thread : threads) {
			thread.join();
		}

		threads.clear();
		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_BODIES-1) : (iStart + chunk_size - 1);
			if (GRAVITY_P2P) threads.emplace_back(&Simulation::CalcAccelRangeP2P, this, iStart, iEnd);
			if (GRAVITY_OCT) threads.emplace_back(&Simulation::CalcAccelRangeOct, this, iStart, iEnd);
		}

		for (auto& thread : threads) {
			thread.join();
		}
	} else {
		PrepareDerivativeDataRange(s, s_d, 0, N_BODIES-1);
		if (GRAVITY_P2P) CalcAccelRangeP2P(0, N_BODIES-1);
		if (GRAVITY_OCT) CalcAccelRangeOct(0, N_BODIES-1);
	}
	
}

void Simulation::BuildOctree()
{
	double p_min[3],p_max[3];
	bool first = true;
	const double MAX_OCT_DST_SQ = 1000000.0;

	for (int i=0; i<N_BODIES; i++)
	{	
		if ((has_gravity[i]) && (pos_sq[i] < MAX_OCT_DST_SQ))
		{
			if (first)
			{
				vcopy(pos[i],p_min);
				vcopy(pos[i],p_max);
				first = false;
			}
			else
			{
				vmin(pos[i],p_min,p_min);
				vmax(pos[i],p_max,p_max);
			}
		}
	}

	vmaxboundcube(p_min,p_max,p_min,p_max);
	Octree->Reset(p_min,p_max);
	
	for (int i=0; i<N_BODIES; i++)
	{
		if ((has_gravity[i]) && (pos_sq[i] < MAX_OCT_DST_SQ))
		{
			Octree->InsertBody(pos[i],mass[i],i);
		}
	}

	Octree->CalcMasses();
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
		threads.clear();
		int chunk_size = N_BODIES / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_BODIES-1) : (iStart + chunk_size - 1);
			threads.emplace_back(&Simulation::CalcOutputsRange, this, iStart, iEnd);
		}

		for (auto& thread : threads) {
			thread.join();
		}
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
		threads.clear();
		int chunk_size = N_BODIES*N_STATES / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_BODIES*N_STATES-1) : (iStart + chunk_size - 1);
			threads.emplace_back(&Simulation::CalcRK4StateEstimateRange, this, s_est, s_curr, s_d, scalar, dt, iStart, iEnd);
		}

		for (auto& thread : threads) {
			thread.join();
		}
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
		threads.clear();
		int chunk_size = N_BODIES / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_BODIES-1) : (iStart + chunk_size - 1);
			threads.emplace_back(&Simulation::CalcLeapFrogPositionsRange, this, iStart, iEnd);
		}

		for (auto& thread : threads) {
			thread.join();
		}
	} else {
		CalcLeapFrogPositionsRange(0, N_BODIES-1);
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
		threads.clear();
		int chunk_size = N_BODIES / numThreads;

		for (int i = 0; i < numThreads; ++i) {
			int iStart = i * chunk_size;
			int iEnd = (i == numThreads - 1) ? (N_BODIES-1) : (iStart + chunk_size - 1);
			threads.emplace_back(&Simulation::CalcLeapFrogVelocitiesRange, this, iStart, iEnd);
		}

		for (auto& thread : threads) {
			thread.join();
		}
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
		CalcLeapFrogVelocities();
	}

	CalcOutputs();

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

void Simulation::InitGL()
{
	//glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);

	glDisable(GL_DEPTH_TEST);
	glDepthMask(false);
	glDisable(GL_CULL_FACE);
	glDisable(GL_LIGHTING);
	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
	glBlendFunc(GL_SRC_ALPHA,GL_ONE);
	glEnable(GL_BLEND);
	
	glColor4d(1.0,0.95,0.92,0.1);	
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	double d = 0.5;	
	dList = glGenLists(2);
	glNewList(dList,GL_COMPILE);
		glPushMatrix();
		glBegin(GL_QUADS);		
			glVertex3d(-d,-d, 0.0);
			glVertex3d(-d, d, 0.0);
			glVertex3d( d, d, 0.0);
			glVertex3d( d,-d, 0.0);
		glEnd();
		glRotated(45,0,0,1);
		glBegin(GL_QUADS);		
			glVertex3d(-d,-d, 0.0);
			glVertex3d(-d, d, 0.0);
			glVertex3d( d, d, 0.0);
			glVertex3d( d,-d, 0.0);
		glEnd();
		glPopMatrix();
	glEndList();

	d = 0.5;
	glNewList(dList+1,GL_COMPILE);
		glColor4d(0.0,0.0,1.0,0.1);
		glDisable(GL_BLEND);
		glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
		glBegin(GL_QUADS);		
			glVertex3d(-d,-d, d);
			glVertex3d(-d, d, d);
			glVertex3d( d, d, d);
			glVertex3d( d,-d, d);
		glEnd();
		glBegin(GL_QUADS);		
			glVertex3d(-d,-d,-d);
			glVertex3d(-d, d,-d);
			glVertex3d( d, d,-d);
			glVertex3d( d,-d,-d);
		glEnd();
		glBegin(GL_QUADS);		
			glVertex3d(-d, d,-d);
			glVertex3d(-d, d, d);
			glVertex3d( d, d, d);
			glVertex3d( d, d,-d);
		glEnd();
		glBegin(GL_QUADS);		
			glVertex3d(-d,-d,-d);
			glVertex3d(-d,-d, d);
			glVertex3d( d,-d, d);
			glVertex3d( d,-d,-d);
		glEnd();
		glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
		glEnable(GL_BLEND);
		glColor4d(1.0,0.95,0.92,0.1);
	glEndList();
}

void Simulation::SetPerspective(double fovY, double aspect, double zNear, double zFar) {
    
	double fH = tan(fovY / 360 * M_PI) * zNear;
    double fW = fH * aspect;
    glFrustum(-fW, fW, -fH, fH, zNear, zFar);
}


void Simulation::ReSizeGL(int width, int height)
{
	if (height==0)
	{
		height=1;
	}

	glViewport(0,0,width,height);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	SetPerspective(45.0,(double)width/(double)height,0.1,10000.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void Simulation::DrawOctant(BHTree *Oct)
{
	double bd[4];
	Oct->GetBounds(bd);
	glPushMatrix();
		glTranslated(bd[0],bd[1],bd[2]);
		glScaled(bd[3],bd[3],bd[3]);
		glCallList(dList+1);
	glPopMatrix();

	for (int i=0; i<8; i++)
	{
		if (Oct->GetOctant(i))
		{
			DrawOctant(Oct->GetOctant(i));
		}
	}
}

void Simulation::DrawGL(GLvoid)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	double phi, theta;

	phi = 180.0/M_PI*Cam.phi;
	theta = 180.0/M_PI*Cam.theta;

	glRotated(-phi,1.0,0.0,0.0);
	glRotated(-theta,0.0,1.0,0.0);
	glTranslated(-Cam.pos[0],-Cam.pos[1],-Cam.pos[2]);

	double r,g,b,a;
	bool sysBody = false;
	int sysIdx = 0;
	for (int i=0; i<N_BODIES; i++)
	{
		glPushMatrix();
			glTranslated(pos[i][0],pos[i][1],pos[i][2]);
			glRotated(theta,0.0f,1.0f,0.0f);
			glRotated(phi,1.0f,0.0f,0.0f);

			sysBody = false;
			sysIdx = 0;
			for (int j=0; j<N_SYSTEMS; j++) {
				if (sysIdx == i) sysBody = true;
				sysIdx += N_SYSTEM_BODIES[j];
			}

			if (sysBody) {
				r = 0.0;
				g = 1.0;
				b = 0.0;
				a = 1.0;
			} else {
				r = pow(acc_sq[i]/1000000.0,1.0/3.0);
				g = 0.3;
				b = 1.0 - r;
					
				r = (r<0.3)?0.3:r;
				g = (g<0.3)?0.3:g;
				b = (b<0.3)?0.3:b;
					
				a = 10000.0/N_BODIES;
				a = (a>0.2)?0.2:a;
				a = (a<0.02)?0.02:a;
			}
			
			glColor4d(r,g,b,a);
			

			glCallList(dList);		
		glPopMatrix();
	}

	if (DrawOctree) DrawOctant(Octree);
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