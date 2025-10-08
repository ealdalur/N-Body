#include <random>
#include "Simulation.h"

Simulation::Simulation()
{
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
		fopen_s(&DataLog, "Simulation.dat", "wb");
		fwrite(&N_BODIES, sizeof(N_BODIES), 1, DataLog);
	}
}

Simulation::~Simulation()
{
	glDeleteLists(dList,2);
	delete Octree;

	if (DATA_LOG) fclose(DataLog);
}

void Simulation::LoadGalaxyDiscState(double M, double Mfrac, double R, double Ri, double Vtol){
	
	double m,r,theta,phi,vm,m_orbit;
	double p[3],v[3];
	double y[3] = {0.0, 1.0, 0.0};

	m = Mfrac*M/(N_SYSTEM_BODIES[0]-1);
	mass[0] = M;
	has_gravity[0] = true;
	vset(0.0,0.0,0.0,pos[0]);
	vset(0.0,0.0,0.0,vel[0]);

	for (int i=1; i<N_SYSTEM_BODIES[0]; i++)
	{
		mass[i] = m;
		has_gravity[i] = true;
		
		r = (R-Ri)*drand() + Ri;
		
		phi = M_PI/2 + M_PI/64*(2*drand()-1);
		theta = 2*M_PI*drand();
		vset(r*sin(phi)*cos(theta),r*cos(phi),r*sin(phi)*sin(theta),p);

		m_orbit = M + m*(N_SYSTEM_BODIES[0]-1)*(r*r/(R*R));
		vm = sqrt(G*m_orbit/r);
		vm = (-1.0+Vtol*(2*drand()-1))*vm;

		vcopy(p,v); vnorm(v);
		vcross(v,y,v); vnorm(v);
		vscale(v,1.0*vm,v);
		vadd(p,pos[0],pos[i]);
		vadd(v,vel[0],vel[i]);
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
	FDE = 0*1.0;
	dt = 0.001;	
	r_soft = 0.5;
	
	LoadGalaxyDiscState(1.0e7, 0.1, 250.0, 25.0, 0.1);
	//LoadSphericalUniverseState(1.0e7, 217.5, 200.0);

	vset(0.0,500.0,500.0,Cam.pos);
	Cam.phi = atan2(-Cam.pos[1],Cam.pos[2]);
	Cam.theta = 0;
}

void Simulation::CalcDerivatives(double *s, double *s_d)
{
	double a[3];
	double r;
		
	int ki;
	for (int i=0; i<N_BODIES; i++)
	{
		ki = i*N_STATES;

		pos_t[i] = s+ki;
		acc_t[i] = s_d+ki+3;

		vcopy(s+ki+3,s_d+ki);
		vset(0.0,0.0,0.0,acc_t[i]);
	}

	for (int i=0; i<N_BODIES; i++)
	{
		
		vscaleadd(pos_t[i],FDE,acc_t[i]);

		if (GRAVITY_P2P)
		{
			for (int j=0; j<N_BODIES; j++)
			{
				if (j > i)
				{
					vsub(pos_t[j],pos_t[i],a);
					r = vmagsoft(a,r_soft);
					vscale(a,G/(r*r*r),a);
					if (has_gravity[j]) vscaleadd(a,mass[j],acc_t[i]);
					if (has_gravity[i]) vscaleadd(a,-1.0*mass[i],acc_t[j]);
				}
			}
		}
		else if (GRAVITY_OCT)
		{
			Octree->CalcAcceleration(pos_t[i],i,G,r_soft,a);
			vadd(acc_t[i],a,acc_t[i]);
		}	
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

void Simulation::CalcOutputs()
{
	for (int i=0; i<N_BODIES; i++)
	{
		pos_sq[i] = vdot(pos[i],pos[i]);
		vel_sq[i] = vdot(vel[i],vel[i]);
		acc_sq[i] = vdot(acc[i],acc[i]);
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
			for (int i=0; i<(N_BODIES*N_STATES); i++)
			{
				states_e[i] = states[i] + dt*c[k]*states_d[k-1][i];
			}
			CalcDerivatives(states_e,states_d[k]);
		}

		for (int k=0; k<4; k++)
		{
			for (int i=0; i<(N_BODIES*N_STATES); i++)
			{
				states[i] = states[i] + dt*b[k]*states_d[k][i];
			}
		}
	}

	if (SOLVER_LEAP_FROG)
	{
		double a[3];
		for (int i=0; i<N_BODIES; i++)
		{
			vcopy(acc[i],acc_prev[i]);

			vscaleadd(vel[i],dt,pos[i]);
			vscaleadd(acc[i],0.5*dt*dt,pos[i]);
		}
		CalcDerivatives(states,states_d[0]);
		for (int i=0; i<N_BODIES; i++)
		{
			vadd(acc[i],acc_prev[i],a);
			vscaleadd(a,0.5*dt,vel[i]);
		}
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
	for (int i=0; i<N_BODIES; i++)
	{
		glPushMatrix();
			glTranslated(pos[i][0],pos[i][1],pos[i][2]);
			glRotated(theta,0.0f,1.0f,0.0f);
			glRotated(phi,1.0f,0.0f,0.0f);


			r = pow(acc_sq[i]/400000.0,1.0/3.0);
			g = 0.3;
			b = 1.0 - r;
				
			r = (r<0.3)?0.3:r;
			g = (g<0.3)?0.3:g;
			b = (b<0.3)?0.3:b;
				
			a = 10000.0/N_BODIES;
			a = (a>0.2)?0.2:a;
			a = (a<0.02)?0.02:a;
			glColor4d(r,g,b,a);
			

			glCallList(dList);		
		glPopMatrix();
	}

	if (DrawOctree) DrawOctant(Octree);
}

void Simulation::SaveState()
{
	FILE *StateFile;

	fopen_s(&StateFile, "State.dat", "wb");
	
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

	fopen_s(&StateFile, "State.dat", "rb");
	if (StateFile == nullptr) 
	{
        return false;
    }
	
	fread(&G,			sizeof(double), 1, StateFile);
	fread(&FDE,			sizeof(double), 1, StateFile);
	fread(&dt,			sizeof(double), 1, StateFile);
	fread(&r_soft,		sizeof(double), 1, StateFile);

	fread(Cam.pos,		sizeof(double), 3, StateFile);
	fread(&Cam.phi,		sizeof(double), 1, StateFile);
	fread(&Cam.theta,	sizeof(double), 1, StateFile);

	fread(states,		sizeof(double), N_BODIES*N_STATES,	StateFile);
	fread(mass,			sizeof(double), N_BODIES,			StateFile);
	fread(has_gravity,	sizeof(bool),	N_BODIES,			StateFile);
	
	fclose(StateFile);

	return true;
}