#include "BHTree.h"

BHTree::BHTree(double *p_min_new, double *p_max_new, BHTree *Parent_new)
{
	Parent = Parent_new;
	for (int i=0; i<8; i++)
	{
		Octants[i] = nullptr;
	}
	
	Initialize(p_min_new,p_max_new);
}

BHTree::~BHTree()
{
	DeleteOctants();
}

void BHTree::Reset(double *p_min_new, double *p_max_new)
{
	DeleteOctants();
	Initialize(p_min_new,p_max_new);
}

void BHTree::Initialize(double *p_min_new, double *p_max_new)
{
	vcopy(p_min_new,p_min);
	vcopy(p_max_new,p_max);
	vmean(p_min,p_max,p_c);
	vcopy(p_c,p_cm);
	Mass = 0;
	Num_Bodies = 0;
	BodyIdx = -1;
}

void BHTree::DeleteOctants()
{
	for (int i=0; i<8; i++)
	{
		delete Octants[i];
		Octants[i] = nullptr;
	}
}

Octant BHTree::DetermineOctant(double *p)
{
	if (p[2] >= p_c[2])
	{
		if (p[1] >= p_c[1])
		{
			if (p[0] >= p_c[0])
			{
				return FNE;
			}
			else
			{
				return FNW;
			}
		}
		else
		{
			if (p[0] >= p_c[0])
			{
				return FSE;
			}
			else
			{
				return FSW;
			}
		}
	}
	else
	{
		if (p[1] >= p_c[1])
		{
			if (p[0] >= p_c[0])
			{
				return BNE;
			}
			else
			{
				return BNW;
			}
		}
		else
		{
			if (p[0] >= p_c[0])
			{
				return BSE;
			}
			else
			{
				return BSW;
			}
		}
	}
}

BHTree *BHTree::CreateNode(Octant Oct)
{
	double p_min_new[3], p_max_new[3];
	switch (Oct)
	{
		case FNE: 
			vset(p_c[0],	p_c[1],		p_c[2],		p_min_new);
			vset(p_max[0],	p_max[1],	p_max[2],	p_max_new);
			break;
		case FNW: 
			vset(p_min[0],	p_c[1],		p_c[2],		p_min_new);
			vset(p_c[0],	p_max[1],	p_max[2],	p_max_new);
			break;
		case FSW: 
			vset(p_min[0],	p_min[1],	p_c[2],		p_min_new);
			vset(p_c[0],	p_c[1],		p_max[2],	p_max_new);
			break;
		case FSE: 
			vset(p_c[0],	p_min[1],	p_c[2],		p_min_new);
			vset(p_max[0],	p_c[1],		p_max[2],	p_max_new);
			break;
		case BNE: 
			vset(p_c[0],	p_c[1],		p_min[2],	p_min_new);
			vset(p_max[0],	p_max[1],	p_c[2],		p_max_new);
			break;
		case BNW: 
			vset(p_min[0],	p_c[1],		p_min[2],	p_min_new);
			vset(p_c[0],	p_max[1],	p_c[2],		p_max_new);
			break;
		case BSW: 
			vset(p_min[0],	p_min[1],	p_min[2],	p_min_new);
			vset(p_c[0],	p_c[1],		p_c[2],		p_max_new);
			break;
		case BSE: 
			vset(p_c[0],	p_min[1],	p_min[2],	p_min_new);
			vset(p_max[0],	p_c[1],		p_c[2],		p_max_new);
			break;
	}

	return new BHTree(p_min_new,p_max_new,this);
}

void BHTree::InsertBody(double *p, double m, int BodyIdxIn)
{
	Octant Oct;

	if (Num_Bodies == 0)
	{
		vcopy(p,p_cm);
		Mass = m;
		BodyIdx = BodyIdxIn;
	}
	else
	{
		if (Num_Bodies == 1)
		{
			Oct = DetermineOctant(p_cm);
			if (!Octants[Oct]) Octants[Oct] = CreateNode(Oct);
			if (vequal(p_cm,p))
			{
				Octants[Oct]->InsertBody(p_cm,Mass+m,BodyIdx);
			}
			else
			{
				Octants[Oct]->InsertBody(p_cm,Mass,BodyIdx);

				Oct = DetermineOctant(p);
				if (!Octants[Oct]) Octants[Oct] = CreateNode(Oct);
				Octants[Oct]->InsertBody(p,m,BodyIdxIn);
			}
		}
		else
		{
			Oct = DetermineOctant(p);
			if (!Octants[Oct]) Octants[Oct] = CreateNode(Oct);
			Octants[Oct]->InsertBody(p,m,BodyIdxIn);
		}
	}

	Num_Bodies++;
}

void BHTree::CalcMasses()
{
	if (Num_Bodies > 1)
	{
		Mass = 0.0;
		vset(0.0,0.0,0.0,p_cm);
		double oct_cm[3];

		for (int i=0; i<8; i++)
		{
			if (Octants[i])
			{
				Octants[i]->CalcMasses();
				Mass += Octants[i]->Mass;
				vscale(Octants[i]->p_cm,Octants[i]->Mass,oct_cm);
				vadd(p_cm,oct_cm,p_cm);
			}
		}
		
		vscale(p_cm,1.0/Mass,p_cm);
	}
}

void BHTree::CalcAcceleration(double *p, int BodyIdxIn, double G, double r_soft, double *a)
{
	vset(0.0,0.0,0.0,a);
	double v[3];
	double r,d;
	
	if (Num_Bodies == 1)
	{
		if (BodyIdx != BodyIdxIn)
		{
			vsub(p_cm,p,v);
			r = vmagsoft(v,r_soft);
			vscale(v,G*Mass/(r*r*r),a);
		}
	}
	else
	{
		vsub(p_cm,p,v);
		r = vmag(v);
		d = p_max[0] - p_min[0];
		if (d/r <= 1.0)
		{
			vscale(v,G*Mass/(r*r*r),a);
		}
		else
		{
			double a_oct[3];
			for (int i=0; i<8; i++)
			{
				if (Octants[i])
				{
					Octants[i]->CalcAcceleration(p,BodyIdxIn,G,r_soft,a_oct);
					vadd(a,a_oct,a);
				}
			}
		}
	}
}

void BHTree::GetBounds(double *p_min_out, double *p_max_out)
{
	vcopy(p_min,p_min_out);
	vcopy(p_max,p_max_out);
}

void BHTree::GetBounds(double *BoundData)
{
	vcopy(p_c,BoundData);
	BoundData[3] = p_max[0] - p_min[0];
}

BHTree *BHTree::GetOctant(int Oct)
{
	return Octants[Oct];
}