#include "VectorMath.h"

enum Octant
{
	FNE = 0,
	FNW,
	FSW,
	FSE,
	BNE,
	BNW,
	BSW,
	BSE,
};

class BHTree
{
private:
	double p_min[3], p_max[3], p_c[3], p_cm[3];
	double Mass;
	int BodyIdx;
	int Num_Bodies;
	BHTree *Parent;
	BHTree *Octants[8];

	void Initialize(double *p_min_new, double *p_max_new);
	void DeleteOctants();
	Octant DetermineOctant(double *p);
	BHTree *CreateNode(Octant Oct);
public:
	BHTree(double *p_min_new, double *p_max_new, BHTree *Parent_new);
	~BHTree();

	void Reset(double *p_min_new, double *p_max_new);
	void InsertBody(double *p, double m, int BodyIdxIn);
	void CalcMasses();
	void CalcAcceleration(double *p, int BodyIdxIn, double G, double r_soft, double *a);
	void GetBounds(double *p_min_out, double *p_max_out);
	void GetBounds(double *BoundData);
	BHTree *GetOctant(int Oct);
};
