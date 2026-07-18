#pragma once

#include "VectorMath.h"
#include <vector>

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

struct BHNode
{
	double p_min[3], p_max[3], p_c[3], p_cm[3];
	double Mass;
	int BodyIdx;
	int Num_Bodies;
	int Octants[8];
};

class BHTree
{
private:
	std::vector<BHNode> nodes;
	int nodeCount;

	void InitNode(int idx, double *p_min_new, double *p_max_new);
	int AllocNode(double *p_min_new, double *p_max_new);
	Octant DetermineOctant(int idx, double *p);
	int CreateChild(int idx, Octant Oct);

public:
	BHTree();

	void Reset(double *p_min_new, double *p_max_new);
	void InsertBody(double *p, double m, int BodyIdxIn);
	void CalcMasses();
	void CalcMassesNode(int idx);
	void CalcAcceleration(double *p, int BodyIdxIn, double G, double r_soft, double *a);
	void GetBounds(int idx, double *BoundData);
	int GetRoot() { return 0; }
	int GetOctant(int idx, int oct) { return nodes[idx].Octants[oct]; }
	int GetNodeCount() { return nodeCount; }
	BHNode& GetNode(int idx) { return nodes[idx]; }
};
