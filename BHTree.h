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
	float p_min[3], p_max[3], p_c[3], p_cm[3];
	float Mass;
	int BodyIdx;
	int Num_Bodies;
	int Octants[8];
};

class BHTree
{
private:
	std::vector<BHNode> nodes;
	int nodeCount;

	void InitNode(int idx, float *p_min_new, float *p_max_new);
	int AllocNode(float *p_min_new, float *p_max_new);
	Octant DetermineOctant(int idx, float *p);
	int CreateChild(int idx, Octant Oct);

public:
	BHTree();

	void Reset(float *p_min_new, float *p_max_new);
	void InsertBody(float *p, float m, int BodyIdxIn);
	void CalcMasses();
	void CalcMassesNode(int idx);
	void CalcAcceleration(float *p, int BodyIdxIn, float G, float r_soft, float theta_sq, double *a);
	void GetBounds(int idx, float *BoundData);
	int GetRoot() { return 0; }
	int GetOctant(int idx, int oct) { return nodes[idx].Octants[oct]; }
	int GetNodeCount() { return nodeCount; }
	BHNode& GetNode(int idx) { return nodes[idx]; }
};
