#include "BHTree.h"
#include <math.h>
#include <cstring>

BHTree::BHTree()
{
	nodes.resize(512000);
	nodeCount = 0;
}

void BHTree::InitNode(int idx, float *p_min_new, float *p_max_new)
{
	BHNode &n = nodes[idx];
	vcopyf(p_min_new, n.p_min);
	vcopyf(p_max_new, n.p_max);
	vmeanf(n.p_min, n.p_max, n.p_c);
	n.Mass = 0;
	n.Num_Bodies = 0;
	n.BodyIdx = -1;
	memset(n.Octants, 0xFF, sizeof(n.Octants));
}

int BHTree::AllocNode(float *p_min_new, float *p_max_new)
{
	int idx = nodeCount++;
	if (idx >= (int)nodes.size()) {
		nodes.resize(nodes.size() * 2);
	}
	InitNode(idx, p_min_new, p_max_new);
	return idx;
}

void BHTree::Reset(float *p_min_new, float *p_max_new)
{
	nodeCount = 0;
	AllocNode(p_min_new, p_max_new);
}

Octant BHTree::DetermineOctant(int idx, float *p)
{
	BHNode &n = nodes[idx];
	int oct = 0;
	if (p[2] >= n.p_c[2]) {
		if (p[1] >= n.p_c[1]) {
			oct = (p[0] >= n.p_c[0]) ? FNE : FNW;
		} else {
			oct = (p[0] >= n.p_c[0]) ? FSE : FSW;
		}
	} else {
		if (p[1] >= n.p_c[1]) {
			oct = (p[0] >= n.p_c[0]) ? BNE : BNW;
		} else {
			oct = (p[0] >= n.p_c[0]) ? BSE : BSW;
		}
	}
	return (Octant)oct;
}

int BHTree::CreateChild(int idx, Octant Oct)
{
	BHNode &n = nodes[idx];
	float p_min_new[3], p_max_new[3];

	switch (Oct)
	{
		case FNE:
			vsetf(n.p_c[0], n.p_c[1], n.p_c[2], p_min_new);
			vsetf(n.p_max[0], n.p_max[1], n.p_max[2], p_max_new);
			break;
		case FNW:
			vsetf(n.p_min[0], n.p_c[1], n.p_c[2], p_min_new);
			vsetf(n.p_c[0], n.p_max[1], n.p_max[2], p_max_new);
			break;
		case FSW:
			vsetf(n.p_min[0], n.p_min[1], n.p_c[2], p_min_new);
			vsetf(n.p_c[0], n.p_c[1], n.p_max[2], p_max_new);
			break;
		case FSE:
			vsetf(n.p_c[0], n.p_min[1], n.p_c[2], p_min_new);
			vsetf(n.p_max[0], n.p_c[1], n.p_max[2], p_max_new);
			break;
		case BNE:
			vsetf(n.p_c[0], n.p_c[1], n.p_min[2], p_min_new);
			vsetf(n.p_max[0], n.p_max[1], n.p_c[2], p_max_new);
			break;
		case BNW:
			vsetf(n.p_min[0], n.p_c[1], n.p_min[2], p_min_new);
			vsetf(n.p_c[0], n.p_max[1], n.p_c[2], p_max_new);
			break;
		case BSW:
			vsetf(n.p_min[0], n.p_min[1], n.p_min[2], p_min_new);
			vsetf(n.p_c[0], n.p_c[1], n.p_c[2], p_max_new);
			break;
		case BSE:
			vsetf(n.p_c[0], n.p_min[1], n.p_min[2], p_min_new);
			vsetf(n.p_max[0], n.p_c[1], n.p_c[2], p_max_new);
			break;
	}

	int child = AllocNode(p_min_new, p_max_new);
	nodes[idx].Octants[Oct] = child;
	return child;
}

void BHTree::InsertBody(float *p, float m, int BodyIdxIn)
{
	int stack[64];
	int top = 0;
	stack[top++] = 0;
	int depth = 0;

	while (top > 0) {
		if (depth >= 60) {
			nodes[stack[top-1]].Mass += m;
			nodes[stack[top-1]].Num_Bodies++;
			break;
		}

		int idx = stack[--top];

		if (nodes[idx].Num_Bodies == 0) {
			vcopyf(p, nodes[idx].p_cm);
			nodes[idx].Mass = m;
			nodes[idx].BodyIdx = BodyIdxIn;
			nodes[idx].Num_Bodies++;
		} else if (nodes[idx].Num_Bodies == 1) {
			float existing_cm[3];
			vcopyf(nodes[idx].p_cm, existing_cm);
			float existing_mass = nodes[idx].Mass;
			int existing_body = nodes[idx].BodyIdx;

			Octant oct_existing = DetermineOctant(idx, existing_cm);
			if (nodes[idx].Octants[oct_existing] == -1) CreateChild(idx, oct_existing);

			if (vequalf(existing_cm, p)) {
				int child_idx = nodes[idx].Octants[oct_existing];
				vcopyf(existing_cm, nodes[child_idx].p_cm);
				nodes[child_idx].Mass = existing_mass + m;
				nodes[child_idx].BodyIdx = existing_body;
				nodes[child_idx].Num_Bodies = 1;
			} else {
				int child_existing = nodes[idx].Octants[oct_existing];
				vcopyf(existing_cm, nodes[child_existing].p_cm);
				nodes[child_existing].Mass = existing_mass;
				nodes[child_existing].BodyIdx = existing_body;
				nodes[child_existing].Num_Bodies = 1;

				Octant oct_new = DetermineOctant(idx, p);
				if (nodes[idx].Octants[oct_new] == -1) CreateChild(idx, oct_new);

				stack[top++] = nodes[idx].Octants[oct_new];
				nodes[idx].Num_Bodies++;
				depth++;
				continue;
			}
			nodes[idx].Num_Bodies++;
		} else {
			Octant oct = DetermineOctant(idx, p);
			if (nodes[idx].Octants[oct] == -1) CreateChild(idx, oct);
			nodes[idx].Num_Bodies++;
			stack[top++] = nodes[idx].Octants[oct];
			depth++;
		}
	}
}

void BHTree::CalcMasses()
{
	for (int idx = nodeCount - 1; idx >= 0; idx--) {
		BHNode &n = nodes[idx];
		if (n.Num_Bodies <= 1) continue;

		n.Mass = 0.0f;
		vsetf(0.0f, 0.0f, 0.0f, n.p_cm);
		float oct_cm[3];

		for (int i = 0; i < 8; i++) {
			if (n.Octants[i] != -1) {
				BHNode &child = nodes[n.Octants[i]];
				n.Mass += child.Mass;
				vscalef(child.p_cm, child.Mass, oct_cm);
				vaddf(n.p_cm, oct_cm, n.p_cm);
			}
		}

		vscalef(n.p_cm, 1.0f / n.Mass, n.p_cm);
	}
}

void BHTree::CalcMassesNode(int idx)
{
	(void)idx;
}

void BHTree::CalcAcceleration(float *p, int BodyIdxIn, float G, float r_soft, double *a)
{
	vset(0.0, 0.0, 0.0, a);

	int stack[512];
	int top = 0;
	stack[top++] = 0;

	float v[3];

	while (top > 0) {
		int idx = stack[--top];
		BHNode &n = nodes[idx];

		if (n.Num_Bodies == 1) {
			if (n.BodyIdx != BodyIdxIn) {
				vsubf(n.p_cm, p, v);
				float dsq = vmagsqsoftf(v, r_soft);
				float r3_inv = fast_r3_inv(dsq);
				vscaleaddf(v, G * n.Mass * r3_inv, a);
			}
		} else {
			vsubf(n.p_cm, p, v);
			float dsq = vmagsqf(v);
			float d = n.p_max[0] - n.p_min[0];
			if (d * d <= dsq) {
				float dsq_soft = dsq + r_soft;
				float r3_inv = fast_r3_inv(dsq_soft);
				vscaleaddf(v, G * n.Mass * r3_inv, a);
			} else {
				for (int i = 0; i < 8; i++) {
					if (n.Octants[i] != -1 && top < 510) {
						stack[top++] = n.Octants[i];
					}
				}
			}
		}
	}
}

void BHTree::GetBounds(int idx, float *BoundData)
{
	BHNode &n = nodes[idx];
	vcopyf(n.p_c, BoundData);
	BoundData[3] = n.p_max[0] - n.p_min[0];
}
