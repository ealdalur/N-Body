#include "BHTree.h"

BHTree::BHTree()
{
	nodes.resize(512000);
	nodeCount = 0;
}

void BHTree::InitNode(int idx, double *p_min_new, double *p_max_new)
{
	BHNode &n = nodes[idx];
	vcopy(p_min_new, n.p_min);
	vcopy(p_max_new, n.p_max);
	vmean(n.p_min, n.p_max, n.p_c);
	vcopy(n.p_c, n.p_cm);
	n.Mass = 0;
	n.Num_Bodies = 0;
	n.BodyIdx = -1;
	for (int i = 0; i < 8; i++) n.Octants[i] = -1;
}

int BHTree::AllocNode(double *p_min_new, double *p_max_new)
{
	int idx = nodeCount++;
	if (idx >= (int)nodes.size()) {
		nodes.resize(nodes.size() * 2);
	}
	InitNode(idx, p_min_new, p_max_new);
	return idx;
}

void BHTree::Reset(double *p_min_new, double *p_max_new)
{
	nodeCount = 0;
	AllocNode(p_min_new, p_max_new);
}

Octant BHTree::DetermineOctant(int idx, double *p)
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
	double p_min_new[3], p_max_new[3];

	switch (Oct)
	{
		case FNE:
			vset(n.p_c[0], n.p_c[1], n.p_c[2], p_min_new);
			vset(n.p_max[0], n.p_max[1], n.p_max[2], p_max_new);
			break;
		case FNW:
			vset(n.p_min[0], n.p_c[1], n.p_c[2], p_min_new);
			vset(n.p_c[0], n.p_max[1], n.p_max[2], p_max_new);
			break;
		case FSW:
			vset(n.p_min[0], n.p_min[1], n.p_c[2], p_min_new);
			vset(n.p_c[0], n.p_c[1], n.p_max[2], p_max_new);
			break;
		case FSE:
			vset(n.p_c[0], n.p_min[1], n.p_c[2], p_min_new);
			vset(n.p_max[0], n.p_c[1], n.p_max[2], p_max_new);
			break;
		case BNE:
			vset(n.p_c[0], n.p_c[1], n.p_min[2], p_min_new);
			vset(n.p_max[0], n.p_max[1], n.p_c[2], p_max_new);
			break;
		case BNW:
			vset(n.p_min[0], n.p_c[1], n.p_min[2], p_min_new);
			vset(n.p_c[0], n.p_max[1], n.p_c[2], p_max_new);
			break;
		case BSW:
			vset(n.p_min[0], n.p_min[1], n.p_min[2], p_min_new);
			vset(n.p_c[0], n.p_c[1], n.p_c[2], p_max_new);
			break;
		case BSE:
			vset(n.p_c[0], n.p_min[1], n.p_min[2], p_min_new);
			vset(n.p_max[0], n.p_c[1], n.p_c[2], p_max_new);
			break;
	}

	int child = AllocNode(p_min_new, p_max_new);
	nodes[idx].Octants[Oct] = child;
	return child;
}

void BHTree::InsertBody(double *p, double m, int BodyIdxIn)
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
			vcopy(p, nodes[idx].p_cm);
			nodes[idx].Mass = m;
			nodes[idx].BodyIdx = BodyIdxIn;
			nodes[idx].Num_Bodies++;
		} else if (nodes[idx].Num_Bodies == 1) {
			double existing_cm[3];
			vcopy(nodes[idx].p_cm, existing_cm);
			double existing_mass = nodes[idx].Mass;
			int existing_body = nodes[idx].BodyIdx;

			Octant oct_existing = DetermineOctant(idx, existing_cm);
			if (nodes[idx].Octants[oct_existing] == -1) CreateChild(idx, oct_existing);

			if (vequal(existing_cm, p)) {
				int child_idx = nodes[idx].Octants[oct_existing];
				vcopy(existing_cm, nodes[child_idx].p_cm);
				nodes[child_idx].Mass = existing_mass + m;
				nodes[child_idx].BodyIdx = existing_body;
				nodes[child_idx].Num_Bodies = 1;
			} else {
				int child_existing = nodes[idx].Octants[oct_existing];
				vcopy(existing_cm, nodes[child_existing].p_cm);
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
	CalcMassesNode(0);
}

void BHTree::CalcMassesNode(int idx)
{
	BHNode &n = nodes[idx];
	if (n.Num_Bodies <= 1) return;

	n.Mass = 0.0;
	vset(0.0, 0.0, 0.0, n.p_cm);
	double oct_cm[3];

	for (int i = 0; i < 8; i++) {
		if (n.Octants[i] != -1) {
			CalcMassesNode(n.Octants[i]);
			BHNode &child = nodes[n.Octants[i]];
			n.Mass += child.Mass;
			vscale(child.p_cm, child.Mass, oct_cm);
			vadd(n.p_cm, oct_cm, n.p_cm);
		}
	}

	vscale(n.p_cm, 1.0 / n.Mass, n.p_cm);
}

void BHTree::CalcAcceleration(double *p, int BodyIdxIn, double G, double r_soft, double *a)
{
	vset(0.0, 0.0, 0.0, a);

	int stack[512];
	int top = 0;
	stack[top++] = 0;

	double v[3];

	while (top > 0) {
		int idx = stack[--top];
		BHNode &n = nodes[idx];

		if (n.Num_Bodies == 1) {
			if (n.BodyIdx != BodyIdxIn) {
				vsub(n.p_cm, p, v);
				double dsq = vmagsqsoft(v, r_soft);
				double r3_inv = 1.0 / (dsq * sqrt(dsq));
				vscaleadd(v, G * n.Mass * r3_inv, a);
			}
		} else {
			vsub(n.p_cm, p, v);
			double dsq = vmagsq(v);
			double d = n.p_max[0] - n.p_min[0];
			if (d * d <= dsq) {
				double dsq_soft = dsq + r_soft;
				double r3_inv = 1.0 / (dsq_soft * sqrt(dsq_soft));
				vscaleadd(v, G * n.Mass * r3_inv, a);
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

void BHTree::GetBounds(int idx, double *BoundData)
{
	BHNode &n = nodes[idx];
	vcopy(n.p_c, BoundData);
	BoundData[3] = n.p_max[0] - n.p_min[0];
}
