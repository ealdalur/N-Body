#define _USE_MATH_DEFINES

#include <stdlib.h>
#include <math.h>

double drand();
void vrand(double *v);
void vset(double x, double y, double z, double *v);
void vcopy(double *v, double *vr);
double vmag(double *v);
double vmagsoft(double *v, double r_soft);
void vscale(double *v, double s, double *vr);
void vnorm(double *v);
void vneg(double *v);
void vadd(double *v1, double *v2, double *vr);
void vscaleadd(double *v, double s, double *vr);
void vsub(double *v1, double *v2, double *vr);
double vdot(double *v1, double *v2);
void vcross(double *v1, double *v2, double *vr);
void vmean(double *v1, double *v2, double *vr);
bool vequal(double *v1, double *v2);
void vmin(double *v1, double *v2, double *vr);
void vmax(double *v1, double *v2, double *vr);
double vmaxel(double *v);
double vmaxsub(double *v1, double *v2);
void vmaxboundcube(double *p_min_in, double *p_max_in, double *p_min_out, double *p_max_out);
void vrot_zy(double *v, double angle, double *vr);