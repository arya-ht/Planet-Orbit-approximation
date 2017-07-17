#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#pragma warning(disable:4996)

#define EPSILON 0.000001
#define MAXT  1.509*(24 * 60 * 60)

#define G 6.67 *pow(10,-11)
#define M1 .08*1.989*pow(10,30)
#define M2 .85*5.9722*pow(10,24)
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
/*
	*Ignore this... I was trying to make a generalized equation, but then I realized it became much more
	*resource demanding than the regular function calling. This is pretty close to being done though.
	*#define R1(u) u
	*#define R2(u) u
	*define R12 (xn,xn2,yn,yn2) sqrt(pow((xn-xn2),2)+pow((yn-yn2),2))
	*#define DR1(t,r1,r2,xn,xn2,yn,yn2) - ((G*M2(R1(xn)-R2(xn2)))/pow(R12(xn,xn2,yn,yn2),3))
	*#define DR2(t,r1,r2,xn,xn2,yn,yn2) - ((G*M1(R2(xn2)-R2(xn)))/pow(R12(xn,xn2,yn,yn2),3))
	*
	*#define US(t,x,r1,r2,rd1,rd2,dt) x + (dt*0.5)*DR1(t,r1,r2,r12,xn,xn2,yn,yn2);
	*#define US(t,x,r1,r2,rd1,rd2,dt) x + (dt*0.5)*DR2(t,r1,r2,r12,xn,xn2,yn,yn2));
	*#define US(t,x,r1,r2,rd1,rd2,dt) x + ((dt*0.5)*R1());
	*
	*#define USS(t,x,r1,r2,rd1,rd2,dt) x + dt*(0.5 * (t + 0.5*dt, US(t,x,r1,r2,rd1,rd2,dt)));
	*#define USSS(t,x,r1,r2,rd1,rd2,dt) x + dt*(t, USS(t,x,r1,r2,rd1,rd2,dt));
	*#define UN1(t,x,r1,r2,rd1,rd2,dt)  x + (dt / 6)*((t,x) + 2 * (t + 0.5*dt, US(t,x,r1,r2,rd1,rd2,dt), + 2 * (t + (1 / 2)*dt, USS(t,x,r1,r2,rd1,rd2,dt)) + (t, USSS(t,x,r1,r2,rd1,rd2,dt)));
	*/
double f1(double);
double f2(double, double);
double f3(double, double);
double r12(double, double, double, double);
void calculateOrbit(double, double, double, double, double, double, double, double, double, double);

int main() {
	double t = 0.0, dt = 60.0, xn = 0, xn2 = 1.9889*pow(10, 8), yn = 0, yn2 = -1.651*pow(10, 9),
		vn = 0, v2 = 79348.8, un = 0, u2 = 9078.256;
	calculateOrbit(t, dt, xn, xn2, yn, yn2, vn, v2, un, u2);
	return 0;
}
void calculateOrbit(double t, double dt, double xn, double xn2, double yn, double yn2, double vn, double v2, double un, double u2) {
	double r = r12(xn, xn2, yn, yn2), r1, r2, r3, xs, xss, xsss, xs2,
		xss2, xsss2, ys, yss, ysss, ys2, yss2, ysss2, vs, vss, vsss, vs2,
		vss2, vsss2, us, uss, usss, us2, uss2, usss2, xn1, yn1, yn21, vn1,
		un1, xn21, v21, u21, TOL, delta;

	FILE *f = fopen("7_3.csv", "w");

	printf("t: %.lf\tx1: %e\ty1: %e\tx2: %e\ty2: %e\n", t, xn, yn, xn2, yn2);
	fprintf(f, "Time, X1, Y1, X2, Y2\n%lf, %.15lf, %.15lf, %.15lf, %.15lf", t, xn, yn, xn2, yn2);

	while (t < MAXT) {
		dt = MIN(dt, MAXT - t);

		r = r12(xn, xn2, yn, yn2);
		xs = xn + .5*dt*f1(vn);
		vs = vn + .5*dt*f2(xn, xn2) / r;
		ys = yn + .5*dt*f1(un);
		us = un + .5*dt*f2(yn, yn2) / r;
		xs2 = xn2 + .5*dt*f1(v2);
		vs2 = v2 + .5*dt*f3(xn, xn2) / r;
		ys2 = yn2 + .5*dt*f1(u2);
		us2 = u2 + .5*dt*f3(yn, yn2) / r;

		r1 = r12(xs, xs2, ys, ys2);
		xss = xn + .5*dt*f1(vs);
		vss = vn + .5*dt*f2(xs, xs2) / r1;
		yss = yn + .5*dt*f1(us);
		uss = un + .5*dt*f2(ys, ys2) / r1;
		xss2 = xn2 + .5*dt*f1(vs2);
		vss2 = v2 + .5*dt*f3(xs, xs2) / r1;
		yss2 = yn2 + .5*dt*f1(us2);
		uss2 = u2 + .5*dt*f3(ys, ys2) / r1;

		r2 = r12(xss, xss2, yss, yss2);
		xsss = xn + dt*f1(vss);
		vsss = vn + dt*f2(xss, xss2) / r2;
		ysss = yn + dt*f1(uss);
		usss = un + dt*f2(yss, yss2) / r2;
		xsss2 = xn2 + dt*f1(vss2);
		vsss2 = v2 + dt*f3(xss, xss2) / r2;
		ysss2 = yn2 + dt*f1(uss2);
		usss2 = u2 + dt*f3(yss, uss2) / r2;

		r3 = r12(xsss, xsss2, ysss, ysss2);
		xn1 = xn + dt / 6 * (f1(vn) + 2 * f1(vs) + 2 * f1(vss) + f1(vsss));
		vn1 = vn + dt / 6 * (f2(xn, xn2) / r + 2 * f2(xs, xs2) / r1 + 2 * f2(xss, xss2) / r2 + f2(xsss, xsss2) / r3);
		yn1 = yn + dt / 6 * (f1(un) + 2 * f1(us) + 2 * f1(uss) + f1(usss));
		un1 = un + dt / 6 * (f2(yn, yn2) / r + 2 * f2(ys, ys2) / r1 + 2 * f2(yss, yss2) / r2 + f2(ysss, ysss2) / r3);
		xn21 = xn2 + dt / 6 * (f1(v2) + 2 * f1(vs2) + 2 * f1(vss2) + f1(vsss2));
		v21 = v2 + dt / 6 * (f3(xn, xn2) / r + 2 * f3(xs, xs2) / r1 + 2 * f3(xss, xss2) / r2 + f3(xsss, xsss2) / r3);
		yn21 = yn2 + dt / 6 * (f1(u2) + 2 * f1(us2) + 2 * f1(uss2) + f1(usss2));
		u21 = u2 + dt / 6 * (f3(yn, yn2) / r + 2 * f3(ys, ys2) / r1 + 2 * f3(yss, yss2) / r2 + f3(ysss, ysss2) / r3);
		
		TOL = fabs(xn1 - xn)*dt;
		delta = 0.84*(pow(EPSILON / TOL, 0.25));

		if (TOL <= EPSILON) {
			xn = xn1; vn = vn1; yn = yn1; un = un1; xn2 = xn21; v2 = v21; u2 = u21; yn2 = yn21;
			dt *= delta;
			printf("Time = %lf, X = %e, Y = %e, X2 = %e, Y2 = %e\n", t, xn1, yn1, xn21, yn21);
			fprintf(f, "\n%lf, %.15lf, %.15lf, %.15lf, %.15lf", t, xn1, yn1, xn21, yn21);
			t += dt;
		}
		else {dt *= delta;}
	}
	fclose(f);
}

double f1(double v) {
	return v;
}
double f2(double r1, double r2) {
	return -1 * G*M2*(r1 - r2);
}
double f3(double r1, double r2) {

	return -1 * G*M1*(r2 - r1);
}
double r12(double xn, double xn2, double yn, double yn2) {
	return pow((pow((xn - xn2), 2) + pow((yn - yn2), 2)), 1.5);
}
