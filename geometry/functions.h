#ifndef __FUNCTION_H__
#define __FUNCTION_H__

#include <math.h>

#include <iostream>
#include <vector>
using namespace std;

#define RADIANS_TO_DEGREES(r) (180.0 * ((r) / M_PI))
#define DEGREES_TO_RADIANS(r) (M_PI/180.0 * (r))

void EigenSort(vector<float> & values, vector<vector<float> > & vectors);
int Jacobi (vector<vector<float> > & q, vector<float> & values, vector<vector<float> > & vectors);

void TRTtransform(vector<vector<float> > & pts, vector<float> & t1, vector<vector<float> > & rot, vector<float> & t2);
void TRTtransform1(vector<float> & pt, vector<float> & t1, vector<vector<float> > & rot, vector<float> & t2);

float findSuperpositionTransform( vector<vector<float> > & from, vector<vector<float> > & onto, vector<float> & t1, vector<vector<float> > & rot, vector<float> & t2 );

void sphereVolSample(vector<float>& cen, float rad, vector<float>& sample);
float dihedDiff(float d0, float d1);
float withinPlusMinus180(float num);

void findTriangleBaseHeight(float a, float b, float c, float &base, float &perp);

float calcDist(vector<float> & A, vector<float> & B);
float calcAngle(vector<float> & A, vector<float> & B, vector<float> & C);
float calcAngle(vector<float> & A, vector<float> & B);
float calcDihed(vector<float> &a, vector<float> &b, vector<float> &c, vector<float> &d);

void find4thPoint(vector<float>& p4,
	vector<float>& p1, vector<float>& p2, vector<float>& p3,
	float dist, float ang, float dihed);

float magnitude(vector<float>& p);

void normalize_point(vector<float>& n, vector<float>& p);

void linear_combination(vector<float>& pt, float a, vector<float>& A, float b, vector<float>& B);

void findYZgivenX( vector<float> & X, vector<float> & Y, vector<float> & Z );
int findSphereSphereAngleIntx(float r1, vector<float> & c1, float r2, vector<float> & c2, vector<float> & q, float desiredAngle, vector<float>& p1, vector<float>& p2);

void randomNormalVector(vector<float> & v);

void rotate(vector<vector<float> > & rotOperator, vector<float>& pt);
void findRotnOperator(vector<vector<float> > & op, vector<float> & axis, float angle);

void cross_product(vector<float>& p, const vector<float>& a, const vector<float>& b);

#endif //__FUNCTION_H__
