#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <cstring>
#include <GL/glut.h>
#include "disney.h"
#include "microfacet.h"
#include "levmar.h"
double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
const int NUMHEMIPOINT = 10000;
const int BLINNSAMPLE = 10000;
const int DRAWSAMPLE = 10000;
vec3 HemiPoint[NUMHEMIPOINT][2];
double Blinn[BLINNSAMPLE];
float Rb=100, Rd, C;
double p[] = { 0.5, 0.2 };
double x[BLINNSAMPLE];
float scale = 0.7;
float eye[] = { 0, 0, 10 };
float center[] = { 0, 0, 0 };
bool bPersp = false;
bool bAnim = false;
bool bWire = false;
bool bDorB = false;
float fRotate = 0;
float fDistance = 0.2f;
bool reOutput = true;
vec3 Wi = UniformSampleHemisphere(random(), random());
void reshape(int width, int height)
{
	if (height == 0)										// Prevent A Divide By Zero By
	{
		height = 1;										// Making Height Equal One
	}
	glViewport(0, 0, width, height);						// Reset The Current Viewport
	glMatrixMode(GL_PROJECTION);						// Select The Projection Matrix
	glLoadIdentity();									// Reset The Projection Matrix

	// Calculate The Aspect Ratio Of The Window
	gluPerspective(45.0f, (GLfloat)width / (GLfloat)height, 2.0f, 100.0f);
	glMatrixMode(GL_MODELVIEW);							// Select The Modelview Matrix
	glLoadIdentity();									// Reset The Modelview Matrix
}
void idle()
{
	glutPostRedisplay();
}
float maxc = 1.0;
void key(unsigned char k, int x, int y)
{
	switch (k)
	{
		case ' ': {bAnim = !bAnim; break; }
		case 'o': {bWire = !bWire; break; }
		case 'p': {bDorB = !bDorB; if (bDorB) printf("Blinn\n"); else printf("Disney\n"); reOutput = true; break; }
		case 'r': { Wi = UniformSampleHemisphere(random(), random()); printf("Wi %f %f %f\n", Wi.x, Wi.y, Wi.z); reOutput = true; maxc = 1.0; break; }
		case 'a': {eye[0] = eye[0] + fDistance; center[0] = center[0] + fDistance;
			break;
		}
		case 'd': {eye[0] = eye[0] - fDistance; center[0] = center[0] - fDistance;
			break;
		}
		case 'w': {eye[1] = eye[1] - fDistance; center[1] = center[1] - fDistance;
			break;
		}
		case 's': {eye[1] = eye[1] + fDistance; center[1] = center[1] + fDistance;
			break;
		}
		case 'z': {eye[2] = eye[2] * 0.95;
			break;
		}
		case 'c': {eye[2] = eye[2] * 1.05;
			break;
		}
		case 27:
		case 'q': {exit(0); break; }
				  //	case 'p': {bPersp = !bPersp; updateView(wWidth, wHeight); break; }
	}
}
//const float F = 0.06 + (1 - 0.06)*std::pow(1 - dot(bRec.wo, lH), 5.0);


void drawContrast(){
	MicrofacetDistribution dis(MicrofacetDistribution::EPhong, Rb);
	float color, max = 0, min = 9999999999999;
	glBegin(GL_LINES);
	glColor3f(1.0f, 1.0f, 1.0f);
	glVertex3f(0.0f, 0.0f, 0.0f);
	glVertex3f(Wi.x, Wi.y, Wi.z);
	glEnd();
	glBegin(GL_POINTS);
	for (int i = 0; i < DRAWSAMPLE;i++)
	{
		vec3 Wo = UniformSampleHemisphere(random(), random());
		vec3 Wh = normalize(Wi + Wo);
		if (Wh.z < 0)continue;
		if (bDorB)
		{
			float F = 0.06 + (1 - 0.06)*pow(1 - dot(Wo, Wh), 5.0f);
			float D = dis.eval(Wh);
			color = D * F;
		}
		else
		{				
			float D = D_GGX(Rd, Wh.z);
			float G = G_Smith(Rd, Wi.z, Wo.z);
			color = C * D * G;
		}
		
		glColor3f(color / maxc*scale + (1 - scale), color / maxc*scale + (1 - scale), color / maxc*scale + (1 - scale));
		Wo = Wo*(color / maxc*scale + (1 - scale));
		glVertex3f(Wo.x, Wo.y, Wo.z);
		if (color > max) max = color;
		if (color < min) min = color;
	}
	glEnd();
	if (max>maxc)maxc = max;
	if (max < 1.0&&abs(maxc-1.0)<1e-5) maxc = max;
	if (reOutput)
	{
		printf("%f %f %f\n", min, max, maxc);
		reOutput = false;
	}
}
void redraw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity();									// Reset The Current Modelview Matrix
	gluLookAt(eye[0], eye[1], eye[2],
		center[0], center[1], center[2],
		0, 1, 0);
	glEnable(GL_DEPTH_TEST);
	if (bWire) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	}
	else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}

	glRotatef(fRotate, 0, 1.0f, 0);			// Rotate around Y axis
	if (bAnim) fRotate += 1.0f;	
	glColor3f(1.0, 1.0, 1.0);
	glutWireCube(2.0);
	drawContrast();
	glutSwapBuffers();
	//system("pause");
	//D_GGX * G_Smith
}

void BlinnF(double *p, double *hx, int m, int n, void *adata){
	MicrofacetDistribution dis(MicrofacetDistribution::EPhong, p[0]);
	for (int i = 0; i < BLINNSAMPLE; i++)
	{
		vec3 H = normalize(HemiPoint[i][0] + HemiPoint[i][1]);
		const float F = 0.06 + (1 - 0.06)*pow(1 - dot(HemiPoint[i][1], H), 5.0f);
		hx[i] = Blinn[i] - dis.eval(H)*F;
	}
}
void DisneyF_1(double *p, double *hx, int m, int n, void *adata){
	for (int i = 0; i < BLINNSAMPLE; i++)
	{
		vec3 H = normalize(HemiPoint[i][0] + HemiPoint[i][1]);
		if (H.z < 0)continue;
		float D = D_GGX(p[0], H.z);
		float G = G_Smith(p[0], HemiPoint[i][1].z, HemiPoint[i][0].z);
		hx[i] = Blinn[i]-D*G;
	}
}
void DisneyF_2(double *p, double *hx, int m, int n, void *adata){
	for (int i = 0; i < BLINNSAMPLE; i++)
	{
		vec3 H = normalize(HemiPoint[i][0] + HemiPoint[i][1]);
		if (H.z < 0)continue;
		float D = D_GGX(p[0], H.z);
		float G = G_Smith(p[0], HemiPoint[i][1].z, HemiPoint[i][0].z);
		hx[i] = Blinn[i] - p[1] * D*G;
	}
}
void fitDfromB(){
	for (int i = 0; i < NUMHEMIPOINT; i++)
	{
		HemiPoint[i][0] = UniformSampleHemisphere(random(), random());
		HemiPoint[i][1] = UniformSampleHemisphere(random(), random());
	}
	int m = 2;
	int n = BLINNSAMPLE;
	opts[0] = LM_INIT_MU; opts[1] = 1E-15; opts[2] = 1E-15; opts[3] = 1E-20;
	opts[4] = LM_DIFF_DELTA;

	MicrofacetDistribution dis(MicrofacetDistribution::EPhong, Rb);
	for (int i = 0; i < n; i++)
	{
		vec3 H = normalize(HemiPoint[i][0] + HemiPoint[i][1]);
		const float F = 0.06 + (1 - 0.06)*pow(1 - dot(HemiPoint[i][1], H), 5.0f);
		if (H.z<0)continue;
		float D = dis.eval(H);
		Blinn[i] = D*F;
	}
	memset(x, 0, sizeof(x));
	int itmax = 1000;
	//if (p[0] < 1e-4) p[0] = 0.01;
	//if (p[1] < 1e-4) p[1] = 0.01;
//	p[0] = 0.5;
//	p[1] = 0.2;
	int ret = dlevmar_dif(DisneyF_2, p, x, m, n, itmax, opts, info, NULL, NULL, NULL);
	printf("Levenberg-Marquardt returned %d in %g iter, reason %g\nSolution: ", ret, info[5], info[6]);
	for (int i = 0; i<m; ++i)
		printf("%.7g ", p[i]);
	printf("\n\nMinimization info:\n");
	for (int i = 0; i<LM_INFO_SZ; ++i)
		printf("%g ", info[i]);
	printf("\n");
	//system("pause");
	Rd = p[0];
	C = p[1];
}
double table[10000][2];
int main(int argc, char *argv[]){
	int seq = 0;
	for (float i = 0.0; i < 100000; i+=10){
		Rb = i;
		fitDfromB();
		table[seq][0] = Rd;
		table[seq][1] = C;
		seq++;
		printf("%f\n", i);
	}
	FILE *f;
	fopen_s(&f,"DUIZHAO.txt", "a+");
	for (int i = 0; i < seq; i++)
	{
		fprintf(f, "%f,%f,", table[i][0],table[i][1]);
	}	
	fclose(f);
	srand((unsigned)time(0));
	glutInit(&argc, argv);
	glutInitWindowSize(800, 800);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	int windowHandle
		= glutCreateWindow("Simple GLUT App");
	glutDisplayFunc(redraw);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(key);
	glutIdleFunc(idle);
	glutMainLoop();
	return 0;
}
/*
void tryBlinee(){
	for (int i = 0; i < NUMHEMIPOINT; i++)
	{
		HemiPoint[i][0] = UniformSampleHemisphere(random(), random());
		HemiPoint[i][1] = UniformSampleHemisphere(random(), random());
	}
	int m = 1;
	int n = BLINNSAMPLE;
	opts[0] = LM_INIT_MU; opts[1] = 1E-15; opts[2] = 1E-15; opts[3] = 1E-20;
	opts[4] = LM_DIFF_DELTA;
	double p[] = { 1 };
	double x[BLINNSAMPLE];
	MicrofacetDistribution dis(MicrofacetDistribution::EPhong, 12.235);
	for (int i = 0; i < n; i++)
	{
		vec3 H = normalize(HemiPoint[i][0] + HemiPoint[i][1]);
		const float F = 0.06 + (1 - 0.06)*pow(1 - dot(HemiPoint[i][1], H), 5.0f);
		if (H.z<0)continue;
		float D = dis.eval(H);
		Blinn[i] = D*F;
	}
	memset(x, 0, sizeof(x));
	int itmax = 1000;
	int ret = dlevmar_dif(BlinnF, p, x, m, n, itmax, opts, info, NULL, NULL, NULL);
	printf("Levenberg-Marquardt returned %d in %g iter, reason %g\nSolution: ", ret, info[5], info[6]);
	for (int i = 0; i<m; ++i)
		printf("%.7g ", p[i]);
	printf("\n\nMinimization info:\n");
	for (int i = 0; i<LM_INFO_SZ; ++i)
		printf("%g ", info[i]);
	printf("\n");
	system("pause");
}
double F[100];
void line(double *p, double *hx, int m, int n, void *adata){
	for (int i = 0; i < n; i++)
	{
		hx[i] = F[i] - p[0] * i * 11 - p[1];
	}
}
void tryLM(){
	int m = 2;
	int n = 100;
	opts[0] = LM_INIT_MU; opts[1] = 1E-15; opts[2] = 1E-15; opts[3] = 1E-20;
	opts[4] = LM_DIFF_DELTA;
	double p[] = { 2, 5 };
	double x[100];
	for (int i = 0; i < n; i++)
	{
		double x = i * 11;
		F[i] = x*3.12 + 6.7 + random();
	}
	memset(x, 0, sizeof(x));
	int itmax = 1000;
	int ret = dlevmar_dif(line, p, x, m, n, itmax, opts, info, NULL, NULL, NULL);
	printf("Levenberg-Marquardt returned %d in %g iter, reason %g\nSolution: ", ret, info[5], info[6]);
	for (int i = 0; i<m; ++i)
		printf("%.7g ", p[i]);
	printf("\n\nMinimization info:\n");
	for (int i = 0; i<LM_INFO_SZ; ++i)
		printf("%g ", info[i]);
	printf("\n");
	system("pause");
}*/
/*
void drawBlinn(){
	MicrofacetDistribution dis(MicrofacetDistribution::EPhong, 20);
	float max = 0, min = 9999999999999;
	glBegin(GL_POINTS);
	float step = 0.05, color;
	int times = 100000;
	while (times--)
	{
		vec3 H = UniformSampleSphere(random(), random());
		if (H.z<0)continue;
		float D = dis.eval(H);
		color = D / 0.4 + 0;
		glColor3f(color, color, color);
		glVertex3f(H.x, H.y, H.z);
		if (D > max) max = D;
		if (D < min) min = D;
	}
	glEnd();
	printf("%f %f\n", min, max);
}
void drawDisney(){
	float max = 0, min = 9999999999999;
	glBegin(GL_POINTS);
	float step = 0.05, color;
	float roughness = 0.5;
	int times = 100000;
	for (float i = 0; i <= 1.0; i += 0.01)
		for (float j = 0; j <= 1.0; j += 0.01)
		{
			//vec3 H = UniformSampleHemisphere(random(), random());
			vec3 H = UniformSampleHemisphere(i, j);
			if (H.z<0)continue;
			float D = D_GGX(roughness, dot(H, vec3(0, 0, 1)));
			//	float G = G_Smith(1.0,)
			color = D / 5 + 0;
			//glColor3f(color, color, color);
			glVertex3f(H.x, H.y, H.z);
			if (D > max) max = D;
			if (D < min) min = D;
		}
	glEnd();
	printf("%f %f\n", min, max);
}*/