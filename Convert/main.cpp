#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <cstring>
#include <GL/glut.h>
#include "disney.h"
#include "microfacet.h"
#include "levmar.h"
float eye[] = { 0, 0, 10 };
float center[] = { 0, 0, 0 };
bool bPersp = false;
bool bAnim = false;
bool bWire = false;
float fRotate = 0;
float fDistance = 0.2f;
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
void key(unsigned char k, int x, int y)
{
	switch (k)
	{
		case ' ': {bAnim = !bAnim; break; }
		case 'o': {bWire = !bWire; break; }

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
	while (times--)
	{
		vec3 H = UniformSampleSphere(random(), random());
		if (H.z<0)continue;
		float D = D_GGX(roughness, dot(H, vec3(0, 0, 1)));
		//	float G = G_Smith(1.0,)
		color = D / 5 + 0;
		glColor3f(color, color, color);
		glVertex3f(H.x, H.y, H.z);
		if (D > max) max = D;
		if (D < min) min = D;
	}
	glEnd();
	printf("%f %f\n", min, max);
}
void redraw()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity();									// Reset The Current Modelview Matrix
	gluLookAt(eye[0], eye[1], eye[2],
		center[0], center[1], center[2],
		0, 1, 0);
	glEnable(GL_DEPTH_TEST);
/*	glEnable(GL_LIGHTING);
	GLfloat white[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat light_pos[] = { 5, 5, 5, 1 };
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, white);
	glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, white);
	glLightfv(GL_LIGHT0, GL_AMBIENT, white);
	glEnable(GL_LIGHT0);*/

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
	drawDisney();
	glutSwapBuffers();
	
	//D_GGX * G_Smith
}
double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
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
	double p[] = { 3,6 };
	double x[100];
	for (int i = 0; i < n; i++)
	{
		double x = i*11;
		F[i] = x*3.12 + 6.7;
	}
	memset(x, sizeof(x), 0);
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
}
int main(int argc, char *argv[]){
	tryLM();
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