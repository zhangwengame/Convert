#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <cstring>
#include <GL/glut.h>
#include "disney.h"
#include "microfacet.h"
#include "levmar.h"
double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
const int NUMHEMIPOINT = 10010;
const int BLINNSAMPLE = 10000;
const int DRAWSAMPLE = 10000;
const int RSAMPLE = 1000;
vec3 HemiPoint[NUMHEMIPOINT][2];
double Importance[NUMHEMIPOINT];
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
bool bforh = true;
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
float maxb = 0.1, maxd = 0.1;
void key(unsigned char k, int x, int y)
{
	switch (k)
	{
		case ' ': {bAnim = !bAnim; break; }
		case 'o': {bWire = !bWire; break; }
		case 'f': {bforh = !bforh; if (bforh) printf("combine\n"); else printf("seperate\n"); reOutput = true; maxc = 1.0; maxb = maxd = 0.1; break;  }
		case 'p': {bDorB = !bDorB; if (bDorB) printf("Blinn\n"); else printf("Disney\n"); reOutput = true; break; }
		case 'r': { Wi = UniformSampleHemisphere(random(), random()); printf("Wi %f %f %f\n", Wi.x, Wi.y, Wi.z); reOutput = true; maxc = 1.0; maxb = maxd = 0.1;  break; }
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
vec3 ImportanceSampleBlinn(double *pdf, double n){
	float phi = 2.0f * PI * random();
	float theta = acos(pow((double)1.0f - random(), (double)1.0f / (n + 1.0f)));
	*pdf = (n + 1) / (2 * PI)*pow((double)cos(theta), (double)n)*sin(theta);
	return vec3(sin(theta) * cos(phi),
		sin(theta) * sin(phi),
		cos(theta));

}

void drawContrast(){
	int ddebug = 0;
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
		double pdf;
		vec3 Wh = ImportanceSampleBlinn(&pdf,Rb);
		vec3 Wo = normalize(Wh*dot(Wi, Wh) * 2 - Wi);
		if (dot(Wo, vec3(0, 0, 1)) < 0) continue;
		if (Wh.z < 0)continue;
		if (bDorB)
		{
			float F = 0.06 + (1 - 0.06)*pow(1 - dot(Wo, Wh), 5.0f);
			float D = dis.eval(Wh);
			color = D * F;
			if (color > maxb) maxb = color;
		}
		else
		{				
			float D = D_GGX(Rd, Wh.z);
			float G = G_Smith(Rd, Wi.z, Wo.z);
			color = C * D * G;
			if (color>1000000)
			{
				ddebug = 1;
				float D = D_GGX(Rd, Wh.z);
				float G = G_Smith(Rd, Wi.z, Wo.z);
			}
			if (color > maxd) maxd = color;
	/*		if (color > maxb)
			{
				ddebug = 1;
				float D = D_GGX(Rd, Wh.z);
				float G = G_Smith(Rd, Wi.z, Wo.z);
			}*/
				
		}
		float cscale;
		if (bforh)
			cscale = maxc;
		else
			cscale = bDorB ? maxb : maxd;
		glColor3f(color / cscale*scale + (1 - scale), color / cscale*scale + (1 - scale), color / cscale*scale + (1 - scale));
		Wo = Wo*(color / cscale*scale + (1 - scale));
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
	double r = p[0] >= 0.0f ? p[0] : 0.0f;
	for (int i = 0; i < BLINNSAMPLE; i++)
	{
		vec3 H = normalize(HemiPoint[i][0] + HemiPoint[i][1]);
		if (H.z < 0)continue;
		float D = D_GGX(r, H.z);
		float G = G_Smith(r, HemiPoint[i][1].z, HemiPoint[i][0].z);
		hx[i] = (Blinn[i] - p[1] * D*G)/Importance[i];
	}
}

void fitDfromB_Force_N(){
	int ddebug = 0;
	for (int i = 0; i < NUMHEMIPOINT; i++)
	{
		vec3 Wo, Wh;
		float theta = 0.50*PI*(0.0001*i + 0.00005);
		HemiPoint[i][0] = vec3(sin(theta),0,cos(theta));
		do{
			Wh = ImportanceSampleBlinn(&Importance[i], Rb);
			Wo = normalize(Wh*dot(HemiPoint[i][0], Wh) * 2 - HemiPoint[i][0]);
			Importance[i] /= (4 * dot(Wo, Wh));
			HemiPoint[i][1] = Wo;
		} while (Wo.z < 0.0f);
	}
	MicrofacetDistribution dis(MicrofacetDistribution::EPhong, Rb);
	for (int i = 0; i < BLINNSAMPLE; i++)
	{
		vec3 H = normalize(HemiPoint[i][0] + HemiPoint[i][1]);
		const float F = 0.06 + (1 - 0.06)*pow(1 - dot(HemiPoint[i][1], H), 5.0f);
		if (H.z<0)continue;
		float D = dis.eval(H);
		Blinn[i] = D*F*HemiPoint[i][0].z;

	}
	float diff[1001];		
	
	float min = FLT_MAX,minRd,minScale;
	for (int i = 1; i <= RSAMPLE; i++){
		float iRd = 1.0/RSAMPLE*i;
		float sum = 0,ascale=1.0;
		double sumfafb = 0.0, sumfb2 = 0.0;
		for (int j = 0; j < BLINNSAMPLE; j++)
		{			
			vec3 H = normalize(HemiPoint[j][0] + HemiPoint[j][1]);
			if (H.z < 0)continue;
			float D = D_GGX(iRd, H.z);
			float G = G_Smith(iRd, HemiPoint[j][1].z, HemiPoint[j][0].z);
			float Disney = D * G* HemiPoint[j][0].z;			
			sumfafb += Blinn[j] * Disney;
			sumfb2 += Disney*Disney;
		}
		ascale = sumfafb / sumfb2;
		for (int j = 0; j < BLINNSAMPLE; j++)
		{
			vec3 H = normalize(HemiPoint[j][0] + HemiPoint[j][1]);
			if (H.z < 0)continue;
			float D = D_GGX(iRd, H.z);
			float G = G_Smith(iRd, HemiPoint[j][1].z, HemiPoint[j][0].z);
			float Disney = D * G* HemiPoint[j][0].z*ascale;
			sum += (Disney - Blinn[j]) * (Disney - Blinn[j]) / (Importance[j] * Importance[j]);
		}		
		diff[i] = sum;
		if (sum < min)
		{
			min = sum;
			minRd = iRd;
			minScale = ascale;
		}
	}
	
	
	Rd = minRd;
	C = minScale;
}
void fitDfromB_Force();
void fitDfromB_LM();
double table[10000][2];
int main(int argc, char *argv[]){
	int seq = 0;
	float begin = 1500000.0f;
	char name[100];
	for (float i =begin; i < begin+5.0; i+=10){
	//for (float i = 0; i <= 100000; i += 10){
		Rb = i;
		fitDfromB_Force_N();
		table[seq][0] = Rd;
		table[seq][1] = C;
		seq++;
		printf("%f Rd %f C %f\n", i,Rd,C);
		if (int(i) % 10000 == 0)
		{
			sprintf_s(name, "DUIZHAO_%d.txt", i);
			FILE *f;
			fopen_s(&f, name, "w");
			for (int i = 0; i < seq; i++)
			{
				fprintf(f, "%f,%f,", table[i][0], table[i][1]);
			}
			fclose(f);
		}

	}

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
int gui = 0, lianguijs = 0;
int liangui = 0;
double l1, l0;
void fitDfromB_Force(){
	int ddebug = 0;
	for (int i = 0; i < NUMHEMIPOINT; i++)
	{
		vec3 Wo, Wh;
		float theta = 0.50*PI*(0.0001*i + 0.00005);
		HemiPoint[i][0] = vec3(sin(theta), 0, cos(theta));
		do{
			Wh = ImportanceSampleBlinn(&Importance[i], Rb);
			Wo = normalize(Wh*dot(HemiPoint[i][0], Wh) * 2 - HemiPoint[i][0]);
			Importance[i] /= (4 * dot(Wo, Wh));
			HemiPoint[i][1] = Wo;
		} while (Wo.z < 0.0f);
	}
	MicrofacetDistribution dis(MicrofacetDistribution::EPhong, Rb);
	float maxb = 0;
	int maxi;
	for (int i = 0; i < BLINNSAMPLE; i++)
	{
		vec3 H = normalize(HemiPoint[i][0] + HemiPoint[i][1]);
		const float F = 0.06 + (1 - 0.06)*pow(1 - dot(HemiPoint[i][1], H), 5.0f);
		if (H.z<0)continue;
		float D = dis.eval(H);
		Blinn[i] = D*F*HemiPoint[i][0].z;
		if (Blinn[i] > maxb)
		{
			maxb = Blinn[i];
			maxi = i;
		}
	}
	float diff[1001];
	float min = FLT_MAX, minRd;
	for (int i = 1; i <= RSAMPLE; i++){
		float iRd = 1.0 / RSAMPLE*i;
		float sum = 0, ascale = 1.0;
		if (i == 300)
			ddebug = 1;
		{
			vec3 H = normalize(HemiPoint[maxi][0] + HemiPoint[maxi][1]);
			float D = D_GGX(iRd, H.z);
			float G = G_Smith(iRd, HemiPoint[maxi][1].z, HemiPoint[maxi][0].z);
			ascale = maxb / (D*G);
		}
		//MicrofacetDistribution dis(MicrofacetDistribution::EPhong, iRd);
		for (int j = 0; j < BLINNSAMPLE; j++)
		{
			/*	vec3 H = normalize(HemiPoint[j][0] + HemiPoint[j][1]);
			const float F = 0.06 + (1 - 0.06)*pow(1 - dot(HemiPoint[j][1], H), 5.0f);
			if (H.z<0)continue;
			float D = dis.eval(H);
			float Blinn_ = D*F*HemiPoint[j][0].z;
			sum += (Blinn_ - Blinn[j]) * (Blinn_ - Blinn[j]);*/
			vec3 H = normalize(HemiPoint[j][0] + HemiPoint[j][1]);
			if (H.z < 0)continue;
			float D = D_GGX(iRd, H.z);
			float G = G_Smith(iRd, HemiPoint[j][1].z, HemiPoint[j][0].z);
			float Disney = D * G* HemiPoint[j][0].z*ascale;
			sum += (Disney - Blinn[j]) * (Disney - Blinn[j]);// / (Importance[j] * Importance[j]);

		}
		diff[i] = sum;
		if (sum < min)
		{
			min = sum;
			minRd = iRd;
		}
	}
	//	minRd = 0.3;
	double sumfafb = 0.0, sumfb2 = 0.0;
	MicrofacetDistribution diss(MicrofacetDistribution::EPhong, minRd);
	for (int i = 0; i < BLINNSAMPLE; i++){
		/*	vec3 H = normalize(HemiPoint[i][0] + HemiPoint[i][1]);
		const float F = 0.06 + (1 - 0.06)*pow(1 - dot(HemiPoint[i][1], H), 5.0f);
		if (H.z<0)continue;
		float D = diss.eval(H);
		float Blinn_ = D*F*HemiPoint[i][0].z;
		sumfafb +=  Blinn[i] * Blinn_;
		sumfb2 += Blinn_*Blinn_;*/
		vec3 H = normalize(HemiPoint[i][0] + HemiPoint[i][1]);
		if (H.z < 0)continue;
		float D = D_GGX(minRd, H.z);
		float G = G_Smith(minRd, HemiPoint[i][1].z, HemiPoint[i][0].z);
		float Disney = D*G*HemiPoint[i][0].z;

		sumfafb += Blinn[i] * Disney;
		sumfb2 += Disney*Disney;
	}
	Rd = minRd;
	C = sumfafb / sumfb2;
}
void fitDfromB_LM(){
F:
	float minE = FLT_MAX, minL0, minL1;
	for (int i = 0; i < NUMHEMIPOINT; i++)
	{
		vec3 Wo, Wh;
		HemiPoint[i][0] = UniformSampleHemisphere(random(), random());
		do{
			Wh = ImportanceSampleBlinn(&Importance[i], Rb);
			Wo = normalize(Wh*dot(HemiPoint[i][0], Wh) * 2 - HemiPoint[i][0]);
			Importance[i] /= (4 * dot(Wo, Wh));
			HemiPoint[i][1] = Wo;
		} while (Wo.z < 0.0f);
	}
	int m = 2;
	int n = BLINNSAMPLE;
	opts[0] = LM_INIT_MU; opts[1] = 1E-15; opts[2] = 1E-15; opts[3] = 1E-20;
	opts[4] = LM_DIFF_DELTA;
	float maxb = 0;
	int maxi = 0;
	MicrofacetDistribution dis(MicrofacetDistribution::EPhong, Rb);
	for (int i = 0; i < n; i++)
	{
		vec3 H = normalize(HemiPoint[i][0] + HemiPoint[i][1]);
		const float F = 0.06 + (1 - 0.06)*pow(1 - dot(HemiPoint[i][1], H), 5.0f);
		if (H.z<0)continue;
		float D = dis.eval(H);
		Blinn[i] = D*F;
		if (Blinn[i] > maxb) { maxb = Blinn[i]; maxi = i; };
	}
	memset(x, 0, sizeof(x));
	int itmax = 1000;
	//if (p[0] < 1e-4) p[0] = 0.01;
	//if (p[1] < 1e-4) p[1] = 0.01;
	{
		p[0] = random()*0.5;
		vec3 H = normalize(HemiPoint[maxi][0] + HemiPoint[maxi][1]);
		float D = D_GGX(p[0], H.z);
		float G = G_Smith(p[0], HemiPoint[maxi][1].z, HemiPoint[maxi][0].z);
		p[1] = maxb / (D*G);
	}
	int ret = dlevmar_dif(DisneyF_2, p, x, m, n, itmax, opts, info, NULL, NULL, NULL);
	printf("Levenberg-Marquardt returned %d in %g iter, reason %g\nSolution: ", ret, info[5], info[6]);
	for (int i = 0; i<m; ++i)
		printf("%.7g ", p[i]);
	printf("\n\nMinimization info:\n");
	for (int i = 0; i<LM_INFO_SZ; ++i)
		printf("%g ", info[i]);
	printf("\n");
	printf("gui le %d ci\n", gui);
	//system("pause");

	if (p[0] < 0 || p[1] < 0 || p[0]>1)
	{
		gui++;
		liangui++;
		goto F;
	}
	if (info[1] < minE){
		minE = info[1];
		minL0 = p[0];
		minL1 = p[1];
	}
	if (info[1] > 5000/*||abs((p[0]-l0)/l0)>0.3||abs((p[1]-l1)/l1)>0.8*/)
	{

		gui++;
		liangui++;
		if (liangui<20)
			goto F;
	}
	printf("lianguijs %d\n", lianguijs);
	liangui = 0;
	Rd = minL0;
	C = minL1;
}