#if !defined(__SHARE_H_)
#define __SHARE_H_
#include <math.h>
#include <float.h>
#define M_PI 3.1415926535897932f
#define INV_TWOPI     0.15915494309189533577f
float max(const float&a, const float &b){
	return a > b ? a : b;
}
struct vec3{
	vec3() :x(0.0), y(0.0), z(0.0){}
	vec3(float x, float y, float z) :x(x), y(y), z(z){}
	float x, y, z;
};
inline float dot(const vec3 &a, const  vec3 &b){
	return a.x*b.x + a.y*b.y + a.z*b.z;
}
inline vec3 operator +(const vec3 &a, const  vec3 &b){
	return vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}
inline vec3 operator -(const vec3 &a, const  vec3 &b){
	return vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}
inline vec3 operator *(const float &a, const  vec3 &b){
	return vec3(a*b.x, a*b.y, a*b.z);
}
inline vec3 cross(const vec3 &a, const vec3 &b){
	return vec3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}
inline vec3 normalize(const vec3 &a){
	float length = sqrtf(dot(a, a));
	return vec3(a.x / length, a.y / length, a.z / length);
}
void buildOrthonormalBasis(vec3& omega_1, vec3& omega_2, const vec3& omega_3)
{
	if (omega_3.z < -0.9999999f)
	{
		omega_1 = vec3(0.0f, -1.0f, 0.0f);
		omega_2 = vec3(-1.0f, 0.0f, 0.0f);
	}
	else {
		const float a = 1.0f / (1.0f + omega_3.z);
		const float b = -omega_3.x*omega_3.y*a;
		omega_1 = vec3(1.0f - omega_3.x*omega_3.x*a, b, -omega_3.x);
		omega_2 = vec3(b, 1.0f - omega_3.y*omega_3.y*a, -omega_3.y);
	}
}
float rcp(float x){
	if (x != 0)
		return 1.0 / x;
	else
		return FLT_MAX;
}
float random(){
	float r = (static_cast <float> (rand())*(static_cast <float> (RAND_MAX)+1) + static_cast <float> (rand())) /
		((static_cast <float> (RAND_MAX)+1)*(static_cast <float> (RAND_MAX)+1) - 1);
//	assert(r >= 0 && r <= 1);
	return r;
}
vec3 UniformSampleSphere(float u1, float u2){
	float z = 1.0f - 2.0f*u1;
	float r = sqrtf(max(0.0f, 1.0f - z*z));
	float phi = 2.0f*M_PI*u2;
	float x = r*cosf(phi);
	float y = r*sinf(phi);
	return vec3(x, y, z);
}
#endif