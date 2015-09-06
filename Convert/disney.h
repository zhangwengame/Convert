#if !defined(__DISNEY_H)
#include "share.h"
#define __DISNEY_H
#define PI 3.1415926535897932
//Diffuse_Lambert + D_GGX * G_Smith * F_None
/*float3 Diffuse_Lambert(float3 DiffuseColor, float NoL)
{
	return NoL * DiffuseColor / PI;
}
*/
// GGX / Trowbridge-Reitz
// [Walter et al. 2007, "Microfacet models for refraction through rough surfaces"]
float D_GGX(float Roughness, float NoH)
{
	float m = Roughness * Roughness;
	float m2 = m * m;
	float d = (NoH * m2 - NoH) * NoH + 1;	// 2 mad
	return m2 / (d*d)/PI;					// 2 mul, 1 rcp
}

// Smith term for GGX modified by Disney to be less "hot" for small roughness values
// [Smith 1967, "Geometrical shadowing of a random rough surface"]
// [Burley 2012, "Physically-Based Shading at Disney"]
float G_Smith(float Roughness, float NoV, float NoL) //V ≥ˆ…‰ L »Î…‰ 
{
	float a = sqrt(Roughness);
	float a2 = a*a;

	float G_SmithV = NoV + sqrt(NoV * (NoV - NoV * a2) + a2);
	float G_SmithL = NoL + sqrt(NoL * (NoL - NoL * a2) + a2);
	return rcp(G_SmithV * G_SmithL);
}

/*float3 F_None(float3 SpecularColor)
{
	return SpecularColor;
}*/
#endif