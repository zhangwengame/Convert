#if !defined(__MICROFACET_H)
#define __MICROFACET_H
#include "share.h"

class MicrofacetDistribution {
public:
	/// Supported distribution types
	enum EType {
		/// Beckmann distribution derived from Gaussian random surfaces
		EBeckmann         = 0,

		/// GGX: Long-tailed distribution for very rough surfaces (aka. Trowbridge-Reitz distr.)
		EGGX              = 1,

		/// Phong distribution (with the anisotropic extension by Ashikhmin and Shirley)
		EPhong            = 2
	};
	inline MicrofacetDistribution(EType type, float glossiness)
		: m_type(type), m_exponentU(glossiness), m_exponentV(glossiness) {
		m_exponentU = max(m_exponentU, 0.0f);
		m_exponentV = max(m_exponentV, 0.0f);
	}
	
	/// Return the distribution type
	inline EType getType() const { return m_type; }	

	/// Return the Phong exponent (isotropic case)
	inline float getExponent() const { return m_exponentU; }	

	inline float eval(const vec3 &m) const {
		if (m.z<= 0)
			return 0.0f;
		float result;
		switch (m_type) {	

			case EPhong: {
					/* Isotropic case: Phong distribution. Anisotropic case: Ashikhmin-Shirley distribution */
					float exponent = m_exponentU;
					result = sqrt((m_exponentU + 2) * (m_exponentV + 2))
						* INV_TWOPI * pow(m.z, exponent);
				}
				break;
		}		
		if (result * m.z < 1e-20f)
			result = 0;
		return result;
	}	


protected:
	EType m_type;
	float m_exponentU, m_exponentV;
};


#endif /* __MICROFACET_H */
