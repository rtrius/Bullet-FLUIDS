
#ifndef VECTOR3DF_DEF
#define VECTOR3DF_DEF

#include <math.h>

///USE_BTVECTOR3_IN_FLUIDS also enables btAlignedObjectArray for Fluid storage,
///(FluidSystem.m_particles) as std::vector does not seem to be able to store 
///btVector3 due to alignment issues. A side effect of this is that the OpenCL
///port breaks, possibly as a result of padding. Additionally, btAlignedObjectArray
///seems to be slightly slower than std::vector.

//#define USE_BTVECTOR3_IN_FLUIDS
#ifndef USE_BTVECTOR3_IN_FLUIDS
	typedef float VTYPE;
	struct Vector3DF 
	{
		VTYPE m_floats[4];

		inline Vector3DF() {}
		inline Vector3DF(const VTYPE x, const VTYPE y, const VTYPE z) {m_floats[0] = x; m_floats[1] = y; m_floats[2] = z;}

		inline Vector3DF &operator+= (const Vector3DF &op) 
			{m_floats[0] += op.m_floats[0]; m_floats[1] += op.m_floats[1]; m_floats[2] += op.m_floats[2]; return *this;}
		inline Vector3DF &operator-= (const Vector3DF &op) 
			{m_floats[0] -= op.m_floats[0]; m_floats[1] -= op.m_floats[1]; m_floats[2] -= op.m_floats[2]; return *this;}
			
		inline Vector3DF &operator*= (const VTYPE op) {m_floats[0] *= op; m_floats[1] *= op; m_floats[2] *= op; return *this;}
		inline Vector3DF &operator/= (const VTYPE op) {m_floats[0] /= op; m_floats[1] /= op; m_floats[2] /= op; return *this;}

		inline const VTYPE &x() const	{ return m_floats[0]; }
		inline const VTYPE &y() const	{ return m_floats[1]; }
		inline const VTYPE &z() const	{ return m_floats[2]; }
				
		inline VTYPE dot(const Vector3DF &vec) const 
			{ return m_floats[0]*vec.m_floats[0] + m_floats[1]*vec.m_floats[1] + m_floats[2]*vec.m_floats[2]; }
		inline VTYPE length() const { return sqrt( length2() ); }
		inline VTYPE length2() const { return m_floats[0]*m_floats[0] + m_floats[1]*m_floats[1] + m_floats[2]*m_floats[2]; }
		
		inline void setValue (const VTYPE x, const VTYPE y, const VTYPE z) {m_floats[0] = x; m_floats[1] = y; m_floats[2] = z;}
		inline Vector3DF &normalize()
		{
			VTYPE n = m_floats[0]*m_floats[0] + m_floats[1]*m_floats[1] + m_floats[2]*m_floats[2];
			if (n != 0.0)
			{
				n = sqrt(n);
				m_floats[0] /= n; 
				m_floats[1] /= n; 
				m_floats[2] /= n;
			}
			return *this;
		}
	};
#else
	#include "LinearMath/btVector3.h"
	typedef btVector3 Vector3DF;
#endif

#endif

