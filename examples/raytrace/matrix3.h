#ifndef __MALLIE_MATRIX3_H__
#define __MALLIE_MATRIX3_H__

#include <cstddef>//std::size_t
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <sstream>
#include <cassert>
#include <cstring>

#include "vector3.h"

namespace mallie {
	
	template<class T>
	class matrix3t {
	public:
        matrix3t() { }
        matrix3t(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22) {
            m[0] = m00;
            m[1] = m01;
            m[2] = m02;
            m[3] = m10;
            m[4] = m11;
            m[5] = m12;
            m[6] = m20;
            m[7] = m21;
            m[8] = m22;
        }
        matrix3t(const matrix3t& rhs){
        	memcpy(m,rhs.m,sizeof(T)*9);
        }
        
              T* operator[](size_t i)     {return m+3*i;}
        const T* operator[](size_t i)const{return m+3*i;}
        
        
        //0,1,2
        //3,4,5
        //6.7.8
        void transpose(){
            std::swap(m[1],m[3]);
            std::swap(m[2],m[6]);
            std::swap(m[5],m[7]);
        }

        T m[3*3];
    };
    
    
    //----------------------------------------------------------
    template<class T, class X>
    inline vector3t<X> operator*(const matrix3t<T>& m, const vector3t<X>& v)
    {
        vector3t<X> r;
        r[0] = m.m[0] * v[0] + m.m[1] * v[1] + m.m[2] * v[2];
        r[1] = m.m[3] * v[0] + m.m[4] * v[1] + m.m[5] * v[2];
        r[2] = m.m[6] * v[0] + m.m[7] * v[1] + m.m[8] * v[2];
        return r;
    }
    
    typedef matrix3t<float> matrix3;
    typedef matrix3t<float> matrix3f;
    typedef matrix3t<double> matrix3d;
    

	
	
	
}


#endif
