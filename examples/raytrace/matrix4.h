#ifndef __MALLIE_MATRIX4_H__
#define __MALLIE_MATRIX4_H__

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
	
	class matrix4 {
	public:
        matrix4() { }
        matrix4(
        	float m00, float m01, float m02, float m03,
        	float m10, float m11, float m12, float m13,
        	float m20, float m21, float m22, float m23,
        	float m30, float m31, float m32, float m33
        ) {
            m[0] = m00;
            m[1] = m01;
            m[2] = m02;
            m[3] = m03;
            
            m[4] = m10;
            m[5] = m11;
            m[6] = m12;
            m[7] = m13;
            
            m[8] = m20;
            m[9] = m21;
            m[10] = m22;
            m[11] = m23;
            
            m[12] = m30;
            m[13] = m31;
            m[14] = m32;
            m[15] = m33;

        }
        matrix4(const matrix4& rhs){
        	memcpy(m,rhs.m,sizeof(float)*16);
        }
        
        float * operator[] (int i){            //pointer what self is const.
	        return element[i];
	    }
	    const float * operator[] (int i)const{//pointer what self is const.
	        return element[i];
	    }
        
        
        // 0, 1, 2, 3
        // 4, 5, 6, 7
        // 8, 9,10,11
        //12,13,14,15
        matrix4& transpose(){
            std::swap(m[1],m[4]);
            std::swap(m[2],m[8]);
            std::swap(m[3],m[12]);
            std::swap(m[6],m[9]);
            std::swap(m[7],m[13]);
            std::swap(m[11],m[14]);
            return *this;
        }

        union{
        	float m[16];
        	float element[4][4];
        };
    };
    
    
    //----------------------------------------------------------
    inline vector3 operator*(const matrix4& m, const vector3& v)
    {
        float r[4];
        r[0] = m.m[0] * v[0] + m.m[1] * v[1] + m.m[2] * v[2] + m.m[3];
        r[1] = m.m[4] * v[0] + m.m[5] * v[1] + m.m[6] * v[2] + m.m[7];
        r[2] = m.m[8] * v[0] + m.m[9] * v[1] + m.m[10] * v[2] + m.m[11];
        r[3] = m.m[12] * v[0] + m.m[13] * v[1] + m.m[14] * v[2] + m.m[15];
        
        float ir = 1.0f/r[3];
        return vector3(r[0]*ir,r[1]*ir,r[2]*ir);
    }
    
    inline matrix4 operator*(const matrix4& lhs, const matrix4& rhs)
    {
    	return matrix4(
            lhs[0][0] * rhs[0][0] + lhs[0][1] * rhs[1][0] + lhs[0][2] * rhs[2][0] + lhs[0][3] * rhs[3][0],
            lhs[0][0] * rhs[0][1] + lhs[0][1] * rhs[1][1] + lhs[0][2] * rhs[2][1] + lhs[0][3] * rhs[3][1],
            lhs[0][0] * rhs[0][2] + lhs[0][1] * rhs[1][2] + lhs[0][2] * rhs[2][2] + lhs[0][3] * rhs[3][2],
            lhs[0][0] * rhs[0][3] + lhs[0][1] * rhs[1][3] + lhs[0][2] * rhs[2][3] + lhs[0][3] * rhs[3][3],
        
            lhs[1][0] * rhs[0][0] + lhs[1][1] * rhs[1][0] + lhs[1][2] * rhs[2][0] + lhs[1][3] * rhs[3][0],
            lhs[1][0] * rhs[0][1] + lhs[1][1] * rhs[1][1] + lhs[1][2] * rhs[2][1] + lhs[1][3] * rhs[3][1],
            lhs[1][0] * rhs[0][2] + lhs[1][1] * rhs[1][2] + lhs[1][2] * rhs[2][2] + lhs[1][3] * rhs[3][2],
            lhs[1][0] * rhs[0][3] + lhs[1][1] * rhs[1][3] + lhs[1][2] * rhs[2][3] + lhs[1][3] * rhs[3][3],
        
            lhs[2][0] * rhs[0][0] + lhs[2][1] * rhs[1][0] + lhs[2][2] * rhs[2][0] + lhs[2][3] * rhs[3][0],
            lhs[2][0] * rhs[0][1] + lhs[2][1] * rhs[1][1] + lhs[2][2] * rhs[2][1] + lhs[2][3] * rhs[3][1],
            lhs[2][0] * rhs[0][2] + lhs[2][1] * rhs[1][2] + lhs[2][2] * rhs[2][2] + lhs[2][3] * rhs[3][2],
            lhs[2][0] * rhs[0][3] + lhs[2][1] * rhs[1][3] + lhs[2][2] * rhs[2][3] + lhs[2][3] * rhs[3][3],
        
            lhs[3][0] * rhs[0][0] + lhs[3][1] * rhs[1][0] + lhs[3][2] * rhs[2][0] + lhs[3][3] * rhs[3][0],
            lhs[3][0] * rhs[0][1] + lhs[3][1] * rhs[1][1] + lhs[3][2] * rhs[2][1] + lhs[3][3] * rhs[3][1],
            lhs[3][0] * rhs[0][2] + lhs[3][1] * rhs[1][2] + lhs[3][2] * rhs[2][2] + lhs[3][3] * rhs[3][2],
            lhs[3][0] * rhs[0][3] + lhs[3][1] * rhs[1][3] + lhs[3][2] * rhs[2][3] + lhs[3][3] * rhs[3][3]        
        );
    }
    
    inline matrix4 transpose(const matrix4& rhs){
    	return (matrix4(rhs).transpose());
    }
    
    class mat4_gen{
    public:
    	inline static matrix4 rotation_z(const float s, const float c){
	        return
	            matrix4(
	                c,-s, 0, 0,
	                s, c, 0, 0,
	                0, 0, 1, 0,
	                0, 0, 0, 1
	            );
	    }
	
	    inline static matrix4 rotation_z(const float radian){
	        float _sin = sin(radian);
	        float _cos = cos(radian);
	        return rotation_z(_sin,_cos);
	    }
	    
	    inline static matrix4 scaling(const float x, const float y, const float z){
	        return
	            matrix4(
	                x, 0, 0, 0,
	                0, y, 0, 0,
	                0, 0, z, 0,
	                0, 0, 0, 1
	            );
	    }

    };
    
}


#endif
