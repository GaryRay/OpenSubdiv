#ifndef __BEZIER_PATCH_HPP__
#define __BEZIER_PATCH_HPP__

#include "vector2.h"
#include "vector3.h"
#include "matrix3.h"
#include "matrix4.h"
#include "bezier.h"
#include <vector>

namespace mallie{
	
	template<class T>
	class bezier_patch
	{
	public:
		typedef T value_type;
		typedef bezier_patch<T> this_type;

		static const int DEFAULT_ORDER = 4;
	protected:
		std::vector<T> default_cp(int nu, int nv){
			return std::vector<T>(nu*nv);
		}
	public:
		bezier_patch()
			:nu_(DEFAULT_ORDER), nv_(DEFAULT_ORDER), cp_(default_cp(DEFAULT_ORDER,DEFAULT_ORDER))
		{
		}

		bezier_patch(int nu, int nv)
			:nu_(nu), nv_(nv), cp_(default_cp(nu,nv))
		{
		}

		bezier_patch(int nu, int nv, const std::vector<T>& p)
			:nu_(nu),nv_(nv),cp_(p)
		{
			assert(nu*nv == p.size());
		}

                template<class S>
                bezier_patch(int nu, int nv, S const *p)
                    :nu_(nu),nv_(nv),cp_(default_cp(nu,nv))
		{
                    for (int i = 0; i < nu*nv; ++i) cp_[i] = p[i];
			assert(nu*nv == p.size());
		}

		bezier_patch(const bezier_patch<T>& rhs)
			:nu_(rhs.nu_),nv_(rhs.nv_),cp_(rhs.cp_)
		{
		}

		void set_cp(int nu, int nv, const std::vector<T>& cp)
		{
			assert(nu*nv == cp.size());
			nu_ = nu;
			nv_ = nv;
			cp_ = cp;
		}

		this_type& operator=(const bezier_patch<T>& rhs)
		{
			set_cp(rhs.nu_,rhs.nv_,rhs.cp_);
			return *this;
		}

		void swap(bezier_patch<T>& rhs)
		{
			std::swap(nu_,rhs.nu_);
			std::swap(nv_,rhs.nv_);
			cp_.swap(rhs.cp_);
		}
	public:
		int get_nu()const{return nu_;}
		int get_nv()const{return nv_;}
		T get_at(int i, int j)const{return cp_[nu_*j+i];}
		void    set_at(int i, int j, const T& p){cp_[nu_*j+i] = p;}

		int get_cp_size_u()const{return get_nu();}
		int get_cp_size_v()const{return get_nv();}
		T get_cp_at(int i, int j)const{return get_at(i,j);}
		void    set_cp_at(int i, int j, const T& p){set_at(i,j,p);}

		const std::vector<T>& get_cp()const{return cp_;}
	public:
		T evaluate     (real u, real v)const
		{
			int nu = get_nu();
			int nv = get_nv();
			return bezier_evaluate( &(cp_[0]), nu, nv, u, v);
		}

        T evaluate_deriv_u(real u, real v)const
		{
            int nu = get_nu();
			int nv = get_nv();
			return bezier_evaluate_deriv_u( &(cp_[0]), nu, nv, u, v);
        }

		T evaluate_dPdU(real u, real v)const
		{
			return evaluate_deriv_u(u, v);
		}

		T evaluate_deriv_v(real u, real v)const
		{
			int nu = get_nu();
			int nv = get_nv();
			return bezier_evaluate_deriv_v( &(cp_[0]), nu, nv, u, v);
		}

        T evaluate_dPdV(real u, real v)const
		{
            return evaluate_deriv_v(u, v);
        }

		T min()const
		{
			return bezier_min( &(cp_[0]), nu_*nv_ );
		}

		T max()const
		{
			return bezier_max( &(cp_[0]), nu_*nv_ );
		}

		void minmax(T& min, T& max)const
		{
			bezier_minmax(min, max,  &(cp_[0]), nu_*nv_ );
		}
	public:
		this_type& swap_u()
		{
			int nu = get_nu();
			int nv = get_nv();

			std::vector<T> cp(nu*nv);
			for(int u = 0;u<nu;u++){
				for(int v = 0;v<nv;v++){
					cp[nu*v+(nu-1-u)] = cp_[nu*v+u];
				}
			}		
			cp_ = cp;

			return *this;
		}

		this_type& swap_v()
		{
			int nu = get_nu();
			int nv = get_nv();

			std::vector<T> cp(nu*nv);
			for(int u = 0;u<nu;u++){
				for(int v = 0;v<nv;v++){
					cp[nu*(nv-1-v)+u] = cp_[nu*v+u];
				}
			}		
			cp_ = cp;

			return *this;
		}

	public:
		bool equal(const bezier_patch<T>& rhs)const
		{
			return (get_nu() == rhs.get_nu())&&(get_nv() == rhs.get_nv())&&(get_cp() == rhs.get_cp());
		}

		this_type& transform(const matrix3& m)
		{
			size_t sz = cp_.size();
			for(size_t i = 0;i<sz;i++){
				cp_[i] = m*cp_[i];
			}
			return *this;
		}

		this_type& transform(const matrix4& m)
		{
			size_t sz = cp_.size();
			for(size_t i = 0;i<sz;i++){
				cp_[i] = m*cp_[i];
			}
			return *this;
		}
	public:
		void split(this_type patches[4], real u, real v)const
		{
			int nu = nu_;
			int nv = nv_;
			int sz = nu*nv;
			for(int i=0;i<4;i++){
				patches[i].nu_ = nu;
				patches[i].nv_ = nv;	
				patches[i].cp_.resize(sz);
			}
			bezier_split(
				&(patches[0].cp_[0]),
				&(patches[1].cp_[0]),
				&(patches[2].cp_[0]),
				&(patches[3].cp_[0]),
				&(cp_[0]),
				nu,nv,
				u,v
			);
		}
        void split_u(this_type patches[2], real u)const
        {
            int nu = nu_;
			int nv = nv_;
			int sz = nu*nv;
			for(int i=0;i<2;i++){
				patches[i].nu_ = nu;
				patches[i].nv_ = nv;	
				patches[i].cp_.resize(sz);
			}
			bezier_split_u(
                &(patches[0].cp_[0]),
                &(patches[1].cp_[0]),
                &(cp_[0]),
                nu,nv,
                u
            );
        }
        void split_v(this_type patches[2], real v)const
        {
            int nu = nu_;
			int nv = nv_;
			int sz = nu*nv;
			for(int i=0;i<2;i++){
				patches[i].nu_ = nu;
				patches[i].nv_ = nv;	
				patches[i].cp_.resize(sz);
			}
			bezier_split_v(
                &(patches[0].cp_[0]),
                &(patches[1].cp_[0]),
                &(cp_[0]),
                nu,nv,
                v
            );
        }

		void crop_u(this_type& patch, real u0, real u1)const
        {
			int nu = nu_;
			int nv = nv_;
			int sz = nu*nv;
			{
				patch.nu_ = nu;
				patch.nv_ = nv;	
				patch.cp_.resize(sz);
			}
			bezier_crop_u(
				&(patch.cp_[0]),
                &(cp_[0]),
                nu,nv,	
				u0,u1
			);
		}

		void crop_v(this_type& patch, real v0, real v1)const
        {
			int nu = nu_;
			int nv = nv_;
			int sz = nu*nv;
			{
				patch.nu_ = nu;
				patch.nv_ = nv;	
				patch.cp_.resize(sz);
			}
			bezier_crop_v(
				&(patch.cp_[0]),
                &(cp_[0]),
                nu,nv,	
				v0,v1
			);

		}
	private:
		int nu_;
		int nv_;
		std::vector<T> cp_;
	};

	template<class T>
	inline bool operator==(const bezier_patch<T>& lhs, const bezier_patch<T>& rhs)
	{
		return lhs.equal(rhs);
	}

	template<class T>
	inline bool operator!=(const bezier_patch<T>& lhs, const bezier_patch<T>& rhs)
	{
		return !lhs.equal(rhs);
	}

	template<class T>
	inline bezier_patch<T> operator*(const matrix3& m, const bezier_patch<T>& patch)
	{
		return bezier_patch<T>(patch).transform(m);
	}

	template<class T>
	inline bezier_patch<T> operator*(const matrix4& m, const bezier_patch<T>& patch)
	{
		return bezier_patch<T>(patch).transform(m);
	}

	typedef bezier_patch<float>    bezier_patch1;
	typedef bezier_patch<vector2> bezier_patch2;
	typedef bezier_patch<vector3> bezier_patch3;
	//typedef bezier_patch<vector4> bezier_patch4;


}

#endif
