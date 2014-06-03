#include "bvh_composite_intersection.h"

#include <functional>
#include <vector>
#include <limits>

//#include "logger.h"
#include "count_ptr.hpp"

#include <memory>

#define print_log printf

#define USE_SSE 1


#define test_AABB IntersectAABB

namespace mallie{

	static
    inline vector3 Conv(const real3& r)
    {
        return vector3(r.x, r.y, r.z); 
    }

    static
    inline real3 Conv(const vector3& r)
    {
        return real3(r[0], r[1], r[2]); 
    }
    
	typedef const intersector* PCINTER;
	typedef PCINTER* INTER_ITERATOR;
	
    struct range_AABB{
        float tmin;
        float tmax;
    };

    static inline bool IntersectAABB(range_AABB* rng, const vector3 & min, const vector3 & max, const Ray & r, float tmin, float tmax)
    {
        //int phase = r.phase();
        int sign[3] = {r.dirSign[0],r.dirSign[1],r.dirSign[2]};
        vector3 box[2] = {min,max};
        vector3 org  = Conv(r.org);
        vector3 idir = Conv(r.invDir);

        for(int i = 0;i<3;i++){
            tmin = std::max<float>(tmin,(box[  sign[i]][i]-org[i])*idir[i]);
            tmax = std::min<float>(tmax,(box[1-sign[i]][i]-org[i])*idir[i]);
        }

        rng->tmin = tmin;
        rng->tmax = tmax;

		return rng->tmin<=rng->tmax;
    }


	//-------------------------------------------------------------------------------
	namespace {
		
		class bvh_node{
		public:
			virtual ~bvh_node(){}
		public:
			virtual bool test(const Ray& r, float tmin, float tmax)const = 0;
			virtual bool test(Intersection* info, const Ray& r, float tmin, float tmax)const = 0;
			virtual vector3 min()const=0;
			virtual vector3 max()const=0;
		};

		class bvh_node_null:public bvh_node{
		public:
			bool test(const Ray& r, float tmin, float tmax)const{return false;}
			bool test(Intersection* info, const Ray& r, float tmin, float tmax)const{return false;}
			vector3 min()const{return vector3(0,0,0);}
			vector3 max()const{return vector3(0,0,0);}
		};

		class bvh_node_branch:public bvh_node{
		public:
			bvh_node_branch(INTER_ITERATOR s, INTER_ITERATOR e);
			~bvh_node_branch();
		protected:
			void initialize(INTER_ITERATOR s, INTER_ITERATOR e);
		public:
			bool test(const Ray& r, float tmin, float tmax)const;
			bool test(Intersection* info, const Ray& r, float tmin, float tmax)const;
			vector3 min()const;
			vector3 max()const;
		protected:
			bool test_internal(const Ray & r, float tmin, float tmax)const;
			bool test_internal(Intersection* info, const Ray & r, float tmin, float tmax)const;
		private:
			int plane_;
			bvh_node* nodes_[2];
			vector3 min_;
			vector3 max_;
		};

		class bvh_node_leaf:public bvh_node{
		public:
			bvh_node_leaf(INTER_ITERATOR s, INTER_ITERATOR e);
			~bvh_node_leaf();
		public:
			bool test(const Ray& r, float tmin, float tmax)const;
			bool test(Intersection* info, const Ray& r, float tmin, float tmax)const;
			vector3 min()const;
			vector3 max()const;
		private:
			PCINTER p_inter_;
		};

		//------------------------------------------------------------
		//------------------------------------------------------------
		inline int eval_phase(int phase,int axis){//rid<0
			return (phase>>axis) & 0x1;
		}

		bvh_node_branch::bvh_node_branch(INTER_ITERATOR s, INTER_ITERATOR e)
		{
			nodes_[0] = 0;
			nodes_[1] = 0;
			initialize(s,e);
		}

		bvh_node_branch::~bvh_node_branch()
		{
			if(nodes_[0])delete nodes_[0];
			if(nodes_[1])delete nodes_[1];
		}

		struct bvh_sorter:public std::binary_function<PCINTER, PCINTER, bool>{
			bvh_sorter(int plane):plane_(plane){}
			bool operator()(PCINTER a, PCINTER b)const{
				int plane = plane_;
				float ac = (a->min()[plane]+a->max()[plane]);
				float bc = (b->min()[plane]+b->max()[plane]);
				return ac<bc;
			}
			int plane_;
		};

		void bvh_node_branch::initialize(INTER_ITERATOR s, INTER_ITERATOR e){
            static const float EPSILON = std::numeric_limits<float>::epsilon()*1024;
            
			vector3 min,max;
			min = (*s)->min();
			max = (*s)->max();
			{
				INTER_ITERATOR it = s+1;
				while(it != e){
					vector3 cmin = (*it)->min();
					vector3 cmax = (*it)->max();
					for(int i=0;i<3;i++){
						if(cmin[i]<min[i])min[i]=cmin[i];
						if(cmax[i]>max[i])max[i]=cmax[i];
					}
					it++;
				}
			}
            
            min -= vector3(EPSILON,EPSILON,EPSILON);
            max += vector3(EPSILON,EPSILON,EPSILON);

			vector3 wid = max-min;

			int plane = 0;
			if(wid[plane]<wid[1])plane=1;
			if(wid[plane]<wid[2])plane=2;

			std::sort(s,e, bvh_sorter(plane));

			size_t sz = e-s;
			size_t msz = sz>>1;

			assert(sz>0);

			INTER_ITERATOR m = s+msz;

			std::auto_ptr<bvh_node> node_l;
			std::auto_ptr<bvh_node> node_r;

			size_t lsz = m-s;
			size_t rsz = e-m;

			if(lsz){
				if(lsz == 1){
					node_l.reset(new bvh_node_leaf(s,m));
				}else{
					node_l.reset(new bvh_node_branch(s,m));
				}
			}

			if(rsz){
				if(rsz == 1){
					node_r.reset(new bvh_node_leaf(m,e));
				}else{
					node_r.reset(new bvh_node_branch(m,e));
				}
			}


			nodes_[0] = node_l.release();
			nodes_[1] = node_r.release();
            
			plane_ = plane;
			min_ = min;
			max_ = max;
		}

		bool bvh_node_branch::test(const Ray& r, float tmin, float tmax)const{
			range_AABB rng;
			if(test_AABB(&rng, min_, max_, r, tmin, tmax)){
				tmin = std::max<float>(tmin, rng.tmin);
				tmax = std::min<float>(tmax, rng.tmax);
				return this->test_internal(r, tmin, tmax);
			}
			return false;
		}

		bool bvh_node_branch::test(Intersection* info, const Ray& r, float tmin, float tmax)const{
			range_AABB rng;
			if(test_AABB(&rng, min_, max_, r, tmin, tmax)){
				tmin = std::max<float>(tmin, rng.tmin);
				tmax = std::min<float>(tmax, rng.tmax);
				return this->test_internal(info, r, tmin, tmax);
			}
			return false; 
		}

		bool bvh_node_branch::test_internal(const Ray & r, float tmin, float tmax)const{
			int nNear = r.dirSign[plane_];
			int nFar  = 1-nNear;

			if(nodes_[nNear] && nodes_[nNear]->test(r, tmin, tmax))return true;
			if(nodes_[nFar ] && nodes_[nFar ]->test(r, tmin, tmax))return true;

			return false;
		}

		bool bvh_node_branch::test_internal(Intersection* info, const Ray & r, float tmin, float tmax)const{
			int nNear = r.dirSign[plane_];
			int nFar  = 1-nNear;

			bool bRet = false; 
			if(nodes_[nNear] && nodes_[nNear]->test(info, r, tmin, tmax)){
				tmax = info->t;
				bRet = true;
			}

			if(nodes_[nFar ] && nodes_[nFar ]->test(info, r, tmin, tmax)){
				tmax = info->t;
				bRet = true;
			}

			return bRet;
		}

		vector3 bvh_node_branch::min()const{
			return min_;
		}

		vector3 bvh_node_branch::max()const{
			return max_;
		}

		//------------------------------------------------------------
		//------------------------------------------------------------


		bvh_node_leaf::bvh_node_leaf(INTER_ITERATOR s, INTER_ITERATOR e){
			assert((e-s)==1);
			p_inter_ = *s;
		}

		bvh_node_leaf::~bvh_node_leaf(){
			//
		}

		bool bvh_node_leaf::test(const Ray& r, float tmin, float tmax)const{
			return p_inter_->test(r, tmax);
		}

		bool bvh_node_leaf::test(Intersection* info, const Ray& r, float tmin, float tmax)const{
			return p_inter_->test(info, r, tmax);
		}
		
		vector3 bvh_node_leaf::min()const{
			return p_inter_->min();
		}

		vector3 bvh_node_leaf::max()const{
			return p_inter_->max();
		}

		//------------------------------------------------------------
		//------------------------------------------------------------



	}
	//-------------------------------------------------------------------------------
	class bvh_composite_intersection_imp{
	public:
		bvh_composite_intersection_imp();
		~bvh_composite_intersection_imp();
	public:
		bool test(const Ray& r, float dist)const;
		bool test(Intersection* info, const Ray& r, float dist)const;
		//void finalize(Intersection* info, const Ray& r, float dist)const;// -- //NO IMPLEMENT!
		
		vector3 min()const;
		vector3 max()const;
	public:
		void construct();
		void add(const count_ptr<intersector>& inter);
	private:
		bvh_node* root_;
		std::vector< count_ptr<intersector> > mv_;
	};

	//-------------------------------------------------------------------------------

	bvh_composite_intersection_imp::bvh_composite_intersection_imp(){
		root_ = new bvh_node_null();
	}

	bvh_composite_intersection_imp::~bvh_composite_intersection_imp(){
		delete root_;
	}

	bool bvh_composite_intersection_imp::test(const Ray& r, float dist)const{
		return root_->test(r, 0, dist);
	}

	bool bvh_composite_intersection_imp::test(Intersection* info, const Ray& r, float dist)const{
		return root_->test(info, r, 0, dist);
	}
		
	vector3 bvh_composite_intersection_imp::min()const{
		return root_->min();
	}

	vector3 bvh_composite_intersection_imp::max()const{
		return root_->max();
	}
    
	void bvh_composite_intersection_imp::add(const count_ptr<intersector>& inter){
		mv_.push_back(inter);
	}

	void bvh_composite_intersection_imp::construct(){
		size_t sz = mv_.size();
        if(sz==0)return;
		std::vector<const intersector*> tmp;
		for(size_t i = 0;i<sz;i++){
			tmp.push_back(mv_[i].get());
		}
		std::auto_ptr<bvh_node_branch> ap( new bvh_node_branch(&tmp[0],&tmp[0]+sz) );

		if(root_)delete root_;
		root_ = ap.release();
	}

	//-------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------
	//-------------------------------------------------------------------------------

	
	bvh_composite_intersection::bvh_composite_intersection(){
		 imp_ = new bvh_composite_intersection_imp();
	}

	bvh_composite_intersection::~bvh_composite_intersection(){
		delete imp_;
	}
	
	bool bvh_composite_intersection::test(const Ray& r, float dist)const{
		return imp_->test(r, dist);
	}

	bool bvh_composite_intersection::test(Intersection* info, const Ray& r, float dist)const{
		return imp_->test(info, r, dist);
	}

	void bvh_composite_intersection::finalize(Intersection* info, const Ray& r, float dist)const{
		// -- //NO IMPLEMENT!
	}
		
	vector3 bvh_composite_intersection::min()const{
		return imp_->min();
	}

	vector3 bvh_composite_intersection::max()const{
		return imp_->max();
	}

	void bvh_composite_intersection::construct(){
		imp_->construct();
	}

	void bvh_composite_intersection::add(const auto_count_ptr<intersector>& inter){
		imp_->add(inter);
	}
}
