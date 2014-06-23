#if !defined( SIMPLE_VERTEX_XYZ )
#define SIMPLE_VERTEX_XYZ

// ====================================
// SimpleVertexXYZ
// Simple VertexPosition Implementation for OSD Hbr/Far
// ====================================
template< typename T >
struct SimpleVertexXYZ{
public:
	//Typedefs
	typedef T ValType;

public:
	//Creators
	inline SimpleVertexXYZ(){ p_[0] = p_[1] = p_[2] = static_cast<ValType>( 0 ); };
	inline SimpleVertexXYZ( int ){ p_[0] = p_[1] = p_[2] = static_cast<ValType>( 0 ); };
	inline SimpleVertexXYZ( const ValType& x,
							const ValType& y,
							const ValType& z )
	{
		p_[0] = x;
		p_[1] = y;
		p_[2] = z;
	};

	inline SimpleVertexXYZ( const SimpleVertexXYZ& rhs ){
		p_[0] = rhs.p_[0];
		p_[1] = rhs.p_[1];
		p_[2] = rhs.p_[2];
	};

	template< typename O >
	inline SimpleVertexXYZ( const typename SimpleVertexXYZ<O>& rhs ){
		p_[0] = static_cast<T>( rhs.p_[0] );
		p_[1] = static_cast<T>( rhs.p_[1] );
		p_[2] = static_cast<T>( rhs.p_[2] );
	};

	virtual ~SimpleVertexXYZ(){};

public:
	//Members
	inline void AddWithWeight( const SimpleVertexXYZ& src, float weight ){
		p_[0] += src.p_[0] * weight;
		p_[1] += src.p_[1] * weight;
		p_[2] += src.p_[2] * weight;
	};

	inline void AddVaryingWithWeight( const SimpleVertexXYZ&, float ){ /*Do nothing*/ };
	
	inline void Clear( void* pArg = 0 ){ p_[0] = p_[1] = p_[2] = static_cast<ValType>(0); };
	
	inline void ApplyVertexEdit( const OpenSubdiv::HbrVertexEdit<SimpleVertexXYZ> & edit) {
		const float *src = edit.GetEdit();
		switch(edit.GetOperation()) {
		 case OpenSubdiv::HbrHierarchicalEdit<SimpleVertexXYZ>::Set:
			p_[0] = src[0];
			p_[1] = src[1];
			p_[2] = src[2];
			break;
		 case OpenSubdiv::HbrHierarchicalEdit<SimpleVertexXYZ>::Add:
			p_[0] += src[0];
			p_[1] += src[1];
			p_[2] += src[2];
			break;
		 case OpenSubdiv::HbrHierarchicalEdit<SimpleVertexXYZ>::Subtract:
			p_[0] -= src[0];
			p_[1] -= src[1];
			p_[2] -= src[2];
			break;
		}
	}

	void ApplyVertexEdit(OpenSubdiv::FarVertexEdit const & edit) {
		const float *src = edit.GetEdit();
		switch(edit.GetOperation()) {
		 case OpenSubdiv::FarVertexEdit::Set:
			p_[0] = src[0];
			p_[1] = src[1];
			p_[2] = src[2];
			break;
		 case OpenSubdiv::FarVertexEdit::Add:
			p_[0] += src[0];
			p_[1] += src[1];
			p_[2] += src[2];
			break;
		}
	}

	void ApplyMovingVertexEdit(const OpenSubdiv::HbrMovingVertexEdit<SimpleVertexXYZ> &) { /*Do nothing*/ }

	const ValType* GetPos() const { return p_; }

	//
	ValType& X( void ){ return p_[0]; };
	const ValType& X( void )const { return p_[0]; };
	ValType& Y( void ){ return p_[1]; };
	const ValType& Y( void )const { return p_[1]; };
	ValType& Z( void ){ return p_[2]; };
	const ValType& Z( void )const { return p_[2]; };

protected:
	ValType p_[3];
	
};//SimpleVertexXYZ

typedef SimpleVertexXYZ<float>	SimpleVertexXYZf;
typedef SimpleVertexXYZ<double>	SimpleVertexXYZd;
typedef SimpleVertexXYZ<int>	SimpleVertexXYZi;

#endif