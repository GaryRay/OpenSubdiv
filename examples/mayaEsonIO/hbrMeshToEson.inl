#include <memory>

// ====================================
// HbrMeshToEson
// @Memo:
// ====================================

// - - - - - - - - - - - - - - - - - -
template< typename T >
void	HbrMeshToEson<T>::operator()( HMesh*	hbrMesh, std::string outputPath, int smoothLevel, bool useAdaptive ){

	using namespace OpenSubdiv;

	FMeshFactory meshFactory( hbrMesh, smoothLevel, useAdaptive );
	std::shared_ptr<FMesh> farMesh = std::shared_ptr<FMesh>( meshFactory.Create() );
	static OpenSubdiv::FarComputeController computeController;
	computeController.Refine( farMesh.get() );

	std::vector<float> vertices;
	{
		for (int i = 0; i < (int)farMesh->GetVertices().size(); ++i) {
			auto v = farMesh->GetVertices()[i];
			vertices.push_back( getVertexPosAsFloat( v, 0 ) );
			vertices.push_back( getVertexPosAsFloat( v, 1 ) );
			vertices.push_back( getVertexPosAsFloat( v, 2 ) );
		}
	}

	FarPatchTables const *patchTables = farMesh->GetPatchTables();
	FarPatchTables::PatchArrayVector const &patchArrays = patchTables->GetPatchArrayVector();
	FarPatchTables::PatchParamTable const &patchParam = patchTables->GetPatchParamTable();

	int numPatches = 0;
	std::vector<float> bezierVertices;
	std::vector<float> bezierBounds;
	std::vector<int> cpIndices; // 16 * numPatches
	// iterate patch types.
	for (FarPatchTables::PatchArrayVector::const_iterator it = patchArrays.begin();
		 it != patchArrays.end(); ++it) {

		switch(it->GetDescriptor().GetType()) {
		 case FarPatchTables::REGULAR:
			numPatches += convertRegular(bezierVertices, bezierBounds,
										 cpIndices, &vertices[0], patchTables, *it);
			break;
		 case FarPatchTables::BOUNDARY:
			numPatches += convertBoundary(bezierVertices, bezierBounds,
										  cpIndices, &vertices[0], patchTables, *it);
			break;
		 case FarPatchTables::CORNER:
			numPatches += convertCorner(bezierVertices, bezierBounds,
										cpIndices, &vertices[0], patchTables, *it);
			break;
		 case FarPatchTables::GREGORY:
			numPatches += convertGregory(bezierVertices, bezierBounds,
										 cpIndices, &vertices[0], patchTables, *it);
			break;
		 case FarPatchTables::GREGORY_BOUNDARY:
			numPatches += convertBoundaryGregory(bezierVertices, bezierBounds,
												 cpIndices, &vertices[0], patchTables, *it);
			break;
		 default:
			break;
		}

	}

	//FarPatch to eson
	eson::Object mesh;
	int nverts = farMesh->GetNumVertices();
	mesh["num_vertices"] =
		eson::Value((int64_t)nverts);
	mesh["vertices"] =
		eson::Value((uint8_t*)&vertices[0], sizeof(float)*nverts*3);
	if( numPatches > 0 ){
		mesh["num_bezier_patches"] = eson::Value((int64_t)numPatches);
	}
	if( patchParam.size() > 0 ){
		mesh["patch_param"] = eson::Value((uint8_t*)&patchParam[0], sizeof(OpenSubdiv::FarPatchParam)*patchParam.size());
	}

	if( bezierVertices.size() > 0 ){
		mesh["bezier_vertices"]  =eson::Value((uint8_t*)&bezierVertices[0], sizeof(float)*bezierVertices.size());
	}

	if( bezierBounds.size() > 0 ){
		mesh["bezier_bounds"] = eson::Value((uint8_t*)&bezierBounds[0], sizeof(float)*bezierBounds.size());
	}

	assert(numPatches*16*3 == (int)bezierVertices.size());

	eson::Value v = eson::Value(mesh);
	int64_t size = v.Size();

	std::vector<uint8_t> buf(size);
	uint8_t* ptr = &buf[0];

	ptr = v.Serialize(ptr);
	assert((ptr-&buf[0]) == size);


	{
		FILE* fp = fopen( outputPath.c_str(), "wb" );
		if( !fp ){
			throw std::exception( "fopen failed" );
		}
		
		fwrite( &buf[0], 1, size, fp );
		fclose( fp );
	}
	
}

// - - - - - - - - - - - - - - - - - -
template< typename T >
float HbrMeshToEson<T>::getVertexPosAsFloat( const VertexType& v, int index ){
	//@Todo: fix
	return v[index];
}

//Specialized
template<>
float HbrMeshToEson<SimpleVertexXYZf>::getVertexPosAsFloat( const VertexType& v, int index ){
	assert( 0 <= index && index < 3 );
	return v.GetPos()[index];
}
