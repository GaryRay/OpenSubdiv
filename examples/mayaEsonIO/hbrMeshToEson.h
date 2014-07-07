//C/C++
#include <string>
#include <stdio.h>

//OSD
#include <osd/vertex.h>
#include <osd/mesh.h>
#include <osd/cpuComputeContext.h>
#include <osd/cpuComputeController.h>
#include <osd/cpuVertexBuffer.h>
#include <far/meshFactory.h>
#include <far/dispatcher.h>

//Local
#include "eson.h"
#include "simpleVertexXYZ.h"
#include "mayaFvarDataDesc.h"

#if !defined( HBR_MESH_TOESON_INCLUDED__ )
#define HBR_MESH_TOESON_INCLUDED__

// ====================================
// HbrMeshToEson
// @Memo:
// ====================================
template< typename T >
class HbrMeshToEson{
public:
	typedef T VertexType;

	typedef OpenSubdiv::HbrMesh<VertexType>               HMesh;
	typedef OpenSubdiv::HbrFace<VertexType>               HFace;
	typedef OpenSubdiv::HbrVertex<VertexType>             HVertex;
	typedef OpenSubdiv::HbrHalfedge<VertexType>           HHalfedge;
	typedef OpenSubdiv::HbrFVarData<VertexType>           HFvarData;
	typedef OpenSubdiv::HbrCatmarkSubdivision<VertexType> HCatmark;

	typedef OpenSubdiv::FarMesh<VertexType>               FMesh;
	typedef OpenSubdiv::FarMeshFactory<VertexType>        FMeshFactory;
	typedef OpenSubdiv::FarSubdivisionTables              FSubdivision;
	typedef OpenSubdiv::FarPatchTables                    FPatches;
	
public:
	void	operator()( HMesh*	hMesh,
						MayaFVarDataDesc& fvarDataDesc,
						std::vector<short>& faceMatIds,
						std::string outputPath,
						int smoothLevel,
						bool useAdaptive = true );

protected:
	virtual float	getVertexPosAsFloat( const VertexType& v, int index );

};

#include "hbrMeshToEson.inl"

#endif
