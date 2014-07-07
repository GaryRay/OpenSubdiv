//Maya
#include <maya/MStringArray.h>
#include<maya/MFnMesh.h>
#include<maya/MItMeshPolygon.h>

//OSD
#include <osd/vertex.h>
#include <osd/mesh.h>
#include <osd/cpuComputeContext.h>
#include <osd/cpuComputeController.h>
#include <osd/cpuVertexBuffer.h>
#include <far/meshFactory.h>
#include <far/dispatcher.h>

//Local
#include "simpleVertexXYZ.h"
#include "mayaFvarDataDesc.h"

#if !defined( MAYA_MESH_TOHBR_MESH_INCLUDED__ )
#define MAYA_MESH_TOHBR_MESH_INCLUDED__

// ====================================
// MayaMeshToHbrMesh
// @Memo:
// ====================================
template< typename T >
class MayaMeshToHbrMesh{
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
	virtual HMesh* operator()(	MFnMesh const & inMeshFn,
								MItMeshPolygon & inMeshItPolygon,
								std::vector<int> & fvarIndices,
								std::vector<int> & fvarWidths,
								MayaFVarDataDesc* fvarDataDesc,
								float * maxCreaseSharpness = 0 );

	virtual MStatus getMayaFvarFieldParams(	MFnMesh const & inMeshFn,
											MStringArray & uvSetNames,
											MStringArray & colorSetNames,
											std::vector<int> & colorSetChannels,
											std::vector<MFnMesh::MColorRepresentation> &colorSetReps,
											int & totalColorSetChannels);

 protected:
	virtual void setHbrVertices( HMesh* hbrMesh,
								 MFnMesh const & inMeshFn );


	virtual float applyCreaseEdges(	MFnMesh const & inMeshFn,
									HMesh * hbrMesh);

	virtual float applyCreaseVertices(	MFnMesh const & inMeshFn,
										HMesh * hbrMesh );

};
#include "mayaMeshToHbrMesh.inl"

#endif
