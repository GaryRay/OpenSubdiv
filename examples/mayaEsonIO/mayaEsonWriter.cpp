//C/C++
#include<fstream>
#include<time.h>

//Maya
#include <maya/MStringArray.h>

#include<maya/MItDag.h>
#include<maya/MItDependencyGraph.h>
#include<maya/MItDependencyNodes.h>
#include<maya/MItGeometry.h>
#include<maya/MItSelectionList.h>
#include<maya/MPlug.h>
#include<maya/MPlugArray.h>
#include<maya/MGlobal.h>

#include<maya/MUintArray.h>
#include<maya/MFloatArray.h>
#include<maya/MDoubleArray.h>
#include<maya/MPointArray.h>
#include<maya/MFloatPointArray.h>
#include<maya/MObjectArray.h>

#include<maya/MVector.h>
#include<maya/MMatrix.h>
#include<maya/MEulerRotation.h>
#include<maya/MTransformationMatrix.h>

#include<maya/MSelectionList.h>
#include<maya/MDagPath.h>
#include<maya/MDagPathArray.h>

#include<maya/MAnimControl.h>

#include<maya/MFnAttribute.h>
#include<maya/MFnNumericAttribute.h>
#include<maya/MFnTypedAttribute.h>
#include<maya/MArgDatabase.h>

#include<maya/MFnMesh.h>
#include<maya/MItMeshPolygon.h>

#include<maya/MFileIO.h>

#include<maya/MSyntax.h>

//OSD
#include <osd/vertex.h>
#include <osd/mesh.h>
#include <osd/cpuComputeContext.h>
#include <osd/cpuComputeController.h>
#include <osd/cpuVertexBuffer.h>
#include <far/meshFactory.h>
#include <far/dispatcher.h>

//Local
#include "../regression/common/shape_utils.h"
#include "eson.h"
#include "../raytrace/convert_bezier.h"
#include "mayaMeshToHbrMesh.h"
#include "mayaEsonWriter.h"
#include "simpleVertexXYZ.h"

struct xyzVV;

#if 0
typedef OpenSubdiv::OsdVertex VertexType;
#else
typedef SimpleVertexXYZf VertexType;
#endif

typedef OpenSubdiv::HbrMesh<VertexType>               HMesh;
typedef OpenSubdiv::HbrFace<VertexType>               HFace;
typedef OpenSubdiv::HbrVertex<VertexType>             HVertex;
typedef OpenSubdiv::HbrHalfedge<VertexType>           HHalfedge;
typedef OpenSubdiv::HbrFVarData<VertexType>           HFvarData;
typedef OpenSubdiv::HbrCatmarkSubdivision<VertexType> HCatmark;

typedef OpenSubdiv::FarMesh<VertexType>               FMesh;
typedef OpenSubdiv::FarMeshFactory<VertexType>        FMeshFactory;
typedef OpenSubdiv::FarSubdivisionTables        FSubdivision;
typedef OpenSubdiv::FarPatchTables              FPatches;


// ====================================
// MayaEsonWriter
// ====================================
MayaEsonWriter::MayaEsonWriter( const MString& options ){

	this->behavior = MayaEsonWriter::generateBehaviorFromOptionString( options );
}
// - - - - - - - - - - - - - - - - - -
MayaEsonWriter::~MayaEsonWriter(){
	
}

// - - - - - - - - - - - - - - - - - -
MayaEsonWriter::Behavior MayaEsonWriter::generateBehaviorFromOptionString( const MString& options ){

	MStringArray	optionList;
	MStringArray	theOption;

	options.split(';', optionList);

	Behavior retBehavior;

	for( unsigned int cnt = 0; cnt < optionList.length(); ++cnt ){
		theOption.clear();
		optionList[cnt].split( '=', theOption );
		if( theOption.length() < 2 ){
			continue;
		}

		if( theOption[0] == MString( "exportSel" ) && theOption[1].asInt() ){
			retBehavior.exportSelOnly = true;
		}
		
	}

	return retBehavior;
}

// - - - - - - - - - - - - - - - - - -
void MayaEsonWriter::doWrite( const MFileObject& file ){

    using namespace OpenSubdiv;
	MStatus status;

	//- - - - - - - - - - - - - - - - - -
	//DagAuxUtl
	//@todo:extract to toplevel class
	//- - - - - - - - - - - - - - - - - -
	struct DagAuxUtl{
		inline static MObject getLastNode( const MDagPath& dg ){
			if( dg.childCount() > 0 ){
				return dg.child( dg.childCount() );
			}else{
				return dg.node();
			}
		};
		
		//- - - - - - - - - - - - - - - - - -
		//
		inline static MStatus getLastDagPath(		MDagPath*			retDagPath,
													MObject*			retObj,
													const MObject&		argObj ){
			MStatus		status;

			if( !retDagPath || !retObj ){
				return MS::kFailure;
			}

			//get Last DagPath
			MFnDagNode	fnTempMeshDagNode( argObj, &status );

			MDagPath	tmpDagPath;
			fnTempMeshDagNode.getPath( tmpDagPath );

			*retObj = getLastNode( tmpDagPath );
			MFnDagNode fnMeshDagNode( *retObj );
			fnMeshDagNode.getPath( *retDagPath );

			return MS::kSuccess;
		};
	};//DagAuxUtl
	
	//Traverse nodes
	MSelectionList	activeSelList;
	MGlobal::getActiveSelectionList( activeSelList );
	
	for( MItDependencyNodes	depIt( MFn::kMesh ); !depIt.isDone(); depIt.next() ){
		MObject		itObj( depIt.item() );

		//need DagPath
		MDagPath dagPath;
		MObject nodeObj;
		
		#if 0
		MFnDagNode	fnTempMeshDagNode( itObj, &status );
		fnTempMeshDagNode.getPath( dagPath );
		#else
		if( DagAuxUtl::getLastDagPath( &dagPath, &nodeObj, itObj ) != MS::kSuccess ){
			continue;
		}
		#endif

		//selective export
		if( behavior.exportSelOnly && !activeSelList.hasItem( dagPath ) ){
			continue;
		}

		//fnMesh
		MFnMesh	fnMesh( dagPath, &status );
		if( status != MS::kSuccess ){
			continue;
		}
		#if 0
		MItMeshPolygon itMeshPoly( itObj, &status );
		#else
		MItMeshPolygon itMeshPoly( nodeObj, &status );
		#endif

		if( status != MS::kSuccess ){
			continue;
		}

		//build osdMesh and write
		{
			std::vector<int> fvarIndices;
			std::vector<int> fvarWidths;
			float maxCreaseSharpness = 0.0;

			MayaMeshToHbrMesh<VertexType >	hbrMeshGenerator;
			HMesh* hbrMesh = hbrMeshGenerator( fnMesh, itMeshPoly, fvarIndices, fvarWidths, &maxCreaseSharpness );

			//set SubDiv params
			HMesh::InterpolateBoundaryMethod vertInterpBoundaryMethod = HMesh::k_InterpolateBoundaryNone;
			HMesh::InterpolateBoundaryMethod fvarInterpBoundaryMethod = HMesh::k_InterpolateBoundaryNone;
			HCatmark::CreaseSubdivision creaseMethod = HCatmark::k_CreaseChaikin;
			HCatmark::TriangleSubdivision triangleSubdivision = HCatmark::k_New;
			HCatmark *catmarkSubdivision = dynamic_cast<HCatmark *>(hbrMesh->GetSubdivision());
			if (catmarkSubdivision) {
				catmarkSubdivision->SetTriangleSubdivisionMethod(triangleSubdivision);
			}
			hbrMesh->Finish();

			int levels = 3;
			bool isAdaptive = true;
			FMeshFactory meshFactory( hbrMesh, levels, isAdaptive );
			FMesh* farMesh = meshFactory.Create();
			static OpenSubdiv::FarComputeController computeController;
			computeController.Refine( farMesh );


			// centering/normalize vertices.
			std::vector<float> vertices;
			{
				#if 0
				float min[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
				float max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
				for (int i = 0; i < (int)farMesh->GetVertices().size(); ++i) {
					auto const &v = farMesh->GetVertices()[i];
					for (int j = 0; j < 3; ++j) {
						min[j] = std::min(min[j], v.GetPos()[j]);
						max[j] = std::max(max[j], v.GetPos()[j]);
					}
				}
				float center[3] = { (max[0]+min[0])*0.5f,
					(max[1]+min[1])*0.5f,
					(max[2]+min[2])*0.5f };
				float radius = std::max(std::max(max[0]-min[0], max[1]-min[1]), max[2]-min[2]);
				for (int i = 0; i < (int)farMesh->GetVertices().size(); ++i) {
					auto v = farMesh->GetVertices()[i];
					float x = v.GetPos()[0];
					float y = v.GetPos()[1];
					float z = v.GetPos()[2];
					vertices.push_back((x-center[0])/radius);
					vertices.push_back((y-center[1])/radius);
					vertices.push_back((z-center[2])/radius);
				}
				#else
				for (int i = 0; i < (int)farMesh->GetVertices().size(); ++i) {
					auto v = farMesh->GetVertices()[i];
					vertices.push_back( v.GetPos()[0] );
					vertices.push_back( v.GetPos()[1] );
					vertices.push_back( v.GetPos()[2] );
				}
				#endif
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

			std::string filename = file.rawFullName().asChar();
			FILE* fp = fopen(filename.c_str(), "wb");
			fwrite(&buf[0], 1, size, fp);
			fclose(fp);

			delete hbrMesh;
			delete farMesh;

		}
	}

}
