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
#include "hbrMeshToEson.h"
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
			std::shared_ptr<HMesh> hbrMesh = std::shared_ptr<HMesh>( hbrMeshGenerator( fnMesh, itMeshPoly, fvarIndices, fvarWidths, &maxCreaseSharpness ) );
			if( !hbrMesh ){
				continue;
			}

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

			//Eson
			std::string filename = file.rawFullName().asChar();
			try{
				HbrMeshToEson<VertexType>()( hbrMesh.get(), filename );
			}catch(...){
				continue;
			}

		}
	}

}
