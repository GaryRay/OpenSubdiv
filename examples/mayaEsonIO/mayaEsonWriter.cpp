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
#include "mayaFvarDataDesc.h"

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



namespace MayaAuxUtl{
	
	// ====================================
	// GetAttr
	// ====================================
	template< typename T >
	T GetAttr(	MFnMesh& fnMesh,
				MString attrName )
	{
		T ret;
		MStatus status;
		MString cmd = MString( "getAttr " ) + fnMesh.fullPathName() + "." + attrName;
		status = MGlobal::executeCommand( cmd, ret, false, false );

		if( status != MS::kSuccess ){
			throw std::exception("getAttr failed");
		}
		
		return ret;
	}

}//namespace MayaAuxUtl

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

	//Traverse nodes
	MSelectionList	activeSelList;
	MGlobal::getActiveSelectionList( activeSelList );

	int exportCnt = 0;

	for( MItDag	dagIt( MItDag::kDepthFirst, MFn::kMesh, &status ); !dagIt.isDone(); dagIt.next() ){

		#if 1
		if( exportCnt > 0 ){
			break;
		}
		#endif

		//need DagPath
		MDagPath dagPath;
		status = dagIt.getPath( dagPath );
		MObject nodeObj = dagPath.node();

		//selective export
		if( behavior.exportSelOnly && !activeSelList.hasItem( nodeObj ) ){
			continue;
		}

		//fnMesh
		MFnMesh	fnMesh( dagPath, &status );
		if( status != MS::kSuccess ){
			continue;
		}
		MItMeshPolygon itMeshPoly( nodeObj, &status );

		if( status != MS::kSuccess ){
			continue;
		}

		//check current mesh is Subdivisition surface or not.
		int dispSmoothMeshVal = 0;
		int smoothLevel = 2;
		try{
			dispSmoothMeshVal = MayaAuxUtl::GetAttr<int>( fnMesh, "displaySmoothMesh" );
			smoothLevel = MayaAuxUtl::GetAttr<int>( fnMesh, "smoothLevel" );
		}catch(...){
			continue;
		}

		if( dispSmoothMeshVal == 0 || smoothLevel == 0 ){
			//@Todo: 'polygon mesh' eson format export here

		}else if( behavior.useOSDEson ){
			//build osdMesh and write
			std::vector<int> fvarIndices;
			std::vector<int> fvarWidths;
			float maxCreaseSharpness = 0.0;
			MayaFVarDataDesc fvarDesc;

			//Generate
			MayaMeshToHbrMesh<VertexType >	hbrMeshGenerator;
			std::shared_ptr<HMesh> hbrMesh = std::shared_ptr<HMesh>( hbrMeshGenerator(
				fnMesh, itMeshPoly, fvarIndices, fvarWidths, &fvarDesc, &maxCreaseSharpness ) );
			if( !hbrMesh ){
				continue;
			}

			//SubDiv params var.
			HMesh::InterpolateBoundaryMethod vertInterpBoundaryMethod = HMesh::k_InterpolateBoundaryEdgeAndCorner;
			HMesh::InterpolateBoundaryMethod fvarInterpBoundaryMethod = HMesh::k_InterpolateBoundaryEdgeAndCorner;
			HCatmark::CreaseSubdivision creaseMethod = HCatmark::k_CreaseNormal;
			HCatmark::TriangleSubdivision triangleSubdivision = HCatmark::k_Normal;
			HCatmark *catmarkSubdivision = dynamic_cast<HCatmark *>( hbrMesh->GetSubdivision() );
			bool fvarPropCorners = false;

			//Maya Attribute
			int smoothDrawType = 0;
			try{
				smoothDrawType = MayaAuxUtl::GetAttr<int>( fnMesh, "smoothDrawType" );
			}catch(...){
				//Maya2014 or below version will here and it it ok. smoothDrawType still 0.
			}
			
			if( smoothDrawType <= 1 ){
				//Maya Catmarll-Clark
				//What is type smoothDrawType == 1? I don't know... but process here.
				
				

			}else if( smoothDrawType == 2 ){
				//Maya2015 or above OSD Catmarll-Clark
				try{
					if( MayaAuxUtl::GetAttr<int>( fnMesh, "osdVertBoundary" ) == 1 ){
						vertInterpBoundaryMethod = HMesh::k_InterpolateBoundaryEdgeOnly;
					}

					switch( MayaAuxUtl::GetAttr<int>( fnMesh, "osdFvarBoundary" ) ){
						case 0:{
							fvarInterpBoundaryMethod = HMesh::k_InterpolateBoundaryNone;
							break;
						}
						case 1:{
							fvarInterpBoundaryMethod = HMesh::k_InterpolateBoundaryEdgeAndCorner;
							break;
						}
						case 2:{
							fvarInterpBoundaryMethod = HMesh::k_InterpolateBoundaryEdgeOnly;
							break;
						}
						case 3:{
							fvarInterpBoundaryMethod = HMesh::k_InterpolateBoundaryAlwaysSharp;
							break;
						}
						default:
						break;
					}
					
					if( MayaAuxUtl::GetAttr<int>( fnMesh, "osdSmoothTriangles" ) == 1 ){
						triangleSubdivision = HCatmark::k_New;
					}
					
					if( MayaAuxUtl::GetAttr<int>( fnMesh, "osdCreaseMethod" ) == 1 ){
						creaseMethod = HCatmark::k_CreaseNormal;
					}

				}catch(...){
				}
			}
			
			//set SubDiv params
			hbrMesh->SetInterpolateBoundaryMethod( vertInterpBoundaryMethod );
			hbrMesh->SetFVarInterpolateBoundaryMethod( fvarInterpBoundaryMethod );
			hbrMesh->SetFVarPropagateCorners( fvarPropCorners );
			hbrMesh->GetSubdivision()->SetCreaseSubdivisionMethod( creaseMethod );
			if( catmarkSubdivision ){
				catmarkSubdivision->SetTriangleSubdivisionMethod( triangleSubdivision );
			}
			hbrMesh->Finish();

			MObjectArray shaders;
			MIntArray shaderIndices;
			fnMesh.getConnectedShaders( 0, shaders, shaderIndices );
			std::vector<short> shaderIndicesVec;

			for( unsigned int j = 0; j < shaderIndices.length(); ++j ){
				shaderIndicesVec.push_back( shaderIndices[j] );
			}

			//Eson
			std::string filename = file.rawFullName().asChar();
			try{
				HbrMeshToEson<VertexType>()( hbrMesh.get(), fvarDesc, shaderIndicesVec, filename, smoothLevel );
			}catch(...){
				continue;
			}
			++exportCnt;
		}else{
			//
		}


	}

}
