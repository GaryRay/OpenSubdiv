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

#include "mayaEsonWriter.h"


struct xyzVV;

#if 0
typedef OpenSubdiv::OsdVertex VertexType;
#else
typedef xyzVV VertexType;
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

// =============================================================================================================================
// from osd2eson

//------------------------------------------------------------------------------
// Vertex class implementation
struct xyzVV {

    xyzVV() { }

    xyzVV( int /*i*/ ) { }

    xyzVV( float x, float z, float y ) { _pos[0]=x; _pos[1]=y; _pos[2]=z; }

    xyzVV( const xyzVV & src ) { _pos[0]=src._pos[0]; _pos[1]=src._pos[1]; _pos[2]=src._pos[2]; }

   ~xyzVV( ) { }

    void AddWithWeight(const xyzVV& src, float weight) { 
        _pos[0]+=weight*src._pos[0]; 
        _pos[1]+=weight*src._pos[1]; 
        _pos[2]+=weight*src._pos[2]; 
    }

    void AddVaryingWithWeight(const xyzVV& , float) { }

    void Clear( void * =0 ) { _pos[0]=_pos[1]=_pos[2]=0.0f; }

    void SetPosition(float x, float z, float y) { _pos[0]=x; _pos[1]=y; _pos[2]=z; }

    void ApplyVertexEdit(const OpenSubdiv::HbrVertexEdit<xyzVV> & edit) {
        const float *src = edit.GetEdit();
        switch(edit.GetOperation()) {
          case OpenSubdiv::HbrHierarchicalEdit<xyzVV>::Set:
            _pos[0] = src[0];
            _pos[1] = src[1];
            _pos[2] = src[2];
            break;
          case OpenSubdiv::HbrHierarchicalEdit<xyzVV>::Add:
            _pos[0] += src[0];
            _pos[1] += src[1];
            _pos[2] += src[2];
            break;
          case OpenSubdiv::HbrHierarchicalEdit<xyzVV>::Subtract:
            _pos[0] -= src[0];
            _pos[1] -= src[1];
            _pos[2] -= src[2];
            break;
        }
    }

    void ApplyVertexEdit(OpenSubdiv::FarVertexEdit const & edit) {
        const float *src = edit.GetEdit();
        switch(edit.GetOperation()) {
          case OpenSubdiv::FarVertexEdit::Set:
            _pos[0] = src[0];
            _pos[1] = src[1];
            _pos[2] = src[2];
            break;
          case OpenSubdiv::FarVertexEdit::Add:
            _pos[0] += src[0];
            _pos[1] += src[1];
            _pos[2] += src[2];
            break;
        }
    }
    
    void ApplyMovingVertexEdit(const OpenSubdiv::HbrMovingVertexEdit<xyzVV> &) { }

    const float * GetPos() const { return _pos; }

private:
    float _pos[3];
};
// =============================================================================================================================

// =============================================================================================================================
// from osdPolySmooth.cpp
// @todo: rewrite and split to independent source
// =============================================================================================================================

// ====================================
// Macros
// ====================================
#define MCHECKERR(status,message)       \
    if( MStatus::kSuccess != status ) {   \
        cerr << "ERROR: " << message << "[" << status << "]" << endl;        \
        return status;                    \
    }

#define MWARNERR(status,message)       \
    if( MStatus::kSuccess != status ) {   \
        cerr << "ERROR: " << message << "[" << status << "]" << endl;        \
    }

//

// Reference: OSD shape_utils.h:: applyTags() "crease"
float
applyCreaseEdges(MFnMesh const & inMeshFn, HMesh * hbrMesh) {

    MStatus returnStatus;
    MUintArray tEdgeIds;
    MDoubleArray tCreaseData;
    float maxCreaseValue = 0.0f;

    if (inMeshFn.getCreaseEdges(tEdgeIds, tCreaseData)) {

        assert( tEdgeIds.length() == tCreaseData.length() );

        // Has crease edges
        int2 edgeVerts;
        for (unsigned int j=0; j < tEdgeIds.length(); j++) {

            // Get vert ids from maya edge
            int edgeId = tEdgeIds[j];
            returnStatus = inMeshFn.getEdgeVertices(edgeId, edgeVerts);

            // Assumption: The OSD vert ids are identical to those of the Maya mesh
            HVertex const * v = hbrMesh->GetVertex( edgeVerts[0] ),
                          * w = hbrMesh->GetVertex( edgeVerts[1] );

            HHalfedge * e = 0;
            if( v and w ) {

                if( (e = v->GetEdge(w)) == 0) {
                    e = w->GetEdge(v);
                }

                if(e) {
                    assert( tCreaseData[j] >= 0.0 );
                    e->SetSharpness( (float)tCreaseData[j] );

                    maxCreaseValue = std::max(float(tCreaseData[j]), maxCreaseValue);
                } else {
                    fprintf(stderr,
                        "warning: cannot find edge for crease tag (%d,%d)\n",
                            edgeVerts[0], edgeVerts[1] );
                }
            }
        }
    }
    return maxCreaseValue;
}


// Reference: OSD shape_utils.h:: applyTags() "corner"
float
applyCreaseVertices( MFnMesh const & inMeshFn, HMesh * hbrMesh ) {

    MUintArray tVertexIds;
    MDoubleArray tCreaseData;
    float maxCreaseValue = 0.0f;

    if ( inMeshFn.getCreaseVertices(tVertexIds, tCreaseData) ) {

        assert( tVertexIds.length() == tCreaseData.length() );

        // Has crease vertices
        for (unsigned int j=0; j < tVertexIds.length(); j++) {

            // Assumption: The OSD vert ids are identical to those of the Maya mesh

            HVertex * v = hbrMesh->GetVertex( tVertexIds[j] );
            if(v) {

                assert( tCreaseData[j] >= 0.0 );

                v->SetSharpness( (float)tCreaseData[j] );

                maxCreaseValue = std::max(float(tCreaseData[j]), maxCreaseValue);
            } else {
                fprintf(stderr,
                    "warning: cannot find vertex for corner tag (%d)\n",
                        tVertexIds[j] );
           }
        }
    }
    return maxCreaseValue;
}

// Collect UVs and ColorSet info to represent them as face-varying in OpenSubdiv
MStatus
getMayaFvarFieldParams(
    MFnMesh const & inMeshFn,
    MStringArray & uvSetNames,
    MStringArray & colorSetNames,
    std::vector<int> & colorSetChannels,
    std::vector<MFnMesh::MColorRepresentation> &colorSetReps,
    int & totalColorSetChannels) {

    MStatus returnStatus;

    returnStatus = inMeshFn.getUVSetNames(uvSetNames);
    MCHECKERR(returnStatus, "Cannot get uvSet names");

    returnStatus = inMeshFn.getColorSetNames(colorSetNames);
    MCHECKERR(returnStatus, "Cannot get colorSet names");

    colorSetChannels.resize(colorSetNames.length());
    colorSetReps.resize(colorSetNames.length());
    totalColorSetChannels = 0;

    for (unsigned int i=0; i < colorSetNames.length(); i++) {

        colorSetReps[i] = inMeshFn.getColorRepresentation(colorSetNames[i], &returnStatus);
        MCHECKERR(returnStatus, "Cannot get colorSet representation");

               if (colorSetReps[i] == MFnMesh::kAlpha) {
            colorSetChannels[i] = 1;
        } else if (colorSetReps[i] == MFnMesh::kRGB) {
            colorSetChannels[i] = 3;
        } else {
            colorSetChannels[i] = 4; // kRGBA
        }
        totalColorSetChannels += colorSetChannels[i];
    }
    return MS::kSuccess;
}

//! Create OSD HBR mesh.
//! Caller is expected to delete the resulting hbrMesh returned
HMesh *
createOsdHbrFromPoly( MFnMesh const & inMeshFn,
                      MItMeshPolygon & inMeshItPolygon,
                      std::vector<int> & fvarIndices,
                      std::vector<int> & fvarWidths,
                      float * maxCreaseSharpness=0)
{
    MStatus returnStatus;

    // == Mesh Properties

    // Gather FVarData
    MStringArray uvSetNames;
    MStringArray colorSetNames;
    std::vector<int> colorSetChannels;
    std::vector<MFnMesh::MColorRepresentation> colorSetReps;
    int totalColorSetChannels = 0;
    returnStatus = getMayaFvarFieldParams(inMeshFn, uvSetNames, colorSetNames, colorSetChannels, colorSetReps, totalColorSetChannels);
    MWARNERR(returnStatus, "Failed to retrieve Maya face-varying parameters");

    // Create face-varying data with independent float channels of dimension 1
    // Note: This FVarData needs to be kept around for the duration of the HBR mesh
    int totalFvarWidth = 2*uvSetNames.length() + totalColorSetChannels;
    fvarIndices.resize(totalFvarWidth);
    fvarWidths.resize(totalFvarWidth);
    for (int i=0; i < totalFvarWidth; ++i) {
        fvarIndices[i] = i;
        fvarWidths[i] = 1;
    }

    // temp storage for UVs and ColorSets for a face
    MIntArray fvArray; // face vertex array
    std::vector<MFloatArray> uvSet_uCoords(uvSetNames.length());
    std::vector<MFloatArray> uvSet_vCoords(uvSetNames.length());
    std::vector<MColorArray> colorSet_colors(colorSetNames.length());

    // =====================================
    // Init HBR
    // =====================================

    // Determine HBR Subdiv Method
    static HCatmark _catmark;

    // Create HBR mesh
    assert(fvarIndices.size() == fvarWidths.size());
    HMesh * hbrMesh = new HMesh( &_catmark,
                                 (int)fvarIndices.size(),
                                 (fvarIndices.size() > 0) ? &fvarIndices[0] : NULL,
                                 (fvarWidths.size()  > 0) ? &fvarWidths[0] : NULL,
                                 totalFvarWidth );

    // Create Stub HBR Vertices
    int numVerts = inMeshFn.numVertices();

	#if 1
	MPointArray points;
	inMeshFn.getPoints( points );
    for(int i=0; i<numVerts; i++ ){
	    VertexType v( (float)points[i][0], (float)points[i][1], (float)points[i][2] );
        hbrMesh->NewVertex( i, v );
    }
	#else
    for(int i=0; i<numVerts; i++ ) {
        hbrMesh->NewVertex( i, v );
    }
	#endif

    // ===================================================
    // Create HBR Topology
    // ===================================================

    assert(totalFvarWidth == hbrMesh->GetTotalFVarWidth());
    unsigned int ptxidx = 0;

    for( inMeshItPolygon.reset(); !inMeshItPolygon.isDone(); inMeshItPolygon.next() ) {

        // Verts for this face
        inMeshItPolygon.getVertices(fvArray);
        unsigned int nv = fvArray.length();

        bool valid = true;

        for(unsigned int j=0;j<fvArray.length(); j++) {

            HVertex const * origin      = hbrMesh->GetVertex( fvArray[j] ),
                          * destination = hbrMesh->GetVertex( fvArray[(j+1)%nv] );
            HHalfedge const * opposite = destination->GetEdge(origin);

            if(origin==NULL || destination==NULL) {
                fprintf(stderr, "Skipping face: An edge was specified that connected a nonexistent vertex\n");
                valid = false;
                break;
            }

            if(origin == destination) {
                fprintf(stderr, "Skipping face: An edge was specified that connected a vertex to itself\n");
                valid = false;
                break;
            }

            if(opposite && opposite->GetOpposite() ) {
                fprintf(stderr, "Skipping face: A non-manifold edge incident to more than 2 faces was found\n");
                valid = false;
                break;
            }

            if(origin->GetEdge(destination)) {
                fprintf(stderr, "Skipping face: An edge connecting two vertices was specified more than once."
                       " It's likely that an incident face was flipped\n");
                valid = false;
                break;
            }
        }

        // Update faces
        const int * fvArrayPtr = &(fvArray[0]); // get the raw int* array from std::vector<int>
        OpenSubdiv::HbrFace<VertexType> * face = hbrMesh->NewFace(nv, fvArrayPtr, 0);

        // Increment ptex index
        face->SetPtexIndex(ptxidx);

        // Add FaceVaryingData (UVSets, ...)
        if (totalFvarWidth > 0) {

            // Retrieve all UV and ColorSet data
            for (unsigned int i=0; i < uvSetNames.length(); ++i) {
                inMeshItPolygon.getUVs(uvSet_uCoords[i], uvSet_vCoords[i], &uvSetNames[i] );
            }
            for (unsigned int i=0; i < colorSetNames.length(); ++i) {
                inMeshItPolygon.getColors(colorSet_colors[i], &colorSetNames[i]);
            }

            std::vector<float> fvarItem(totalFvarWidth); // storage for all the face-varying channels for this face-vertex

            // loop over each uvSet and the uvs within
            for (unsigned int fvid=0; fvid < fvArray.length(); ++fvid) {

                int fvarItemIndex = 0;
                // Handle uvSets
                for( unsigned int uvSetIt=0; uvSetIt < uvSetNames.length(); ++uvSetIt ) {
                    if (fvid < uvSet_uCoords[uvSetIt].length()) {
                        fvarItem[fvarItemIndex  ] = uvSet_uCoords[uvSetIt][fvid];
                        fvarItem[fvarItemIndex+1] = uvSet_vCoords[uvSetIt][fvid];
                    } else {
                        // getUVs() can return incomplete or empty arrays
                        fvarItem[fvarItemIndex  ] = 0.0f;
                        fvarItem[fvarItemIndex+1] = 0.0f;
                    }
                    fvarItemIndex += 2;
                }
                // Handle colorSets
                for( unsigned int colorSetIt=0; colorSetIt < colorSetNames.length(); ++colorSetIt ) {

                    int nchannels = colorSetChannels[colorSetIt];
                    for (int channel=0; channel < nchannels; ++channel) {
                        if (fvid < colorSet_colors[colorSetIt].length()) {
                            fvarItem[fvarItemIndex+channel] = colorSet_colors[colorSetIt][fvid][channel];
                        } else {
                            // getColors() can return incomplete or empty arrays
                            fvarItem[fvarItemIndex+channel] = 0.0f;
                        }
                    }
                    fvarItemIndex += nchannels;
                }
                assert((fvarItemIndex) == totalFvarWidth); // For UVs, sanity check the resulting value

                // Insert into the HBR structure for that face
                HVertex * hbrVertex = hbrMesh->GetVertex( fvArray[fvid] );
                HFvarData & fvarData = hbrVertex->GetFVarData(face);

                if (not fvarData.IsInitialized()) {
                    fvarData.SetAllData(totalFvarWidth, &fvarItem[0]);
                } else if (not fvarData.CompareAll(totalFvarWidth, &fvarItem[0])) {

                    // If data exists for this face vertex, but is different
                    // (e.g. we're on a UV seam) create another fvar datum
                    OpenSubdiv::HbrFVarData<VertexType> &fvarData_new = hbrVertex->NewFVarData(face);
                    fvarData_new.SetAllData( totalFvarWidth, &fvarItem[0] );
                }
            }
        }

        // Increment ptxidx and store off ptex index values
        //   The number of ptexIds needed is 1 if a quad.  Otherwise it is the number of
        //   vertices for the face.
        int numPtexIdsForFace;
        if (valid) {
            numPtexIdsForFace = ( nv != 4 ) ? nv : 1 ;
        }
        else {
            numPtexIdsForFace = 0;
        }
        ptxidx += numPtexIdsForFace;
    }

    // Apply Creases
    float maxEdgeCrease = applyCreaseEdges( inMeshFn, hbrMesh );
    float maxVertexCrease = applyCreaseVertices( inMeshFn, hbrMesh );

    if (maxCreaseSharpness) {
        *maxCreaseSharpness = std::max(maxEdgeCrease, maxVertexCrease);
    }

    // Return the resulting HBR Mesh
    // Note that boundaryMethods and hbrMesh->Finish() still need to be called
    return hbrMesh;
}

// =============================================================================================================================


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

			HMesh *hbrMesh = createOsdHbrFromPoly( fnMesh, itMeshPoly, fvarIndices, fvarWidths, &maxCreaseSharpness );
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
