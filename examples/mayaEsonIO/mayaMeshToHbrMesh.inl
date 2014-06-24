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
#include<maya/MItMeshEdge.h>

//OSD
#include <osd/vertex.h>
#include <osd/mesh.h>
#include <osd/cpuComputeContext.h>
#include <osd/cpuComputeController.h>
#include <osd/cpuVertexBuffer.h>
#include <far/meshFactory.h>
#include <far/dispatcher.h>


// ====================================
// Macros
// ====================================
//
#define MCHECKERR(status,message)\
	if( MStatus::kSuccess != status ) {\
		cerr << "ERROR: " << message << "[" << status << "]" << endl;\
		return status;\
	}

#define MWARNERR(status,message)\
	if( MStatus::kSuccess != status ) {\
		cerr << "ERROR: " << message << "[" << status << "]" << endl;\
	}
// ====================================
// MayaMeshToHbrMesh
// ====================================
template< typename T >
typename MayaMeshToHbrMesh<T>::HMesh* MayaMeshToHbrMesh<T>::operator()(
	MFnMesh const & inMeshFn,
	MItMeshPolygon & inMeshItPolygon,
	std::vector<int> & fvarIndices,
	std::vector<int> & fvarWidths,
	float * maxCreaseSharpness )
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
	setHbrVertices( hbrMesh, inMeshFn );

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

// - - - - - - - - - - - - - - - - - -
template< typename T >
void MayaMeshToHbrMesh<T>::setHbrVertices(	HMesh* hbrMesh,
											 MFnMesh const & inMeshFn )
{
	int numVerts = inMeshFn.numVertices();
	VertexType v;
	for(int i=0; i<numVerts; i++ ) {
		hbrMesh->NewVertex( i, v );
	}
}

//Specialize for simpleVertexXYZ
template<>
void MayaMeshToHbrMesh< SimpleVertexXYZf >::setHbrVertices( HMesh* hbrMesh,
											 MFnMesh const & inMeshFn )
{
	int numVerts = inMeshFn.numVertices();
	MPointArray points;
	inMeshFn.getPoints( points );
	for(int i=0; i<numVerts; i++ ){
		VertexType v( (float)points[i][0], (float)points[i][1], (float)points[i][2] );
		hbrMesh->NewVertex( i, v );
	}
}

template<>
void MayaMeshToHbrMesh< SimpleVertexXYZd >::setHbrVertices( HMesh* hbrMesh,
											 MFnMesh const & inMeshFn )
{
	int numVerts = inMeshFn.numVertices();
	MPointArray points;
	inMeshFn.getPoints( points );
	for(int i=0; i<numVerts; i++ ){
		VertexType v( (double)points[i][0], (double)points[i][1], (double)points[i][2] );
		hbrMesh->NewVertex( i, v );
	}
}

// - - - - - - - - - - - - - - - - - -
template< typename T >
MStatus MayaMeshToHbrMesh<T>::getMayaFvarFieldParams(
	MFnMesh const & inMeshFn,
	MStringArray & uvSetNames,
	MStringArray & colorSetNames,
	std::vector<int> & colorSetChannels,
	std::vector<MFnMesh::MColorRepresentation> &colorSetReps,
	int & totalColorSetChannels)
{
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

// - - - - - - - - - - - - - - - - - -
template< typename T >
float MayaMeshToHbrMesh<T>::applyCreaseEdges(
	MFnMesh const & inMeshFn,
	HMesh * hbrMesh)
{

	// - - - - - - - - - - - - - - - - - -
	struct setCreaseFromMayaEdgeId{
		inline void operator()( int edgeId,
								MFnMesh const & inMeshFn,
								HMesh * hbrMesh,
								float* maxCreaseValue,
								float sharpVal
								)
		{
			int2 edgeVerts;
			
			MStatus returnStatus = inMeshFn.getEdgeVertices( edgeId, edgeVerts );

			// Assumption: The OSD vert ids are identical to those of the Maya mesh
			HVertex const *v = hbrMesh->GetVertex( edgeVerts[0] );
			HVertex const *w = hbrMesh->GetVertex( edgeVerts[1] );

			HHalfedge *e = 0;
			if( v and w ) {

				if( (e = v->GetEdge(w)) == 0) {
					e = w->GetEdge(v);
				}

				if( e ){
					e->SetSharpness( sharpVal );
					*maxCreaseValue = std::max<float>( sharpVal, *maxCreaseValue );
				} else {
					fprintf(stderr,
							"warning: cannot find edge for crease tag (%d,%d)\n",
							edgeVerts[0], edgeVerts[1] );
				}
			}
			
		};
	};
	// - - - - - - - - - - - - - - - - - -

	MStatus returnStatus;
	MUintArray tEdgeIds;
	MDoubleArray tCreaseData;
	float maxCreaseValue = 0.0f;

	if (inMeshFn.getCreaseEdges(tEdgeIds, tCreaseData)) {
		assert( tEdgeIds.length() == tCreaseData.length() );
		// Has crease edges
		for (unsigned int j=0; j < tEdgeIds.length(); j++) {
			setCreaseFromMayaEdgeId()( tEdgeIds[j], inMeshFn, hbrMesh, &maxCreaseValue, (float)( tCreaseData[j] ) );
		}
	}

	MItMeshEdge itMeshEdge( inMeshFn.object() );
	for( ;!itMeshEdge.isDone(); itMeshEdge.next() ){
		if( !itMeshEdge.onBoundary() ){
			continue;
		}
		int index = itMeshEdge.index();
		setCreaseFromMayaEdgeId()( index, inMeshFn, hbrMesh, &maxCreaseValue, 3.0f );
	}

	
	return maxCreaseValue;
}

// - - - - - - - - - - - - - - - - - -
template< typename T >
float MayaMeshToHbrMesh<T>::applyCreaseVertices(
	MFnMesh const & inMeshFn,
	HMesh * hbrMesh )
{

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

#undef MCHECKERR
#undef MWARNERR
