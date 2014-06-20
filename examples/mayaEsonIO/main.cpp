#include<maya/MFnPlugin.h>
#include<maya/MPxFileTranslator.h>
#include<maya/MPxCommand.h>

#include "mayaEsonIO.h"

// ====================================
// initializePlugin
MStatus initializePlugin( MObject obj ){
	MFnPlugin plugin(obj, "mayaEsonIO", "1.0", "any");
	
	plugin.registerFileTranslator( "mayaEsonIO",
									NULL,
									MayaEsonIO::creator,
									nullptr,
									nullptr,
									false);
	
	return MS::kSuccess;
}

// ====================================
// uninitializePlugin
MStatus uninitializePlugin( MObject obj ){
	MFnPlugin plugin(obj);
	plugin.deregisterFileTranslator("mayaEsonIO");
	
	return MS::kSuccess;
}

