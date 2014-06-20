//Local
#include "mayaEsonWriter.h"
#include "mayaEsonIO.h"

// ====================================
// MayaEsonIO
// ====================================
MayaEsonIO::MayaEsonIO(){
}

// - - - - - - - - - - - - - - - - - -
MayaEsonIO::~MayaEsonIO(){
}

// - - - - - - - - - - - - - - - - - -
void* MayaEsonIO::creator(){
	return new MayaEsonIO();
}

// - - - - - - - - - - - - - - - - - -
bool MayaEsonIO::haveWriteMethod() const{
	return true;
}

// - - - - - - - - - - - - - - - - - -
MString MayaEsonIO::filter() const{
	return "*.eson";
}

// - - - - - - - - - - - - - - - - - -
MString MayaEsonIO::defaultExtension() const{
	return "eson";
}

// - - - - - - - - - - - - - - - - - -
MPxFileTranslator::MFileKind MayaEsonIO::identifyFile(const MFileObject& file, const char* buffer, short size) const
{
	MFileKind ret;
	ret = kIsMyFileType;

	//@TODO: check filetype
	
	return ret;
}

// - - - - - - - - - - - - - - - - - -
MStatus MayaEsonIO::writer( const MFileObject& file, const MString& options, FileAccessMode mode){

	MayaEsonWriter	writerImp( options );
	try{
		writerImp.doWrite( file );
	}catch(...){
		return MS::kFailure;
	}
	
	return MS::kSuccess;
}

// - - - - - - - - - - - - - - - - - -
MStatus MayaEsonIO::reader( const MFileObject& file, const MString& options, FileAccessMode mode){
	
	return MS::kFailure;
	
}


