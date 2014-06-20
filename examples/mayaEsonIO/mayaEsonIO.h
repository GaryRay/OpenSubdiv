//Maya
#include<maya/MPxFileTranslator.h>

#if !defined( MAYA_ESON_IO_INCLUDED__ )
#define MAYA_ESON_IO_INCLUDED__

//Forward declarations
class MayaEsonIO;


// ====================================
// MayaEsonIO class
// ====================================
class MayaEsonIO: public MPxFileTranslator{
public:
	//creators
	MayaEsonIO();
	virtual ~MayaEsonIO();
	
public:
	//manipulators
	bool haveWriteMethod( void ) const;
	MString filter() const;
	MString defaultExtension() const;
	
	MFileKind identifyFile( const MFileObject& file, const char* buffer, short size) const;
	MStatus writer( const MFileObject& file, const MString& options, FileAccessMode mode);
	MStatus reader( const MFileObject& file, const MString& options, FileAccessMode mode);

	static void* creator();

protected:
	//members
	MString					m_melModeStr;
	FileAccessMode			m_mode;
	MString					m_fileStr;
	
};//MayaEsonIO class



#endif
