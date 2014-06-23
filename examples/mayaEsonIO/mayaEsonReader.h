//Maya
#include <maya/MString.h>

//Local

#if !defined( MAYA_ESON_READER_INCLUDED__ )
#define MAYA_ESON_READER_INCLUDED__

// ====================================
// MayaEsonReader class
// ====================================
class MayaEsonReader{
public:
	MayaEsonReader( const MString& options );
	virtual ~MayaEsonReader();

protected:
	

};//MayaEsonReader

#endif
