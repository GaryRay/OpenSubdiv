//Maya
#include <maya/MString.h>
#include <maya/MFileObject.h>
//Local

#if !defined( MAYA_ESON_WRITER_INCLUDED__ )
#define MAYA_ESON_WRITER_INCLUDED__

//Forward declarations

// ====================================
// MayaEsonWriter class
// ====================================
class MayaEsonWriter{
public:
	struct Behavior{
		inline Behavior(): exportSelOnly( false ){};
		bool		exportSelOnly;
		MString		outputBasePath;
	};
	
public:
	//creators
	MayaEsonWriter( const MString& options );
	virtual ~MayaEsonWriter();

public:
	//manipulators
	void doWrite( const MFileObject& file );
	static Behavior generateBehaviorFromOptionString( const MString& options );

protected:
	Behavior	behavior;

};//MayaEsonWriter

#endif
