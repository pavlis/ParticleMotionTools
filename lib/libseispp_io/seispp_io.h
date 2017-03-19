#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
/*! \brief Abstract base class for generic object reader.

This is the base class for a collection of C++ classes designed to
abstract and simplify the process of reading and writing of seispp data
objects.   The purpose of this base class is to provide the foundation for
a set of polymorphic objects for reading data in various generic forms.
*/
class BasicObjectReader
{
  public:
    /*! Pure virtual method used to allow an abstract base class.*/
    virtual bool good()=0;
    virtual long number_available();
};
/*! \brief Abstract base class for generic object writer.

This is the base class for a collection of C++ classes designed to
abstract and simplify the process of reading and writing of seispp data
objects.   The purpose of this base class is to provide the foundation for
a set of polymorphic objects for saving data in varous forms.
*/
class BasicObjectWriter
{
  public:
    virtual long number_already_written();
};
/*! We write the number of objects in set of concatenated serial objects
   at the end of the file.  We seek back this many bytes to read the
   number of objects written */
const int TextIOStreamEOFOffset(64);
/*! This string is written to the very end of the file and tested on input to
validate file integrity.*/
const string magic_tag("SEISPP_FIlE_IS_VALID");
/*! \brief Generic object reader to read serialized test data from a sequential file.

This object abstracts reading of objects based on two standard concepts:
(1) access is sequential as in reading from a tape, and (2) data are
serialized and written as ascii text.   The implementation uses the
boost text archive library, but the design of the interface is intended to
insulate the application from this implementation detail. */
class TextIOStreamReader : BasicObjectReader
{
  public:
    /*! \brief Default constructor.

      Default constructor uses boost text serialization to stdin*/
    TextIOStreamReader();
    /*! \brief Create handle to read from file.

      Creates an input handle to read from a file.
      \param fname - file name opened as ifs

      \exception - throws a SeisppError object if operation fails.  Boost
        constructors may also throw special error object (needs research).
        */
    TextIOStreamReader(string fname);
    /*! Standard copy constructor.  Seems impossible with boost text_iarchive*/
    //TextIOStreamReader(const TextIOStreamReader& parent);
    /*! Destructor - has to close io channel */
     ~TextIOStreamReader();
     /*! Standard assignment operator. Appears impossible for boost serialization*/
    //TextIOStreamReader& operator=(const TextIOStreamReader& parent);
    template <class T> T read();
    /*! Returns number of objects in the file being read. */
    long number_available(){return nobjects;};
    /*! \brief Return the number of objects already read.

    In loops going through a sequential file it can be useful to know the
    position.   This effectively returns a count of position as an integer
    of how many objects are read.   Since most objects are not a fixed size
    this is not directly linked to file position.  It is of use in a loop
    of for information on how a job is progressing. */
    long number_already_read(){return n_previously_read;};
    /*! \brief Test if ok to read more.

    This method can be used to drive a while loop.  Returns true as long
    as the count of object read is less than the number in the file */
    bool good();
    /*! Test for end of file condition.  */
    bool eof();
    /*! \brief position to beginning of file.

    Since a serial file is basd on the concept of sequential access it is
    useful to have the concept of rewind as in a tape. */
    void rewind();
  private:
    boost::archive::text_iarchive *ar;
    /* input from stdio is special and is flagged by this boolean */
    bool input_is_stdio;
    /* This string defines the parent filename opened as ifstream ifs*/
    string parent_filename;
    /* stream linked to ar */
    ifstream ifs;
    /* To support multiple objects in one serial file we need to
       cache the number of objects expected and the number already
       read */
    long nobjects;
    long n_previously_read;
};
/*! \brief Generic object writer saving data in a sequential text file.

This object abstracts reading of objects based on two standard concepts:
(1) access is sequential as in reading from a tape, and (2) data are
serialized and written as ascii text.   The implementation uses the
boost text archive library, but the design of the interface is intended to
insulate the application from this implementation detail. */
class TextIOStreamWriter : BasicObjectWriter
{
  public:
    /*! \brief Default constructor.

      Default constructor uses boost text serialization to stdout*/
    TextIOStreamWriter();
    /*! \brief Create handle to write to file.

      Creates an input handle to write to a file.
      \param file - file to open for output

      \exception - throws a SeisppError object if operation fails.  Boost
        constructors may also throw special error object (needs research).
        */
    TextIOStreamWriter(string fname);
    /*! Destructor - has to close io channel */
     ~TextIOStreamWriter();
    /*! \brief write one object.

    This is the primary method of this object.  Writes a single object to
    the output.  Assumes d has a serialization defined.

    \param d - object to be written
    */
    template <class T> void write(T& d);
    long number_already_written(){return nobjects;};
  private:
    boost::archive::text_oarchive *ar;
    /* input from stdio is special and is flagged by this boolean */
    bool output_is_stdio;
    /* This is the name of the file linked to ofs.*/
    string parent_filename;
    /* stream linked to ar */
    ofstream ofs;
    /* To support multiple objects in one serial file we need to
       count the number of objects written.   This number is written to
       the end of the file.   It is incremted by each write. */
    long nobjects;
};
/***** Should be able to add: (1) read/write binary sequential, (2) read/write
xml sequential, (3) forms of indexed files with keys generated seperately, and
(4) database access through a nonsql db */

/*! Legacy writer for archive connected to stdout - original seispp_filters */
template <class OutputObject> void write_object(OutputObject& d,
        boost::archive::text_oarchive& oa)
{
    try {
        oa << d;
    }catch(...)
    {
        throw SeisppError(string("write_object failed\n")
                +"Is serialization defined for this object type?\n"
                +"Do you have write permission for output directory?");
    }
}
/*! Legacy reader for archive connected to stdin- original seispp_filters */

template <class InputObject> InputObject
    read_object(boost::archive::text_iarchive& ia)
{
    InputObject d;
    try{
        ia>>d;
    }catch(...)
    {
        /* This template ignores errors thrown by boost
           and converts to a standard error for the seispp
           library - perhaps not ideal but what is done here. */
        throw SeisppError(string("read_object failed:  ")
                    + "Check that input file is a boost text archive file");
    }
    return d;
}
