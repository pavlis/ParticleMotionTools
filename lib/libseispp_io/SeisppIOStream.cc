#include <string>
#include <iostream>
#include <fstream>
#include "seispp.h"
#include "seispp_io.h"
using namespace std;
TextIOStreamReader::TextIOStreamReader()
{
  input_is_stdio=true;
  /* Only one object is allowed for stdin*/
  nobjects=1;
  parent_filename="STDIN";
  n_previously_read=0;
  ar=new boost::archive::text_iarchive(std::cin);
}
TextIOStreamReader::TextIOStreamReader(string fname)
{
  try{
    const string base_error("TextIOStreamReader file constructor:  ");
    ifs.open(fname.c_str(),ios::in);
    if(ifs.fail())
    {
      throw SeisppError(base_error+"cannot open file "+fname+" for input");
    }
    parent_filename=fname;
    input_is_stdio=false;
    n_previously_read=0;
    ifs.seekg(ifs.end-TextIOStreamEOFOffset);
    string magic_test;
    ifs >> magic_test;
    if(magic_test!=magic_tag) throw SeisppError(base_error + "File "
        + fname + " does not appear to be a valid seispp boost serialization file");
    ifs >> nobjects;
    ifs.seekg(0,ios_base::beg);
    ar=new boost::archive::text_iarchive(ifs);
  }catch(...){throw;};
}
/*
TextIOStreamReader::TextIOStreamReader(const TextIOStreamReader& parent) : ar(parent.ar)
{
  parent_filename=parent.parent_filename;
  input_is_stdio=parent.input_is_stdio;
  nobjects=parent.nobjects;
  n_previously_read=parent.n_previously_read;
  ifs=parent.ifs;
  ar=parent.ar;
}
TextIOStreamReader& TextIOStreamReader::operator=(const TextIOStreamReader& parent)
{
  if(this!=&parent)
  {
    parent_filename=parent.parent_filename;
    input_is_stdio=parent.input_is_stdio;
    nobjects=parent.nobjects;
    n_previously_read=parent.n_previously_read;
    ifs=parent.ifs;
    ar=parent.ar;
  }
  return *this;
}
*/
TextIOStreamReader::~TextIOStreamReader()
{
  delete ar;
  if(!input_is_stdio) ifs.close();
}
template <class InputObject> InputObject TextIOStreamReader::read()
{
  const string base_error("TextIOStreamReader read method:  ");
  InputObject d;
  try{

    if(n_previously_read>=(nobjects-1)) throw SeisppError(base_error
        + "Trying to read past end of file - code should test for this condition with at_eof method");
    (*ar)>>d;
    ++n_previously_read;
    return d;
  }catch(...)
  {
    throw SeisppError(base_error
      + "boost text serialization read failed\nCheck that input is a valid boost text serialization file");
  }
}
void TextIOStreamReader::set_method(string name)
{
  /* This does nothing.  Purely a place holder to link to base until I
  discover something more appropriate.*/
}
bool TextIOStreamReader::good()
{
  if(n_previously_read<(nobjects-1))
    return true;
  else
    return true;
}
bool TextIOStreamReader::eof()
{
  if(n_previously_read>=nobjects)
    return true;
  else
    return false;
}
TextIOStreamWriter::TextIOStreamWriter()
{
  output_is_stdio=true;
  nobjects=0;
  parent_filename="STDIN";
  ar=new boost::archive::text_oarchive(std::cout);
}
TextIOStreamWriter::TextIOStreamWriter(string fname)
{
  try{
    const string base_error("TextIOStreamWriter file constructor:  ");
    ofs.open(fname.c_str(),ios::out | ios::trunc);
    if(ofs.fail())
    {
      throw SeisppError(base_error+"open failed on file "+fname+" for output");
    }
    parent_filename=fname;
    output_is_stdio=false;
    nobjects=0;
    ar=new boost::archive::text_oarchive(ofs);
  }catch(...){throw;};
}
TextIOStreamWriter::~TextIOStreamWriter()
{
  char *buf=new char [TextIOStreamEOFOffset];
  sprintf(buf,"%s %ld",magic_tag.c_str(),nobjects);
  int i;
  for(i=0;i<TextIOStreamEOFOffset;++i) ofs<<buf[i];
  ofs.close();
  delete ar;
}
template <class OutputObject> void TextIOStreamWriter::write(OutputObject& d)
{
    try {
        (*ar) << d;
        ++nobjects;
    }catch(...)
    {
        throw SeisppError(string("TextIOStreamWriter write method failed\n")
                +"Is serialization defined for this object type?\n"
                +"Do you have write permission for output directory?");
    }
}
/* Binary versions go next - silly to do them until we at least get
the painfully similar text versions above to at least compile */
