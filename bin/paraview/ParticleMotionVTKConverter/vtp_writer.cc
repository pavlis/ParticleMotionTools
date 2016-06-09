#include <fstream>
#include <sstream>
#include "perf.h"
#include "ParticleMotionData.h"
using namespace std;
using namespace SEISPP;
void WriteTimeWindow(vector<ParticleMotionData>& d,TimeWindow tw,ofstream& out)
{
    try{
        int i;
        vector<ParticleMotionData>::iterator dptr;
        vector<dmatrix> points;
        vector<dmatrix>::iterator pptr;
        points.reserve(d.size());
        int TotalNumberPoints(0);
        int NumberPMCurves(0);
        double rmin(-1.0),rmax(0.0);  // These may not be necessary 
        for(dptr=d.begin();dptr!=d.end();++dptr)
        {
            try{
                dmatrix pm(dptr->particle_motion(tw));
                points.push_back(pm);
                TotalNumberPoints+=pm.columns();
                ++NumberPMCurves;
                if(rmin<0.0) rmin=dnrm2(3,pm.get_address(0,0),1);
                for(i=0;i<pm.columns();++i)
                {
                    double r=dnrm2(3,pm.get_address(0,i),1);
                    if(r>rmax) rmax=r;
                    if(r<rmin) rmin=r;
                }
            }catch(SeisppError& serr)
            {
                /* this call to get_string is a bit dangerous if
                   this is used for nonstandard data */
                cerr << "Error building particle motion data for "
                    << " station "<< dptr->get_string("sta")<<endl;
                serr.log_error();
                cerr << "Data not saved"<<endl;
            }
        }
        /* Write this preamble.  This is a serious maintenance issue
           as the file format could easily change. */
        out << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"
            <<endl
            <<"<PolyData>"<<endl
            << "<Piece NumberOfPoints=\""<<TotalNumberPoints<<"\" "
            << "NumberOfVerts=\"0\" NumberOfLines=\"" << NumberPMCurves
            << "\" " 
            << "NumberOfStrips=\"0\" NumberOfPolys=\"0\">" <<endl;
        /* Not sure these are necessary, but best leave them empty than
           not present I presume */
        out << "<PointData>\n</PointData>"<<endl;
        out << "<CellData>\n</CellData>"<<endl;
        out << "<Points>"<<endl;
        out << "<DataArray type=\"Float32\" Name=\"Points\""
                <<" NumberOfComponents=\"3\" format=\"ascii\""
                <<" RangeMin=\""<<rmin
                <<"\" RangeMax=\""<<rmax<<"\">" << endl;
        for(pptr=points.begin();pptr!=points.end();++pptr)
        {
            dmatrix pmt;
            pmt=tr(*pptr);
            out << pmt;
        }
        out << "</DataArray>"<<endl;
        out << "</Points>"<<endl;
        /* A Verts attribute may be required here, but I think we can 
           ignore it for pure lines */
        out << "<Lines>"<<endl;
        long int indexmin(0),indexmax;
        indexmax=TotalNumberPoints-1;
        long index(0);
        out << "<DataArray type=\"Int64\" Name=\"connectivity\" "
                << "format=\"ascii\" "
                << "RangeMin=\""<<indexmin<<"\" RangeMax=\""<<indexmax<<"\">"
                <<endl;
        for(pptr=points.begin();pptr!=points.end();++pptr)
        {
            indexmax=indexmin + pptr->columns() - 1;
            for(i=0;i<pptr->columns();++i)
            {
                out << index <<" ";
                if( ((i+1)%10)==0) out << endl;
                ++index;
            }
        }
        out << endl <<"</DataArray>"<<endl;
        /* This section defines the beginning of each line.  The count
           must match the NumberOfLines definition earlier. */
        out << "<DataArray type=\"Int64\" Name=\"offsets\" "
                << "format=\"ascii\" "
                << "RangeMin=\""<<indexmin<<"\" RangeMax=\""<<indexmax<<"\">"
                <<endl;
        /* This seems to demand to start at 1 instead of 0 */
        index=0;
        for(pptr=points.begin();pptr!=points.end();++pptr)
        {
            if(index==0)
                out << index+1 <<endl;
            else
                out << index <<endl;
            index+=pptr->columns();

        }
        out << "</DataArray>"<<endl;
        out << "</Lines>"<<endl;
        /* My example file has an empty Strips and Polys section.  Hopefully
           this can be avoided*/
        out << "</Piece>"<<endl
            << "</PolyData>"<<endl
            << "</VTKFile>"<<endl;
    }catch(...){throw;};
}


