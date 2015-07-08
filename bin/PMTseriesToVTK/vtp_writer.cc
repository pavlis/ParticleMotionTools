#include <fstream>
#include <sstream>
#include "dmatrix.h"
#include "perf.h"
#include "ParticleMotionEllipse.h"
#include "PMTimeSeries.h"
using namespace std;
using namespace SEISPP;
/* ParticleMotionEllipse object has a points method that returns
    a matrix of point defining an ellipse figure.   To close that figure
    we have to repeat the first point.  This is the return of this
    procedure - i.e. copy of dm but with one extra row repeating row 0.
 
 This procedure intentionally does to try to catch indexing errors in
 the matrix as this routine is assumed locked to this code where that
 cannot happen.*/
dmatrix close_ellipse(dmatrix& dm)
{
    int nr=dm.rows();
    dmatrix result(nr+1,3);
    /* We could do this more efficiently with memcp, but using the
       indexing is clearer. */
    int i,j;
    for(j=0;j<3;++j)
        for(i=0;i<nr;++i)
            result(i,j)=dm(i,j);
    for(j=0;j<3;++j) result(nr,j)=dm(0,j);
    return result;
}
void WriteEllipses(vector<PMTimeSeries>& d,
        double t, double scale, int np_per_ellipse,
        ofstream& out)
{
    try{
        int i;
        vector<PMTimeSeries>::iterator dptr;
        vector<dmatrix> points;
        vector<dmatrix>::iterator pptr;
        points.reserve(d.size());
        int TotalNumberPoints(0);
        int NumberPMCurves(0);
        double rmin(-1.0),rmax(0.0);  // These may not be necessary 
        for(dptr=d.begin();dptr!=d.end();++dptr)
        {
            try{
                i=dptr->sample_number(t);
                ParticleMotionEllipse pme;  // assume initializes to 0;
                try {
                    pme=dptr->ellipse(i);
                }catch(SeisppError& serr)
                {
                    cerr << "Warning:  time="<<t<<" is outside range of data for station="
                        << dptr->get_string("sta")<<endl
                        << "Message from ParticleMotionEllipse::ellipse method"
                        <<endl;
                    serr.log_error();
                    cerr << "Setting ellipse to zero size"<<endl;
                }
                dmatrix pm=pme.points(np_per_ellipse);
                /* Need to add a point to close the ellipse figure */
                pm=close_ellipse(pm);
                /* A simple way to scale the ellipse */
                pm=scale*pm;
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
            /* points are in a 3xN matrix so we can just write them
               using ostream operator. */
            out << *pptr;
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


