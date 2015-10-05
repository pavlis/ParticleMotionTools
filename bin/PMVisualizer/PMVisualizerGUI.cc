#include "PMVisualizerGUI.h"
#include "filter++.h"
PMVisualizerGUI::PMVisualizerGUI() : TimeWindowPicker(), guimd()
{
    filtername="DEMEAN";   //frozen default filter 
} 
PMVisualizerGUI::PMVisualizerGUI(Metadata& md) : TimeWindowPicker(md),guimd()
{
    filtername="DEMEAN";   //frozen default filter 
    winmd=md;
}
void PMVisualizerGUI::load_pfdata(PfStyleMetadata& md)
{
    guimd=md;
}
void PMVisualizerGUI::filter_and_plot(TimeSeriesEnsemble& din)
{
    try {
        list<string> filterlist;
        filterlist=guimd.get_tbl("BRTT_filter_definition");
        /* the filtername is a comma separated list so we have to reformat*/
        list<string>::iterator sptr;
        filtername.clear();
        for(sptr=filterlist.begin();sptr!=filterlist.end();++sptr)
        {
            if(sptr==filterlist.begin()) 
                filtername=(*sptr);
            else
                filtername=filtername+","+(*sptr);
        }
        TimeInvariantFilter f(filtername);
        TimeSeriesEnsemble d(din);
        FilterEnsemble(d,f);
        // Do not block - later call too pick time window is needed
        this->plot(d);
        //DEBUG - see if we can set display marker 
        /*
        this->setpick(TimeWindow(5.0,20.0));
        this->repick();
        */
    }catch(...){throw;};
}

    
