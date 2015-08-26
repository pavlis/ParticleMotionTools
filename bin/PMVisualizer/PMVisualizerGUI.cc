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
void PMVisualizerGUI::load_pfdata(Metadata& md)
{
    guimd=md;
}
void PMVisualizerGUI::filter_and_plot(TimeSeriesEnsemble d)
{
    try {
        filtername=guimd.get_string("BRTT_filter_definition");
        TimeInvariantFilter f(filtername);
        FilterEnsemble(d,f);
        this->plot(d);
    }catch(...){throw;};
}

    
