#include "TimeWindowPicker.h"
using namespace SEISPP;
class PMVisualizerGUI : public TimeWindowPicker
{
    public:
        PMVisualizerGUI();
        PMVisualizerGUI(Metadata& md);
        void load_pfdata(PfStyleMetadata& md);
        /* Intentionally not a reference because copy is altered*/
        void filter_and_plot(TimeSeriesEnsemble& d);
        Metadata get_parameters()
        {
            return guimd;
        };
        string get_filter()
        {
            return filtername;
        };
    private:
        Metadata winmd;
        PfStyleMetadata guimd;
        string filtername;
};

