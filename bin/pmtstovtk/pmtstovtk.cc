/* This program takes the output of mwpm, which computes particle motion
   esimates as a function of frequency and time, and spits out a vtk
   file that can be used with paraview to visualize the motion in 3D. */
void usage()
{
    cerr << "pmtstovtk [-pf pffile ] file1 file2 ... filen "<<endl;
}
int main(int argc, char **argv)
{
}
