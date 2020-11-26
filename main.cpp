
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        BL_PROFILE("main");

        IntVect n_cell;
        int max_grid_size;
        int nghost;
        {
            ParmParse pp;

            Vector<int> n_cell_v;
            pp.getarr("n_cell", n_cell_v);
            n_cell = IntVect(AMREX_D_DECL(n_cell_v[0],n_cell_v[1],n_cell_v[2]));

            pp.get("max_grid_size", max_grid_size);
            pp.get("nghost", nghost);
        }

        Box domain(IntVect(0), n_cell-1);
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);
    
        MultiFab xmf(amrex::convert(ba,IntVect::TheDimensionVector(0)), dm, 1, nghost);
        MultiFab ymf(amrex::convert(ba,IntVect::TheDimensionVector(1)), dm, 1, nghost);
        MultiFab zmf(amrex::convert(ba,IntVect::TheDimensionVector(2)), dm, 1, nghost);
        MultiFab mf3(amrex::convert(ba,IntVect::TheDimensionVector(0)), dm, 3, nghost);

        xmf.setVal(0.0);
        ymf.setVal(1.0);
        zmf.setVal(2.0);
        mf3.setVal(3.0);

        Periodicity period(n_cell);

        xmf.FillBoundary(period);
        ymf.FillBoundary(period);
        zmf.FillBoundary(period);
        mf3.FillBoundary(period);

        ParallelDescriptor::Barrier();
        {
            BL_PROFILE("FB-vector");
            FillBoundary<MultiFab>({&xmf,&ymf,&zmf}, period);
        }

        ParallelDescriptor::Barrier();
        {
            BL_PROFILE("FB-111");
            xmf.FillBoundary(period);
            ymf.FillBoundary(period);
            zmf.FillBoundary(period);
        }

        ParallelDescriptor::Barrier();
        {
            BL_PROFILE("FB-3");
            mf3.FillBoundary(period);
        }
    }
    amrex::Finalize();
}
