#include <openbabel/obconversion.h>
#include "FileParser.h"
#include "Containers.h"
#include "UnitChain.h"
#include "SearchTypes/RandomRotorSearch.h"
#include "Chain.h"
#include "SearchTypes/MonteCarloRotorSearch.h"

using namespace PNAB;
using namespace OpenBabel;

int main(int argc, char *argv[]) {
    if (argc < 1) {
        std::cerr << "Too few input arguments. Include input file." << std::endl;
        exit(1);
    }
    // Read input file
    FileParser f(argv[1]);
    f.readFile();
    f.print();
    RuntimeParameters rp(f);
    Backbone backbone(f);
    Bases bases(f);
    //UnitChain uc({"ADA", "ADA"},bases, backbone);
    OBConversion conv;
    std::filebuf fb;
    fb.open("out.pdb", std::ios::out);
    std::ostream fileStream(&fb);
    conv.SetOutFormat("PDB");
    conv.SetOutStream(&fileStream);
    OBMol mol;
    std::array<double, 3> translate_arr = {0, 0, 3.4},
                          rotation_arr  = {0, 0,  30};
    //uc.updateHelicalParameters(translate_arr, rotation_arr);
    //mol = uc.getUnit(0,0).getMol();
    PNAB::RuntimeParameters runtime_params(f);
    //RandomRotorSearch rrs(runtime_params,uc);
    //rrs.speak();
    //mol = rrs.search();

    HelicalParameters hp(f);
    MonteCarloRotorSearch mcrs(runtime_params, bases.bases[0], bases.bases[0], backbone, hp);
    mcrs.run();
    auto mol_w_conf = mcrs.getMolecule();
    auto coord_vec = mcrs.getCoordinates();
    if (!coord_vec.empty())
        mol_w_conf.SetCoordinates(coord_vec[0]);

//    Chain c(BaseUnit(PNAB::Base(), PNAB::Backbone()));
//    mol = c.getChain();
//    conv.Write(&mol);

    std::filebuf fb_1;
    fb_1.open("out_mc.pdb", std::ios::out);
    std::ostream fileStream_1(&fb_1);
    conv.SetOutFormat("PDB");
    conv.SetOutStream(&fileStream_1);

    conv.Write(&mol_w_conf);

    // Create UnitChain
    // Perform Search

    return 0;
}