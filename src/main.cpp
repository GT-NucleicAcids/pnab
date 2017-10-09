#include <openbabel/obconversion.h>
#include "FileParser.h"
#include "Containers.h"
#include "UnitChain.h"

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
    UnitChain uc({"ada,ada,ada"},bases, backbone);
    OBConversion conv;
    std::filebuf fb;
    fb.open("out.cml", std::ios::out);
    std::ostream fileStream(&fb);
    conv.SetOutFormat("CML");
    conv.SetOutStream(&fileStream);
    OBMol mol;
    mol = *uc.getUnit(0,0).getMol();
    conv.Write(&mol);

    // Create UnitChain
    // Perform Search

    return 0;
}