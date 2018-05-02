#include <openbabel/obconversion.h>
#include "FileParser.h"
#include "Containers.h"
#include "Chain.h"
#include "SearchTypes/MonteCarloRotorSearch.h"

using namespace PNAB;
using namespace OpenBabel;
using namespace std;

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Too few input arguments. Include input file." << std::endl;
        throw 1;
    }
    // Read input file
    FileParser f(argv[1]);
    f.readFile();
    f.print();
    RuntimeParameters rp(f);
    Backbone backbone(f);
    Bases bases(f);
    PNAB::RuntimeParameters runtime_params(f);
    HelicalParameters hp(f);
    MonteCarloRotorSearch mcrs(runtime_params, backbone, hp,  bases);
    mcrs.run();
    return 0;
}