#include "FileParser.h"

int main(int argc, char *argv[]) {
    if (argc < 1) {
        std::cerr << "Too few input arguments. Include input file." << std::endl;
        exit(1);
    }
    // Read input file
    FileParser f(argv[1]);
    f.readFile();
    f.print();
    // Create UnitChain
    // Perform Search

    return 0;
}