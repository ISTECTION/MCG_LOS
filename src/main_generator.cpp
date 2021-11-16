#include "utils/generator.hpp"

int main(int argc, char* argv[]) {
    // Generator gen("file/generator/random", 10, 0.00000001, 100000, Chance { 1, 100 });

    for (size_t i = 1; i < 21; i++)
        GenerateGilbert gen("file/generator/gilbert", i, 1e-14, 100000);

    // Generate_Ak gen("file/generator/Ak", 10, 1e-14, 100000, 20);

    // Generate_diagDomination gen("file/generator/diagonal", 10, 1e-14, 100000);

    return 0;
}
