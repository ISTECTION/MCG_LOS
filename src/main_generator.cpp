#include "utils/generator.hpp"

int main(int argc, char* argv[]) {
    Generator gen("file/generator/random", 1000, 0.00000001, 100000, Chance { 1, 500 });


    // for (size_t i = 1; i < 13; i++)
    //     GenerateGilbert gen("file/generator/gilbert", i, 0.000000000001, 10000);

    return 0;
}