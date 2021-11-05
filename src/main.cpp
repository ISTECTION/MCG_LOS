#include <iostream>

#include "MCG.hpp"
#include "LOS.hpp"
#include "utils/Timer.hpp"

int main(int argc, char* argv[]) {

    Timer::Timer timer;
    MCG<double> a("file/generator/Ak/10/10");
    a.solve(Conditional::NONE, false);
    timer.setTimeEnd();
    std::cout << timer << std::endl;

    a.printX();

    return 0;
}
