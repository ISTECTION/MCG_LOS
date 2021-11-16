#include <iostream>

#include "MCG.hpp"
#include "LOS.hpp"
#include "utils/Timer.hpp"

int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);

    // Timer::Timer timer;
    // MCG<double> a("file/generator/diagonal/10/p");
    // a.solve(Conditional::NONE, true);
    // timer.setTimeEnd();
    // std::cout << std::fixed << timer << std::endl;
    // a.printX();

    Timer::Timer timer;
    MCG<double> a("file/5");
    a.solve(Conditional::HOLLESKY, true);
    timer.setTimeEnd();
    std::cout << std::fixed << timer << std::endl;
    a.printX();


    return 0;
}
