#include "utils/Automation.hpp"

int main(int argc, char* argv[]) {

    param.metod  = Method::LOS;
    param.matrix = Matrix::DIAGONAL_DOMINATION;
    param.cond   = Conditional::HOLLESKY;

    Automation _auto2(param);
    _auto2.start(10);

    return 0;
}
