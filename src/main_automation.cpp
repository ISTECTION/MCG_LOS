#include "utils/Automation.hpp"

int main(int argc, char* argv[]) {

    param.metod  = Method::MCG;
    param.matrix = Matrix::DIAGONAL_DOMINATION;
    param.cond   = Conditional::NONE;

    Automation _auto2(param);
    _auto2.start(10);

    return 0;
}
