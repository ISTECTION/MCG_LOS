#include "utils/Automation.hpp"

int main(int argc, char* argv[]) {

    Automation _auto(Method::LOS, Matrix::GILBERT, symmetric::Conditional::NONE);
    _auto.start(10);

    return 0;
}
