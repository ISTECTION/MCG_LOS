#ifndef _AUTOMATION_HPP_
#define _AUTOMATION_HPP_
#include "MCG.hpp"
#include "LOS.hpp"
#include "utils/Timer.hpp"

#include <locale>

using namespace symmetric;

using std::filesystem::is_directory;

class Punct : public std::numpunct<char>
{
protected:
    char do_decimal_point() const { return ','; }
};

enum class Method {
    MCG,
    LOS
};

enum class Matrix {
    DIAGONAL_DOMINATION,
    GILBERT,
    AK
};

struct Params {
    Method metod;
    Matrix matrix;
    Conditional cond;
} param { Method::MCG, Matrix::GILBERT, Conditional::NONE };

class Automation
{
private:
    Method method;
    Matrix matrix;
    Conditional cond;

    void writeFile(path _path, std::ostringstream &ostream) const;
    void step(path _path, std::ostringstream &ostream, size_t _size);

public:
    Automation(Params _param = param) :
        method(_param.metod ),
        matrix(_param.matrix),
        cond  (_param.cond  ) { }
    ~Automation() { }

    void start(size_t _size = 0);
};

void Automation::start(size_t _size) {

    std::filesystem::path _path;

    switch(matrix)
    {
        case Matrix::DIAGONAL_DOMINATION:
            _path = "file/generator/diagonal/" + std::to_string(_size);
            break;
        case Matrix::GILBERT:
            _path = "file/generator/gilbert";
            break;
        case Matrix::AK:
            _path = "file/generator/Ak/" + std::to_string(_size);
            break;
    }

    std::ostringstream ostream;
    ostream.imbue(std::locale(std::locale(), new Punct));

    ostream << std::setw(5)  << "Size:"
            << std::setw(10) << "Iter:"
            << std::setw(12) << "Time:"
            << std::setw(30) << "x:"
            << std::setw(18) << "x` - x:"
            << '\n';

    switch(matrix)
    {
        case Matrix::DIAGONAL_DOMINATION:
            ostream << '\n' << "negativ" << '\n';
            step(_path / "n", ostream, 1);
            ostream << "positiv" << '\n';
            step(_path / "p", ostream, 2);
            break;

        case Matrix::GILBERT:
        case Matrix::AK:
            for (size_t _size = 1; is_directory(_path / std::to_string(_size)); _size++)
                step(_path / std::to_string(_size), ostream, _size);
            break;
    }

    writeFile("file/automation", ostream);
}

void Automation::step(path _path, std::ostringstream &ostream, size_t _size) {

    std::vector<double> x, _expected;
    Timer::Timer timer; size_t _iter;

    switch (method)
    {
        case Method::MCG:
        {
            timer.setTimeStart();
            MCG<double> a(_path);
            a.solve(cond, false);
            timer.setTimeEnd();

            _iter = a.getIteration();
            x     = a.getX();
            break;
        }
        case Method::LOS:
        {
            timer.setTimeStart();
            LOS<double> a(_path);
            a.solve(cond, false);
            timer.setTimeEnd();

            _iter = a.getIteration();
            x     = a.getX();

            break;
        }
    }

    _expected.resize(x.size());
    std::iota(_expected.begin(), _expected.end(), 1);

    ostream << std::setw(5)  << std::fixed << _size
            << std::setw(10) << std::fixed << _iter;

    ostream.precision(6);
    ostream << std::setw(12) << timer.getElapsed();

    for (size_t i = 0; i < x.size(); i++)
    {
        ostream.precision(14);
        ostream << (i == 0 ? std::setw(30) : std::setw(57))
                << std::fixed << x[i];

        ostream.precision(3);
        ostream << std::setw(18) << std::scientific
                << _expected[i] - x[i] << '\n';
    }
    ostream << std::setw(76) << '\n';
}

void Automation::writeFile(path _path, std::ostringstream &ostream) const {

    switch(method)
    {
        case Method::MCG:
            _path /= "MCG";
            break;
        case Method::LOS:
            _path /= "LOS";
            break;
    }

    switch(cond)
    {
        case symmetric::Conditional::NONE:
            _path /= "none";
            break;
        case symmetric::Conditional::DIAGONAL:
            _path /= "diag";
            break;
        case symmetric::Conditional::HOLLESKY:
            _path /= "hollesky";
            break;
    }

    std::filesystem::create_directories(_path);
    switch(matrix)
    {
        case Matrix::DIAGONAL_DOMINATION:
            _path /= "diagonal.txt";
            break;
        case Matrix::GILBERT:
            _path /= "gilbert.txt";
            break;
        case Matrix::AK:
            _path /= "ak.txt";
            break;
    }

    std::ofstream fout(_path);
    fout << ostream.str();
    fout.close();
}
#endif // _AUTOMATION_HPP_