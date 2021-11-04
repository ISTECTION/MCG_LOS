#ifndef _GENERATOR_HPP_
#define _GENERATOR_HPP_
#include "../Data.hpp"
#include <filesystem>
#include <numeric>
#include <random>

const int MIN = 1;
const int MAX = 9;

using std::filesystem::path;

int getRandom(int _minValue, int _maxValue) {
    std::random_device rd;
    std::mt19937 engine(rd());

    return (engine() % (_maxValue - _minValue) + _minValue);
}
bool getBool() { return getRandom(0, 2); }


struct Chance {
    uint32_t here;
    uint32_t from;
};

struct Vector {
    std::vector<double> di, gg, pr;
    std::vector<size_t> ig, jg;
};

void writeFile(path, Vector, symmetric::Param, uint8_t _comma = 0);
void write(path, std::vector<double> &, uint8_t);
void write(path, std::vector<size_t> &, uint8_t);


/// Класс генерации положительно определённых матриц
class Generator
{
private:
    symmetric::Param _param;
    Chance _chance;
    Vector _v;
    void generateFile(path _path);

public:
    Generator(path _path, size_t n, double eps, size_t max_iter, Chance chance = { 5, 100 })
        : _param  { n, eps, max_iter },
          _chance { chance.here, chance.from } { generateFile(_path); }
    ~Generator() { }
};

void Generator::generateFile(path _path) {

    _v.di.resize(_param.n);
    _v.pr.resize(_param.n);
    _v.ig.resize(_param.n + 1);

    _v.ig[0] = 0;

    size_t _Count = 0;
    for (size_t i = 0; i < _param.n; i++) {

        _v.di[i] = getRandom(MIN, MAX);

        for(size_t j = 0; j < i; j++) {

            if (getRandom(0, _chance.from) < _chance.here) {
                _Count++;
                _v.gg.push_back(getRandom(MIN, MAX));
                _v.jg.push_back(j);
            }
        }
        _v.ig[i + 1] = _Count;
    }

    std::vector<std::vector<uint32_t>> _none_Posi, _posi;
    _none_Posi.resize(_param.n, std::vector<uint32_t>(_param.n));
    _posi.resize(_param.n, std::vector<uint32_t>(_param.n));

    size_t jj = 0;
    for (size_t i = 0; i < _param.n; i++) {
        _none_Posi[i][i] = _v.di[i];

        for (size_t j = _v.ig[i]; j < _v.ig[i + 1]; j++, jj++) {
            _none_Posi[i][_v.jg[jj]] = _v.gg[jj];
            _none_Posi[_v.jg[jj]][i] = _v.gg[jj];
        }
    }

    for (size_t i = 0; i < _param.n; i++)
        for (size_t j = 0; j < _param.n; j++)
        {
            _posi[i][j] = 0;
            for (size_t k = 0; k < _param.n; k++)
                _posi[i][j] += _none_Posi[i][k] * _none_Posi[k][j];
        }

    std::vector<uint32_t> x(_param.n);
    std::iota(x.begin(), x.end(), 1);

    for (size_t i = 0; i < _param.n; i++) {
        uint32_t _sum = 0;
        for (size_t j = 0; j < _param.n; j++)
            _sum += _posi[i][j] * x[j];

        _v.pr[i] = _sum;
    }

    _v.gg.clear();
    _v.jg.clear();

    _Count = 0;
    for (size_t i = 0; i < _param.n; i++) {

        _v.di[i] = _posi[i][i];

        for(size_t j = 0; j < i; j++) {

            if (_posi[i][j] != 0) {
                _Count++;
                _v.gg.push_back(_posi[i][j]);
                _v.jg.push_back(j);
            }
        }
        _v.ig[i + 1] = _Count;
    }

    _path /= std::to_string(_param.n);
    _path /= "gg-" + std::to_string(_Count);
    std::filesystem::create_directories(_path);

    writeFile(_path, _v, _param);
}



class GenerateGilbert
{
private:
    symmetric::Param _param;
    Vector _v;

    void generateFile(path _path);

public:
    GenerateGilbert(path _path, size_t n, double eps, size_t max_iter)
        : _param { n, eps, max_iter } { generateFile(_path); };
    ~GenerateGilbert() { };

};

void GenerateGilbert::generateFile(path _path) {

    _path /= std::to_string(_param.n);
    std::filesystem::create_directories(_path);

    _v.di.resize(_param.n);
    _v.pr.resize(_param.n);
    _v.ig.resize(_param.n + 1);

    _v.ig[0] = 0;

    size_t _count_L = _param.n * (_param.n - 1) / 2;
    _v.gg.resize(_count_L);
    _v.jg.resize(_count_L);

    size_t _Count = 0;
    for (size_t i = 0; i < _param.n; i++) {

        _v.di[i] = 1. / ( (i + 1) * 2 - 1);

        for(size_t j = 0; j < i; j++) {

            _v.gg[_Count] = 1. / ( (i + 1) + (j + 1) - 1 );
            _v.jg[_Count] = j;

            _Count++;
        }
        _v.ig[i + 1] = _Count;
    }

    std::vector<uint32_t> x(_param.n);
    std::iota(x.begin(), x.end(), 1);

    for (size_t i = 0, jj = 0; i < x.size(); i++) {
        _v.pr[i] = _v.di[i] * x[i];

        for (size_t j = _v.ig[i]; j < _v.ig[i + 1]; j++, jj++) {
            _v.pr[i] += _v.gg[jj] * x[_v.jg[jj]];
            _v.pr[_v.jg[jj]] += _v.gg[jj] * x[i];
        }
    }

    writeFile(_path, _v, _param, 14);
}


void writeFile(path _path, Vector _v, symmetric::Param _param, uint8_t _comma) {
    std::ofstream fout(_path / "kuslau.txt");
    fout << _param.n
         << ' ' << _param.max_iter
         << ' ' << _param.epsilon;
    fout.close();

    write(_path / "di.txt", _v.di, _comma);
    write(_path / "gg.txt", _v.gg, _comma);
    write(_path / "ig.txt", _v.ig, _comma);
    write(_path / "jg.txt", _v.jg, _comma);
    write(_path / "pr.txt", _v.pr, _comma);
}

void write(path _path, std::vector<double> &_v, uint8_t _comma) {
    std::ofstream fout(_path);
    if(_comma) fout.precision(_comma);
    for (size_t i = 0; i < _v.size(); i++)
        fout << _v[i] << '\n';
    fout.close();
}

void write(path _path, std::vector<size_t> &_v, uint8_t _comma) {
    std::ofstream fout(_path);
    if(_comma) fout.precision(_comma);
    for (size_t i = 0; i < _v.size(); i++)
        fout << _v[i] << '\n';
    fout.close();
}
#endif // _GENERATOR_HPP_