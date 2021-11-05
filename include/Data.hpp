#ifndef _DATA_HPP_
#define _DATA_HPP_
#include "utils/namespace.hpp"

#include <filesystem>
#include <iostream>
#include <cassert>
#include <fstream>
#include <sstream>
#include <vector>

_SYMMETRIC_BEG                  /// Пространство имён симметричной матрицы

using std::filesystem::path;    /// Тип path используется для хранения путей

struct Param {
    size_t n;                   /// n        - Размерность матрицы
    double epsilon;             /// epsilon  - Точность решения СЛАУ
    size_t max_iter;            /// max_iter - MAX количество итераций
};

enum class Conditional
{
    NONE,                       /// Без предобусловливания
    DIAGONAL,                   /// С диагональным предобусловливанием
    HOLLESKY                    /// С предобусловливанием Холесского
};

template <class T>
class Data
{
protected:
    Param param;
    std::vector<T> di;          /// Диагональные элементы
    std::vector<T> gg;          /// Внедиагональные элементы
    std::vector<size_t> ig;     /// Указатели начала строк
    std::vector<size_t> jg;     /// Номера столбцов
    std::vector<T> pr;          /// Вектор правой части
    std::vector<T> x;           /// Вектор решения
    size_t iter = 0;            /// Количество итераций


public:
    Data(path _path = "file/default") { assert(loadData(_path)); }
    ~Data() { }

    void reset();

    std::vector<T>& getX() { return x; }
    size_t getIteration() const { return iter; }

    std::vector<T> mult(std::vector<T> _V);
    void printX(unsigned int count = 0) const;
private:
    bool loadData(path _Path);
    void resize(size_t _Mem);
    bool read(path _Path, std::vector<T>& _Vec);
    bool read(path _Path, std::vector<size_t>& _Vec);
};

template <class T>
void Data<T>::reset() {
    std::fill(x.begin(), x.end(), 0);
    iter = 0;
}

template <class T>
std::vector<T> Data<T>::mult(std::vector<T> _Vec) {
    std::vector<T> pr(_Vec.size());

    int jj = 0;
    for (size_t i = 0; i < _Vec.size(); i++) {
        pr[i] = di[i] * _Vec[i];

        for (size_t j = ig[i]; j < ig[i + 1]; j++, jj++) {
            pr[i] += gg[jj] * _Vec[jg[jj]];
            pr[jg[jj]] += gg[jj] * _Vec[i];
        }
    }
    return pr;
}

template <class T>
void Data<T>::printX(unsigned int count) const {
    std::ostringstream ostream;
    if (count)
        ostream.precision(count);

    ostream << "[ ";
    for (size_t i = 0; i < x.size(); i++)
        ostream << x[i] << " ";
    ostream << "]";
    std::cout << ostream.str();
}

template <class T>
bool Data<T>::loadData(path _Path) {
    bool isCor { true };

    std::ifstream fin(_Path / "kuslau.txt");
    if(!fin) return false;
        fin >> param.n
            >> param.max_iter
            >> param.epsilon;
    fin.close();

    ig.resize(param.n + 1);
    isCor &= read(_Path / "ig.txt", ig);

    resize(ig.back());

    isCor &= read(_Path / "gg.txt", gg);
    isCor &= read(_Path / "jg.txt", jg);
    isCor &= read(_Path / "di.txt", di);
    isCor &= read(_Path / "pr.txt", pr);

    return isCor;
}

template <class T>
bool Data<T>::read(path _Path, std::vector<T>& _Vec) {
    std::ifstream fin(_Path);
    if(!fin) return false;
    for (size_t i = 0; i < _Vec.size(); i++)
            fin >> _Vec[i];
    fin.close();
    return true;
}

template <class T>
bool Data<T>::read(path _Path, std::vector<size_t>& _Vec) {
    std::ifstream fin(_Path);
    if(!fin) return false;
    for (size_t i = 0; i < _Vec.size(); i++)
            fin >> _Vec[i];
    fin.close();
    return true;
}

template <class T>
void Data<T>::resize(size_t _Mem) {
    gg.resize(_Mem);
    jg.resize(_Mem);

    di.resize(param.n);
    pr.resize(param.n);
    x.resize (param.n);
}

_SYMMETRIC_END
#endif // _DATA_HPP_