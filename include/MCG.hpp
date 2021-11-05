#ifndef _MCG_HPP_
#define _MCG_HPP_
#include "Data.hpp"
#include "utils/Function.hpp"

using namespace symmetric;

template <class T>
class MCG : public Data<T>
{
public:
     MCG(path _Path) : Data<T>(_Path) { }
    ~MCG() { }

    void solve(Conditional _cond, bool isLog = true);
private:
    void none    (bool);
    void diagonal(bool);
    void hollesky(bool);
};

template <class T>
void MCG<T>::solve(Conditional _cond, bool isLog) {
    switch (_cond) {
        case Conditional::NONE:     none    (isLog); break;
        case Conditional::DIAGONAL: diagonal(isLog); break;
        case Conditional::HOLLESKY: hollesky(isLog); break;
    }
}

template <class T>
void MCG<T>::none(bool isLog) {
    std::vector<T> r (this->param.n),
                   z (this->param.n),
                   Az(this->param.n);

    T alpha, betta, eps;

    r = this->pr - this->mult(this->x);
    z = r;

    do {

        Az      = this->mult(z);
        betta   = scalar(r, r);
        alpha   = betta / scalar(Az, z);
        this->x = this->x + alpha * z;
        r       = r - alpha * Az;
        betta   = scalar(r, r) / betta;
        z       = r + betta * z;
        eps     = norm(r) / norm(this->pr);

        this->iter++;
        if (isLog)
            printLog(this->iter, eps);

    } while (this->iter < this->param.max_iter
                && eps  > this->param.epsilon );
}

template <class T>
void MCG<T>::diagonal(bool isLog) {
    std::vector<T> r(this->param.n),
                   z(this->param.n),
                   Az(this->param.n),
                   Mr(this->param.n);

    std::vector<T> M(this->param.n, 1);
    for (size_t i = 0; i < M.size(); i++)
        M[i] /= this->di[i];

    T alpha, betta, eps;

    r = this->pr - this->mult(this->x);
    z = M * r;

    do {
        Az      = this->mult(z);
        betta   = scalar(M * r, r);
        alpha   = betta / scalar(Az, z);
        this->x = this->x + alpha * z;
        r       = r - alpha * Az;
        Mr      = M * r;
        betta   = scalar(Mr, r) / betta;
        z       = Mr + betta * z;
        eps     = norm(r) / norm(this->pr);

        this->iter++;
        if (isLog)
            printLog(this->iter, eps);

    } while (this->iter < this->param.max_iter
                && eps  > this->param.epsilon );
}

template <class T>
void MCG<T>::hollesky(bool isLog) { }

#endif // _MCG_HPP_