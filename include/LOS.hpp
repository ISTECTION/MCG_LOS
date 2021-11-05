#ifndef _LOS_HPP_
#define _LOS_HPP_
#include "Data.hpp"
#include "utils/Function.hpp"

using namespace symmetric;

template <class T>
class LOS : public Data<T>
{
public:
     LOS(path _Path) : Data<T>(_Path) { }
    ~LOS() { }

    void solve(Conditional _cond, bool isLog = true);

private:
    void none    (bool);
    void diagonal(bool);
    void hollesky(bool);
};

template <class T>
void LOS<T>::solve(Conditional _cond, bool isLog) {
    switch (_cond) {
        case Conditional::NONE:     none    (isLog); break;
        case Conditional::DIAGONAL: diagonal(isLog); break;
        case Conditional::HOLLESKY: hollesky(isLog); break;
    }
}

template <class T>
void LOS<T>::none(bool isLog) {
    std::vector<T> r (this->param.n),
                   z (this->param.n),
                   p (this->param.n),
                   Ar(this->param.n);

    T alpha, betta, eps;

    r = this->pr - this->mult(this->x);
    z = r;
    p = this->mult(z);

    do {
        betta   = scalar(p, p);
        alpha   = scalar(p, r) / betta;
        this->x = this->x + alpha * z;
        r       = r - alpha * p;
        Ar      = this->mult(r);
        betta   = -scalar(p, Ar) / betta;
        z       = r + betta * z;
        p       = Ar + betta * p;
        eps     = scalar(r, r);

        this->iter++;
        if (isLog)
            printLog(this->iter, eps);

    } while (this->iter < this->param.max_iter
                && eps  > this->param.epsilon );
}

template <class T>
void LOS<T>::diagonal(bool isLog) {
    std::vector<T>  r(this->param.n),
                    z(this->param.n),
                    Ar(this->param.n),
                    p(this->param.n);

    std::vector<T> L(this->param.n, 1);
    for (size_t i = 0; i < L.size(); i++)
        L[i] /= sqrt(this->di[i]);

    T alpha, betta, eps;

    r = L * (this->pr - this->mult(this->x));
    z = L * r;
    p = L * this->mult(z);
    do {
        betta = scalar(p, p);
        alpha = scalar(p, r) / betta;
        this->x = this->x + alpha * z;
        r = r - alpha * p;
        betta = -scalar(p, L * this->mult(L * r)) / betta;
        z = L * r + betta * z;
        p = L * this->mult(L * r) + betta * p;
        eps = scalar(r, r);

        this->iter++;
        if (isLog)
            printLog(this->iter, eps);

    } while (this->iter < this->param.max_iter
                && eps  > this->param.epsilon );

}

template <class T>
void LOS<T>::hollesky(bool isLog) { }

#endif // _LOS_HPP_