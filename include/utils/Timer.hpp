#ifndef _TIMER_HPP_
#define _TIMER_HPP_
#include <iostream>
#include <chrono>

namespace Timer {
    using std::chrono::high_resolution_clock;
    using std::chrono::time_point;
    using std::chrono::duration;

    class Timer
    {
    private:
        time_point<high_resolution_clock> timeStart;
        time_point<high_resolution_clock> timeEnd;
    public:
        Timer() : timeStart(high_resolution_clock::now()){ }
        ~Timer() { };

        void setTimeEnd  ();
        void setTimeStart();
        float getElapsed () const;

        friend std::ostream& operator<< (std::ostream& out, Timer& point) {
            return out << "SECONDS: " << point.getElapsed() << std::endl;
        }
    };
    void  Timer::setTimeEnd()       { timeEnd = high_resolution_clock::now();   }
    void  Timer::setTimeStart()     { timeStart = high_resolution_clock::now(); };
    float Timer::getElapsed() const { return duration<float>(timeEnd - timeStart).count(); }
}
#endif // _TIMER_HPP_