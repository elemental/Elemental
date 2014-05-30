/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

Timer::Timer( const std::string& name )
: name_(name)
{ }

void Timer::Start()
{
    DEBUG_ONLY(
        if( running_ )
            LogicError("Forgot to stop timer before restarting.");
    )
#ifdef EL_HAVE_STEADYCLOCK
    lastTime_ = steady_clock::now();
#else
    lastTime_ = high_resolution_clock::now();
#endif
    running_ = true;
}

double Timer::Stop()
{
    DEBUG_ONLY(
        if( !running_ )
            LogicError("Tried to stop a timer before starting it.");
    )
    lastPartialTime_ = Partial();
    running_ = false;
    totalTime_ += lastPartialTime_;
    return lastPartialTime_;
}

void Timer::Reset( const std::string& name )
{ 
    name_ = name;
    running_ = false;
    totalTime_ = 0; 
    lastPartialTime_ = 0;
}

const std::string& Timer::Name() const { return name_; }

double Timer::Partial() const
{ 
    if( running_ )
    {
#ifdef EL_HAVE_STEADYCLOCK
        auto now = steady_clock::now();
#else
        auto now = high_resolution_clock::now();
#endif
        auto timeSpan = duration_cast<duration<double>>(now-lastTime_);
        return timeSpan.count();
    }
    else
        return lastPartialTime_; 
}

double Timer::Total() const
{
    if( running_ )
        return totalTime_ + Partial();
    else
        return totalTime_;
}

} // namespace El
