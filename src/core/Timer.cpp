/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

Timer::Timer( const string& name )
: name_(name)
{ }

void Timer::Start()
{
    DEBUG_ONLY(
      if( running_ )
          LogicError("Forgot to stop timer before restarting.");
    )
    lastTime_ = Clock::now();
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

void Timer::Reset( const string& name )
{ 
    name_ = name;
    running_ = false;
    totalTime_ = 0; 
    lastPartialTime_ = 0;
}

const string& Timer::Name() const { return name_; }

double Timer::Partial() const
{ 
    if( running_ )
    {
        auto now = Clock::now();
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
