/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef ELEMENTAL_TIMER_HPP
#define ELEMENTAL_TIMER_HPP 1

namespace elem {

class Timer
{
public:
    Timer();
    Timer( const std::string name );

    void Start();
    void Stop();
    void Reset();

    const std::string Name() const;
    double Time() const;
    
private:
    bool running_;
    double lastStartTime_;
    double time_;
    const std::string name_;
};

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline Timer::Timer()
: running_(false), time_(0), name_("[blank]")
{ }

inline Timer::Timer( const std::string name )
: running_(false), time_(0), name_(name)
{ }

inline void Timer::Start()
{
#ifndef RELEASE
    PushCallStack("Timer::Start");
    if( running_ )
        throw std::logic_error("Forgot to stop timer before restarting");
#endif
    lastStartTime_ = mpi::Time();
    running_ = true;
    running_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void Timer::Stop()
{
#ifndef RELEASE
    PushCallStack("Timer::Stop");
    if( !running_ )
        throw std::logic_error("Tried to stop a timer before starting it");
#endif
    time_ += mpi::Time()-lastStartTime_;
    running_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void Timer::Reset()
{ time_ = 0; }

inline const std::string Timer::Name() const
{ return name_; }

inline double Timer::Time() const
{
#ifndef RELEASE
    PushCallStack("Timer::Time");
    if( running_ )
        throw std::logic_error("Asked for time while still timing");
    PopCallStack();
#endif
    return time_;
}

} // namespace elem

#endif /* ELEMENTAL_TIMER_HPP */

