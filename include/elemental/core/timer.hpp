/*
   Copyright (c) 2009-2011, Jack Poulson
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

namespace elemental {

class Timer
{
    bool _running;
    double _lastStartTime;
    double _time;
    const std::string _name;
public:
    Timer();
    Timer( const std::string name );

    void Start();
    void Stop();
    void Reset();

    const std::string Name() const;
    double Time() const;
};

} // elemental

// Implementations
inline elemental::Timer::Timer()
: _running(false), _time(0), _name("[blank]")
{ }

inline elemental::Timer::Timer( const std::string name )
: _running(false), _time(0), _name(name)
{ }

inline void
elemental::Timer::Start()
{
#ifndef RELEASE
    PushCallStack("Timer::Start");
    if( _running )
        throw std::logic_error("Forgot to stop timer before restarting.");
#endif
    _lastStartTime = mpi::Time();
    _running = true;
    _running = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::Timer::Stop()
{
#ifndef RELEASE
    PushCallStack("Timer::Stop");
    if( !_running )
        throw std::logic_error("Tried to stop a timer before starting it.");
#endif
    _time += mpi::Time()-_lastStartTime;
    _running = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::Timer::Reset()
{ _time = 0; }

inline const std::string
elemental::Timer::Name() const
{ return _name; }

inline double
elemental::Timer::Time() const
{
#ifndef RELEASE
    PushCallStack("Timer::Time");
    if( _running )
        throw std::logic_error("Asked for time while still timing.");
    PopCallStack();
#endif
    return _time;
}

#endif /* ELEMENTAL_TIMER_HPP */

