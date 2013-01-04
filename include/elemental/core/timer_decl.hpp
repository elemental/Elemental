/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

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

} // namespace elem
