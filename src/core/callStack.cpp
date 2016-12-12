/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <stack>

namespace {

// Debugging
EL_DEBUG_ONLY(
  std::stack<std::string> callStack;
  bool tracingEnabled = false;
)

}

namespace El {

// If we are not in RELEASE mode, then implement wrappers for a call stack
EL_DEBUG_ONLY(

  void EnableTracing() { ::tracingEnabled = true; }
  void DisableTracing() { ::tracingEnabled = false; }

  void PushCallStack( string s )
  { 
      // [1]:
      // It was discovered that a global instantiation of a BigInt
      // (::bigIntZero) led to pushing to the call stack in global scope
      // before entering main, and this could possibly have been before the
      // call stack was constructed, leading to PopCallStack() causing an
      // exception due to the stack being empty. But this isn't completely
      // verified.
      if( !Initialized() )
          return;
#ifdef EL_HYBRID
      if( omp_get_thread_num() != 0 )
          return;
#endif
      const size_t maxStackSize = 300;
      if( ::callStack.size() > maxStackSize )
      {
          DumpCallStack();
          return;
      }
      ::callStack.push(s); 
      if( ::tracingEnabled )
      {
          const int stackSize = ::callStack.size();
          ostringstream os;
          for( int j=0; j<stackSize; ++j )
              os << " "; 
          os << s << endl;
          cout <<  os.str();
      }
  }

  void PopCallStack()
  { 
      // See note [1] above.
      if( !Initialized() )
          return;
#ifdef EL_HYBRID
      if( omp_get_thread_num() != 0 )
          return;
#endif
      if( ::callStack.empty() )
          LogicError("Attempted to pop an empty call stack");
      ::callStack.pop(); 
  }

  void DumpCallStack( ostream& os )
  {
      ostringstream msg;
      while( ! ::callStack.empty() )
      {
          msg << "[" << ::callStack.size() << "]: " << ::callStack.top() 
              << "\n";
          ::callStack.pop();
      }
      os << msg.str();
      os.flush();
  }

) // EL_DEBUG_ONLY

} // namespace El
