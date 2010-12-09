/*
 * This file exists simply as a means to force CMake into appending link flags
 * by creating a phony target.
 */
void CMakeDummyFunction( int* a )
{
    *a = 42;  
}

