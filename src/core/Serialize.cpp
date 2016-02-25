/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#ifdef EL_HAVE_MPC

namespace El {

byte* Serialize( Int n, const BigInt* x, byte* buf )
{
    DEBUG_ONLY(CSE cse("Serialize"))
    for( Int j=0; j<n; ++j )
        buf = x[j].Serialize( buf );
    return buf;
}

byte* Serialize( Int n, const BigFloat* x, byte* buf )
{
    DEBUG_ONLY(CSE cse("Serialize"))
    for( Int j=0; j<n; ++j )
        buf = x[j].Serialize( buf );
    return buf;
}

byte* Serialize( Int n, const ValueInt<BigInt>* x, byte* buf )
{
    DEBUG_ONLY(CSE cse("Serialize"))
    for( Int j=0; j<n; ++j )
    {
        buf = x[j].value.Serialize( buf );
        std::memcpy( buf, &x[j].index, sizeof(Int) );
        buf += sizeof(Int);
    }
    return buf;
}

byte* Serialize( Int n, const ValueInt<BigFloat>* x, byte* buf )
{
    DEBUG_ONLY(CSE cse("Serialize"))
    for( Int j=0; j<n; ++j )
    {
        buf = x[j].value.Serialize( buf );
        std::memcpy( buf, &x[j].index, sizeof(Int) );
        buf += sizeof(Int);
    }
    return buf;
}

byte* Serialize( Int n, const Entry<BigInt>* x, byte* buf )
{
    DEBUG_ONLY(CSE cse("Serialize"))
    for( Int k=0; k<n; ++k )
    {
        std::memcpy( buf, &x[k].i, sizeof(Int) );
        buf += sizeof(Int);
        std::memcpy( buf, &x[k].j, sizeof(Int) );
        buf += sizeof(Int);
        buf = x[k].value.Serialize( buf );
    }
    return buf;
}

byte* Serialize( Int n, const Entry<BigFloat>* x, byte* buf )
{
    DEBUG_ONLY(CSE cse("Serialize"))
    for( Int k=0; k<n; ++k )
    {
        std::memcpy( buf, &x[k].i, sizeof(Int) );
        buf += sizeof(Int);
        std::memcpy( buf, &x[k].j, sizeof(Int) );
        buf += sizeof(Int);
        buf = x[k].value.Serialize( buf );
    }
    return buf;
}

void ReserveSerialized( Int n, const BigInt* x, std::vector<byte>& buf )
{
    DEBUG_ONLY(CSE cse("ReserveSerialized"))
    if( n == 0 )
    {
        buf.resize(0);
        return;
    }
    const auto packedSize = x[0].SerializedSize();
    buf.resize( n*packedSize );
}

void ReserveSerialized( Int n, const BigFloat* x, std::vector<byte>& buf )
{
    DEBUG_ONLY(CSE cse("ReserveSerialized"))
    if( n == 0 )
    {
        buf.resize(0);
        return;
    }
    const auto packedSize = x[0].SerializedSize();
    buf.resize( n*packedSize );
}

void ReserveSerialized
( Int n, const ValueInt<BigInt>* x, std::vector<byte>& buf )
{
    DEBUG_ONLY(CSE cse("ReserveSerialized"))
    if( n == 0 )
    {
        buf.resize(0);
        return;
    }
    const auto packedSize = x[0].value.SerializedSize() + sizeof(Int);
    buf.resize( n*packedSize );
}

void ReserveSerialized
( Int n, const ValueInt<BigFloat>* x, std::vector<byte>& buf )
{
    DEBUG_ONLY(CSE cse("ReserveSerialized"))
    if( n == 0 )
    {
        buf.resize(0);
        return;
    }
    const auto packedSize = x[0].value.SerializedSize() + sizeof(Int);
    buf.resize( n*packedSize );
}

void ReserveSerialized
( Int n, const Entry<BigInt>* x, std::vector<byte>& buf )
{
    DEBUG_ONLY(CSE cse("ReserveSerialized"))
    if( n == 0 )
    {
        buf.resize(0);
        return;
    }
    const auto packedSize = x[0].value.SerializedSize() + 2*sizeof(Int);
    buf.resize( n*packedSize );
}

void ReserveSerialized
( Int n, const Entry<BigFloat>* x, std::vector<byte>& buf )
{
    DEBUG_ONLY(CSE cse("ReserveSerialized"))
    if( n == 0 )
    {
        buf.resize(0);
        return;
    }
    const auto packedSize = x[0].value.SerializedSize() + 2*sizeof(Int);
    buf.resize( n*packedSize );
}

void Serialize( Int n, const BigInt* x, std::vector<byte>& buf )
{
    DEBUG_ONLY(CSE cse("Serialize"))
    ReserveSerialized( n, x, buf );
    Serialize( n, x, buf.data() );
}

void Serialize( Int n, const BigFloat* x, std::vector<byte>& buf )
{
    DEBUG_ONLY(CSE cse("Serialize"))
    ReserveSerialized( n, x, buf );
    Serialize( n, x, buf.data() );
}

void Serialize( Int n, const ValueInt<BigInt>* x, std::vector<byte>& buf )
{
    DEBUG_ONLY(CSE cse("Serialize"))
    ReserveSerialized( n, x, buf );
    Serialize( n, x, buf.data() );
}

void Serialize( Int n, const ValueInt<BigFloat>* x, std::vector<byte>& buf )
{
    DEBUG_ONLY(CSE cse("Serialize"))
    ReserveSerialized( n, x, buf );
    Serialize( n, x, buf.data() );
}

void Serialize( Int n, const Entry<BigInt>* x, std::vector<byte>& buf )
{
    DEBUG_ONLY(CSE cse("Serialize"))
    ReserveSerialized( n, x, buf );
    Serialize( n, x, buf.data() );
}

void Serialize( Int n, const Entry<BigFloat>* x, std::vector<byte>& buf )
{
    DEBUG_ONLY(CSE cse("Serialize"))
    ReserveSerialized( n, x, buf );
    Serialize( n, x, buf.data() );
}

byte* Deserialize( Int n, byte* buf, BigInt* x )
{
    DEBUG_ONLY(CSE cse("Deserialize [BigInt]"))
    for( Int j=0; j<n; ++j )
        buf = x[j].Deserialize( buf );
    return buf;
}

byte* Deserialize( Int n, byte* buf, BigFloat* x )
{
    DEBUG_ONLY(CSE cse("Deserialize [BigFloat]"))
    for( Int j=0; j<n; ++j )
        buf = x[j].Deserialize( buf );
    return buf;
}

byte* Deserialize( Int n, byte* buf, ValueInt<BigInt>* x )
{
    DEBUG_ONLY(CSE cse("Deserialize [ValueInt<BigInt>]"))
    for( Int j=0; j<n; ++j )
    {
        buf = x[j].value.Deserialize( buf );
        std::memcpy( &x[j].index, buf, sizeof(Int) );
        buf += sizeof(Int);
    }
    return buf;
}

byte* Deserialize( Int n, byte* buf, ValueInt<BigFloat>* x )
{
    DEBUG_ONLY(CSE cse("Deserialize [ValueInt<BigFloat>]"))
    for( Int j=0; j<n; ++j )
    {
        buf = x[j].value.Deserialize( buf );
        std::memcpy( &x[j].index, buf, sizeof(Int) );
        buf += sizeof(Int);
    }
    return buf;
}

byte* Deserialize( Int n, byte* buf, Entry<BigInt>* x )
{
    DEBUG_ONLY(CSE cse("Deserialize [Entry<BigInt>]"))
    for( Int k=0; k<n; ++k )
    {
        std::memcpy( &x[k].i, buf, sizeof(Int) );
        buf += sizeof(Int);
        std::memcpy( &x[k].j, buf, sizeof(Int) );
        buf += sizeof(Int);
        buf = x[k].value.Deserialize( buf );
    }
    return buf;
}

byte* Deserialize( Int n, byte* buf, Entry<BigFloat>* x )
{
    DEBUG_ONLY(CSE cse("Deserialize [Entry<BigFloat>]"))
    for( Int k=0; k<n; ++k )
    {
        std::memcpy( &x[k].i, buf, sizeof(Int) );
        buf += sizeof(Int);
        std::memcpy( &x[k].j, buf, sizeof(Int) );
        buf += sizeof(Int);
        buf = x[k].value.Deserialize( buf );
    }
    return buf;
}

const byte* Deserialize( Int n, const byte* buf, BigInt* x )
{
    DEBUG_ONLY(CSE cse("Deserialize [BigInt]"))
    for( Int j=0; j<n; ++j )
        buf = x[j].Deserialize( buf );
    return buf;
}

const byte* Deserialize( Int n, const byte* buf, BigFloat* x )
{
    DEBUG_ONLY(CSE cse("Deserialize [BigFloat]"))
    for( Int j=0; j<n; ++j )
        buf = x[j].Deserialize( buf );
    return buf;
}

const byte* Deserialize( Int n, const byte* buf, ValueInt<BigInt>* x )
{
    DEBUG_ONLY(CSE cse("Deserialize [ValueInt<BigInt>]"))
    for( Int j=0; j<n; ++j )
    {
        buf = x[j].value.Deserialize( buf );
        std::memcpy( &x[j].index, buf, sizeof(Int) );
        buf += sizeof(Int);
    }
    return buf;
}

const byte* Deserialize( Int n, const byte* buf, ValueInt<BigFloat>* x )
{
    DEBUG_ONLY(CSE cse("Deserialize [ValueInt<BigFloat>]"))
    for( Int j=0; j<n; ++j )
    {
        buf = x[j].value.Deserialize( buf );
        std::memcpy( &x[j].index, buf, sizeof(Int) );
        buf += sizeof(Int);
    }
    return buf;
}

const byte* Deserialize( Int n, const byte* buf, Entry<BigInt>* x )
{
    DEBUG_ONLY(CSE cse("Deserialize [Entry<BigInt>]"))
    for( Int k=0; k<n; ++k )
    {
        std::memcpy( &x[k].i, buf, sizeof(Int) );
        buf += sizeof(Int);
        std::memcpy( &x[k].j, buf, sizeof(Int) );
        buf += sizeof(Int);
        buf = x[k].value.Deserialize( buf );
    }
    return buf;
}

const byte* Deserialize( Int n, const byte* buf, Entry<BigFloat>* x )
{
    DEBUG_ONLY(CSE cse("Deserialize [Entry<BigFloat>]"))
    for( Int k=0; k<n; ++k )
    {
        std::memcpy( &x[k].i, buf, sizeof(Int) );
        buf += sizeof(Int);
        std::memcpy( &x[k].j, buf, sizeof(Int) );
        buf += sizeof(Int);
        buf = x[k].value.Deserialize( buf );
    }
    return buf;
}

void Deserialize( Int n, const std::vector<byte>& buf, BigInt* x )
{
    Deserialize( n, buf.data(), x ); 
}
void Deserialize( Int n, const std::vector<byte>& buf, ValueInt<BigInt>* x )
{
    Deserialize( n, buf.data(), x );
}
void Deserialize( Int n, const std::vector<byte>& buf, Entry<BigInt>* x )
{
    Deserialize( n, buf.data(), x );
}

void Deserialize( Int n, const std::vector<byte>& buf, BigFloat* x )
{
    Deserialize( n, buf.data(), x ); 
}
void Deserialize( Int n, const std::vector<byte>& buf, ValueInt<BigFloat>* x )
{
    Deserialize( n, buf.data(), x );
}
void Deserialize( Int n, const std::vector<byte>& buf, Entry<BigFloat>* x )
{
    Deserialize( n, buf.data(), x );
}

} // namespace El

#endif // ifdef EL_HAVE_MPC
