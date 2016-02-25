/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SERIALIZE_HPP
#define EL_SERIALIZE_HPP

namespace El {

#ifdef EL_HAVE_MPC

byte* Serialize( Int n, const BigInt* x, byte* xPacked );
byte* Serialize( Int n, const ValueInt<BigInt>* x, byte* xPacked );
byte* Serialize( Int n, const Entry<BigInt>* x, byte* xPacked );

byte* Serialize( Int n, const BigFloat* x, byte* xPacked );
byte* Serialize( Int n, const ValueInt<BigFloat>* x, byte* xPacked );
byte* Serialize( Int n, const Entry<BigFloat>* x, byte* xPacked );

byte* Deserialize( Int n, byte* xPacked, BigInt* x );
byte* Deserialize( Int n, byte* xPacked, ValueInt<BigInt>* x );
byte* Deserialize( Int n, byte* xPacked, Entry<BigInt>* x );

byte* Deserialize( Int n, byte* xPacked, BigFloat* x );
byte* Deserialize( Int n, byte* xPacked, ValueInt<BigFloat>* x );
byte* Deserialize( Int n, byte* xPacked, Entry<BigFloat>* x );

const byte* Deserialize( Int n, const byte* xPacked, BigInt* x );
const byte* Deserialize( Int n, const byte* xPacked, ValueInt<BigInt>* x );
const byte* Deserialize( Int n, const byte* xPacked, Entry<BigInt>* x );

const byte* Deserialize( Int n, const byte* xPacked, BigFloat* x );
const byte* Deserialize( Int n, const byte* xPacked, ValueInt<BigFloat>* x );
const byte* Deserialize( Int n, const byte* xPacked, Entry<BigFloat>* x );

void ReserveSerialized
( Int n, const BigInt* x, std::vector<byte>& xPacked );
void ReserveSerialized
( Int n, const ValueInt<BigInt>* x, std::vector<byte>& xPacked );
void ReserveSerialized
( Int n, const Entry<BigInt>* x, std::vector<byte>& xPacked );

void ReserveSerialized
( Int n, const BigFloat* x, std::vector<byte>& xPacked );
void ReserveSerialized
( Int n, const ValueInt<BigFloat>* x, std::vector<byte>& xPacked );
void ReserveSerialized
( Int n, const Entry<BigFloat>* x, std::vector<byte>& xPacked );

void Serialize
( Int n, const BigInt* x, std::vector<byte>& xPacked );
void Serialize
( Int n, const ValueInt<BigInt>* x, std::vector<byte>& xPacked );
void Serialize
( Int n, const Entry<BigInt>* x, std::vector<byte>& xPacked );

void Serialize
( Int n, const BigFloat* x, std::vector<byte>& xPacked );
void Serialize
( Int n, const ValueInt<BigFloat>* x, std::vector<byte>& xPacked );
void Serialize
( Int n, const Entry<BigFloat>* x, std::vector<byte>& xPacked );

void Deserialize
( Int n, const std::vector<byte>& xPacked, BigInt* x );
void Deserialize
( Int n, const std::vector<byte>& xPacked, ValueInt<BigInt>* x );
void Deserialize
( Int n, const std::vector<byte>& xPacked, Entry<BigInt>* x );

void Deserialize
( Int n, const std::vector<byte>& xPacked, BigFloat* x );
void Deserialize
( Int n, const std::vector<byte>& xPacked, ValueInt<BigFloat>* x );
void Deserialize
( Int n, const std::vector<byte>& xPacked, Entry<BigFloat>* x );

#endif // ifdef EL_HAVE_MPC

} // namespace El

#endif // ifndef EL_SERIALIZE_HPP
