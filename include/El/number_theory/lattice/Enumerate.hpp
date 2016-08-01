/*
   Copyright (c) 2009-2016, Jack Poulson,
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LATTICE_ENUMERATE_HPP
#define EL_LATTICE_ENUMERATE_HPP

namespace El {

enum EnumType {
  FULL_ENUM,
  GNR_ENUM,
  YSPARSE_ENUM
};

template<typename Real>
struct EnumCtrl
{
    EnumType enumType=FULL_ENUM;

    bool disablePrecDrop=false;
    Real fudge=Real(1.5); // fudge factor for number of bits of precision

    bool time=false;
    bool progress=false;

    // For monitoring the core (bounded) enumeration procedure
    bool innerProgress=false;

    // Explicitly transpose 'N' to encourage unit-stride access
    bool explicitTranspose=true;

    // GNR_ENUM
    // --------
    // TODO: Add ability to further tune the bounding function
    bool linearBounding=false;
    Int numTrials=1000;

    // YSPARSE_ENUM
    // ------------
    Int phaseLength=10;
    double enqueueProb=1.;

    bool customStartIndex=false;
    Int startIndex;

    bool customMinInfNorms=false;
    vector<Int> minInfNorms;

    bool customMaxInfNorms=false;
    vector<Int> maxInfNorms;

    bool customMinOneNorms=false;
    vector<Int> minOneNorms;

    bool customMaxOneNorms=false;
    vector<Int> maxOneNorms;

    Int progressLevel=0;

    template<typename OtherReal>
    EnumCtrl<Real>& operator=( const EnumCtrl<OtherReal>& ctrl )
    {
        enumType = ctrl.enumType;
        disablePrecDrop = ctrl.disablePrecDrop;
        fudge = Real(ctrl.fudge);
        time = ctrl.time;
        progress = ctrl.progress;
        innerProgress = ctrl.innerProgress;
        explicitTranspose = ctrl.explicitTranspose;

        // GNR_ENUM
        // --------
        linearBounding = ctrl.linearBounding;
        numTrials = ctrl.numTrials;

        // YSPARSE_ENUM
        // ------------
        phaseLength = ctrl.phaseLength;
        enqueueProb = ctrl.enqueueProb;
        customStartIndex = ctrl.customStartIndex;
        startIndex = ctrl.startIndex;

        customMinInfNorms = ctrl.customMinInfNorms;
        minInfNorms = ctrl.minInfNorms;
        customMaxInfNorms = ctrl.customMaxInfNorms;
        maxInfNorms = ctrl.maxInfNorms;

        customMinOneNorms = ctrl.customMinOneNorms;
        minOneNorms = ctrl.minOneNorms;
        customMaxOneNorms = ctrl.customMaxOneNorms;
        maxOneNorms = ctrl.maxOneNorms;

        progressLevel = ctrl.progressLevel;

        return *this;
    }

    EnumCtrl() { }
    EnumCtrl( const EnumCtrl<Real>& ctrl ) { *this = ctrl; }
    template<typename OtherReal>
    EnumCtrl( const EnumCtrl<OtherReal>& ctrl ) { *this = ctrl; }
};

namespace svp {

// The "GNR" enumeration algorithm has been extended to support complex
// arithmetic below by reformulating the original algorithm in terms of 
// a vector of states, each of which allows for a simple traversal of a 
// discrete spiral in the complex plane that is centered about a particular
// point (with the initial direction and orientation of the spiral chosen
// based upon the difference between the center and its rounding).
//
// While it is possible to easily analytically solve for the "next closest"
// integer to a given real number relative to a given integer with O(1) work,
// and such a  solution can easily drive a walk "outward" from a center point, 
// the author is not aware of an O(1) analogue in the complex plane.
//
// The easier part of generalizing GNR enumeration is in generalizing the notion
// of traversing Z modulo multiplication by -1 (which can be represented by Z+)
// to Z^2 modulo multiplication by i. In particular, the latter can be 
// represented by 
//
//     { (a,b) in Z^2 : a > 0, |b| < a } U {(0,0)},
//
// which can be easily traversed by trying each admissible b after incrementing
// a.

// An example of traversing a clockwise Z^2 spiral which begins rightward
// (which makes the most sense if the target center was rounded up and left)
// is given here:
//
//           20- ...
//           |
//           19  6 - 7 - 8 - 9
//           |   |           |
//           18  5   0 - 1   10
//           |   |       |   |
//           17  4 - 3 - 2   11
//           |               |
//           16- 15- 14- 13- 12
//
// where the root node, marked as 0, was the complex rounding of the
// (non-integer) center. Since the legs of the spiral are of length 
//
//   1, 1, 2, 2, 3, 3, etc.
//
// we can easily traverse the spiral with O(1) state. And, while such an
// ordering of Z^2 does not precisely equate with an ordering based upon the
// distance to the center node, it is a reasonable approximation and there are
// no tie-breaking subtleties.
//

template<typename Real>
struct SpiralState
{
private:
    bool constrained_;
    Real position_;

    // Only relevant if unconstrained
    bool forward_;
    Real jump_;

public:

    void Initialize( bool constrain=false )
    {
        constrained_ = constrain;
        position_ = 0;
        if( constrained_ )
            return;

        forward_ = true;
        jump_ = 1;
    }

    void Initialize( const Real& center )
    {
        constrained_ = false;
        position_ = Round(center);
        forward_ = (position_ <= center);
        jump_ = 1;
    }

    Real Step()
    {
        if( constrained_ )
        {
            position_ += 1;
        }
        else
        {
            if( forward_ )
                position_ += jump_;
            else
                position_ -= jump_;
            forward_ = !forward_;
            jump_ += 1;
        }
        return position_;
    }
};

template<typename Real>
class SpiralState<Complex<Real>>
{
private:
    bool constrained_;
    Complex<Real> position_;

    // Only relevant if unconstrained
    Int legLength_;
    Int numSteps_;
    Int direction_; // 0=right, 1=down, 2=left, 3=up
    bool clockwise_;
    bool firstLeg_;

    Complex<Real> ConstrainedStep()
    {
        Real realPos = position_.real();
        Real imagPos = position_.imag();

        // The coset Z^2 modulo multiplication by i can be represented by the
        // set
        // 
        //    { (a,b) in Z^2 : a > 0, |b| < a } U {(0,0)}, 
        //
        // and the origin can be ignored for our purposes. We traverse this set
        // by searching each admissible value of b after incrementing a.

        if( imagPos > Real(0) )
        {
            if( imagPos == realPos-1 )
                imagPos = -1;
            else
                imagPos += 1;
        }
        else if( imagPos < Real(0) )
        {
            if( imagPos == -(realPos-1) )
            {
                realPos += 1;
                imagPos = 0;
            }
            else
            {
                imagPos -= 1;
            }
        }
        else
        {
            imagPos += 1;
        }
        position_ = Complex<Real>(realPos,imagPos);
        return position_;
    }

    Complex<Real> UnconstrainedStep()
    {
        Real realPos = position_.real();
        Real imagPos = position_.imag();
        if( direction_ == 0 )
            realPos += 1;
        else if( direction_ == 1 )
            imagPos -= 1;
        else if( direction_ == 2 )
            realPos -= 1;
        else if( direction_ == 3 )
            imagPos += 1;

        ++numSteps_;
        if( numSteps_ == legLength_ )
        {
            numSteps_ = 0;

            if( clockwise_ )
                direction_ = Mod( direction_+1, 4 );
            else
                direction_ = Mod( direction_-1, 4 );

            if( firstLeg_ )
            {
                firstLeg_ = false;
            }
            else
            {
                ++legLength_;
                firstLeg_ = true;
            }
        }

        position_ = Complex<Real>(realPos,imagPos);
        return position_;
    }

public:

    void Initialize( bool constrain=false )
    {
        constrained_ = constrain;
        position_ = 0;
        if( constrained_ )
            return;    

        legLength_ = 1;
        numSteps_ = 0;
        firstLeg_ = true;

        direction_ = 1; // right
        clockwise_ = true;
    }

    void Initialize( const Complex<Real>& center )
    {
        constrained_ = false;

        legLength_ = 1;
        numSteps_ = 0;
        firstLeg_ = true;

        const Real realCenter = center.real();
        const Real imagCenter = center.imag();
        const Real realStart = Round(realCenter);
        const Real imagStart = Round(imagCenter);
        position_ = Complex<Real>(realStart,imagStart);
        const Real realDist = Abs(realCenter-realStart);
        const Real imagDist = Abs(imagCenter-imagStart);
        if( realDist < imagDist )
        {
            if( realCenter > realStart )
            {
                direction_ = 0; // right
                clockwise_ = (imagCenter < imagStart);
            }
            else
            {
                direction_ = 2; // left
                clockwise_ = (imagCenter > imagStart);
            }
        }
        else
        {
            if( imagCenter > imagStart )
            {
                direction_ = 3; // up
                clockwise_ = (realCenter > realStart);
            }
            else
            {
                direction_ = 2; // left
                clockwise_ = (realCenter < realStart);
            }
        }
    }

    Complex<Real> Step()
    {
        if( constrained_ )
            return ConstrainedStep();
        else
            return UnconstrainedStep();
    }
};

// If successful, fills 'v' with the integer coordinates of the columns of 
// an m x n matrix B (represented by its n x n upper-triangular Gaussian Normal
// Form; the 'R' from the QR factorization) which had a norm profile
// underneath the vector 'u' of upper bounds (|| (B v)(0:j) ||_2 < u(j)).
// Notice that the inequalities are strict.
//
// If not successful, the return value is a value greater than u(n-1) and 
// the contents of 'v' should be ignored.
//
// NOTE: There is not currently a complex implementation, though algorithms
//       exist.
template<typename F>
Base<F> GNREnumeration
( const Matrix<Base<F>>& d,
  const Matrix<F>& N,
  const Matrix<Base<F>>& u,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

// Convert to/from the so-called "y-sparse" representation of
//
//   Dan Ding, Guizhen Zhu, Yang Yu, and Zhongxiang Zheng,
//   "A fast phase-based enumeration algorithm for SVP challenge through
//    y-sparse representations of short lattice vectors"

template<typename F>
void CoordinatesToSparse
( const Matrix<F>& N, const Matrix<F>& v, Matrix<F>& y );
template<typename F>
void TransposedCoordinatesToSparse
( const Matrix<F>& NTrans, const Matrix<F>& v, Matrix<F>& y );

template<typename F>
void SparseToCoordinates
( const Matrix<F>& N, const Matrix<F>& y, Matrix<F>& v );
template<typename F>
void TransposedSparseToCoordinates
( const Matrix<F>& NTrans, const Matrix<F>& y, Matrix<F>& v );

template<typename F>
Base<F> CoordinatesToNorm
( const Matrix<Base<F>>& d, const Matrix<F>& N, const Matrix<F>& v );
template<typename F>
Base<F> TransposedCoordinatesToNorm
( const Matrix<Base<F>>& d, const Matrix<F>& NTrans, const Matrix<F>& v );

template<typename F>
Base<F> SparseToNorm
( const Matrix<Base<F>>& d, const Matrix<F>& N, const Matrix<F>& y );
template<typename F>
Base<F> TransposedSparseToNorm
( const Matrix<Base<F>>& d, const Matrix<F>& NTrans, const Matrix<F>& y );

template<typename F>
Base<F> PhaseEnumeration
( const Matrix<F>& B,
  const Matrix<Base<F>>& d,
  const Matrix<F>& N,
        Base<F> normUpperBound,
        Int startIndex,
        Int phaseLength,
        double enqueueProb,
  const vector<Int>& minInfNorms,
  const vector<Int>& maxInfNorms,
  const vector<Int>& minOneNorms,
  const vector<Int>& maxOneNorms,
        Matrix<F>& v,
        Int progressLevel=0 );
template<typename F>
std::pair<Base<F>,Int>
PhaseEnumeration
( const Matrix<F>& B,
  const Matrix<Base<F>>& d,
  const Matrix<F>& N,
  const Matrix<Base<F>>& normUpperBounds,
        Int startIndex,
        Int phaseLength,
        double enqueueProb,
  const vector<Int>& minInfNorms,
  const vector<Int>& maxInfNorms,
  const vector<Int>& minOneNorms,
  const vector<Int>& maxOneNorms,
        Matrix<F>& v,
        Int progressLevel=0 );

} // namespace svp

// Given a reduced lattice B and its Gaussian Normal Form, R, either find a
// member of the lattice (given by B v, with v the output) with norm less than
// the upper bound (and return its norm), or return a value greater than the
// upper bound.
template<typename F>
Base<F> ShortVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
        Base<F> normUpperBound,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

template<typename F>
std::pair<Base<F>,Int>
MultiShortVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
  const Matrix<Base<F>>& normUpperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

// Given a reduced lattice B and its Gaussian Normal Form, R, find the shortest
// member of the lattice (with the shortest vector given by B v).
//
// The return value is the norm of the (approximately) shortest vector.
template<typename F>
Base<F> ShortestVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

// If an upper-bound on the shortest vector which is better than || b_0 ||_2 is
// available
template<typename F>
Base<F> ShortestVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
        Base<F> normUpperBound,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

template<typename F>
std::pair<Base<F>,Int>
MultiShortestVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
  const Matrix<Base<F>>& normUpperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

// If a shorter vector is found, insert it into the first position
// ---------------------------------------------------------------

// The return value is the norm of the (approximately) shortest vector.
template<typename F>
Base<F> ShortestVectorEnrichment
(       Matrix<F>& B,
  const Matrix<F>& R,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );
template<typename F>
Base<F> ShortestVectorEnrichment
(       Matrix<F>& B,
        Matrix<F>& U,
  const Matrix<F>& R,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

// If an upper-bound on the shortest vector which is better than || b_0 ||_2 is
// available
template<typename F>
Base<F> ShortestVectorEnrichment
(       Matrix<F>& B,
  const Matrix<F>& R,
        Base<F> normUpperBound,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );
template<typename F>
Base<F> ShortestVectorEnrichment
(       Matrix<F>& B,
        Matrix<F>& U,
  const Matrix<F>& R,
        Base<F> normUpperBound,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

template<typename F>
std::pair<Base<F>,Int>
MultiShortestVectorEnrichment
(       Matrix<F>& B,
  const Matrix<F>& R,
  const Matrix<Base<F>>& normUpperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );
template<typename F>
std::pair<Base<F>,Int>
MultiShortestVectorEnrichment
(       Matrix<F>& B,
        Matrix<F>& U,
  const Matrix<F>& R,
  const Matrix<Base<F>>& normUpperBounds,
        Matrix<F>& v,
  const EnumCtrl<Base<F>>& ctrl=EnumCtrl<Base<F>>() );

} // namespace El

#endif // ifndef EL_LATTICE_ENUMERATE_HPP
