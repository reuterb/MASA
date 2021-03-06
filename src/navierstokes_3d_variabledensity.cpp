// -*-c++-*-
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010,2011,2012,2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#include <masa_internal.h>

#ifdef HAVE_METAPHYSICL

#include <ad_masa.h>

// Private methods declarations
template <typename Scalar, typename Scalar2>
Scalar helper_zeta(Scalar2 kx, Scalar2 kz, Scalar x, Scalar z);

template <typename Scalar, typename Scalar2>
Scalar helper_psi(Scalar2 kx, Scalar2 kz, Scalar y, Scalar2 noSlip);
  
template <typename Scalar, typename Scalar2>
Scalar helper_phi(Scalar2 kx, Scalar2 kz, Scalar y, Scalar2 noSlip, Scalar2 addMean);

template <typename Scalar, typename Scalar2>
Scalar helper_zetaAlpha(Scalar2 kx1, Scalar2 kz1, Scalar2 kx2, Scalar2 kz2,
                        Scalar2 kx3, Scalar2 kz3, Scalar2 kx4, Scalar2 kz4,
                        Scalar2 kx5, Scalar2 kz5, Scalar2 numModes,
                        Scalar x, Scalar y, Scalar z);

template <typename Scalar, typename Scalar2>
Scalar helper_zetaDelta(Scalar2 kx1, Scalar2 kz1, Scalar2 kx2, Scalar2 kz2,
                        Scalar2 kx3, Scalar2 kz3, Scalar2 kx4, Scalar2 kz4,
                        Scalar2 kx5, Scalar2 kz5, Scalar2 numModes,
                        Scalar x, Scalar y, Scalar z);

template <typename Scalar, typename Scalar2>
Scalar helper_zetaPhi(Scalar2 kx1, Scalar2 kz1, Scalar2 kx2, Scalar2 kz2,
                      Scalar2 kx3, Scalar2 kz3, Scalar2 kx4, Scalar2 kz4,
                      Scalar2 kx5, Scalar2 kz5, Scalar2 numModes, 
                      Scalar2 noSlip, Scalar2 addMean,
                      Scalar x, Scalar y, Scalar z);

template <typename Scalar, typename Scalar2>
Scalar helper_zetaPsi(Scalar2 kx1, Scalar2 kz1, Scalar2 kx2, Scalar2 kz2,
                      Scalar2 kx3, Scalar2 kz3, Scalar2 kx4, Scalar2 kz4,
                      Scalar2 kx5, Scalar2 kz5, Scalar2 numModes,
                      Scalar2 noSlip,
                      Scalar x, Scalar y, Scalar z);

template <typename Scalar, typename Scalar2>
Scalar helper_zetaGammaPlus(Scalar2 kx1, Scalar2 kz1, Scalar2 kx2, Scalar2 kz2,
                            Scalar2 kx3, Scalar2 kz3, Scalar2 kx4, Scalar2 kz4,
                            Scalar2 kx5, Scalar2 kz5, Scalar2 numModes,
                            Scalar x, Scalar y, Scalar z);


template <typename Scalar, typename Scalar2>
Scalar helper_zetaPsiCases(Scalar2 kx1, Scalar2 kz1, Scalar2 kx2, Scalar2 kz2,
                           Scalar2 kx3, Scalar2 kz3, Scalar2 kx4, Scalar2 kz4,
                           Scalar2 kx5, Scalar2 kz5, Scalar2 numModes,
                           Scalar2 noSlip,
                           Scalar x, Scalar y, Scalar z, bool m3Flag);


template <typename Scalar, typename Scalar2>
Scalar helper_zetaPsiTimesK(Scalar2 kx1, Scalar2 kz1, Scalar2 kx2, Scalar2 kz2,
                            Scalar2 kx3, Scalar2 kz3, Scalar2 kx4, Scalar2 kz4,
                            Scalar2 kx5, Scalar2 kz5, Scalar2 numModes,
                            Scalar2 noSlip,
                            Scalar x, Scalar y, Scalar z, int BCswitch);

template <typename Scalar>
Scalar helper_alpha(Scalar y);

template <typename Scalar>
Scalar helper_beta(Scalar y);

template <typename Scalar>
Scalar helper_gamma(Scalar y);

template <typename Scalar>
Scalar helper_delta(Scalar y);

template <typename Scalar>
Scalar helper_T(Scalar t);

template <typename Scalar>
Scalar rhoMapVD(Scalar Z, Scalar At);

template <typename Scalar>
Scalar drhoMapVD(Scalar Z, Scalar At);

//template <typename Scalar>
//NumberVector DivTauij(Scalar x1, Scalar y1, Scalar z1, Scalar t)

typedef ShadowNumber<double, long double> RawScalar;
const unsigned int NDIM = 3;  // Need time dependence!! But put in manually!
typedef DualNumber<RawScalar, NumberVector<NDIM, RawScalar> > D1Type;
typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
typedef D2Type ADType;

using namespace MASA;

template <typename Scalar>
MASA::navierstokes_3d_variabledensity<Scalar>::navierstokes_3d_variabledensity()
{
  this->mmsname = "navierstokes_3d_variabledensity";
  this->dimension = 4;

  this->register_var("re",&re);
  this->register_var("sc",&sc);
  this->register_var("at",&at);
  this->register_var("kx1",&kx1);
  this->register_var("kz1",&kz1);
  this->register_var("kx2",&kx2);
  this->register_var("kz2",&kz2);
  this->register_var("kx3",&kx3);
  this->register_var("kz3",&kz3);
  this->register_var("kx4",&kx4);
  this->register_var("kz4",&kz4);
  this->register_var("kx5",&kx5);
  this->register_var("kz5",&kz5);
  this->register_var("zero",&zero);
  this->register_var("one",&one);
 
  this->register_var("numModes",&numModes);

  this->register_var("noSlip",&noSlip);
  this->register_var("addMean",&addMean);

  this->init_var();

} // done with constructor

template <typename Scalar>
int MASA::navierstokes_3d_variabledensity<Scalar>::init_var()
{
  int err = 0;

  double kx,kz,kMag;

  kx = 1.0;
  kz = 1.0;

  err += this->set_var("re",1.0);
  err += this->set_var("sc",1.0);
  err += this->set_var("at",0.0);
  err += this->set_var("kx1",kx);
  err += this->set_var("kz1",kz);
  err += this->set_var("kx2",kx);
  err += this->set_var("kz2",kz); 
  err += this->set_var("kx3",kx);
  err += this->set_var("kz3",kz); 
  err += this->set_var("kx4",kx);
  err += this->set_var("kz4",kz); 
  err += this->set_var("kx5",kx);
  err += this->set_var("kz5",kz); 
  err += this->set_var("numModes",1); 
  err += this->set_var("noSlip",0);
  err += this->set_var("addMean",0);
  err += this->set_var("zero",0.);
  err += this->set_var("one",1.);

  return err;

} // done with init_var

// ----------------------------------------
// Source Terms
// ----------------------------------------

// public static method, that can be called from eval_q_t
template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_q_omega(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef DualNumber<D3Type, NumberVector<NDIM, D3Type> > D4Type;
  typedef D4Type ADScalar;

  typedef DualNumber<Scalar, NumberVector<NDIM+1, Scalar> > D1TimeType;
  typedef DualNumber<D1TimeType, NumberVector<NDIM+1, D1TimeType> > D2TimeType;
  typedef DualNumber<D2TimeType, NumberVector<NDIM+1, D2TimeType> > D3TimeType;
  typedef DualNumber<D3TimeType, NumberVector<NDIM+1, D3TimeType> > D4TimeType;
  typedef D4TimeType ADTimeScalar;

  // Treat velocity, momentum as a vector
  NumberVector<NDIM, D3Type> U,mD,mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  
  ADScalar zetaPhi = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,noSlip,
                   addMean,x,y,z);
  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mD[0] =    zetaPhi.derivatives()[1]   * helper_T(t1);
  mD[1] = (  zetaPhi.derivatives()[2] 
           - zetaPhi.derivatives()[0] ) * helper_T(t1);
  mD[2] =  - zetaPhi.derivatives()[1]   * helper_T(t1);

  mC = gradient(zetaPsi)*helper_T(t1);

  D3Type x3 = D3Type(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  D3Type y3 = D3Type(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  D3Type z3 = D3Type(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  D3Type rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                 x3,y3,z3) * helper_T(t1); 
  D3Type mu  = helper_beta(y3); 

  U = (mD + mC) / rho;

  // Get omega
  ADTimeScalar xt = 
    ADTimeScalar(x1,NumberVectorUnitVector<NDIM+1, 0, Scalar>::value());
  ADTimeScalar yt = 
    ADTimeScalar(y1,NumberVectorUnitVector<NDIM+1, 1, Scalar>::value());
  ADTimeScalar zt = 
    ADTimeScalar(z1,NumberVectorUnitVector<NDIM+1, 2, Scalar>::value());
  ADTimeScalar tt = 
    ADTimeScalar(t1,NumberVectorUnitVector<NDIM+1, 3, Scalar>::value());

  ADTimeScalar zetaPhiTime = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,noSlip,
                   addMean,xt,yt,zt)*helper_T(tt);
 
  // Omega = d m_1^d / dz - d m_3^d / dx

  //ADTimeScalar omega = 
  D2TimeType omega =
        // d m_1^d / dz
        zetaPhiTime.derivatives()[1].derivatives()[2]
    +
        // - d m_3^d / dz
        zetaPhiTime.derivatives()[1].derivatives()[0];

  // get DivCij

  NumberVector<NDIM, D2Type> DivC = -divergence(rho*U.outerproduct(U));

  // get DivTauij

  NumberVector<NDIM, NumberVector<NDIM, Scalar> > I = 
    NumberVector<NDIM, Scalar>::identity();
  NumberVector<NDIM, D2Type> DivTau = 
    divergence(
        mu * (gradient(U) + transpose(gradient(U)) 
              - 2./3. * divergence(U)*I) );

  // Omega equation residuals
  Scalar Q_omega = 
    raw_value(  
                // Temporal part
                omega.derivatives()[3]

        ) +
    raw_value(
                // Convective part
              - DivC[0].derivatives()[2] 
              + DivC[2].derivatives()[0]
                // Viscous part
              - 1./re*DivTau[0].derivatives()[2]
              + 1./re*DivTau[2].derivatives()[0]
              
              );

  return Q_omega;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_q_phi(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef DualNumber<D3Type, NumberVector<NDIM, D3Type> > D4Type;
  typedef DualNumber<D4Type, NumberVector<NDIM, D4Type> > D5Type;
  typedef D5Type ADScalar;

  typedef DualNumber<Scalar, NumberVector<NDIM+1, Scalar> > D1TimeType;
  typedef DualNumber<D1TimeType, NumberVector<NDIM+1, D1TimeType> > D2TimeType;
  typedef DualNumber<D2TimeType, NumberVector<NDIM+1, D2TimeType> > D3TimeType;
  typedef DualNumber<D3TimeType, NumberVector<NDIM+1, D3TimeType> > D4TimeType;
  typedef DualNumber<D4TimeType, NumberVector<NDIM+1, D4TimeType> > D5TimeType;
  typedef D5TimeType ADTimeScalar;

  // Treat velocity, momentum as a vector
  NumberVector<NDIM, D4Type> U,mD,mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPhi = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,noSlip,
                   addMean,x,y,z);

  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mD[0] =    zetaPhi.derivatives()[1]   * helper_T(t1);
  mD[1] = (  zetaPhi.derivatives()[2] 
           - zetaPhi.derivatives()[0] ) * helper_T(t1);
  mD[2] =  - zetaPhi.derivatives()[1]   * helper_T(t1);

  mC = gradient(zetaPsi)*helper_T(t1);

  D4Type x4 = D4Type(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  D4Type y4 = D4Type(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  D4Type z4 = D4Type(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  D4Type rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                 x4,y4,z4) * helper_T(t1);
  D4Type mu  = helper_beta(y4);

  U = (mD + mC) / rho;

  // Get phi
  ADTimeScalar xt = 
    ADTimeScalar(x1,NumberVectorUnitVector<NDIM+1, 0, Scalar>::value());
  ADTimeScalar yt = 
    ADTimeScalar(y1,NumberVectorUnitVector<NDIM+1, 1, Scalar>::value());
  ADTimeScalar zt = 
    ADTimeScalar(z1,NumberVectorUnitVector<NDIM+1, 2, Scalar>::value());
  ADTimeScalar tt = 
    ADTimeScalar(t1,NumberVectorUnitVector<NDIM+1, 3, Scalar>::value());

  // phi = lap( m_2^d ) 
   
  ADTimeScalar zetaPhiTime = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,noSlip,
                   addMean,xt,yt,zt)*helper_T(tt);

  //ADTimeScalar m2_D = 
  D4TimeType m2_D =
      zetaPhiTime.derivatives()[2] 
    - 
      zetaPhiTime.derivatives()[0];

  //ADTimeScalar phi = 
  D2TimeType phi = 
          m2_D.derivatives()[0].derivatives()[0]
        + m2_D.derivatives()[1].derivatives()[1]
        + m2_D.derivatives()[2].derivatives()[2];

  // get DivCij

  NumberVector<NDIM, D3Type > DivC = -divergence(rho*U.outerproduct(U));

  // get DivTauij

  NumberVector<NDIM, NumberVector<NDIM, Scalar> > I = 
    NumberVector<NDIM, Scalar>::identity();
  NumberVector<NDIM, D3Type> DivTau = 
    divergence(
        mu * (gradient(U) + transpose(gradient(U)) 
              - 2./3. * divergence(U)*I) );

  // get divergences
  
  D2Type DivDivC   = divergence(DivC);
  D2Type DivDivTau = divergence(DivTau);

  // phi equation residual
  Scalar Q_phi = 
    raw_value(  
                // Temporal part
                phi.derivatives()[3]

        ) +
    raw_value(
                // Convective part
              - divergence(gradient(DivC[1]))
              + DivDivC.derivatives()[1]
                // Viscous part
              - 1./re*divergence(gradient(DivTau[1]))
              + 1./re*DivDivTau.derivatives()[1]
                // gravity
              + ( rho.derivatives()[0].derivatives()[0] + 
                  rho.derivatives()[2].derivatives()[2] )

              );

  return Q_phi;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_q_rho(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef D2Type ADScalar;

  typedef DualNumber<Scalar, NumberVector<NDIM+1, Scalar> > D1TimeType;
  typedef D1TimeType ADTimeScalar;

  // Treat momentum as a vector
  NumberVector<NDIM, D1Type> mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mC = gradient(zetaPsi)*helper_T(t1);

  // Time part

  ADTimeScalar xt = 
    ADTimeScalar(x1,NumberVectorUnitVector<NDIM+1, 0, Scalar>::value());
  ADTimeScalar yt = 
    ADTimeScalar(y1,NumberVectorUnitVector<NDIM+1, 1, Scalar>::value());
  ADTimeScalar zt = 
    ADTimeScalar(z1,NumberVectorUnitVector<NDIM+1, 2, Scalar>::value());
  ADTimeScalar tt = 
    ADTimeScalar(t1,NumberVectorUnitVector<NDIM+1, 3, Scalar>::value());

  ADTimeScalar rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                 xt,yt,zt) * helper_T(tt);

  // Rho equation residuals
  Scalar Q_rho = 
                // Temporal part
                rho.derivatives()[3]

        +
                // Div mC
                divergence(mC);

  return Q_rho;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_q_z(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef D2Type ADScalar;

  typedef DualNumber<Scalar, NumberVector<NDIM+1, Scalar> > D1TimeType;
  typedef D1TimeType ADTimeScalar;

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U,mC,mD;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPhi = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,noSlip,
                   addMean,x,y,z);
 
  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mD[0] =    zetaPhi.derivatives()[1]   * helper_T(t1);
  mD[1] = (  zetaPhi.derivatives()[2] 
           - zetaPhi.derivatives()[0] ) * helper_T(t1);
  mD[2] =  - zetaPhi.derivatives()[1]   * helper_T(t1);

  mC = gradient(zetaPsi)*helper_T(t1);

  ADScalar rho  = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                   x,y,z) * helper_T(t1);
  ADScalar mu   = helper_beta(y);
  ADScalar zVar = 
    helper_zetaGammaPlus(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                         x,y,z) * helper_T(t1);

  U = (mD + mC) / rho;

  // Time part

  ADTimeScalar xt = 
    ADTimeScalar(x1,NumberVectorUnitVector<NDIM+1, 0, Scalar>::value());
  ADTimeScalar yt = 
    ADTimeScalar(y1,NumberVectorUnitVector<NDIM+1, 1, Scalar>::value());
  ADTimeScalar zt = 
    ADTimeScalar(z1,NumberVectorUnitVector<NDIM+1, 2, Scalar>::value());
  ADTimeScalar tt = 
    ADTimeScalar(t1,NumberVectorUnitVector<NDIM+1, 3, Scalar>::value());

  ADTimeScalar zTime = 
    helper_zetaGammaPlus(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                         xt,yt,zt)*helper_T(tt);

  // Z equation residuals
  Scalar Q_z = 
                // Temporal part
                zTime.derivatives()[3]

        +
                // U dot grad(Z)
                // For new iterative method, we will 
                // need to be sure that the continuity solve
                // gives us back m^c (scaled appropriately)
                raw_value(U.dot(gradient(zVar)))
        
        -
                // viscous
       raw_value(   1.0/re*1.0/sc*(
                    1.0/rho*gradient(mu).dot(gradient(zVar))
                  +  mu/rho*divergence(gradient(zVar))));

  return Q_z;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_q_m1(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef D3Type ADScalar;

  typedef DualNumber<Scalar, NumberVector<NDIM+1, Scalar> > D1TimeType;
  typedef DualNumber<D1TimeType, NumberVector<NDIM+1, D1TimeType> > D2TimeType;
  typedef D2TimeType ADTimeScalar;

  // Treat velocity,momentum as a vector
  NumberVector<NDIM, D2Type> mC,mD, U;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mC = gradient(zetaPsi)*helper_T(t1);

  ADScalar zetaPhi = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,noSlip,
                   addMean,x,y,z);

  mD[0] =    zetaPhi.derivatives()[1]   * helper_T(t1);
  mD[1] = (  zetaPhi.derivatives()[2] 
           - zetaPhi.derivatives()[0] ) * helper_T(t1);
  mD[2] =  - zetaPhi.derivatives()[1]   * helper_T(t1);

  D2Type x2 = D2Type(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  D2Type y2 = D2Type(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  D2Type z2 = D2Type(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  D2Type rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                 x2,y2,z2) * helper_T(t1); 
  D2Type mu  = helper_beta(y2); 

  U = (mD + mC) / rho;

  // Time part

  ADTimeScalar xt = 
    ADTimeScalar(x1,NumberVectorUnitVector<NDIM+1, 0, Scalar>::value());
  ADTimeScalar yt = 
    ADTimeScalar(y1,NumberVectorUnitVector<NDIM+1, 1, Scalar>::value());
  ADTimeScalar zt = 
    ADTimeScalar(z1,NumberVectorUnitVector<NDIM+1, 2, Scalar>::value());
  ADTimeScalar tt = 
    ADTimeScalar(t1,NumberVectorUnitVector<NDIM+1, 3, Scalar>::value());

  ADTimeScalar zetaPhiTime = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,noSlip,
                   addMean,xt,yt,zt)*helper_T(tt);

  ADTimeScalar m1_D = zetaPhiTime.derivatives()[1];

  ADTimeScalar zetaPsiTime = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,xt,yt,zt)*helper_T(tt);

  ADTimeScalar m1_C = zetaPsiTime.derivatives()[0];

  NumberVector<NDIM, NumberVector<NDIM, Scalar> > I = 
    NumberVector<NDIM, Scalar>::identity();
 
  NumberVector<NDIM, NumberVector<NDIM,D2Type> > Tau = 
   mu * (gradient(U) + transpose(gradient(U)) - 2./3.*divergence(U)*I);

  D2Type C12 = -rho*U[0]*U[1];

  // m1 equation residuals
  Scalar Q_m1 = 
                // Temporal part
                raw_value(m1_D.derivatives()[3] + m1_C.derivatives()[3])

        -
                // RHS
                raw_value(C12.derivatives()[1]
                + 0.*1./re*Tau[0][1].derivatives()[1]);

  return Q_m1;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_q_m3(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef D3Type ADScalar;

  typedef DualNumber<Scalar, NumberVector<NDIM+1, Scalar> > D1TimeType;
  typedef DualNumber<D1TimeType, NumberVector<NDIM+1, D1TimeType> > D2TimeType;
  typedef D2TimeType ADTimeScalar;

  // Treat velocity,momentum as a vector
  NumberVector<NDIM, D2Type> mC,mD, U;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mC = gradient(zetaPsi)*helper_T(t1);

  ADScalar zetaPhi = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,noSlip,
                   addMean,x,y,z);

  mD[0] =    zetaPhi.derivatives()[1]   * helper_T(t1);
  mD[1] = (  zetaPhi.derivatives()[2] 
           - zetaPhi.derivatives()[0] ) * helper_T(t1);
  mD[2] =  - zetaPhi.derivatives()[1]   * helper_T(t1);

  D2Type x2 = D2Type(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  D2Type y2 = D2Type(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  D2Type z2 = D2Type(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  D2Type rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                 x2,y2,z2) * helper_T(t1); 
  D2Type mu  = helper_beta(y2); 

  U = (mD + mC) / rho;

  // Time part

  ADTimeScalar xt = 
    ADTimeScalar(x1,NumberVectorUnitVector<NDIM+1, 0, Scalar>::value());
  ADTimeScalar yt = 
    ADTimeScalar(y1,NumberVectorUnitVector<NDIM+1, 1, Scalar>::value());
  ADTimeScalar zt = 
    ADTimeScalar(z1,NumberVectorUnitVector<NDIM+1, 2, Scalar>::value());
  ADTimeScalar tt = 
    ADTimeScalar(t1,NumberVectorUnitVector<NDIM+1, 3, Scalar>::value());

  ADTimeScalar zetaPhiTime = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,noSlip,
                   addMean,xt,yt,zt)*helper_T(tt);

  ADTimeScalar m3_D = -zetaPhiTime.derivatives()[1];

  ADTimeScalar zetaPsiTime = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,xt,yt,zt)*helper_T(tt);

  ADTimeScalar m3_C = zetaPsiTime.derivatives()[2];

  NumberVector<NDIM, NumberVector<NDIM, Scalar> > I = 
    NumberVector<NDIM, Scalar>::identity();
 
  NumberVector<NDIM, NumberVector<NDIM,D2Type> > Tau = 
   mu * (gradient(U) + transpose(gradient(U)) - 2./3.*divergence(U)*I);

  D2Type C23 = -rho*U[1]*U[2];

  // m3 equation residuals
  Scalar Q_m3 = 
                // Temporal part
                raw_value( m3_D.derivatives()[3] + m3_C.derivatives()[3] )

        -
                // RHS
                raw_value(C23.derivatives()[1] 
                + 1./re*Tau[1][2].derivatives()[1]);

  return Q_m3;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_q_p(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef DualNumber<D3Type, NumberVector<NDIM, D3Type> > D4Type;
  typedef D4Type ADScalar;

  typedef DualNumber<Scalar, NumberVector<NDIM+1, Scalar> > D1TimeType;
  typedef DualNumber<D1TimeType, NumberVector<NDIM+1, D1TimeType> > D2TimeType;
  typedef DualNumber<D2TimeType, NumberVector<NDIM+1, D2TimeType> > D3TimeType;
  typedef D3TimeType ADTimeScalar;

  // Treat velocity,momentum as a vector
  NumberVector<NDIM, D3Type> mC,mD, U;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mC = gradient(zetaPsi)*helper_T(t1);

  ADScalar zetaPhi = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,noSlip,
                   addMean,x,y,z);

  mD[0] =    zetaPhi.derivatives()[1]   * helper_T(t1);
  mD[1] = (  zetaPhi.derivatives()[2] 
           - zetaPhi.derivatives()[0] ) * helper_T(t1);
  mD[2] =  - zetaPhi.derivatives()[1]   * helper_T(t1);

  D3Type x3 = D3Type(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  D3Type y3 = D3Type(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  D3Type z3 = D3Type(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  D3Type rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                 x3,y3,z3) * helper_T(t1); 
  D3Type mu  = helper_beta(y3); 

  U = (mD + mC) / rho;

  D3Type p = helper_zetaDelta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                              x3,y3,z3) * helper_T(t1); 

  // get DivCij

  NumberVector<NDIM, D3Type > DivC = -divergence(rho*U.outerproduct(U));

  // Time part

  ADTimeScalar xt = 
    ADTimeScalar(x1,NumberVectorUnitVector<NDIM+1, 0, Scalar>::value());
  ADTimeScalar yt = 
    ADTimeScalar(y1,NumberVectorUnitVector<NDIM+1, 1, Scalar>::value());
  ADTimeScalar zt = 
    ADTimeScalar(z1,NumberVectorUnitVector<NDIM+1, 2, Scalar>::value());
  ADTimeScalar tt = 
    ADTimeScalar(t1,NumberVectorUnitVector<NDIM+1, 3, Scalar>::value());

//  ADTimeScalar rhoTime = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
//                                 xt,yt,zt) * helper_T(tt); 

  ADTimeScalar zetaPsiTime = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,xt,yt,zt)*helper_T(tt);

  ADTimeScalar zetaPhiTime = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,addMean,xt,yt,zt)*helper_T(tt);

  NumberVector<NDIM, NumberVector<NDIM, Scalar> > I = 
    NumberVector<NDIM, Scalar>::identity();
 
  NumberVector<NDIM, NumberVector<NDIM,D3Type> > Tau = 
   mu * (gradient(U) + transpose(gradient(U)) - 2./3.*divergence(U)*I);

  Scalar Q_p = 

    raw_value(
    // laplacian of pressure
    
    divergence(gradient(p)) +
    
    // temporal part
    
    (zetaPsiTime.derivatives()[0].derivatives()[0].derivatives()[3] + 
     zetaPsiTime.derivatives()[1].derivatives()[1].derivatives()[3] +
     zetaPsiTime.derivatives()[2].derivatives()[2].derivatives()[3] +
     zetaPhiTime.derivatives()[1].derivatives()[0].derivatives()[3] +
     zetaPhiTime.derivatives()[2].derivatives()[1].derivatives()[3] -
     zetaPhiTime.derivatives()[0].derivatives()[1].derivatives()[3] -
     zetaPhiTime.derivatives()[1].derivatives()[2].derivatives()[3]) -
    
    // convection

    divergence(DivC) -

    // viscous
    
    divergence(divergence(Tau)) );

  return Q_p;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_q_pBound(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef D3Type ADScalar;

  typedef DualNumber<Scalar, NumberVector<NDIM+1, Scalar> > D1TimeType;
  typedef DualNumber<D1TimeType, NumberVector<NDIM+1, D1TimeType> > D2TimeType;
  typedef D2TimeType ADTimeScalar;

  // Treat velocity,momentum as a vector
  NumberVector<NDIM, D2Type> mC,mD, U;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mC = gradient(zetaPsi)*helper_T(t1);

  ADScalar zetaPhi = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,noSlip,
                   addMean,x,y,z);

  mD[0] =    zetaPhi.derivatives()[1]   * helper_T(t1);
  mD[1] = (  zetaPhi.derivatives()[2] 
           - zetaPhi.derivatives()[0] ) * helper_T(t1);
  mD[2] =  - zetaPhi.derivatives()[1]   * helper_T(t1);

  D2Type x2 = D2Type(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  D2Type y2 = D2Type(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  D2Type z2 = D2Type(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  D2Type rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                 x2,y2,z2) * helper_T(t1); 
  D2Type mu  = helper_beta(y2); 

  U = (mD + mC) / rho;

  D2Type p = helper_zetaDelta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                              x2,y2,z2) * helper_T(t1); 

  // get DivCij

  NumberVector<NDIM, D2Type > DivC = -divergence(rho*U.outerproduct(U));

  // Time part

  ADTimeScalar xt = 
    ADTimeScalar(x1,NumberVectorUnitVector<NDIM+1, 0, Scalar>::value());
  ADTimeScalar yt = 
    ADTimeScalar(y1,NumberVectorUnitVector<NDIM+1, 1, Scalar>::value());
  ADTimeScalar zt = 
    ADTimeScalar(z1,NumberVectorUnitVector<NDIM+1, 2, Scalar>::value());
  ADTimeScalar tt = 
    ADTimeScalar(t1,NumberVectorUnitVector<NDIM+1, 3, Scalar>::value());

  ADTimeScalar zetaPhiTime = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,noSlip,
                   addMean,xt,yt,zt)*helper_T(tt);

  ADTimeScalar m2_D = (zetaPhiTime.derivatives()[2] -
                       zetaPhiTime.derivatives()[0]);

  ADTimeScalar zetaPsiTime = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,xt,yt,zt)*helper_T(tt);

  ADTimeScalar m2_C = zetaPsiTime.derivatives()[1];


  NumberVector<NDIM, NumberVector<NDIM, Scalar> > I = 
    NumberVector<NDIM, Scalar>::identity();
 
  NumberVector<NDIM, NumberVector<NDIM,D2Type> > Tau = 
   mu * (gradient(U) + transpose(gradient(U)) - 2./3.*divergence(U)*I);

  Scalar Q_pBound = 

    raw_value(
    // pressure gradient
    
    p.derivatives()[1] )  +
    
    // temporal part
    raw_value(
    m2_D.derivatives()[3] + m2_C.derivatives()[3] ) -

    // convection

    raw_value(

    DivC[1] + 

    // viscous

    divergence(Tau)[1] ) ;

  return Q_pBound;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_q_top(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef D2Type ADScalar;

  // Treat momentum as a vector
  NumberVector<NDIM, D1Type> mC_m1, mC_m3;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPsi_m1 = 
    helper_zetaPsiCases(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                        noSlip,x,y,z,false);

  ADScalar zetaPsi_m3 = 
    helper_zetaPsiCases(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                        noSlip,x,y,z,true);

  ADScalar zetaPsiTimesK_m3 = 
    helper_zetaPsiTimesK(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                         noSlip,x,y,z,1);

  ADScalar zetaPsiTimesK_m1 = 
    helper_zetaPsiTimesK(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                         noSlip,x,y,z,2);

  mC_m1 = gradient(zetaPsi_m1)*helper_T(t1);
  mC_m3 = gradient(zetaPsi_m3)*helper_T(t1);

  // BC residuals
  Scalar Q_top = raw_value(zetaPsiTimesK_m3.derivatives()[2])*helper_T(t1)
               + raw_value(zetaPsiTimesK_m1.derivatives()[0])*helper_T(t1)
       + mC_m1[1].derivatives()[0] 
       + mC_m3[1].derivatives()[2]; 

  return Q_top;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_q_bottom(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef D2Type ADScalar;

  // Treat momentum as a vector
  NumberVector<NDIM, D1Type> mC_m1, mC_m3;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPsi_m1 = 
    helper_zetaPsiCases(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                        noSlip,x,y,z,false);

  ADScalar zetaPsi_m3 = 
    helper_zetaPsiCases(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                        noSlip,x,y,z,true);

  ADScalar zetaPsiTimesK_m3 = 
    helper_zetaPsiTimesK(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                         noSlip,x,y,z,1);

  ADScalar zetaPsiTimesK_m1 = 
    helper_zetaPsiTimesK(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                         noSlip,x,y,z,2);

  mC_m1 = gradient(zetaPsi_m1)*helper_T(t1);
  mC_m3 = gradient(zetaPsi_m3)*helper_T(t1);

  // BC residuals
  Scalar Q_bottom = raw_value(zetaPsiTimesK_m3.derivatives()[2])*helper_T(t1)
                  + raw_value(zetaPsiTimesK_m1.derivatives()[0])*helper_T(t1)
       - mC_m1[1].derivatives()[0] 
       - mC_m3[1].derivatives()[2]; 

  return Q_bottom;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_q_phiTop(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef D2Type ADScalar;

  // Treat momentum as a vector
  NumberVector<NDIM, D1Type> mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  ADScalar zetaPsiTimesK = 
    helper_zetaPsiTimesK(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                         noSlip,x,y,z,3);

  mC = gradient(zetaPsi)*helper_T(t1);

  // BC residuals
  Scalar Q_top = raw_value(zetaPsiTimesK.derivatives()[1])*helper_T(t1)
       + mC[1].derivatives()[1];

  return Q_top;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_q_phiBottom(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef D2Type ADScalar;

  // Treat momentum as a vector
  NumberVector<NDIM, D1Type> mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  ADScalar zetaPsiTimesK = 
    helper_zetaPsiTimesK(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                         noSlip,x,y,z,3);

  mC = gradient(zetaPsi)*helper_T(t1);

  // BC residuals
  Scalar Q_bottom = raw_value(zetaPsiTimesK.derivatives()[1])*helper_T(t1)
       - mC[1].derivatives()[1];

  return Q_bottom;

}

// ----------------------------------------
// Analytical Terms
// ----------------------------------------
//
template <typename Scalar, typename Scalar2>
Scalar helper_zetaPsiTimesK(Scalar2 kx1, Scalar2 kz1, Scalar2 kx2, Scalar2 kz2,
                           Scalar2 kx3, Scalar2 kz3, Scalar2 kx4, Scalar2 kz4,
                           Scalar2 kx5, Scalar2 kz5, Scalar2 numModes,
                           Scalar2 noSlip,
                           Scalar x, Scalar y, Scalar z, int BCswitch)
{
  Scalar func = 0.0;

  int modes = int(numModes);

  if (BCswitch == 1) {

    // If were nabbing m3, we want to ignore the case with kz = 0 and
    // kx = anything but zero

    if ( ! ( (kz1 == 0.0) && (std::abs(kx1) > 1.0e-13) ) ) {
  
          func += (std::sqrt(kx1*kx1 + kz1*kz1))
                  *helper_zeta(kx1,kz1,x,z)*helper_psi(kx1,kz1,y,noSlip);
        }

    switch (modes) {
  
      case 1:

        // done!

        break;
  
      case 5:

         if ( ! ( (kz2 == 0.0) && (std::abs(kx2) > 1.0e-13) ) ) {
  
          func += (std::sqrt(kx2*kx2 + kz2*kz2))
                  *helper_zeta(kx2,kz2,x,z)*helper_psi(kx2,kz2,y,noSlip);
        }   

         if ( ! ( (kz3 == 0.0) && (std::abs(kx3) > 1.0e-13) ) ) {
  
          func += (std::sqrt(kx3*kx3 + kz3*kz3))
                  *helper_zeta(kx3,kz3,x,z)*helper_psi(kx3,kz3,y,noSlip);
        }   

        if ( ! ( (kz4 == 0.0) && (std::abs(kx4) > 1.0e-13) ) ) {
  
          func += (std::sqrt(kx4*kx4 + kz1*kz4))
                  *helper_zeta(kx4,kz4,x,z)*helper_psi(kx4,kz4,y,noSlip);
        }   

        if ( ! ( (kz5 == 0.0) && (std::abs(kx5) > 1.0e-13) ) ) {
  
          func += (std::sqrt(kx5*kx5 + kz5*kz5))
                  *helper_zeta(kx5,kz5,x,z)*helper_psi(kx5,kz5,y,noSlip);
        }   

        break;
  
    }
  
  } else if (BCswitch == 2) {
    
    // ONLY if kz = 0 and kx is non-zero do we solve for m1 instead
    
    if ( ( (kz1 == 0.0) && (std::abs(kx1) > 1.0e-13) ) ) {
  
          func += (std::sqrt(kx1*kx1 + kz1*kz1))
                  *helper_zeta(kx1,kz1,x,z)*helper_psi(kx1,kz1,y,noSlip);
        }

    switch (modes) {
  
      case 1:

        // done!

        break;
  
      case 5:

         if ( ( (kz2 == 0.0) && (std::abs(kx2) > 1.0e-13) ) ) {
  
          func += (std::sqrt(kx2*kx2 + kz2*kz2))
                  *helper_zeta(kx2,kz2,x,z)*helper_psi(kx2,kz2,y,noSlip);
        }   

         if ( ( (kz3 == 0.0) && (std::abs(kx3) > 1.0e-13) ) ) {
  
          func += (std::sqrt(kx3*kx3 + kz3*kz3))
                  *helper_zeta(kx3,kz3,x,z)*helper_psi(kx3,kz3,y,noSlip);
        }   

        if ( ( (kz4 == 0.0) && (std::abs(kx4) > 1.0e-13) ) ) {
  
          func += (std::sqrt(kx4*kx4 + kz1*kz4))
                  *helper_zeta(kx4,kz4,x,z)*helper_psi(kx4,kz4,y,noSlip);
        }   

        if ( ( (kz5 == 0.0) && (std::abs(kx5) > 1.0e-13) ) ) {
  
          func += (std::sqrt(kx5*kx5 + kz5*kz5))
                  *helper_zeta(kx5,kz5,x,z)*helper_psi(kx5,kz5,y,noSlip);
        }   

        break;
  
    }

  } else {

    
    func += (std::sqrt(kx1*kx1 + kz1*kz1))
              *helper_zeta(kx1,kz1,x,z)*helper_psi(kx1,kz1,y,noSlip);

    switch (modes) {
    
      case 1:

        // done!

        break;
    
      case 5:
    
        func += (std::sqrt(kx2*kx2 + kz2*kz2))
                *helper_zeta(kx2,kz2,x,z)*helper_psi(kx2,kz2,y,noSlip);
        func += (std::sqrt(kx3*kx3 + kz3*kz3))
                *helper_zeta(kx3,kz3,x,z)*helper_psi(kx3,kz3,y,noSlip);
        func += (std::sqrt(kx4*kx4 + kz1*kz4))
                *helper_zeta(kx4,kz4,x,z)*helper_psi(kx4,kz4,y,noSlip);
        func += (std::sqrt(kx5*kx5 + kz5*kz5))
                *helper_zeta(kx5,kz5,x,z)*helper_psi(kx5,kz5,y,noSlip);

        break;
    
      }

  }

  return func;

}

template <typename Scalar, typename Scalar2>
Scalar helper_zetaPsiCases(Scalar2 kx1, Scalar2 kz1, Scalar2 kx2, Scalar2 kz2,
                      Scalar2 kx3, Scalar2 kz3, Scalar2 kx4, Scalar2 kz4,
                      Scalar2 kx5, Scalar2 kz5, Scalar2 numModes,
                      Scalar2 noSlip,
                      Scalar x, Scalar y, Scalar z, bool m3Flag)
{

  Scalar func = 0.0;

  int modes = int(numModes);
  if (m3Flag == true) {

    // If were nabbing m3, we want to ignore the case with kz = 0 and
    // kx = anything but zero

    if ( ! ( (kz1 == 0.0) && (std::abs(kx1) > 1.0e-13) ) ) {
  
          func += helper_zeta(kx1,kz1,x,z)*helper_psi(kx1,kz1,y,noSlip);
        }

    switch (modes) {
  
      case 1:

        // done!

        break;
  
      case 5:

         if ( ! ( (kz2 == 0.0) && (std::abs(kx2) > 1.0e-13) ) ) {
  
          func += helper_zeta(kx2,kz2,x,z)*helper_psi(kx2,kz2,y,noSlip);
        }   

         if ( ! ( (kz3 == 0.0) && (std::abs(kx3) > 1.0e-13) ) ) {
  
          func += helper_zeta(kx3,kz3,x,z)*helper_psi(kx3,kz3,y,noSlip);
        }   

        if ( ! ( (kz4 == 0.0) && (std::abs(kx4) > 1.0e-13) ) ) {
  
          func += helper_zeta(kx4,kz4,x,z)*helper_psi(kx4,kz4,y,noSlip);
        }   

        if ( ! ( (kz5 == 0.0) && (std::abs(kx5) > 1.0e-13) ) ) {
  
          func += helper_zeta(kx5,kz5,x,z)*helper_psi(kx5,kz5,y,noSlip);
        }   

        break;
  
    }
  
  } else {
    
    // ONLY if kz = 0 and kx is non-zero do we solve for m1 instead
    
    if ( ( (kz1 == 0.0) && (std::abs(kx1) > 1.0e-13) ) ) {
  
          func += helper_zeta(kx1,kz1,x,z)*helper_psi(kx1,kz1,y,noSlip);
        }

    switch (modes) {
  
      case 1:

        // done!

        break;
  
      case 5:

         if ( ( (kz2 == 0.0) && (std::abs(kx2) > 1.0e-13) ) ) {
  
          func += helper_zeta(kx2,kz2,x,z)*helper_psi(kx2,kz2,y,noSlip);
        }   

         if ( ( (kz3 == 0.0) && (std::abs(kx3) > 1.0e-13) ) ) {
  
          func += helper_zeta(kx3,kz3,x,z)*helper_psi(kx3,kz3,y,noSlip);
        }   

        if ( ( (kz4 == 0.0) && (std::abs(kx4) > 1.0e-13) ) ) {
  
          func += helper_zeta(kx4,kz4,x,z)*helper_psi(kx4,kz4,y,noSlip);
        }   

        if ( ( (kz5 == 0.0) && (std::abs(kx5) > 1.0e-13) ) ) {
  
          func += helper_zeta(kx5,kz5,x,z)*helper_psi(kx5,kz5,y,noSlip);
        }   

        break;
  
    }

  }

  return func;

}


template <typename Scalar, typename Scalar2>
Scalar helper_zetaPsi(Scalar2 kx1, Scalar2 kz1, Scalar2 kx2, Scalar2 kz2,
                      Scalar2 kx3, Scalar2 kz3, Scalar2 kx4, Scalar2 kz4,
                      Scalar2 kx5, Scalar2 kz5, Scalar2 numModes,
                      Scalar2 noSlip,
                      Scalar x, Scalar y, Scalar z)
{

  Scalar func = 0.0;

  int modes = int(numModes);

  switch (modes) {

    case 1:

      func += helper_zeta(kx1,kz1,x,z)*helper_psi(kx1,kz1,y,noSlip);
      break;

    case 5:

      func += helper_zeta(kx1,kz1,x,z)*helper_psi(kx1,kz1,y,noSlip)
           +  helper_zeta(kx2,kz2,x,z)*helper_psi(kx2,kz2,y,noSlip)
           +  helper_zeta(kx3,kz3,x,z)*helper_psi(kx3,kz3,y,noSlip)
           +  helper_zeta(kx4,kz4,x,z)*helper_psi(kx4,kz4,y,noSlip)
           +  helper_zeta(kx5,kz5,x,z)*helper_psi(kx5,kz5,y,noSlip);
      break;

  }

  return func;

}

template <typename Scalar, typename Scalar2>
Scalar helper_zetaPhi(Scalar2 kx1, Scalar2 kz1, Scalar2 kx2, Scalar2 kz2,
                      Scalar2 kx3, Scalar2 kz3, Scalar2 kx4, Scalar2 kz4,
                      Scalar2 kx5, Scalar2 kz5, Scalar2 numModes, 
                      Scalar2 noSlip, Scalar2 addMean,
                      Scalar x, Scalar y, Scalar z)
{

  Scalar func = 0.0;

  int modes = int(numModes);

  switch (modes) {

    case 1:

      func += helper_zeta(kx1,kz1,x,z)*helper_phi(kx1,kz1,y,noSlip,addMean);
      break;

    case 5:

      func += helper_zeta(kx1,kz1,x,z)*helper_phi(kx1,kz1,y,noSlip,addMean)
           +  helper_zeta(kx2,kz2,x,z)*helper_phi(kx2,kz2,y,noSlip,addMean)
           +  helper_zeta(kx3,kz3,x,z)*helper_phi(kx3,kz3,y,noSlip,addMean)
           +  helper_zeta(kx4,kz4,x,z)*helper_phi(kx4,kz4,y,noSlip,addMean)
           +  helper_zeta(kx5,kz5,x,z)*helper_phi(kx5,kz5,y,noSlip,addMean);
      break;

  }

  return func;

}

template <typename Scalar, typename Scalar2>
Scalar helper_zetaGammaPlus(Scalar2 kx1, Scalar2 kz1, Scalar2 kx2, Scalar2 kz2,
                            Scalar2 kx3, Scalar2 kz3, Scalar2 kx4, Scalar2 kz4,
                            Scalar2 kx5, Scalar2 kz5, Scalar2 numModes, 
                            Scalar x, Scalar y, Scalar z)
{

  Scalar func = 0.5;

  int modes = int(numModes);

  switch (modes) {

    case 1:

      func += helper_zeta(kx1,kz1,x,z)*helper_gamma(y);
      break;

    case 5:

      func += (helper_zeta(kx1,kz1,x,z)*helper_gamma(y)
           +   helper_zeta(kx2,kz2,x,z)*helper_gamma(y)
           +   helper_zeta(kx3,kz3,x,z)*helper_gamma(y)
           +   helper_zeta(kx4,kz4,x,z)*helper_gamma(y)
           +   helper_zeta(kx5,kz5,x,z)*helper_gamma(y))/5.0;
      break;

  }

  return func;

}

template <typename Scalar, typename Scalar2>
Scalar helper_zetaAlpha(Scalar2 kx1, Scalar2 kz1, Scalar2 kx2, Scalar2 kz2,
                        Scalar2 kx3, Scalar2 kz3, Scalar2 kx4, Scalar2 kz4,
                        Scalar2 kx5, Scalar2 kz5, Scalar2 numModes, 
                        Scalar x, Scalar y, Scalar z)
{

  Scalar func = 6.; 

  return func; 

  return func * helper_alpha(y);

  int modes = int(numModes);

  switch (modes) {

    case 1:

      func += .6*helper_zeta(kx1,kz1,x,z)*helper_alpha(y);
      //func += .6*helper_alpha(y);
      //func += .1*helper_zeta(kx1,kz1,x,z)*helper_alpha(y);
      break;

    case 5:

      func += (helper_zeta(kx1,kz1,x,z)*helper_alpha(y)
           +   helper_zeta(kx2,kz2,x,z)*helper_alpha(y)
           +   helper_zeta(kx3,kz3,x,z)*helper_alpha(y)
           +   helper_zeta(kx4,kz4,x,z)*helper_alpha(y)
           +   helper_zeta(kx5,kz5,x,z)*helper_alpha(y))/5.0;
      break;

  }

  return func;

}

template <typename Scalar, typename Scalar2>
Scalar helper_zetaDelta(Scalar2 kx1, Scalar2 kz1, Scalar2 kx2, Scalar2 kz2,
                        Scalar2 kx3, Scalar2 kz3, Scalar2 kx4, Scalar2 kz4,
                        Scalar2 kx5, Scalar2 kz5, Scalar2 numModes, 
                        Scalar x, Scalar y, Scalar z)
{

  Scalar func = 0.;

  int modes = int(numModes);

  switch (modes) {

    case 1:

      func += helper_zeta(kx1,kz1,x,z)*helper_delta(y);
      break;

    case 5:

      func += (helper_zeta(kx1,kz1,x,z)*helper_delta(y)
           +   helper_zeta(kx2,kz2,x,z)*helper_delta(y)
           +   helper_zeta(kx3,kz3,x,z)*helper_delta(y)
           +   helper_zeta(kx4,kz4,x,z)*helper_delta(y)
           +   helper_zeta(kx5,kz5,x,z)*helper_delta(y));
      break;

  }

  return func;

}

// example of a private method, called from exact_t
template <typename Scalar, typename Scalar2>
Scalar helper_zeta(Scalar2 kx, Scalar2 kz, Scalar x, Scalar z)
{
  Scalar func;
  func = 2.*( std::cos(kx*x)*std::cos(kz*z) - std::sin(kx*x)*std::sin(kz*z)
             -std::sin(kx*x)*std::cos(kz*z) - std::cos(kx*x)*std::sin(kz*z) );
  return func;
}

template <typename Scalar, typename Scalar2>
Scalar helper_psi(Scalar2 kx, Scalar2 kz, Scalar y, Scalar2 noSlip)
{
  Scalar func;

  Scalar2 kMag = std::sqrt(kx*kx + kz*kz);

  Scalar2 a_psi, b_psi, c_psi, d_psi, e_psi;

  if ( kMag == 0. ) {
    
    a_psi = -1.0;
    b_psi =  0.0;
    c_psi =  3.0;
    d_psi =  0.0;
    e_psi =  0.0;

  } else {

    if ( noSlip == 0.) {
  
      a_psi = - (kMag + 1.0)/2.0;
      b_psi = 0.0;
      c_psi = 1.0 - a_psi;
      d_psi = 0.0;
      e_psi = 1.0/2.0*(-2.0/kMag - 2.0/a_psi - c_psi); 

    } else {

      a_psi =  4.0;
      b_psi =  0.0;
      c_psi = -4.0;
      d_psi = 0.0;
      e_psi = 1.0;

    } 

  }

  func =   1.0/4.0*a_psi*std::pow(y,Scalar(4.0))
         + 1.0/3.0*b_psi*std::pow(y,Scalar(3.0))
         + 1.0/2.0*c_psi*std::pow(y,Scalar(2.0))
         + d_psi*y
         + e_psi;

  return func;
}
 
template <typename Scalar, typename Scalar2>
Scalar helper_phi(Scalar2 kx, Scalar2 kz, Scalar y, Scalar2 noSlip, 
                                                    Scalar2 addMean)
{
  Scalar func;

  Scalar2 kMag = std::sqrt(kx*kx + kz*kz);

  Scalar2 a_phi, b_phi, c_phi, d_phi;

  if ( kMag == 0. ) {

    a_phi = 1.0;
    b_phi = -10.0/3.0;
    c_phi = -6.0*a_phi;
    d_phi = 0.0;

    if ( addMean == 0.0 ) return 0.0;

  } else {

    if ( noSlip == 0.0 ) {
      a_phi = -(2.0 + kMag)/(4.0 + kMag);
      b_phi = -(8.0 + 3.0*kMag)/(3.0 + kMag) - a_phi*(kMag + 4.0)/(kMag + 3.0);
      c_phi = 1.0;
      d_phi = 1.0;

    } else {

      a_phi = 0.0;
      b_phi = -2.0;
      c_phi = 0.0;
      d_phi = 1.0;

    }

  }

  func =  std::pow(y,Scalar(5.0))
         + a_phi*std::pow(y,Scalar(4.0))
         + b_phi*std::pow(y,Scalar(3.0))
         + c_phi*std::pow(y,Scalar(2.0))
         + d_phi*y;

  return func;
} 

template <typename Scalar>
Scalar helper_alpha(Scalar y)
{
  Scalar func;
  // make sure drho/dy = 0 at the boundary
  // TODO doesn't necessary force U (unless all mean...)
  func = 2.0 + std::pow(y,Scalar(2.0)) - .5*std::pow(y,Scalar(4.0));

  return func;
}

template <typename Scalar>
Scalar helper_beta(Scalar y)
{
  Scalar func;
  func = 2.0 - std::pow(y,Scalar(2.0));

  return func;
}

template <typename Scalar>
Scalar helper_gamma(Scalar y)
{
  Scalar func;

  func = std::sqrt(Scalar(2.0))/16.0*(std::pow(y,Scalar(3.0)) + 1.0);

  return func;
}

template <typename Scalar>
Scalar helper_delta(Scalar y)
{

  Scalar func;

  func = 1. - std::pow(y,Scalar(4.0));

  return func;

}

template <typename Scalar>
Scalar helper_T(Scalar t)
{
  Scalar func;

  func = 1.;
  //func = (1. + 10.*t + 20.*t*t + 30.*t*t*t); // + 1.e13*t*t + 0.*30.*t*t*t);//std::exp(t);
  //func = std::exp(t/1e-4);
  //func = std::exp(t);
  return func;
}

template <typename Scalar>
Scalar rhoMapVD(Scalar Z, Scalar at)
{
  Scalar rho;

  rho = ( 1.0 + at ) / ( 1.0 + at - 2.0*at*Z );

  return rho;

}

template <typename Scalar>
Scalar drhoMapVD(Scalar Z, Scalar at)
{
  Scalar drho;

  drho = ( 1.0 + at ) / ( 1.0 + at - 2.0*at*Z );

  return drho;

}

//
// main functions
//
// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_mD_1(Scalar x1, 
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef D1Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  Scalar exact_1;
  exact_1 = helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                           noSlip,addMean,x,y,z).derivatives()[1];

  return exact_1 * helper_T(t1);
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_mD_2(Scalar x1, 
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef D1Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  Scalar exact_2;
  exact_2 = ( helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                             noSlip,addMean,x,y,z).derivatives()[2]
            - helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                             noSlip,addMean,x,y,z).derivatives()[0] );

  return exact_2 * helper_T(t1);
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_mD_3(Scalar x1, 
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef D1Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  Scalar exact_3;

  exact_3 = -helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                           noSlip,addMean,x,y,z).derivatives()[1];

  return exact_3 * helper_T(t1);
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_mC_1(Scalar x1, 
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef D1Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  Scalar exact_1;

  exact_1 = helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                           noSlip,x,y,z).derivatives()[0];

  return exact_1 * helper_T(t1);
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_mC_2(Scalar x1, 
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef D1Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  Scalar exact_2;

  exact_2 = helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                           noSlip,x,y,z).derivatives()[1];

  return exact_2 * helper_T(t1);
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_mean_mC_2(Scalar x1, 
                                                                           Scalar y1, 
                                                                           Scalar z1,
                                                                           Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef D1Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  Scalar exact_2;

  exact_2 = helper_zetaPsi(zero,zero,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,one,
                           noSlip,x,y,z).derivatives()[1];

  return exact_2 * helper_T(t1);
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_mC_3(Scalar x1, 
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef D1Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  Scalar exact_3;

  exact_3 = helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                           noSlip,x,y,z).derivatives()[2];

  return exact_3 * helper_T(t1);
}

// public method

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_rho(Scalar x1,
                                                                     Scalar y1, 
                                                                     Scalar z1,
                                                                     Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef D1Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  Scalar exact_rho;

  exact_rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                 x1,y1,z1) * helper_T(t1);

  return exact_rho;
}

// public method

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_p(Scalar x1,
                                                                   Scalar y1, 
                                                                   Scalar z1,
                                                                   Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef D1Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  Scalar exact_p;

  exact_p = helper_zetaDelta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                             x1,y1,z1) * helper_T(t1);

  return exact_p;
}

// public method

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_drho(Scalar x1,
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef D1Type ADScalar;

  typedef DualNumber<Scalar, NumberVector<NDIM+1, Scalar> > D1TimeType;
  typedef D1TimeType ADTimeScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  ADTimeScalar xt = 
    ADTimeScalar(x1,NumberVectorUnitVector<NDIM+1, 0, Scalar>::value());
  ADTimeScalar yt = 
    ADTimeScalar(y1,NumberVectorUnitVector<NDIM+1, 1, Scalar>::value());
  ADTimeScalar zt = 
    ADTimeScalar(z1,NumberVectorUnitVector<NDIM+1, 2, Scalar>::value());
  ADTimeScalar tt = 
    ADTimeScalar(t1,NumberVectorUnitVector<NDIM+1, 3, Scalar>::value());

  ADTimeScalar exact_rho;

  exact_rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                 x1,y1,z1) * helper_T(tt);

  return exact_rho.derivatives()[3];
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_z(Scalar x1,
                                                                   Scalar y1, 
                                                                   Scalar z1,
                                                                   Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef D1Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  Scalar exact_z;

  exact_z = helper_zetaGammaPlus(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                 x1,y1,z1) * helper_T(t1);

  return exact_z;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_mu(Scalar x1,
                                                                    Scalar y1, 
                                                                    Scalar z1,
                                                                    Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef D1Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  Scalar exact_mu;

  exact_mu = helper_beta(y1); //helper_zeta(kx,kz,x1,z1) * helper_beta(y1);// * helper_T(t1);

  return exact_mu;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_omega(Scalar x1, 
                                                                       Scalar y1, 
                                                                       Scalar z1, 
                                                                       Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef D2Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // Omega = d m_1^d / dz - d m_3^d / dx

  ADScalar zetaPhi = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,addMean,x,y,z);

                  // d m_1^d / dz
  Scalar omega =  (zetaPhi.derivatives()[1]).derivatives()[2] 
                                                      *helper_T(t1)
                - // d m_3^d / dx
                 -(zetaPhi.derivatives()[1]).derivatives()[0]
                                                      *helper_T(t1);

  return omega;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_phi(Scalar x1,
                                                                     Scalar y1, 
                                                                     Scalar z1, 
                                                                     Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef DualNumber<D3Type, NumberVector<NDIM, D3Type> > D4Type;
  typedef D3Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  // phi = lap( m_2^d ) 
 
  ADScalar m2_D = 
    ( helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                     noSlip,addMean,x,y,z).derivatives()[2] 
    - 
      helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                     noSlip,addMean,x,y,z).derivatives()[0] ) * helper_T(t1);

  ADScalar phi = divergence(gradient(m2_D)); 

  return raw_value(phi);

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_div_mC(Scalar x1,
                                                                        Scalar y1, 
                                                                        Scalar z1, 
                                                                        Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef D2Type ADScalar;

  // Treat momentum as a vector
  NumberVector<NDIM, D1Type> mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mC = gradient(zetaPsi)*helper_T(t1);

  return divergence(mC);

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_RHSomega(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef DualNumber<D3Type, NumberVector<NDIM, D3Type> > D4Type;
  typedef D4Type ADScalar;

  // Treat velocity, momentum as a vector
  NumberVector<NDIM, D3Type> U,mD,mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPhi = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,addMean,x,y,z);
  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mD[0] =    zetaPhi.derivatives()[1]   * helper_T(t1);
  mD[1] = (  zetaPhi.derivatives()[2] 
           - zetaPhi.derivatives()[0] ) * helper_T(t1);
  mD[2] =  - zetaPhi.derivatives()[1]   * helper_T(t1);

  mC = gradient(zetaPsi)*helper_T(t1);

  D3Type x3  = D3Type(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  D3Type y3  = D3Type(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  D3Type z3  = D3Type(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  D3Type rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                 x3,y3,z3) * helper_T(t1); 
  D3Type mu  = helper_beta(y3); //helper_zeta(kx,kz,x,z) * helper_beta(y);// * helper_T(t1); 

  U = (mD + mC) / rho;

  // get DivCij

  NumberVector<NDIM, D2Type > DivC = -divergence(rho*U.outerproduct(U));

  // get DivTauij

  NumberVector<NDIM, NumberVector<NDIM, Scalar> > I = 
    NumberVector<NDIM, Scalar>::identity();
  NumberVector<NDIM, D2Type> DivTau = 
    divergence(
        mu * (gradient(U) + transpose(gradient(U)) 
              - 2./3. * divergence(U)*I) );

  Scalar RHS_omega =
    raw_value(
               // Convective part
                DivC[0].derivatives()[2]
              - DivC[2].derivatives()[0]
               // Viscous part
              + 1./re*DivTau[0].derivatives()[2]
              - 1./re*DivTau[2].derivatives()[0]

             );

  return RHS_omega;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_RHSphi(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef DualNumber<D3Type, NumberVector<NDIM, D3Type> > D4Type;
  typedef DualNumber<D4Type, NumberVector<NDIM, D4Type> > D5Type;
  typedef D5Type ADScalar;

  // Treat velocity, momentum as a vector
  NumberVector<NDIM, D4Type> U,mD,mC;
 
  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPhi = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,addMean,x,y,z);
  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mD[0] =    zetaPhi.derivatives()[1]   * helper_T(t1);
  mD[1] = (  zetaPhi.derivatives()[2] 
           - zetaPhi.derivatives()[0] ) * helper_T(t1);
  mD[2] =  - zetaPhi.derivatives()[1]   * helper_T(t1);

  mC = gradient(zetaPsi)*helper_T(t1);

  D4Type x4  = D4Type(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  D4Type y4  = D4Type(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  D4Type z4  = D4Type(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  D4Type rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                 x4,y4,z4) * helper_T(t1); 

  D4Type mu  = helper_beta(y4); //helper_zeta(kx,kz,x,z) * helper_beta(y);// * helper_T(t1); 

  U = (mD + mC) / rho;

  // get DivCij

  NumberVector<NDIM, D3Type > DivC = -divergence(rho*U.outerproduct(U));

  // get DivTauij

  NumberVector<NDIM, NumberVector<NDIM, Scalar> > I = 
    NumberVector<NDIM, Scalar>::identity();
  NumberVector<NDIM, D3Type> DivTau = 
    divergence(
        mu * (gradient(U) + transpose(gradient(U)) 
              - 2./3. * divergence(U)*I) );

  // get divergences
  
  D2Type DivDivC   = divergence(DivC);
  D2Type DivDivTau = divergence(DivTau);

  // phi RHS
  Scalar RHS_phi = 
    raw_value(
                // Convective part
                divergence(gradient(DivC[1]))
              - DivDivC.derivatives()[1]
              +
                // Viscous part
                1./re*divergence(gradient(DivTau[1]))
              - 1./re*DivDivTau.derivatives()[1]

                // gravity
              - ( rho.derivatives()[0].derivatives()[0] + 
                  rho.derivatives()[2].derivatives()[2] )

             );

  return RHS_phi;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_g_DivTau(Scalar x1, Scalar y1, Scalar z1, Scalar t1, int i)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef DualNumber<D3Type, NumberVector<NDIM, D3Type> > D4Type;
  typedef D3Type ADScalar;

  // Treat velocity, momentum as a vector
  NumberVector<NDIM, D2Type> U,mD,mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPhi = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,addMean,x,y,z);
  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mD[0] =    zetaPhi.derivatives()[1]   * helper_T(t1);
  mD[1] = (  zetaPhi.derivatives()[2] 
           - zetaPhi.derivatives()[0] ) * helper_T(t1);
  mD[2] =  - zetaPhi.derivatives()[1]   * helper_T(t1);

  mC = gradient(zetaPsi)*helper_T(t1);

  D2Type x2 = D2Type(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  D2Type y2 = D2Type(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  D2Type z2 = D2Type(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  D2Type rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                 x2,y2,z2) * helper_T(t1); 
  D2Type mu  = helper_beta(y2); //helper_zeta(kx,kz,x,z) * helper_beta(y);// * helper_T(t1); 

  U = (mD + mC) / rho;

  // get DivTauij

  NumberVector<NDIM, NumberVector<NDIM, Scalar> > I = 
    NumberVector<NDIM, Scalar>::identity();
  NumberVector<NDIM, D1Type> DivTau = 
    divergence(
        mu * (gradient(U) + transpose(gradient(U)) 
              - 2./3. * divergence(U)*I) );

  return raw_value(DivTau[i]);
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_g_C(Scalar x1, Scalar y1, Scalar z1, Scalar t1, int i)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef DualNumber<D3Type, NumberVector<NDIM, D3Type> > D4Type;
  typedef D1Type ADScalar;

  // Treat velocity, momentum as a vector
  NumberVector<NDIM, ADScalar> U,mD,mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPhi = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,addMean,x,y,z);
  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mD[0] =    zetaPhi.derivatives()[1]   * helper_T(t1);
  mD[1] = (  zetaPhi.derivatives()[2] 
           - zetaPhi.derivatives()[0] ) * helper_T(t1);
  mD[2] =  - zetaPhi.derivatives()[1]   * helper_T(t1);

  mC = gradient(zetaPsi)*helper_T(t1);

  ADScalar rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                  x,y,z) * helper_T(t1); 
 
  U = (mD + mC) / rho;

  // get Cij
  
  NumberVector<NDIM, NumberVector<NDIM,ADScalar> > C = -rho*U.outerproduct(U);

  Scalar Cval = 0.0;

  switch (i)
  {
    case 1:
      Cval = raw_value(C[0][0] - C[1][1]);
      break;
    case 2:
      Cval = raw_value(C[2][2] - C[1][1]);
      break;
    case 3:
      Cval = raw_value(C[0][1]);
      break;
    case 4:
      Cval = raw_value(C[0][2]);
      break;
    case 5:
      Cval = raw_value(C[1][2]);
      break;
  }

  return Cval;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_g_Tau(Scalar x1, Scalar y1, Scalar z1, Scalar t1, int i, int j)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef DualNumber<D3Type, NumberVector<NDIM, D3Type> > D4Type;
  typedef D3Type ADScalar;

  // Treat velocity, momentum as a vector
  NumberVector<NDIM, D2Type> U,mD,mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPhi = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,addMean,x,y,z);
  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mD[0] =    zetaPhi.derivatives()[1]   * helper_T(t1);
  mD[1] = (  zetaPhi.derivatives()[2] 
           - zetaPhi.derivatives()[0] ) * helper_T(t1);
  mD[2] =  - zetaPhi.derivatives()[1]   * helper_T(t1);

  mC = gradient(zetaPsi)*helper_T(t1);

  D2Type x2 = D2Type(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  D2Type y2 = D2Type(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  D2Type z2 = D2Type(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  D2Type rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                 x2,y2,z2) * helper_T(t1); 
  D2Type mu  = helper_beta(y2); //helper_zeta(kx,kz,x,z) * helper_beta(y);// * helper_T(t1); 

  U = (mD + mC) / rho;

  // get Tauij

  NumberVector<NDIM, NumberVector<NDIM, Scalar> > I = 
    NumberVector<NDIM, Scalar>::identity();
  NumberVector<NDIM, NumberVector<NDIM,D2Type> > Tau = 
        mu * (gradient(U) + transpose(gradient(U)) 
              - 2./3. * divergence(U)*I) ;

  return raw_value(Tau[i][j]);

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_g_u(Scalar x1, Scalar y1, Scalar z1, Scalar t1, int i)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef DualNumber<D3Type, NumberVector<NDIM, D3Type> > D4Type;
  typedef D3Type ADScalar;

  // Treat velocity, momentum as a vector
  NumberVector<NDIM, D2Type> U,mD,mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPhi = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,addMean,x,y,z);
  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mD[0] =    zetaPhi.derivatives()[1]   * helper_T(t1);
  mD[1] = (  zetaPhi.derivatives()[2] 
           - zetaPhi.derivatives()[0] ) * helper_T(t1);
  mD[2] =  - zetaPhi.derivatives()[1]   * helper_T(t1);

  mC = gradient(zetaPsi)*helper_T(t1);

  D2Type x2 = D2Type(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  D2Type y2 = D2Type(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  D2Type z2 = D2Type(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  D2Type rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                 x2,y2,z2) * helper_T(t1); 
  D2Type mu  = helper_beta(y2); //helper_zeta(kx,kz,x,z) * helper_beta(y);// * helper_T(t1); 

  U = (mD + mC) / rho;

  // get grad(u)

  NumberVector<NDIM, D2Type> du = gradient(U[0]);

  return raw_value(du[i]);

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_g_v(Scalar x1, Scalar y1, Scalar z1, Scalar t1, int i)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef DualNumber<D3Type, NumberVector<NDIM, D3Type> > D4Type;
  typedef D3Type ADScalar;

  // Treat velocity, momentum as a vector
  NumberVector<NDIM, D2Type> U,mD,mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPhi = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,addMean,x,y,z);
  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mD[0] =    zetaPhi.derivatives()[1]   * helper_T(t1);
  mD[1] = (  zetaPhi.derivatives()[2] 
           - zetaPhi.derivatives()[0] ) * helper_T(t1);
  mD[2] =  - zetaPhi.derivatives()[1]   * helper_T(t1);

  mC = gradient(zetaPsi)*helper_T(t1);

  D2Type x2 = D2Type(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  D2Type y2 = D2Type(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  D2Type z2 = D2Type(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  D2Type rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                 x2,y2,z2) * helper_T(t1); 
  D2Type mu  = helper_beta(y2); //helper_zeta(kx,kz,x,z) * helper_beta(y);// * helper_T(t1); 

  U = (mD + mC) / rho;

  // get grad(U)

  NumberVector<NDIM, D2Type> dv = gradient(U[1]);

  return raw_value(dv[i]);

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_g_w(Scalar x1, Scalar y1, Scalar z1, Scalar t1, int i)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef DualNumber<D3Type, NumberVector<NDIM, D3Type> > D4Type;
  typedef D3Type ADScalar;

  // Treat velocity, momentum as a vector
  NumberVector<NDIM, D2Type> U,mD,mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPhi = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,addMean,x,y,z);
  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mD[0] =    zetaPhi.derivatives()[1]   * helper_T(t1);
  mD[1] = (  zetaPhi.derivatives()[2] 
           - zetaPhi.derivatives()[0] ) * helper_T(t1);
  mD[2] =  - zetaPhi.derivatives()[1]   * helper_T(t1);

  mC = gradient(zetaPsi)*helper_T(t1);

  D2Type x2 = D2Type(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  D2Type y2 = D2Type(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  D2Type z2 = D2Type(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  D2Type rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                 x2,y2,z2) * helper_T(t1); 
  D2Type mu  = helper_beta(y2); //helper_zeta(kx,kz,x,z) * helper_beta(y);// * helper_T(t1); 

  U = (mD + mC) / rho;

  // get grad(w)

  NumberVector<NDIM, D2Type> dw = gradient(U[2]);

  return raw_value(dw[i]);

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_g_Z(Scalar x1, Scalar y1, Scalar z1, Scalar t1, int i)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef D1Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  NumberVector<NDIM, D1Type> dZ = gradient(
      helper_zetaGammaPlus(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
        x,y,z) ) * helper_T(t1);

  return raw_value(dZ[i]);

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_RHSz(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef D2Type ADScalar;

  // Treat velocity, momentum as a vector
  NumberVector<NDIM, ADScalar> U,mD,mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zetaPhi = 
    helper_zetaPhi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,addMean,x,y,z);
  ADScalar zetaPsi = 
    helper_zetaPsi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                   noSlip,x,y,z);

  mD[0] =    zetaPhi.derivatives()[1]   * helper_T(t1);
  mD[1] = (  zetaPhi.derivatives()[2] 
           - zetaPhi.derivatives()[0] ) * helper_T(t1);
  mD[2] =  - zetaPhi.derivatives()[1]   * helper_T(t1);

  mC = gradient(zetaPsi)*helper_T(t1);

  ADScalar rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                  x,y,z) * helper_T(t1); 
 
  ADScalar mu  = helper_beta(y); //helper_zeta(kx,kz,x,z) * helper_beta(y);// * helper_T(t1); 

  // set mC to zero for new method
  
  U = (mD + 0.*mC) / rho;
 
  ADScalar zVar = helper_zetaGammaPlus(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                       x,y,z) * helper_T(t1);

  Scalar RHS_z = raw_value(-U.dot(gradient(zVar)))
               + raw_value(1.0/(re*sc)*(1.0/rho*gradient(mu).dot(gradient(zVar))
                   + mu/rho*divergence(gradient(zVar))));

  return RHS_z;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_g_jacF(Scalar x1, Scalar y1, Scalar z1, Scalar t1, int i)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef DualNumber<D3Type, NumberVector<NDIM, D3Type> > D4Type;
  typedef D3Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zVar = helper_zetaGammaPlus(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                       x,y,z);

  D2Type x2 = D2Type(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  D2Type y2 = D2Type(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  D2Type z2 = D2Type(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  D2Type rho = helper_zetaAlpha(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,
                                  x2,y2,z2) * helper_T(t1); 

  return raw_value(-1./rho*gradient(zVar)[i]);
}

// ----------------------------------------
// Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::navierstokes_3d_variabledensity);


#endif // HAVE_METAPHYSICL
