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
Scalar helper_zeta(Scalar2 kx1, Scalar2 kz1, 
                   Scalar2 kx2, Scalar2 kz2, 
                   Scalar2 kx3, Scalar2 kz3, 
                   Scalar2 kx4, Scalar2 kz4, 
                   Scalar2 kx5, Scalar2 kz5,
                   Scalar2 numModes, 
                   Scalar x, Scalar z);

template <typename Scalar, typename Scalar2>
Scalar helper_meanOnly(Scalar2 kx, Scalar2 kz, Scalar x, Scalar z);

template <typename Scalar, typename Scalar2>
Scalar helper_psi(Scalar2 kx1, Scalar2 kz1, 
                  Scalar2 kx2, Scalar2 kz2, 
                  Scalar2 kx3, Scalar2 kz3, 
                  Scalar2 kx4, Scalar2 kz4, 
                  Scalar2 kx5, Scalar2 kz5,
                  Scalar2 numModes, 
                  Scalar y);
  
template <typename Scalar, typename Scalar2>
Scalar helper_phi(Scalar2 kx1, Scalar2 kz1, 
                  Scalar2 kx2, Scalar2 kz2, 
                  Scalar2 kx3, Scalar2 kz3, 
                  Scalar2 kx4, Scalar2 kz4, 
                  Scalar2 kx5, Scalar2 kz5,
                  Scalar2 numModes, 
                  Scalar y);

template <typename Scalar>
Scalar helper_alpha(Scalar y);

template <typename Scalar>
Scalar helper_beta(Scalar y);

template <typename Scalar>
Scalar helper_T(Scalar t);

//template <typename Scalar>
//NumberVector DivTauij(Scalar x1, Scalar y1, Scalar z1, Scalar t)

typedef ShadowNumber<double, long double> RawScalar;
const unsigned int NDIM = 3;  // Need time dependence!! But put in manually!
typedef DualNumber<RawScalar, NumberVector<NDIM, RawScalar> > FirstDerivType;
typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
typedef SecondDerivType ADType;

using namespace MASA;

template <typename Scalar>
MASA::navierstokes_3d_variabledensity<Scalar>::navierstokes_3d_variabledensity()
{
  this->mmsname = "navierstokes_3d_variabledensity";
  this->dimension = 4;

  this->register_var("re",&re);
  this->register_var("sc",&sc);
  this->register_var("numModes",&numModes);
  this->register_var("kx1",&kx1);
  this->register_var("kx2",&kx2);
  this->register_var("kx3",&kx3);
  this->register_var("kx4",&kx4);
  this->register_var("kx5",&kx5);
  this->register_var("kz1",&kz1);
  this->register_var("kz2",&kz2);
  this->register_var("kz3",&kz3);
  this->register_var("kz4",&kz4);
  this->register_var("kz5",&kz5);
  this->register_var("kMag",&kMag);

  this->init_var();

} // done with constructor

template <typename Scalar>
int MASA::navierstokes_3d_variabledensity<Scalar>::init_var()
{
  int err = 0;

  double kx,kz,kMag;

  kx = 1.0;
  kz = 1.0;
  kMag = std::sqrt(kx*kx + kz*kz);

  err += this->set_var("re",1);
  err += this->set_var("sc",1);
  err += this->set_var("numModes",1);
  err += this->set_var("kx1",kx);
  err += this->set_var("kx2",kx);
  err += this->set_var("kx3",kx);
  err += this->set_var("kx4",kx);
  err += this->set_var("kx5",kx);
  err += this->set_var("kz1",kx);
  err += this->set_var("kz2",kz); 
  err += this->set_var("kz3",kz); 
  err += this->set_var("kz4",kz); 
  err += this->set_var("kz5",kz); 
  err += this->set_var("kMag",kMag); 

  return err;

} // done with init_var

// ----------------------------------------
// Source Terms
// ----------------------------------------

// public static method, that can be called from eval_q_t
template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_q_omega(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  typedef DualNumber<Scalar, NumberVector<NDIM+1, Scalar> > FirstDerivTimeType;
  typedef DualNumber<FirstDerivTimeType, NumberVector<NDIM+1, FirstDerivTimeType> > SecondDerivTimeType;
  typedef DualNumber<SecondDerivTimeType, NumberVector<NDIM+1, SecondDerivTimeType> > ThirdDerivTimeType;
  typedef ThirdDerivTimeType ADTimeScalar;


  // Treat velocity, momentum as a vector
  NumberVector<NDIM, ADScalar> U,mD,mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  //ADScalar t = ADScalar(t1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());

  mD[0] =  helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
           helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1);
  mD[1] = (  helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2] 
           - helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0] ) * 
          helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1);
  mD[2] = -helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
           helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1);

  mC[0] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0] *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1); 
  mC[1] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1); 
  mC[2] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2] *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1); 

  ADScalar rho = helper_alpha(y);//*helper_meanOnly(kx,kz,x1,z1);// * helper_alpha(y);// * helper_T(t1); 
  ADScalar mu  = helper_beta(y);//*helper_meanOnly(kx,kz,x1,z1);// * helper_beta(y);// * helper_T(t1); 

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

  // Omega = d m_1^d / dz - d m_3^d / dx

  ADTimeScalar omega = 
        // d m_1^d / dz
    ( helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,xt,zt).derivatives()[2] *
           helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,yt).derivatives()[1] * helper_T(tt) )
    -   // d m_3^d / dz
    (-helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,xt,zt).derivatives()[0] *
           helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,yt).derivatives()[1] * helper_T(tt) );

  // get DivCij

  NumberVector<NDIM, ADScalar> DivC = -divergence(rho*U.outerproduct(U));

  // get DivTauij

  NumberVector<NDIM, NumberVector<NDIM, Scalar> > I = 
    NumberVector<NDIM, Scalar>::identity();
  NumberVector<NDIM, ADScalar> DivTau = 
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
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  typedef DualNumber<Scalar, NumberVector<NDIM+1, Scalar> > FirstDerivTimeType;
  typedef DualNumber<FirstDerivTimeType, NumberVector<NDIM+1, FirstDerivTimeType> > SecondDerivTimeType;
  typedef DualNumber<SecondDerivTimeType, NumberVector<NDIM+1, SecondDerivTimeType> > ThirdDerivTimeType;
  typedef ThirdDerivTimeType ADTimeScalar;


  // Treat velocity, momentum as a vector
  NumberVector<NDIM, ADScalar> U,mD,mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  //ADScalar t = ADScalar(t1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());

  mD[0] =  helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
           helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1);
  mD[1] = (  helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2] 
           - helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0] ) * 
          helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1);
  mD[2] = -helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
           helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1);

  mC[0] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0] *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1); 
  mC[1] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1); 
  mC[2] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2] *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1); 

  ADScalar rho = helper_alpha(y);//*helper_meanOnly(kx,kz,x1,z1);// * helper_alpha(y);// * helper_T(t1); 
  ADScalar mu  = helper_beta(y);//*helper_meanOnly(kx,kz,x1,z1);// * helper_beta(y);// * helper_T(t1); 

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

  ADTimeScalar phi = 
        // d^2 m_2^d / dx^2
    (  ((helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,xt,zt).derivatives()[2])
                                 .derivatives()[0]).derivatives()[0] 
      -((helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,xt,zt).derivatives()[0])
                                 .derivatives()[0]).derivatives()[0]) * 
          helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,yt) * helper_T(tt) +
        // d^2 m_2^d / dy^2
    (   helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,xt,zt).derivatives()[2]
           - helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,xt,zt).derivatives()[0] ) * 
       (helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,yt).derivatives()[1]).derivatives()[1] * 
        helper_T(tt) +
        // d^2 m_2^d / dz^2
    (  ((helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,xt,zt).derivatives()[2])
                                 .derivatives()[2]).derivatives()[2] 
      -((helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,xt,zt).derivatives()[0])
                                 .derivatives()[2]).derivatives()[2]) * 
          helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,yt) * helper_T(tt);

  // get DivCij

  NumberVector<NDIM, ADScalar > DivC = -divergence(rho*U.outerproduct(U));

  // get DivTauij

  NumberVector<NDIM, NumberVector<NDIM, Scalar> > I = 
    NumberVector<NDIM, Scalar>::identity();
  NumberVector<NDIM, ADScalar> DivTau = 
    divergence(
        mu * (gradient(U) + transpose(gradient(U)) 
              - 2./3. * divergence(U)*I) );

  // get divergences

  ADScalar DivDivC   = divergence(DivC);
  ADScalar DivDivTau = divergence(DivTau);

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

              );

  return Q_phi;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_q_rho(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  typedef DualNumber<Scalar, NumberVector<NDIM+1, Scalar> > FirstDerivTimeType;
  typedef DualNumber<FirstDerivTimeType, NumberVector<NDIM+1, FirstDerivTimeType> > SecondDerivTimeType;
  typedef DualNumber<SecondDerivTimeType, NumberVector<NDIM+1, SecondDerivTimeType> > ThirdDerivTimeType;
  typedef ThirdDerivTimeType ADTimeScalar;


  // Treat momentum as a vector
  NumberVector<NDIM, ADScalar> mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  //ADScalar t = ADScalar(t1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());

  mC[0] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0] *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1); 
  mC[1] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1); 
  mC[2] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2] *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1); 

  // Time part

  ADTimeScalar xt = 
    ADTimeScalar(x1,NumberVectorUnitVector<NDIM+1, 0, Scalar>::value());
  ADTimeScalar yt = 
    ADTimeScalar(y1,NumberVectorUnitVector<NDIM+1, 1, Scalar>::value());
  ADTimeScalar zt = 
    ADTimeScalar(z1,NumberVectorUnitVector<NDIM+1, 2, Scalar>::value());
  ADTimeScalar tt = 
    ADTimeScalar(t1,NumberVectorUnitVector<NDIM+1, 3, Scalar>::value());

  ADTimeScalar rho = helper_alpha(yt);//*helper_meanOnly(kx,kz,x1,z1);//* helper_alpha(yt);// * helper_T(tt); 

  // Rho equation residuals
  Scalar Q_rho = 
    raw_value(  
                // Temporal part
                rho.derivatives()[3]

        ) +
    raw_value(
                // Div mC
                divergence(mC)

              );

  return Q_rho;
}
//
//
//// public static method, that can be called from eval_q_t
//template <typename Scalar>
//Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_q_v(Scalar x1, Scalar y1, Scalar z1)
//{
//  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
//  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
//  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
//  typedef ThirdDerivType ADScalar;
//
//  // Treat velocity as a vector
//  NumberVector<NDIM, ADScalar> U;
//
//  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
//  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
//  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
//
//  // Arbitrary manufactured solutions
//  U[0]       = a * helper_f(beta,kx,x)                  * helper_g(y).derivatives()[1] * helper_h(gamma,kz,z).derivatives()[2];
//  U[1]       = b * helper_f(beta,kx,x).derivatives()[0] * helper_g(y)                  * helper_h(gamma,kz,z).derivatives()[2];
//  U[2]       = c * helper_f(beta,kx,x).derivatives()[0] * helper_g(y).derivatives()[1] * helper_h(gamma,kz,z);
//  ADScalar P = d * helper_f(beta,kx,x)                  * helper_gt(y)                 * helper_h(gamma,kz,z);
//
//  // NS equation residuals
//  NumberVector<NDIM, Scalar> Q_rho_u = 
//    raw_value(
//
//	      // convective term
//	      - divergence(U.outerproduct(U))
//
//	      // pressure
//	      - P.derivatives()
//
//	      + nu * divergence(gradient(U)));
//
//  return -Q_rho_u[1];
//}
//
//
//// public static method, that can be called from eval_q_t
//template <typename Scalar>
//Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_q_w(Scalar x1, Scalar y1, Scalar z1)
//{
//  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
//  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
//  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
//  typedef ThirdDerivType ADScalar;
//
//  // Treat velocity as a vector
//  NumberVector<NDIM, ADScalar> U;
//
//  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
//  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
//  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
//
//  // Arbitrary manufactured solutions
//  U[0]       = a * helper_f(beta,kx,x)                  * helper_g(y).derivatives()[1] * helper_h(gamma,kz,z).derivatives()[2];
//  U[1]       = b * helper_f(beta,kx,x).derivatives()[0] * helper_g(y)                  * helper_h(gamma,kz,z).derivatives()[2];
//  U[2]       = c * helper_f(beta,kx,x).derivatives()[0] * helper_g(y).derivatives()[1] * helper_h(gamma,kz,z);
//  ADScalar P = d * helper_f(beta,kx,x)                  * helper_gt(y)                 * helper_h(gamma,kz,z);
//
//  // NS equation residuals
//  NumberVector<NDIM, Scalar> Q_rho_u = 
//    raw_value(
//
//	      // convective term
//	      - divergence(U.outerproduct(U))
//
//	      // pressure
//	      - P.derivatives()
//
//	      // dissipation
//	      + nu * divergence(gradient(U)));
//
//  return -Q_rho_u[2];
//}

// ----------------------------------------
// Analytical Terms
// ----------------------------------------

// example of a private method, called from exact_t
template <typename Scalar, typename Scalar2>
Scalar helper_zeta(Scalar2 kx1, Scalar2 kz1, 
                   Scalar2 kx2, Scalar2 kz2, 
                   Scalar2 kx3, Scalar2 kz3, 
                   Scalar2 kx4, Scalar2 kz4, 
                   Scalar2 kx5, Scalar2 kz5,
                   Scalar2 numModes, 
                   Scalar x, Scalar z)
{
  Scalar func = 0.0;
  
  Scalar2 kx,kz;

  int modes = int(numModes);

  for (int i = 1; i <= modes; i++) {

    switch(i) {
      case 1:
        kx = kx1;
        kz = kz1;
        break;
      case 2:
        kx = kx2;
        kz = kz2;
        break;
      case 3:
        kx = kx3;
        kz = kz3;
        break;
      case 4:
        kx = kx4;
        kz = kz4;
        break;
      case 5:
        kx = kx5;
        kz = kz5;
        break;
    } 

  func += 2*( std::cos(kx*x)*std::cos(kz*z) - std::sin(kx*x)*std::sin(kz*z)
             -std::sin(kx*x)*std::cos(kz*z) - std::cos(kx*x)*std::sin(kz*z) );
  }

  return func;
}

template <typename Scalar, typename Scalar2>
Scalar helper_psi(Scalar2 kx1, Scalar2 kz1, 
                  Scalar2 kx2, Scalar2 kz2, 
                  Scalar2 kx3, Scalar2 kz3, 
                  Scalar2 kx4, Scalar2 kz4, 
                  Scalar2 kx5, Scalar2 kz5,
                  Scalar2 numModes, 
                  Scalar y)
{
  Scalar func = 0.0;

  Scalar2 kMag, kx, kz;

  Scalar2 a_psi, c_psi, e_psi;
  
  int modes = int(numModes);

  for (int i = 1; i <= modes; i++) {

    switch(i) {
      case 1:
        kx = kx1;
        kz = kz1;
        break;
      case 2:
        kx = kx2;
        kz = kz2;
        break;
      case 3:
        kx = kx3;
        kz = kz3;
        break;
      case 4:
        kx = kx4;
        kz = kz4;
        break;
      case 5:
        kx = kx5;
        kz = kz5;
        break;
    } 

    kMag = sqrt(kx*kx+kz*kz);
 
    if ( kMag == 0. ) {
      
      a_psi = -1.0;
      c_psi = 3.0;
      e_psi = 0.0;
  
    } else {
    
      a_psi = - (kMag + 1.0)/2.0;
      c_psi = 1.0 - a_psi;
      e_psi = 1.0/2.0*(-2.0/kMag - 2.0/a_psi - c_psi);
  
    }
  
    func +=  1/4*a_psi*std::pow(y,Scalar(4.0))
           + 1/2*c_psi*std::pow(y,Scalar(2.0))
           + e_psi;
  
  }

  return func;
}
 
template <typename Scalar, typename Scalar2>
Scalar helper_phi(Scalar2 kx1, Scalar2 kz1, 
                  Scalar2 kx2, Scalar2 kz2, 
                  Scalar2 kx3, Scalar2 kz3, 
                  Scalar2 kx4, Scalar2 kz4, 
                  Scalar2 kx5, Scalar2 kz5,
                  Scalar2 numModes, 
                  Scalar y)
{
  Scalar func;

  Scalar2 kx,kz,kMag;

  Scalar2 a_phi, b_phi;

  int modes = int(numModes);

  for (int i = 1; i <= modes; i++) {

    switch(i) {
      case 1:
        kx = kx1;
        kz = kz1;
        break;
      case 2:
        kx = kx2;
        kz = kz2;
        break;
      case 3:
        kx = kx3;
        kz = kz3;
        break;
      case 4:
        kx = kx4;
        kz = kz4;
        break;
      case 5:
        kx = kx5;
        kz = kz5;
        break;
    } 

    kMag = sqrt(kx*kx+kz*kz);
    if ( kMag == 0. ) {

      func += 0.0;

    } else {

      a_phi = -(2.0 + kMag)/(4.0 + kMag);
      b_phi = -(8.0 + 3.0*kMag)/(3.0 + kMag) - a_phi*(kMag + 4.0)/(kMag + 3.0);

    }

    func +=   std::pow(y,Scalar(5.0))
            + a_phi*std::pow(y,Scalar(4.0))
            + b_phi*std::pow(y,Scalar(3.0))
            + std::pow(y,Scalar(2.0))
            + y;

  }

  return func;
} 

template <typename Scalar>
Scalar helper_alpha(Scalar y)
{
  Scalar func;
  func = 2.0 - std::pow(y,Scalar(4.0));

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
Scalar helper_T(Scalar t)
{
  Scalar func;
  func = std::exp(t);

  return func;
}

template <typename Scalar, typename Scalar2>
Scalar helper_meanOnly(Scalar2 kx, Scalar2 kz, Scalar x, Scalar z)
{
  Scalar func;

  if ( (kx == 0.0) and (kz == 0.0) ) {
    func = 1.0;
  } else {
    func = 0.0;
  }

  return func;
}
//template <typename Scalar>
//Scalar helper_gt(Scalar y)
////Scalar MASA::navierstokes_3d_variabledensity<Scalar>::helper_gt(Scalar y)
//{
//  Scalar func;
//  func = std::pow(y,Scalar(5.0))+std::pow(y,Scalar(4.0))+std::pow(y,Scalar(3.0))+std::pow(y,Scalar(2.0))+y;
//  return func;
//}
//
//template <typename Scalar, typename Scalar2>
////Scalar MASA::navierstokes_3d_variabledensity<Scalar>::helper_h(Scalar z)
//Scalar helper_h(Scalar2 gamma, Scalar2 kz, Scalar z)
//{
//  Scalar func;
//  func = 1/(gamma+std::sin(kz*z));
//  return func;
//}

//
// main functions
// 
//
// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_mD_1(Scalar x1, 
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef FirstDerivType ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  //ADScalar t = ADScalar(t1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());

  Scalar exact_1;
  exact_1 = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x1,z1) * helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1];
  return exact_1 * helper_T(t1);
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_mD_2(Scalar x1, 
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef FirstDerivType ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  //ADScalar t = ADScalar(t1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());

  Scalar exact_2;
  exact_2 = (  helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2] 
             - helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0]) * helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y1);
  return exact_2 * helper_T(t1);
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_mD_3(Scalar x1, 
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef FirstDerivType ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  //ADScalar t = ADScalar(t1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());

  Scalar exact_3;
  exact_3 = -helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x1,z1) * helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1];
  return exact_3 * helper_T(t1);
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_mC_1(Scalar x1, 
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef FirstDerivType ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  //ADScalar t = ADScalar(t1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());

  Scalar exact_1;
  exact_1 = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0] 
            * helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y1);
  return exact_1 * helper_T(t1);
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_mC_2(Scalar x1, 
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef FirstDerivType ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  //ADScalar t = ADScalar(t1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());

  Scalar exact_2;
  exact_2 = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x1,z1) 
            * helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1];
  return exact_2 * helper_T(t1);
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_mC_3(Scalar x1, 
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef FirstDerivType ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  //ADScalar t = ADScalar(t1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());

  Scalar exact_3;
  exact_3 = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2]
            * helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y1);
  return exact_3 * helper_T(t1);
}

//// public method
//
//template <typename Scalar>
//Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_om(Scalar x1,
//                                                                    Scalar y1, 
//                                                                    Scalar z1,
//                                                                    Scalar t1)
//{
//  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
//  typedef FirstDerivType ADScalar;
//
//  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
//  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
//  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
//  ADScalar t = ADScalar(z1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());
//
//  NumberVector<NDIM, ADScalar> mD;
//
//  mD[0] = eval_exact_mD_1(x,y,z,t);
//  mD[1] = eval_exact_mD_2(x,y,z,t);
//  mD[2] = eval_exact_mD_3(x,y,z,t);
//
//  NumberVector<NDIM, NumberVector<NDIM, Scalar> > gradmD = 
//    raw_value(gradient(mD));
//
//  return gradmD[0][2] - gradmD[2][0];
//}
//
//// public method
//
//template <typename Scalar>
//Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_phi(Scalar x1,
//                                                                     Scalar y1, 
//                                                                     Scalar z1,
//                                                                     Scalar t1)
//{
//  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
//  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
//  typedef SecondDerivType ADScalar;
//
//  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
//  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
//  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
//  ADScalar t = ADScalar(z1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());
//
//  NumberVector<NDIM, ADScalar> mD;
//
//  mD[0] = eval_exact_mD_1(x,y,z,t);
//  mD[1] = eval_exact_mD_2(x,y,z,t);
//  mD[2] = eval_exact_mD_3(x,y,z,t);
//
//  NumberVector<NDIM, Scalar> phi = 
//    raw_value(divergence(gradient(mD)));
//
//  return phi[1];
//}

// public method

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_rho(Scalar x1,
                                                                     Scalar y1, 
                                                                     Scalar z1,
                                                                     Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef FirstDerivType ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  //ADScalar t = ADScalar(t1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());

  Scalar exact_rho;

  exact_rho = helper_alpha(y1);//*helper_meanOnly(kx,kz,x1,z1);//helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x1,z1) * helper_alpha(y1);// * helper_T(t1);

  return exact_rho;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_mu(Scalar x1,
                                                                    Scalar y1, 
                                                                    Scalar z1,
                                                                    Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef FirstDerivType ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  //ADScalar t = ADScalar(t1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());

  Scalar exact_mu;

  exact_mu = helper_beta(y1);//*helper_meanOnly(kx,kz,x1,z1); //helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x1,z1) * helper_beta(y1);// * helper_T(t1);

  return exact_mu;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_omega(Scalar x1, 
                                                                       Scalar y1, 
                                                                       Scalar z1, 
                                                                       Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  //ADScalar t = ADScalar(t1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());

  // Omega = d m_1^d / dz - d m_3^d / dx

  Scalar omega = 
        // d m_1^d / dz
    raw_value( helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2] *
           helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1) )
    -   // d m_3^d / dz
    raw_value(-helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0] *
           helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1) );

  return omega;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_phi(Scalar x1,
                                                                     Scalar y1, 
                                                                     Scalar z1, 
                                                                     Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  //ADScalar t = ADScalar(t1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());

  // phi = lap( m_2^d ) 

  ADScalar phi = 
        // d^2 m_2^d / dx^2
    (  ((helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2])
                                .derivatives()[0]).derivatives()[0] 
      -((helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0])
                                .derivatives()[0]).derivatives()[0]) * 
          helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1) +
        // d^2 m_2^d / dy^2
    (   helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2]
           - helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0] ) * 
       (helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1]).derivatives()[1] * 
        helper_T(t1) +
        // d^2 m_2^d / dz^2
    (  ((helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2])
                               .derivatives()[2]).derivatives()[2] 
      -((helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0])
                                .derivatives()[2]).derivatives()[2]) * 
          helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1);

  return raw_value(phi);

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_RHSomega(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  // Treat velocity, momentum as a vector
  NumberVector<NDIM, ADScalar> U,mD,mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  //ADScalar t = ADScalar(t1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());

  mD[0] =  helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
           helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1);
  mD[1] = (  helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2] 
           - helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0] ) * 
          helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1);
  mD[2] = -helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
           helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1);

  mC[0] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0] *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1); 
  mC[1] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1); 
  mC[2] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2] *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1); 

  ADScalar rho = helper_alpha(y);//*helper_meanOnly(kx,kz,x1,z1);// * helper_alpha(y);// * helper_T(t1); 
  ADScalar mu  = helper_beta(y);//*helper_meanOnly(kx,kz,x1,z1);// * helper_beta(y);// * helper_T(t1); 

  U = (mD + mC) / rho;

  // get DivCij

  NumberVector<NDIM, ADScalar> DivC = -divergence(rho*U.outerproduct(U));

  // get DivTauij

  NumberVector<NDIM, NumberVector<NDIM, Scalar> > I = 
    NumberVector<NDIM, Scalar>::identity();
  NumberVector<NDIM, ADScalar> DivTau = 
    divergence(
        mu * (gradient(U) + transpose(gradient(U)) 
              - 2./3. * divergence(U)*I) );

  // Omega equation residuals
  Scalar RHS_omega = 
    raw_value( 0.
                // Convective part
  //              DivC[0].derivatives()[2] 
  //            - DivC[2].derivatives()[0]
                // Viscous part
              + 1./re*DivTau[0].derivatives()[2]
              - 1./re*DivTau[2].derivatives()[0]
              
              );

  return RHS_omega;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_RHSphi(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  // Treat velocity, momentum as a vector
  NumberVector<NDIM, ADScalar> U,mD,mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  //ADScalar t = ADScalar(t1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());

  mD[0] =  helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
           helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1);
  mD[1] = (  helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2] 
           - helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0] ) * 
          helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1);
  mD[2] = -helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
           helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1);

  mC[0] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0] *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1); 
  mC[1] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1); 
  mC[2] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2] *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1); 

  ADScalar rho = helper_alpha(y);//*helper_meanOnly(kx,kz,x1,z1);// * helper_alpha(y);// * helper_T(t1); 
  ADScalar mu  = helper_beta(y);//*helper_meanOnly(kx,kz,x1,z1);// * helper_beta(y);// * helper_T(t1); 

  U = (mD + mC) / rho;

  // get DivCij

  NumberVector<NDIM, ADScalar > DivC = -divergence(rho*U.outerproduct(U));

  // get DivTauij

  NumberVector<NDIM, NumberVector<NDIM, Scalar> > I = 
    NumberVector<NDIM, Scalar>::identity();
  NumberVector<NDIM, ADScalar> DivTau = 
    divergence(
        mu * (gradient(U) + transpose(gradient(U)) 
              - 2./3. * divergence(U)*I) );

  // get divergences

  ADScalar DivDivC   = divergence(DivC);
  ADScalar DivDivTau = divergence(DivTau);

  // phi equation residual
  Scalar RHS_phi = 
    raw_value(
                // Convective part
                divergence(gradient(DivC[1]))
              - DivDivC.derivatives()[1]
              +  
                // Viscous part
                1./re*divergence(gradient(DivTau[1]))
              - 1./re*DivDivTau.derivatives()[1]

              );

  return RHS_phi;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_g_DivTau(Scalar x1, Scalar y1, Scalar z1, Scalar t1, int i)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  // Treat velocity, momentum as a vector
  NumberVector<NDIM, ADScalar> U,mD,mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  //ADScalar t = ADScalar(t1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());

  mD[0] =  helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
           helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1);
  mD[1] = (  helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2] 
           - helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0] ) * 
          helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1);
  mD[2] = -helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
           helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1);

  mC[0] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0] *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1); 
  mC[1] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1); 
  mC[2] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2] *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1); 

  ADScalar rho = helper_alpha(y);//*helper_meanOnly(kx,kz,x1,z1);// * helper_alpha(y);// * helper_T(t1); 
  ADScalar mu  = helper_beta(y);//*helper_meanOnly(kx,kz,x1,z1); //helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) * helper_beta(y);// * helper_T(t1); 

  U = (mD + mC) / rho;

  // get DivTauij

  NumberVector<NDIM, NumberVector<NDIM, Scalar> > I = 
    NumberVector<NDIM, Scalar>::identity();
  NumberVector<NDIM, ADScalar> DivTau = 
    divergence(
        mu * (gradient(U) + transpose(gradient(U)) 
              - 2./3. * divergence(U)*I) );

  return raw_value(DivTau[i]);
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_g_C(Scalar x1, Scalar y1, Scalar z1, Scalar t1, int i)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
  typedef DualNumber<SecondDerivType, NumberVector<NDIM, SecondDerivType> > ThirdDerivType;
  typedef ThirdDerivType ADScalar;

  // Treat velocity, momentum as a vector
  NumberVector<NDIM, ADScalar> U,mD,mC;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
  //ADScalar t = ADScalar(t1,NumberVectorUnitVector<NDIM, 3, Scalar>::value());

  mD[0] =  helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
           helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1);
  mD[1] = (  helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2] 
           - helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0] ) * 
          helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1);
  mD[2] = -helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
           helper_phi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1);

  mC[0] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[0] *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1); 
  mC[1] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y).derivatives()[1] * helper_T(t1); 
  mC[2] = helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z).derivatives()[2] *
          helper_psi(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,y) * helper_T(t1); 

  ADScalar rho = helper_alpha(y);//*helper_meanOnly(kx,kz,x1,z1);// * helper_alpha(y);// * helper_T(t1); 
  ADScalar mu  = helper_beta(y);//*helper_meanOnly(kx,kz,x1,z1); //helper_zeta(kx1,kz1,kx2,kz2,kx3,kz3,kx4,kz4,kx5,kz5,numModes,x,z) * helper_beta(y);// * helper_T(t1); 

  U = (mD + mC) / rho;

  // get Cij

  NumberVector<NDIM, NumberVector<NDIM,ADScalar> > C = -rho*U.outerproduct(U);

  Scalar Cval = 0.;

  switch(i)
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
//template <typename Scalar>
//NumberVector<NDIM-1, NumberVector<NDIM-1,Scalar> > 
//  MASA::navierstokes_3d_variabledensity<Scalar>::Cij(Scalar x, 
//                                                     Scalar y, 
//                                                     Scalar z,
//                                                     Scalar t)
//{ 
//  NumberVector<NDIM-1, NumberVector<NDIM-1,Scalar> > Cij;
//
//  Scalar m1,m2,m3,rho;
//
//  rho = eval_exact_rho(x,y,z,t);
//  m1  = eval_exact_mD_1(x,y,z,t) + eval_exact_mC_1(x,y,z,t);
//  m2  = eval_exact_mD_2(x,y,z,t) + eval_exact_mC_2(x,y,z,t);
//  m3  = eval_exact_mD_3(x,y,z,t) + eval_exact_mC_3(x,y,z,t);
//
//  Cij[0][0] = -m1*m1/rho;
//  Cij[1][0] = -m2*m1/rho;
//  Cij[2][0] = -m3*m1/rho;
//  Cij[1][1] = -m2*m2/rho;
//  Cij[1][2] = -m2*m3/rho;
//  Cij[2][2] = -m3*m3/rho;
//
//  Cij[2][1] = Cij[1][2];
//  Cij[0][1] = Cij[1][0]
//  Cij[0][2] = Cij[2][0];
//
//  return Cij;
//
//}

//template <typename Scalar>
//NumberVector<NDIM-1, Scalar> 
//  MASA::navierstokes_3d_variabledensity<Scalar>::DivTauij(Scalar x1, 
//                                                          Scalar y1, 
//                                                          Scalar z1,
//                                                          Scalar t)
//{ 
//  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > FirstDerivType;
//  typedef DualNumber<FirstDerivType, NumberVector<NDIM, FirstDerivType> > SecondDerivType;
//  typedef SecondDerivType ADScalar;
//
//  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
//  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
//  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());
//
//  NumberVector<NDIM-1, NumberVector<NDIM-1,ADScalar> > Tij;
//
//  NumberVector<NDIM-1, NumberVector<NDIM-1,Scalar> > I = 
//    NumberVector<NDIM-1, Scalar>::identity;
//
//  NumberVector<NDIM-1, ADScalar> U;
//
//  ADScalar mu, rho;
//
//  mu    = eval_exact_mu(x,y,z,t);
//  rho   = eval_exact_rho(x,y,z,t);
//  U[0]  = (eval_exact_mD_1(x,y,z,t) + eval_exact_mC_1(x,y,z,t))/rho;
//  U[1]  = (eval_exact_mD_2(x,y,z,t) + eval_exact_mC_2(x,y,z,t))/rho;
//  U[2]  = (eval_exact_mD_3(x,y,z,t) + eval_exact_mC_3(x,y,z,t))/rho;
//
//  Tij = mu*(gradient(U) + transpose(gradient(U))) - 2./3.*mu*divergence(U)*I;
//  
//  return raw_value(divergence(Tij));
//
//}
// example of a public method called from eval_exact_t
//template <typename Scalar>
//Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_u(Scalar x, Scalar y1, Scalar z1)
//{
//  typedef DualNumber<Scalar, Scalar> OneDDerivType;
//  OneDDerivType y = OneDDerivType(y1,1);
//  OneDDerivType z = OneDDerivType(z1,1);
//
//  Scalar exact_u;
//  exact_u =   a *  helper_f(beta,kx,x) * helper_g(y).derivatives() *  helper_h(gamma,kz,z).derivatives();
//  return exact_u;
//}
//
//// public method
//template <typename Scalar>
//Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_v(Scalar x1, Scalar y, Scalar z1)
//{
//  typedef DualNumber<Scalar, Scalar> OneDDerivType;
//  OneDDerivType x = OneDDerivType(x1,1);
//  OneDDerivType z = OneDDerivType(z1,1);
//
//  Scalar exact_v;
//  exact_v = b * helper_f(beta,kx,x).derivatives() *  helper_g(y) * helper_h(gamma,kz,z).derivatives();
//  return exact_v;
//}
//
//// public method
//template <typename Scalar>
//Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_w(Scalar x1, Scalar y1, Scalar z)
//{
//  typedef DualNumber<Scalar, Scalar> OneDDerivType;
//  OneDDerivType x = OneDDerivType(x1,1);
//  OneDDerivType y = OneDDerivType(y1,1);
//
//  Scalar exact_w;
//  exact_w = c * helper_f(beta,kx,x).derivatives() * helper_g(y).derivatives() *  helper_h(gamma,kz,z);
//  return exact_w;
//}
//
//// public method
//template <typename Scalar>
//Scalar MASA::navierstokes_3d_variabledensity<Scalar>::eval_exact_p(Scalar x, Scalar y, Scalar z)
//{
//  Scalar P = d *  helper_f(beta,kx,x) * helper_gt(y) *  helper_h(gamma,kz,z);
//  return P;
//}


// ----------------------------------------
// Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::navierstokes_3d_variabledensity);


#endif // HAVE_METAPHYSICL
