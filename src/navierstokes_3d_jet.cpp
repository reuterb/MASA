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
template <typename Scalar>
Scalar helper_U(Scalar y);

template <typename Scalar>
Scalar helper_Z(Scalar y);

template <typename Scalar>
Scalar helper_T(Scalar t);

template <typename Scalar, typename Scalar2>
Scalar rhoMap(Scalar Z, Scalar2 at);

//template <typename Scalar>
//NumberVector DivTauij(Scalar x1, Scalar y1, Scalar z1, Scalar t)

typedef ShadowNumber<double, long double> RawScalar;
const unsigned int NDIM = 3;  // Need time dependence!! But put in manually!
typedef DualNumber<RawScalar, NumberVector<NDIM, RawScalar> > D1Type;
typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
typedef D2Type ADType;

using namespace MASA;

double off = .25;
double delta = .25;

template <typename Scalar>
MASA::navierstokes_3d_jet<Scalar>::navierstokes_3d_jet()
{
  this->mmsname = "navierstokes_3d_jet";
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
int MASA::navierstokes_3d_jet<Scalar>::init_var()
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
  err += this->set_var("one",1);

  return err;

} // done with init_var

// ----------------------------------------
// Source Terms
// ----------------------------------------

// public static method, that can be called from eval_q_t
template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_q_omega(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  return 0.;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_q_phi(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  return 0.;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_q_rho(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef D2Type ADScalar;

  typedef DualNumber<Scalar, NumberVector<NDIM+1, Scalar> > D1TimeType;
  typedef D1TimeType ADTimeScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar U = helper_U(y) * helper_T(t1);

  ADScalar zVar = helper_U(y) * helper_T(t1);
  ADScalar rho  = rhoMap(zVar,at);

  ADScalar V = -rho.derivatives()[1];

  // Time part

  ADTimeScalar xt = 
    ADTimeScalar(x1,NumberVectorUnitVector<NDIM+1, 0, Scalar>::value());
  ADTimeScalar yt = 
    ADTimeScalar(y1,NumberVectorUnitVector<NDIM+1, 1, Scalar>::value());
  ADTimeScalar zt = 
    ADTimeScalar(z1,NumberVectorUnitVector<NDIM+1, 2, Scalar>::value());
  ADTimeScalar tt = 
    ADTimeScalar(t1,NumberVectorUnitVector<NDIM+1, 3, Scalar>::value());

  ADTimeScalar zVarTime = helper_U(y1) * helper_T(tt);

  ADTimeScalar rhoTime = rhoMap(zVarTime,at);

  ADScalar mC = rho*V;

  // Rho equation residuals
  Scalar Q_rho = 
                // Temporal part
                rhoTime.derivatives()[3]

        +
                // Div mC
                raw_value(mC.derivatives()[1]);

  return Q_rho;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_q_z(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef D2Type ADScalar;

  typedef DualNumber<Scalar, NumberVector<NDIM+1, Scalar> > D1TimeType;
  typedef D1TimeType ADTimeScalar;

  // Treat velocity as a vector
  NumberVector<NDIM, ADScalar> U;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zVar = helper_U(y) * helper_T(t1);

  ADScalar rho = rhoMap(zVar,at);
  ADScalar mu  = rho; 

  U[0] =  helper_U(y) * helper_T(t1);
  U[1] = -rho.derivatives()[1];
  U[2] = 0.;

  // Time part

  ADTimeScalar xt = 
    ADTimeScalar(x1,NumberVectorUnitVector<NDIM+1, 0, Scalar>::value());
  ADTimeScalar yt = 
    ADTimeScalar(y1,NumberVectorUnitVector<NDIM+1, 1, Scalar>::value());
  ADTimeScalar zt = 
    ADTimeScalar(z1,NumberVectorUnitVector<NDIM+1, 2, Scalar>::value());
  ADTimeScalar tt = 
    ADTimeScalar(t1,NumberVectorUnitVector<NDIM+1, 3, Scalar>::value());

  ADTimeScalar zTime = helper_U(yt) * helper_T(tt);

  // Z equation residuals
  Scalar Q_z = 
                // Temporal part
                zTime.derivatives()[3]

        +
                // U dot grad(Z)
                raw_value(U.dot(gradient(zVar)))
        
        -
                // viscous
       raw_value(   1.0/re*1.0/sc*(
                    1.0/rho*gradient(mu).dot(gradient(zVar))
                  +  mu/rho*divergence(gradient(zVar))));

  return Q_z;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_q_m1(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef DualNumber<D2Type, NumberVector<NDIM, D2Type> > D3Type;
  typedef D3Type ADScalar;

  typedef DualNumber<Scalar, NumberVector<NDIM+1, Scalar> > D1TimeType;
  typedef DualNumber<D1TimeType, NumberVector<NDIM+1, D1TimeType> > D2TimeType;
  typedef D2TimeType ADTimeScalar;

  // Treat velocity,momentum as a vector
  NumberVector<NDIM, D2Type> U;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  D2Type x2 = D2Type(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  D2Type y2 = D2Type(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  D2Type z2 = D2Type(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  D2Type rho = rhoMap(helper_U(y2) * helper_T(t1), at);

  D2Type mu  = rho;

  U[0] = helper_U(y2) * helper_T(t1);
  U[1] = -rhoMap(helper_U(y) * helper_T(t1), at).derivatives()[1];

  // Time part

  ADTimeScalar xt = 
    ADTimeScalar(x1,NumberVectorUnitVector<NDIM+1, 0, Scalar>::value());
  ADTimeScalar yt = 
    ADTimeScalar(y1,NumberVectorUnitVector<NDIM+1, 1, Scalar>::value());
  ADTimeScalar zt = 
    ADTimeScalar(z1,NumberVectorUnitVector<NDIM+1, 2, Scalar>::value());
  ADTimeScalar tt = 
    ADTimeScalar(t1,NumberVectorUnitVector<NDIM+1, 3, Scalar>::value());

  ADTimeScalar m1Time = raw_value(rho) * helper_U(y1) * helper_T(tt);

  NumberVector<NDIM, NumberVector<NDIM, Scalar> > I = 
    NumberVector<NDIM, Scalar>::identity();
 
  NumberVector<NDIM, NumberVector<NDIM,D2Type> > Tau = 
   mu * (gradient(U) + transpose(gradient(U)) - 2./3.*divergence(U)*I);

  D2Type C12 = -rho*U[0]*U[1];

  // m1 equation residuals
  Scalar Q_m1 = 
                // Temporal part
                raw_value(m1Time.derivatives()[3])

        -
                // RHS
                raw_value(C12.derivatives()[1]
                + 1./re*Tau[0][1].derivatives()[1]);

  return Q_m1;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_q_m3(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{

  return 0.;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_q_top(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{

  return 0.;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_q_bottom(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{

  return 0.;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_q_phiTop(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{

  return 0.;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_q_phiBottom(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{

  return 0.;

}

// ----------------------------------------
// Analytical Terms
// ----------------------------------------
//
template <typename Scalar>
Scalar helper_U(Scalar y)
{
  Scalar func;
  func = .5 * (std::tanh( (y + off)/delta ) - std::tanh( (y - off)/delta ));

  return func;
}

template <typename Scalar>
Scalar helper_T(Scalar t)
{
  Scalar func;

  //func = 1.;
  //func = (1. + 10.*t + 20.*t*t + 30.*t*t*t); // + 1.e13*t*t + 0.*30.*t*t*t);//std::exp(t);
  //func = std::exp(t/1e-3);
  func = std::exp(-t);

  return func;
}

template <typename Scalar, typename Scalar2>
Scalar rhoMap(Scalar Z, Scalar2 at)
{
  Scalar rho;

  rho = ( 1.0 + at ) / ( 1.0 + at - 2.0*at*Z );

  return rho;

}

//
// main functions
//
// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_exact_mD_1(Scalar x1, 
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

  Scalar zVal, rho;

  zVal = raw_value(helper_U(y) * helper_T(t1));

  rho = rhoMap(zVal,at);

  exact_1 = raw_value(rho * helper_U(y) * helper_T(t1));

  return exact_1;
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_exact_mD_2(Scalar x1, 
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  return 0.;
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_exact_mD_3(Scalar x1, 
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  return 0.;
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_exact_mC_1(Scalar x1, 
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  return 0.;
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_exact_mC_2(Scalar x1, 
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef D1Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar rho,zVal;

  Scalar exact_2;

  zVal = helper_U(y)*helper_T(t1);
  rho  = rhoMap(zVal,at);

  exact_2 = -rho.derivatives()[1];

  return exact_2;
}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_exact_mean_mC_2(Scalar x1, 
                                                                           Scalar y1, 
                                                                           Scalar z1,
                                                                           Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef D1Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar rho,zVal;

  Scalar exact_2;

  zVal = helper_U(y)*helper_T(t1);
  rho  = rhoMap(zVal,at);

  exact_2 = -rho.derivatives()[1];

  return exact_2;

}

// public method
template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_exact_mC_3(Scalar x1, 
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  return 0.;
}

// public method

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_exact_rho(Scalar x1,
                                                                     Scalar y1, 
                                                                     Scalar z1,
                                                                     Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef D1Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  Scalar zVar = raw_value(helper_U(y) * helper_T(t1));

  Scalar exact_rho = rhoMap(zVar,at);
  
  return exact_rho;
}

// public method

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_exact_drho(Scalar x1,
                                                                      Scalar y1, 
                                                                      Scalar z1,
                                                                      Scalar t1)
{
  return 0.;  // Nothing!
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_exact_z(Scalar x1,
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

  exact_z = raw_value(helper_U(y) * helper_T(t1));

  return exact_z;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_exact_mu(Scalar x1,
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

  Scalar zVar = raw_value(helper_U(y) * helper_T(t1));

  Scalar exact_rho = rhoMap(zVar,at);

  exact_mu = exact_rho;

  return exact_mu;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_exact_omega(Scalar x1, 
                                                                       Scalar y1, 
                                                                       Scalar z1, 
                                                                       Scalar t1)
{
  return 0.;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_exact_phi(Scalar x1,
                                                                     Scalar y1, 
                                                                     Scalar z1, 
                                                                     Scalar t1)
{
  return 0.;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_exact_div_mC(Scalar x1,
                                                                        Scalar y1, 
                                                                        Scalar z1, 
                                                                        Scalar t1)
{

  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef D2Type ADScalar;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar rho,zVal,mC;

  zVal = helper_U(y)*helper_T(t1);
  rho  = rhoMap(zVal,at);

  mC = -rho.derivatives()[1];

  return raw_value(mC.derivatives()[1]);

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_exact_RHSomega(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  return 0.;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_exact_RHSphi(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  return 0.;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_g_DivTau(Scalar x1, Scalar y1, Scalar z1, Scalar t1, int i)
{
  return 0.;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_g_C(Scalar x1, Scalar y1, Scalar z1, Scalar t1, int i)
{
  return 0.;
}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_exact_RHSz(Scalar x1, Scalar y1, Scalar z1, Scalar t1)
{
  typedef DualNumber<Scalar, NumberVector<NDIM, Scalar> > D1Type;
  typedef DualNumber<D1Type, NumberVector<NDIM, D1Type> > D2Type;
  typedef D2Type ADScalar;

  // Treat velocity, momentum as a vector
  NumberVector<NDIM, ADScalar> U;

  ADScalar x = ADScalar(x1,NumberVectorUnitVector<NDIM, 0, Scalar>::value());
  ADScalar y = ADScalar(y1,NumberVectorUnitVector<NDIM, 1, Scalar>::value());
  ADScalar z = ADScalar(z1,NumberVectorUnitVector<NDIM, 2, Scalar>::value());

  ADScalar zVar = helper_U(y) * helper_T(t1);

  ADScalar rho = rhoMap(zVar,at);

  ADScalar mu  = rho; 

  U[0] =  helper_U(y) * helper_T(t1);
  U[1] = -rho.derivatives()[1];
  U[2] = 0.;

  Scalar RHS_z = raw_value(-U.dot(gradient(zVar)))
               + raw_value(1.0/(re*sc)*(1.0/rho*gradient(mu).dot(gradient(zVar))
                   + mu/rho*divergence(gradient(zVar))));

  return RHS_z;

}

template <typename Scalar>
Scalar MASA::navierstokes_3d_jet<Scalar>::eval_g_jacF(Scalar x1, Scalar y1, Scalar z1, Scalar t1, int i)
{
  return 0.;
}

// ----------------------------------------
// Template Instantiation(s)
// ----------------------------------------

MASA_INSTANTIATE_ALL(MASA::navierstokes_3d_jet);


#endif // HAVE_METAPHYSICL
