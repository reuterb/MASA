// -*-c++-*-
//
//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MASA - Manufactured Analytical Solutions Abstraction Library
//
// Copyright (C) 2010 The PECOS Development Team
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
// $Author$
// $Id$
//
// axi_euler.cpp: These are the MASA class member functions and constructors
// for the Axisymmetric Euler Equations
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//

#include <masa_internal.h> 

using namespace MASA;

/* ------------------------------------------------
 *
 *         Axisymmetric EULER EQUATION 1D
 *
 *
 *
 * -----------------------------------------------
 */ 

MASA::axi_euler::axi_euler()
{
  mmsname = "axisymmetric_euler";
  dimension=2;

  register_var("R",&R);
  //register_var("k",&k);

  register_var("p_0",&p_0);
  register_var("p_1",&p_1);
  register_var("rho_0",&rho_0);
  register_var("rho_1",&rho_1);
  register_var("u_1",&u_1);
  register_var("w_0",&w_0);
  register_var("w_1",&w_1);
  register_var("a_pr",&a_pr);
  register_var("a_pz",&a_pz);
  register_var("a_rhor",&a_rhor);
  register_var("a_rhoz",&a_rhoz);
  register_var("a_ur",&a_ur);
  register_var("a_uz",&a_uz);
  register_var("a_wr",&a_wr);
  register_var("a_wz",&a_wz);
  register_var("L",&L);

  register_var("Gamma",&Gamma);
  register_var("mu",&mu);

}//done with constructor

int MASA::axi_euler::init_var()
{

  int err = 0;

  // currently randomly generated
  err += set_var("R",1.01);

  //err += set_var("k",.1791);
  err += set_var("p_0",.381);
  err += set_var("p_1",.1791);
  err += set_var("rho_0",.1791);
  err += set_var("rho_1",.1791);
  err += set_var("u_1",.1791);
  err += set_var("w_0",.1791);
  err += set_var("w_1",.1791);
  err += set_var("a_pr",.1791);
  err += set_var("a_pz",.1791);
  err += set_var("a_rhor",.1791);
  err += set_var("a_rhoz",.1791);
  err += set_var("a_ur",.1791);
  err += set_var("a_uz",.1791);
  err += set_var("a_wr",.1791);
  err += set_var("a_wz",.1791);
  err += set_var("L",1);

  err += set_var("Gamma",.1791);
  err += set_var("mu",.1791);

  return err;

} // done with variable initializer


// ----------------------------------------
//   Source Terms
// ----------------------------------------

double MASA::axi_euler::eval_q_u(double r,double z)
{
  double Q_u;
  Q_u = p_1 * cos(a_pr * PI * r / L) * cos(a_pz * PI * z / L) * a_pr * PI / L - u_1 * u_1 * pow(sin(a_uz * PI * z / L), 0.2e1) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * rho_1 * sin(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L) * a_rhor * PI / L + u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * rho_1 * cos(a_rhor * PI * r / L) * cos(a_rhoz * PI * z / L) * sin(a_uz * PI * z / L) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * a_rhoz * PI / L - 0.2e1 * u_1 * u_1 * pow(sin(a_uz * PI * z / L), 0.2e1) * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_ur * PI * r / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * a_ur * PI / L + u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * cos(a_uz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * a_uz * PI / L + u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * w_1 * cos(a_wr * PI * r / L) * cos(a_wz * PI * z / L) * sin(a_uz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * a_wz * PI / L + u_1 * u_1 * pow(sin(a_uz * PI * z / L), 0.2e1) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) / r;
  return(Q_u);
}

double MASA::axi_euler::eval_q_w(double r,double z)
{
  double Q_w;
  Q_w = -p_1 * sin(a_pr * PI * r / L) * sin(a_pz * PI * z / L) * a_pz * PI / L - u_1 * sin(a_uz * PI * z / L) * rho_1 * sin(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * (cos(a_ur * PI * r / L) - 0.1e1) * a_rhor * PI / L + pow(w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L), 0.2e1) * rho_1 * cos(a_rhor * PI * r / L) * cos(a_rhoz * PI * z / L) * a_rhoz * PI / L - u_1 * sin(a_uz * PI * z / L) * sin(a_ur * PI * r / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * a_ur * PI / L - u_1 * sin(a_uz * PI * z / L) * w_1 * sin(a_wr * PI * r / L) * sin(a_wz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (cos(a_ur * PI * r / L) - 0.1e1) * a_wr * PI / L + (0.2e1 * w_0 + 0.2e1 * w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * w_1 * cos(a_wr * PI * r / L) * cos(a_wz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * a_wz * PI / L + u_1 * sin(a_uz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * (cos(a_ur * PI * r / L) - 0.1e1) / r;
  return(Q_w);
}

double MASA::axi_euler::eval_q_rho(double r,double z)
{
  double Q_rho;
  Q_rho = -(cos(a_ur * PI * r / L) - 0.1e1) * a_rhor * PI * rho_1 * u_1 * sin(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L) * sin(a_uz * PI * z / L) / L + (w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L) + w_0) * a_rhoz * PI * rho_1 * cos(a_rhor * PI * r / L) * cos(a_rhoz * PI * z / L) / L - (rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L) + rho_0) * a_ur * PI * u_1 * sin(a_ur * PI * r / L) * sin(a_uz * PI * z / L) / L + (rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L) + rho_0) * a_wz * PI * w_1 * cos(a_wr * PI * r / L) * cos(a_wz * PI * z / L) / L + (rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L) + rho_0) * (cos(a_ur * PI * r / L) - 0.1e1) * u_1 * sin(a_uz * PI * z / L) / r;
  return(Q_rho);
}

double MASA::axi_euler::eval_q_e(double r,double z)
{
  double Q_e;
  Q_e = Gamma * cos(a_pz * PI * z / L) * cos(a_pr * PI * r / L) * p_1 * u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_uz * PI * z / L) * a_pr * PI / L / (Gamma - 0.1e1) - Gamma * sin(a_pz * PI * z / L) * sin(a_pr * PI * r / L) * p_1 * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * a_pz * PI / L / (Gamma - 0.1e1) - sin(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L) * sin(a_uz * PI * z / L) * u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * (pow(w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L), 0.2e1) + pow(sin(a_uz * PI * z / L), 0.2e1) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * u_1 * u_1) * a_rhor * PI * rho_1 / L / 0.2e1 + cos(a_rhor * PI * r / L) * cos(a_rhoz * PI * z / L) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * (pow(w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L), 0.2e1) + pow(sin(a_uz * PI * z / L), 0.2e1) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * u_1 * u_1) * a_rhoz * PI * rho_1 / L / 0.2e1 - sin(a_uz * PI * z / L) * sin(a_ur * PI * r / L) * (p_0 + p_1 * sin(a_pr * PI * r / L) * cos(a_pz * PI * z / L)) * a_ur * PI * u_1 * Gamma / L / (Gamma - 0.1e1) - sin(a_uz * PI * z / L) * sin(a_ur * PI * r / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (pow(w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L), 0.2e1) + 0.3e1 * pow(sin(a_uz * PI * z / L), 0.2e1) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * u_1 * u_1) * a_ur * PI * u_1 / L / 0.2e1 + (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * u_1 * u_1 * cos(a_uz * PI * z / L) * sin(a_uz * PI * z / L) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * a_uz * PI / L - sin(a_uz * PI * z / L) * (cos(a_ur * PI * r / L) - 0.1e1) * u_1 * w_1 * sin(a_wr * PI * r / L) * sin(a_wz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L)) * a_wr * PI / L + cos(a_wr * PI * r / L) * cos(a_wz * PI * z / L) * (p_0 + p_1 * sin(a_pr * PI * r / L) * cos(a_pz * PI * z / L)) * a_wz * PI * w_1 * Gamma / L / (Gamma - 0.1e1) + cos(a_wr * PI * r / L) * cos(a_wz * PI * z / L) * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (0.3e1 * pow(w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L), 0.2e1) + pow(sin(a_uz * PI * z / L), 0.2e1) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * u_1 * u_1) * a_wz * PI * w_1 / L / 0.2e1 + sin(a_uz * PI * z / L) * (cos(a_ur * PI * r / L) - 0.1e1) * u_1 * (p_0 + p_1 * sin(a_pr * PI * r / L) * cos(a_pz * PI * z / L)) * Gamma / (Gamma - 0.1e1) / r + sin(a_uz * PI * z / L) * (cos(a_ur * PI * r / L) - 0.1e1) * u_1 * (rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L)) * (pow(w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L), 0.2e1) + pow(sin(a_uz * PI * z / L), 0.2e1) * pow(cos(a_ur * PI * r / L) - 0.1e1, 0.2e1) * u_1 * u_1) / r / 0.2e1;
  return(Q_e);
}

// ----------------------------------------
//   Analytical Solutions
// ----------------------------------------

double MASA::axi_euler::eval_an_u(double r,double z)
{
  double u_an;
  u_an = u_1 * (cos(a_ur * PI * r / L) - 0.1e1) * sin(a_uz * PI * z / L);
  return u_an;
}

double MASA::axi_euler::eval_an_w(double r,double z)
{
  double w_an;
  w_an = w_0 + w_1 * cos(a_wr * PI * r / L) * sin(a_wz * PI * z / L);
  return w_an;
}

double MASA::axi_euler::eval_an_p(double r,double z)
{
  double p_an;
  p_an = p_0 + p_1 * sin(a_pr * PI * r / L) * cos(a_pz * PI * z / L);
  return p_an;
}

double MASA::axi_euler::eval_an_rho(double r,double z)
{
  double rho_an;
  rho_an = rho_0 + rho_1 * cos(a_rhor * PI * r / L) * sin(a_rhoz * PI * z / L);
  return rho_an;
}
