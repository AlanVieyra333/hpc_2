// -*- C++ -*-
// $Id$

/**
 * Code generated by the The ACE ORB (TAO) IDL Compiler v2.2.0
 * TAO and the TAO IDL Compiler have been developed by:
 *       Center for Distributed Object Computing
 *       Washington University
 *       St. Louis, MO
 *       USA
 *       http://www.cs.wustl.edu/~schmidt/doc-center.html
 * and
 *       Distributed Object Computing Laboratory
 *       University of California at Irvine
 *       Irvine, CA
 *       USA
 * and
 *       Institute for Software Integrated Systems
 *       Vanderbilt University
 *       Nashville, TN
 *       USA
 *       http://www.isis.vanderbilt.edu/
 *
 * Information about TAO is available at:
 *     http://www.cs.wustl.edu/~schmidt/TAO.html
 **/

// TAO_IDL - Generated from
// be/be_codegen.cpp:1616

#ifndef MULTMATRIZI_TQU0FF_H_
#define MULTMATRIZI_TQU0FF_H_

#include "MultMatrizS.h"

#if !defined (ACE_LACKS_PRAGMA_ONCE)
#pragma once
#endif /* ACE_LACKS_PRAGMA_ONCE */

class  MultMatriz_i
  : public virtual POA_MultMatriz
{
public:
  // Constructor
  MultMatriz_i (void);

  // Destructor
  virtual ~MultMatriz_i (void);

  virtual
  ::CORBA::Float mult (
    ::CORBA::Float a,
    ::CORBA::Float b);
};


#endif /* MULTMATRIZI_H_  */
