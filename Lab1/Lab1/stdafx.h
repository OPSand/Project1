// stdafx.h : fichier Include pour les fichiers Include système standard,
// ou les fichiers Include spécifiques aux projets qui sont utilisés fréquemment,
// et sont rarement modifiés
//

#pragma once

#include "targetver.h"
#include <time.h>
#include <stdio.h>
#include <tchar.h>
#include "armadillo"
#include "lib.h"

// We only have to modify this line to change the type of all of our variables
#define TYPE double

// Matrix / vector types
typedef arma::Col<TYPE> arr; // it's good to be a pirate! (typed array/vector)
typedef arma::Mat<TYPE> matr; // typed matrix

// Program flow control variables - uncomment/comment to include/exclude a certain exercise
//#define EXA
#define EXB
//#define EXC
//#define EXD
//#define EXE

// Preprocessing instruction ...
/*template<typename T1>
extern "C" bool lu
  (
         Mat<typename T1::elem_type>&    L,
         Mat<typename T1::elem_type>&    U,
  const Base<typename T1::elem_type,T1>& X,
  const typename arma_blas_type_only<typename T1::elem_type>::result* junk = 0
  )
  {
  arma_extra_debug_sigprint();
  arma_ignore(junk);
  
  arma_debug_check( (&L == &U), "lu(): L and U are the same object");
  
  const bool status = auxlib::lu(L, U, X);
  
  if(status == false)
    {
    L.reset();
    U.reset();
    arma_bad("lu(): failed to converge", false);
    }
  
  return status;
  };*/