// stdafx.h�: fichier Include pour les fichiers Include syst�me standard,
// ou les fichiers Include sp�cifiques aux projets qui sont utilis�s fr�quemment,
// et sont rarement modifi�s
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
//#define EXA
//#define EXB
//#define EXC
#define EXD
typedef arma::vec arr; // it's good to be a pirate! (this is our untyped array/vector type)
typedef arma::mat matr; // matrix class (now untyped)

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

// TODO: faites r�f�rence ici aux en-t�tes suppl�mentaires n�cessaires au programme