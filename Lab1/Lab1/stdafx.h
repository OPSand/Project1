// stdafx.h�: fichier Include pour les fichiers Include syst�me standard,
// ou les fichiers Include sp�cifiques aux projets qui sont utilis�s fr�quemment,
// et sont rarement modifi�s
//

#pragma once

#include "targetver.h"

#include <stdio.h>
#include <tchar.h>
#include "armadillo"

// We only have to modify this line to change the type of all of our variables
#define TYPE double
//#define EXA
typedef arma::vec arr; // it's good to be a pirate! (this is our untyped array/vector type)
typedef arma::mat matr; // matrix class (now untyped)

// TODO: faites r�f�rence ici aux en-t�tes suppl�mentaires n�cessaires au programme