// stdafx.h : fichier Include pour les fichiers Include système standard,
// ou les fichiers Include spécifiques aux projets qui sont utilisés fréquemment,
// et sont rarement modifiés
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

// TODO: faites référence ici aux en-têtes supplémentaires nécessaires au programme