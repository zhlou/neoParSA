/*
 *       Filename:  staticLam.h
 *        Created:  09/24/2015 03:47:26 PM
 *         Author:  Zhihao Lou
 *
 *           Note:  This cooling schedule relies heavily on the properties
 *                  of Rastrigin function so it has to be here.
 */

#ifndef STATICLAM_H
#define STATICLAM_H

#include <vector>
#include <libxml/tree.h>

#include "dynDebug.h"
#include "aState.h"

class staticLam {
private:
    std::vector<double> energy;
    std::vector<double> variance;
}



#endif
