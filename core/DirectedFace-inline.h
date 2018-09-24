#ifndef DIRECTEDFACE_INLINE_H
#define DIRECTEDFACE_INLINE_H

#include "DirectedFace.h"

void DirectedFace::updateAdditionalInterfacialTension() { 
  _additionalInterfacialTension = cell()->type()->additionalInterfacialTensionWith(otherCell()->type()); 
}

#endif /* DIRECTEDFACE_INLINE_H */

