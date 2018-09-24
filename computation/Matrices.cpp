#include "computations/EigenDenseBackend.h"
#include "computations/EigenSparseBackend.h"
#include "Matrices.h"
#include "Matrices-inline.h"

const std::vector<PeriodicBoxDof> DefaultBoxDofs = {};
const std::vector<SpringType> DefaultSprings = {SpringType::Surface, SpringType::Volume};

template class Matrices<EigenDenseBackend>;
//template class Matrices<EigenSparseBackend>;
