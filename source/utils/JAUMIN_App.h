#ifndef included_JAUMIN_App
#define included_JAUMIN_App

// NOTE:
//      The following headers cover most (not every) headers that
//      a JAUMIN application can use. They are included here just
//      for convenience. The downside is that it may include some
//      headers that you do not need at all, so the compile time
//      is longer. It's up to you.

#include "Pointer.h"
#include "Array.h"
#include "IEEE.h"
#include "MPI.h"
#include "PIO.h"
#include "Utilities.h"

#include "GridGeometry.h"
#include "GridTopology.h"
#include "PatchHierarchy.h"
#include "HierarchyTimeIntegrator.h"
#include "JaVisDataWriter.h"
#include "JAUMINManager.h"
#include "InputManager.h"
#include "RestartManager.h"
#include "TimerManager.h"
#include "VariableDatabase.h"
#include "PatchHierarchy.h"
#include "PatchLevel.h"
#include "PatchLevelLocalIterator.h"

#include "StandardComponentPatchStrategy.h"
#include "TimeIntegratorLevelStrategy.h"
#include "NumericalIntegratorComponent.h"
#include "DtIntegratorComponent.h"
#include "InitializeIntegratorComponent.h"
#include "CopyIntegratorComponent.h"
#include "MemoryIntegratorComponent.h"
#include "ReductionIntegratorComponent.h"

#include "Patch.h"
#include "PatchTopology.h"
#include "PatchGeometry.h"
#include "EntityUtilities.h"
#include "EntitySet.h"
#include "PatchBoundary.h"

#include "CellData.h"
#include "NodeData.h"
#include "EdgeData.h"
#if NDIM > 2
#include "FaceData.h"
#endif
#include "NodeVariable.h"
#include "CellVariable.h"
#include "EdgeVariable.h"
#if NDIM > 2
#include "FaceVariable.h"
#endif

#include "OuterNodeData.h"
#include "OuterEdgeData.h"
#if NDIM > 2
#include "OuterFaceData.h"
#endif

#include "VectorVariable.h"
#include "CSRMatrixVariable.h"
#include "VectorData.h"
#include "CSRMatrixData.h"
#include "LinearSolverManager.h"
#include "BaseLinearSolver.h"

#endif // included_JAUMIN_App
