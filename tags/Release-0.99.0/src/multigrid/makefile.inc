LDIR = src/multigrid

# Sources and objects
SOURCES += $(LDIR)/MGLevels.f90
SOURCES += $(LDIR)/GaussSeidel.f90
SOURCES += $(LDIR)/GaussSeidelARRAY.f90
SOURCES += $(LDIR)/GaussSeidelARRAYback.f90
SOURCES += $(LDIR)/GaussSeidelReverse.f90
SOURCES += $(LDIR)/GaussSeidelSpecialPoints.f90
SOURCES += $(LDIR)/GetSparseLinearLaplaceMatrix.f90
SOURCES += $(LDIR)/MGPreProcess.f90
SOURCES += $(LDIR)/MGProlongation3DExplicit.f90
SOURCES += $(LDIR)/MGRestriction3DExplicit.f90
SOURCES += $(LDIR)/MultiGridSolveLinearSystem.f90
SOURCES += $(LDIR)/MultiGridSolver.f90
