LDIR = src/curvilinear

# Sources and objects
SOURCES += $(LDIR)/BuildLinearSystemTransformedCurvilinear.f90
SOURCES += $(LDIR)/BuildLinearSystemTransformedMatrixCurvilinear.f90
SOURCES += $(LDIR)/ComputeNormalVectors.f90
SOURCES += $(LDIR)/DetermineCurvilinearTransform2D.f90
SOURCES += $(LDIR)/DetermineGenericStencils.f90
SOURCES += $(LDIR)/DetermineGenericStencilsUniform.f90
SOURCES += $(LDIR)/DiffXuniform2D.f90
SOURCES += $(LDIR)/DiffXuniform3D.f90
SOURCES += $(LDIR)/DiffYuniform2D.f90
SOURCES += $(LDIR)/DiffYuniform3D.f90
SOURCES += $(LDIR)/DiffZuniform3D.f90
