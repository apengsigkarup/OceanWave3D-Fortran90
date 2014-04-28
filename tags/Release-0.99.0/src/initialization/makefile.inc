LDIR = src/initialization

# Sources and objects
SOURCES += $(LDIR)/BuildStencilsGridX.f90
SOURCES += $(LDIR)/BuildStencilsGridY.f90
SOURCES += $(LDIR)/CloseVariables.f90
SOURCES += $(LDIR)/GammaFunctions.f90
SOURCES += $(LDIR)/Initialize.f90
SOURCES += $(LDIR)/InitializeVariables.f90
SOURCES += $(LDIR)/PreProcessDiffStencils.f90
SOURCES += $(LDIR)/PreprocessRelaxationZones.f90
SOURCES += $(LDIR)/PreprocessPDampingZones.f90
SOURCES += $(LDIR)/SetupCompDomain.f90
SOURCES += $(LDIR)/SetupInitialConditions.f90
SOURCES += $(LDIR)/setupWavePaddle.f90
