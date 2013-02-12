LDIR  = src/OpenFoam
IODIR = src/OpenFoam/IO
OFDIR = src/OpenFoam/OpenFoam

# Sources and objects
SOURCES += $(OFDIR)/calculateKinematics.f90
SOURCES += $(OFDIR)/interpolationInitialize.f90
SOURCES += $(OFDIR)/makeStencils.f90
SOURCES += $(OFDIR)/openFoamNearestNeighbour.f90
SOURCES += $(OFDIR)/setDeltaT.f90
#
SOURCES += $(LDIR)/interpolation_functions.f90
SOURCES += $(LDIR)/interpolation.f90			
#
SOURCES += $(IODIR)/openFoamInterfaceEta.f90
SOURCES += $(IODIR)/openFoamInterfaceU.f90
SOURCES += $(IODIR)/writeOceanwave3D.f90



