LDIR = src/IO

# Sources and objects
SOURCES += $(LDIR)/CheckError.f90			
SOURCES += $(LDIR)/StoreData.f90
# SOURCES += $(LDIR)/StoreDataVTK.f90
SOURCES += $(LDIR)/StoreRealArray3D.f90
SOURCES += $(LDIR)/CloseIOFiles.f90
SOURCES += $(LDIR)/StoreDataAscii.f90
SOURCES += $(LDIR)/StoreSparseMatrix.f90
SOURCES += $(LDIR)/GatherData.f90
SOURCES += $(LDIR)/StoreKinematicData.f90		
SOURCES += $(LDIR)/ReadInputFileParameters.f90	
SOURCES += $(LDIR)/StoreRealArray.f90
SOURCES += $(LDIR)/StoreWaveGauges.f90
