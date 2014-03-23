LDIR = src/functions

# Sources and objects
SOURCES += $(LDIR)/ALLOCATE_Wavefield_Type.f90
SOURCES += $(LDIR)/ASSIGN_IncidentWavefield.f90
SOURCES += $(LDIR)/AnalyticWaveMaker2D.f90
SOURCES += $(LDIR)/BESSH.f90
SOURCES += $(LDIR)/BESSHP.f90
SOURCES += $(LDIR)/BESSJ0.f90
SOURCES += $(LDIR)/BESSJ1.f90
SOURCES += $(LDIR)/BESSJP.f90
SOURCES += $(LDIR)/BESSY.f90
SOURCES += $(LDIR)/BESSY0.f90
SOURCES += $(LDIR)/BESSY1.f90
SOURCES += $(LDIR)/BESSYP.f90
SOURCES += $(LDIR)/BuildStencil_1D_Uneven.f90
SOURCES += $(LDIR)/BuildLinearSystem.f90
SOURCES += $(LDIR)/BuildPDampMatrices.f90
SOURCES += $(LDIR)/BuildLinearSystemMatrix.f90
SOURCES += $(LDIR)/CSRdiaREORDER.f90
SOURCES += $(LDIR)/CleanSparseMatrixCOO.f90
SOURCES += $(LDIR)/CleanSparseMatrixCOOold.f90
SOURCES += $(LDIR)/ConstructTableCrossDerivatives.f90
SOURCES += $(LDIR)/ConstructTableCrossDerivatives_Curvilinear.f90
SOURCES += $(LDIR)/ConvertCOOtoCSR.f90
SOURCES += $(LDIR)/ConvertCSRtoCOO.f90
SOURCES += $(LDIR)/DEALLOCATE_Wavefield_Type.f90
SOURCES += $(LDIR)/DetermineTransformationConstants.f90
SOURCES += $(LDIR)/DetermineTransformationConstantsArray.f90
SOURCES += $(LDIR)/DfDx_1D_Uneven.f90
SOURCES += $(LDIR)/DiffStencilsZ_modif.f90
SOURCES += $(LDIR)/DiffXEven.f90
SOURCES += $(LDIR)/DiffXEvenCurvilinear.f90
SOURCES += $(LDIR)/DiffXEven_CD.f90
SOURCES += $(LDIR)/DiffXXEvenCurvilinear.f90
SOURCES += $(LDIR)/DiffXXEvenCurvilinear_CD.f90
SOURCES += $(LDIR)/DiffXuniform3D_CD.f90
SOURCES += $(LDIR)/DiffYEven.f90
SOURCES += $(LDIR)/DiffYEvenCurvilinear.f90
SOURCES += $(LDIR)/DiffYEven_CD.f90
SOURCES += $(LDIR)/DiffYYEvenCurvilinear.f90
SOURCES += $(LDIR)/DiffYYEvenCurvilinear_CD.f90
SOURCES += $(LDIR)/DiffYuniform3D_CD.f90
SOURCES += $(LDIR)/DiffZArbitrary.f90
SOURCES += $(LDIR)/DiffZArbitrary_CD.f90
SOURCES += $(LDIR)/DiffZuniform3D_CD.f90
SOURCES += $(LDIR)/DiscardGhostLayer.f90
SOURCES += $(LDIR)/DiscardSpaces.f90
SOURCES += $(LDIR)/FILTERING.f90
SOURCES += $(LDIR)/FILTERING_NEW.f90
SOURCES += $(LDIR)/FILTERING_SWENSE.f90
SOURCES += $(LDIR)/FactorPreconditioner.f90
SOURCES += $(LDIR)/FilterInit.f90
SOURCES += $(LDIR)/McCamyFuchs54LocalXY.f90
SOURCES += $(LDIR)/PreProcessDiffStencilsZ.f90
SOURCES += $(LDIR)/PrepareFullOperatorStencils.f90
SOURCES += $(LDIR)/PrepareFullOperatorStencilsForGaussSeidel.f90
SOURCES += $(LDIR)/PrepareFullOperatorStencilsForGaussSeidelARRAY.f90
SOURCES += $(LDIR)/PreparePreconditioner.f90
SOURCES += $(LDIR)/PreparePreconditionerOLD.f90
SOURCES += $(LDIR)/RelaxationModule.f90
SOURCES += $(LDIR)/RelaxationModule_new.f90
SOURCES += $(LDIR)/SetupCompDomain_test.f90
SOURCES += $(LDIR)/TaylorFDStencils1DArbitrary.f90
SOURCES += $(LDIR)/UpdateGhostLayerCurvilinear.f90
SOURCES += $(LDIR)/UpdateGhostLayerECurvilinear.f90
SOURCES += $(LDIR)/UpdateGhostLayerNCurvilinear.f90
SOURCES += $(LDIR)/UpdateGhostLayerX.f90
SOURCES += $(LDIR)/UpdateGhostLayerXCurvilinear.f90
SOURCES += $(LDIR)/UpdateGhostLayerY.f90
SOURCES += $(LDIR)/UpdateGhostLayerYCurvilinear.f90
SOURCES += $(LDIR)/VerticalFreeSurfaceVelocity.f90
SOURCES += $(LDIR)/Zero_Corners.f90
SOURCES += $(LDIR)/bessj.f90
SOURCES += $(LDIR)/eta_linear.f90
SOURCES += $(LDIR)/eta_linear_3D.f90
SOURCES += $(LDIR)/eta_linear_3D_bottom.f90
SOURCES += $(LDIR)/eta_linear_3D_standing.f90
SOURCES += $(LDIR)/eta_linear_bottom.f90
SOURCES += $(LDIR)/eta_linear_standing.f90
SOURCES += $(LDIR)/eta_nonlinear.f90
SOURCES += $(LDIR)/eta_nonlinear_3D.f90
SOURCES += $(LDIR)/eta_nonlinear_3D_standing.f90
SOURCES += $(LDIR)/eta_nonlinear_standing.f90
SOURCES += $(LDIR)/incident_linear_wf_finite.f90
SOURCES += $(LDIR)/incident_linear_wf_finite_bottom.f90
SOURCES += $(LDIR)/incident_linear_wf_finite_standing.f90
SOURCES += $(LDIR)/incident_nonlinear_wf_finite_standing.f90
SOURCES += $(LDIR)/incident_wf_finite.f90
SOURCES += $(LDIR)/maxnorm.f90
SOURCES += $(LDIR)/phi_linear.f90
SOURCES += $(LDIR)/phi_linear_3D.f90
SOURCES += $(LDIR)/phi_linear_3D_bottom.f90
SOURCES += $(LDIR)/phi_linear_3D_standing.f90
SOURCES += $(LDIR)/phi_linear_bottom.f90
SOURCES += $(LDIR)/phi_linear_standing.f90
SOURCES += $(LDIR)/phi_nonlinear.f90
SOURCES += $(LDIR)/phi_nonlinear_3D.f90
SOURCES += $(LDIR)/phi_nonlinear_3D_standing.f90
SOURCES += $(LDIR)/phi_nonlinear_standing.f90
SOURCES += $(LDIR)/rhsFreeSurface3D_SWENSE.f90
SOURCES += $(LDIR)/rhsLinearFreeSurface3D_SWENSE.f90