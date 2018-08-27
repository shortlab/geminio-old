#include "GeminioApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

// Include kernel header files
#include "DefectSink.h"
#include "DefectSource.h"
#include "DefectRecombination.h"
#include "MatPropDiffusion.h"
#include "FuncCoefVariable.h"
//for cluster
#include "VariableProduct.h"
#include "CoeffDiffusion.h"
#include "SingleVariable.h"
#include "DislocationSink.h"
#include "UserObjectVariableProduct.h"
#include "UserObjectSingleVariable.h"
#include "UserObjectDiffusion.h"
#include "MobileDefects.h"
#include "ImmobileDefects.h"

// Include auxkernel header files
#include "DislocationSinkRate.h"
#include "DefectRecombinationRateConstant.h"
#include "VoidSinkRate.h"
#include "GVoidSwelling.h"
#include "GSumSIAClusterDensity.h"
#include "ClusterDensity.h"


// Next come the materials
#include "RadiationMaterial.h"
#include "ArcMaterial.h"

//*************Postprocessors**********************//
#include "NodalConservationCheck.h"
#include "TotalDefectLoss.h"

//*************UserObjects**************************//
#include "GroupConstant.h"
#include "MaterialConstants.h"
#include "TestProperty.h"
#include "GroupingTest.h"
#include "BCCIronProperty.h"




/***************grouping method start*********************/
//###################kernels################//
#include "GImmobileL0.h"
#include "GImmobileL1.h"
#include "GMobile.h"
#include "ConstantKernel.h"
//#####################Actions##############//
#include "AddGVariable.h"
#include "AddGImmobile.h"
#include "AddGMobile.h"
#include "AddGTimeDerivative.h"
#include "AddGConstantKernels.h"
//#####################user objects##########//
#include "GGroup.h"
#include "GMaterialConstants.h"
#include "GGroupingTest.h"
#include "GIron.h"
#include "GTungsten.h"
/***************grouping method end*********************/


template<>
InputParameters validParams<GeminioApp>()
{
  InputParameters params = validParams<MooseApp>();
  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  return params;
}

GeminioApp::GeminioApp(const InputParameters & parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  GeminioApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  GeminioApp::associateSyntax(_syntax, _action_factory);

  Moose::registerExecFlags(_factory);
  GeminioApp::registerExecFlags(_factory);
}

GeminioApp::~GeminioApp() {}

// External entry point for dynamic application loading
extern "C" void
GeminioApp__registerApps()
{
  GeminioApp::registerApps();
}
void
GeminioApp::registerApps()
{
  registerApp(GeminioApp);
}

// External entry point for dynamic object registration
extern "C" void
GeminioApp__registerObjects(Factory & factory)
{
  GeminioApp::registerObjects(factory);
}
void
GeminioApp::registerObjects(Factory & factory)
{
  Registry::registerObjectsTo(factory, {"GeminioApp"});
}

// External entry point for dynamic syntax association
extern "C" void
GeminioApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  GeminioApp::associateSyntax(syntax, action_factory);
}
void
GeminioApp::associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  Registry::registerActionsTo(action_factory, {"GeminioApp"});

  // Action syntax
  syntax.registerActionSyntax("AddLotsOfVariableProduct", "LotsOfVariableProduct/*");
  syntax.registerActionSyntax("AddUserObjectVariableProduct", "LotsOfUserObjectVariableProduct/*");
  syntax.registerActionSyntax("AddLotsOfTimeDerivative", "LotsOfTimeDerivative/*");
  syntax.registerActionSyntax("AddLotsOfCoeffDiffusion", "LotsOfCoeffDiffusion/*");
  syntax.registerActionSyntax("AddLotsOfSingleVariable", "LotsOfSingleVariable/*");
  syntax.registerActionSyntax("AddUserObjectSingleVariable", "LotsOfUserObjectSingleVariable/*");
  syntax.registerActionSyntax("AddUserObjectDiffusion", "LotsOfUserObjectDiffusion/*");
  syntax.registerActionSyntax("AddUserObjectDislocationSink", "LotsOfUserObjectDislocationSink/*");
  syntax.registerActionSyntax("AddLotsOfSource", "LotsOfSource/*");
  syntax.registerActionSyntax("AddClusterICAction", "ClusterIC/*");
  syntax.registerActionSyntax("AddLotsOfDislocationSinks", "LotsOfSink_disl/*"); // deprecated
  syntax.registerActionSyntax("AddLotsOfDislocationSinks", "LotsOfDislocationSinks/*");
  syntax.registerActionSyntax("AddLotsOfVariableAction", "LotsOfVariables/*");
  syntax.registerActionSyntax("AddMobileDefects", "MobileDefects/*");
  syntax.registerActionSyntax("AddImmobileDefects", "ImmobileDefects/*");
  syntax.registerActionSyntax("AddClusterDensity", "ClusterDensity/*");

  // Grouping method actions
  syntax.registerActionSyntax("AddGVariable","GVariable/*");
  syntax.registerActionSyntax("AddGMobile","GMobile/*");
  syntax.registerActionSyntax("AddGImmobile","GImmobile/*");
  syntax.registerActionSyntax("AddGTimeDerivative", "GTimeDerivative/*");
  syntax.registerActionSyntax("AddGConstantKernels", "Sources/*");
  syntax.registerActionSyntax("AddGVoidSwelling", "GVoidSwelling/*");
  syntax.registerActionSyntax("AddGSumSIAClusterDensity", "GSumSIAClusterDensity/*");
}

// External entry point for dynamic execute flag registration
extern "C" void
GeminioApp__registerExecFlags(Factory & factory)
{
  GeminioApp::registerExecFlags(factory);
}
void
GeminioApp::registerExecFlags(Factory & /*factory*/)
{
}
