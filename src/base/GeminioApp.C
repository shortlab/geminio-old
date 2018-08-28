#include "GeminioApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

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
