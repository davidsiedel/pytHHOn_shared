/*!
* \file   IsotropicLinearKinematicHardeningPlasticity-generic.cxx
* \brief  This file implements the umat interface for the IsotropicLinearKinematicHardeningPlasticity behaviour law
* \author Thomas Helfer
* \date   14 / 10 / 2016
*/

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif /* NOMINMAX */
#include <windows.h>
#ifdef small
#undef small
#endif /* small */
#endif /* _WIN32 */

#ifndef MFRONT_SHAREDOBJ
#define MFRONT_SHAREDOBJ TFEL_VISIBILITY_EXPORT
#endif /* MFRONT_SHAREDOBJ */

#include<iostream>
#include<cstdlib>
#include"TFEL/Material/OutOfBoundsPolicy.hxx"
#include"TFEL/Math/t2tot2.hxx"
#include"TFEL/Math/t2tost2.hxx"
#include"TFEL/Material/IsotropicLinearKinematicHardeningPlasticity.hxx"
#include"MFront/GenericBehaviour/Integrate.hxx"

#include"MFront/GenericBehaviour/IsotropicLinearKinematicHardeningPlasticity-generic.hxx"

static tfel::material::OutOfBoundsPolicy&
IsotropicLinearKinematicHardeningPlasticity_getOutOfBoundsPolicy(){
using namespace tfel::material;
static OutOfBoundsPolicy policy = None;
return policy;
}

#ifdef __cplusplus
extern "C"{
#endif /* __cplusplus */

MFRONT_SHAREDOBJ const char* 
IsotropicLinearKinematicHardeningPlasticity_build_id = "";

MFRONT_SHAREDOBJ const char* 
IsotropicLinearKinematicHardeningPlasticity_mfront_ept = "IsotropicLinearKinematicHardeningPlasticity";

MFRONT_SHAREDOBJ const char* 
IsotropicLinearKinematicHardeningPlasticity_tfel_version = "4.0.0-dev";

MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_mfront_mkt = 1u;

MFRONT_SHAREDOBJ const char *
IsotropicLinearKinematicHardeningPlasticity_mfront_interface = "Generic";

MFRONT_SHAREDOBJ const char *
IsotropicLinearKinematicHardeningPlasticity_src = "IsotropicLinearKinematicHardeningPlasticity.mfront";

MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_nModellingHypotheses = 5u;

MFRONT_SHAREDOBJ const char * 
IsotropicLinearKinematicHardeningPlasticity_ModellingHypotheses[5u] = {"AxisymmetricalGeneralisedPlaneStrain",
"Axisymmetrical",
"PlaneStrain",
"GeneralisedPlaneStrain",
"Tridimensional"};

MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_nMainVariables = 1;
MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_nGradients = 1;

MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_GradientsTypes[1] = {1};
MFRONT_SHAREDOBJ const char * IsotropicLinearKinematicHardeningPlasticity_Gradients[1] = {"Strain"};
MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_nThermodynamicForces = 1;

MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_ThermodynamicForcesTypes[1] = {1};
MFRONT_SHAREDOBJ const char * IsotropicLinearKinematicHardeningPlasticity_ThermodynamicForces[1] = {"Stress"};
MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_nTangentOperatorBlocks = 2;

MFRONT_SHAREDOBJ const char * IsotropicLinearKinematicHardeningPlasticity_TangentOperatorBlocks[2] = {"Stress",
"Strain"};
MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_BehaviourType = 1u;

MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_BehaviourKinematic = 1u;

MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_SymmetryType = 0u;

MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_ElasticSymmetryType = 0u;

MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_UsableInPurelyImplicitResolution = 1;

MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_nMaterialProperties = 0u;

MFRONT_SHAREDOBJ const char * const *IsotropicLinearKinematicHardeningPlasticity_MaterialProperties = nullptr;

MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_nInternalStateVariables = 2;
MFRONT_SHAREDOBJ const char * IsotropicLinearKinematicHardeningPlasticity_InternalStateVariables[2] = {"ElasticStrain",
"a"};
MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_InternalStateVariablesTypes [] = {1,1};

MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_nExternalStateVariables = 0;
MFRONT_SHAREDOBJ const char * const * IsotropicLinearKinematicHardeningPlasticity_ExternalStateVariables = nullptr;

MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_nParameters = 6;
MFRONT_SHAREDOBJ const char * IsotropicLinearKinematicHardeningPlasticity_Parameters[6] = {"YoungModulus",
"PoissonRatio","YieldStrength","HardeningSlope","minimal_time_step_scaling_factor","maximal_time_step_scaling_factor"};
MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_ParametersTypes [] = {0,0,0,0,0,0};

MFRONT_SHAREDOBJ double IsotropicLinearKinematicHardeningPlasticity_YoungModulus_ParameterDefaultValue = 70000000000;

MFRONT_SHAREDOBJ double IsotropicLinearKinematicHardeningPlasticity_PoissonRatio_ParameterDefaultValue = 0.34;

MFRONT_SHAREDOBJ double IsotropicLinearKinematicHardeningPlasticity_YieldStrength_ParameterDefaultValue = 300000000;

MFRONT_SHAREDOBJ double IsotropicLinearKinematicHardeningPlasticity_HardeningSlope_ParameterDefaultValue = 10000000000;

MFRONT_SHAREDOBJ double IsotropicLinearKinematicHardeningPlasticity_minimal_time_step_scaling_factor_ParameterDefaultValue = 0.1;

MFRONT_SHAREDOBJ double IsotropicLinearKinematicHardeningPlasticity_maximal_time_step_scaling_factor_ParameterDefaultValue = 1.7976931348623e+308;

MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_requiresStiffnessTensor = 0;
MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_requiresThermalExpansionCoefficientTensor = 0;
MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_nInitializeFunctions= 0;

MFRONT_SHAREDOBJ const char * const * IsotropicLinearKinematicHardeningPlasticity_InitializeFunctions = nullptr;


MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_ComputesInternalEnergy = 1;

MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_ComputesDissipatedEnergy = 1;

MFRONT_SHAREDOBJ void
IsotropicLinearKinematicHardeningPlasticity_setOutOfBoundsPolicy(const int p){
if(p==0){
IsotropicLinearKinematicHardeningPlasticity_getOutOfBoundsPolicy() = tfel::material::None;
} else if(p==1){
IsotropicLinearKinematicHardeningPlasticity_getOutOfBoundsPolicy() = tfel::material::Warning;
} else if(p==2){
IsotropicLinearKinematicHardeningPlasticity_getOutOfBoundsPolicy() = tfel::material::Strict;
} else {
std::cerr << "IsotropicLinearKinematicHardeningPlasticity_setOutOfBoundsPolicy: invalid argument\n";
}
}

MFRONT_SHAREDOBJ int
IsotropicLinearKinematicHardeningPlasticity_setParameter(const char *const key,const double value){
using tfel::material::IsotropicLinearKinematicHardeningPlasticityParametersInitializer;
auto& i = IsotropicLinearKinematicHardeningPlasticityParametersInitializer::get();
try{
i.set(key,value);
} catch(std::runtime_error& e){
std::cerr << e.what() << std::endl;
return 0;
}
return 1;
}

MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_AxisymmetricalGeneralisedPlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN;
using Behaviour = IsotropicLinearKinematicHardeningPlasticity<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, IsotropicLinearKinematicHardeningPlasticity_getOutOfBoundsPolicy());
return r;
} // end of IsotropicLinearKinematicHardeningPlasticity_AxisymmetricalGeneralisedPlaneStrain

MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_Axisymmetrical(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::AXISYMMETRICAL;
using Behaviour = IsotropicLinearKinematicHardeningPlasticity<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, IsotropicLinearKinematicHardeningPlasticity_getOutOfBoundsPolicy());
return r;
} // end of IsotropicLinearKinematicHardeningPlasticity_Axisymmetrical

MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_PlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::PLANESTRAIN;
using Behaviour = IsotropicLinearKinematicHardeningPlasticity<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, IsotropicLinearKinematicHardeningPlasticity_getOutOfBoundsPolicy());
return r;
} // end of IsotropicLinearKinematicHardeningPlasticity_PlaneStrain

MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_GeneralisedPlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::GENERALISEDPLANESTRAIN;
using Behaviour = IsotropicLinearKinematicHardeningPlasticity<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, IsotropicLinearKinematicHardeningPlasticity_getOutOfBoundsPolicy());
return r;
} // end of IsotropicLinearKinematicHardeningPlasticity_GeneralisedPlaneStrain

MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_Tridimensional(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::TRIDIMENSIONAL;
using Behaviour = IsotropicLinearKinematicHardeningPlasticity<h,real,false>;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, IsotropicLinearKinematicHardeningPlasticity_getOutOfBoundsPolicy());
return r;
} // end of IsotropicLinearKinematicHardeningPlasticity_Tridimensional

#ifdef __cplusplus
}
#endif /* __cplusplus */

