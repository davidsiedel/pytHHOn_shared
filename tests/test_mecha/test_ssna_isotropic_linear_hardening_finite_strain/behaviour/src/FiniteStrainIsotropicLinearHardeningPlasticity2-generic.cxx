/*!
* \file   FiniteStrainIsotropicLinearHardeningPlasticity2-generic.cxx
* \brief  This file implements the umat interface for the FiniteStrainIsotropicLinearHardeningPlasticity2 behaviour law
* \author Thomas Helfer
* \date   09 / 03 / 2021
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
#include"TFEL/Material/FiniteStrainBehaviourTangentOperator.hxx"
#include"TFEL/Material/FiniteStrainIsotropicLinearHardeningPlasticity2.hxx"
#include"MFront/GenericBehaviour/Integrate.hxx"

#include"MFront/GenericBehaviour/FiniteStrainIsotropicLinearHardeningPlasticity2-generic.hxx"

static tfel::material::OutOfBoundsPolicy&
FiniteStrainIsotropicLinearHardeningPlasticity2_getOutOfBoundsPolicy(){
using namespace tfel::material;
static OutOfBoundsPolicy policy = None;
return policy;
}

#ifdef __cplusplus
extern "C"{
#endif /* __cplusplus */

MFRONT_SHAREDOBJ const char* 
FiniteStrainIsotropicLinearHardeningPlasticity2_build_id = "";

MFRONT_SHAREDOBJ const char* 
FiniteStrainIsotropicLinearHardeningPlasticity2_mfront_ept = "FiniteStrainIsotropicLinearHardeningPlasticity2";

MFRONT_SHAREDOBJ const char* 
FiniteStrainIsotropicLinearHardeningPlasticity2_tfel_version = "4.0.0-dev";

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_mfront_mkt = 1u;

MFRONT_SHAREDOBJ const char *
FiniteStrainIsotropicLinearHardeningPlasticity2_mfront_interface = "Generic";

MFRONT_SHAREDOBJ const char *
FiniteStrainIsotropicLinearHardeningPlasticity2_src = "isotropic_linear_hardening_finite_strain_3.mfront";

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_nModellingHypotheses = 7u;

MFRONT_SHAREDOBJ const char * 
FiniteStrainIsotropicLinearHardeningPlasticity2_ModellingHypotheses[7u] = {"AxisymmetricalGeneralisedPlaneStrain",
"AxisymmetricalGeneralisedPlaneStress",
"Axisymmetrical",
"PlaneStress",
"PlaneStrain",
"GeneralisedPlaneStrain",
"Tridimensional"};

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_nMainVariables = 1;
MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_nGradients = 1;

MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_GradientsTypes[1] = {3};
MFRONT_SHAREDOBJ const char * FiniteStrainIsotropicLinearHardeningPlasticity2_Gradients[1] = {"DeformationGradient"};
MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_nThermodynamicForces = 1;

MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_ThermodynamicForcesTypes[1] = {1};
MFRONT_SHAREDOBJ const char * FiniteStrainIsotropicLinearHardeningPlasticity2_ThermodynamicForces[1] = {"Stress"};
MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_nTangentOperatorBlocks = 2u;

MFRONT_SHAREDOBJ const char * FiniteStrainIsotropicLinearHardeningPlasticity2_TangentOperatorBlocks[2] = {"Stress",
"DeformationGradient"};
MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_BehaviourType = 2u;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_BehaviourKinematic = 3u;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_SymmetryType = 0u;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_ElasticSymmetryType = 0u;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_UsableInPurelyImplicitResolution = 1;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_nMaterialProperties = 0u;

MFRONT_SHAREDOBJ const char * const *FiniteStrainIsotropicLinearHardeningPlasticity2_MaterialProperties = nullptr;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_nInternalStateVariables = 2;
MFRONT_SHAREDOBJ const char * FiniteStrainIsotropicLinearHardeningPlasticity2_InternalStateVariables[2] = {"ElasticStrain",
"EquivalentPlasticStrain"};
MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_InternalStateVariablesTypes [] = {1,0};

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_nExternalStateVariables = 0;
MFRONT_SHAREDOBJ const char * const * FiniteStrainIsotropicLinearHardeningPlasticity2_ExternalStateVariables = nullptr;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_nParameters = 11;
MFRONT_SHAREDOBJ const char * FiniteStrainIsotropicLinearHardeningPlasticity2_Parameters[11] = {"epsilon",
"theta","YoungModulus","PoissonRatio","RelativeValueForTheEquivalentStressLowerBoundDefinition","ihr_R0_",
"ihr_H_","minimal_time_step_scaling_factor","maximal_time_step_scaling_factor","numerical_jacobian_epsilon","iterMax"};
MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_ParametersTypes [] = {0,0,0,0,0,0,0,0,0,0,2};

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_epsilon_ParameterDefaultValue = 1e-14;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_theta_ParameterDefaultValue = 1;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_YoungModulus_ParameterDefaultValue = 70000000000;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_PoissonRatio_ParameterDefaultValue = 0.34;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_RelativeValueForTheEquivalentStressLowerBoundDefinition_ParameterDefaultValue = 1e-12;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_ihr_R0__ParameterDefaultValue = 300000000;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_ihr_H__ParameterDefaultValue = 10000000000;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_minimal_time_step_scaling_factor_ParameterDefaultValue = 0.1;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_maximal_time_step_scaling_factor_ParameterDefaultValue = 1.7976931348623e+308;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_numerical_jacobian_epsilon_ParameterDefaultValue = 1e-15;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_iterMax_ParameterDefaultValue  = 100;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_requiresStiffnessTensor = 0;
MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_requiresThermalExpansionCoefficientTensor = 0;
MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_nInitializeFunctions= 0;

MFRONT_SHAREDOBJ const char * const * FiniteStrainIsotropicLinearHardeningPlasticity2_InitializeFunctions = nullptr;


MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_ComputesInternalEnergy = 0;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_ComputesDissipatedEnergy = 0;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_UsableInPurelyImplicitResolution = 1;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_nMaterialProperties = 0u;

MFRONT_SHAREDOBJ const char * const *FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_MaterialProperties = nullptr;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_nInternalStateVariables = 3;
MFRONT_SHAREDOBJ const char * FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_InternalStateVariables[3] = {"ElasticStrain",
"EquivalentPlasticStrain","AxialStrain"};
MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_InternalStateVariablesTypes [] = {1,0,0};

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_nExternalStateVariables = 1;
MFRONT_SHAREDOBJ const char * FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_ExternalStateVariables[1] = {"AxialStress"};
MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_nParameters = 11;
MFRONT_SHAREDOBJ const char * FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_Parameters[11] = {"epsilon",
"theta","YoungModulus","PoissonRatio","RelativeValueForTheEquivalentStressLowerBoundDefinition","ihr_R0_",
"ihr_H_","minimal_time_step_scaling_factor","maximal_time_step_scaling_factor","numerical_jacobian_epsilon","iterMax"};
MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_ParametersTypes [] = {0,0,0,0,0,0,0,0,0,0,2};

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_epsilon_ParameterDefaultValue = 1e-14;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_theta_ParameterDefaultValue = 1;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_YoungModulus_ParameterDefaultValue = 70000000000;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_PoissonRatio_ParameterDefaultValue = 0.34;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_RelativeValueForTheEquivalentStressLowerBoundDefinition_ParameterDefaultValue = 1e-12;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_ihr_R0__ParameterDefaultValue = 300000000;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_ihr_H__ParameterDefaultValue = 10000000000;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_minimal_time_step_scaling_factor_ParameterDefaultValue = 0.1;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_maximal_time_step_scaling_factor_ParameterDefaultValue = 1.7976931348623e+308;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_numerical_jacobian_epsilon_ParameterDefaultValue = 1e-15;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_iterMax_ParameterDefaultValue  = 100;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_requiresStiffnessTensor = 0;
MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_requiresThermalExpansionCoefficientTensor = 0;
MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_nInitializeFunctions= 0;

MFRONT_SHAREDOBJ const char * const * FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_InitializeFunctions = nullptr;


MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_ComputesInternalEnergy = 0;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_ComputesDissipatedEnergy = 0;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_UsableInPurelyImplicitResolution = 1;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_nMaterialProperties = 0u;

MFRONT_SHAREDOBJ const char * const *FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_MaterialProperties = nullptr;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_nInternalStateVariables = 3;
MFRONT_SHAREDOBJ const char * FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_InternalStateVariables[3] = {"ElasticStrain",
"EquivalentPlasticStrain","AxialStrain"};
MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_InternalStateVariablesTypes [] = {1,0,0};

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_nExternalStateVariables = 0;
MFRONT_SHAREDOBJ const char * const * FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_ExternalStateVariables = nullptr;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_nParameters = 11;
MFRONT_SHAREDOBJ const char * FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_Parameters[11] = {"epsilon",
"theta","YoungModulus","PoissonRatio","RelativeValueForTheEquivalentStressLowerBoundDefinition","ihr_R0_",
"ihr_H_","minimal_time_step_scaling_factor","maximal_time_step_scaling_factor","numerical_jacobian_epsilon","iterMax"};
MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_ParametersTypes [] = {0,0,0,0,0,0,0,0,0,0,2};

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_epsilon_ParameterDefaultValue = 1e-14;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_theta_ParameterDefaultValue = 1;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_YoungModulus_ParameterDefaultValue = 70000000000;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_PoissonRatio_ParameterDefaultValue = 0.34;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_RelativeValueForTheEquivalentStressLowerBoundDefinition_ParameterDefaultValue = 1e-12;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_ihr_R0__ParameterDefaultValue = 300000000;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_ihr_H__ParameterDefaultValue = 10000000000;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_minimal_time_step_scaling_factor_ParameterDefaultValue = 0.1;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_maximal_time_step_scaling_factor_ParameterDefaultValue = 1.7976931348623e+308;

MFRONT_SHAREDOBJ double FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_numerical_jacobian_epsilon_ParameterDefaultValue = 1e-15;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_iterMax_ParameterDefaultValue  = 100;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_requiresStiffnessTensor = 0;
MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_requiresThermalExpansionCoefficientTensor = 0;
MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_nInitializeFunctions= 0;

MFRONT_SHAREDOBJ const char * const * FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_InitializeFunctions = nullptr;


MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_ComputesInternalEnergy = 0;

MFRONT_SHAREDOBJ unsigned short FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_ComputesDissipatedEnergy = 0;

MFRONT_SHAREDOBJ void
FiniteStrainIsotropicLinearHardeningPlasticity2_setOutOfBoundsPolicy(const int p){
if(p==0){
FiniteStrainIsotropicLinearHardeningPlasticity2_getOutOfBoundsPolicy() = tfel::material::None;
} else if(p==1){
FiniteStrainIsotropicLinearHardeningPlasticity2_getOutOfBoundsPolicy() = tfel::material::Warning;
} else if(p==2){
FiniteStrainIsotropicLinearHardeningPlasticity2_getOutOfBoundsPolicy() = tfel::material::Strict;
} else {
std::cerr << "FiniteStrainIsotropicLinearHardeningPlasticity2_setOutOfBoundsPolicy: invalid argument\n";
}
}

MFRONT_SHAREDOBJ int
FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_setParameter(const char *const key,const double value){
using tfel::material::FiniteStrainIsotropicLinearHardeningPlasticity2AxisymmetricalGeneralisedPlaneStressParametersInitializer;
auto& i = FiniteStrainIsotropicLinearHardeningPlasticity2AxisymmetricalGeneralisedPlaneStressParametersInitializer::get();
try{
i.set(key,value);
} catch(std::runtime_error& e){
std::cerr << e.what() << std::endl;
return 0;
}
return 1;
}

MFRONT_SHAREDOBJ int
FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_setUnsignedShortParameter(const char *const key,const unsigned short value){
using tfel::material::FiniteStrainIsotropicLinearHardeningPlasticity2AxisymmetricalGeneralisedPlaneStressParametersInitializer;
auto& i = FiniteStrainIsotropicLinearHardeningPlasticity2AxisymmetricalGeneralisedPlaneStressParametersInitializer::get();
try{
i.set(key,value);
} catch(std::runtime_error& e){
std::cerr << e.what() << std::endl;
return 0;
}
return 1;
}

MFRONT_SHAREDOBJ int
FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_setParameter(const char *const key,const double value){
using tfel::material::FiniteStrainIsotropicLinearHardeningPlasticity2PlaneStressParametersInitializer;
auto& i = FiniteStrainIsotropicLinearHardeningPlasticity2PlaneStressParametersInitializer::get();
try{
i.set(key,value);
} catch(std::runtime_error& e){
std::cerr << e.what() << std::endl;
return 0;
}
return 1;
}

MFRONT_SHAREDOBJ int
FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_setUnsignedShortParameter(const char *const key,const unsigned short value){
using tfel::material::FiniteStrainIsotropicLinearHardeningPlasticity2PlaneStressParametersInitializer;
auto& i = FiniteStrainIsotropicLinearHardeningPlasticity2PlaneStressParametersInitializer::get();
try{
i.set(key,value);
} catch(std::runtime_error& e){
std::cerr << e.what() << std::endl;
return 0;
}
return 1;
}

MFRONT_SHAREDOBJ int
FiniteStrainIsotropicLinearHardeningPlasticity2_setParameter(const char *const key,const double value){
using tfel::material::FiniteStrainIsotropicLinearHardeningPlasticity2ParametersInitializer;
auto& i = FiniteStrainIsotropicLinearHardeningPlasticity2ParametersInitializer::get();
try{
i.set(key,value);
} catch(std::runtime_error& e){
std::cerr << e.what() << std::endl;
return 0;
}
return 1;
}

MFRONT_SHAREDOBJ int
FiniteStrainIsotropicLinearHardeningPlasticity2_setUnsignedShortParameter(const char *const key,const unsigned short value){
using tfel::material::FiniteStrainIsotropicLinearHardeningPlasticity2ParametersInitializer;
auto& i = FiniteStrainIsotropicLinearHardeningPlasticity2ParametersInitializer::get();
try{
i.set(key,value);
} catch(std::runtime_error& e){
std::cerr << e.what() << std::endl;
return 0;
}
return 1;
}

MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using TangentOperator = FiniteStrainBehaviourTangentOperatorBase;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN;
using Behaviour = FiniteStrainIsotropicLinearHardeningPlasticity2<h,real,false>;
// stress measure 
enum struct StressMeasure { PK1, PK2, CAUCHY };
const auto sm = [&d]{
  if(d->K[1]<0.5){
    return StressMeasure::CAUCHY;
  } else if (d->K[1]<1.5){
    return StressMeasure::PK2;
  } else if (d->K[1]<2.5){
    return StressMeasure::PK1;
  } else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
  }
}();
// stiffness type
const auto smf = [&d]{
  if((d->K[0]>-0.5)&&(d->K[0]<0.5)){
    // no stiffness requested, 
    // returned value is meaningless
    return TangentOperator::DSIG_DF;
  }
  if(d->K[2]<0.5){
    return TangentOperator::DSIG_DF;
  } else if (d->K[2]<1.5){
    return TangentOperator::DS_DEGL;
  } else if (d->K[2]<2.5){
    return TangentOperator::DPK1_DF;
  } else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
  }
}();
tfel::math::st2tost2<1,real> K;
tfel::math::tensor<1,real> F0;
tfel::math::tensor<1,real> F1;
tfel::math::stensor<1,real> s0;
tfel::fsalgo::copy<3>::exe(d->s0.gradients,F0.begin());
tfel::fsalgo::copy<3>::exe(d->s1.gradients,F1.begin());
const auto setting = (smf==TangentOperator::DSIG_DF) ? 
LogarithmicStrainHandler<1,real>::EULERIAN :
LogarithmicStrainHandler<1,real>::LAGRANGIAN;
LogarithmicStrainHandler<1,real> lgh0(setting,F0);
LogarithmicStrainHandler<1,real> lgh1(setting,F1);
auto e0 = lgh0.getHenckyLogarithmicStrain();
auto e1 = lgh1.getHenckyLogarithmicStrain();
auto T0 = tfel::math::stensor<1,real>{};
auto T1 = tfel::math::stensor<1,real>{};
if (sm == StressMeasure::CAUCHY) {
tfel::fsalgo::copy<3>::exe(d->s0.thermodynamic_forces,s0.begin());
T0 = lgh0.convertFromCauchyStress(s0);
} else if (sm == StressMeasure::PK1) {
auto pk0 = tfel::math::tensor<1,real>{};
tfel::fsalgo::copy<3>::exe(d->s0.thermodynamic_forces,pk0.begin());
s0 = tfel::math::convertFirstPiolaKirchhoffStressToCauchyStress(pk0,F0);
T0 = lgh0.convertFromCauchyStress(s0);
} else if (sm == StressMeasure::PK2) {
auto S0 = tfel::math::stensor<1,real>{};
tfel::fsalgo::copy<3>::exe(d->s0.thermodynamic_forces,S0.begin());
s0 = convertSecondPiolaKirchhoffStressToCauchyStress(S0,F0);
T0 = lgh0.convertFromCauchyStress(s0);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
auto *const gradients0_old = d->s0.gradients;
auto *const gradients1_old = d->s1.gradients;
auto *const thermodynamic_forces0_old = d->s0.thermodynamic_forces;
auto *const thermodynamic_forces1_old = d->s1.thermodynamic_forces;
auto *const K_old = d->K;
K(0,0) = d->K[0];
d->s0.gradients = e0.begin();
d->s1.gradients = e1.begin();
d->s0.thermodynamic_forces = T0.begin();
d->s1.thermodynamic_forces = T1.begin();
d->K = K.begin();
const auto bp = K(0,0) < -0.5;
const auto bk = K(0,0) > 0.5;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, FiniteStrainIsotropicLinearHardeningPlasticity2_getOutOfBoundsPolicy());
d->s0.gradients = gradients0_old;
d->s1.gradients = gradients1_old;
d->s0.thermodynamic_forces = thermodynamic_forces0_old;
d->s1.thermodynamic_forces = thermodynamic_forces1_old;
d->K = K_old;
if(r){
if(bp){
if(smf==TangentOperator::DSIG_DF){
const auto Cs = lgh0.convertToSpatialTangentModuli(K,T0);
const auto Dt = convert<TangentOperator::DTAU_DF,TangentOperator::SPATIAL_MODULI>(Cs,F0,F0,s0);
tfel::math::T2toST2View<1,real>(d->K) = convert<TangentOperator::DSIG_DF,        TangentOperator::DTAU_DF>(Dt,F0,F0,s0);
} else if(smf==TangentOperator::DS_DEGL){
tfel::math::ST2toST2View<1,real>(d->K) = lgh0.convertToMaterialTangentModuli(K,T0);
} else if(smf==TangentOperator::DPK1_DF){
const auto Cse = lgh0.convertToMaterialTangentModuli(K,T0);
tfel::math::T2toT2View<1,real>(d->K) = convert<TangentOperator::DPK1_DF,        TangentOperator::DS_DEGL>(Cse,F0,F0,s0);
} else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
}
} else { // if(bp)
const auto s1 = lgh1.convertToCauchyStress(T1);
if(sm==StressMeasure::CAUCHY){
tfel::fsalgo::copy<3>::exe(s1.begin(),d->s1.thermodynamic_forces);
} else if(sm==StressMeasure::PK1){
tfel::math::TensorView<1,real> pk1(d->s1.thermodynamic_forces);
pk1 = tfel::math::convertCauchyStressToFirstPiolaKirchhoffStress(s1,F1);
} else if(sm==StressMeasure::PK2){
tfel::math::StensorView<1,real> S1(d->s1.thermodynamic_forces);
S1 = tfel::math::convertCauchyStressToSecondPiolaKirchhoffStress(s1,F1);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
if(bk){
if(smf==TangentOperator::DSIG_DF){
const auto Cs = lgh1.convertToSpatialTangentModuli(K,T1);
const auto Dt = convert<TangentOperator::DTAU_DF,                        TangentOperator::SPATIAL_MODULI>(Cs,F0,F1,s1);
tfel::math::T2toST2View<1,real>(d->K) = convert<TangentOperator::DSIG_DF,        TangentOperator::DTAU_DF>(Dt,F0,F1,s1);
} else if(smf==TangentOperator::DS_DEGL){
tfel::math::ST2toST2View<1,real>(d->K) = lgh1.convertToMaterialTangentModuli(K,T1);
} else if(smf==TangentOperator::DPK1_DF){
const auto Cse = lgh1.convertToMaterialTangentModuli(K,T1);
tfel::math::T2toT2View<1,real>(d->K) = convert<TangentOperator::DPK1_DF,        TangentOperator::DS_DEGL>(Cse,F0,F1,s1);
} else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
}
} // end of if(bk)
} // end of if(bp)
}
return r;
} // end of FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStrain

MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using TangentOperator = FiniteStrainBehaviourTangentOperatorBase;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS;
using Behaviour = FiniteStrainIsotropicLinearHardeningPlasticity2<h,real,false>;
// stress measure 
enum struct StressMeasure { PK1, PK2, CAUCHY };
const auto sm = [&d]{
  if(d->K[1]<0.5){
    return StressMeasure::CAUCHY;
  } else if (d->K[1]<1.5){
    return StressMeasure::PK2;
  } else if (d->K[1]<2.5){
    return StressMeasure::PK1;
  } else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
  }
}();
// stiffness type
const auto smf = [&d]{
  if((d->K[0]>-0.5)&&(d->K[0]<0.5)){
    // no stiffness requested, 
    // returned value is meaningless
    return TangentOperator::DSIG_DF;
  }
  if(d->K[2]<0.5){
    return TangentOperator::DSIG_DF;
  } else if (d->K[2]<1.5){
    return TangentOperator::DS_DEGL;
  } else if (d->K[2]<2.5){
    return TangentOperator::DPK1_DF;
  } else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
  }
}();
tfel::math::st2tost2<1,real> K;
tfel::math::tensor<1,real> F0;
tfel::math::tensor<1,real> F1;
tfel::math::stensor<1,real> s0;
tfel::fsalgo::copy<3>::exe(d->s0.gradients,F0.begin());
tfel::fsalgo::copy<3>::exe(d->s1.gradients,F1.begin());
const auto setting = (smf==TangentOperator::DSIG_DF) ? 
LogarithmicStrainHandler<1,real>::EULERIAN :
LogarithmicStrainHandler<1,real>::LAGRANGIAN;
LogarithmicStrainHandler<1,real> lgh0(setting,F0);
LogarithmicStrainHandler<1,real> lgh1(setting,F1);
auto e0 = lgh0.getHenckyLogarithmicStrain();
auto e1 = lgh1.getHenckyLogarithmicStrain();
lgh0.updateAxialDeformationGradient(std::exp(d->s0.internal_state_variables[4]));
auto T0 = tfel::math::stensor<1,real>{};
auto T1 = tfel::math::stensor<1,real>{};
if (sm == StressMeasure::CAUCHY) {
tfel::fsalgo::copy<3>::exe(d->s0.thermodynamic_forces,s0.begin());
T0 = lgh0.convertFromCauchyStress(s0);
} else if (sm == StressMeasure::PK1) {
auto pk0 = tfel::math::tensor<1,real>{};
tfel::fsalgo::copy<3>::exe(d->s0.thermodynamic_forces,pk0.begin());
s0 = tfel::math::convertFirstPiolaKirchhoffStressToCauchyStress(pk0,F0);
T0 = lgh0.convertFromCauchyStress(s0);
} else if (sm == StressMeasure::PK2) {
auto S0 = tfel::math::stensor<1,real>{};
tfel::fsalgo::copy<3>::exe(d->s0.thermodynamic_forces,S0.begin());
s0 = convertSecondPiolaKirchhoffStressToCauchyStress(S0,F0);
T0 = lgh0.convertFromCauchyStress(s0);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
auto *const gradients0_old = d->s0.gradients;
auto *const gradients1_old = d->s1.gradients;
auto *const thermodynamic_forces0_old = d->s0.thermodynamic_forces;
auto *const thermodynamic_forces1_old = d->s1.thermodynamic_forces;
auto *const K_old = d->K;
K(0,0) = d->K[0];
d->s0.gradients = e0.begin();
d->s1.gradients = e1.begin();
d->s0.thermodynamic_forces = T0.begin();
d->s1.thermodynamic_forces = T1.begin();
d->K = K.begin();
const auto bp = K(0,0) < -0.5;
const auto bk = K(0,0) > 0.5;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, FiniteStrainIsotropicLinearHardeningPlasticity2_getOutOfBoundsPolicy());
d->s0.gradients = gradients0_old;
d->s1.gradients = gradients1_old;
d->s0.thermodynamic_forces = thermodynamic_forces0_old;
d->s1.thermodynamic_forces = thermodynamic_forces1_old;
d->K = K_old;
if(r){
if(bp){
if(smf==TangentOperator::DSIG_DF){
const auto Cs = lgh0.convertToSpatialTangentModuli(K,T0);
const auto Dt = convert<TangentOperator::DTAU_DF,TangentOperator::SPATIAL_MODULI>(Cs,F0,F0,s0);
tfel::math::T2toST2View<1,real>(d->K) = convert<TangentOperator::DSIG_DF,        TangentOperator::DTAU_DF>(Dt,F0,F0,s0);
} else if(smf==TangentOperator::DS_DEGL){
tfel::math::ST2toST2View<1,real>(d->K) = lgh0.convertToMaterialTangentModuli(K,T0);
} else if(smf==TangentOperator::DPK1_DF){
const auto Cse = lgh0.convertToMaterialTangentModuli(K,T0);
tfel::math::T2toT2View<1,real>(d->K) = convert<TangentOperator::DPK1_DF,        TangentOperator::DS_DEGL>(Cse,F0,F0,s0);
} else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
}
} else { // if(bp)
lgh1.updateAxialDeformationGradient(std::exp(d->s1.internal_state_variables[4]));
const auto s1 = lgh1.convertToCauchyStress(T1);
if(sm==StressMeasure::CAUCHY){
tfel::fsalgo::copy<3>::exe(s1.begin(),d->s1.thermodynamic_forces);
} else if(sm==StressMeasure::PK1){
tfel::math::TensorView<1,real> pk1(d->s1.thermodynamic_forces);
pk1 = tfel::math::convertCauchyStressToFirstPiolaKirchhoffStress(s1,F1);
} else if(sm==StressMeasure::PK2){
tfel::math::StensorView<1,real> S1(d->s1.thermodynamic_forces);
S1 = tfel::math::convertCauchyStressToSecondPiolaKirchhoffStress(s1,F1);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
if(bk){
if(smf==TangentOperator::DSIG_DF){
const auto Cs = lgh1.convertToSpatialTangentModuli(K,T1);
const auto Dt = convert<TangentOperator::DTAU_DF,                        TangentOperator::SPATIAL_MODULI>(Cs,F0,F1,s1);
tfel::math::T2toST2View<1,real>(d->K) = convert<TangentOperator::DSIG_DF,        TangentOperator::DTAU_DF>(Dt,F0,F1,s1);
} else if(smf==TangentOperator::DS_DEGL){
tfel::math::ST2toST2View<1,real>(d->K) = lgh1.convertToMaterialTangentModuli(K,T1);
} else if(smf==TangentOperator::DPK1_DF){
const auto Cse = lgh1.convertToMaterialTangentModuli(K,T1);
tfel::math::T2toT2View<1,real>(d->K) = convert<TangentOperator::DPK1_DF,        TangentOperator::DS_DEGL>(Cse,F0,F1,s1);
} else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
}
} // end of if(bk)
} // end of if(bp)
}
return r;
} // end of FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress

MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_Axisymmetrical(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using TangentOperator = FiniteStrainBehaviourTangentOperatorBase;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::AXISYMMETRICAL;
using Behaviour = FiniteStrainIsotropicLinearHardeningPlasticity2<h,real,false>;
// stress measure 
enum struct StressMeasure { PK1, PK2, CAUCHY };
const auto sm = [&d]{
  if(d->K[1]<0.5){
    return StressMeasure::CAUCHY;
  } else if (d->K[1]<1.5){
    return StressMeasure::PK2;
  } else if (d->K[1]<2.5){
    return StressMeasure::PK1;
  } else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
  }
}();
// stiffness type
const auto smf = [&d]{
  if((d->K[0]>-0.5)&&(d->K[0]<0.5)){
    // no stiffness requested, 
    // returned value is meaningless
    return TangentOperator::DSIG_DF;
  }
  if(d->K[2]<0.5){
    return TangentOperator::DSIG_DF;
  } else if (d->K[2]<1.5){
    return TangentOperator::DS_DEGL;
  } else if (d->K[2]<2.5){
    return TangentOperator::DPK1_DF;
  } else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
  }
}();
tfel::math::st2tost2<2,real> K;
tfel::math::tensor<2,real> F0;
tfel::math::tensor<2,real> F1;
tfel::math::stensor<2,real> s0;
tfel::fsalgo::copy<5>::exe(d->s0.gradients,F0.begin());
tfel::fsalgo::copy<5>::exe(d->s1.gradients,F1.begin());
const auto setting = (smf==TangentOperator::DSIG_DF) ? 
LogarithmicStrainHandler<2,real>::EULERIAN :
LogarithmicStrainHandler<2,real>::LAGRANGIAN;
LogarithmicStrainHandler<2,real> lgh0(setting,F0);
LogarithmicStrainHandler<2,real> lgh1(setting,F1);
auto e0 = lgh0.getHenckyLogarithmicStrain();
auto e1 = lgh1.getHenckyLogarithmicStrain();
auto T0 = tfel::math::stensor<2,real>{};
auto T1 = tfel::math::stensor<2,real>{};
if (sm == StressMeasure::CAUCHY) {
tfel::fsalgo::copy<4>::exe(d->s0.thermodynamic_forces,s0.begin());
T0 = lgh0.convertFromCauchyStress(s0);
} else if (sm == StressMeasure::PK1) {
auto pk0 = tfel::math::tensor<2,real>{};
tfel::fsalgo::copy<5>::exe(d->s0.thermodynamic_forces,pk0.begin());
s0 = tfel::math::convertFirstPiolaKirchhoffStressToCauchyStress(pk0,F0);
T0 = lgh0.convertFromCauchyStress(s0);
} else if (sm == StressMeasure::PK2) {
auto S0 = tfel::math::stensor<2,real>{};
tfel::fsalgo::copy<4>::exe(d->s0.thermodynamic_forces,S0.begin());
s0 = convertSecondPiolaKirchhoffStressToCauchyStress(S0,F0);
T0 = lgh0.convertFromCauchyStress(s0);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
auto *const gradients0_old = d->s0.gradients;
auto *const gradients1_old = d->s1.gradients;
auto *const thermodynamic_forces0_old = d->s0.thermodynamic_forces;
auto *const thermodynamic_forces1_old = d->s1.thermodynamic_forces;
auto *const K_old = d->K;
K(0,0) = d->K[0];
d->s0.gradients = e0.begin();
d->s1.gradients = e1.begin();
d->s0.thermodynamic_forces = T0.begin();
d->s1.thermodynamic_forces = T1.begin();
d->K = K.begin();
const auto bp = K(0,0) < -0.5;
const auto bk = K(0,0) > 0.5;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, FiniteStrainIsotropicLinearHardeningPlasticity2_getOutOfBoundsPolicy());
d->s0.gradients = gradients0_old;
d->s1.gradients = gradients1_old;
d->s0.thermodynamic_forces = thermodynamic_forces0_old;
d->s1.thermodynamic_forces = thermodynamic_forces1_old;
d->K = K_old;
if(r){
if(bp){
if(smf==TangentOperator::DSIG_DF){
const auto Cs = lgh0.convertToSpatialTangentModuli(K,T0);
const auto Dt = convert<TangentOperator::DTAU_DF,TangentOperator::SPATIAL_MODULI>(Cs,F0,F0,s0);
tfel::math::T2toST2View<2,real>(d->K) = convert<TangentOperator::DSIG_DF,        TangentOperator::DTAU_DF>(Dt,F0,F0,s0);
} else if(smf==TangentOperator::DS_DEGL){
tfel::math::ST2toST2View<2,real>(d->K) = lgh0.convertToMaterialTangentModuli(K,T0);
} else if(smf==TangentOperator::DPK1_DF){
const auto Cse = lgh0.convertToMaterialTangentModuli(K,T0);
tfel::math::T2toT2View<2,real>(d->K) = convert<TangentOperator::DPK1_DF,        TangentOperator::DS_DEGL>(Cse,F0,F0,s0);
} else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
}
} else { // if(bp)
const auto s1 = lgh1.convertToCauchyStress(T1);
if(sm==StressMeasure::CAUCHY){
tfel::fsalgo::copy<4>::exe(s1.begin(),d->s1.thermodynamic_forces);
} else if(sm==StressMeasure::PK1){
tfel::math::TensorView<2,real> pk1(d->s1.thermodynamic_forces);
pk1 = tfel::math::convertCauchyStressToFirstPiolaKirchhoffStress(s1,F1);
} else if(sm==StressMeasure::PK2){
tfel::math::StensorView<2,real> S1(d->s1.thermodynamic_forces);
S1 = tfel::math::convertCauchyStressToSecondPiolaKirchhoffStress(s1,F1);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
if(bk){
if(smf==TangentOperator::DSIG_DF){
const auto Cs = lgh1.convertToSpatialTangentModuli(K,T1);
const auto Dt = convert<TangentOperator::DTAU_DF,                        TangentOperator::SPATIAL_MODULI>(Cs,F0,F1,s1);
tfel::math::T2toST2View<2,real>(d->K) = convert<TangentOperator::DSIG_DF,        TangentOperator::DTAU_DF>(Dt,F0,F1,s1);
} else if(smf==TangentOperator::DS_DEGL){
tfel::math::ST2toST2View<2,real>(d->K) = lgh1.convertToMaterialTangentModuli(K,T1);
} else if(smf==TangentOperator::DPK1_DF){
const auto Cse = lgh1.convertToMaterialTangentModuli(K,T1);
tfel::math::T2toT2View<2,real>(d->K) = convert<TangentOperator::DPK1_DF,        TangentOperator::DS_DEGL>(Cse,F0,F1,s1);
} else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
}
} // end of if(bk)
} // end of if(bp)
}
return r;
} // end of FiniteStrainIsotropicLinearHardeningPlasticity2_Axisymmetrical

MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using TangentOperator = FiniteStrainBehaviourTangentOperatorBase;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::PLANESTRESS;
using Behaviour = FiniteStrainIsotropicLinearHardeningPlasticity2<h,real,false>;
// stress measure 
enum struct StressMeasure { PK1, PK2, CAUCHY };
const auto sm = [&d]{
  if(d->K[1]<0.5){
    return StressMeasure::CAUCHY;
  } else if (d->K[1]<1.5){
    return StressMeasure::PK2;
  } else if (d->K[1]<2.5){
    return StressMeasure::PK1;
  } else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
  }
}();
// stiffness type
const auto smf = [&d]{
  if((d->K[0]>-0.5)&&(d->K[0]<0.5)){
    // no stiffness requested, 
    // returned value is meaningless
    return TangentOperator::DSIG_DF;
  }
  if(d->K[2]<0.5){
    return TangentOperator::DSIG_DF;
  } else if (d->K[2]<1.5){
    return TangentOperator::DS_DEGL;
  } else if (d->K[2]<2.5){
    return TangentOperator::DPK1_DF;
  } else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
  }
}();
tfel::math::st2tost2<2,real> K;
tfel::math::tensor<2,real> F0;
tfel::math::tensor<2,real> F1;
tfel::math::stensor<2,real> s0;
tfel::fsalgo::copy<5>::exe(d->s0.gradients,F0.begin());
tfel::fsalgo::copy<5>::exe(d->s1.gradients,F1.begin());
const auto setting = (smf==TangentOperator::DSIG_DF) ? 
LogarithmicStrainHandler<2,real>::EULERIAN :
LogarithmicStrainHandler<2,real>::LAGRANGIAN;
LogarithmicStrainHandler<2,real> lgh0(setting,F0);
LogarithmicStrainHandler<2,real> lgh1(setting,F1);
auto e0 = lgh0.getHenckyLogarithmicStrain();
auto e1 = lgh1.getHenckyLogarithmicStrain();
lgh0.updateAxialDeformationGradient(std::exp(d->s0.internal_state_variables[5]));
auto T0 = tfel::math::stensor<2,real>{};
auto T1 = tfel::math::stensor<2,real>{};
if (sm == StressMeasure::CAUCHY) {
tfel::fsalgo::copy<4>::exe(d->s0.thermodynamic_forces,s0.begin());
T0 = lgh0.convertFromCauchyStress(s0);
} else if (sm == StressMeasure::PK1) {
auto pk0 = tfel::math::tensor<2,real>{};
tfel::fsalgo::copy<5>::exe(d->s0.thermodynamic_forces,pk0.begin());
s0 = tfel::math::convertFirstPiolaKirchhoffStressToCauchyStress(pk0,F0);
T0 = lgh0.convertFromCauchyStress(s0);
} else if (sm == StressMeasure::PK2) {
auto S0 = tfel::math::stensor<2,real>{};
tfel::fsalgo::copy<4>::exe(d->s0.thermodynamic_forces,S0.begin());
s0 = convertSecondPiolaKirchhoffStressToCauchyStress(S0,F0);
T0 = lgh0.convertFromCauchyStress(s0);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
auto *const gradients0_old = d->s0.gradients;
auto *const gradients1_old = d->s1.gradients;
auto *const thermodynamic_forces0_old = d->s0.thermodynamic_forces;
auto *const thermodynamic_forces1_old = d->s1.thermodynamic_forces;
auto *const K_old = d->K;
K(0,0) = d->K[0];
d->s0.gradients = e0.begin();
d->s1.gradients = e1.begin();
d->s0.thermodynamic_forces = T0.begin();
d->s1.thermodynamic_forces = T1.begin();
d->K = K.begin();
const auto bp = K(0,0) < -0.5;
const auto bk = K(0,0) > 0.5;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, FiniteStrainIsotropicLinearHardeningPlasticity2_getOutOfBoundsPolicy());
d->s0.gradients = gradients0_old;
d->s1.gradients = gradients1_old;
d->s0.thermodynamic_forces = thermodynamic_forces0_old;
d->s1.thermodynamic_forces = thermodynamic_forces1_old;
d->K = K_old;
if(r){
if(bp){
if(smf==TangentOperator::DSIG_DF){
const auto Cs = lgh0.convertToSpatialTangentModuli(K,T0);
const auto Dt = convert<TangentOperator::DTAU_DF,TangentOperator::SPATIAL_MODULI>(Cs,F0,F0,s0);
tfel::math::T2toST2View<2,real>(d->K) = convert<TangentOperator::DSIG_DF,        TangentOperator::DTAU_DF>(Dt,F0,F0,s0);
} else if(smf==TangentOperator::DS_DEGL){
tfel::math::ST2toST2View<2,real>(d->K) = lgh0.convertToMaterialTangentModuli(K,T0);
} else if(smf==TangentOperator::DPK1_DF){
const auto Cse = lgh0.convertToMaterialTangentModuli(K,T0);
tfel::math::T2toT2View<2,real>(d->K) = convert<TangentOperator::DPK1_DF,        TangentOperator::DS_DEGL>(Cse,F0,F0,s0);
} else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
}
} else { // if(bp)
lgh1.updateAxialDeformationGradient(std::exp(d->s1.internal_state_variables[5]));
const auto s1 = lgh1.convertToCauchyStress(T1);
if(sm==StressMeasure::CAUCHY){
tfel::fsalgo::copy<4>::exe(s1.begin(),d->s1.thermodynamic_forces);
} else if(sm==StressMeasure::PK1){
tfel::math::TensorView<2,real> pk1(d->s1.thermodynamic_forces);
pk1 = tfel::math::convertCauchyStressToFirstPiolaKirchhoffStress(s1,F1);
} else if(sm==StressMeasure::PK2){
tfel::math::StensorView<2,real> S1(d->s1.thermodynamic_forces);
S1 = tfel::math::convertCauchyStressToSecondPiolaKirchhoffStress(s1,F1);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
if(bk){
if(smf==TangentOperator::DSIG_DF){
const auto Cs = lgh1.convertToSpatialTangentModuli(K,T1);
const auto Dt = convert<TangentOperator::DTAU_DF,                        TangentOperator::SPATIAL_MODULI>(Cs,F0,F1,s1);
tfel::math::T2toST2View<2,real>(d->K) = convert<TangentOperator::DSIG_DF,        TangentOperator::DTAU_DF>(Dt,F0,F1,s1);
} else if(smf==TangentOperator::DS_DEGL){
tfel::math::ST2toST2View<2,real>(d->K) = lgh1.convertToMaterialTangentModuli(K,T1);
} else if(smf==TangentOperator::DPK1_DF){
const auto Cse = lgh1.convertToMaterialTangentModuli(K,T1);
tfel::math::T2toT2View<2,real>(d->K) = convert<TangentOperator::DPK1_DF,        TangentOperator::DS_DEGL>(Cse,F0,F1,s1);
} else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
}
} // end of if(bk)
} // end of if(bp)
}
return r;
} // end of FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress

MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using TangentOperator = FiniteStrainBehaviourTangentOperatorBase;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::PLANESTRAIN;
using Behaviour = FiniteStrainIsotropicLinearHardeningPlasticity2<h,real,false>;
// stress measure 
enum struct StressMeasure { PK1, PK2, CAUCHY };
const auto sm = [&d]{
  if(d->K[1]<0.5){
    return StressMeasure::CAUCHY;
  } else if (d->K[1]<1.5){
    return StressMeasure::PK2;
  } else if (d->K[1]<2.5){
    return StressMeasure::PK1;
  } else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
  }
}();
// stiffness type
const auto smf = [&d]{
  if((d->K[0]>-0.5)&&(d->K[0]<0.5)){
    // no stiffness requested, 
    // returned value is meaningless
    return TangentOperator::DSIG_DF;
  }
  if(d->K[2]<0.5){
    return TangentOperator::DSIG_DF;
  } else if (d->K[2]<1.5){
    return TangentOperator::DS_DEGL;
  } else if (d->K[2]<2.5){
    return TangentOperator::DPK1_DF;
  } else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
  }
}();
tfel::math::st2tost2<2,real> K;
tfel::math::tensor<2,real> F0;
tfel::math::tensor<2,real> F1;
tfel::math::stensor<2,real> s0;
tfel::fsalgo::copy<5>::exe(d->s0.gradients,F0.begin());
tfel::fsalgo::copy<5>::exe(d->s1.gradients,F1.begin());
const auto setting = (smf==TangentOperator::DSIG_DF) ? 
LogarithmicStrainHandler<2,real>::EULERIAN :
LogarithmicStrainHandler<2,real>::LAGRANGIAN;
LogarithmicStrainHandler<2,real> lgh0(setting,F0);
LogarithmicStrainHandler<2,real> lgh1(setting,F1);
auto e0 = lgh0.getHenckyLogarithmicStrain();
auto e1 = lgh1.getHenckyLogarithmicStrain();
auto T0 = tfel::math::stensor<2,real>{};
auto T1 = tfel::math::stensor<2,real>{};
if (sm == StressMeasure::CAUCHY) {
tfel::fsalgo::copy<4>::exe(d->s0.thermodynamic_forces,s0.begin());
T0 = lgh0.convertFromCauchyStress(s0);
} else if (sm == StressMeasure::PK1) {
auto pk0 = tfel::math::tensor<2,real>{};
tfel::fsalgo::copy<5>::exe(d->s0.thermodynamic_forces,pk0.begin());
s0 = tfel::math::convertFirstPiolaKirchhoffStressToCauchyStress(pk0,F0);
T0 = lgh0.convertFromCauchyStress(s0);
} else if (sm == StressMeasure::PK2) {
auto S0 = tfel::math::stensor<2,real>{};
tfel::fsalgo::copy<4>::exe(d->s0.thermodynamic_forces,S0.begin());
s0 = convertSecondPiolaKirchhoffStressToCauchyStress(S0,F0);
T0 = lgh0.convertFromCauchyStress(s0);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
auto *const gradients0_old = d->s0.gradients;
auto *const gradients1_old = d->s1.gradients;
auto *const thermodynamic_forces0_old = d->s0.thermodynamic_forces;
auto *const thermodynamic_forces1_old = d->s1.thermodynamic_forces;
auto *const K_old = d->K;
K(0,0) = d->K[0];
d->s0.gradients = e0.begin();
d->s1.gradients = e1.begin();
d->s0.thermodynamic_forces = T0.begin();
d->s1.thermodynamic_forces = T1.begin();
d->K = K.begin();
const auto bp = K(0,0) < -0.5;
const auto bk = K(0,0) > 0.5;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, FiniteStrainIsotropicLinearHardeningPlasticity2_getOutOfBoundsPolicy());
d->s0.gradients = gradients0_old;
d->s1.gradients = gradients1_old;
d->s0.thermodynamic_forces = thermodynamic_forces0_old;
d->s1.thermodynamic_forces = thermodynamic_forces1_old;
d->K = K_old;
if(r){
if(bp){
if(smf==TangentOperator::DSIG_DF){
const auto Cs = lgh0.convertToSpatialTangentModuli(K,T0);
const auto Dt = convert<TangentOperator::DTAU_DF,TangentOperator::SPATIAL_MODULI>(Cs,F0,F0,s0);
tfel::math::T2toST2View<2,real>(d->K) = convert<TangentOperator::DSIG_DF,        TangentOperator::DTAU_DF>(Dt,F0,F0,s0);
} else if(smf==TangentOperator::DS_DEGL){
tfel::math::ST2toST2View<2,real>(d->K) = lgh0.convertToMaterialTangentModuli(K,T0);
} else if(smf==TangentOperator::DPK1_DF){
const auto Cse = lgh0.convertToMaterialTangentModuli(K,T0);
tfel::math::T2toT2View<2,real>(d->K) = convert<TangentOperator::DPK1_DF,        TangentOperator::DS_DEGL>(Cse,F0,F0,s0);
} else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
}
} else { // if(bp)
const auto s1 = lgh1.convertToCauchyStress(T1);
if(sm==StressMeasure::CAUCHY){
tfel::fsalgo::copy<4>::exe(s1.begin(),d->s1.thermodynamic_forces);
} else if(sm==StressMeasure::PK1){
tfel::math::TensorView<2,real> pk1(d->s1.thermodynamic_forces);
pk1 = tfel::math::convertCauchyStressToFirstPiolaKirchhoffStress(s1,F1);
} else if(sm==StressMeasure::PK2){
tfel::math::StensorView<2,real> S1(d->s1.thermodynamic_forces);
S1 = tfel::math::convertCauchyStressToSecondPiolaKirchhoffStress(s1,F1);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
if(bk){
if(smf==TangentOperator::DSIG_DF){
const auto Cs = lgh1.convertToSpatialTangentModuli(K,T1);
const auto Dt = convert<TangentOperator::DTAU_DF,                        TangentOperator::SPATIAL_MODULI>(Cs,F0,F1,s1);
tfel::math::T2toST2View<2,real>(d->K) = convert<TangentOperator::DSIG_DF,        TangentOperator::DTAU_DF>(Dt,F0,F1,s1);
} else if(smf==TangentOperator::DS_DEGL){
tfel::math::ST2toST2View<2,real>(d->K) = lgh1.convertToMaterialTangentModuli(K,T1);
} else if(smf==TangentOperator::DPK1_DF){
const auto Cse = lgh1.convertToMaterialTangentModuli(K,T1);
tfel::math::T2toT2View<2,real>(d->K) = convert<TangentOperator::DPK1_DF,        TangentOperator::DS_DEGL>(Cse,F0,F1,s1);
} else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
}
} // end of if(bk)
} // end of if(bp)
}
return r;
} // end of FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStrain

MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_GeneralisedPlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using TangentOperator = FiniteStrainBehaviourTangentOperatorBase;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::GENERALISEDPLANESTRAIN;
using Behaviour = FiniteStrainIsotropicLinearHardeningPlasticity2<h,real,false>;
// stress measure 
enum struct StressMeasure { PK1, PK2, CAUCHY };
const auto sm = [&d]{
  if(d->K[1]<0.5){
    return StressMeasure::CAUCHY;
  } else if (d->K[1]<1.5){
    return StressMeasure::PK2;
  } else if (d->K[1]<2.5){
    return StressMeasure::PK1;
  } else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
  }
}();
// stiffness type
const auto smf = [&d]{
  if((d->K[0]>-0.5)&&(d->K[0]<0.5)){
    // no stiffness requested, 
    // returned value is meaningless
    return TangentOperator::DSIG_DF;
  }
  if(d->K[2]<0.5){
    return TangentOperator::DSIG_DF;
  } else if (d->K[2]<1.5){
    return TangentOperator::DS_DEGL;
  } else if (d->K[2]<2.5){
    return TangentOperator::DPK1_DF;
  } else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
  }
}();
tfel::math::st2tost2<2,real> K;
tfel::math::tensor<2,real> F0;
tfel::math::tensor<2,real> F1;
tfel::math::stensor<2,real> s0;
tfel::fsalgo::copy<5>::exe(d->s0.gradients,F0.begin());
tfel::fsalgo::copy<5>::exe(d->s1.gradients,F1.begin());
const auto setting = (smf==TangentOperator::DSIG_DF) ? 
LogarithmicStrainHandler<2,real>::EULERIAN :
LogarithmicStrainHandler<2,real>::LAGRANGIAN;
LogarithmicStrainHandler<2,real> lgh0(setting,F0);
LogarithmicStrainHandler<2,real> lgh1(setting,F1);
auto e0 = lgh0.getHenckyLogarithmicStrain();
auto e1 = lgh1.getHenckyLogarithmicStrain();
auto T0 = tfel::math::stensor<2,real>{};
auto T1 = tfel::math::stensor<2,real>{};
if (sm == StressMeasure::CAUCHY) {
tfel::fsalgo::copy<4>::exe(d->s0.thermodynamic_forces,s0.begin());
T0 = lgh0.convertFromCauchyStress(s0);
} else if (sm == StressMeasure::PK1) {
auto pk0 = tfel::math::tensor<2,real>{};
tfel::fsalgo::copy<5>::exe(d->s0.thermodynamic_forces,pk0.begin());
s0 = tfel::math::convertFirstPiolaKirchhoffStressToCauchyStress(pk0,F0);
T0 = lgh0.convertFromCauchyStress(s0);
} else if (sm == StressMeasure::PK2) {
auto S0 = tfel::math::stensor<2,real>{};
tfel::fsalgo::copy<4>::exe(d->s0.thermodynamic_forces,S0.begin());
s0 = convertSecondPiolaKirchhoffStressToCauchyStress(S0,F0);
T0 = lgh0.convertFromCauchyStress(s0);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
auto *const gradients0_old = d->s0.gradients;
auto *const gradients1_old = d->s1.gradients;
auto *const thermodynamic_forces0_old = d->s0.thermodynamic_forces;
auto *const thermodynamic_forces1_old = d->s1.thermodynamic_forces;
auto *const K_old = d->K;
K(0,0) = d->K[0];
d->s0.gradients = e0.begin();
d->s1.gradients = e1.begin();
d->s0.thermodynamic_forces = T0.begin();
d->s1.thermodynamic_forces = T1.begin();
d->K = K.begin();
const auto bp = K(0,0) < -0.5;
const auto bk = K(0,0) > 0.5;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, FiniteStrainIsotropicLinearHardeningPlasticity2_getOutOfBoundsPolicy());
d->s0.gradients = gradients0_old;
d->s1.gradients = gradients1_old;
d->s0.thermodynamic_forces = thermodynamic_forces0_old;
d->s1.thermodynamic_forces = thermodynamic_forces1_old;
d->K = K_old;
if(r){
if(bp){
if(smf==TangentOperator::DSIG_DF){
const auto Cs = lgh0.convertToSpatialTangentModuli(K,T0);
const auto Dt = convert<TangentOperator::DTAU_DF,TangentOperator::SPATIAL_MODULI>(Cs,F0,F0,s0);
tfel::math::T2toST2View<2,real>(d->K) = convert<TangentOperator::DSIG_DF,        TangentOperator::DTAU_DF>(Dt,F0,F0,s0);
} else if(smf==TangentOperator::DS_DEGL){
tfel::math::ST2toST2View<2,real>(d->K) = lgh0.convertToMaterialTangentModuli(K,T0);
} else if(smf==TangentOperator::DPK1_DF){
const auto Cse = lgh0.convertToMaterialTangentModuli(K,T0);
tfel::math::T2toT2View<2,real>(d->K) = convert<TangentOperator::DPK1_DF,        TangentOperator::DS_DEGL>(Cse,F0,F0,s0);
} else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
}
} else { // if(bp)
const auto s1 = lgh1.convertToCauchyStress(T1);
if(sm==StressMeasure::CAUCHY){
tfel::fsalgo::copy<4>::exe(s1.begin(),d->s1.thermodynamic_forces);
} else if(sm==StressMeasure::PK1){
tfel::math::TensorView<2,real> pk1(d->s1.thermodynamic_forces);
pk1 = tfel::math::convertCauchyStressToFirstPiolaKirchhoffStress(s1,F1);
} else if(sm==StressMeasure::PK2){
tfel::math::StensorView<2,real> S1(d->s1.thermodynamic_forces);
S1 = tfel::math::convertCauchyStressToSecondPiolaKirchhoffStress(s1,F1);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
if(bk){
if(smf==TangentOperator::DSIG_DF){
const auto Cs = lgh1.convertToSpatialTangentModuli(K,T1);
const auto Dt = convert<TangentOperator::DTAU_DF,                        TangentOperator::SPATIAL_MODULI>(Cs,F0,F1,s1);
tfel::math::T2toST2View<2,real>(d->K) = convert<TangentOperator::DSIG_DF,        TangentOperator::DTAU_DF>(Dt,F0,F1,s1);
} else if(smf==TangentOperator::DS_DEGL){
tfel::math::ST2toST2View<2,real>(d->K) = lgh1.convertToMaterialTangentModuli(K,T1);
} else if(smf==TangentOperator::DPK1_DF){
const auto Cse = lgh1.convertToMaterialTangentModuli(K,T1);
tfel::math::T2toT2View<2,real>(d->K) = convert<TangentOperator::DPK1_DF,        TangentOperator::DS_DEGL>(Cse,F0,F1,s1);
} else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
}
} // end of if(bk)
} // end of if(bp)
}
return r;
} // end of FiniteStrainIsotropicLinearHardeningPlasticity2_GeneralisedPlaneStrain

MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_Tridimensional(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using TangentOperator = FiniteStrainBehaviourTangentOperatorBase;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::TRIDIMENSIONAL;
using Behaviour = FiniteStrainIsotropicLinearHardeningPlasticity2<h,real,false>;
// stress measure 
enum struct StressMeasure { PK1, PK2, CAUCHY };
const auto sm = [&d]{
  if(d->K[1]<0.5){
    return StressMeasure::CAUCHY;
  } else if (d->K[1]<1.5){
    return StressMeasure::PK2;
  } else if (d->K[1]<2.5){
    return StressMeasure::PK1;
  } else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
  }
}();
// stiffness type
const auto smf = [&d]{
  if((d->K[0]>-0.5)&&(d->K[0]<0.5)){
    // no stiffness requested, 
    // returned value is meaningless
    return TangentOperator::DSIG_DF;
  }
  if(d->K[2]<0.5){
    return TangentOperator::DSIG_DF;
  } else if (d->K[2]<1.5){
    return TangentOperator::DS_DEGL;
  } else if (d->K[2]<2.5){
    return TangentOperator::DPK1_DF;
  } else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
  }
}();
tfel::math::st2tost2<3,real> K;
tfel::math::tensor<3,real> F0;
tfel::math::tensor<3,real> F1;
tfel::math::stensor<3,real> s0;
tfel::fsalgo::copy<9>::exe(d->s0.gradients,F0.begin());
tfel::fsalgo::copy<9>::exe(d->s1.gradients,F1.begin());
const auto setting = (smf==TangentOperator::DSIG_DF) ? 
LogarithmicStrainHandler<3,real>::EULERIAN :
LogarithmicStrainHandler<3,real>::LAGRANGIAN;
LogarithmicStrainHandler<3,real> lgh0(setting,F0);
LogarithmicStrainHandler<3,real> lgh1(setting,F1);
auto e0 = lgh0.getHenckyLogarithmicStrain();
auto e1 = lgh1.getHenckyLogarithmicStrain();
auto T0 = tfel::math::stensor<3,real>{};
auto T1 = tfel::math::stensor<3,real>{};
if (sm == StressMeasure::CAUCHY) {
tfel::fsalgo::copy<6>::exe(d->s0.thermodynamic_forces,s0.begin());
T0 = lgh0.convertFromCauchyStress(s0);
} else if (sm == StressMeasure::PK1) {
auto pk0 = tfel::math::tensor<3,real>{};
tfel::fsalgo::copy<9>::exe(d->s0.thermodynamic_forces,pk0.begin());
s0 = tfel::math::convertFirstPiolaKirchhoffStressToCauchyStress(pk0,F0);
T0 = lgh0.convertFromCauchyStress(s0);
} else if (sm == StressMeasure::PK2) {
auto S0 = tfel::math::stensor<3,real>{};
tfel::fsalgo::copy<6>::exe(d->s0.thermodynamic_forces,S0.begin());
s0 = convertSecondPiolaKirchhoffStressToCauchyStress(S0,F0);
T0 = lgh0.convertFromCauchyStress(s0);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
auto *const gradients0_old = d->s0.gradients;
auto *const gradients1_old = d->s1.gradients;
auto *const thermodynamic_forces0_old = d->s0.thermodynamic_forces;
auto *const thermodynamic_forces1_old = d->s1.thermodynamic_forces;
auto *const K_old = d->K;
K(0,0) = d->K[0];
d->s0.gradients = e0.begin();
d->s1.gradients = e1.begin();
d->s0.thermodynamic_forces = T0.begin();
d->s1.thermodynamic_forces = T1.begin();
d->K = K.begin();
const auto bp = K(0,0) < -0.5;
const auto bk = K(0,0) > 0.5;
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, FiniteStrainIsotropicLinearHardeningPlasticity2_getOutOfBoundsPolicy());
d->s0.gradients = gradients0_old;
d->s1.gradients = gradients1_old;
d->s0.thermodynamic_forces = thermodynamic_forces0_old;
d->s1.thermodynamic_forces = thermodynamic_forces1_old;
d->K = K_old;
if(r){
if(bp){
if(smf==TangentOperator::DSIG_DF){
const auto Cs = lgh0.convertToSpatialTangentModuli(K,T0);
const auto Dt = convert<TangentOperator::DTAU_DF,TangentOperator::SPATIAL_MODULI>(Cs,F0,F0,s0);
tfel::math::T2toST2View<3,real>(d->K) = convert<TangentOperator::DSIG_DF,        TangentOperator::DTAU_DF>(Dt,F0,F0,s0);
} else if(smf==TangentOperator::DS_DEGL){
tfel::math::ST2toST2View<3,real>(d->K) = lgh0.convertToMaterialTangentModuli(K,T0);
} else if(smf==TangentOperator::DPK1_DF){
const auto Cse = lgh0.convertToMaterialTangentModuli(K,T0);
tfel::math::T2toT2View<3,real>(d->K) = convert<TangentOperator::DPK1_DF,        TangentOperator::DS_DEGL>(Cse,F0,F0,s0);
} else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
}
} else { // if(bp)
const auto s1 = lgh1.convertToCauchyStress(T1);
if(sm==StressMeasure::CAUCHY){
tfel::fsalgo::copy<6>::exe(s1.begin(),d->s1.thermodynamic_forces);
} else if(sm==StressMeasure::PK1){
tfel::math::TensorView<3,real> pk1(d->s1.thermodynamic_forces);
pk1 = tfel::math::convertCauchyStressToFirstPiolaKirchhoffStress(s1,F1);
} else if(sm==StressMeasure::PK2){
tfel::math::StensorView<3,real> S1(d->s1.thermodynamic_forces);
S1 = tfel::math::convertCauchyStressToSecondPiolaKirchhoffStress(s1,F1);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
if(bk){
if(smf==TangentOperator::DSIG_DF){
const auto Cs = lgh1.convertToSpatialTangentModuli(K,T1);
const auto Dt = convert<TangentOperator::DTAU_DF,                        TangentOperator::SPATIAL_MODULI>(Cs,F0,F1,s1);
tfel::math::T2toST2View<3,real>(d->K) = convert<TangentOperator::DSIG_DF,        TangentOperator::DTAU_DF>(Dt,F0,F1,s1);
} else if(smf==TangentOperator::DS_DEGL){
tfel::math::ST2toST2View<3,real>(d->K) = lgh1.convertToMaterialTangentModuli(K,T1);
} else if(smf==TangentOperator::DPK1_DF){
const auto Cse = lgh1.convertToMaterialTangentModuli(K,T1);
tfel::math::T2toT2View<3,real>(d->K) = convert<TangentOperator::DPK1_DF,        TangentOperator::DS_DEGL>(Cse,F0,F1,s1);
} else {
  std::cerr << "invalid choice for consistent tangent "
               "operator\n";
  std::exit(-1);
}
} // end of if(bk)
} // end of if(bp)
}
return r;
} // end of FiniteStrainIsotropicLinearHardeningPlasticity2_Tridimensional

#ifdef __cplusplus
}
#endif /* __cplusplus */

