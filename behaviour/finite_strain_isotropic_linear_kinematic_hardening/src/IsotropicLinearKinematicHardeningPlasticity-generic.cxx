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
#include"TFEL/Material/FiniteStrainBehaviourTangentOperator.hxx"
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
IsotropicLinearKinematicHardeningPlasticity_src = "FiniteStrainIsotropicLinearKinematicHardeningPlasticity.mfront";

MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_nModellingHypotheses = 5u;

MFRONT_SHAREDOBJ const char * 
IsotropicLinearKinematicHardeningPlasticity_ModellingHypotheses[5u] = {"AxisymmetricalGeneralisedPlaneStrain",
"Axisymmetrical",
"PlaneStrain",
"GeneralisedPlaneStrain",
"Tridimensional"};

MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_nMainVariables = 1;
MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_nGradients = 1;

MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_GradientsTypes[1] = {3};
MFRONT_SHAREDOBJ const char * IsotropicLinearKinematicHardeningPlasticity_Gradients[1] = {"DeformationGradient"};
MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_nThermodynamicForces = 1;

MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_ThermodynamicForcesTypes[1] = {1};
MFRONT_SHAREDOBJ const char * IsotropicLinearKinematicHardeningPlasticity_ThermodynamicForces[1] = {"Stress"};
MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_nTangentOperatorBlocks = 2u;

MFRONT_SHAREDOBJ const char * IsotropicLinearKinematicHardeningPlasticity_TangentOperatorBlocks[2] = {"Stress",
"DeformationGradient"};
MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_BehaviourType = 2u;

MFRONT_SHAREDOBJ unsigned short IsotropicLinearKinematicHardeningPlasticity_BehaviourKinematic = 3u;

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
using TangentOperator = FiniteStrainBehaviourTangentOperatorBase;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN;
using Behaviour = IsotropicLinearKinematicHardeningPlasticity<h,real,false>;
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
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, IsotropicLinearKinematicHardeningPlasticity_getOutOfBoundsPolicy());
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
} // end of IsotropicLinearKinematicHardeningPlasticity_AxisymmetricalGeneralisedPlaneStrain

MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_Axisymmetrical(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using TangentOperator = FiniteStrainBehaviourTangentOperatorBase;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::AXISYMMETRICAL;
using Behaviour = IsotropicLinearKinematicHardeningPlasticity<h,real,false>;
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
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, IsotropicLinearKinematicHardeningPlasticity_getOutOfBoundsPolicy());
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
} // end of IsotropicLinearKinematicHardeningPlasticity_Axisymmetrical

MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_PlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using TangentOperator = FiniteStrainBehaviourTangentOperatorBase;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::PLANESTRAIN;
using Behaviour = IsotropicLinearKinematicHardeningPlasticity<h,real,false>;
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
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, IsotropicLinearKinematicHardeningPlasticity_getOutOfBoundsPolicy());
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
} // end of IsotropicLinearKinematicHardeningPlasticity_PlaneStrain

MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_GeneralisedPlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using TangentOperator = FiniteStrainBehaviourTangentOperatorBase;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::GENERALISEDPLANESTRAIN;
using Behaviour = IsotropicLinearKinematicHardeningPlasticity<h,real,false>;
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
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, IsotropicLinearKinematicHardeningPlasticity_getOutOfBoundsPolicy());
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
} // end of IsotropicLinearKinematicHardeningPlasticity_GeneralisedPlaneStrain

MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_Tridimensional(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using TangentOperator = FiniteStrainBehaviourTangentOperatorBase;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::TRIDIMENSIONAL;
using Behaviour = IsotropicLinearKinematicHardeningPlasticity<h,real,false>;
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
const auto r = mfront::gb::integrate<Behaviour>(*d,Behaviour::STANDARDTANGENTOPERATOR, IsotropicLinearKinematicHardeningPlasticity_getOutOfBoundsPolicy());
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
} // end of IsotropicLinearKinematicHardeningPlasticity_Tridimensional

#ifdef __cplusplus
}
#endif /* __cplusplus */

