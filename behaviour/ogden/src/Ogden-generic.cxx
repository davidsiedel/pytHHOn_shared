/*!
* \file   Ogden-generic.cxx
* \brief  This file implements the umat interface for the Ogden behaviour law
* \author Thomas Helfer
* \date   20 / 12 / 2016
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
#include"TFEL/Material/Ogden.hxx"
#include"MFront/GenericBehaviour/Integrate.hxx"

#include"MFront/GenericBehaviour/Ogden-generic.hxx"

static tfel::material::OutOfBoundsPolicy&
Ogden_getOutOfBoundsPolicy(){
using namespace tfel::material;
static OutOfBoundsPolicy policy = None;
return policy;
}

#ifdef __cplusplus
extern "C"{
#endif /* __cplusplus */

MFRONT_SHAREDOBJ const char* 
Ogden_build_id = "";

MFRONT_SHAREDOBJ const char* 
Ogden_mfront_ept = "Ogden";

MFRONT_SHAREDOBJ const char* 
Ogden_tfel_version = "4.0.0-dev";

MFRONT_SHAREDOBJ unsigned short Ogden_mfront_mkt = 1u;

MFRONT_SHAREDOBJ const char *
Ogden_mfront_interface = "Generic";

MFRONT_SHAREDOBJ const char *
Ogden_src = "ogden.mfront";

MFRONT_SHAREDOBJ unsigned short Ogden_nModellingHypotheses = 5u;

MFRONT_SHAREDOBJ const char * 
Ogden_ModellingHypotheses[5u] = {"AxisymmetricalGeneralisedPlaneStrain",
"Axisymmetrical",
"PlaneStrain",
"GeneralisedPlaneStrain",
"Tridimensional"};

MFRONT_SHAREDOBJ unsigned short Ogden_nMainVariables = 1;
MFRONT_SHAREDOBJ unsigned short Ogden_nGradients = 1;

MFRONT_SHAREDOBJ int Ogden_GradientsTypes[1] = {3};
MFRONT_SHAREDOBJ const char * Ogden_Gradients[1] = {"DeformationGradient"};
MFRONT_SHAREDOBJ unsigned short Ogden_nThermodynamicForces = 1;

MFRONT_SHAREDOBJ int Ogden_ThermodynamicForcesTypes[1] = {1};
MFRONT_SHAREDOBJ const char * Ogden_ThermodynamicForces[1] = {"Stress"};
MFRONT_SHAREDOBJ unsigned short Ogden_nTangentOperatorBlocks = 2u;

MFRONT_SHAREDOBJ const char * Ogden_TangentOperatorBlocks[2] = {"Stress",
"DeformationGradient"};
MFRONT_SHAREDOBJ unsigned short Ogden_BehaviourType = 2u;

MFRONT_SHAREDOBJ unsigned short Ogden_BehaviourKinematic = 3u;

MFRONT_SHAREDOBJ unsigned short Ogden_SymmetryType = 0u;

MFRONT_SHAREDOBJ unsigned short Ogden_ElasticSymmetryType = 0u;

MFRONT_SHAREDOBJ unsigned short Ogden_UsableInPurelyImplicitResolution = 1;

MFRONT_SHAREDOBJ unsigned short Ogden_nMaterialProperties = 0u;

MFRONT_SHAREDOBJ const char * const *Ogden_MaterialProperties = nullptr;

MFRONT_SHAREDOBJ unsigned short Ogden_nInternalStateVariables = 0;
MFRONT_SHAREDOBJ const char * const * Ogden_InternalStateVariables = nullptr;

MFRONT_SHAREDOBJ const int * Ogden_InternalStateVariablesTypes = nullptr;

MFRONT_SHAREDOBJ unsigned short Ogden_nExternalStateVariables = 0;
MFRONT_SHAREDOBJ const char * const * Ogden_ExternalStateVariables = nullptr;

MFRONT_SHAREDOBJ unsigned short Ogden_nParameters = 5;
MFRONT_SHAREDOBJ const char * Ogden_Parameters[5] = {"alpha",
"mu","K","minimal_time_step_scaling_factor","maximal_time_step_scaling_factor"};
MFRONT_SHAREDOBJ int Ogden_ParametersTypes [] = {0,0,0,0,0};

MFRONT_SHAREDOBJ double Ogden_alpha_ParameterDefaultValue = 28.8;

MFRONT_SHAREDOBJ double Ogden_mu_ParameterDefaultValue = 27778;

MFRONT_SHAREDOBJ double Ogden_K_ParameterDefaultValue = 69444444;

MFRONT_SHAREDOBJ double Ogden_minimal_time_step_scaling_factor_ParameterDefaultValue = 0.1;

MFRONT_SHAREDOBJ double Ogden_maximal_time_step_scaling_factor_ParameterDefaultValue = 1.7976931348623e+308;

MFRONT_SHAREDOBJ unsigned short Ogden_requiresStiffnessTensor = 0;
MFRONT_SHAREDOBJ unsigned short Ogden_requiresThermalExpansionCoefficientTensor = 0;
MFRONT_SHAREDOBJ unsigned short Ogden_nInitializeFunctions= 0;

MFRONT_SHAREDOBJ const char * const * Ogden_InitializeFunctions = nullptr;


MFRONT_SHAREDOBJ unsigned short Ogden_ComputesInternalEnergy = 0;

MFRONT_SHAREDOBJ unsigned short Ogden_ComputesDissipatedEnergy = 0;

MFRONT_SHAREDOBJ void
Ogden_setOutOfBoundsPolicy(const int p){
if(p==0){
Ogden_getOutOfBoundsPolicy() = tfel::material::None;
} else if(p==1){
Ogden_getOutOfBoundsPolicy() = tfel::material::Warning;
} else if(p==2){
Ogden_getOutOfBoundsPolicy() = tfel::material::Strict;
} else {
std::cerr << "Ogden_setOutOfBoundsPolicy: invalid argument\n";
}
}

MFRONT_SHAREDOBJ int
Ogden_setParameter(const char *const key,const double value){
using tfel::material::OgdenParametersInitializer;
auto& i = OgdenParametersInitializer::get();
try{
i.set(key,value);
} catch(std::runtime_error& e){
std::cerr << e.what() << std::endl;
return 0;
}
return 1;
}

MFRONT_SHAREDOBJ int Ogden_AxisymmetricalGeneralisedPlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using TangentOperator = FiniteStrainBehaviourTangentOperatorBase;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRAIN;
using Behaviour = Ogden<h,real,false>;
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
tfel::math::stensor<1,real> s0;
tfel::math::stensor<1,real> s1;
auto *const thermodynamic_forces0_old = d->s0.thermodynamic_forces;
auto *const thermodynamic_forces1_old = d->s1.thermodynamic_forces;
if(sm!=StressMeasure::CAUCHY){
tfel::math::tensor<1,real> F0;
tfel::fsalgo::copy<3>::exe(d->s0.gradients,F0.begin());
if(sm==StressMeasure::PK1){
tfel::math::TensorView<1,real> pk0(d->s0.thermodynamic_forces);
s0 = tfel::math::convertFirstPiolaKirchhoffStressToCauchyStress(pk0,F0);
} else if(sm==StressMeasure::PK2){
tfel::math::StensorView<1,real> S0(d->s0.thermodynamic_forces);
s0 = tfel::math::convertSecondPiolaKirchhoffStressToCauchyStress(S0,F0);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
d->s0.thermodynamic_forces = s0.begin();
d->s1.thermodynamic_forces = s1.begin();
}
const auto r = mfront::gb::integrate<Behaviour>(*d,smf, Ogden_getOutOfBoundsPolicy());
if((r) && (sm != StressMeasure::CAUCHY)){
tfel::math::tensor<1,real> F1;
tfel::fsalgo::copy<3>::exe(d->s1.gradients,F1.begin());
d->s0.thermodynamic_forces = thermodynamic_forces0_old;
d->s1.thermodynamic_forces = thermodynamic_forces1_old;
if(sm==StressMeasure::PK1){
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
}
if((!r) && (sm!=StressMeasure::CAUCHY)){
d->s0.thermodynamic_forces = thermodynamic_forces0_old;
d->s1.thermodynamic_forces = thermodynamic_forces1_old;
}
return r;
} // end of Ogden_AxisymmetricalGeneralisedPlaneStrain

MFRONT_SHAREDOBJ int Ogden_Axisymmetrical(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using TangentOperator = FiniteStrainBehaviourTangentOperatorBase;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::AXISYMMETRICAL;
using Behaviour = Ogden<h,real,false>;
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
tfel::math::stensor<2,real> s0;
tfel::math::stensor<2,real> s1;
auto *const thermodynamic_forces0_old = d->s0.thermodynamic_forces;
auto *const thermodynamic_forces1_old = d->s1.thermodynamic_forces;
if(sm!=StressMeasure::CAUCHY){
tfel::math::tensor<2,real> F0;
tfel::fsalgo::copy<5>::exe(d->s0.gradients,F0.begin());
if(sm==StressMeasure::PK1){
tfel::math::TensorView<2,real> pk0(d->s0.thermodynamic_forces);
s0 = tfel::math::convertFirstPiolaKirchhoffStressToCauchyStress(pk0,F0);
} else if(sm==StressMeasure::PK2){
tfel::math::StensorView<2,real> S0(d->s0.thermodynamic_forces);
s0 = tfel::math::convertSecondPiolaKirchhoffStressToCauchyStress(S0,F0);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
d->s0.thermodynamic_forces = s0.begin();
d->s1.thermodynamic_forces = s1.begin();
}
const auto r = mfront::gb::integrate<Behaviour>(*d,smf, Ogden_getOutOfBoundsPolicy());
if((r) && (sm != StressMeasure::CAUCHY)){
tfel::math::tensor<2,real> F1;
tfel::fsalgo::copy<5>::exe(d->s1.gradients,F1.begin());
d->s0.thermodynamic_forces = thermodynamic_forces0_old;
d->s1.thermodynamic_forces = thermodynamic_forces1_old;
if(sm==StressMeasure::PK1){
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
}
if((!r) && (sm!=StressMeasure::CAUCHY)){
d->s0.thermodynamic_forces = thermodynamic_forces0_old;
d->s1.thermodynamic_forces = thermodynamic_forces1_old;
}
return r;
} // end of Ogden_Axisymmetrical

MFRONT_SHAREDOBJ int Ogden_PlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using TangentOperator = FiniteStrainBehaviourTangentOperatorBase;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::PLANESTRAIN;
using Behaviour = Ogden<h,real,false>;
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
tfel::math::stensor<2,real> s0;
tfel::math::stensor<2,real> s1;
auto *const thermodynamic_forces0_old = d->s0.thermodynamic_forces;
auto *const thermodynamic_forces1_old = d->s1.thermodynamic_forces;
if(sm!=StressMeasure::CAUCHY){
tfel::math::tensor<2,real> F0;
tfel::fsalgo::copy<5>::exe(d->s0.gradients,F0.begin());
if(sm==StressMeasure::PK1){
tfel::math::TensorView<2,real> pk0(d->s0.thermodynamic_forces);
s0 = tfel::math::convertFirstPiolaKirchhoffStressToCauchyStress(pk0,F0);
} else if(sm==StressMeasure::PK2){
tfel::math::StensorView<2,real> S0(d->s0.thermodynamic_forces);
s0 = tfel::math::convertSecondPiolaKirchhoffStressToCauchyStress(S0,F0);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
d->s0.thermodynamic_forces = s0.begin();
d->s1.thermodynamic_forces = s1.begin();
}
const auto r = mfront::gb::integrate<Behaviour>(*d,smf, Ogden_getOutOfBoundsPolicy());
if((r) && (sm != StressMeasure::CAUCHY)){
tfel::math::tensor<2,real> F1;
tfel::fsalgo::copy<5>::exe(d->s1.gradients,F1.begin());
d->s0.thermodynamic_forces = thermodynamic_forces0_old;
d->s1.thermodynamic_forces = thermodynamic_forces1_old;
if(sm==StressMeasure::PK1){
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
}
if((!r) && (sm!=StressMeasure::CAUCHY)){
d->s0.thermodynamic_forces = thermodynamic_forces0_old;
d->s1.thermodynamic_forces = thermodynamic_forces1_old;
}
return r;
} // end of Ogden_PlaneStrain

MFRONT_SHAREDOBJ int Ogden_GeneralisedPlaneStrain(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using TangentOperator = FiniteStrainBehaviourTangentOperatorBase;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::GENERALISEDPLANESTRAIN;
using Behaviour = Ogden<h,real,false>;
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
tfel::math::stensor<2,real> s0;
tfel::math::stensor<2,real> s1;
auto *const thermodynamic_forces0_old = d->s0.thermodynamic_forces;
auto *const thermodynamic_forces1_old = d->s1.thermodynamic_forces;
if(sm!=StressMeasure::CAUCHY){
tfel::math::tensor<2,real> F0;
tfel::fsalgo::copy<5>::exe(d->s0.gradients,F0.begin());
if(sm==StressMeasure::PK1){
tfel::math::TensorView<2,real> pk0(d->s0.thermodynamic_forces);
s0 = tfel::math::convertFirstPiolaKirchhoffStressToCauchyStress(pk0,F0);
} else if(sm==StressMeasure::PK2){
tfel::math::StensorView<2,real> S0(d->s0.thermodynamic_forces);
s0 = tfel::math::convertSecondPiolaKirchhoffStressToCauchyStress(S0,F0);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
d->s0.thermodynamic_forces = s0.begin();
d->s1.thermodynamic_forces = s1.begin();
}
const auto r = mfront::gb::integrate<Behaviour>(*d,smf, Ogden_getOutOfBoundsPolicy());
if((r) && (sm != StressMeasure::CAUCHY)){
tfel::math::tensor<2,real> F1;
tfel::fsalgo::copy<5>::exe(d->s1.gradients,F1.begin());
d->s0.thermodynamic_forces = thermodynamic_forces0_old;
d->s1.thermodynamic_forces = thermodynamic_forces1_old;
if(sm==StressMeasure::PK1){
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
}
if((!r) && (sm!=StressMeasure::CAUCHY)){
d->s0.thermodynamic_forces = thermodynamic_forces0_old;
d->s1.thermodynamic_forces = thermodynamic_forces1_old;
}
return r;
} // end of Ogden_GeneralisedPlaneStrain

MFRONT_SHAREDOBJ int Ogden_Tridimensional(mfront_gb_BehaviourData* const d){
using namespace tfel::material;
using TangentOperator = FiniteStrainBehaviourTangentOperatorBase;
using real = mfront::gb::real;
constexpr auto h = ModellingHypothesis::TRIDIMENSIONAL;
using Behaviour = Ogden<h,real,false>;
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
tfel::math::stensor<3,real> s0;
tfel::math::stensor<3,real> s1;
auto *const thermodynamic_forces0_old = d->s0.thermodynamic_forces;
auto *const thermodynamic_forces1_old = d->s1.thermodynamic_forces;
if(sm!=StressMeasure::CAUCHY){
tfel::math::tensor<3,real> F0;
tfel::fsalgo::copy<9>::exe(d->s0.gradients,F0.begin());
if(sm==StressMeasure::PK1){
tfel::math::TensorView<3,real> pk0(d->s0.thermodynamic_forces);
s0 = tfel::math::convertFirstPiolaKirchhoffStressToCauchyStress(pk0,F0);
} else if(sm==StressMeasure::PK2){
tfel::math::StensorView<3,real> S0(d->s0.thermodynamic_forces);
s0 = tfel::math::convertSecondPiolaKirchhoffStressToCauchyStress(S0,F0);
} else {
  std::cerr << "invalid choice for the "
               "stress measure";
  std::exit(-1);
}
d->s0.thermodynamic_forces = s0.begin();
d->s1.thermodynamic_forces = s1.begin();
}
const auto r = mfront::gb::integrate<Behaviour>(*d,smf, Ogden_getOutOfBoundsPolicy());
if((r) && (sm != StressMeasure::CAUCHY)){
tfel::math::tensor<3,real> F1;
tfel::fsalgo::copy<9>::exe(d->s1.gradients,F1.begin());
d->s0.thermodynamic_forces = thermodynamic_forces0_old;
d->s1.thermodynamic_forces = thermodynamic_forces1_old;
if(sm==StressMeasure::PK1){
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
}
if((!r) && (sm!=StressMeasure::CAUCHY)){
d->s0.thermodynamic_forces = thermodynamic_forces0_old;
d->s1.thermodynamic_forces = thermodynamic_forces1_old;
}
return r;
} // end of Ogden_Tridimensional

#ifdef __cplusplus
}
#endif /* __cplusplus */

