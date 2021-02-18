/*!
* \file   TFEL/Material/IsotropicLinearHardeningPlasticity.hxx
* \brief  this file implements the IsotropicLinearHardeningPlasticity Behaviour.
*         File generated by tfel version 4.0.0-dev
* \author Thomas Helfer
* \date   14 / 10 / 2016
 */

#ifndef LIB_TFELMATERIAL_ISOTROPICLINEARHARDENINGPLASTICITY_HXX
#define LIB_TFELMATERIAL_ISOTROPICLINEARHARDENINGPLASTICITY_HXX

#include<string>
#include<iostream>
#include<limits>
#include<stdexcept>
#include<algorithm>

#include"TFEL/Raise.hxx"
#include"TFEL/PhysicalConstants.hxx"
#include"TFEL/Config/TFELConfig.hxx"
#include"TFEL/Config/TFELTypes.hxx"
#include"TFEL/TypeTraits/IsFundamentalNumericType.hxx"
#include"TFEL/TypeTraits/IsReal.hxx"
#include"TFEL/Math/General/IEEE754.hxx"
#include"TFEL/Material/MaterialException.hxx"
#include"TFEL/Material/MechanicalBehaviour.hxx"
#include"TFEL/Material/MechanicalBehaviourTraits.hxx"
#include"TFEL/Material/OutOfBoundsPolicy.hxx"
#include"TFEL/Material/BoundsCheck.hxx"
#include"TFEL/Material/IsotropicPlasticity.hxx"
#include"TFEL/Material/Lame.hxx"
#include"TFEL/Material/Hosford1972YieldCriterion.hxx"
#include"TFEL/Material/LogarithmicStrainComputeAxialStrainIncrementElasticPrediction.hxx"
#include"TFEL/Material/IsotropicLinearHardeningPlasticityBehaviourData.hxx"
#include"TFEL/Material/IsotropicLinearHardeningPlasticityIntegrationData.hxx"

#include "MFront/GenericBehaviour/State.hxx"
#include "MFront/GenericBehaviour/BehaviourData.hxx"
namespace tfel{

namespace material{

struct IsotropicLinearHardeningPlasticityParametersInitializer
{
static IsotropicLinearHardeningPlasticityParametersInitializer&
get();

double young;
double nu;
double H;
double s0;
double minimal_time_step_scaling_factor;
double maximal_time_step_scaling_factor;

void set(const char* const,const double);

/*!
 * \brief convert a string to double
 * \param[in] p : parameter
 * \param[in] v : value
 */
static double getDouble(const std::string&,const std::string&);
private :

IsotropicLinearHardeningPlasticityParametersInitializer();

IsotropicLinearHardeningPlasticityParametersInitializer(const IsotropicLinearHardeningPlasticityParametersInitializer&);

IsotropicLinearHardeningPlasticityParametersInitializer&
operator=(const IsotropicLinearHardeningPlasticityParametersInitializer&);
/*!
 * \brief read the parameters from the given file
 * \param[out] pi : parameters initializer
 * \param[in]  fn : file name
 */
static void readParameters(IsotropicLinearHardeningPlasticityParametersInitializer&,const char* const);
};

//! \brief forward declaration
template<ModellingHypothesis::Hypothesis,typename Type,bool use_qt>
class IsotropicLinearHardeningPlasticity;

//! \brief forward declaration
template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
std::ostream&
 operator <<(std::ostream&,const IsotropicLinearHardeningPlasticity<hypothesis,Type,false>&);

/*!
* \class IsotropicLinearHardeningPlasticity
* \brief This class implements the IsotropicLinearHardeningPlasticity behaviour.
* \param hypothesis, modelling hypothesis.
* \param Type, numerical type.
* \author Thomas Helfer
* \date   14 / 10 / 2016
* An implicit implementation of a simple 
* isotropic plasticity behaviour with 
* isotropic linear hardening . 
* 
* The yield surface is defined by : 
* \[ 
*   f(\sigmaeq,p) = \sigmaeq-s_{0}-H\,p 
* \] 
*/
template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
class IsotropicLinearHardeningPlasticity<hypothesis,Type,false> final
: public MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>,
public IsotropicLinearHardeningPlasticityBehaviourData<hypothesis,Type,false>,
public IsotropicLinearHardeningPlasticityIntegrationData<hypothesis,Type,false>
{

static constexpr unsigned short N = ModellingHypothesisToSpaceDimension<hypothesis>::value;

static_assert(N==1||N==2||N==3);
static_assert(tfel::typetraits::IsFundamentalNumericType<Type>::cond);
static_assert(tfel::typetraits::IsReal<Type>::cond);

friend std::ostream& operator<< <>(std::ostream&,const IsotropicLinearHardeningPlasticity&);

static constexpr unsigned short TVectorSize = N;
typedef tfel::math::StensorDimeToSize<N> StensorDimeToSize;
static constexpr unsigned short StensorSize = StensorDimeToSize::value;
typedef tfel::math::TensorDimeToSize<N> TensorDimeToSize;
static constexpr unsigned short TensorSize = TensorDimeToSize::value;

using ushort =  unsigned short;
using Types = tfel::config::Types<N,Type,false>;
using real                = typename Types::real;
using time                = typename Types::time;
using length              = typename Types::length;
using frequency           = typename Types::frequency;
using stress              = typename Types::stress;
using strain              = typename Types::strain;
using strainrate          = typename Types::strainrate;
using stressrate          = typename Types::stressrate;
using temperature         = typename Types::temperature;
using thermalexpansion    = typename Types::thermalexpansion;
using thermalconductivity = typename Types::thermalconductivity;
using massdensity         = typename Types::massdensity;
using energydensity         = typename Types::energydensity;
using TVector             = typename Types::TVector;
using Stensor             = typename Types::Stensor;
using Stensor4            = typename Types::Stensor4;
using FrequencyStensor    = typename Types::FrequencyStensor;
using ForceTVector        = typename Types::ForceTVector;
using StressStensor       = typename Types::StressStensor;
using StressRateStensor   = typename Types::StressRateStensor;
using DisplacementTVector = typename Types::DisplacementTVector;
using StrainStensor       = typename Types::StrainStensor;
using StrainRateStensor   = typename Types::StrainRateStensor;
using StiffnessTensor     = typename Types::StiffnessTensor;
using Tensor              = typename Types::Tensor;
using FrequencyTensor     = typename Types::FrequencyTensor;
using StressTensor        = typename Types::StressTensor;
using ThermalExpansionCoefficientTensor = typename Types::ThermalExpansionCoefficientTensor;
using DeformationGradientTensor         = typename Types::DeformationGradientTensor;
using DeformationGradientRateTensor     = typename Types::DeformationGradientRateTensor;
using TemperatureGradient = typename Types::TemperatureGradient;
using HeatFlux = typename Types::HeatFlux;
using TangentOperator   = StiffnessTensor;
using PhysicalConstants = tfel::PhysicalConstants<real>;

public :

typedef IsotropicLinearHardeningPlasticityBehaviourData<hypothesis,Type,false> BehaviourData;
typedef IsotropicLinearHardeningPlasticityIntegrationData<hypothesis,Type,false> IntegrationData;
typedef typename MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>::SMFlag SMFlag;
typedef typename MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>::SMType SMType;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>::ELASTIC;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>::SECANTOPERATOR;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>::TANGENTOPERATOR;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>::CONSISTENTTANGENTOPERATOR;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>::NOSTIFFNESSREQUESTED;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>::STANDARDTANGENTOPERATOR;
using IntegrationResult = typename MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>::IntegrationResult;

using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>::SUCCESS;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>::FAILURE;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>::UNRELIABLE_RESULTS;

using StressFreeExpansionType = StrainStensor;

private :



#line 18 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
StrainStensor deel;
#line 20 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
strain dp;


#line 23 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
real young;
#line 25 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
real nu;
#line 27 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
real H;
#line 29 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
real s0;
real minimal_time_step_scaling_factor;
real maximal_time_step_scaling_factor;

//! Tangent operator;
TangentOperator Dt;
//! alias to the tangent operator;
TangentOperator& dsig_ddeto;
/*!
* \brief Update internal variables at end of integration
*/
void updateIntegrationVariables(){
}

/*!
* \brief Update internal variables at end of integration
*/
void updateStateVariables(){
this->eel += this->deel;
this->p += this->dp;
}

/*!
* \brief Update auxiliary state variables at end of integration
*/
void updateAuxiliaryStateVariables()
{}

//! \brief Default constructor (disabled)
IsotropicLinearHardeningPlasticity() =delete ;
//! \brief Copy constructor (disabled)
IsotropicLinearHardeningPlasticity(const IsotropicLinearHardeningPlasticity&) = delete;
//! \brief Assignement operator (disabled)
IsotropicLinearHardeningPlasticity& operator = (const IsotropicLinearHardeningPlasticity&) = delete;

public:

/*!
* \brief Constructor
*/
IsotropicLinearHardeningPlasticity(const IsotropicLinearHardeningPlasticityBehaviourData<hypothesis,Type,false>& src1,
const IsotropicLinearHardeningPlasticityIntegrationData<hypothesis,Type,false>& src2)
: IsotropicLinearHardeningPlasticityBehaviourData<hypothesis,Type,false>(src1),
IsotropicLinearHardeningPlasticityIntegrationData<hypothesis,Type,false>(src2),
deel(typename tfel::math::MathObjectTraits<StrainStensor>::NumType(0)),
dp(strain(0)),
dsig_ddeto(Dt)
{
using namespace std;
using namespace tfel::math;
using std::vector;
this->young = IsotropicLinearHardeningPlasticityParametersInitializer::get().young;
this->nu = IsotropicLinearHardeningPlasticityParametersInitializer::get().nu;
this->H = IsotropicLinearHardeningPlasticityParametersInitializer::get().H;
this->s0 = IsotropicLinearHardeningPlasticityParametersInitializer::get().s0;
this->minimal_time_step_scaling_factor = IsotropicLinearHardeningPlasticityParametersInitializer::get().minimal_time_step_scaling_factor;
this->maximal_time_step_scaling_factor = IsotropicLinearHardeningPlasticityParametersInitializer::get().maximal_time_step_scaling_factor;
}

/*
 * \brief constructor for the Generic interface
 * \param[in] mgb_d: behaviour data
 */
IsotropicLinearHardeningPlasticity(const mfront::gb::BehaviourData& mgb_d)
: IsotropicLinearHardeningPlasticityBehaviourData<hypothesis,Type,false>(mgb_d),
IsotropicLinearHardeningPlasticityIntegrationData<hypothesis,Type,false>(mgb_d),
deel(typename tfel::math::MathObjectTraits<StrainStensor>::NumType(0)),
dp(strain(0)),
dsig_ddeto(Dt)
{
using namespace std;
using namespace tfel::math;
using std::vector;
this->young = IsotropicLinearHardeningPlasticityParametersInitializer::get().young;
this->nu = IsotropicLinearHardeningPlasticityParametersInitializer::get().nu;
this->H = IsotropicLinearHardeningPlasticityParametersInitializer::get().H;
this->s0 = IsotropicLinearHardeningPlasticityParametersInitializer::get().s0;
this->minimal_time_step_scaling_factor = IsotropicLinearHardeningPlasticityParametersInitializer::get().minimal_time_step_scaling_factor;
this->maximal_time_step_scaling_factor = IsotropicLinearHardeningPlasticityParametersInitializer::get().maximal_time_step_scaling_factor;
tfel::fsalgo::copy<StensorSize >::exe(mgb_d.s0.gradients,this->eto.begin());
tfel::fsalgo::transform<StensorSize>::exe(mgb_d.s1.gradients,mgb_d.s0.gradients,this->deto.begin(),std::minus<real>());
tfel::fsalgo::copy<StensorSize >::exe(mgb_d.s0.thermodynamic_forces,this->sig.begin());
}

/*!
 * \ brief initialize the behaviour with user code
 */
void initialize(){
using namespace std;
using namespace tfel::math;
using std::vector;
}

/*!
* \brief set the policy for "out of bounds" conditions
*/
void
setOutOfBoundsPolicy(const OutOfBoundsPolicy policy_value){
this->policy = policy_value;
} // end of setOutOfBoundsPolicy

/*!
* \return the modelling hypothesis
*/
constexpr ModellingHypothesis::Hypothesis
getModellingHypothesis() const{
return hypothesis;
} // end of getModellingHypothesis

/*!
* \brief check bounds
*/
void checkBounds() const{
} // end of checkBounds

IntegrationResult
computePredictionOperator(const SMFlag smflag,const SMType smt) override{
using namespace std;
using namespace tfel::math;
using std::vector;
tfel::raise_if(smflag!=MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>::STANDARDTANGENTOPERATOR,
"invalid prediction operator flag");
#line 41 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
const auto lambda = computeLambda((this->young),(this->nu));
#line 42 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
const auto mu     = computeMu((this->young),(this->nu));
#line 43 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
(this->Dt) = lambda*Stensor4::IxI()+2*mu*Stensor4::Id();return SUCCESS;
}

real getMinimalTimeStepScalingFactor() const override{
  return this->minimal_time_step_scaling_factor;
}

std::pair<bool,real>
computeAPrioriTimeStepScalingFactor(const real current_time_step_scaling_factor) const override{
const auto time_scaling_factor = this->computeAPrioriTimeStepScalingFactorII();
return {time_scaling_factor.first,
        std::min(std::min(std::max(time_scaling_factor.second,
                                   this->minimal_time_step_scaling_factor),
                          this->maximal_time_step_scaling_factor),
                  current_time_step_scaling_factor)};
}

/*!
* \brief Integrate behaviour  over the time step
*/
IntegrationResult
integrate(const SMFlag smflag, const SMType smt) override{
using namespace std;
using namespace tfel::math;
raise_if(smflag!=MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>::STANDARDTANGENTOPERATOR,
"invalid tangent operator flag");
bool computeTangentOperator_ = smt!=NOSTIFFNESSREQUESTED;
#line 51 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
const auto lambda = computeLambda(this->young,this->nu);
#line 52 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
const auto mu     = computeMu(this->young,this->nu);
#line 53 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
this->eel += this->deto;
#line 54 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
const auto se     = 2*mu*deviator(this->eel);
#line 55 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
const auto seq_e  = sigmaeq(se);
#line 56 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
const auto b      = seq_e-this->s0-this->H*this->p>stress{0};
#line 57 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
if(b){
#line 58 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
const auto iseq_e = 1/seq_e;
#line 59 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
const auto n      = eval(3*se/(2*seq_e));
#line 60 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
const auto cste   = 1/(this->H+3*mu);
#line 61 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
this->dp   = (seq_e-this->s0-this->H*this->p)*cste;
#line 62 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
this->eel -= this->dp*n;
#line 63 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
if (computeTangentOperator_) {
#line 64 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
if (smt == CONSISTENTTANGENTOPERATOR) {
#line 65 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
this->Dt = (lambda * Stensor4::IxI() + 2 * mu * Stensor4::Id() -
#line 66 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
4 * mu * mu *
#line 67 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
(this->dp * iseq_e * (Stensor4::M() - (n ^ n)) + cste * (n ^ n)));
#line 68 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
} else {
#line 69 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
this->Dt = lambda * Stensor4::IxI() + 2 * mu * Stensor4::Id();
#line 70 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
}
#line 71 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
}
#line 72 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
} else {
#line 73 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
if (computeTangentOperator_) {
#line 74 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
this->Dt = lambda * Stensor4::IxI() + 2 * mu * Stensor4::Id();
#line 75 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
}
#line 76 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
}
#line 77 "FiniteStrainIsotropicLinearHardeningPlasticity.mfront"
this->sig = lambda * trace(this->eel) * Stensor::Id() + 2 * mu * this->eel;
this->updateIntegrationVariables();
this->updateStateVariables();
this->updateAuxiliaryStateVariables();
return MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>::SUCCESS;
}

std::pair<bool,real>
computeAPosterioriTimeStepScalingFactor(const real current_time_step_scaling_factor) const override{
const auto time_scaling_factor = this->computeAPosterioriTimeStepScalingFactorII();
return {time_scaling_factor.first,
        std::min(std::min(std::max(time_scaling_factor.second,
                                   this->minimal_time_step_scaling_factor),
                          this->maximal_time_step_scaling_factor),
                 current_time_step_scaling_factor)};
}

/*!
* \brief Update the internal energy at end of the time step
* \param[in] Psi_s: internal energy at end of the time step
*/
void computeInternalEnergy(real& Psi_s) const
{
Psi_s=0;
}

/*!
* \brief Update the dissipated energy at end of the time step
* \param[in] Psi_d: dissipated energy at end of the time step
*/
void computeDissipatedEnergy(real& Psi_d) const
{
Psi_d=0;
}

const TangentOperator& getTangentOperator() const{
return this->Dt;
}

void updateExternalStateVariables(){
this->eto  += this->deto;
this->T += this->dT;
}

//!
~IsotropicLinearHardeningPlasticity()
 override = default;

private:

std::pair<bool,real> computeAPrioriTimeStepScalingFactorII() const{
return {true,this->maximal_time_step_scaling_factor};
}

std::pair<bool,real> computeAPosterioriTimeStepScalingFactorII() const{
return {true,this->maximal_time_step_scaling_factor};
}

//! policy for treating out of bounds conditions
OutOfBoundsPolicy policy = None;
}; // end of IsotropicLinearHardeningPlasticity class

template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
std::ostream&
operator <<(std::ostream& os,const IsotropicLinearHardeningPlasticity<hypothesis,Type,false>& b)
{
os << "εᵗᵒ : " << b.eto << '\n';
os << "Δεᵗᵒ : " << b.deto << '\n';
os << "σ : " << b.sig << '\n';
os << "Δt : " << b.dt << '\n';
os << "eel : " << b.eel << '\n';
os << "Δeel : " << b.deel << '\n';
os << "p : " << b.p << '\n';
os << "Δp : " << b.dp << '\n';
os << "T : " << b.T << '\n';
os << "ΔT : " << b.dT << '\n';
os << "young : " << b.young << '\n';
os << "nu : " << b.nu << '\n';
os << "H : " << b.H << '\n';
os << "s0 : " << b.s0 << '\n';
os << "minimal_time_step_scaling_factor : " << b.minimal_time_step_scaling_factor << '\n';
os << "maximal_time_step_scaling_factor : " << b.maximal_time_step_scaling_factor << '\n';
return os;
}

/*!
* Partial specialisation for IsotropicLinearHardeningPlasticity.
*/
template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
class MechanicalBehaviourTraits<IsotropicLinearHardeningPlasticity<hypothesis,Type,false> >
{
static constexpr unsigned short N = ModellingHypothesisToSpaceDimension<hypothesis>::value;
static constexpr unsigned short TVectorSize = N;
typedef tfel::math::StensorDimeToSize<N> StensorDimeToSize;
static constexpr unsigned short StensorSize = StensorDimeToSize::value;
typedef tfel::math::TensorDimeToSize<N> TensorDimeToSize;
static constexpr unsigned short TensorSize = TensorDimeToSize::value;
public:
static constexpr bool is_defined = true;
static constexpr bool use_quantities = false;
static constexpr bool hasStressFreeExpansion = false;
static constexpr bool handlesThermalExpansion = false;
static constexpr unsigned short dimension = N;
typedef Type NumType;
static constexpr unsigned short material_properties_nb = 0;
static constexpr unsigned short internal_variables_nb  = 1+StensorSize;
static constexpr unsigned short external_variables_nb  = 1;
static constexpr unsigned short external_variables_nb2 = 0;
static constexpr bool hasConsistentTangentOperator = true;
static constexpr bool isConsistentTangentOperatorSymmetric = true;
static constexpr bool hasPredictionOperator = true;
static constexpr bool hasAPrioriTimeStepScalingFactor = false;
static constexpr bool hasComputeInternalEnergy = false;
static constexpr bool hasComputeDissipatedEnergy = false;
/*!
* \return the name of the class.
*/
static const char* getName(){
return "IsotropicLinearHardeningPlasticity";
}

};

/*!
* Partial specialisation for IsotropicLinearHardeningPlasticity.
*/
template<typename Type>
class MechanicalBehaviourTraits<IsotropicLinearHardeningPlasticity<ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS,Type,false> >
{
public:
static constexpr bool is_defined = false;
static constexpr bool use_quantities = false;
static constexpr bool hasStressFreeExpansion = false;
static constexpr bool handlesThermalExpansion = false;
static constexpr unsigned short dimension = 0u;
typedef Type NumType;
static constexpr unsigned short material_properties_nb = 0;
static constexpr unsigned short internal_variables_nb  = 0;
static constexpr unsigned short external_variables_nb  = 0;
static constexpr unsigned short external_variables_nb2 = 0;
static constexpr bool hasConsistentTangentOperator = false;
static constexpr bool isConsistentTangentOperatorSymmetric = false;
static constexpr bool hasPredictionOperator = false;
static constexpr bool hasAPrioriTimeStepScalingFactor = false;
static constexpr bool hasComputeInternalEnergy = false;
static constexpr bool hasComputeDissipatedEnergy = false;
/*!
* \return the name of the class.
*/
static const char* getName(){
return "IsotropicLinearHardeningPlasticity";
}

};

/*!
* Partial specialisation for IsotropicLinearHardeningPlasticity.
*/
template<typename Type>
class MechanicalBehaviourTraits<IsotropicLinearHardeningPlasticity<ModellingHypothesis::PLANESTRESS,Type,false> >
{
public:
static constexpr bool is_defined = false;
static constexpr bool use_quantities = false;
static constexpr bool hasStressFreeExpansion = false;
static constexpr bool handlesThermalExpansion = false;
static constexpr unsigned short dimension = 0u;
typedef Type NumType;
static constexpr unsigned short material_properties_nb = 0;
static constexpr unsigned short internal_variables_nb  = 0;
static constexpr unsigned short external_variables_nb  = 0;
static constexpr unsigned short external_variables_nb2 = 0;
static constexpr bool hasConsistentTangentOperator = false;
static constexpr bool isConsistentTangentOperatorSymmetric = false;
static constexpr bool hasPredictionOperator = false;
static constexpr bool hasAPrioriTimeStepScalingFactor = false;
static constexpr bool hasComputeInternalEnergy = false;
static constexpr bool hasComputeDissipatedEnergy = false;
/*!
* \return the name of the class.
*/
static const char* getName(){
return "IsotropicLinearHardeningPlasticity";
}

};

} // end of namespace material

} // end of namespace tfel

#endif /* LIB_TFELMATERIAL_ISOTROPICLINEARHARDENINGPLASTICITY_HXX */
