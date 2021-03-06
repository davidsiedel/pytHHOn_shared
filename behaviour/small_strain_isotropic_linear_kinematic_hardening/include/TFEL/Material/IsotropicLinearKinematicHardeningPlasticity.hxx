/*!
* \file   TFEL/Material/IsotropicLinearKinematicHardeningPlasticity.hxx
* \brief  this file implements the IsotropicLinearKinematicHardeningPlasticity Behaviour.
*         File generated by tfel version 4.0.0-dev
* \author Thomas Helfer
* \date   14 / 10 / 2016
 */

#ifndef LIB_TFELMATERIAL_ISOTROPICLINEARKINEMATICHARDENINGPLASTICITY_HXX
#define LIB_TFELMATERIAL_ISOTROPICLINEARKINEMATICHARDENINGPLASTICITY_HXX

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
#include"TFEL/Material/IsotropicLinearKinematicHardeningPlasticityBehaviourData.hxx"
#include"TFEL/Material/IsotropicLinearKinematicHardeningPlasticityIntegrationData.hxx"

#include "MFront/GenericBehaviour/State.hxx"
#include "MFront/GenericBehaviour/BehaviourData.hxx"
namespace tfel{

namespace material{

struct IsotropicLinearKinematicHardeningPlasticityParametersInitializer
{
static IsotropicLinearKinematicHardeningPlasticityParametersInitializer&
get();

double young;
double nu;
double s0;
double C;
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

IsotropicLinearKinematicHardeningPlasticityParametersInitializer();

IsotropicLinearKinematicHardeningPlasticityParametersInitializer(const IsotropicLinearKinematicHardeningPlasticityParametersInitializer&);

IsotropicLinearKinematicHardeningPlasticityParametersInitializer&
operator=(const IsotropicLinearKinematicHardeningPlasticityParametersInitializer&);
/*!
 * \brief read the parameters from the given file
 * \param[out] pi : parameters initializer
 * \param[in]  fn : file name
 */
static void readParameters(IsotropicLinearKinematicHardeningPlasticityParametersInitializer&,const char* const);
};

//! \brief forward declaration
template<ModellingHypothesis::Hypothesis,typename Type,bool use_qt>
class IsotropicLinearKinematicHardeningPlasticity;

//! \brief forward declaration
template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
std::ostream&
 operator <<(std::ostream&,const IsotropicLinearKinematicHardeningPlasticity<hypothesis,Type,false>&);

/*!
* \class IsotropicLinearKinematicHardeningPlasticity
* \brief This class implements the IsotropicLinearKinematicHardeningPlasticity behaviour.
* \param hypothesis, modelling hypothesis.
* \param Type, numerical type.
* \author Thomas Helfer
* \date   14 / 10 / 2016
* An explicit implementation of a simple 
* isotropic plasticity behaviour . 
*/
template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
class IsotropicLinearKinematicHardeningPlasticity<hypothesis,Type,false> final
: public MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>,
public IsotropicLinearKinematicHardeningPlasticityBehaviourData<hypothesis,Type,false>,
public IsotropicLinearKinematicHardeningPlasticityIntegrationData<hypothesis,Type,false>
{

static constexpr unsigned short N = ModellingHypothesisToSpaceDimension<hypothesis>::value;

static_assert(N==1||N==2||N==3);
static_assert(tfel::typetraits::IsFundamentalNumericType<Type>::cond);
static_assert(tfel::typetraits::IsReal<Type>::cond);

friend std::ostream& operator<< <>(std::ostream&,const IsotropicLinearKinematicHardeningPlasticity&);

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

typedef IsotropicLinearKinematicHardeningPlasticityBehaviourData<hypothesis,Type,false> BehaviourData;
typedef IsotropicLinearKinematicHardeningPlasticityIntegrationData<hypothesis,Type,false> IntegrationData;
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



#line 10 "IsotropicLinearKinematicHardeningPlasticity.mfront"
StrainStensor deel;
#line 12 "IsotropicLinearKinematicHardeningPlasticity.mfront"
StrainStensor da;

#line 23 "IsotropicLinearKinematicHardeningPlasticity.mfront"
StressStensor sig0;
#line 24 "IsotropicLinearKinematicHardeningPlasticity.mfront"
Stensor n;
#line 25 "IsotropicLinearKinematicHardeningPlasticity.mfront"
stress lambda;
#line 26 "IsotropicLinearKinematicHardeningPlasticity.mfront"
stress mu;
#line 27 "IsotropicLinearKinematicHardeningPlasticity.mfront"
stress seq;
#line 28 "IsotropicLinearKinematicHardeningPlasticity.mfront"
strain dp;
#line 29 "IsotropicLinearKinematicHardeningPlasticity.mfront"
bool b;

#line 14 "IsotropicLinearKinematicHardeningPlasticity.mfront"
real young;
#line 16 "IsotropicLinearKinematicHardeningPlasticity.mfront"
real nu;
#line 18 "IsotropicLinearKinematicHardeningPlasticity.mfront"
real s0;
#line 20 "IsotropicLinearKinematicHardeningPlasticity.mfront"
real C;
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
this->a += this->da;
}

/*!
* \brief Update auxiliary state variables at end of integration
*/
void updateAuxiliaryStateVariables()
{}

//! \brief Default constructor (disabled)
IsotropicLinearKinematicHardeningPlasticity() =delete ;
//! \brief Copy constructor (disabled)
IsotropicLinearKinematicHardeningPlasticity(const IsotropicLinearKinematicHardeningPlasticity&) = delete;
//! \brief Assignement operator (disabled)
IsotropicLinearKinematicHardeningPlasticity& operator = (const IsotropicLinearKinematicHardeningPlasticity&) = delete;

public:

/*!
* \brief Constructor
*/
IsotropicLinearKinematicHardeningPlasticity(const IsotropicLinearKinematicHardeningPlasticityBehaviourData<hypothesis,Type,false>& src1,
const IsotropicLinearKinematicHardeningPlasticityIntegrationData<hypothesis,Type,false>& src2)
: IsotropicLinearKinematicHardeningPlasticityBehaviourData<hypothesis,Type,false>(src1),
IsotropicLinearKinematicHardeningPlasticityIntegrationData<hypothesis,Type,false>(src2),
deel(typename tfel::math::MathObjectTraits<StrainStensor>::NumType(0)),
da(typename tfel::math::MathObjectTraits<StrainStensor>::NumType(0)),
dsig_ddeto(Dt)
{
using namespace std;
using namespace tfel::math;
using std::vector;
this->young = IsotropicLinearKinematicHardeningPlasticityParametersInitializer::get().young;
this->nu = IsotropicLinearKinematicHardeningPlasticityParametersInitializer::get().nu;
this->s0 = IsotropicLinearKinematicHardeningPlasticityParametersInitializer::get().s0;
this->C = IsotropicLinearKinematicHardeningPlasticityParametersInitializer::get().C;
this->minimal_time_step_scaling_factor = IsotropicLinearKinematicHardeningPlasticityParametersInitializer::get().minimal_time_step_scaling_factor;
this->maximal_time_step_scaling_factor = IsotropicLinearKinematicHardeningPlasticityParametersInitializer::get().maximal_time_step_scaling_factor;
}

/*
 * \brief constructor for the Generic interface
 * \param[in] mgb_d: behaviour data
 */
IsotropicLinearKinematicHardeningPlasticity(const mfront::gb::BehaviourData& mgb_d)
: IsotropicLinearKinematicHardeningPlasticityBehaviourData<hypothesis,Type,false>(mgb_d),
IsotropicLinearKinematicHardeningPlasticityIntegrationData<hypothesis,Type,false>(mgb_d),
deel(typename tfel::math::MathObjectTraits<StrainStensor>::NumType(0)),
da(typename tfel::math::MathObjectTraits<StrainStensor>::NumType(0)),
dsig_ddeto(Dt)
{
using namespace std;
using namespace tfel::math;
using std::vector;
this->young = IsotropicLinearKinematicHardeningPlasticityParametersInitializer::get().young;
this->nu = IsotropicLinearKinematicHardeningPlasticityParametersInitializer::get().nu;
this->s0 = IsotropicLinearKinematicHardeningPlasticityParametersInitializer::get().s0;
this->C = IsotropicLinearKinematicHardeningPlasticityParametersInitializer::get().C;
this->minimal_time_step_scaling_factor = IsotropicLinearKinematicHardeningPlasticityParametersInitializer::get().minimal_time_step_scaling_factor;
this->maximal_time_step_scaling_factor = IsotropicLinearKinematicHardeningPlasticityParametersInitializer::get().maximal_time_step_scaling_factor;
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
#line 32 "IsotropicLinearKinematicHardeningPlasticity.mfront"
this->lambda = computeLambda(this->young,this->nu);
#line 33 "IsotropicLinearKinematicHardeningPlasticity.mfront"
this->mu     = computeMu(this->young,this->nu);
#line 34 "IsotropicLinearKinematicHardeningPlasticity.mfront"
this->sig0   = this->sig;
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
#line 38 "IsotropicLinearKinematicHardeningPlasticity.mfront"
(this->Dt) = (this->lambda)*Stensor4::IxI()+2*(this->mu)*Stensor4::Id();return SUCCESS;
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
#line 42 "IsotropicLinearKinematicHardeningPlasticity.mfront"
this->eel += this->deto;
#line 43 "IsotropicLinearKinematicHardeningPlasticity.mfront"
const auto X  = 2*this->C*this->a/3;
#line 44 "IsotropicLinearKinematicHardeningPlasticity.mfront"
const auto s  = 2*this->mu*deviator(this->eel)-X;
#line 45 "IsotropicLinearKinematicHardeningPlasticity.mfront"
this->seq = sigmaeq(s);
#line 46 "IsotropicLinearKinematicHardeningPlasticity.mfront"
this->b   = this->seq-this->s0>stress{0};
#line 47 "IsotropicLinearKinematicHardeningPlasticity.mfront"
if(this->b){
#line 48 "IsotropicLinearKinematicHardeningPlasticity.mfront"
this->n    = (3*s)/(2*this->seq);
#line 49 "IsotropicLinearKinematicHardeningPlasticity.mfront"
this->dp   = (this->seq-this->s0)/(3*this->mu+this->C);
#line 50 "IsotropicLinearKinematicHardeningPlasticity.mfront"
this->a   += this->dp*this->n;
#line 51 "IsotropicLinearKinematicHardeningPlasticity.mfront"
this->eel -= this->dp*this->n;
#line 52 "IsotropicLinearKinematicHardeningPlasticity.mfront"
}
#line 53 "IsotropicLinearKinematicHardeningPlasticity.mfront"
this->sig = this->lambda*trace(this->eel)*Stensor::Id()+2*this->mu*this->eel;
this->updateIntegrationVariables();
this->updateStateVariables();
this->updateAuxiliaryStateVariables();
if(computeTangentOperator_){
if(!this->computeConsistentTangentOperator(smt)){
return MechanicalBehaviour<MechanicalBehaviourBase::STANDARDSTRAINBASEDBEHAVIOUR,hypothesis,Type,false>::FAILURE;
}
}
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
void computeInternalEnergy(real& Psi_s) const{
using namespace std;
using namespace tfel::math;
#line 66 "IsotropicLinearKinematicHardeningPlasticity.mfront"
const auto tr = trace(this->eel);
#line 67 "IsotropicLinearKinematicHardeningPlasticity.mfront"
Psi_s = this->lambda/2*tr*tr+this->mu*(this->eel|this->eel);
}

/*!
* \brief Update the dissipated energy at end of the time step
* \param[in] Psi_d: dissipated energy at end of the time step
*/
void computeDissipatedEnergy(real& Psi_d) const{
using namespace std;
using namespace tfel::math;
#line 71 "IsotropicLinearKinematicHardeningPlasticity.mfront"
Psi_d += ((this->sig+this->sig0)|this->n)*this->dp/2;
}

bool computeConsistentTangentOperator(const SMType smt){
using namespace std;
using namespace tfel::math;
using std::vector;
#line 57 "IsotropicLinearKinematicHardeningPlasticity.mfront"
if((smt==CONSISTENTTANGENTOPERATOR)&&(this->b)){
#line 58 "IsotropicLinearKinematicHardeningPlasticity.mfront"
this->Dt = (this->lambda*Stensor4::IxI()+2*this->mu*Stensor4::Id()
#line 59 "IsotropicLinearKinematicHardeningPlasticity.mfront"
-4*this->mu*this->mu*((this->dp/this->seq)*(Stensor4::M()-(this->n^this->n))+(this->n^this->n)/(3*this->mu+this->C)));
#line 60 "IsotropicLinearKinematicHardeningPlasticity.mfront"
} else {
#line 61 "IsotropicLinearKinematicHardeningPlasticity.mfront"
this->Dt = this->lambda*Stensor4::IxI()+2*this->mu*Stensor4::Id();
#line 62 "IsotropicLinearKinematicHardeningPlasticity.mfront"
}
return true;
}

const TangentOperator& getTangentOperator() const{
return this->Dt;
}

void updateExternalStateVariables(){
this->eto  += this->deto;
this->T += this->dT;
}

//!
~IsotropicLinearKinematicHardeningPlasticity()
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
}; // end of IsotropicLinearKinematicHardeningPlasticity class

template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
std::ostream&
operator <<(std::ostream& os,const IsotropicLinearKinematicHardeningPlasticity<hypothesis,Type,false>& b)
{
os << "εᵗᵒ : " << b.eto << '\n';
os << "Δεᵗᵒ : " << b.deto << '\n';
os << "σ : " << b.sig << '\n';
os << "Δt : " << b.dt << '\n';
os << "eel : " << b.eel << '\n';
os << "Δeel : " << b.deel << '\n';
os << "a : " << b.a << '\n';
os << "Δa : " << b.da << '\n';
os << "T : " << b.T << '\n';
os << "ΔT : " << b.dT << '\n';
os << "sig0 : " << b.sig0 << '\n';
os << "n : " << b.n << '\n';
os << "young : " << b.young << '\n';
os << "nu : " << b.nu << '\n';
os << "s0 : " << b.s0 << '\n';
os << "C : " << b.C << '\n';
os << "minimal_time_step_scaling_factor : " << b.minimal_time_step_scaling_factor << '\n';
os << "maximal_time_step_scaling_factor : " << b.maximal_time_step_scaling_factor << '\n';
return os;
}

/*!
* Partial specialisation for IsotropicLinearKinematicHardeningPlasticity.
*/
template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
class MechanicalBehaviourTraits<IsotropicLinearKinematicHardeningPlasticity<hypothesis,Type,false> >
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
static constexpr unsigned short internal_variables_nb  = 2*StensorSize;
static constexpr unsigned short external_variables_nb  = 1;
static constexpr unsigned short external_variables_nb2 = 0;
static constexpr bool hasConsistentTangentOperator = true;
static constexpr bool isConsistentTangentOperatorSymmetric = false;
static constexpr bool hasPredictionOperator = true;
static constexpr bool hasAPrioriTimeStepScalingFactor = false;
static constexpr bool hasComputeInternalEnergy = true;
static constexpr bool hasComputeDissipatedEnergy = true;
/*!
* \return the name of the class.
*/
static const char* getName(){
return "IsotropicLinearKinematicHardeningPlasticity";
}

};

/*!
* Partial specialisation for IsotropicLinearKinematicHardeningPlasticity.
*/
template<typename Type>
class MechanicalBehaviourTraits<IsotropicLinearKinematicHardeningPlasticity<ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS,Type,false> >
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
return "IsotropicLinearKinematicHardeningPlasticity";
}

};

/*!
* Partial specialisation for IsotropicLinearKinematicHardeningPlasticity.
*/
template<typename Type>
class MechanicalBehaviourTraits<IsotropicLinearKinematicHardeningPlasticity<ModellingHypothesis::PLANESTRESS,Type,false> >
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
return "IsotropicLinearKinematicHardeningPlasticity";
}

};

} // end of namespace material

} // end of namespace tfel

#endif /* LIB_TFELMATERIAL_ISOTROPICLINEARKINEMATICHARDENINGPLASTICITY_HXX */
