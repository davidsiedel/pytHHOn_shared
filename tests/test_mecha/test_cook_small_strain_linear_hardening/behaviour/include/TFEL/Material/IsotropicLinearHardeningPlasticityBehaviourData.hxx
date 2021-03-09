/*!
* \file   TFEL/Material/IsotropicLinearHardeningPlasticityBehaviourData.hxx
* \brief  this file implements the IsotropicLinearHardeningPlasticityBehaviourData class.
*         File generated by tfel version 4.0.0-dev
* \author Thomas Helfer
* \date   14 / 10 / 2016
 */

#ifndef LIB_TFELMATERIAL_ISOTROPICLINEARHARDENINGPLASTICITY_BEHAVIOUR_DATA_HXX
#define LIB_TFELMATERIAL_ISOTROPICLINEARHARDENINGPLASTICITY_BEHAVIOUR_DATA_HXX

#include<limits>
#include<string>
#include<sstream>
#include<iostream>
#include<stdexcept>
#include<algorithm>

#include"TFEL/Raise.hxx"
#include"TFEL/PhysicalConstants.hxx"
#include"TFEL/Config/TFELConfig.hxx"
#include"TFEL/Config/TFELTypes.hxx"
#include"TFEL/TypeTraits/IsFundamentalNumericType.hxx"
#include"TFEL/TypeTraits/IsReal.hxx"
#include"TFEL/Math/General/IEEE754.hxx"
#include"TFEL/Math/stensor.hxx"
#include"TFEL/Math/Stensor/StensorConceptIO.hxx"
#include"TFEL/Math/tmatrix.hxx"
#include"TFEL/Math/Matrix/tmatrixIO.hxx"
#include"TFEL/Math/st2tost2.hxx"
#include"TFEL/Math/ST2toST2/ST2toST2ConceptIO.hxx"
#include"TFEL/Material/ModellingHypothesis.hxx"

#include "MFront/GenericBehaviour/State.hxx"
#include "MFront/GenericBehaviour/BehaviourData.hxx"
namespace tfel{

namespace material{

//! \brief forward declaration
template<ModellingHypothesis::Hypothesis hypothesis,typename,bool>
class IsotropicLinearHardeningPlasticityBehaviourData;

//! \brief forward declaration
template<ModellingHypothesis::Hypothesis hypothesis,typename Type,bool use_qt>
class IsotropicLinearHardeningPlasticityIntegrationData;

//! \brief forward declaration
template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
std::ostream&
 operator <<(std::ostream&,const IsotropicLinearHardeningPlasticityBehaviourData<hypothesis,Type,false>&);

template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
class IsotropicLinearHardeningPlasticityBehaviourData<hypothesis,Type,false>
{

static constexpr unsigned short N = ModellingHypothesisToSpaceDimension<hypothesis>::value;
static_assert(N==1||N==2||N==3);
static_assert(tfel::typetraits::IsFundamentalNumericType<Type>::cond);
static_assert(tfel::typetraits::IsReal<Type>::cond);

friend std::ostream& operator<< <>(std::ostream&,const IsotropicLinearHardeningPlasticityBehaviourData&);

/* integration data is declared friend to access   driving variables at the beginning of the time step */
friend class IsotropicLinearHardeningPlasticityIntegrationData<hypothesis,Type,false>;

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

protected:

StrainStensor eto;

StressStensor sig;


#line 16 "IsotropicLinearHardeningPlasticity.mfront"
StrainStensor eel;
#line 18 "IsotropicLinearHardeningPlasticity.mfront"
strain p;
temperature T;

public:

/*!
* \brief Default constructor
*/
IsotropicLinearHardeningPlasticityBehaviourData()
{}

/*!
* \brief copy constructor
*/
IsotropicLinearHardeningPlasticityBehaviourData(const IsotropicLinearHardeningPlasticityBehaviourData& src)
: eto(src.eto),
sig(src.sig),
eel(src.eel),
p(src.p),
T(src.T)
{}

/*
 * \brief constructor for the Generic interface
 * \param[in] mgb_d: behaviour data
 */
IsotropicLinearHardeningPlasticityBehaviourData(const mfront::gb::BehaviourData& mgb_d)
: eel(&mgb_d.s0.internal_state_variables[0]),
p(mgb_d.s0.internal_state_variables[StensorSize]),
T(mgb_d.s0.external_state_variables[0])
{
}


/*
* \brief Assignement operator
*/
IsotropicLinearHardeningPlasticityBehaviourData&
operator=(const IsotropicLinearHardeningPlasticityBehaviourData& src){
this->eto = src.eto;
this->sig = src.sig;
this->eel = src.eel;
this->p = src.p;
this->T = src.T;
return *this;
}

void exportStateData(mfront::gb::State& mbg_s1) const
{
using namespace tfel::math;
tfel::fsalgo::copy<StensorSize>::exe(this->sig.begin(), mbg_s1.thermodynamic_forces);
tfel::fsalgo::copy<StensorSize>::exe(this->eel.begin(), mbg_s1.internal_state_variables);
mbg_s1.internal_state_variables[StensorSize] = this->p;
} // end of exportStateData

}; // end of IsotropicLinearHardeningPlasticityBehaviourDataclass

template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
std::ostream&
operator <<(std::ostream& os,const IsotropicLinearHardeningPlasticityBehaviourData<hypothesis,Type,false>& b)
{
os << "εᵗᵒ : " << b.eto << '\n';
os << "σ : " << b.sig << '\n';
os << "eel : " << b.eel << '\n';
os << "p : " << b.p << '\n';
os << "T : " << b.T << '\n';
return os;
}

} // end of namespace material

} // end of namespace tfel

#endif /* LIB_TFELMATERIAL_ISOTROPICLINEARHARDENINGPLASTICITY_BEHAVIOUR_DATA_HXX */