/*!
* \file   TFEL/Material/OgdenBehaviourData.hxx
* \brief  this file implements the OgdenBehaviourData class.
*         File generated by tfel version 4.0.0-dev
* \author Thomas Helfer
* \date   20 / 12 / 2016
 */

#ifndef LIB_TFELMATERIAL_OGDEN_BEHAVIOUR_DATA_HXX
#define LIB_TFELMATERIAL_OGDEN_BEHAVIOUR_DATA_HXX

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
#include"TFEL/Math/tensor.hxx"
#include"TFEL/Math/Tensor/TensorConceptIO.hxx"
#include"TFEL/Math/t2tot2.hxx"
#include"TFEL/Math/T2toT2/T2toT2ConceptIO.hxx"
#include"TFEL/Math/t2tost2.hxx"
#include"TFEL/Math/T2toST2/T2toST2ConceptIO.hxx"
#include"TFEL/Math/st2tot2.hxx"
#include"TFEL/Math/ST2toT2/ST2toT2ConceptIO.hxx"
#include"TFEL/Math/ST2toST2/ConvertToTangentModuli.hxx"
#include"TFEL/Math/ST2toST2/ConvertSpatialModuliToKirchhoffJaumanRateModuli.hxx"
#include"TFEL/Material/FiniteStrainBehaviourTangentOperator.hxx"
#include"TFEL/Material/ModellingHypothesis.hxx"

#include "MFront/GenericBehaviour/State.hxx"
#include "MFront/GenericBehaviour/BehaviourData.hxx"
namespace tfel{

namespace material{

//! \brief forward declaration
template<ModellingHypothesis::Hypothesis hypothesis,typename,bool>
class OgdenBehaviourData;

//! \brief forward declaration
template<ModellingHypothesis::Hypothesis hypothesis,typename Type,bool use_qt>
class OgdenIntegrationData;

//! \brief forward declaration
template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
std::ostream&
 operator <<(std::ostream&,const OgdenBehaviourData<hypothesis,Type,false>&);

template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
class OgdenBehaviourData<hypothesis,Type,false>
{

static constexpr unsigned short N = ModellingHypothesisToSpaceDimension<hypothesis>::value;
static_assert(N==1||N==2||N==3);
static_assert(tfel::typetraits::IsFundamentalNumericType<Type>::cond);
static_assert(tfel::typetraits::IsReal<Type>::cond);

friend std::ostream& operator<< <>(std::ostream&,const OgdenBehaviourData&);

/* integration data is declared friend to access   driving variables at the beginning of the time step */
friend class OgdenIntegrationData<hypothesis,Type,false>;

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
using TangentOperator   = FiniteStrainBehaviourTangentOperator<N,stress>;
using PhysicalConstants = tfel::PhysicalConstants<real>;

protected:

DeformationGradientTensor F0;

StressStensor sig;


temperature T;

public:

/*!
* \brief Default constructor
*/
OgdenBehaviourData()
{}

/*!
* \brief copy constructor
*/
OgdenBehaviourData(const OgdenBehaviourData& src)
: F0(src.F0),
sig(src.sig),
T(src.T)
{}

/*
 * \brief constructor for the Generic interface
 * \param[in] mgb_d: behaviour data
 */
OgdenBehaviourData(const mfront::gb::BehaviourData& mgb_d)
: T(mgb_d.s0.external_state_variables[0])
{
}


/*
* \brief Assignement operator
*/
OgdenBehaviourData&
operator=(const OgdenBehaviourData& src){
this->F0 = src.F0;
this->sig = src.sig;
this->T = src.T;
return *this;
}

void exportStateData(mfront::gb::State& mbg_s1) const
{
using namespace tfel::math;
tfel::fsalgo::copy<StensorSize>::exe(this->sig.begin(), mbg_s1.thermodynamic_forces);
} // end of exportStateData

}; // end of OgdenBehaviourDataclass

template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
std::ostream&
operator <<(std::ostream& os,const OgdenBehaviourData<hypothesis,Type,false>& b)
{
os << "F₀ : " << b.F0 << '\n';
os << "σ : " << b.sig << '\n';
os << "T : " << b.T << '\n';
return os;
}

} // end of namespace material

} // end of namespace tfel

#endif /* LIB_TFELMATERIAL_OGDEN_BEHAVIOUR_DATA_HXX */
