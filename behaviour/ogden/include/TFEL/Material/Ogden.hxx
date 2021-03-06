/*!
* \file   TFEL/Material/Ogden.hxx
* \brief  this file implements the Ogden Behaviour.
*         File generated by tfel version 4.0.0-dev
* \author Thomas Helfer
* \date   20 / 12 / 2016
 */

#ifndef LIB_TFELMATERIAL_OGDEN_HXX
#define LIB_TFELMATERIAL_OGDEN_HXX

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
#include"TFEL/Material/OgdenBehaviourData.hxx"
#include"TFEL/Material/OgdenIntegrationData.hxx"

#include"TFEL/Math/tensor.hxx"
#include "MFront/GenericBehaviour/State.hxx"
#include "MFront/GenericBehaviour/BehaviourData.hxx"
namespace tfel{

namespace material{

struct OgdenParametersInitializer
{
static OgdenParametersInitializer&
get();

double alpha;
double mu;
double K;
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

OgdenParametersInitializer();

OgdenParametersInitializer(const OgdenParametersInitializer&);

OgdenParametersInitializer&
operator=(const OgdenParametersInitializer&);
/*!
 * \brief read the parameters from the given file
 * \param[out] pi : parameters initializer
 * \param[in]  fn : file name
 */
static void readParameters(OgdenParametersInitializer&,const char* const);
};

//! \brief forward declaration
template<ModellingHypothesis::Hypothesis,typename Type,bool use_qt>
class Ogden;

//! \brief forward declaration
template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
std::ostream&
 operator <<(std::ostream&,const Ogden<hypothesis,Type,false>&);

/*!
* \class Ogden
* \brief This class implements the Ogden behaviour.
* \param hypothesis, modelling hypothesis.
* \param Type, numerical type.
* \author Thomas Helfer
* \date   20 / 12 / 2016
* 
*/
template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
class Ogden<hypothesis,Type,false> final
: public MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>,
public OgdenBehaviourData<hypothesis,Type,false>,
public OgdenIntegrationData<hypothesis,Type,false>
{

static constexpr unsigned short N = ModellingHypothesisToSpaceDimension<hypothesis>::value;

static_assert(N==1||N==2||N==3);
static_assert(tfel::typetraits::IsFundamentalNumericType<Type>::cond);
static_assert(tfel::typetraits::IsReal<Type>::cond);

friend std::ostream& operator<< <>(std::ostream&,const Ogden&);

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

public :

typedef OgdenBehaviourData<hypothesis,Type,false> BehaviourData;
typedef OgdenIntegrationData<hypothesis,Type,false> IntegrationData;
typedef typename MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::SMFlag SMFlag;
typedef typename MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::SMType SMType;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::ELASTIC;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::SECANTOPERATOR;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::TANGENTOPERATOR;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::CONSISTENTTANGENTOPERATOR;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::NOSTIFFNESSREQUESTED;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::DSIG_DF;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::DSIG_DDF;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::C_TRUESDELL;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::SPATIAL_MODULI;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::C_TAU_JAUMANN;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::ABAQUS;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::DSIG_DDE;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::DTAU_DF;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::DTAU_DDF;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::DS_DF;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::DS_DDF;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::DS_DC;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::DS_DEGL;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::DT_DELOG;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::DPK1_DF;
using IntegrationResult = typename MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::IntegrationResult;

using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::SUCCESS;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::FAILURE;
using MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::UNRELIABLE_RESULTS;

using StressFreeExpansionType = StrainStensor;

private :




#line 14 "ogden.mfront"
StiffnessTensor dS_dC;

#line 8 "ogden.mfront"
real alpha;
#line 10 "ogden.mfront"
real mu;
#line 12 "ogden.mfront"
real K;
real minimal_time_step_scaling_factor;
real maximal_time_step_scaling_factor;

//! Tangent operator;
TangentOperator Dt;
/*!
* \brief Update internal variables at end of integration
*/
void updateIntegrationVariables()
{}

/*!
* \brief Update internal variables at end of integration
*/
void updateStateVariables()
{}

/*!
* \brief Update auxiliary state variables at end of integration
*/
void updateAuxiliaryStateVariables()
{}

//! \brief Default constructor (disabled)
Ogden() =delete ;
//! \brief Copy constructor (disabled)
Ogden(const Ogden&) = delete;
//! \brief Assignement operator (disabled)
Ogden& operator = (const Ogden&) = delete;

public:

/*!
* \brief Constructor
*/
Ogden(const OgdenBehaviourData<hypothesis,Type,false>& src1,
const OgdenIntegrationData<hypothesis,Type,false>& src2)
: OgdenBehaviourData<hypothesis,Type,false>(src1),
OgdenIntegrationData<hypothesis,Type,false>(src2)
{
using namespace std;
using namespace tfel::math;
using std::vector;
this->alpha = OgdenParametersInitializer::get().alpha;
this->mu = OgdenParametersInitializer::get().mu;
this->K = OgdenParametersInitializer::get().K;
this->minimal_time_step_scaling_factor = OgdenParametersInitializer::get().minimal_time_step_scaling_factor;
this->maximal_time_step_scaling_factor = OgdenParametersInitializer::get().maximal_time_step_scaling_factor;
}

/*
 * \brief constructor for the Generic interface
 * \param[in] mgb_d: behaviour data
 */
Ogden(const mfront::gb::BehaviourData& mgb_d)
: OgdenBehaviourData<hypothesis,Type,false>(mgb_d),
OgdenIntegrationData<hypothesis,Type,false>(mgb_d)
{
using namespace std;
using namespace tfel::math;
using std::vector;
this->alpha = OgdenParametersInitializer::get().alpha;
this->mu = OgdenParametersInitializer::get().mu;
this->K = OgdenParametersInitializer::get().K;
this->minimal_time_step_scaling_factor = OgdenParametersInitializer::get().minimal_time_step_scaling_factor;
this->maximal_time_step_scaling_factor = OgdenParametersInitializer::get().maximal_time_step_scaling_factor;
tfel::fsalgo::copy<TensorSize >::exe(mgb_d.s0.gradients,this->F0.begin());
tfel::fsalgo::copy<TensorSize >::exe(mgb_d.s1.gradients,this->F1.begin());
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

IntegrationResult computePredictionOperator(const SMFlag,const SMType) override{
tfel::raise("Ogden::computePredictionOperator: "
"unsupported prediction operator flag");
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
bool computeTangentOperator_ = smt!=NOSTIFFNESSREQUESTED;
#line 17 "ogden.mfront"
const auto id  = Stensor::Id();
#line 18 "ogden.mfront"
const auto J   = det(this->F1);
#line 19 "ogden.mfront"
const auto C   = computeRightCauchyGreenTensor(this->F1);
#line 21 "ogden.mfront"
const auto C2       = square(C);
#line 22 "ogden.mfront"
const auto I1       = trace(C);
#line 23 "ogden.mfront"
const auto I2       = (I1*I1-trace(C2))/2;
#line 24 "ogden.mfront"
const auto I3       = J*J;
#line 25 "ogden.mfront"
const auto dI3_dC   = C2-I1*C+I2*id;
#line 28 "ogden.mfront"
const auto dPv_dJ   = this->K*(J-1);
#line 29 "ogden.mfront"
const StressStensor Sv = dPv_dJ/J*dI3_dC;
#line 33 "ogden.mfront"
const auto iJb        =  1/cbrt(I3);
#line 34 "ogden.mfront"
const auto iJb2       =  power<2>(iJb);
#line 35 "ogden.mfront"
const auto iJb4       =  iJb2*iJb2;
#line 36 "ogden.mfront"
const auto iJb7       =  iJb4*power<3>(iJb);
#line 37 "ogden.mfront"
const auto diJb_dI3   = -iJb4/3;
#line 38 "ogden.mfront"
const auto d2iJb_dI32 = 4*iJb7/9;
#line 39 "ogden.mfront"
const auto diJb_dC    = diJb_dI3*dI3_dC;
#line 41 "ogden.mfront"
Stensor n0,n1,n2;
#line 42 "ogden.mfront"
tvector<3u,real> vp;
#line 43 "ogden.mfront"
tmatrix<3u,3u,real> m;
#line 44 "ogden.mfront"
C.computeEigenVectors(vp,m);
#line 45 "ogden.mfront"
Stensor::computeEigenTensors(n0,n1,n2,m);
#line 46 "ogden.mfront"
const auto a = this->alpha/2;
#line 47 "ogden.mfront"
const tvector<3u,real> pwv = {pow(vp(0),a-2),pow(vp(1),a-2),pow(vp(2),a-2)};
#line 48 "ogden.mfront"
const tvector<3u,real> dfv = {a*vp(0)*pwv(0),a*vp(1)*pwv(1),a*vp(2)*pwv(2)};
#line 49 "ogden.mfront"
const auto fv  = vp(0)*vp(0)*pwv(0)+vp(1)*vp(1)*pwv(1)+vp(2)*vp(2)*pwv(2);
#line 50 "ogden.mfront"
const auto df_dC = dfv(0)*n0+dfv(1)*n1+dfv(2)*n2;
#line 51 "ogden.mfront"
const auto c     = pow(iJb,a-2);
#line 52 "ogden.mfront"
const StressStensor Si = this->mu*c*iJb*((fv*diJb_dC+(iJb/a)*df_dC));
#line 54 "ogden.mfront"
this->sig = convertSecondPiolaKirchhoffStressToCauchyStress(Sv+Si,this->F1);
#line 55 "ogden.mfront"
if(computeTangentOperator_){
#line 56 "ogden.mfront"
const auto d2I3_dC2 = computeDeterminantSecondDerivative(C);
#line 58 "ogden.mfront"
const auto d2Pv_dJ2 = this->K;
#line 59 "ogden.mfront"
this->dS_dC = ((d2Pv_dJ2-dPv_dJ/J)/(2*I3)*(dI3_dC^dI3_dC)+
#line 60 "ogden.mfront"
dPv_dJ/J*d2I3_dC2);
#line 62 "ogden.mfront"
auto df = [&a](const real x){
#line 63 "ogden.mfront"
return a*pow(x,a-1);
#line 64 "ogden.mfront"
};
#line 65 "ogden.mfront"
auto d2f = [&a](const real x){
#line 66 "ogden.mfront"
return a*(a-1)*pow(x,a-2);
#line 67 "ogden.mfront"
};
#line 68 "ogden.mfront"
Stensor4 d2f_dC2;
#line 69 "ogden.mfront"
Stensor::computeIsotropicFunctionDerivative(d2f_dC2,df,d2f,
#line 70 "ogden.mfront"
vp,m,1.e-12);
#line 71 "ogden.mfront"
const auto d2iJb_dC2  =
#line 72 "ogden.mfront"
d2iJb_dI32*(dI3_dC^dI3_dC)+ diJb_dI3*d2I3_dC2;
#line 73 "ogden.mfront"
this->dS_dC += this->mu*c*((a-1)*fv*(diJb_dC^diJb_dC)+
#line 74 "ogden.mfront"
iJb*(fv*d2iJb_dC2+
#line 75 "ogden.mfront"
((diJb_dC^df_dC)+(df_dC^diJb_dC))+
#line 76 "ogden.mfront"
iJb/a*d2f_dC2));
#line 77 "ogden.mfront"
}
this->updateIntegrationVariables();
this->updateStateVariables();
this->updateAuxiliaryStateVariables();
if(computeTangentOperator_){
if(!this->computeConsistentTangentOperator(smflag,smt)){
return MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::FAILURE;
}
}
return MechanicalBehaviour<MechanicalBehaviourBase::STANDARDFINITESTRAINBEHAVIOUR,hypothesis,Type,false>::SUCCESS;
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

bool computeConsistentTangentOperator_DSIG_DF(const SMType smt){
using namespace tfel::math;
// computing DS_DC
this->computeConsistentTangentOperator_DS_DC(smt);
const st2tost2<N,stress> tangentOperator_DS_DC = this->Dt.template get<st2tost2<N,stress> >();
const auto tangentOperator_DS_DF = convert<DS_DF,DS_DC>(tangentOperator_DS_DC,this->F0,this->F1,this->sig);
const auto tangentOperator_DTAU_DF = convert<DTAU_DF,DS_DF>(tangentOperator_DS_DF,this->F0,this->F1,this->sig);
this->Dt = convert<DSIG_DF,DTAU_DF>(tangentOperator_DTAU_DF,this->F0,this->F1,this->sig);
return true;
}

bool computeConsistentTangentOperator_DSIG_DDF(const SMType smt){
using namespace tfel::math;
// computing DS_DC
this->computeConsistentTangentOperator_DS_DC(smt);
const st2tost2<N,stress> tangentOperator_DS_DC = this->Dt.template get<st2tost2<N,stress> >();
const auto tangentOperator_DS_DF = convert<DS_DF,DS_DC>(tangentOperator_DS_DC,this->F0,this->F1,this->sig);
const auto tangentOperator_DTAU_DF = convert<DTAU_DF,DS_DF>(tangentOperator_DS_DF,this->F0,this->F1,this->sig);
const auto tangentOperator_DSIG_DF = convert<DSIG_DF,DTAU_DF>(tangentOperator_DTAU_DF,this->F0,this->F1,this->sig);
this->Dt = convert<DSIG_DDF,DSIG_DF>(tangentOperator_DSIG_DF,this->F0,this->F1,this->sig);
return true;
}

bool computeConsistentTangentOperator_C_TRUESDELL(const SMType smt){
using namespace tfel::math;
// computing DS_DC
this->computeConsistentTangentOperator_DS_DC(smt);
const st2tost2<N,stress> tangentOperator_DS_DC = this->Dt.template get<st2tost2<N,stress> >();
const auto tangentOperator_DS_DEGL = convert<DS_DEGL,DS_DC>(tangentOperator_DS_DC,this->F0,this->F1,this->sig);
const auto tangentOperator_SPATIAL_MODULI = convert<SPATIAL_MODULI,DS_DEGL>(tangentOperator_DS_DEGL,this->F0,this->F1,this->sig);
this->Dt = convert<C_TRUESDELL,SPATIAL_MODULI>(tangentOperator_SPATIAL_MODULI,this->F0,this->F1,this->sig);
return true;
}

bool computeConsistentTangentOperator_SPATIAL_MODULI(const SMType smt){
using namespace tfel::math;
// computing DS_DC
this->computeConsistentTangentOperator_DS_DC(smt);
const st2tost2<N,stress> tangentOperator_DS_DC = this->Dt.template get<st2tost2<N,stress> >();
const auto tangentOperator_DS_DEGL = convert<DS_DEGL,DS_DC>(tangentOperator_DS_DC,this->F0,this->F1,this->sig);
this->Dt = convert<SPATIAL_MODULI,DS_DEGL>(tangentOperator_DS_DEGL,this->F0,this->F1,this->sig);
return true;
}

bool computeConsistentTangentOperator_C_TAU_JAUMANN(const SMType smt){
using namespace tfel::math;
// computing DS_DC
this->computeConsistentTangentOperator_DS_DC(smt);
const st2tost2<N,stress> tangentOperator_DS_DC = this->Dt.template get<st2tost2<N,stress> >();
const auto tangentOperator_DS_DF = convert<DS_DF,DS_DC>(tangentOperator_DS_DC,this->F0,this->F1,this->sig);
const auto tangentOperator_DTAU_DF = convert<DTAU_DF,DS_DF>(tangentOperator_DS_DF,this->F0,this->F1,this->sig);
this->Dt = convert<C_TAU_JAUMANN,DTAU_DF>(tangentOperator_DTAU_DF,this->F0,this->F1,this->sig);
return true;
}

bool computeConsistentTangentOperator_ABAQUS(const SMType smt){
using namespace tfel::math;
// computing DS_DC
this->computeConsistentTangentOperator_DS_DC(smt);
const st2tost2<N,stress> tangentOperator_DS_DC = this->Dt.template get<st2tost2<N,stress> >();
const auto tangentOperator_DS_DEGL = convert<DS_DEGL,DS_DC>(tangentOperator_DS_DC,this->F0,this->F1,this->sig);
this->Dt = convert<ABAQUS,DS_DEGL>(tangentOperator_DS_DEGL,this->F0,this->F1,this->sig);
return true;
}

bool computeConsistentTangentOperator_DSIG_DDE(const SMType){
tfel::raise("Ogden::computeConsistentTangentOperator_DSIG_DDE: "
"computing the tangent operator 'DSIG_DDE' is not supported");
}

bool computeConsistentTangentOperator_DTAU_DF(const SMType smt){
using namespace tfel::math;
// computing DS_DC
this->computeConsistentTangentOperator_DS_DC(smt);
const st2tost2<N,stress> tangentOperator_DS_DC = this->Dt.template get<st2tost2<N,stress> >();
const auto tangentOperator_DS_DF = convert<DS_DF,DS_DC>(tangentOperator_DS_DC,this->F0,this->F1,this->sig);
this->Dt = convert<DTAU_DF,DS_DF>(tangentOperator_DS_DF,this->F0,this->F1,this->sig);
return true;
}

bool computeConsistentTangentOperator_DTAU_DDF(const SMType smt){
using namespace tfel::math;
// computing DS_DC
this->computeConsistentTangentOperator_DS_DC(smt);
const st2tost2<N,stress> tangentOperator_DS_DC = this->Dt.template get<st2tost2<N,stress> >();
const auto tangentOperator_DS_DF = convert<DS_DF,DS_DC>(tangentOperator_DS_DC,this->F0,this->F1,this->sig);
const auto tangentOperator_DTAU_DF = convert<DTAU_DF,DS_DF>(tangentOperator_DS_DF,this->F0,this->F1,this->sig);
this->Dt = convert<DTAU_DDF,DTAU_DF>(tangentOperator_DTAU_DF,this->F0,this->F1,this->sig);
return true;
}

bool computeConsistentTangentOperator_DS_DF(const SMType smt){
using namespace tfel::math;
// computing DS_DC
this->computeConsistentTangentOperator_DS_DC(smt);
const st2tost2<N,stress> tangentOperator_DS_DC = this->Dt.template get<st2tost2<N,stress> >();
this->Dt = convert<DS_DF,DS_DC>(tangentOperator_DS_DC,this->F0,this->F1,this->sig);
return true;
}

bool computeConsistentTangentOperator_DS_DDF(const SMType){
tfel::raise("Ogden::computeConsistentTangentOperator_DS_DDF: "
"computing the tangent operator 'DS_DDF' is not supported");
}

bool computeConsistentTangentOperator_DS_DC(const SMType smt){
using namespace std;
using namespace tfel::math;
using std::vector;
#line 81 "ogden.mfront"
this->Dt = this->dS_dC;
return true;
}

bool computeConsistentTangentOperator_DS_DEGL(const SMType smt){
using namespace tfel::math;
// computing DS_DC
this->computeConsistentTangentOperator_DS_DC(smt);
const st2tost2<N,stress> tangentOperator_DS_DC = this->Dt.template get<st2tost2<N,stress> >();
this->Dt = convert<DS_DEGL,DS_DC>(tangentOperator_DS_DC,this->F0,this->F1,this->sig);
return true;
}

bool computeConsistentTangentOperator_DT_DELOG(const SMType){
tfel::raise("Ogden::computeConsistentTangentOperator_DT_DELOG: "
"computing the tangent operator 'DT_DELOG' is not supported");
}

bool computeConsistentTangentOperator_DPK1_DF(const SMType smt){
using namespace tfel::math;
// computing DS_DC
this->computeConsistentTangentOperator_DS_DC(smt);
const st2tost2<N,stress> tangentOperator_DS_DC = this->Dt.template get<st2tost2<N,stress> >();
const auto tangentOperator_DS_DEGL = convert<DS_DEGL,DS_DC>(tangentOperator_DS_DC,this->F0,this->F1,this->sig);
this->Dt = convert<DPK1_DF,DS_DEGL>(tangentOperator_DS_DEGL,this->F0,this->F1,this->sig);
return true;
}

bool computeConsistentTangentOperator(const SMFlag smflag,const SMType smt){
switch(smflag){
case DSIG_DF:
return this->computeConsistentTangentOperator_DSIG_DF(smt);
case DSIG_DDF:
return this->computeConsistentTangentOperator_DSIG_DDF(smt);
case C_TRUESDELL:
return this->computeConsistentTangentOperator_C_TRUESDELL(smt);
case SPATIAL_MODULI:
return this->computeConsistentTangentOperator_SPATIAL_MODULI(smt);
case C_TAU_JAUMANN:
return this->computeConsistentTangentOperator_C_TAU_JAUMANN(smt);
case ABAQUS:
return this->computeConsistentTangentOperator_ABAQUS(smt);
case DSIG_DDE:
return this->computeConsistentTangentOperator_DSIG_DDE(smt);
case DTAU_DF:
return this->computeConsistentTangentOperator_DTAU_DF(smt);
case DTAU_DDF:
return this->computeConsistentTangentOperator_DTAU_DDF(smt);
case DS_DF:
return this->computeConsistentTangentOperator_DS_DF(smt);
case DS_DDF:
return this->computeConsistentTangentOperator_DS_DDF(smt);
case DS_DC:
return this->computeConsistentTangentOperator_DS_DC(smt);
case DS_DEGL:
return this->computeConsistentTangentOperator_DS_DEGL(smt);
case DT_DELOG:
return this->computeConsistentTangentOperator_DT_DELOG(smt);
case DPK1_DF:
return this->computeConsistentTangentOperator_DPK1_DF(smt);
}
tfel::raise("Ogden::computeConsistentTangentOperator: "
"unsupported tangent operator flag");
}

const TangentOperator& getTangentOperator() const{
return this->Dt;
}

void updateExternalStateVariables(){
this->F0  = this->F1;
this->T += this->dT;
}

//!
~Ogden()
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
}; // end of Ogden class

template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
std::ostream&
operator <<(std::ostream& os,const Ogden<hypothesis,Type,false>& b)
{
os << "F₀ : " << b.F0 << '\n';
os << "F₁ : " << b.F1 << '\n';
os << "σ : " << b.sig << '\n';
os << "Δt : " << b.dt << '\n';
os << "T : " << b.T << '\n';
os << "ΔT : " << b.dT << '\n';
os << "dS_dC : " << b.dS_dC << '\n';
os << "alpha : " << b.alpha << '\n';
os << "mu : " << b.mu << '\n';
os << "K : " << b.K << '\n';
os << "minimal_time_step_scaling_factor : " << b.minimal_time_step_scaling_factor << '\n';
os << "maximal_time_step_scaling_factor : " << b.maximal_time_step_scaling_factor << '\n';
return os;
}

/*!
* Partial specialisation for Ogden.
*/
template<ModellingHypothesis::Hypothesis hypothesis,typename Type>
class MechanicalBehaviourTraits<Ogden<hypothesis,Type,false> >
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
static constexpr unsigned short internal_variables_nb  = 0;
static constexpr unsigned short external_variables_nb  = 1;
static constexpr unsigned short external_variables_nb2 = 0;
static constexpr bool hasConsistentTangentOperator = true;
static constexpr bool isConsistentTangentOperatorSymmetric = false;
static constexpr bool hasPredictionOperator = false;
static constexpr bool hasAPrioriTimeStepScalingFactor = false;
static constexpr bool hasComputeInternalEnergy = false;
static constexpr bool hasComputeDissipatedEnergy = false;
/*!
* \return the name of the class.
*/
static const char* getName(){
return "Ogden";
}

};

/*!
* Partial specialisation for Ogden.
*/
template<typename Type>
class MechanicalBehaviourTraits<Ogden<ModellingHypothesis::AXISYMMETRICALGENERALISEDPLANESTRESS,Type,false> >
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
return "Ogden";
}

};

/*!
* Partial specialisation for Ogden.
*/
template<typename Type>
class MechanicalBehaviourTraits<Ogden<ModellingHypothesis::PLANESTRESS,Type,false> >
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
return "Ogden";
}

};

} // end of namespace material

} // end of namespace tfel

#endif /* LIB_TFELMATERIAL_OGDEN_HXX */
