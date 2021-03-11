/*!
* \file   FiniteStrainIsotropicLinearHardeningPlasticity2-generic.hxx
* \brief  This file declares the umat interface for the FiniteStrainIsotropicLinearHardeningPlasticity2 behaviour law
* \author Thomas Helfer
* \date   09 / 03 / 2021
*/

#ifndef LIB_GENERIC_FINITESTRAINISOTROPICLINEARHARDENINGPLASTICITY2_HXX
#define LIB_GENERIC_FINITESTRAINISOTROPICLINEARHARDENINGPLASTICITY2_HXX

#include"TFEL/Config/TFELConfig.hxx"
#include"MFront/GenericBehaviour/BehaviourData.h"

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

#ifdef __cplusplus
extern "C"{
#endif /* __cplusplus */

MFRONT_SHAREDOBJ void
FiniteStrainIsotropicLinearHardeningPlasticity2_setOutOfBoundsPolicy(const int);

MFRONT_SHAREDOBJ int
FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_setParameter(const char *const,const double);

MFRONT_SHAREDOBJ int
FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress_setUnsignedShortParameter(const char *const,const unsigned short);

MFRONT_SHAREDOBJ int
FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_setParameter(const char *const,const double);

MFRONT_SHAREDOBJ int
FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress_setUnsignedShortParameter(const char *const,const unsigned short);

MFRONT_SHAREDOBJ int
FiniteStrainIsotropicLinearHardeningPlasticity2_setParameter(const char *const,const double);

MFRONT_SHAREDOBJ int
FiniteStrainIsotropicLinearHardeningPlasticity2_setUnsignedShortParameter(const char *const,const unsigned short);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStrain(mfront_gb_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_AxisymmetricalGeneralisedPlaneStress(mfront_gb_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_Axisymmetrical(mfront_gb_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStress(mfront_gb_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_PlaneStrain(mfront_gb_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_GeneralisedPlaneStrain(mfront_gb_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int FiniteStrainIsotropicLinearHardeningPlasticity2_Tridimensional(mfront_gb_BehaviourData* const);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* LIB_GENERIC_FINITESTRAINISOTROPICLINEARHARDENINGPLASTICITY2_HXX */
