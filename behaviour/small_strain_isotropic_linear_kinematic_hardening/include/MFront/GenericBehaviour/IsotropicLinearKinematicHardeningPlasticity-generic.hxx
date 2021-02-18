/*!
* \file   IsotropicLinearKinematicHardeningPlasticity-generic.hxx
* \brief  This file declares the umat interface for the IsotropicLinearKinematicHardeningPlasticity behaviour law
* \author Thomas Helfer
* \date   14 / 10 / 2016
*/

#ifndef LIB_GENERIC_ISOTROPICLINEARKINEMATICHARDENINGPLASTICITY_HXX
#define LIB_GENERIC_ISOTROPICLINEARKINEMATICHARDENINGPLASTICITY_HXX

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
IsotropicLinearKinematicHardeningPlasticity_setOutOfBoundsPolicy(const int);

MFRONT_SHAREDOBJ int
IsotropicLinearKinematicHardeningPlasticity_setParameter(const char *const,const double);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_AxisymmetricalGeneralisedPlaneStrain(mfront_gb_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_Axisymmetrical(mfront_gb_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_PlaneStrain(mfront_gb_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_GeneralisedPlaneStrain(mfront_gb_BehaviourData* const);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int IsotropicLinearKinematicHardeningPlasticity_Tridimensional(mfront_gb_BehaviourData* const);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* LIB_GENERIC_ISOTROPICLINEARKINEMATICHARDENINGPLASTICITY_HXX */
