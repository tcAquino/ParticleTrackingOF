/**
* \file PTOF/Steppers.h
* \author Tomás Aquino
* \date 09/03/2022
*/

#ifndef PTOF_STEPPERS_H
#define PTOF_STEPPERS_H

#include <cstddef>
#include <type_traits>
#include <utility>
#include "CTRW/JumpGenerator.h"
#include "CTRW/TimeGenerator.h"

namespace ptof
{
  /** \struct Steppers_Advection_RK4_Diffusion_Euler PTOF/Steppers.h "PTOF/Steppers.h"
   * \brief  Time steppers for RK4 advection
   * and stochastic forward Euler diffusion. */
  struct Steppers_Advection_RK4_Diffusion_Euler
  {
    /** Make constant time step TimeGenerator. */
    template <typename Parameters>
    static auto makeTimeGenerator(Parameters const& params)
    {
      if constexpr (has_time_step<Parameters>::value)
        return ctrw::TimeGenerator_Step{ params.time_step };
      return ctrw::TimeGenerator_Step{ 0. };
    }
    
    /** Make advection--diffusion JumpGenerator.
     * Given:
     * - velocity field;
     * - boundary enforcer;
     * - transport parameters; 
     * - solver parameters;
     * - spatial dimension. */
    template
    <typename VelocityField,
    typename Boundary,
    typename TransportParameters,
    typename SolverParameters>
    static auto makeJumpGenerator
    (VelocityField&& velocity_field,
     Boundary&& boundary,
     TransportParameters const& params_transport,
     SolverParameters const& params_solvers,
     std::size_t dim)
    {
      if constexpr (has_time_step<SolverParameters>::value)
        return
          ctrw::JumpGenerator_Add{
            ctrw::JumpGenerator_Velocity_RK4{
              std::forward<VelocityField>(velocity_field),
              params_solvers.time_step,
              std::forward<Boundary>(boundary) },
            ctrw::JumpGenerator_Diffusion{
              params_transport.diff_coeff,
              params_solvers.time_step,
              dim } };
      return
        ctrw::JumpGenerator_Add{
          ctrw::JumpGenerator_Velocity_RK4{
            std::forward<VelocityField>(velocity_field),
            0.,
            std::forward<Boundary>(boundary) },
          ctrw::JumpGenerator_Diffusion{
            params_transport.diff_coeff,
            0.,
            dim } };
    }
    
    /** Make purely-advective JumpGenerator.
     * Given:
     * - velocity field; 
     * - boundary enforcer;
     * - transport parameters;
     * - solver parameters;
     * - spatial dimension. */
    template
    <typename VelocityField,
    typename Boundary,
    typename TransportParameters,
    typename SolverParameters>
    static auto makeJumpGenerator_Advection
    (VelocityField&& velocity_field,
     Boundary&& boundary,
     TransportParameters const& params_transport,
     SolverParameters const& params_solvers,
     std::size_t dim)
    {
      if constexpr (has_time_step<SolverParameters>::value)
        return
          ctrw::JumpGenerator_Velocity_RK4{
            std::forward<VelocityField>(velocity_field),
            params_solvers.time_step,
            std::forward<Boundary>(boundary) };
      return
        ctrw::JumpGenerator_Velocity_RK4{
          std::forward<VelocityField>(velocity_field),
          0.,
          std::forward<Boundary>(boundary) };
    }
    
  private:
    // Check if type T has member double time_step
    // Adapted from kispaljr's answer here:
    // https://stackoverflow.com/questions/257288/templated-check-for-the-existence-of-a-class-member-function
    template <typename T> struct has_time_step
    {
        typedef char (&Yes)[1];
        typedef char (&No)[2];

        template<class U>
        static Yes test(U* data,
                        typename std::enable_if<std::is_same<
                          double,
                          decltype(data->time_step)>::value>::type* = 0);
        static No test(...);
        static const bool value = sizeof(Yes) == sizeof(has_time_step::test((typename std::remove_reference<T>::type*)0));
    };
  };
  
  /** \struct Steppers_Advection_Euler_Diffusion_Euler PTOF/Steppers.h "PTOF/Steppers.h"
   * \brief  Time steppers for forward Euler advection
   * and stochastic forward Euler diffusion. */
  struct Steppers_Advection_Euler_Diffusion_Euler
  {
    // Make constant time step TimeGenerator
    template <typename Parameters>
    static auto makeTimeGenerator(Parameters const& params)
    {
//      if constexpr (std::is_member_pointer_v<decltype(&Parameters::time_step)>)
//        return ctrw::TimeGenerator_Step{ params.time_step };
//      else
        return ctrw::TimeGenerator_Step{ 0. };
    }
    
    /** Make advection--diffusion JumpGenerator.
     * Given:
     * - velocity field;
     * - boundary enforcer;
     * - transport parameters; 
     * - solver parameters;
     * - spatial dimension. */
    template
    <typename VelocityField,
    typename Boundary,
    typename TransportParameters,
    typename SolverParameters>
    static auto makeJumpGenerator
    (VelocityField&& velocity_field,
     Boundary&& boundary,
     TransportParameters const& params_transport,
     SolverParameters const& params_solvers,
     std::size_t dim)
    {
      if constexpr (has_time_step<SolverParameters>::value)
        return
          ctrw::JumpGenerator_Add{
            ctrw::JumpGenerator_Velocity{
              std::forward<VelocityField>(velocity_field),
              params_solvers.time_step },
            ctrw::JumpGenerator_Diffusion{
              params_transport.diff_coeff,
              params_solvers.time_step,
              dim } };
      return
        ctrw::JumpGenerator_Add{
          ctrw::JumpGenerator_Velocity{
            std::forward<VelocityField>(velocity_field),
            0. },
          ctrw::JumpGenerator_Diffusion{
            params_transport.diff_coeff,
            0.,
            dim } };
    }
    
    /** Make advection JumpGenerator.
     * Given:
     * - velocity field;
     * - boundary enforcer;
     * - transport parameters;
     * - solver parameters;
     * - spatial dimension. */
    template
    <typename VelocityField,
    typename Boundary,
    typename TransportParameters,
    typename SolverParameters>
    static auto makeJumpGenerator_Advection
    (VelocityField&& velocity_field,
     Boundary&& boundary,
     TransportParameters const& params_transport,
     SolverParameters const& params_solvers,
     std::size_t dim)
    {
      if constexpr (has_time_step<SolverParameters>::value)
        return
          ctrw::JumpGenerator_Velocity{
            std::forward<VelocityField>(velocity_field),
            params_solvers.time_step };
      return
        ctrw::JumpGenerator_Velocity{
          std::forward<VelocityField>(velocity_field),
          0. };
    }
    
  private:
    // Check if type T has member double time_step
    // Adapted from kispaljr's answer here:
    // https://stackoverflow.com/questions/257288/templated-check-for-the-existence-of-a-class-member-function
    template <typename T> struct has_time_step
    {
        typedef char (&Yes)[1];
        typedef char (&No)[2];

        template<class U>
        static Yes test(U* data,
                        typename std::enable_if<std::is_same<
                          double,
                          decltype(data->time_step)>::value>::type* = 0);
        static No test(...);
        static const bool value = sizeof(Yes) == sizeof(has_time_step::test((typename std::remove_reference<T>::type*)0));
    };
  };
}


#endif /* PTOF_STEPPERS_H */
