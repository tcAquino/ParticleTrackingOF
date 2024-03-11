/**
 \file CTRW/JumpGenerator.h
 \author Tomás Aquino
 \date 20/01/2018
 
 \brief Jump generators for particle tracking.
 
  A JumpGenerator should implement the following minimal functionality:
 
 \code{.cpp}
  class JumpGenerator
  {
    // Return the jump increment
    // with type of state element to be incremented
    //   template <typename State>
    //   auto operator() (State const& state)
    {
      // Return the jump increment
    }
  };
 \endcode
*/

#ifndef CTRW_JUMPGENERATOR_H
#define CTRW_JUMPGENERATOR_H

#include <algorithm>
#include <cmath>
#include <random>
#include <type_traits>
#include <utility>
#include "General/Constants.h"
#include "General/Operations.h"
#include "General/Useful.h"
#include "Geometry/Boundary.h"

namespace ctrw
{
  /** \class JumpGenerator_Step CTRW/JumpGenerator.h "CTRW/JumpGenerator.h"
   * \brief Jump a fixed step. */
  template <typename Step = double>
  class JumpGenerator_Step
  {
  public:
    /** Constructor.
     \param step Fixed jump step.
    */
    JumpGenerator_Step(Step step)
    : step{ step }
    {}
    
    /**
     \param state Particle state (unused).
     \return Jump increment.
    */
    template <typename State = useful::Empty>
    auto operator() (State const& state = {})
    { return step; }
    
    const Step step;  /**< The fixed step to jump. */
  };
  
  /** \class JumpGenerator_Add CTRW/JumpGenerator.h "CTRW/JumpGenerator.h"
   * \brief Jump according to the sum of two jump generators. */
  template <typename JumpGenerator_1, typename JumpGenerator_2>
  class JumpGenerator_Add
  {
  public:
    /** Constructor.
     \param jump_generator_1 First jump generator.
     \param jump_generator_2 Second jump generator.
    */
    JumpGenerator_Add
    (JumpGenerator_1 jump_generator_1,
     JumpGenerator_2 jump_generator_2)
    : _jump_generator_1{ jump_generator_1 }
    , _jump_generator_2{ jump_generator_2 }
    {}
    
    /**
     \param state Particle state (unused).
     \return Jump increment.
    */
    template <typename State = useful::Empty>
    auto operator() (State const& state = {})
    {
      return operation::plus(_jump_generator_1(state),
                             _jump_generator_2(state));
    }
    
    /** \return First jump generator. */
    auto const& jump_generator_1() const
    { return _jump_generator_1; }
    
    /** \return Second jump generator. */
    auto const& jump_generator_2() const
    { return _jump_generator_2; }
    
    /** \param time_step Time step to set. */
    void time_step(double time_step)
    {
      if constexpr (useful::has_time_step_setter<JumpGenerator_1>::value)
        _jump_generator_1.time_step(time_step);
      if constexpr (useful::has_time_step_setter<JumpGenerator_2>::value)
        _jump_generator_2.time_step(time_step);
    }
    
  private:
    JumpGenerator_1 _jump_generator_1;    /**< First jump generator. */
    JumpGenerator_2 _jump_generator_2;    /**< Second jump generator. */
  };
  
  /** \class JumpGenerator_Velocity CTRW/JumpGenerator.h "CTRW/JumpGenerator.h"
    \brief Jump according to a velocity field and time step (forward Euler).
  */
  template <typename VelocityField>
  class JumpGenerator_Velocity
  {
  public:
    /** Constructor.
      \param velocity_field Velocity field (as a function of state)
      \param time_step Time step.
    */
    JumpGenerator_Velocity
    (VelocityField&& velocity_field, double time_step)
    : _velocity_field{ std::forward<VelocityField>(velocity_field) }
    , _time_step{ time_step }
    {}

    /** \param time_step Time step to set. */
    void time_step(double time_step)
    { _time_step = time_step; }

    /** \return Time step. */
    double time_step() const
    { return _time_step; }

    /**
     \param state Particle state.
     \return Jump increment.
    */
    template <typename State>
    auto operator() (State const& state)
    {
      return operation::times_scalar(_time_step, velocity(state));
    }
    
    /**
     \param state Particle state.
     \return Velocity value.
    */
    template <typename State>
    auto velocity(State const& state) const
    {
      return _velocity_field(state);
    }
    
    /** \return Velocity field. */
    auto const& velocity_field() const
    { return _velocity_field; }
    
  private:
    VelocityField _velocity_field;  /**< Velocity as a function of state. */
    double _time_step;              /**< Time step.                       */
  };
  template <typename VelocityField>
  JumpGenerator_Velocity(VelocityField&&, double)
  -> JumpGenerator_Velocity<VelocityField>;
  
  /** \class JumpGenerator_Velocity_RK4 CTRW/JumpGenerator.h "CTRW/JumpGenerator.h"
    \brief Jump according to a velocity field and time step, using RK4 scheme.
    \details Boundary conditions are enforced in predictor-corrector steps.
    \note Particle state must define:
      - position */
  template <typename VelocityField, typename Boundary = geometry::Boundary_DoNothing>
  class JumpGenerator_Velocity_RK4
  {
  public:
    /** Constructor.
     \param velocity_field the time step, and the boundary to enforce boundary conditions.
     \param time_step Time step.
     \param boundary Object to apply boundary conditions to current state given current state and previous state.
     */
    JumpGenerator_Velocity_RK4
    (VelocityField&& velocity_field, double time_step, Boundary&& boundary = {})
    : _velocity_field{ std::forward<VelocityField>(velocity_field) }
    , _time_step{ time_step }
    , _boundary{ std::forward<Boundary>(boundary) }
    {}

    /** \param time_step Time step to set. */
    void time_step(double time_step)
    { _time_step = time_step; }

    /** \return Time step. */
    double time_step() const
    { return _time_step; }

    /**
     \param state Particle state.
     \return Jump increment.
    */
    template <typename State>
    auto operator() (State const& state)
    {
      auto state_intermediate = state;

      auto k1 = velocity(state);
      operation::linearOp(_time_step/2., k1, state.position,
                          state_intermediate.position);
      boundary(state_intermediate, state);
      
      auto k2 = velocity(state_intermediate);
      operation::linearOp(_time_step/2., k2, state.position,
                          state_intermediate.position);
      boundary(state_intermediate, state);
      
      auto k3 = velocity(state_intermediate);
      operation::linearOp(_time_step, k3, state.position,
                          state_intermediate.position);
      boundary(state_intermediate, state);

      auto k4 = velocity(state_intermediate);
      auto jump = operation::plus(k1, k4);
      operation::linearOp(2., k2, jump, jump);
      operation::linearOp(2., k3, jump, jump);
      operation::times_scalar_InPlace(_time_step/6., jump);
      
      return jump;
    }
    
    /**
     \param state Particle state.
     \return Velocity value.
    */
    template <typename State>
    auto velocity(State const& state) const
    {
      return _velocity_field(state);
    }
    
    /** \return Velocity field. */
    auto const& get_velocity_field() const
    { return _velocity_field; }
    
  private:
    VelocityField _velocity_field;  /**< Velocity as a function of state.*/
    double _time_step;              /**< Time step.                      */
    Boundary _boundary;             /**< Boundary conditions.             */
  };
  template <typename VelocityField, typename Boundary>
  JumpGenerator_Velocity_RK4
  (VelocityField&&, double, Boundary&&)
  -> JumpGenerator_Velocity_RK4<VelocityField, Boundary>;
  
  /** \class JumpGenerator_Diffusion_1d CTRW/JumpGenerator.h "CTRW/JumpGenerator.h"
   * \brief One dimensional diffusion jumps using Gaussian distribution. */
  template <typename RNG = std::mt19937>
  class JumpGenerator_Diffusion_1d
  {
    double _time_step;          /**< Time step.            */
    double _diff_coeff;        /**< Diffusion coefficient.*/
    
    double _diff_aux{ std::sqrt(2.*_diff_coeff*_time_step) };   /**< Typicall jump size. */
    RNG _rng{ std::random_device{}() };                         /**< Random number generator. */
    std::normal_distribution<double> _normal_dist{ 0., 1. };    /**< Unit normal distribution.*/

  public:

    /** Constructor.
     \param diff_coeff Diffusion coefficient.
     \param time_step Diffusion coefficient.
    */
    JumpGenerator_Diffusion_1d(double diff_coeff, double time_step)
    : _time_step{ time_step }
    , _diff_coeff{ diff_coeff }
    {}

    /** \param time_step Time step to set. */
    void time_step(double time_step)
    {
      _time_step = time_step;
      _diff_aux = std::sqrt(2.*_diff_coeff*_time_step);
    }

    /** \param diff_coeff Diffusion coefficient to set. */
    void diff(double diff_coeff)
    {
      _diff_coeff = diff_coeff;
      _diff_aux = std::sqrt(2.*_diff_coeff*_time_step);
    }

    /** \return Time step. */
    double time_step() const
    { return _time_step; }
    
    /** \return Diffusion coefficient. */
    double diff() const
    { return _diff_coeff; }

    /**
     \param state Particle state (unused).
     \return Jump increment.
    */
    template <typename State = useful::Empty>
    double operator() (State const& state = {})
    {
      return _diff_aux*_normal_dist(_rng);
    }
  };
  
  /** \class JumpGenerator_Diffusion CTRW/JumpGenerator.h "CTRW/JumpGenerator.h"
   * \brief Diffusion jumps along each dimension using Gaussian distributions. */
  class JumpGenerator_Diffusion
  {
  public:
    /** Constructor.
     \param diff_coeffs Diffusion coefficients for each dimension.
     \param time_step Diffusion coefficient.
    */
    JumpGenerator_Diffusion
    (std::vector<double> const& diff_coeffs, double time_step)
    {
      _jump_generator.reserve(diff_coeffs.size());
      for (auto diff_coeff : diff_coeffs)
        _jump_generator.emplace_back(diff_coeff, time_step);
    }
    
    /** Constructor.
     \param diff_coeff Diffusion coefficient (same for all dimensions).
     \param time_step Diffusion coefficient.
     \param dim Number of dimensions.
    */
    JumpGenerator_Diffusion
    (double diff_coeff, double time_step, std::size_t dim)
    {
      _jump_generator.reserve(dim);
      for (std::size_t dd = 0; dd < dim; ++dd)
        _jump_generator.emplace_back(diff_coeff, time_step);
    }
    
    /** \param time_step Time step to set. */
    void time_step(double time_step)
    {
      for (std::size_t dd = 0; dd < dim(); ++dd)
        _jump_generator[dd].time_step(time_step);
    }

    /** \param diff_coeff Diffusion coefficient to set (same for all dimensions). */
    void diff(double diff_coeff)
    {
      for (std::size_t dd = 0; dd < dim(); ++dd)
        _jump_generator[dd].diff(diff_coeff);
    }
    
    /** \param diff_coeffs Diffusion coefficients to set for each dimension. */
    void diff(std::vector<double> const& diff_coeffs)
    {
      for (std::size_t dd = 0; dd < dim(); ++dd)
        _jump_generator[dd].diff(diff_coeffs[dd]);
    }
    
    /** \return Time step. */
    double time_step() const
    { return _jump_generator[0].time_step(); }
    
    /**
     \param dd Dimension.
     \return Diffusion coefficient along dimension. */
    double diff(std::size_t dd) const
    { return _jump_generator[dd].diff(); }
    
    /** \return Number of spatial dimensions. */
    std::size_t dim() const
    { return _jump_generator.size(); }
    
    /**
     \param state Particle state (unused).
     \return Jump increment.
    */
    template <typename State = useful::Empty>
    auto operator() (State const& state = {})
    {
      auto jump = state.position;
      for (std::size_t dd = 0; dd < dim(); ++dd)
        jump[dd] = _jump_generator[dd](state);
      
      return jump;
    }
    
  private:
    std::vector<JumpGenerator_Diffusion_1d<>> _jump_generator;  /**< Diffusion jump generators along each dimension.*/
  };
    
  /** \class JumpGenerator_RW_1d CTRW/JumpGenerator.h "CTRW/JumpGenerator.h"
   * \brief One-dimensional random walk jumps with fixed step sizes to left or right. */
  template <typename Distance_type = double, typename RNG = std::mt19937>
  class JumpGenerator_RW_1d
  {
    const Distance_type _jump_size;                             /**< Fixed jump size. */
    const double _prob_right;                                   /**< Probability to jump right. */
    RNG _rng{ std::random_device{}() };                         /**< Random number generator. */
    std::bernoulli_distribution _bernoulli_dist{ _prob_right }; /**< Bernoulli distribution. */

  public:
    /** Constructor.
     \param jump_size Fixed jump size.
     \param prob_right Probability of jumping to the right.
    */
    JumpGenerator_RW_1d
    (Distance_type jump_size = 1, double prob_right = 0.5)
    : _jump_size{ jump_size }
    , _prob_right{ prob_right }
    {}

    /**
     \param state Particle state (unused).
     \return Jump increment.
    */
    template <typename State = useful::Empty>
    double operator() (State const& state = {})
    { return _jump_size*(2*_bernoulli_dist(_rng)-1); }
  };
  
  /** \class JumpGenerator_Uniform_1d CTRW/JumpGenerator.h "CTRW/JumpGenerator.h"
   * \brief One-dimensional random walk jumps with uniformly-distributed step sizes. */
  template <typename RNG = std::mt19937>
  class JumpGenerator_Uniform_1d
  {
  public:
    const double max_jump_size;  /**< Jumps uniformly distributed between this maximum size and its negative. */
    
    /** Constructor.
     \param max_jump_size Maximum jump size.
     */
    JumpGenerator_Uniform_1d(double max_jump_size)
    : max_jump_size{ max_jump_size }
    {}
    
    /**
     \param state Particle state (unused).
     \return Jump increment.
    */
    template <typename State = useful::Empty>
    double operator() (State const& state = {})
    { return max_jump_size*(2.*_dist(_rng) - 1.); }
    
  private:
    RNG _rng{ std::random_device{}() };              /**< Random number generator.   */
    std::uniform_real_distribution<double> _dist{};  /**< Unit uniform distribution. */
  };
  
  /** \class JumpGenerator_JumpAngle_2d CTRW/JumpGenerator.h "CTRW/JumpGenerator.h"
    \brief Random walk jumps along current orientation with fixed step size.
    \note Particle state must define:
    -  orientation [scalar type]
   */
  class JumpGenerator_JumpAngle_2d
  {
  public:
    const double jump_size;  /**< Fixed jump size. */
    
    /** Constructor.
     \param jump_size Fixed jump size
     */
    JumpGenerator_JumpAngle_2d(double jump_size)
    : jump_size{ jump_size }
    {}
    
    /**
     \param state Particle state.
     \return Orientation increment.
    */
    template <typename State>
    std::vector<double> operator()(State const& state)
    {
      return { jump_size*std::cos(state.orientation),
        jump_size*std::sin(state.orientation) };
    }
  };
  
  /** \class OrientationGenerator_Gradient_1d CTRW/JumpGenerator.h "CTRW/JumpGenerator.h"
    \brief Gaussian-distributed angle jumps around local gradient towards preferred concentration value.
    \note State must define:
    - orientation [scalar type] */
  template
  <typename Concentration, typename Gradient, typename RNG = std::mt19937>
  class OrientationGenerator_Gradient_1d
  {
  public:
    double variance;                 /**< Variance of jumps around local gradient.*/
    double preferred_concentration;  /**< Concentration to jump towards.          */
    
    /** Constructor
     \param concentration Concentration as a function of state,
     \param gradient Concentration gradient as a function of state,
     \param variance Variance of jumps around local gradient.
     \param preferred_concentration Preferred concentration to seek.
    */
    OrientationGenerator_Gradient_1d
    (Concentration&& concentration,
     Gradient&& gradient,
     double variance, double preferred_concentration)
    : variance{ variance }
    , preferred_concentration{ preferred_concentration }
    , _concentration{ std::forward<Concentration>(concentration) }
    , _gradient{ std::forward<Gradient>(gradient) }
    {}
    
    /**
     \param state Particle state.
     \return Orientation increment.
    */
    template <typename State>
    double operator() (State const& state)
    {
      double concentration_val = _concentration(state);
      std::vector<double> gradient_val = _gradient(state);
      
      // If local gradient is zero jump in a random direction
      if (operation::abs(gradient_val) == 0.)
        return constants::pi*(2.*_dist(_rng) - 1.);
      
      // Orient mean jump direction along gradient if below preferred concentration,
      // or along opposite direction of gradient if above
      double reference_angle =
        concentration_val < preferred_concentration
        ? std::atan2(gradient_val[1], gradient_val[0])
        : std::atan2(-gradient_val[1], -gradient_val[0]);
      
      // Compute angle jump and bound it to ]-pi,pi]
      double angle = variance*_dist(_rng) + reference_angle;
      angle = std::abs(angle) > constants::pi
      ? constants::pi
      : angle;
      
      return angle - state.orientation;
    }
    
  private:
    RNG _rng{ std::random_device{}() };                /**< Random number generator. */
    std::normal_distribution<double> _dist{ 0., 1. };  /**< Unit normal distribution.  */
    Concentration _concentration;                      /**< Concentration as a function of state. */
    Gradient _gradient;                                /**< Gradient as a function of state. */
  };
  template
  <typename Concentration, typename Gradient>
  OrientationGenerator_Gradient_1d(Concentration, Gradient, double, double)
  ->OrientationGenerator_Gradient_1d<Concentration, Gradient, std::mt19937>;
    
  /** \class OrientationGenerator_Flip CTRW/JumpGenerator.h "CTRW/JumpGenerator.h"
   * \brief Flip orientation angle. */
  class OrientationGenerator_Flip
  {
  public:
    /**
     \param state Particle state (unused).
     \return Orientation increment.
    */
    template <typename State = useful::Empty>
    double operator() (State const& state = {})
    { return constants::pi; }
  };
}


#endif /* CTRW_JUMPGENERATOR_H */
