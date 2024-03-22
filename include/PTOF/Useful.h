/**
 \file PTOF/Useful.h
 \author Tomás Aquino
 \date 2022/02/17
*/

#ifndef PTOF_USEFUL_H
#define PTOF_USEFUL_H

#include <iostream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <fieldTypes.H>
#include "CTRW/StateGetter.h"
#include "General/Operations.h"
#include "General/Useful.h"

namespace ptof
{
  /** \struct CheckOptions PTOF/Useful.h "PTOF/Useful.h"
   \brief Options for checking conditions. */
  struct CheckOptions
  {
    /** \struct CheckOptions::NoChecl PTOF/Useful.h "PTOF/Useful.h"
      \brief No bounds checking. */
    struct NoCheck{};
    
    /** \struct CheckOptions::Check PTOF/Useful.h "PTOF/Useful.h"
      \brief Bounds checking, no warning. */
    struct Check{};
    
    /** \struct CheckOptions::Warn PTOF/Useful.h "PTOF/Useful.h"
      \brief Bounds checking and warning. */
    struct Warn{};
  };
  
  /** \struct ParallelOptions PTOF/Useful.h "PTOF/Useful.h"
   \brief Options to choose between serial and parallel implementations. */
  struct ParallelOptions
  {
    /** \struct ParallelOptions::Serial PTOF/Useful.h "PTOF/Useful.h"
    \brief Serial implementations. */
    struct Serial{};
    
    /** \struct ParallelOptions::Parallel PTOF/Useful.h "PTOF/Useful.h"
    \brief Parallel implementations. */
    struct Parallel{};
  };
  
  /** \struct SearchOptions PTOF/Useful.h "PTOF/Useful.h"
   \brief Options for mesh searching. */
  struct SearchOptions
  {
    /** \struct SearchOptions::NeighborPrecheck PTOF/Useful.h "PTOF/Useful.h"
    \brief Level of neighbors to check before searching mesh for position.. */
    template <int level>
    struct NeighborPrecheck
    { static constexpr int neighbor_check_level{ level }; };
    
    /** \brief Dot not check neighbors, even self, before searching mesh for position. */
    using NoNeighborPrecheck = NeighborPrecheck<-1>;
    
    /** \brief Check if position is in given cell before searching mesh for position */
    using ZeroNeighborPrecheck = NeighborPrecheck<0>;
    
    /** \brief Check if position is in given cell and in first orthogonal neighbors before searching mesh for position */
    using FirstNeighborPrecheck = NeighborPrecheck<1>;
    
    /** \brief Check if position is in given cell and in first and second orthogonal neighbors before searching mesh for position */
    using SecondNeighborPrecheck = NeighborPrecheck<2>;
  };
  
  /** \struct Dynamics PTOF/Useful.h "PTOF/Useful.h"
     \brief Keep track of names of dynamics types for different types of simulations. */
  struct Dynamics
  {
    /** \enum Dynamics::Type
     *  \brief Implemented dynamics types. */
    enum class Type
    {
      transport,       /**< No parameters, custom boundary does nothing.   */
      firstpassage,    /**< Parameters should hold InitialCondition object.*/
    };
    
    /** \brief Get type from name. */
    static auto type(std::string const& name)
    { return name_to_type.at(name); }
    
    /** \brief Get name from type. */
    static auto name(Type type)
    { return type_to_name.at(type); }
    
    /** \brief Check if name exists. */
    static bool contains(std::string const& name)
    { return name_to_type.count(name); }
    
    /** \brief Map names to types. */
    inline static const
    std::unordered_map<std::string, Type> name_to_type
    {
      { "transport", Type::transport },
      { "firstpassage", Type::firstpassage },
    };
    
    /** \brief Map types to names. */
    inline static const
    std::unordered_map<Type, std::string> type_to_name
    {
      { Type::transport, "transport" },
      { Type::firstpassage, "firstpassage" },
    };
  };
  
  /** \struct Periodicity PTOF/Useful.h "PTOF/Useful.h"
     \brief Keep track of names of periodicity types for boundary conditions. */
  struct Periodicity
  {
    /** \enum Periodicity::Type
     *  \brief Implemented periodicity types. */
    enum class Type
    {
      cartesian,       /**< Cartesian periodicity.                         */
      symmetryplanes   /**< Periodicity according to symmetry planes.      */
    };
    
    /** \brief Get type from name. */
    static auto type(std::string const& name)
    { return name_to_type.at(name); }
    
    /** \brief Get name from type. */
    static auto name(Type type)
    { return type_to_name.at(type); }
    
    /** \brief Check if name exists. */
    static bool contains(std::string const& name)
    { return name_to_type.count(name); }
    
    /** \brief Map names to types. */
    inline static const
    std::unordered_map<std::string, Type> name_to_type
    {
      { "cartesian", Type::cartesian },
      { "symmetryplanes", Type::symmetryplanes }
    };
    
    /** \brief Map types to names. */
    inline static const
    std::unordered_map<Type, std::string> type_to_name
    {
      { Type::cartesian, "cartesian" },
      { Type::symmetryplanes, "symmetryplanes" }
    };
  };
  
  /** \brief Return whether cell index indicates outside mesh
  \tparam warn_if_outside Output warning if true, do not if false
  \param position Position to check.
  \param cell Mesh index of cell position is in.
  \return true if outside, false otherwise.
   \note \p cell is the actual index of the cell position is in, determined elsewhere, not a hint.
  */
  bool outside
  (Foam::label cell)
  {
    if (cell < 0)
      return 1;
    return 0;
  }
  
  /** \brief Return whether cell index indicates outside mesh, and warn and print position if warning is enabled
  \tparam warn_if_outside Output warning if true, do not if false
  \param position Position to check.
  \param cell Mesh index of cell position is in.
  \param extra_warning_info Additional information to output if warning is issued.
  \return true if outside, false otherwise.
   \note \p cell is the actual index of the cell position is in, determined elsewhere, not a hint.
  */
  template <bool warn_if_outside>
  bool outside
  (Foam::label cell, Foam::point const& position, std::string const& extra_warning_info = {})
  {
    if (cell < 0)
    {
      if constexpr (warn_if_outside)
      {
        std::cerr << "Warning: Requested cell side at position "
                  << "("
                  << position[0] << ", "
                  << position[1] << ", "
                  << position[2] << ")"
                  << " outside mesh.";
        if (!extra_warning_info.empty())
          std::cerr << " ";
        std::cerr << extra_warning_info << "\n";
      }
      return 1;
    }
    return 0;
  }
  
  /** \brief Make 3D point from 2D point. */
  auto make_point(Foam::Vector2D<Foam::scalar> const& point)
  {
    return Foam::point{ point[0], point[1], 0. };
  }

  /** \brief Make 3D point from 1D point. */
  auto make_point(Foam::scalar point)
  {
    return Foam::point{ point, 0., 0. };
  }

  /** \brief Make 3D point from 3D point (simple copy). */
  auto make_point(Foam::point const& point)
  { return point; }
  
  /** \brief Make point from vector (simple copy). */
  Foam::point make_point(std::vector<double> const& point)
  {
    if (point.size() >= 3)
      return { point[0], point[1], point[2] };
    if (point.size() == 2)
      return Foam::point{ point[0], point[1], 0. };
    if (point.size() == 1)
      return Foam::point{ point[0], 0., 0. };
    return { 0., 0., 0. };
  }
  
  /** \brief Unit normal at boundary face, pointing outward. */
  template <typename Mesh>
  auto unit_normal_outward
  (Foam::label face, Mesh const& mesh)
  {
    return mesh.faces()[face].unitNormal(mesh.points());
  }
  
  /** \brief Unit normal at boundary face, pointing inward. */
  template <typename Mesh>
  auto unit_normal_inward
  (Foam::label face, Mesh const& mesh)
  {
    return -mesh.faces()[face].unitNormal(mesh.points());
  }
  
  /** \brief Area (magnitude) of a face. */
  template <typename Mesh>
  auto face_area
  (Foam::label face, Mesh const& mesh)
  {
    return mesh.faces()[face].mag(mesh.points());
  }
  
  /** \brief Center of a face. */
  template <typename Mesh>
  auto face_center
  (Foam::label face, Mesh const& mesh)
  {
    return mesh.faces()[face].centre(mesh.points());
  }
  
  /** \brief Volume of a cell. */
  template <typename Mesh>
  auto cell_volume
  (Foam::label cell, Mesh const& mesh)
  {
    return mesh.cellVolumes()[cell];
  }
  
  /** \brief Center of a cell. */
  template <typename Mesh>
  auto cell_center
  (Foam::label cell, Mesh const& mesh)
  {
    return mesh.cellCentres()[cell];
  }
  
  /** Sometimes the face center associated with a mesh face is not considered within the cell. Verify this */
  template <typename Mesh, typename Locator>
  bool face_center_is_in_cell
  (Foam::label face, Mesh const& mesh, Locator const& locator)
  {
    return locator.mesh_search().findCell(face_center(face, locator.mesh_search().mesh()))
      == mesh.owner()[face];
  }
  
  /** \brief Small offset forward given current face and direction.*/
  template <typename Locator>
  Foam::vector offset_face
  (Foam::point const& begin,
  Foam::label face,
  Foam::vector const& direction,
  Locator const& locator)
  {
    Foam::label owner_cell = locator.mesh_search().mesh().faceOwner()[face];
    Foam::point const& cell_center = locator.mesh_search().mesh().
      cellCentres()[owner_cell];
    Foam::scalar typ_dim = Foam::mag(cell_center
                                     - face_center(face,
                                                   locator.mesh_search().mesh()));
   
    return locator.mesh_search().tol_*typ_dim*
      direction/Foam::mag(direction);
  }
  
  /**
   \brief Small offset along face normal (outward).
   \param face Mesh face index to use face center as beginning point.
   \param locator Object to locate positions in mesh.
   \return Offset vector.
  */
  template <typename Locator>
  Foam::vector offset_face
  (Foam::label face,
   Locator const& locator)
  {
    return offset_face(face_center(face, locator.mesh_search().mesh()),
                       face,
                       unit_normal_outward(face, locator.mesh_search().mesh()),
                       locator);
  }
  
  /**
   \brief Small offset forward from a point on a face.
   \param begin Starting point (on a face).
   \param face Mesh face index where \p begin is located.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
  */
  template <typename Locator>
  Foam::point offset_forward_face
  (Foam::point const& begin,
   Foam::label face,
   Foam::vector const& direction,
   Locator const& locator)
  {
    return begin + offset_face(begin, face, direction, locator);;
  }
  
  /**
   \brief Small offset forward from a point on a face.
   \details If begin point is in mesh, guarantee offset point in mesh.
   \param begin Starting point (on a face).
   \param face Mesh face index where \p begin is located.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
  */
  template <typename Locator>
  Foam::point offset_forward_face_keep_inside
  (Foam::point const& begin,
   Foam::label face,
   Foam::vector const& direction,
   Locator const& locator)
  {
    auto offset = offset_face(begin, face, direction, locator);
    auto point = begin + offset;
    Foam::label owner_cell = locator.mesh_search().mesh().faceOwner()[face];
    if (locator(begin, owner_cell) < 0)
      return point;
    while (locator(point, owner_cell) < 0)
    {
      offset /= 2.;
      point = begin + offset;
    }
    return point;
  }
  
  /**
   \brief Small offset backward from a point on a face.
   \param begin Starting point (on a face).
   \param face Mesh face index where \p begin is located.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
  */
  template <typename Locator>
  Foam::point offset_backward_face
  (Foam::point const& begin,
   Foam::label face,
   Foam::vector const& direction,
   Locator const& locator)
  {
    return begin - offset_face(begin, face, direction, locator);
  }
  
  /**
   \brief Small offset backward from a point on a face.
   \details If begin point is in mesh, guarantee offset point in mesh.
   \param begin Starting point (on a face).
   \param face Mesh face index where \p begin is located.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
  */
  template <typename Locator>
  Foam::point offset_backward_face_keep_inside
  (Foam::point const& begin,
   Foam::label face,
   Foam::vector const& direction,
   Locator const& locator)
  {
    auto offset = offset_face(begin, face, direction, locator);
    auto point = begin - offset;
    Foam::label owner_cell = locator.mesh_search().mesh().faceOwner()[face];
    if (locator(begin, owner_cell) < 0)
      return point;
    while (locator(point, owner_cell) < 0)
    {
      offset /= 2.;
      point = begin - offset;
    }
    return point;
  }
  
  /**
   \brief Small offset outward along face normal starting from face center.
   \details If begin point is in mesh, guarantee offset point in mesh.
   \param face Mesh face index.
   \param locator Object to locate positions in mesh.
   \return Offset point.
  */
  template <typename Locator>
  Foam::point offset_outward_face
  (Foam::label face,
   Locator const& locator)
  {
    return face_center(face, locator.mesh_search().mesh())
      + offset_face(face, locator);
  }
  
   /**
    \brief Small offset inward along face normal starting from face center.
    \details If begin point is in mesh, guarantee offset point in mesh.
    \param face Mesh face index.
    \param locator Object to locate positions in mesh.
    \return Offset point.
   */
   template <typename Locator>
   Foam::point offset_inward_face
   (Foam::label face,
    Locator const& locator)
   {
     return face_center(face, locator.mesh_search().mesh())
      - offset_face(face, locator);
   }
  
  /**
   \brief Small offset outward along face normal starting from face center.
   \details If begin point is in mesh, guarantee offset point in mesh.
   \param begin Starting point in cell.
   \param cell Mesh cell index.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
  */
  template <typename Locator>
  Foam::vector offset_cell
  (Foam::point const& begin,
   Foam::label cell,
   Foam::vector const& direction,
   Locator const& locator)
  {
    Foam::point const& center = locator.mesh_search().mesh().cellCentres()[cell];
    Foam::scalar typ_dim = Foam::mag(center - begin);
    
    return locator.mesh_search().tol_*typ_dim*
      direction/Foam::mag(direction);
  }
  
  /**
   \brief Small offset forward.
   \param begin Starting point in cell.
   \param cell Hint for cell index where \p begin is located.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
  */
  template <typename Locator>
  Foam::point offset_forward_cell
  (Foam::point const& begin,
   Foam::label cell,
   Foam::vector const& direction,
   Locator const& locator)
  {
    return begin
      + offset_cell(begin, cell, direction, locator);
  }
  
  /** Small offset forward
   \brief Small offset forward.
   \details If begin point is in mesh, guarantee offset point in mesh.
   \param begin Starting point in cell.
   \param cell Hint for cell index where \p begin is located.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
  */
  template <typename Locator>
  Foam::point offset_forward_cell_keep_inside
  (Foam::point const& begin,
   Foam::label cell,
   Foam::vector const& direction,
   Locator const& locator)
  {
    auto offset = offset_cell(begin, cell, direction, locator);
    auto point = begin + offset;
    if (locator(begin, cell) < 0)
      return point;
    while (locator(point, cell) < 0)
    {
      offset /= 2.;
      point = begin + offset;
    }
    return point;
  }
  
  /**
   \brief Small offset backward.
   \param begin Starting point in cell.
   \param cell Hint for cell index where \p begin is located.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
  */
  template <typename Locator>
  Foam::point offset_backward_cell
  (Foam::point const& begin,
   Foam::label cell,
   Foam::vector const& direction,
   Locator const& locator)
  {
    return begin
      - offset_cell(begin, cell, direction, locator);
  }
  
  /**
   \brief Small offset backward.
   \details If begin point is in mesh, guarantee offset point in mesh.
   \param begin Starting point in cell.
   \param cell Hint for cell index where \p begin is located.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
  */
  template <typename Locator>
  Foam::point offset_backward_cell_keep_inside
  (Foam::point const& begin,
   Foam::label cell,
   Foam::vector const& direction,
   Locator const& locator)
  {
    auto offset = offset_cell(begin, cell, direction, locator);
    auto point = begin - offset;
    if (locator(begin, cell) < 0)
      return point;
    while (locator(point, cell) < 0)
    {
      offset /= 2.;
      point = begin - offset;
    }
    return point;
  }
  
  /** \class Has_periodicity
   \brief Check for existence of periodicity info (to use with State objects).
  */
  template <typename State, typename = int>
  struct Has_periodicity : std::false_type { };
  template <typename State>
  struct Has_periodicity
  <State, decltype((void) State::periodicity, 0)>
  : std::true_type { };
  
  /**
   \param subject CTRW object.
   \param time Current time.
   \return Number of absorbed particles.
   \note Particle states must implement:
   - info.asorbed [std::size_t]
   - time */
  template <typename Subject>
  std::size_t nr_absorbed
  (Subject const& subject, double time)
  {
    std::size_t absorbed = 0;
    for (auto const& part : subject.particles())
      
      if (part.state_new().info.absorbed
          && part.state_old().time <= time)
        ++absorbed;
    return absorbed;
  };
  
  /**
  \param subject CTRW object.
  \param time Current time.
  \return Total mass.
  \note Particle states must implement:
   - mass [double]
   - info.asorbed [std::size_t]
   - time */
  template <typename Subject>
  auto mass
  (Subject const& subject, double time)
  {
    double mass = 0.;
    for (auto const& part : subject.particles())
      if (!part.state_new().info.absorbed
          && part.state_old().time <= time)
      {
        mass += part.state_new().mass;
      }
    return mass;
  };
  
  /**
   \param subject CTRW object.
   \param time Current time.
   \param locator Object to locate positions in mesh.
   \param masks Container of mask reference wrappers. Masks are scalar fields assigned values to mesh cells through operator[].
   \param tolerances Vector of tolerances for each mask, such that cells where a mask is above the tolerance are considered.
   \return Total masses in regions specified by masks.
   \note -Particle states must implement:
   -# mass [double]
   -# info.asorbed [std::size_t]
   -# time
   
   = \c tolerances must have at least the same size as \c masks  */
  template <typename Subject, typename Locator, typename Mask>
  auto mass
  (Subject const& subject, double time,
   Locator const& locator,
   std::vector<std::reference_wrapper<const Mask>> masks,
   std::vector<double> tolerances)
  {
    std::vector<double> masses(masks.size(), 0.);
    for (auto const& part : subject.particles())
    {
      if (!part.state_new().info.absorbed
          && part.state_old().time <= time)
      {
        auto cell = locator(part.state_new());
        for (std::size_t ii = 0; ii < masks.size(); ++ii)
        {
          if (cell != -1 && masks[ii].get()[cell] > tolerances[ii])
            masses[ii] += part.state_new().mass;
        }
      }
    }
    return masses;
  };
  
  /**
    \param subject CTRW object.
    \param time Current time.
    \param getter_position Get position from state, gets position directly by default.
    \return Mean position (weighted by mass).
    \note Particle states must implement:
    - position [for default positition getter]
    - info.asorbed [std::size_t]
    - mass
    - time*/
  template
  <typename Subject, typename GetterPosition = ctrw::Get_position>
  auto position_mean
  (Subject const& subject, double time,
   GetterPosition getter_position = {})
  {
    decltype(getter_position(subject.particles(0).state_new()))
      position_mean = Foam::zero{};
    for (auto const& part : subject.particles())
    {
      auto const& state = part.state_new();
      if (!state.info.absorbed
          && part.state_old().time <= time)
      {
        position_mean += state.mass*getter_position(state);
      }
    }
    return position_mean/mass(subject, time);
  };
  
  /**
   \param subject CTRW object.
   \param time Current time.
   \param getter_position Get position from state, gets position directly by default.
   \return Second moment of position (weighted by mass).
   \note Particle states must implement:
   - position [for default positition getter]
   - info.asorbed [std::size_t]
   - mass */
  template
  <typename Subject, typename GetterPosition = ctrw::Get_position>
  auto position_second_moment
  (Subject const& subject, double time,
   GetterPosition getter_position = {})
  {
    decltype(getter_position(subject.particles(0).state_new()))
      second_moment = Foam::zero{};
    for (auto const& part : subject.particles())
    {
      auto const& state = part.state_new();
      if (!state.info.absorbed
          && part.state_old().time <= time)
      {
        auto position = getter_position(state);
        for (std::size_t dd = 0; dd < second_moment.size(); ++dd)
          second_moment[dd] += position[dd]*position[dd];
        second_moment *= state.mass;
      }
    }
    return second_moment/mass(subject, time);
  };
  
  /**
  \param subject CTRW object.
  \param time Current time.
  \param getter_position Get position from state, gets position directly by default.
  \return Position variance (weighted by mass).
  \note Particle states must implement:
  - position [for default positition getter]
  - info.asorbed [std::size_t]
  - mass */
  template
  <typename Subject, typename GetterPosition = ctrw::Get_position>
  auto position_variance
  (Subject const& subject, double time,
   GetterPosition getter_position = {})
  {
    auto position_mean_sq = position_mean(subject, time, getter_position);
    operation::times_InPlace(position_mean_sq, position_mean_sq);
    return position_second_moment(subject, time, getter_position)
      - position_mean_sq;
  };
  
  /** \brief Output time information.
   \param output Output stream.
   \param params Output parameters, holding time_unit_factor to rescale \p time.
   \param time Current time.
   */
  template <typename OStream, typename ParametersOutput>
  void info_time
  (OStream& output, ParametersOutput const& params, double time)
  {
    output << "Time "
           << "[" << params.time_units << " time units]: "
           << time/params.time_unit_factor
           << "\n";
  }
  
  /** \brief Output information about fraction of particles that have not been absorbed.
   \param output Output stream.
   \param subject CTRW object.
   \param time Current time.
   \note Particle states must implement:
   - info.asorbed [std::size_t]
  */
  template <typename OStream, typename Subject>
  void info_fraction_not_absorbed
  (OStream& output, Subject const& subject, double time)
  {
    output << "Fraction not absorbed: "
          << 1. - double(ptof::nr_absorbed(subject, time))/subject.size()
          << "\n";
  }
  
  /** \return Identifier string for output file names. */
  std::string identifier
  (std::string const& model_name,
   std::string const& case_name,
   std::string const& of_case_name,
   std::string const& params_transport_name,
   std::string const& params_reaction_name,
   std::string const& params_solvers_name,
   std::string const& params_initial_condition_name,
   std::string const& params_output_name)
  {
    return "M_" + model_name
      + "_C_" + case_name
      + "_OF_" + of_case_name
      + "_T_" + params_transport_name
      + "_R_" + params_reaction_name
      + "_S_" + params_solvers_name
      + "_I_" + params_initial_condition_name
      + "_O_" + params_output_name;
  }
  
  /** \return Identifier string for output file names. */
  std::string identifier
  (std::string const& model_name,
   std::string const& case_name,
   std::string const& of_case_name,
   std::string const& params_transport_name,
   std::string const& params_reaction_name,
   std::string const& params_solvers_name,
   std::string const& params_initial_condition_name,
   std::string const& params_output_name,
   std::size_t run_nr)
  {
    return identifier(model_name,
                      case_name,
                      of_case_name,
                      params_transport_name,
                      params_reaction_name,
                      params_solvers_name,
                      params_initial_condition_name,
                      params_output_name)
      + "_RUN_" + std::to_string(run_nr);
  }
  
  /** \return Identifier string for output file names. */
  std::string identifier
  (std::string const& model_name,
   std::string const& case_name,
   std::string const& of_case_name,
   std::string const& params_transport_name,
   std::string const& params_phase_name,
   std::string const& params_reaction_name,
   std::string const& params_solvers_name,
   std::string const& params_initial_condition_name,
   std::string const& params_output_name)
  {
    return "M_" + model_name
      + "_C_" + case_name
      + "_OF_" + of_case_name
      + "_T_" + params_transport_name
      + "_P_" + params_phase_name
      + "_R_" + params_reaction_name
      + "_S_" + params_solvers_name
      + "_I_" + params_initial_condition_name
      + "_O_" + params_output_name;
  }
  
  /** \return Identifier string for output file names. */
  std::string identifier
  (std::string const& model_name,
   std::string const& case_name,
   std::string const& of_case_name,
   std::string const& params_transport_name,
   std::string const& params_phase_name,
   std::string const& params_reaction_name,
   std::string const& params_solvers_name,
   std::string const& params_initial_condition_name,
   std::string const& params_output_name,
   std::size_t run_nr)
  {
    return identifier(model_name,
                      case_name,
                      of_case_name,
                      params_transport_name,
                      params_phase_name,
                      params_reaction_name,
                      params_solvers_name,
                      params_initial_condition_name,
                      params_output_name)
      + "_RUN_" + std::to_string(run_nr);
  }
}

#endif /* PTOF_USEFUL_H */
