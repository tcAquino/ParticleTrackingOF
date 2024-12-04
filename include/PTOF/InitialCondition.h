/**
   \file PTOF/InitialCondition.h
   \author Tomás Aquino
   \date 29/09/2024
   \brief Objects to create particles according to initial condition.
*/

#ifndef PTOF_INITIALCONDITION_H
#define PTOF_INITIALCONDITION_H

#include "CTRW/Meta.h"
#include "General/IO.h"
#include "PTOF/Useful.h"
#include "Stochastic/Random.h"
#include "fvMesh.H"
#include <cmath>
#include <cstddef>
#include <fieldTypes.H>
#include <fstream>
#include <map>
#include <point.H>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace ptof {

/**
   \class InitialCondition PTOF/InitialCondition.h "CTRW/InitialCondition.h"
   \brief Make particles and positions according to prescribed rules.
   \tparam ParticleMaker Abstract object to make a particle given a position.
   \tparam Geometry Object handling domain geometry information.
*/
template <typename ParticleMaker, typename Geometry> struct InitialCondition {
  using Particle = typename std::remove_reference_t<ParticleMaker>::Particle;
  using Position = typename Particle::State::Position;
  struct PositionAndHint {
    Position position;
    Foam::label hint{-1};
  };
  using ParticleContainer = std::vector<Particle>;
  using PositionContainer = std::vector<Position>;

  InitialCondition(ParticleMaker &&particle_maker, Geometry &&geometry)
      : _particle_maker{std::forward<ParticleMaker>(particle_maker)},
        _geometry{std::forward<Geometry>(geometry)} {}

  virtual ~InitialCondition() = default;

  /**
     \brief Make particles.
     \param nr_particles Number of particles to make.
     \return Container with particles.
  */
  ParticleContainer operator()(std::size_t nr_particles) {
    return make_particles(nr_particles);
  }

  /**
     \brief Make a single particle.
     \return Particle.
  */
  virtual Particle make_particle() {
    auto [position, hint] = make_position_and_hint();
    return _particle_maker(position, hint);
  }

  /**
     \brief Make a single position.
     \return Position.
  */
  virtual Position make_position() { return make_position_and_hint().position; }

  /**
     \brief Make a single position along with location hint.
     \param nr_positions Number of positions to make
     \return Position.
  */
  virtual PositionAndHint make_position_and_hint() = 0;

  virtual ParticleContainer make_particles(std::size_t nr_particles) {
    ParticleContainer particles;
    particles.reserve(nr_particles);
    for (std::size_t pp = 0; pp < nr_particles; ++pp)
      particles.emplace_back(make_particle());

    return particles;
  };

  virtual PositionContainer make_positions(std::size_t nr_positions) {
    PositionContainer positions;
    positions.reserve(nr_positions);
    for (std::size_t pp = 0; pp < nr_positions; ++pp)
      positions.emplace_back(make_position());

    return positions;
  };

protected:
  ParticleMaker _particle_maker;
  Geometry _geometry;
};
template <typename ParticleMaker, typename Geometry>
InitialCondition(ParticleMaker &&, Geometry &&)
    -> InitialCondition<ParticleMaker, Geometry>;

/**
   \class InitialCondition_Continuous PTOF/InitialCondition.h
   "CTRW/InitialCondition.h" \brief Continuous injection from a series of
   pulses. \tparam InitialCondition InitalCondition object to handle each pulse.
*/
template <typename InitialCondition>
struct InitialCondition_Continuous : public InitialCondition {
private:
  using IC = InitialCondition;

public:
  InitialCondition_Continuous(InitialCondition &&initial_condition,
                              double injection_time, double injection_time_end,
                              double injection_time_step)
      : IC{std::forward<IC>(initial_condition)},
        _injection_time{injection_time},
        _injection_time_end{injection_time_end}, _injection_time_step{
                                                     injection_time_step} {}

  virtual typename IC::ParticleContainer
  make_particles(std::size_t nr_particles) override {
    typename IC::ParticleContainer particles;
    std::size_t nr_injections = std::ceil(
        (_injection_time_end - _injection_time) / _injection_time_step);
    particles.reserve(nr_injections * nr_particles);

    for (double time = _injection_time; time < _injection_time_end;
         time += _injection_time_step) {
      this->_particle_maker.time = time;
      particles.push_back(IC::make_particle());
    }

    return particles;
  };

  typename IC::PositionAndHint make_position_and_hint() override {
    throw std::runtime_error{"Making positions is not allowed for continuous "
                             "injections"};
  }

  virtual typename IC::Particle make_particle() override {
    throw std::runtime_error{
        "Making single particles is not allowed for continuous "
        "injections"};
  }

  virtual typename IC::Position make_position() override {
    throw std::runtime_error{"Making positions is not allowed for continuous "
                             "injections"};
  }

protected:
  double _injection_time;
  double _injection_time_end;
  double _injection_time_step;
};
template <typename InitialCondition>
InitialCondition_Continuous(InitialCondition &&, double, double, double)
    -> InitialCondition_Continuous<InitialCondition>;

/**
   \class InitialCondition_Point PTOF/InitialCondition.h
   "CTRW/InitialCondition.h"
   \brief Initial condition at a point.
*/
template <typename ParticleMaker, typename Geometry>
struct InitialCondition_Point
    : public InitialCondition<ParticleMaker, Geometry> {
private:
  using IC = InitialCondition<ParticleMaker, Geometry>;

public:
  InitialCondition_Point(ParticleMaker &&particle_maker, Geometry &&geometry,
                         typename IC::Position position)
      : IC{std::forward<ParticleMaker>(particle_maker),
           std::forward<Geometry>(geometry)},
        _position{position} {
    if (outside(_cell_id))
      throw std::runtime_error{
          "Point initial condition position is outside mesh"};
  }

  typename IC::PositionAndHint make_position_and_hint() override {
    return {_position, _cell_id};
  }

protected:
  stochastic::RNGThreaded<
      typename std::remove_reference_t<Geometry>::ParallelOption, std::mt19937>
      _rng;
  typename IC::Position _position;
  Foam::label _cell_id{this->_geometry.locator(_position)};
};
template <typename ParticleMaker, typename Geometry>
InitialCondition_Point(
    ParticleMaker &&, Geometry &&,
    typename InitialCondition<ParticleMaker, Geometry>::Position)
    -> InitialCondition_Point<ParticleMaker, Geometry>;

/**
   \class InitialCondition_DistributedPosition PTOF/InitialCondition.h
   "CTRW/InitialCondition.h"
   \brief Initial condition with randomly distributed positions.
   \tparam Distribution Object that takes an RNG and produces a (random)
   position.
*/
template <typename ParticleMaker, typename Geometry, typename Distribution>
struct InitialCondition_DistributedPosition
    : public InitialCondition<ParticleMaker, Geometry> {
private:
  using IC = InitialCondition<ParticleMaker, Geometry>;

public:
  InitialCondition_DistributedPosition(ParticleMaker &&particle_maker,
                                       Geometry &&geometry, Distribution &&dist,
                                       std::size_t nr_tries)
      : IC{std::forward<ParticleMaker>(particle_maker),
           std::forward<Geometry>(geometry)},
        _dist{std::forward<Distribution>(dist)}, _nr_tries{nr_tries} {}

  /** \brief Make a single position according to prescribed initial condition.
  \param nr_positions Number of positions to make
  \return Position.  */
  typename IC::PositionAndHint make_position_and_hint() override {
    for (std::size_t ii = 0; ii < _nr_tries; ++ii) {
      auto position = IC::Particle::State::make_position(_dist(_rng));
      auto cell_id = this->_geometry.locator(position);
      if (!ptof::outside(cell_id))
        return {position, cell_id};
    }

    throw std::runtime_error{"Failed to make position inside mesh after " +
                             std::to_string(_nr_tries) + " tries"};
  }

protected:
  stochastic::RNGThreaded<
      typename std::remove_reference_t<Geometry>::ParallelOption, std::mt19937>
      _rng;
  Distribution _dist;
  std::size_t _nr_tries;
};
template <typename ParticleMaker, typename Geometry, typename Distribution>
InitialCondition_DistributedPosition(ParticleMaker &&, Geometry &&,
                                     Distribution &&, std::size_t)
    -> InitialCondition_DistributedPosition<ParticleMaker, Geometry,
                                            Distribution>;

/**
   \class InitialCondition_DistributedCellCenters PTOF/InitialCondition.h
   "CTRW/InitialCondition.h"
   \brief Initial condition at randomly distributed cell centers.
   \tparam CellIdContainer Container of cell indices.
   \tparam Distribution Object that takes an RNG and produces a (random)
   cell index.
*/
template <typename ParticleMaker, typename Geometry, typename CellIdContainer,
          typename Distribution>
struct InitialCondition_DistributedCellCenters
    : public InitialCondition<ParticleMaker, Geometry> {
private:
  using IC = InitialCondition<ParticleMaker, Geometry>;

public:
  InitialCondition_DistributedCellCenters(ParticleMaker &&particle_maker,
                                          Geometry &&geometry,
                                          CellIdContainer &&cell_ids,
                                          Distribution &&dist)
      : IC{std::forward<ParticleMaker>(particle_maker),
           std::forward<Geometry>(geometry)},
        _cell_ids{std::forward<CellIdContainer>(cell_ids)},
        _dist{std::forward<Distribution>(dist)} {
    if (_cell_ids.size() == 0)
      throw std::runtime_error{"No cells provided for initial condition"};
  }

  typename IC::PositionAndHint make_position_and_hint() override {
    auto cell_id = _cell_ids[_dist(_rng)];
    return {IC::Particle::State::make_position(
                cell_center(cell_id, this->_geometry.mesh())),
            cell_id};
  }

protected:
  stochastic::RNGThreaded<
      typename std::remove_reference_t<Geometry>::ParallelOption, std::mt19937>
      _rng;
  CellIdContainer _cell_ids;
  Distribution _dist;
};
template <typename ParticleMaker, typename Geometry, typename CellIdContainer,
          typename Distribution>
InitialCondition_DistributedCellCenters(ParticleMaker &&, Geometry &&,
                                        CellIdContainer &&, Distribution &&)
    -> InitialCondition_DistributedCellCenters<ParticleMaker, Geometry,
                                               CellIdContainer, Distribution>;

/**
   \class InitialCondition_DistributedFaceCenters PTOF/InitialCondition.h
   "CTRW/InitialCondition.h"
   \brief Initial condition at randomly distributed face centers.
   \tparam CellIdContainer Container of face indices.
   \tparam Distribution Object that takes an RNG and produces a (random)
   face index.
*/
template <typename ParticleMaker, typename Geometry, typename FaceIdContainer,
          typename Distribution>
struct InitialCondition_DistributedFaceCenters
    : public InitialCondition<ParticleMaker, Geometry> {
private:
  using IC = InitialCondition<ParticleMaker, Geometry>;

public:
  InitialCondition_DistributedFaceCenters(ParticleMaker &&particle_maker,
                                          Geometry &&geometry,
                                          FaceIdContainer &&face_ids,
                                          Distribution &&dist)
      : IC{std::forward<ParticleMaker>(particle_maker),
           std::forward<Geometry>(geometry)},
        _face_ids{std::forward<FaceIdContainer>(face_ids)},
        _dist{std::forward<Distribution>(dist)} {
    if (_face_ids.size() == 0)
      throw std::runtime_error{"No faces provided for initial condition"};
  }

  /** \note If face center is not in owner cell, position is placed at cell
   * center.*/
  typename IC::PositionAndHint make_position_and_hint() override {
    auto const &mesh = this->_geometry.mesh();
    auto face_id = _face_ids[_dist(_rng)];
    return {IC::Particle::State::make_position(
                adjusted_face_center(face_id, this->_geometry.locator)),
            mesh.faceOwner()[face_id]};
  }

protected:
  stochastic::RNGThreaded<
      typename std::remove_reference_t<Geometry>::ParallelOption, std::mt19937>
      _rng;
  FaceIdContainer _face_ids;
  Distribution _dist;
};
template <typename ParticleMaker, typename Geometry, typename FaceIdContainer,
          typename Distribution>
InitialCondition_DistributedFaceCenters(ParticleMaker &&, Geometry &&,
                                        FaceIdContainer &&, Distribution &&)
    -> InitialCondition_DistributedFaceCenters<ParticleMaker, Geometry,
                                               FaceIdContainer, Distribution>;

/**
   \class InitialCondition_UniformCellCenters PTOF/InitialCondition.h
   "CTRW/InitialCondition.h"
   \brief Initial condition at uniformly randomly distributed cell centers
   (weighted by cell volume).
   \tparam CellIdContainer Container of cell indices.
*/
template <typename ParticleMaker, typename Geometry, typename CellIdContainer>
struct InitialCondition_UniformCellCenters
    : public InitialCondition_DistributedCellCenters<
          ParticleMaker, Geometry, CellIdContainer,
          std::discrete_distribution<std::size_t>> {
private:
  using IC = InitialCondition_DistributedCellCenters<
      ParticleMaker, Geometry, CellIdContainer,
      std::discrete_distribution<std::size_t>>;

public:
  InitialCondition_UniformCellCenters(ParticleMaker &&particle_maker,
                                      Geometry &&geometry,
                                      CellIdContainer &&cell_ids)
      : IC{std::forward<ParticleMaker>(particle_maker),
           std::forward<Geometry>(geometry),
           std::forward<CellIdContainer>(cell_ids),
           uniform_cell_distribution(cell_ids, geometry.mesh())} {}
};
template <typename ParticleMaker, typename Geometry, typename CellIdContainer>
InitialCondition_UniformCellCenters(ParticleMaker &&, Geometry &&,
                                    CellIdContainer &&)
    -> InitialCondition_UniformCellCenters<ParticleMaker, Geometry,
                                           CellIdContainer>;

/**
   \class InitialCondition_UniformFaceCenters PTOF/InitialCondition.h
   "CTRW/InitialCondition.h"
   \brief Initial condition at uniformly randomly distributed face centers
   (weighted by face area).
   \tparam FaceIdContainer Container of face indices.
*/
template <typename ParticleMaker, typename Geometry, typename FaceIdContainer>
struct InitialCondition_UniformFaceCenters
    : public InitialCondition_DistributedFaceCenters<
          ParticleMaker, Geometry, FaceIdContainer,
          std::discrete_distribution<std::size_t>> {
private:
  using IC = InitialCondition_DistributedFaceCenters<
      ParticleMaker, Geometry, FaceIdContainer,
      std::discrete_distribution<std::size_t>>;

public:
  InitialCondition_UniformFaceCenters(ParticleMaker &&particle_maker,
                                      Geometry &&geometry,
                                      FaceIdContainer &&face_ids)
      : IC{std::forward<ParticleMaker>(particle_maker),
           std::forward<Geometry>(geometry),
           std::forward<FaceIdContainer>(face_ids),
           uniform_face_distribution(face_ids, geometry.mesh())} {}
};
template <typename ParticleMaker, typename Geometry, typename CellIdContainer>
InitialCondition_UniformFaceCenters(ParticleMaker &&, Geometry &&,
                                    CellIdContainer &&)
    -> InitialCondition_UniformFaceCenters<ParticleMaker, Geometry,
                                           CellIdContainer>;

/**
   \class InitialCondition_FluxweightedCellCenters PTOF/InitialCondition.h
   "CTRW/InitialCondition.h"
   \brief Initial condition at flux-weighted randomly distributed cell centers.
   \tparam CellIdContainer Container of cell indices.
*/
template <typename ParticleMaker, typename Geometry, typename CellIdContainer>
struct InitialCondition_FluxweightedCellCenters
    : public InitialCondition_DistributedCellCenters<
          ParticleMaker, Geometry, CellIdContainer,
          std::discrete_distribution<std::size_t>> {
private:
  using IC = InitialCondition_DistributedCellCenters<
      ParticleMaker, Geometry, CellIdContainer,
      std::discrete_distribution<std::size_t>>;

public:
  template <typename VelocityField>
  InitialCondition_FluxweightedCellCenters(ParticleMaker &&particle_maker,
                                           Geometry &&geometry,
                                           CellIdContainer &&cell_ids,
                                           VelocityField const &velocity_field)
      : IC{std::forward<ParticleMaker>(particle_maker),
           std::forward<Geometry>(geometry),
           std::forward<CellIdContainer>(cell_ids),
           fluxweighted_cell_distribution(cell_ids, velocity_field,
                                          geometry.mesh())} {}
};
template <typename ParticleMaker, typename Geometry, typename CellIdContainer,
          typename VelocityField>
InitialCondition_FluxweightedCellCenters(ParticleMaker &&, Geometry &&,
                                         CellIdContainer &&,
                                         VelocityField const &)
    -> InitialCondition_FluxweightedCellCenters<ParticleMaker, Geometry,
                                                CellIdContainer>;

/**
   \class InitialCondition_FluxweightedCellCenters PTOF/InitialCondition.h
   "CTRW/InitialCondition.h"
   \brief Initial condition at flux-weighted randomly distributed face centers.
   \tparam FaceIdContainer Container of face indices.
*/
template <typename ParticleMaker, typename Geometry, typename FaceIdContainer>
struct InitialCondition_FluxweightedFaceCenters
    : public InitialCondition_DistributedFaceCenters<
          ParticleMaker, Geometry, FaceIdContainer,
          std::discrete_distribution<std::size_t>> {
private:
  using IC = InitialCondition_DistributedFaceCenters<
      ParticleMaker, Geometry, FaceIdContainer,
      std::discrete_distribution<std::size_t>>;

public:
  template <typename VelocityField>
  InitialCondition_FluxweightedFaceCenters(ParticleMaker &&particle_maker,
                                           Geometry &&geometry,
                                           FaceIdContainer &&face_ids,
                                           VelocityField const &velocity_field)
      : IC{std::forward<ParticleMaker>(particle_maker),
           std::forward<Geometry>(geometry),
           std::forward<FaceIdContainer>(face_ids),
           fluxweighted_face_distribution(face_ids, velocity_field,
                                          geometry.locator)} {}
};
template <typename ParticleMaker, typename Geometry, typename FaceIdContainer,
          typename VelocityField>
InitialCondition_FluxweightedFaceCenters(ParticleMaker &&, Geometry &&,
                                         FaceIdContainer &&,
                                         VelocityField const &)
    -> InitialCondition_FluxweightedFaceCenters<ParticleMaker, Geometry,
                                                FaceIdContainer>;

/**
   \class InitialCondition_PrescribedPositions PTOF/InitialCondition.h
   "CTRW/InitialCondition.h"
   \brief Initial condition at prescribed positions read from a file.
   \details File should contain particle positions, one particle per line.
*/
template <typename ParticleMaker, typename Geometry>
struct InitialCondition_PrescribedPositions
    : public InitialCondition<ParticleMaker, Geometry> {
private:
  using IC = InitialCondition<ParticleMaker, Geometry>;

public:
  InitialCondition_PrescribedPositions(ParticleMaker &&particle_maker,
                                       Geometry &&geometry,
                                       std::string filename)
      : IC{std::forward<ParticleMaker>(particle_maker),
           std::forward<Geometry>(geometry)},
        filename{filename} {}

  typename IC::PositionAndHint make_position_and_hint() override {
    auto split_line = io::split_line(input);
    if (split_line.size() < std::remove_reference_t<Geometry>::dim)
      throw std::runtime_error{"Could not parse particle position in file " +
                               filename};
    Foam::point position;
    for (std::size_t dd = 0; dd < std::remove_reference_t<Geometry>::dim; ++dd)
      position[dd] = std::stod(split_line[dd]);

    return {IC::Particle::State::make_position(position), -1};
  }

protected:
  std::string filename;
  std::ifstream input{io::open_read(filename)};
};
template <typename ParticleMaker, typename Geometry>
InitialCondition_PrescribedPositions(ParticleMaker &&, Geometry &&, std::string)
    -> InitialCondition_PrescribedPositions<ParticleMaker, Geometry>;

/**
   \class InitialCondition_PrescribedPositionsMasses PTOF/InitialCondition.h
   "CTRW/InitialCondition.h"
   \brief Initial condition at prescribed positions and masses read from a file.
   \details File should contain particle positions and masses, one particle per
   line.
*/
template <typename ParticleMaker, typename Geometry>
struct InitialCondition_PrescribedPositionsMasses
    : public InitialCondition<ParticleMaker, Geometry> {
private:
  using IC = InitialCondition<ParticleMaker, Geometry>;

public:
  InitialCondition_PrescribedPositionsMasses(ParticleMaker &&particle_maker,
                                             Geometry &&geometry,
                                             std::string filename)
      : IC{std::forward<ParticleMaker>(particle_maker),
           std::forward<Geometry>(geometry)},
        filename{filename} {}

  typename IC::PositionAndHint make_position_and_hint() override {
    auto split_line = io::split_line(input);
    if (split_line.size() < std::remove_reference_t<Geometry>::dim + 1)
      throw std::runtime_error{
          "Could not parse particle position and mass in file " + filename};
    Foam::point position;
    for (std::size_t dd = 0; dd < std::remove_reference_t<Geometry>::dim; ++dd)
      position[dd] = std::stod(split_line[dd]);
    this->_particle_maker.mass =
        std::stod(split_line[std::remove_reference_t<Geometry>::dim]);

    return {IC::Particle::State::make_position(position), -1};
  }

protected:
  std::string filename;
  std::ifstream input{io::open_read(filename)};
};
template <typename ParticleMaker, typename Geometry>
InitialCondition_PrescribedPositionsMasses(ParticleMaker &&, Geometry &&,
                                           std::string)
    -> InitialCondition_PrescribedPositionsMasses<ParticleMaker, Geometry>;

/**
   \class InitialCondition_PrescribedPositionsMassesTags PTOF/InitialCondition.h
   "CTRW/InitialCondition.h"
   \brief Initial condition at prescribed positions, masses, and tags, read from
   a file.
   \details File should contain particle positions, masses, and tags, one
   particle per line.
*/
template <typename ParticleMaker, typename Geometry>
struct InitialCondition_PrescribedPositionsMassesTags
    : public InitialCondition<ParticleMaker, Geometry> {
private:
  using IC = InitialCondition<ParticleMaker, Geometry>;

public:
  InitialCondition_PrescribedPositionsMassesTags(ParticleMaker &&particle_maker,
                                                 Geometry &&geometry,
                                                 std::string filename)
      : IC{std::forward<ParticleMaker>(particle_maker),
           std::forward<Geometry>(geometry)},
        filename{filename} {}

  typename IC::PositionAndHint make_position_and_hint() override {
    auto split_line = io::split_line(input);
    if (split_line.size() < std::remove_reference_t<Geometry>::dim + 2)
      throw std::runtime_error{
          "Could not parse particle position and mass in file " + filename};
    Foam::point position;
    for (std::size_t dd = 0; dd < std::remove_reference_t<Geometry>::dim; ++dd)
      position[dd] = std::stod(split_line[dd]);
    this->_particle_maker.mass =
        std::stod(split_line[std::remove_reference_t<Geometry>::dim]);
    this->_particle_maker.tag =
        std::stoul(split_line[std::remove_reference_t<Geometry>::dim + 1]);

    return {IC::Particle::State::make_position(position), -1};
  }

protected:
  std::string filename;
  std::ifstream input{io::open_read(filename)};
};
template <typename ParticleMaker, typename Geometry>
InitialCondition_PrescribedPositionsMassesTags(ParticleMaker &&, Geometry &&,
                                               std::string)
    -> InitialCondition_PrescribedPositionsMassesTags<ParticleMaker, Geometry>;

/**
   \class InitialCondition_PrescribedPositionsMassesTagsTimes
   PTOF/InitialCondition.h "CTRW/InitialCondition.h"
   \brief Initial condition at
   prescribed positions, masses, tags, and times, read from a file.
   \details File should contain particle positions, masses, tags, and times, one
   particle per line.
*/
template <typename ParticleMaker, typename Geometry>
struct InitialCondition_PrescribedPositionsMassesTagsTimes
    : public InitialCondition<ParticleMaker, Geometry> {
private:
  using IC = InitialCondition<ParticleMaker, Geometry>;

public:
  InitialCondition_PrescribedPositionsMassesTagsTimes(
      ParticleMaker &&particle_maker, Geometry &&geometry, std::string filename)
      : IC{std::forward<ParticleMaker>(particle_maker),
           std::forward<Geometry>(geometry)},
        filename{filename} {}

  typename IC::PositionAndHint make_position_and_hint() override {
    auto split_line = io::split_line(input);
    if (split_line.size() < std::remove_reference_t<Geometry>::dim + 3)
      throw std::runtime_error{
          "Could not parse particle position and mass in file " + filename};
    Foam::point position;
    for (std::size_t dd = 0; dd < std::remove_reference_t<Geometry>::dim; ++dd)
      position[dd] = std::stod(split_line[dd]);
    this->_particle_maker.mass =
        std::stod(split_line[std::remove_reference_t<Geometry>::dim]);
    this->_particle_maker.tag =
        std::stoul(split_line[std::remove_reference_t<Geometry>::dim + 1]);
    this->_particle_maker.time =
        std::stod(split_line[std::remove_reference_t<Geometry>::dim + 2]);

    return {IC::Particle::State::make_position(position), -1};
  }

protected:
  std::string filename;
  std::ifstream input{io::open_read(filename)};
};
template <typename ParticleMaker, typename Geometry>
InitialCondition_PrescribedPositionsMassesTagsTimes(ParticleMaker &&,
                                                    Geometry &&, std::string)
    -> InitialCondition_PrescribedPositionsMassesTagsTimes<ParticleMaker,
                                                           Geometry>;

/**
   \class InitialCondition_DistributedNearFaces PTOF/InitialCondition.h
   "CTRW/InitialCondition.h"
   \brief Initial condition randomly distributed near boundary faces.
   \tparam FaceIdContainer Container of face indices.
*/
template <typename ParticleMaker, typename Geometry, typename FaceIdContainer,
          typename Distribution>
struct InitialCondition_DistributedNearFaces
    : public InitialCondition<ParticleMaker, Geometry> {
private:
  using IC = InitialCondition<ParticleMaker, Geometry>;

public:
  InitialCondition_DistributedNearFaces(ParticleMaker &&particle_maker,
                                        Geometry &&geometry,
                                        FaceIdContainer &&face_ids,
                                        Distribution &&dist, double distance,
                                        std::size_t nr_tries)
      : IC{std::forward<ParticleMaker>(particle_maker),
           std::forward<Geometry>(geometry)},
        _face_ids{std::forward<FaceIdContainer>(face_ids)},
        _dist{std::forward<Distribution>(dist)}, _distance{distance},
        _nr_tries{nr_tries} {
    if (_face_ids.size() == 0)
      throw std::runtime_error{"No faces provided for initial condition"};
  }

  typename IC::PositionAndHint make_position_and_hint() override {
    for (std::size_t ii = 0; ii < _nr_tries; ++ii) {
      auto const &locator = this->_geometry.locator;
      auto const &mesh = locator.mesh();
      auto face_id = _face_ids[_dist(_rng)];
      auto position = face_center(face_id, mesh) +
                      _distance * unit_normal_inward(face_id, mesh);
      auto owner_cell = mesh.faceOwner()[face_id];
      auto cell_id = locator(position, owner_cell);
      if (outside(cell_id))
        continue;
      auto nearest_face_id =
          locator.mesh_search().findNearestBoundaryFace(position, face_id);
      if (nearest_face_id == face_id)
        return {IC::Particle::State::make_position(position), cell_id};
    }

    throw std::runtime_error{"Failed to make position inside mesh and nearest "
                             "to a selected face after " +
                             std::to_string(_nr_tries) + " tries"};
  }

protected:
  stochastic::RNGThreaded<
      typename std::remove_reference_t<Geometry>::ParallelOption, std::mt19937>
      _rng;
  FaceIdContainer _face_ids;
  Distribution _dist;
  double _distance;
  std::size_t _nr_tries;
};
template <typename ParticleMaker, typename Geometry, typename FaceIdContainer,
          typename Distribution>
InitialCondition_DistributedNearFaces(ParticleMaker &&, Geometry &&,
                                      FaceIdContainer &&, Distribution &&,
                                      double, std::size_t)
    -> InitialCondition_DistributedNearFaces<ParticleMaker, Geometry,
                                             FaceIdContainer, Distribution>;

/**
   \class InitialCondition_DistributedNearFaces PTOF/InitialCondition.h
   "CTRW/InitialCondition.h"
   \brief Initial condition uniformly randomly distributed near boundary faces.
   \tparam FaceIdContainer Container of face indices.
*/
template <typename ParticleMaker, typename Geometry, typename FaceIdContainer>
struct InitialCondition_UniformNearFaces
    : public InitialCondition_DistributedNearFaces<
          ParticleMaker, Geometry, FaceIdContainer,
          std::discrete_distribution<std::size_t>> {
private:
  using IC = InitialCondition_DistributedNearFaces<
      ParticleMaker, Geometry, FaceIdContainer,
      std::discrete_distribution<std::size_t>>;

public:
  InitialCondition_UniformNearFaces(ParticleMaker &&particle_maker,
                                    Geometry &&geometry,
                                    FaceIdContainer &&face_ids, double distance,
                                    std::size_t nr_tries)
      : IC{std::forward<ParticleMaker>(particle_maker),
           std::forward<Geometry>(geometry),
           std::forward<FaceIdContainer>(face_ids),
           uniform_face_distribution(face_ids, geometry.mesh()),
           distance,
           nr_tries} {}
};
template <typename ParticleMaker, typename Geometry, typename CellIdContainer>
InitialCondition_UniformNearFaces(ParticleMaker &&, Geometry &&,
                                  CellIdContainer &&)
    -> InitialCondition_UniformNearFaces<ParticleMaker, Geometry,
                                         CellIdContainer>;

/**
   \class InitialCondition_DistributedNearFaces PTOF/InitialCondition.h
   "CTRW/InitialCondition.h"
   \brief Initial condition flux-weighted randomly distributed near boundary
   faces.
   \tparam FaceIdContainer Container of face indices.
*/
template <typename ParticleMaker, typename Geometry, typename FaceIdContainer>
struct InitialCondition_FluxweightedNearFaces
    : public InitialCondition_DistributedNearFaces<
          ParticleMaker, Geometry, FaceIdContainer,
          std::discrete_distribution<std::size_t>> {
private:
  using IC = InitialCondition_DistributedNearFaces<
      ParticleMaker, Geometry, FaceIdContainer,
      std::discrete_distribution<std::size_t>>;

public:
  template <typename VelocityField>
  InitialCondition_FluxweightedNearFaces(ParticleMaker &&particle_maker,
                                         Geometry &&geometry,
                                         FaceIdContainer &&face_ids,
                                         VelocityField const &velocity_field,
                                         double distance, std::size_t nr_tries)
      : IC{std::forward<ParticleMaker>(particle_maker),
           std::forward<Geometry>(geometry),
           std::forward<FaceIdContainer>(face_ids),
           fluxweighted_face_distribution(face_ids, velocity_field,
                                          geometry.locator),
           distance,
           nr_tries} {}
};
template <typename ParticleMaker, typename Geometry, typename FaceIdContainer,
          typename VelocityField>
InitialCondition_FluxweightedNearFaces(ParticleMaker &&, Geometry &&,
                                       FaceIdContainer &&,
                                       VelocityField const &, double,
                                       std::size_t)
    -> InitialCondition_FluxweightedNearFaces<ParticleMaker, Geometry,
                                              FaceIdContainer>;
} // namespace ptof
#endif /* PTOF_INITIALCONDITION_H */
