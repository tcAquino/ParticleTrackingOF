/**
   \file PTOF/ReactionHandler.h
   \author Tomas Aquino
   \date 10/03/2022
   \brief Objects for handling surface and bulk reactions.
*/

#ifndef PTOF_REACTIONHANDLER_H
#define PTOF_REACTIONHANDLER_H

#include "CTRW/TimeGenerator.h"
#include "General/IO.h"
#include "General/Operation.h"
#include "PTOF/BulkReaction.h"
#include "PTOF/Directories.h"
#include "PTOF/SurfaceReaction.h"
#include <limits>
#include <ostream>
#include <random>
#include <string>
#include <type_traits>
#include <vector>

namespace ptof {
struct ReactionHandler_NoBulk_NoSurface {
  ReactionHandler_NoBulk_NoSurface() = delete;

  using BulkReaction = BulkReaction_DoNothing;
  using SurfaceReaction = SurfaceReaction_DoNothing;

  struct Parameters {
    double damkohler{0.};
    double rate_constant{0.};
    double reaction_time{std::numeric_limits<double>::infinity()};

    Parameters() {}

    template <typename Geometry, typename TransportParameters>
    Parameters(Directories const &dir, std::string const &parameter_set_name,
               Geometry const &geometry,
               TransportParameters const &params_transport) {}

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output << io::line() << "Reaction parameters\n"
             << io::line() << "None\n"
             << io::line();
      return output;
    }
  };

  template <typename Geometry, typename TransportParameters,
            typename SolverParameters>
  static auto makeBulkReaction(Geometry const &geometry,
                               Parameters const &params_reaction,
                               TransportParameters const &params_transport,
                               SolverParameters const &params_solvers) {
    return BulkReaction{};
  }

  template <typename Geometry, typename TransportParameters,
            typename SolverParameters>
  static auto makeSurfaceReaction(Geometry const &geometry,
                                  Parameters const &params_reaction,
                                  TransportParameters const &params_transport,
                                  SolverParameters const &params_solvers) {
    return SurfaceReaction{};
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    BulkReaction::info(output) << "\n";
    SurfaceReaction::info(output);
    return output;
  }

  template <typename CTRW>
  static void update(BulkReaction &, SurfaceReaction &, CTRW &,
                     typename CTRW::State::Time, typename CTRW::State::Time) {}
};

template <typename Geometry, bool solid_decay, typename ParallelOption>
struct ReactionHandler_NoBulk_SurfaceDecay {
  ReactionHandler_NoBulk_SurfaceDecay() = delete;

  using BulkReaction = BulkReaction_DoNothing;
  using Locator = typename Geometry::Locator;
  using SurfaceReaction =
      std::conditional_t<solid_decay,
                         SurfaceReaction_AFluidPlusASolidtoNothing<
                             Locator const &, ParallelOption>,
                         SurfaceReaction_AFluidPlusASolidtoASolid>;

  struct Parameters {
  public:
    double damkohler;
    std::string initial_distribution;
    std::vector<std::string> patch_names;
    std::vector<double> surface_concentrations;
    double rate_constant;
    double reaction_time;

    double rate_constant_ratio_solid_to_fluid = 0.;

    template <typename TransportParameters>
    Parameters(Directories const &directories,
               std::string const &parameter_set_name, Geometry const &geometry,
               TransportParameters const &params_transport) {
      std::string filename = directories.dir_parameters + "/reaction_" +
                             parameter_set_name + ".param";
      auto input = io::open_read(filename);
      std::string in_file = std::string{"In file "} + filename + " : ";

      double average_surface_concentration = 0.;
      auto split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
      std::size_t param_index = 0;
      io::read(split_line, param_index,
               in_file +
                   "Could not parse initial surface reactant distribution type",
               initial_distribution);
      std::string for_initial_distribution =
          std::string{"Initial surface reactant distribution type "} +
          initial_distribution + " : ";
      if (initial_distribution == "uniform_patch") {
        std::size_t nr_patches = (split_line.size() - param_index) / 2;
        if (nr_patches == 0) {
          throw std::runtime_error{
              in_file + for_initial_distribution +
              "Could not parse patch names(s) and surface concentration(s): "
              "At least one patch and surface concentration required"};
        }
        patch_names.reserve(nr_patches);
        surface_concentrations.reserve(nr_patches);
        for (std::size_t ii = 0; ii < nr_patches; ++ii) {
          patch_names.push_back(
              io::read<std::string>(split_line, param_index,
                                    in_file + for_initial_distribution +
                                        "Could not parse patch name"));
          surface_concentrations.push_back(
              io::read<double>(split_line, param_index,
                               in_file + for_initial_distribution +
                                   "Could not parse initial surface "
                                   "concentration for patch " +
                                   patch_names[ii]));
        }
        auto areas = patch_areas(patch_names, geometry.mesh());
        average_surface_concentration =
            op::sum(op::times(surface_concentrations, areas)) / op::sum(areas);
      } else {
        throw std::runtime_error{in_file + for_initial_distribution +
                                 "Not supported"};
      }

      split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
      param_index = 0;
      auto rate_distribution = io::read<std::string>(
          split_line, param_index,
          in_file + "Could not parse rate specification type");
      std::string for_rate_distribution =
          std::string{"Rate distibution type "} + rate_distribution + " : ";
      if (rate_distribution == "uniform") {
        auto rate_option =
            io::read<std::string>(split_line, param_index,
                                  in_file + for_rate_distribution +
                                      "Could not parse rate format option");
        std::string for_rate_option =
            std::string{"Rate specification type "} + rate_option + " : ";
        if (rate_option == "damkohler") {
          io::read(split_line, param_index,
                   in_file + for_rate_distribution + for_rate_option +
                       "Could not parse Damkohler number",
                   damkohler);
          rate_constant =
              params_transport.lengthscale * damkohler /
              (average_surface_concentration * params_transport.diffusion_time);
        } else if (rate_option == "rate_constant") {
          io::read(split_line, param_index,
                   in_file + for_rate_distribution + for_rate_option +
                       "Could not parse rate constant",
                   rate_constant);
          damkohler = average_surface_concentration *
                      params_transport.diffusion_time * rate_constant /
                      params_transport.lengthscale;
        } else {
          throw std::runtime_error{in_file + for_rate_distribution +
                                   for_rate_option + "Not supported"};
        }

        if constexpr (solid_decay) {
          io::read(split_line, param_index,
                   in_file + for_rate_distribution +
                       "Could not parse ratio of solid to fluid rate constants",
                   rate_constant_ratio_solid_to_fluid);
        }

        reaction_time = params_transport.lengthscale /
                        (average_surface_concentration * rate_constant);
      } else {
        throw std::runtime_error{in_file + for_rate_distribution +
                                 "Not supported"};
      }
    }

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << io::line() << "Reaction parameters:\n"
          << io::line()
          << "- Solid reactant initial distribution type:\n"
             "  - uniform_patch\n"
             "    - Homogeneous in each specified boundary patch\n"
             "    - Pass on same line:\n"
             "      - Pairs of boundary patch names and surface\n"
             "        concentrations\n"
             "- Rate specification type:\n"
             "  - uniform\n"
             "    - Homogeneous and the same in all boundary patches\n"
             "    - Pass on same line:\n"
             "      - Rate specification format:\n"
             "        - damkohler\n"
             "          - Pass on same line:\n"
             "            - Damkohler number ([diffusion time] * [average\n"
             "              surface concentration] * [surface rate constant]\n"
             "              / [lengthscale])\n";
      if constexpr (solid_decay) {
        output << "            - Ratio of solid to fluid rate constants\n";
      }
      output
          << "            (Note:\n"
             "              - Diffusion time is as defined in transport\n"
             "                parameters)\n"
             "        - rate_constant\n"
             "          - Pass on same line:\n"
             "            - Surface rate constant (units length per surface\n"
             "              concentration per time)\n";
      if constexpr (solid_decay) {
        output << "            - Ratio of solid to fluid rate constants\n";
      }
      output << io::line();
      return output;
    }
  };

  template <typename TransportParameters, typename SolverParameters>
  static auto makeBulkReaction(Geometry const &geometry,
                               Parameters const &params_reactions,
                               TransportParameters const &params_transport,
                               SolverParameters const &params_solvers) {
    return BulkReaction{};
  }

  template <typename TransportParameters, typename SolverParameters>
  static auto makeSurfaceReaction(Geometry const &geometry,
                                  Parameters const &params_reaction,
                                  TransportParameters const &params_transport,
                                  SolverParameters const &params_solvers) {
    if constexpr (solid_decay) {
      return SurfaceReaction{
          params_reaction.rate_constant,
          params_reaction.rate_constant_ratio_solid_to_fluid,
          params_transport.diff_coeff,
          uniform_solid_reactant(params_reaction.patch_names,
                                 params_reaction.surface_concentrations,
                                 geometry.mesh()),
          geometry.locator,
          ParallelOption{}};
    } else {
      return SurfaceReaction{
          params_reaction.rate_constant, params_transport.diff_coeff,
          uniform_solid_reactant(params_reaction.patch_names,
                                 params_reaction.surface_concentrations,
                                 geometry.mesh())};
    }
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    BulkReaction::info(output) << "\n";
    SurfaceReaction::info(output);
    return output;
  }

  template <typename CTRW>
  static void update(BulkReaction &bulk_reaction,
                     SurfaceReaction &surface_reaction, CTRW &ctrw,
                     typename CTRW::State::Time time_new,
                     typename CTRW::State::Time time_old) {
    if constexpr (solid_decay) {
      surface_reaction.update(ctrw, time_new, time_old);
    }
  }
};

template <typename Geometry, typename ParallelOption>
struct ReactionHandler_NoBulk_SurfaceAdsorption {
  ReactionHandler_NoBulk_SurfaceAdsorption() = delete;

  using BulkReaction = BulkReaction_DoNothing;
  using Locator = typename Geometry::Locator;
  using SurfaceReaction = SurfaceReaction_LinearSorption<
      ParallelOption,
      ctrw::TimeGenerator_Dist<std::exponential_distribution<double>,
                               ParallelOption>>;

  struct Parameters {
  public:
    double damkohler;
    std::string initial_distribution;
    std::vector<std::string> patch_names;
    std::vector<double> surface_concentrations;
    double rate_constant;
    double desorption_rate;
    double reaction_time;

    template <typename TransportParameters>
    Parameters(Directories const &directories,
               std::string const &parameter_set_name, Geometry const &geometry,
               TransportParameters const &params_transport) {
      std::string filename = directories.dir_parameters + "/reaction_" +
                             parameter_set_name + ".param";
      auto input = io::open_read(filename);
      std::string in_file = std::string{"In file "} + filename + " : ";

      double average_surface_concentration = 0.;
      auto split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
      std::size_t param_index = 0;
      io::read(split_line, param_index,
               in_file +
                   "Could not parse initial surface reactant distribution type",
               initial_distribution);
      std::string for_initial_distribution =
          std::string{"Initial surface reactant distribution type "} +
          initial_distribution + " : ";
      if (initial_distribution == "uniform_patch") {
        std::size_t nr_patches = (split_line.size() - param_index) / 2;
        if (nr_patches == 0) {
          throw std::runtime_error{
              in_file + for_initial_distribution +
              "Could not parse patch names(s) and surface concentration(s): "
              "At least one patch and surface concentration required"};
        }
        patch_names.reserve(nr_patches);
        surface_concentrations.reserve(nr_patches);
        for (std::size_t ii = 0; ii < nr_patches; ++ii) {
          patch_names.push_back(
              io::read<std::string>(split_line, param_index,
                                    in_file + for_initial_distribution +
                                        "Could not parse patch name"));
          surface_concentrations.push_back(
              io::read<double>(split_line, param_index,
                               in_file + for_initial_distribution +
                                   "Could not parse initial surface "
                                   "concentration for patch " +
                                   patch_names[ii]));
        }
        auto areas = patch_areas(patch_names, geometry.mesh());
        average_surface_concentration =
            op::sum(op::times(surface_concentrations, areas)) / op::sum(areas);
      } else {
        throw std::runtime_error{in_file + for_initial_distribution +
                                 "Not supported"};
      }

      split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
      param_index = 0;
      auto rate_distribution = io::read<std::string>(
          split_line, param_index,
          in_file + "Could not parse rate specification type");
      std::string for_rate_distribution =
          std::string{"Rate specification type "} + rate_distribution + " : ";
      if (rate_distribution == "uniform") {
        auto rate_option =
            io::read<std::string>(split_line, param_index,
                                  in_file + for_rate_distribution +
                                      "Could not parse rate format option");
        std::string for_rate_option =
            std::string{"Rate specification type "} + rate_option + " : ";
        if (rate_option == "damkohler") {
          io::read(split_line, param_index,
                   in_file + for_rate_distribution + for_rate_option +
                       "Could not parse Damkohler number",
                   damkohler);
          rate_constant =
              params_transport.lengthscale * damkohler /
              (average_surface_concentration * params_transport.diffusion_time);
          io::read(split_line, param_index,
                   in_file + for_rate_distribution + for_rate_option +
                       "Could not parse desorption rate",
                   desorption_rate);
        } else if (rate_option == "rate_constant") {
          io::read(split_line, param_index,
                   in_file + for_rate_distribution + for_rate_option +
                       "Could not parse surface rate constant",
                   rate_constant);
          damkohler = average_surface_concentration *
                      params_transport.diffusion_time * rate_constant /
                      params_transport.lengthscale;
          auto damkohler_desorption = io::read<double>(
              split_line, param_index,
              in_file + for_rate_distribution + for_rate_option +
                  "Could not parse desorption Damkohler number");
          desorption_rate =
              damkohler_desorption / params_transport.diffusion_time;
        } else {
          throw std::runtime_error{in_file + for_rate_distribution +
                                   for_rate_option + "Not supported"};
        }

        reaction_time = params_transport.lengthscale /
                        (average_surface_concentration * rate_constant);
      } else {
        throw std::runtime_error{in_file + for_rate_distribution +
                                 "Not supported"};
      }
    }

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << io::line() << "Reaction parameters:\n"
          << io::line()
          << "- Solid reactant initial distribution type:\n"
             "  - uniform_patch\n"
             "    - Homogeneous in each specified boundary patch\n"
             "    - Pass on same line:\n"
             "      - Pairs of boundary patch names and surface\n"
             "        concentrations\n"
             "- Rate distribution type:\n"
             "  - uniform\n"
             "    - Homogeneous and the same in all boundary patches\n"
             "    - Pass on same line:\n"
             "      - Rate specification format:\n"
             "        - damkohler\n"
             "          - Pass on same line:\n"
             "            - Damkohler number ([diffusion time] * [average\n"
             "              surface concentration] * [surface rate constant]\n"
             "              / [lengthscale])\n"
             "            - Desorption Damkohler number ([diffusion time] * \n"
             "              [desorption rate])\n"
             "            (Note:\n"
             "              - Diffusion time is as defined in transport\n"
             "                parameters)\n"
             "        - rate_constant\n"
             "          - Pass on same line:\n"
             "            - Surface rate constant (units length per surface\n"
             "              concentration per time)\n"
             "            - Desorption rate\n"
          << io::line();
      return output;
    }
  };

  template <typename TransportParameters, typename SolverParameters>
  static auto makeBulkReaction(Geometry const &geometry,
                               Parameters const &params_reactions,
                               TransportParameters const &params_transport,
                               SolverParameters const &params_solvers) {
    return BulkReaction{};
  }

  template <typename TransportParameters, typename SolverParameters>
  static auto makeSurfaceReaction(Geometry const &geometry,
                                  Parameters const &params_reaction,
                                  TransportParameters const &params_transport,
                                  SolverParameters const &params_solvers) {
    return SurfaceReaction{
        params_reaction.rate_constant, params_transport.diff_coeff,
        uniform_solid_reactant(params_reaction.patch_names,
                               params_reaction.surface_concentrations,
                               geometry.mesh()),
        params_reaction.desorption_rate};
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    BulkReaction::info(output) << "\n";
    SurfaceReaction::info(output);
    return output;
  }

  template <typename CTRW>
  static void update(BulkReaction &bulk_reaction,
                     SurfaceReaction &surface_reaction, CTRW &ctrw,
                     typename CTRW::State::Time time_new,
                     typename CTRW::State::Time time_old) {}
};
} // namespace ptof

#endif /* PTOF_REACTIONHANDLER_H */
