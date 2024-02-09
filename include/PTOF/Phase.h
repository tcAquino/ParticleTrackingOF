/**
* \file PTOF/Phase.h
* \author Tomás Aquino
* \date 08/02/2024
*/

#ifndef PTOF_PHASE_H
#define PTOF_PHASE_H

#include <cmath>
#include <exception>
#include <fieldTypes.H>
#include <fvcGrad.H>
#include <string>
#include <volFieldsFwd.H>
#include "General/Useful.h"
#include "PTOF/Directories.h"

namespace ptof
{
  struct Phase
  {
    using PhaseField = Foam::volScalarField;
    using GradPhaseField = Foam::volVectorField;
    
    struct Parameters
    {
      std::string excluded_phase_name;
      double leakage_coefficient;
      bool excluded_phase;
      bool compute_gradient;
      double phase_threshold;
      
      Parameters
      (Directories const& directories, std::string const& name)
      {
        auto input = useful::open_read(directories.dir_parameters
                                       + "/parameters_phase_"
                                       + name + ".dat");
        useful::read(input, excluded_phase_name);
        double leakage_tolerance = 1e-3;
        std::string phase_transport;
        useful::read(input, phase_transport);
        if (phase_transport == "excluded")
        {
          excluded_phase = true;
        }
        else if (phase_transport == "transport")
        {
          excluded_phase = false;
        }
        else
          throw std::invalid_argument{
            std::string("Unknown option ") + phase_transport +
            " for phase type" };
        std::string compute;
        useful::read(input, compute);
        if (compute == "compute")
        {
          compute_gradient = true;
        }
        else if (compute == "read")
        {
          compute_gradient = false;
        }
        else
          throw std::invalid_argument{
            std::string("Unknown option ") + compute +
            " for whether to compute or read excluded phase gradient" };
        useful::read(input, leakage_tolerance);
        leakage_coefficient = std::log(leakage_tolerance);
        useful::read(input, phase_threshold);
        input.close();
      }
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Phase parameters\n"
          "--------------------------------------------------\n"
          "- Name of phase\n"
          "- Excluded or transport phase\n"
          "\texcluded: No transport in this phase\n"
          "\ttransport: Transport in this phase\n"
          "- Wheteher to read or compute gradient of excluded phase\n"
          "  (ignored if transport phase is given)\n"
          "\tread: Read gradient\n"
          "\tcompute: Compute gradient\n"
          "- Leakage tolerance\n"
          "- Phase field value threshold to consider phase is present\n"
          "--------------------------------------------------\n";
      }
    };
    
    template <typename Mesh>
    static auto get_excluded_phase_data
    (Mesh const& mesh, Parameters const& parameters)
    {
      if (parameters.excluded_phase)
        return PhaseField{
          Foam::IOobject{
            std::string("alpha.") + parameters.excluded_phase_name,
            mesh.time().timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE
          },
          mesh };
      else
        return PhaseField(1. - Foam::volScalarField{
          Foam::IOobject{
            std::string("alpha.") + parameters.excluded_phase_name,
            mesh.time().timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE
          },
          mesh });
    }
    
    template <typename Mesh>
    static auto get_grad_excluded_phase_data
    (Mesh const& mesh, Parameters const& parameters, Foam::volScalarField const& excluded_phase_field)
    {
      if (parameters.excluded_phase && !parameters.compute_gradient)
        return GradPhaseField{
          Foam::IOobject{
            std::string("gradAlpha.") + parameters.excluded_phase_name,
            mesh.time().timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE
          },
          mesh };
      else
        return Foam::volVectorField(Foam::fvc::grad(excluded_phase_field));
    }
    
    template <typename Mesh>
    static auto get_grad_excluded_phase_data
    (Mesh const& mesh, Parameters const& parameters)
    {
      return GradPhaseField{
        Foam::IOobject{
          std::string("alpha.") + parameters.excluded_phase_name,
          mesh.time().timeName(),
          mesh,
          Foam::IOobject::MUST_READ,
          Foam::IOobject::NO_WRITE
        },
        mesh };
    }
  };
}

#endif /* PTOF_PHASE_H */
