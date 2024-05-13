/**
 \file PTOF/Phase.h
 \author Tomás Aquino
 \date 08/02/2024
*/

#ifndef PTOF_PHASE_H
#define PTOF_PHASE_H

#include <cmath>
#include <stdexcept>
#include <fieldTypes.H>
#include <fvcGrad.H>
#include <string>
#include <volFieldsFwd.H>
#include "General/Useful.h"
#include "PTOF/Directories.h"

namespace ptof
{
  /** \struct Phase PTOF/Phase.h "PTOF/Phase.h"
   \brief Phase field information. */
  struct Phase
  {
    using PhaseField = Foam::volScalarField;      /**< Saturation scalar field as a function of mesh cell index. */
    using GradPhaseField = Foam::volVectorField;  /**< Gradient of saturation vector field as a function of mesh cell index. */
    
    /** \struct Phase::Parameters PTOF/Phase.h "PTOF/Phase.h"
     \brief Parameters for phase-related quantities. */
    struct Parameters
    {
    private:
      std::string comment_sequence = "#"; /**< Sequence of characters marking comment for file parsing. */
      
    public:
      std::string phase_name;             /**< Name of phase saturation field to be read from file. */
      bool excluded_phase;                /**< Weather named phase is excluded phase or carrier phase from file. */
      bool compute_gradient;              /**< Weather to compute or read from file gradient of excluded phase. */
      double leakage_coefficient;         /**< Log of ratio of concentrations between excluded phase and carrier phase. */
      double phase_threshold;             /**< Tolerance in local phase saturation to consider pure phase. */
      
      /** Constructor.
       \param directories Current case directory information.
       \param name Name of phase parameters set.
      */
      template <typename Geometry>
      Parameters
      (Directories const& directories, std::string const& name,
       Geometry const& geometry)
      {
        auto input = useful::open_read(directories.dir_parameters
                                       + "/parameters_phase_"
                                       + name + ".dat");
        useful::read_first_from_line(input, phase_name, comment_sequence);
        auto phase_transport = useful::read_first_from_line<std::string>(input, comment_sequence);
        if (phase_transport == "excluded")
        {
          excluded_phase = true;
        }
        else if (phase_transport == "carrier")
        {
          excluded_phase = false;
        }
        else
          throw std::invalid_argument{
            std::string("Unknown option ") + phase_transport +
            " for phase type" };
        auto compute = useful::read_first_from_line<std::string>(input, comment_sequence);
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
        auto leakage_tolerance = useful::read_first_from_line<double>(input, comment_sequence);
        leakage_coefficient = std::log(leakage_tolerance);
        useful::read_first_from_line(input, phase_threshold, comment_sequence);
        input.close();
      }
      
      /** \brief Output general information about object. */
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Phase parameters\n"
          "--------------------------------------------------\n"
          "- Name of phase to be read from file\n"
          "- Whether named phase is excluded or carrier phase\n"
          "\texcluded: No transport in this phase\n"
          "\tcarrier: Transport in this phase\n"
          "- Wheteher to read or compute gradient of excluded phase\n"
          "  (ignored if carrier phase is given)\n"
          "\tread: Read gradient\n"
          "\tcompute: Compute gradient\n"
          "- Leakage tolerance\n"
          "- Phase field value threshold to consider phase is present\n"
          "--------------------------------------------------\n";
      }
    };
    
    /** \return Excluded phase saturation field.
     \details Obtained from OpenFOAM file data for either the excluded phase or the carrier phase. */
    template <typename Mesh>
    static auto get_excluded_phase_data
    (Mesh const& mesh, Parameters const& parameters)
    {
      if (parameters.excluded_phase)
        return PhaseField{
          Foam::IOobject{
            std::string("alpha.") + parameters.phase_name,
            mesh.time().timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE },
          mesh };
      else
        return static_cast<PhaseField>(1. - PhaseField{
          Foam::IOobject{
            std::string("alpha.") + parameters.phase_name,
            mesh.time().timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE },
          mesh });
    }
    
    /** \return Excluded phase field gradient.
    \details Based on OpenFOAM file data. */
    template <typename Mesh>
    static auto get_grad_excluded_phase_data
    (Mesh const& mesh, Parameters const& parameters)
    {
      if (parameters.excluded_phase)
        return GradPhaseField{
          Foam::IOobject{
            std::string("gradAlpha.") + parameters.phase_name,
            mesh.time().timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE },
          mesh };
      else
        return static_cast<GradPhaseField>(-GradPhaseField{
          Foam::IOobject{
            std::string("gradAlpha.") + parameters.phase_name,
            mesh.time().timeName(),
            mesh,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE },
          mesh });
    }
    
    /** \return Excluded phase field gradient.
     \details Based on OpenFOAM file data or computed from excluded phase saturation field. */
    template <typename Mesh>
    static auto grad_excluded_phase
    (Mesh const& mesh, Parameters const& parameters, PhaseField const& excluded_phase_field)
    {
      if (parameters.compute_gradient)
        return static_cast<GradPhaseField>(Foam::fvc::grad(excluded_phase_field));
      else
        return get_grad_excluded_phase_data(mesh, parameters);
    }
  };
}

#endif /* PTOF_PHASE_H */
