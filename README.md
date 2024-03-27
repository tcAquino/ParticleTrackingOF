# ParticleTrackingOF

Particle tracking based on OpenFOAM-generated meshes and flow fields, including advection, diffusion, and reaction processes. This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software and owner of the OPENFOAM®  and OpenCFD® trade marks.

## Dependencies

Required external libraries (free and open-source):
- OpenFOAM : CFD software (https://www.openfoam.com/)

A version of the g++ or clang++ compiler compatible with C++17 and OpenMP v4.5 are required. Note that the same compiler used for OpenFOAM should be used to avoid (dynamic) linking issues. For compilation of parallel implementations with clang++, llvm's OpenMP must be installed (to install from the terminal, run "brew install libomp" on MacOS or "sudo apt install libomp-dev" on Linux).

## Documentation

Basic information on code structure, compilation, and execution is given in man/Manual.pdf. Doxygen documentation can be found in man/Doxygen/.

## Citations

Please link to this repository (https://github.com/tcAquino/ParticleTrackingOF).

## License

License information can be found in the file LICENSE.txt.

## About the author

Tomás Aquino is a researcher at IDAEA-CSIC, Barcelona, Spain. He can be contacted at tomas.aquino@idaea.csic.es.

## Collaborators

- Nolwenn Delouche
