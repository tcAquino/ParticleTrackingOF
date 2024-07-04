# ParticleTrackingOF (PTOF)

Particle tracking based on OpenFOAM-generated meshes and flow fields, including advection, diffusion, and reaction processes. This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software and owner of the OPENFOAM®  and OpenCFD® trade marks.

## Dependencies

Required external libraries (free and open-source):
- OpenFOAM : CFD software (https://www.openfoam.com/)

A version of the g++ or clang++ compiler compatible with C++17 and OpenMP v4.5 are required. Note that the same compiler used for OpenFOAM should be used to avoid (dynamic) linking issues. For compilation of parallel implementations with clang++, llvm's OpenMP must be installed (to install from the terminal, run "brew install libomp" on MacOS or "sudo apt install libomp-dev" on Linux).

## Documentation

Basic information on code structure, compilation, and execution is given in doc/Manual.pdf. Doxygen documentation can be found in doc/Doxygen.

## Notes and known issues

- If the underlying OpenFOAM case uses dynamicCode functionalities, they are compiled at runtime by OpenFOAM when running PTOF executables on that case for the first time and placed on a dynamicCode folder in the same folder as the executable. If more than one instance of PTOF is run simultaneously during this first run, this can cause a conflict and lead to a crash. A simple workaround is to wait for dynamic compilation, which happens early on, before launching a second instance on the first run.
- Parall simulations employ parallelization over particles across the whole spatial domain. This means that parallel OpenFOAM cases must be reconstructed before using PTOF, and that each thread must have access to the full mesh. Unfortunatlely, OpenFOAM is not designed with this type of parallelization in mind, and even mesh searching methods that are marked const can make changes to mutable members of the mesh object. At least for the moment, this requires each thread to have its own dedicated copy of the full mesh. Because of this, you should be careful of RAM consumption when running parallel simulations.

## Citations

Please link to this repository (https://github.com/tcAquino/ParticleTrackingOF).

## License

License information can be found in the file LICENSE.txt.

## About the author

Tomás Aquino is a researcher at IDAEA-CSIC, Barcelona, Spain. He can be contacted at tomas.aquino@idaea.csic.es.

## Contributors

- Tomás Aquino (main developer)
- Nolwenn Delouche (testing and original Doxygen documentation setup)
- Yang Liu (testing)
