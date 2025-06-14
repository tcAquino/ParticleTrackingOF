# ParticleTrackingOF (PTOF)

Particle tracking based on OpenFOAM-generated meshes and flow fields, including advection, diffusion, and reaction processes. This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software and owner of the OPENFOAM®  and OpenCFD® trade marks.

## Dependencies

PTOF is meant to be used in a Linux system. It has mostly been tested under Ubuntu v22.04. It is also expected to work on MacOS, but is less tested there.

Required external libraries (free and open-source):
- OpenFOAM : CFD software (https://www.openfoam.com/)

A version of the g++ or clang++ compiler compatible with C++17 and OpenMP v4.5 are required. Note that the same compiler used for OpenFOAM should be used to avoid (dynamic) linking issues. For compilation of parallel implementations with clang++, llvm's OpenMP must be installed (to install from the terminal, run "brew install libomp" on MacOS or "sudo apt install libomp-dev" on Linux).

## Documentation

Basic information on code structure, compilation, and execution is given in doc/Manual.pdf. Doxygen documentation can be found in doc/Doxygen. The compiled doxygen documentation is not currently online. To compile and view it, make sure you have doxygen, graphviz, and latex (e.g., texlive-full) installed (all free and available in common repositories like apt), go to the doc folder, and execute doxygen Doxyfile. Then open the doc/html/index.html file to view it in a browser.

## Notes and known issues

- If the underlying OpenFOAM case uses dynamicCode functionalities, they are compiled at runtime by OpenFOAM when running PTOF executables on that case for the first time and placed on a dynamicCode folder in the same folder as the executable. If more than one instance of PTOF is run simultaneously during this first run, this can cause a conflict and lead to a crash. A simple workaround is to wait for dynamic compilation, which happens early on, before launching a second instance on the first run.
- Parallel simulations employ parallelization over particles across the whole spatial domain. This means that parallel OpenFOAM cases must be reconstructed before using PTOF, and that each thread must have access to the full mesh. OpenFOAM is not designed with this type of parallelization in mind. Demand-driven mesh and mesh search tool data must be precomputed in serial before being used in parallel. When implementing extensions, be careful to verify that all necessary data is precomputed.
- Attempting to compile different models or parallel settings corresponding to the same cpp file at the same time is not safe, as it can lead to concurrent modifications to the same file.
- Due to floating point errors in reading and writing, positions that are inside the mesh may be outside after writing and subsequent reading. It does not appear to be possible to fix this in general by increasing output precision. This means that, in general, a position that is read from a file cannot be assumed to be inside the mesh.

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
