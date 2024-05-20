# cmaq-container

Build immutable images for running the CMAQ adjoint for OpenMethane.

OpenMethane uses a customised CH4-only version of CMAQ based on CMAQ v5.0.1.
The adjunct code and it's customisations are not publically available,
so have been vendored into this repository.

## Requirements

[Docker](https://www.docker.com/) is required.

    
### Local development


To build the docker image, run the following command:

```
    docker build . -t cmaq
```

If building on a Arm-based Mac, the following command should be used to target the correct platform.
Otherwise, the build will emulate the amd64 platform and take much longer.

```
    docker build --build_arg PLATFORM=linux/arm64 . -t cmaq
```
