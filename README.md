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

This docker container supports building on both x86 and Arm-based architectures.
The CI only builds the x86 version, but the above command will build the correct version for your architecture.
