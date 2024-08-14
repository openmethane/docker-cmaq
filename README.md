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
    make build
```

This docker container supports building on both x86 and Arm-based architectures.
The CI only builds the x86 version, but the above command will build the correct version for your architecture.
The Arm version is only used for local testing and may not work.


### Testing

A very simple test suite is present to run the adjoint forward and backwards with
minimal inputs.

These tests should be run in the docker container 
to ensure the correct environment is present. 
A `make run` target is provided to start the container with the repository mounted to `/opt/project`.

Some additional test-data is required to run the tests that is too large to check into the repository.
This can be fetched using `make tests/test-data/mcip` (locally or in the docker container).

```
    make run
    cd /opt/project
    pytest tests
```
