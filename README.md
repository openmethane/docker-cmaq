# cmaq-container

Build immutable images for running WRF.


## Requirements

[Docker](https://www.docker.com/) is required.

    
### Local development


To build the docker image, run the following command:

```
    docker build . -t wrf
```

If building on a Arm-based Mac, the following command should be used to target the correct platform.
Otherwise, the build will emulate the amd64 platform and take much longer.

```
    docker build --build_arg PLATFORM=linux/arm64 . -t wrf
```
