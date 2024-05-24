FROM ubuntu:20.04

MAINTAINER Jared Lewis <jared.lewis@climate-resource.com>

ARG DEBIAN_FRONTEND=noninteractive
ENV CMAQ_VERSION="5.0.2"
ENV TZ=Etc/UTC

RUN apt-get update && \
    apt-get install -y build-essential gfortran m4 csh git jq wget libopenmpi-dev unzip && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/cmaq

COPY --from=ghcr.io/climate-resource/wrf:4.5.1 /opt/wrf/libs /opt/wrf/libs

# Build IOAPI
COPY templates/ioapi /opt/cmaq/templates/ioapi
COPY scripts/common.sh /opt/cmaq/scripts/common.sh
COPY scripts/build_00_ioapi.sh /opt/cmaq/scripts/build_00_ioapi.sh
RUN bash /opt/cmaq/scripts/build_00_ioapi.sh

# Build common CMAQ dependencies
COPY src/CMAQv5.0.2_notpollen /opt/cmaq/CMAQv5.0.2_notpollen
COPY scripts/build_10_cmaq.sh /opt/cmaq/scripts/build_10_cmaq.sh
RUN bash /opt/cmaq/scripts/build_10_cmaq.sh

# Build the CMAQ adjoint
COPY scripts/build_11_cmaq_adj.sh /opt/cmaq/scripts/build_11_cmaq_adj.sh
COPY src/cmaq_adj /opt/cmaq/cmaq_adj
RUN bash /opt/cmaq/scripts/build_11_cmaq_adj.sh

ENTRYPOINT ["/bin/bash"]