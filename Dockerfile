FROM debian:bookworm as build

ARG DEBIAN_FRONTEND=noninteractive
ENV CMAQ_VERSION="5.0.2"

RUN apt-get update && \
    apt-get install -y build-essential m4 csh wget libhdf5-dev libmpich-dev libnetcdff-dev gfortran python3-pytest python3-dotenv

RUN wget -nv https://snapshot.debian.org/archive/debian-security/20190413T190428Z/pool/updates/main/j/jasper/libjasper1_1.900.1-debian1-2.4%2Bdeb8u6_amd64.deb && \
    dpkg --force-all -i libjasper1_1.900.1-debian1-2.4+deb8u6_amd64.deb && \
    wget -nv https://snapshot.debian.org/archive/debian-security/20190413T190428Z/pool/updates/main/j/jasper/libjasper-dev_1.900.1-debian1-2.4%2Bdeb8u6_amd64.deb && \
    dpkg -i libjasper-dev_1.900.1-debian1-2.4+deb8u6_amd64.deb

WORKDIR /opt/cmaq

# Build ioapi
COPY templates/ioapi /opt/cmaq/templates/ioapi
COPY scripts/common.sh /opt/cmaq/scripts/common.sh
COPY scripts/build_00_ioapi.sh /opt/cmaq/scripts/build_00_ioapi.sh
RUN bash /opt/cmaq/scripts/build_00_ioapi.sh

# Build a modified version of CMAQ in ch4 only mode
COPY templates/cmaq /opt/cmaq/templates/cmaq
COPY src/CMAQv5.0.2_notpollen /opt/cmaq/CMAQv5.0.2_notpollen
COPY scripts/build_10_cmaq.sh /opt/cmaq/scripts/build_10_cmaq.sh
RUN bash /opt/cmaq/scripts/build_10_cmaq.sh

# Build the CMAQ adjoint
COPY scripts/build_11_cmaq_adj.sh /opt/cmaq/scripts/build_11_cmaq_adj.sh
COPY src/cmaq_adj /opt/cmaq/cmaq_adj
RUN bash /opt/cmaq/scripts/build_11_cmaq_adj.sh


FROM debian:bookworm AS runtime

MAINTAINER Jared Lewis <jared.lewis@climate-resource.com>

ENV TZ=Etc/UTC
ENV CMAQ_VERSION="5.0.2"
ENV DEBUGINFOD_URLS="https://debuginfod.debian.net"

WORKDIR /opt/cmaq
COPY --from=build /opt/cmaq /opt/cmaq

# todo we don't need the -dev libraries for the runtime, could be slimmed down
RUN apt-get update && \
    apt-get install -y make csh wget libhdf5-dev libmpich-dev libnetcdff-dev python3-pytest python3-dotenv && \
    rm -rf /var/lib/apt/lists/*

# todo we probably don't need the -dev library for the runtime
RUN wget -nv https://snapshot.debian.org/archive/debian-security/20190413T190428Z/pool/updates/main/j/jasper/libjasper1_1.900.1-debian1-2.4%2Bdeb8u6_amd64.deb && \
    dpkg --force-all -i libjasper1_1.900.1-debian1-2.4+deb8u6_amd64.deb && \
    wget -nv https://snapshot.debian.org/archive/debian-security/20190413T190428Z/pool/updates/main/j/jasper/libjasper-dev_1.900.1-debian1-2.4%2Bdeb8u6_amd64.deb && \
    dpkg -i libjasper-dev_1.900.1-debian1-2.4+deb8u6_amd64.deb && \
    rm libjasper1_1.900.1-debian1-2.4+deb8u6_amd64.deb libjasper-dev_1.900.1-debian1-2.4+deb8u6_amd64.deb

ENTRYPOINT ["/bin/bash"]
