FROM debian:buster-slim as base

###
# Install certs
###

RUN apt-get update \
      && apt-get -y install \
      ca-certificates \
      && rm -rf /var/lib/apt/lists/*

FROM base as builder

###
# Build the binary. These intermediate layers will be discarded.
###

LABEL description="Build layers - crux-toolkit"

# Required system packages
RUN apt-get update && apt-get -y install \
  build-essential \
  cmake \
  curl \
  git \
  subversion \
  wget \
  libcurl4-openssl-dev \
  libssl-dev \
  uuid-dev \
  zlib1g-dev \
  libpulse-dev

RUN mkdir /app

WORKDIR /src

RUN git config --global user.email "you@example.com" && \
    git config --global user.name "Your Name"

COPY . /src/crux-toolkit

# Next build crux
###

WORKDIR /src/crux-toolkit

RUN cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX:PATH=/app .

RUN make && make install

# Gather binary and dynamic dependencies to copy over to base layer
###

RUN mkdir /gathered

RUN ldd /app/bin/crux | grep "=> /" | awk '{print $3}' |  xargs -I '{}' sh -c "cp -v --parents {} /gathered"

RUN cp --parents /app/bin/crux /gathered

FROM base as runtime

###
# Final image contains binary, dynamic dependencies, and SSL/TLS certs
###

LABEL description="crux-toolkit"

COPY --from=builder /gathered /

ENV PATH="${PATH}:/app/bin"

WORKDIR /app

