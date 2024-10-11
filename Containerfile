# SPDX-FileCopyrightText: 2024 Sascha Brawer <sascha@brawer.ch>
# SPDX-License-Identifier: MIT
#
# Configuration for building the Cadaref tool inside a Linux container.
# Used in the continuous integration pipeline to make sure that Cadaref
# can be compiled and run on Alpine Linux.

FROM alpine:3.20.3
WORKDIR /home/builder
COPY . .
RUN apk add --no-cache cargo gdal-dev rust
RUN cargo build && cargo test
RUN cargo build --release && cargo test --release
