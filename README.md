# Robust Performance Analysis for Time-Varying Multi-Agent Systems with Stochastic Packet Loss

## General

This repository contains an implementation of the algorithms described in the paper

> C. Hespe, H. Werner, *Robust Performance Analysis for Time-Varying Multi-Agent Systems with Stochastic Packet Loss*, submitted to IFAC World Congress 2023.

It may be used to recreate and validate the figures from the paper.
To do so, run either of the two main entry points in the repository, the scripts `eigenvalue_contour.m` and `varying_graphs.m`.
Be advised that the latter of these scripts has a runtime of at least one day.
The raw data used in the figures in the paper is available in the subdirectory `figures`.

## Simulation of Multi-Agent Systems

Both of the main entry points will run simulations of multi-agent systems with packet loss.
For that purpose, an open source library which can be found [on Github](https://github.com/TUHH-ICS/MAS-Simulation) is utilized.
The repository is included as a git submodule and needs to be initialized using `git submodule update --init`.

## Prerequisites

To run the scripts in this repository, you will need a working copy of [*Yalmip*](https://yalmip.github.io/) together with a suitable SDP solver in your *Matlab* path.

The code in this repository was tested in the following environment:

* *Windows 10* Version 21H2
* *Matlab* 2021b
* *Yalmip* 31-March-2021

The *Matlab* [`parfor`](https://de.mathworks.com/help/parallel-computing/parfor.html) feature from the *Parallel Computing Toolbox* is used to speed up the calculations.
*Matlab* should automatically detect if that toolbox is not available and run the iterations sequentially in that case.
However, this will drastically prolong the runtime of the scripts!
You may want to reduce the number of sampling points for the figures or run the calculations for smaller networks.
