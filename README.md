# Autosomal Recessive Variant Filtering Pipeline

This repository contains an R-based workflow for prioritizing candidate variants from annotated whole-exome sequencing datasets under an autosomal recessive inheritance model.

## Overview

The pipeline applies a series of filtering steps commonly used in rare disease and family-based sequencing studies:

1. Quality filtering (QUAL, GQ, PASS).
2. Population frequency filtering.
3. Functional consequence filtering.
4. Family-based segregation filtering under an autosomal recessive model.

Optional analyses include exploration of alternative inheritance models and relaxed segregation criteria.

## Requirements

* R (>= 4.5.2)
* tidyverse

## Input

The pipeline expects an annotated variant table containing:

- Variant information
- Functional annotations
- Population frequency data
- Family genotype information

## Output

The workflow generates progressively filtered variant tables and a final set of candidate variants consistent with the selected inheritance model.

## Data Availability

Example patient data are not included in this repository. Users should provide their own annotated variant datasets.

## Authors
Mayela Fosado-Mendoza

Adrián Guevara-Maqueda

