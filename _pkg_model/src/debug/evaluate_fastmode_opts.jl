## Imports
using Helios
using BenchmarkTools

## Description
# Script to test "fast-mode" features and their performance consequences

# Test with second power-law branch disabled
@btime sim_gfc2023_dict(zeros(5000), 1000.0; powerlaw_include_fast=false);  # => 84.489 ms

# Test with second power-law branch enabled
@btime sim_gfc2023_dict(zeros(5000), 1000.0; powerlaw_include_fast=true, powerlaw_len_memory=100000);  # => 168.918 ms

# Test with second power-law branch enabled but time limited
@btime sim_gfc2023_dict(zeros(5000), 1000.0; powerlaw_include_fast=true, powerlaw_len_memory=5000);  # => 168.918 ms

# Test with second power-law branch enabled but time limited
@btime sim_gfc2023_dict(zeros(5000), 1000.0; powerlaw_include_fast=true, powerlaw_len_memory=500);  # => 168.918 ms