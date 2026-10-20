# Squall line case

This case is the deep convection experiment from Tijhuis et al. (2025, https://arxiv.org/abs/2506.19435) in which a cold pool is initiated, which forces the formation of a squall line. The case serves as a test case for deep convection with the full S&B ice microphysics, but can also be run with S&B warm microphysics, nsw6 microphysics, or without any microphysics scheme.

## Instructions
1. Run `python3 squall_line_input.py`
1. Run `./microhh init squall_line`
1. Run `python3 line_theta.py`
1. Run `./microhh run squall_line`
