# Weisman - Klemp case

This case is the deep convection experiment from Weisman & Klemp (1982, MWR) in which a warm bubble is released, which splits and forms a left and right moving storm. The case has now warm micro enabled and will serve as a benchmark case for deep convection once the full S&B ice microphysics is completed.

## Instructions
1. Run `python3 weisman_klemp_input.py`
1. Run `./microhh init weisman_klemp`
1. Run `python3 bubble_theta.py`
1. Run `./microhh run weisman_klemp`
