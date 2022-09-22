import os
import sys
import argparse
import traceback
import socket
import re
import datetime

MICROHH_HOME = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(MICROHH_HOME, "external/kernel_launcher/python"))
import kernel_launcher as kl


def host_name():
    return socket.gethostname()


def device_name(device=0):
    import pycuda.driver as drv
    drv.init()
    return drv.Device(device).name()


def tune_kernel(filename, args):
    print(f"parsing file: {filename}")

    problem = kl.load_tuning_problem(filename, data_dir=args.data_dir)
    options = dict(
            iterations=args.iterations,
            verify=args.verify,
            atol=args.atol,
    )
    strategy_options = dict(
        time_limit=args.time_limit,
        max_fevals=1e99,
    )

    if all(a is None for a in problem.answers):
        for i, arg in enumerate(problem.args):
            if len(arg) > 512 * 512:
                problem.answers[i] = arg
                break

    print(f"host name: {socket.gethostname()}")
    print(f"device name: {device_name()}")
    print(f"kernel name: {problem.key}")
    print(f"problem size: {problem.problem_size}")

    if args.strategy == "block":
        block_params = dict()

        for k, v in problem.space.default_config().items():
            block_params[k] = [v]

        for expr in problem.kernel.block_size:
            for var in expr.free_variables():
                block_params[var] = problem.space.params[var]

        before = datetime.datetime.now()

        results, env = problem.tune(
            block_params,
            strategy="brute_force",
            strategy_options=dict(max_fevals=1e99),
            **options)

        after = datetime.datetime.now()
        time_remaining = args.time_limit - (after - before).total_seconds()

        if time_remaining > 0:
            best_result = min(results, key=lambda p: p["time"])

            more_params = dict(problem.space.params)
            for expr in problem.kernel.block_size:
                for var in expr.free_variables():
                    more_params[var] = [best_result[var]]

            strategy_options["time_limit"] = time_remaining
            more_results, _ = problem.tune(
                more_params,
                strategy="bayes_opt",
                strategy_options=strategy_options,
                **options)

            results += more_results

    elif args.strategy == "random":
        results, env = problem.tune(
            strategy="random",
            strategy_options=strategy_options,
            **options)

    elif args.strategy == "bayes":
        results, env = problem.tune(
            strategy="bayes_opt",
            strategy_options=strategy_options,
            **options)

    else:
        raise ValueError(f"unknown strategy: {args.strategy}")

    best_result = min(results, key=lambda p: p["time"])

    print(f"finished tuning {problem.key}")
    print(f"best configuration: {best_result!r}")

    print("writing wisdom file")
    kl.write_wisdom_for_problem(args.output, problem, results, env,
                                max_results=1, merge_existing_results=args.append)


def parse_time(input):
    match = re.match("^([0-9]+([.][0-9]*)?)$", input)
    if match:
        return float(match[1])

    match = re.match("^([0-9]+):([0-9]+)$", input)
    if match:
        return int(match[1]) * 60 + int(match[2])

    match = re.match("^([0-9]+):([0-9]+):([0-9]+)$", input)
    if match:
        return int(match[1]) * 60 * 60 + int(match[1]) * 60 + int(match[2])

    raise ValueError(f"failed to parse time: {input}")


def main():
    parser = argparse.ArgumentParser(
            description="Tune given kernel files and store the results in wisdom files")
    parser.add_argument("--iterations", "-i", type=int, default=5,
                        help="Number of benchmark iterations for each kernel")
    parser.add_argument("--time", "-t", type=parse_time, default="15:00", dest="time_limit",
                        help="Maximum time in seconds spend on tuning each kernel.")
    parser.add_argument("--output", "-o", default=os.path.join(MICROHH_HOME, "wisdom"),
                        help="Directory where to store resulting wisdom files.")
    parser.add_argument("--data-dir", "-d", default=None,
                        help="Directory where data files (.bin) are located.")
    parser.add_argument("--no-verify", action="store_false", default=True, dest="verify",
                        help="Verify if the output of each kernel launch is correct.")
    parser.add_argument("--tolerance", "--atol", dest="atol", type=float,
                        help="Absolute tolerance used for verification as interpreted by numpy.isclose.")
    parser.add_argument("--append", "-a", default=False, action="store_true",
                        help="Append new results to existing wisdom instead of overwriting them.")
    parser.add_argument("--strategy", "-s", default="bayes", choices=["block", "bayes", "random"],
                        help="The strategy to use for tuning:\n"
                             " - random: try random configurations until time runs out.\n"
                             " - bayes: use Bayesian optimization to try configurations until time runs out.\n"
                             " - block: brute-force search block sizes and then optimize the remaining parameters.\n")
    parser.add_argument("files", nargs="*")

    args = parser.parse_args()

    if not os.path.isdir(args.output):
        print(f"error: not a valid directory: {args.output}")
        return

    if not args.files:
        print(f"error: no files given")

    for file in args.files:
        try:
            tune_kernel(file, args)
        except Exception as e:
            print(f"error: exception occurred while tuning {file}:")
            traceback.print_exception(type(e), e, e.__traceback__)
            print()



if __name__ == "__main__":
    main()
