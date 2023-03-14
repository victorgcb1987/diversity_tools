
def check_results(program, results):
    if results["return_code"] == 0:
        print("{} successfully run".format(program))
    else:
        print("{} failed".format(program))
        raise RuntimeError(results["log_messages"])
    