
def check_results(program, results):
    if results["return_code"] == 0:
        print("{} successfully run".format(program))
    else:
        print("{} failed".format(program))
        raise RuntimeError(results["log_messages"])
    
def file_exists(filepath):
    if filepath.exists():
        if not filepath.is_dir() and filepath.stat().st_size > 0:
            return True
    return False

def store_results(out_fpath, run_results):
    return {"out_fpath": out_fpath,
            "return_code": run_results.returncode,
            "log_messages": run_results.stderr.decode()}
    