def logging(log_file, string):
    with open(log_file, "a") as f:
        f.write(string)
