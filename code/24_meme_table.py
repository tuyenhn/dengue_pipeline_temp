import sys
import os


def main(fpath):
    print(fpath)

    with open(fpath, mode="r") as f:
        dat = f.readlines()
        dat = list(filter(lambda x: x.startswith("|"), dat))
        print("".join(dat))


if __name__ == "__main__":
    log_fpath = os.path.join(os.getcwd(), sys.argv[1])
    main(log_fpath)
