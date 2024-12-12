import argparse

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("add", type=str,nargs='?')

    args = parser.parse_args()
    print(args.add)
    pass


if __name__ == "__main__":
    main()


