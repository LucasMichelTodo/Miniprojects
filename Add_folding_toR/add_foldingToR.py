import sys

def add_folding(rfile):
    with open(rfile, "r+") as infile:
        for line in infile:
            if line.startswith("##"):
                print("##"+line.strip()+" ####")
            else:
                print(line.strip())

if __name__ == "__main__":
        filename = sys.argv[1]
        add_folding(filename)
